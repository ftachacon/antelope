/// SBEs with 'L'ength gauge, 'W'annier approach, 'M'oving frame 
/**
 * @author Dasol Kim
*/

#pragma once

#include <iostream>

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <libconfig.h++>

#include "momaxis.h"
#include "laser.h"
#include "constant.h"
#include "utility.h"

#include "system_list.h"

using namespace std;

class SBEsLWM
{
public:
    complex **dmatrix; 

    WannierMaterial *material;  ///< dipole&Hamiltonian generator

    // allocated memory with no value - to reduce costs for making temporoary variables
    complex *hamiltonian;
    complex **dipoleMatrix, **jMatrix;
    complex *newdMatrix;
    complex *temp1Matrix, *temp2Matrix;
    complex **rkMatrix;

    complex *uMatrix, *ctransuMatrix;
    double *dephasingMatrix;

    laser *fpulses;

    int Nband;

    //bool isWannier90;       ///< if true, it means target material is written by wannier90 input.
    //std::array<std::array<double, Ndim>, Ndim> vec_lattice; // only for wannier90 case, to project on reciprocal vectors

    //array<double, Ndim*Ndim> BzAxes;
    //array<double, Ndim> BzOrigin;

    // variables for Runge-Kutta method
    int RKorder, NumRK; // order of accumulated error in RK, and number of stages in RK algorithm
    double *tRK;
    double *bRK;
    double **aRK;

    momaxis *kmesh;
    string targetMaterial;
    array<double, Ndim> ksfactor;
    int shotNumber;
    int diagnostic;
private:
    double dt;
    array<int, Ndim> Nk;

public:
    double dephasing_time_inter, dephasing_time_intra;          // dephasing time, intrabnad dephasing has to be changed later
    double dephasing_factor_inter, dephasing_factor_intra;      // dephasing factor(1/dephasing time), to avoid expensive float division operation and remove if statement

    bool isDipoleZero;                                  ///< default is true. Wannier dipole is non-zero in very rare cases.

    SBEsLWM(const libconfig::Setting * _cfg);
    ~SBEsLWM();
    
    void RunSBEs(int _kindex, double _time);

    array<double, Ndim> GenCurrent(int _kindex, double _time);

    void GenDifferentialDM(complex *out, complex *input, int _kindex, double _time);

    void InitializeGeneral(const libconfig::Setting *_calc);     ///< Initialize general stuff inside constructor
    void InitializeLaser(const libconfig::Setting *_laser);     ///< Initialize laser inside constructor

    array<double, Ndim> GenKpulsA(array<double, Ndim> _kpoint, double time);

};

SBEsLWM::SBEsLWM(const libconfig::Setting * _cfg)
{
    const libconfig::Setting &cfg = (*_cfg);
    InitializeGeneral( &cfg["calc"]);
    InitializeLaser( &cfg["laser"]);
    // Check target, get Nband
    try
    {
        targetMaterial = cfg.lookup("target").c_str();
        //Nband = cfg[targetMaterial.c_str()].lookup("Nband");
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        cerr << "No 'target' in inputParam.cfg or no 'Nband' parameter" << endl; 
        exit(EXIT_FAILURE);
    }
    // cout << "Target :  " << targetMaterial << endl;

    // Initialize specific material
    //isWannier90 = false;
    if (targetMaterial == "Haldane")
    {
        material = new Haldane( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "WilsonMass")
    {
        material = new WilsonMass( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "Bi2Se3surf")
    {
        material = new BieSe3surf( &cfg );
    }
    else if (targetMaterial == "TMDC")
    {
        material = new TMDC( &cfg[targetMaterial.c_str()] );
    }
    else
    {
        cerr << "Undefined Material\n";
        exit(EXIT_FAILURE);
        //material = new Wannier90( &(cfg[targetMaterial.c_str()]) );
        //isDipoleZero = dynamic_cast<Wannier90*>(material)->isDipoleZero;
        //vec_lattice = dynamic_cast<Wannier90*>(material)->vec_lattice;
        //isWannier90 = true;
    }
    
    Nband = material->Nband;

    auto [bzaxes, bzori] = material->GenBrillouinzone();

    // generate k-grid
    kmesh = new momaxis( Nk, bzaxes, bzori );

    /*if (isWannier90)
    {
        dynamic_cast<Wannier90*>(material)->CalculateKMesh(kmesh);
    }*/
    

    dmatrix = Create2D<complex>(this->kmesh->Ntotal, Nband*Nband);

    hamiltonian = new complex[Nband * Nband];
    temp1Matrix = new complex[Nband * Nband];
    temp2Matrix = new complex[Nband * Nband];
    newdMatrix = new complex[Nband * Nband];

    jMatrix = Create2D<complex>(Ndim, Nband*Nband);
    dipoleMatrix = Create2D<complex>(Ndim, Nband*Nband);
    uMatrix = new complex[Nband * Nband];
    ctransuMatrix = new complex[Nband * Nband];
    dephasingMatrix = new double[Nband * Nband];

    for (int k = 0; k < this->kmesh->Ntotal; ++k)
    {
        material->GenInitialValue(&dmatrix[k][0], kmesh->kgrid[k]);
    }

    // allocation for rk realted variables
    tRK = new double[NumRK-1];
    bRK = new double[NumRK];
    aRK = new double*[NumRK-1];

    rkMatrix = Create2D<complex>(NumRK, Nband * Nband);

    for (int i = 0; i < NumRK-1; ++i)
    {
        aRK[i] = new double[i+1];
    }

    // Explicit Runge-Kutta method:
    // y_{n+1} = y_{n} + dt * sum_{i=0}^{N-1} b_{i}*k_{i}
    // k0 = f(tn, yn)
    // k1 = f(tn + tRK[0]*dt, yn + dt*aRK[0][0]*k0)
    // ...
    // km = f(tn + tRK[m-1]*dt, yn + dt* sum_{i=0}^{m-1}aRK[m-1][i] * ki)
    switch (RKorder)
    {
        case 4:
            tRK[0] = 0.5;   tRK[1] = 0.5;   tRK[2] = 1.0;
            bRK[0] = 1.0/6.0;   bRK[1] = 2.0/6.0;   bRK[2] = 2.0/6.0;   bRK[3] = 1.0/6.0;
            aRK[0][0] = 0.5;
            aRK[1][0] = 0.0;    aRK[1][1] = 0.5;
            aRK[2][0] = 0.0;    aRK[2][1] = 0.0;    aRK[2][2] = 1.0;
            break;
        case 5:
            tRK[0] = 0.25;  tRK[1] = 3./8.; tRK[2] = 12./13.;   tRK[3] = 1.0;   tRK[4] = 0.5;
            bRK[0] = 16./135.;  bRK[1] = 0; bRK[2] = 6656./12825.;  bRK[3] = 28561./56430.; bRK[4] = -9./50.;    bRK[5] = 2./55.;
            aRK[0][0] = 1.0/4.0;
            aRK[1][0] = 3./32.; aRK[1][1] = 9./32.;
            aRK[2][0] = 1932./2197.;    aRK[2][1] = -7200./2197.;   aRK[2][2] = 7296./2197.;
            aRK[3][0] = 439./216.;  aRK[3][1] = -8.;    aRK[3][2] = 3680./513.; aRK[3][3] = -845./4104.;
            aRK[4][0] = -8./27.;    aRK[4][1] = 2.; aRK[4][2] = -3544./2565; aRK[4][3] = 1859./4104; aRK[4][4] = -11./40.;
            break;
        default:
            throw out_of_range("Undefined Runge-Kutta order");
            break;
    }

    // dephasing matrix
    for (int i = 0; i < Nband; ++i)
    {
        for (int j = 0; j < Nband; ++j)
        {
            if (i == j)
                //dephasingMatrix[i*Nband + j] = dephasing_factor_intra;
                dephasingMatrix[i*Nband + j] = 0.;
            else
                dephasingMatrix[i*Nband + j] = dephasing_factor_inter;
        }
    }
}

void SBEsLWM::InitializeGeneral(const libconfig::Setting *_calc)
{
    // default values of parameters
    isDipoleZero = true;
    fill(ksfactor.begin(), ksfactor.end(), 1.0);
    RKorder = 5;
    dephasing_time_inter = -1.;
    dephasing_time_intra = -1.;

    const libconfig::Setting &calc = (*_calc);
    fill(Nk.begin(), Nk.end(), 1);
    // Get calculation parameters
    try
    {
        const libconfig::Setting &npoints_setting = calc["Npoints"];
        for (int i = 0; i < npoints_setting.getLength(); ++i)
        {
            Nk[i] = npoints_setting[i];
        }
        dt = calc.lookup("dt");

        // read parameters which have default values - no error when they are omitted
        calc.lookupValue("T1", dephasing_time_intra );
        calc.lookupValue("T2", dephasing_time_inter );
        calc.lookupValue("shotNumber", shotNumber);
        calc.lookupValue("diagnostic", diagnostic);
        calc.lookupValue("RKorder", RKorder);
        // ksfactor - value & array
        const libconfig::Setting &ksfactor_setting = calc["ksfactor"];
        for (int i = 0; i < ksfactor_setting.getLength(); ++i)
        {
            ksfactor[i] = ksfactor_setting[i];
        }
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        cerr << "Setting not found while search calc group \n";
        exit(EXIT_FAILURE);
    }
    catch(const libconfig::SettingTypeException &tex )
    {
        cerr << "Type exception in calc group\n";
        exit(EXIT_FAILURE);
    }

    // some post-process
    double eps = 1.0e-16;
    if (dephasing_time_inter <= eps)
        dephasing_factor_inter = 0.;
    else
        dephasing_factor_inter = 1.0 / dephasing_time_inter;
    
    if (dephasing_time_intra <= eps)
        dephasing_factor_intra = 0.;
    else
        dephasing_factor_intra = 1.0 / dephasing_time_intra;
    
    // Runge-Kutta order matching
    switch (RKorder)
    {
        case 4:
            NumRK = 4;
            break;
        case 5:
            NumRK = 6;
            break;
        default:
            cout << "Not-implemented RK order: " << RKorder << "\n";
            exit(1);
            break;
    }
}

void SBEsLWM::InitializeLaser(const libconfig::Setting *_laser)
{
    const libconfig::Setting &laserCfg = (*_laser);

    double offset_before = 0.;
    double offset_after = 0.;
    try
    {
        const libconfig::Setting &pulses = laserCfg["pulses"];
        int pulseNum = pulses.getLength();
        fpulses = new laser();

        double t_E0, t_ellip, t_w0, t_ncycles, t_cep, t_t0, t_phix, t_thetaz, t_phiz;
        string t_env_name;
        for (int i = 0; i < pulseNum; ++i)
        {
            if ( pulses[i].lookupValue("E0", t_E0) && pulses[i].lookupValue("ellip", t_ellip) && pulses[i].lookupValue("w0", t_w0) 
                && pulses[i].lookupValue("ncycles", t_ncycles) && pulses[i].lookupValue("cep", t_cep) && pulses[i].lookupValue("t0", t_t0) 
                && pulses[i].lookupValue("phix", t_phix) && pulses[i].lookupValue("env_name", t_env_name) )
            {
                t_thetaz = 0.;  t_phiz = 0.;    // this parameters can be omitted from config file
                pulses[i].lookupValue("thetaz", t_thetaz);
                pulses[i].lookupValue("phiz", t_phiz);
                t_phix *= pi/180.0;
                t_thetaz *= pi/180.0;
                t_phiz *= pi/180.0;
                fpulses->pulses.push_back(Pulse(t_E0, t_w0, t_ellip, t_ncycles, t_cep, t_t0, t_phix, t_env_name, t_thetaz, t_phiz));
            }
            else
            {
                cerr << "Error while initilizng laser pulses. Check config file\n";
                exit(EXIT_FAILURE);
            }
        }
        if ( laserCfg.exists("offset") )
        {
            offset_before = laserCfg["offset"][0];
            offset_after = laserCfg["offset"][1];
        }
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        cerr << "No laser, check config file \n";
        exit(EXIT_FAILURE);
    }
    fpulses->Initialize(dt, offset_before, offset_after);
}

SBEsLWM::~SBEsLWM()
{
    Delete2D<complex>(dmatrix, this->kmesh->Ntotal, Nband*Nband);

    delete[] hamiltonian;
    delete[] newdMatrix;
    delete[] temp1Matrix, temp2Matrix;
    delete[] dephasingMatrix;
    delete[] uMatrix, ctransuMatrix;
    Delete2D<complex>(rkMatrix, NumRK, Nband*Nband);
    Delete2D<complex>(jMatrix, Ndim, Nband*Nband);
    Delete2D<complex>(dipoleMatrix, Ndim, Nband*Nband);

    delete[] tRK;
    delete[] bRK;
    for (int i = 0; i < NumRK-1; ++i)
    {
        delete[] aRK[i];
    }
    delete[] aRK;

    delete material;
    delete fpulses;

    delete kmesh;
}

void SBEsLWM::RunSBEs(int _kindex, double  _time)
{
    GenDifferentialDM( rkMatrix[0], dmatrix[_kindex], _kindex, _time);
    for (int irk = 0; irk < NumRK-1; ++irk)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                newdMatrix[m*Nband + n] = dmatrix[_kindex][m*Nband + n];
                for (int jrk = 0; jrk < irk+1; ++jrk)
                {
                    newdMatrix[m*Nband + n] += fpulses->dt*aRK[irk][jrk] * rkMatrix[jrk][m*Nband + n];
                }
            }
        }
        GenDifferentialDM(rkMatrix[irk+1], newdMatrix, _kindex, _time + tRK[irk]*fpulses->dt);
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            for (int irk = 0; irk < NumRK; ++irk)
            {
                dmatrix[_kindex][m*Nband + n] += fpulses->dt*bRK[irk] * rkMatrix[irk][m*Nband + n];
            }
        }
    }
}
array<double, Ndim> SBEsLWM::GenCurrent(int _kindex, double _time)
{
    auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    material->GenJMatrix(jMatrix, _tkp);

    // J = Tr{densitymatrix j}
    array<double, Ndim> outCurrent;
    fill(outCurrent.begin(), outCurrent.end(), 0.);
    for (int iaxis = 0; iaxis < Ndim; ++iaxis)
    {
        for (int i = 0; i < Nband; ++i)
        {
            for (int j = 0; j < Nband; ++j)
            {
                outCurrent[iaxis] += real( dmatrix[_kindex][Nband*i + j] * jMatrix[iaxis][Nband*j + i] );
            }
        }
    }
    return outCurrent;
}


void SBEsLWM::GenDifferentialDM(complex *out, complex *input, int _kindex, double _time)
{
    array<double, Ndim> evector = fpulses->elaser(_time);
    auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    material->GenHamiltonian(hamiltonian, _tkp);

    //fill(out, out+Nband*Nband, 0.);
    MatrixMult(out, hamiltonian, input, Nband);
    /*for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            for (int l = 0; l < Nband; ++l)
            {
                out[m*Nband + n] += hamiltonian[m*Nband + l] * input[l*Nband + n];
            }
        }
    }*/
    if (!isDipoleZero)
    {
        material->GenDipole(dipoleMatrix, _tkp);
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int l = 0; l < Nband; ++l)
                {
                    out[m*Nband + n] += (evector[0]*dipoleMatrix[0][m*Nband + l] 
                                        + evector[1]*dipoleMatrix[1][m*Nband + l]
                                        + evector[2]*dipoleMatrix[2][m*Nband + l]) * input[l*Nband + n];
                }
            }
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            hamiltonian[m*Nband + n] = out[m*Nband + n];
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            out[m*Nband + n] -= conj(hamiltonian[n*Nband + m]);
            out[m*Nband + n] *= -I;
        }
    }
    // Add dephainsg process here!
    material->GenUMatrix(uMatrix, _tkp);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            ctransuMatrix[m*Nband + n] = conj(uMatrix[n*Nband + m]);
        }
    }
    MatrixMult(temp1Matrix, input, uMatrix, Nband);
    MatrixMult(temp2Matrix, ctransuMatrix, temp1Matrix, Nband);

    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            temp2Matrix[m*Nband + n] *= dephasingMatrix[m*Nband + n];
        }
    }

    MatrixMult(temp1Matrix, uMatrix, temp2Matrix, Nband);
    MatrixMult(temp2Matrix, temp1Matrix, ctransuMatrix, Nband);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            out[m*Nband + n] -= temp2Matrix[m*Nband + n];
        }
    }
}

array<double, Ndim> SBEsLWM::GenKpulsA(array<double, Ndim> _kpoint, double time)
{
    auto avector = fpulses->avlaser(time);
    array<double, Ndim> tkp = {_kpoint[0]+avector[0], _kpoint[1]+avector[1], _kpoint[2]+avector[2]};
    /*if (isWannier90)
    {
        // k = sum_i g_i bvec_i
        // k.avec_i = 2pi g_i
        array<double, Ndim> returnVal;
        fill(returnVal.begin(), returnVal.end(), 0.0);
        for (int i = 0; i < Ndim; ++i)
        {
            for (int j =0; j < Ndim; ++j)
            {
                returnVal[i] += vec_lattice[i][j] * tkp[j];
            }
            returnVal[i] /= (2*pi);
        }
        return returnVal;
    }
    else
    {
        return tkp;
    }*/
    return tkp;
}