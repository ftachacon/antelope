// cout// SBEs with 'L'ength gauge, 'W'annier approach, 'M'oving frame
/**
 * @author Dasol Kim
*/

#pragma once

#include <iostream>
#include <fstream>
#include <tuple>

#include <string>
#include <list>
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
typedef std::tuple<array<double, Ndim>,  array< array<double, Ndim>,Nband>,  array<array<double, Ndim>, Ndim> > tupleSpinCharge;

// public static class Tuple;
class SBEs
{
public:
    const GaugeType gauge;
    InitialValueType initType;

    complex **dmatrix;
    complex **initMatrix;

    BaseMaterial *material;  ///< dipole&Hamiltonian generator

    // allocated memory with no value - to reduce costs for making temporoary variables
    complex *hamiltonian;
    complex **dipoleMatrix, **jMatrix;
    complex *newdMatrix;
    complex *temp1Matrix, *temp2Matrix;
    complex **rkMatrix;

    complex *uMatrix, *ctransuMatrix;
    double *dephasingMatrix;

    complex ***pMatrix;
    double **edispersion;
    complex ***sdMatrix;

    laser *fpulses;

    int Nband;
    array<int, Ndim> Nk;

    double Ef;          ///< Fermi energy
    double thermalE;    ///< Thermal energy (kT) in a.u.

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
    double sigma[Ndim][2];
    int shotNumber;
    int diagnostic;

    double rx1;
    double ry1;
    double rz1;


    double rx2;
    double ry2;
    double rz2;

    double rx3;
    double ry3;
    double rz3;

    double rx4;
    double ry4;
    double rz4;

    double rmax;

private:
    double dt;

    // used for lapack routines
    int num_of_eig;
    complex *tempHmat;
    double *tempEval;
    int *tempIsuppz;
    double eps;


public:
    double dephasing_time_inter, dephasing_time_intra;          // dephasing time, intrabnad dephasing has to be changed later
    double dephasing_factor_inter, dephasing_factor_intra;      // dephasing factor(1/dephasing time), to avoid expensive float division operation and remove if statement

    bool isDipoleZero;                                  ///< default is true. Wannier dipole is non-zero in very rare cases.

    SBEs(const libconfig::Setting * _cfg, GaugeType _gauge);
    ~SBEs();

    void RunSBEs(int _kindex, double _time);

    // array<double, Ndim> GenCurrent(int _kindex, double _time);
    tupleSpinCharge  GenCurrent(int _kindex, double _time);

    void GenDifferentialDM(complex *out, complex *input, int _kindex, double _time);

    void InitializeGeneral(const libconfig::Setting *_calc);     ///< Initialize general stuff inside constructor
    void InitializeLaser(const libconfig::Setting *_laser);     ///< Initialize laser inside constructor

    void GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint);    ///< Generate U matrix

    array<double, Ndim> GenKpulsA(array<double, Ndim> _kpoint, double time);

    void WannierToHamiltonian(complex *out, complex *input, int _kindex, double _time); ///< Convert Wannier representation to Hamiltonian representation

    void PrintInfo();   ///< Print calculation information

};

SBEs::SBEs(const libconfig::Setting * _cfg, GaugeType _gauge) : gauge(_gauge)
{
    const libconfig::Setting &cfg = (*_cfg);
    InitializeGeneral( &cfg["calc"]);
    InitializeLaser( &cfg["laser"]);
    // Check target, get Nband
    // cout << "debug.inside.sbe1\n";
    try
    {
        //targetMaterial = cfg.lookup("target").c_str();
        if (!cfg.lookupValue("target", targetMaterial))
        {
            cerr << "No 'target' in input file\n";
            exit(EXIT_FAILURE);
        }
        //Nband = cfg[targetMaterial.c_str()].lookup("Nband");
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        cerr << "No 'target' in inputParam.cfg or no 'Nband' parameter" << endl;
        exit(EXIT_FAILURE);
    }
    //cout << "Target :  " << targetMaterial << endl;

    bool usingCustomInitial = false;
    if (cfg.exists(targetMaterial))
        cfg[targetMaterial.c_str()].lookupValue("UsingCustomInitial", usingCustomInitial);
    if ( usingCustomInitial )
    {
        initType = InitialValueType::Custom;
    }
    else if ( cfg.exists(targetMaterial) && cfg[targetMaterial.c_str()].lookupValue("Ef", Ef) )
    {
        initType = InitialValueType::FermiDirac;
    }
    else
    {
        initType = InitialValueType::UniformValence;
    }

    // Initialize specific material
    //isWannier90 = false;
    // cout << "debug.inside.sbe2\n";
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
    else if (targetMaterial == "TMDCs")
    {
        material = new TMDCs( &cfg[targetMaterial.c_str()] );
    }
    else if (targetMaterial == "WDSM_TaAs")
    {
        // cout << "debug.wdsm_taas\n";
        material = new WDSM_TaAs( &cfg[targetMaterial.c_str()] );
        // cout << "debug.inside.sbe3a\n";
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
    // cout << "debug.inside.sbe4\n";
    Nband = material->Nband;

    auto [bzaxes, bzori] = material->GenBrillouinzone();
    // cout << "debug.inside.sbe5\n";
    // cout << bzaxes[0] ;
    // apply ksfactor
    bzaxes[0] *= ksfactor[0];   bzaxes[1] *= ksfactor[0];   bzaxes[2] *= ksfactor[0];
    bzaxes[3] *= ksfactor[1];   bzaxes[4] *= ksfactor[1];   bzaxes[5] *= ksfactor[1];
    bzaxes[6] *= ksfactor[2];   bzaxes[7] *= ksfactor[2];   bzaxes[8] *= ksfactor[2];
    for (int tmpi =0; tmpi<Ndim;++tmpi)
    {
    sigma[tmpi][0]=0;
    sigma[tmpi][1]=1;
    sigma[tmpi][2]=1;
    sigma[tmpi][3]=0;
    };

    // generate k-grid
    kmesh = new momaxis( Nk, bzaxes, bzori );
    // for (int i = 0; i < 3; i++)
    // {
    //   for (int j = 0; j < 3; j++)
    //   {
    //   cout << kmesh[i][j] << " ";
    //   }
    // }

   // // Newline for new row
   // cout << endl;
//    cout << "debug.inside.sbe4\n";

    /*if (isWannier90)
    {
        dynamic_cast<Wannier90*>(material)->CalculateKMesh(kmesh);
    }*/


    dmatrix = Create2D<complex>(this->kmesh->Ntotal, Nband*Nband);
    initMatrix = Create2D<complex>(this->kmesh->Ntotal, Nband*Nband);

    hamiltonian = new complex[Nband * Nband];
    temp1Matrix = new complex[Nband * Nband];
    temp2Matrix = new complex[Nband * Nband];
    newdMatrix = new complex[Nband * Nband];

    jMatrix = Create2D<complex>(Ndim, Nband*Nband);
    dipoleMatrix = Create2D<complex>(Ndim, Nband*Nband);
    uMatrix = new complex[Nband * Nband];
    ctransuMatrix = new complex[Nband * Nband];
    dephasingMatrix = new double[Nband * Nband];

    pMatrix = Create3D<complex>(this->kmesh->Ntotal, Ndim, Nband*Nband);
    sdMatrix = Create3D<complex>(this->kmesh->Ntotal, Ndim, Nband*Nband);
    edispersion = Create2D<double>(this->kmesh->Ntotal, Nband);

    tempIsuppz = new int[Nband];
    tempEval = new double[Nband];
    tempHmat = new complex[Nband*Nband];

    thermalE = 300.;
    thermalE /= 3.15775024804e5;     // hartree energy (4.3597447222071×10−18) / Boltzman constant(1.380649×10−23)

    eps = 1.0e-18;

    // Temporary variables
    double *tempInitVal = new double[Nband];
    double tempInitSum = 0;

    fill(&dmatrix[0][0], &dmatrix[0][0] + kmesh->Ntotal*Nband*Nband, 0.);
    ofstream bands ("bandstructure.dat");
    for (int k = 0; k < kmesh->Ntotal; ++k)
    {
        // edispersion generation & temporary umatrix generation
        int num_of_eig;
        material->GenHamiltonian(hamiltonian, kmesh->kgrid[k]);
        int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, hamiltonian, Nband, 0., 0., 0., 0.,
                                eps, &num_of_eig, edispersion[k], uMatrix, Nband, tempIsuppz );

        bands << kmesh->kgrid[k][0] << "  " << kmesh->kgrid[k][1] << "  " << kmesh->kgrid[k][2];
        bands << "  " << edispersion[k][0] << "  " << edispersion[k][1];
	      bands << "  " << edispersion[k][2] << "  " << edispersion[k][3] <<"\n";


        if (info != 0)
        {
            std::cerr << "Problem in lapack, info = " << info <<  std::endl;
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++ n)
                {
                    std::cout << hamiltonian[m*Nband + n] << "     ";
                }
                std::cout << endl;
            }
            for (int m = 0; m < Nband; ++m)
            {
                std::cout << edispersion[k][m] << "      ";
            }
            std::cout << endl;
            exit(EXIT_FAILURE);
        }
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                ctransuMatrix[m*Nband + n] = conj( uMatrix[n*Nband + m] );
            }
        }

        // Generate intial value in Hamiltonian gauge
        switch (initType)
        {
        case InitialValueType::UniformValence:
            for (int m = 0; m < Nband; ++m)
            {
                if (m < material->Nval)
                    dmatrix[k][m*Nband + m] = 1.0;
            }
            break;

        case InitialValueType::FermiDirac:
            tempInitSum = 0;
            for (int m = 0; m < Nband; ++m)
            {
                tempInitVal[m] = 1.0 / ( exp( (edispersion[k][m] - Ef) / thermalE ) + 1.0 );
                tempInitSum += tempInitVal[m];
            }
            for (int m = 0; m < Nband; ++m)
            {
                dmatrix[k][m*Nband + m] = tempInitVal[m] / tempInitSum;
            }
            break;

        case InitialValueType::Custom:
            material->GenInitialValue(&dmatrix[k][0], kmesh->kgrid[k]);
            break;

        default:
            cerr<<"Undefined Initial type\n";
            exit(EXIT_FAILURE);
            break;
        }
/////////////////////////////////////////////////////////////////////////////////here spin sdmatrix ////////
        if (gauge == GaugeType::LengthWannier)
        {
            MatrixMult(temp1Matrix, &dmatrix[k][0], ctransuMatrix, Nband);
            MatrixMult(&dmatrix[k][0], uMatrix, temp1Matrix, Nband);
        }

        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                initMatrix[k][m*Nband + n] = dmatrix[k][m*Nband + n];
            }
        }
/////////////////////////////////////////////////////////////////////////////////here spin sdmatrix ////////
        // pmatrix generation
        material->GenJMatrix(jMatrix, kmesh->kgrid[k]);
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            MatrixMult(temp1Matrix, jMatrix[iaxis], uMatrix, Nband);
            MatrixMult(pMatrix[k][iaxis], ctransuMatrix, temp1Matrix, Nband);
        }
    }
    bands.close();

    delete[] tempInitVal;

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
                dephasingMatrix[i*Nband + j] = dephasing_factor_intra;
                //dephasingMatrix[i*Nband + j] = 0.;
            else
                dephasingMatrix[i*Nband + j] = dephasing_factor_inter;
        }
    }
}

//checked

void SBEs::InitializeGeneral(const libconfig::Setting *_calc)
{
    // default values of parameters
    isDipoleZero = true;
    fill(ksfactor.begin(), ksfactor.end(), 1.0);
    RKorder = 5;
    dephasing_time_inter = -1.;
    dephasing_time_intra = -1.;

    const libconfig::Setting &calc = (*_calc);
    fill(Nk.begin(), Nk.end(), 1);
    // cout<<Nk[1]<<"\n";
    // Get calculation parameters
    try
    {
        const libconfig::Setting &npoints_setting = calc["Npoints"];
        for (int i = 0; i < npoints_setting.getLength(); ++i)
        {
            Nk[i] = npoints_setting[i];
        }
        //dt = calc.lookup("dt");
        dt = 0.1;
        calc.lookupValue("dt", dt);

        // read parameters which have default values - no error when they are omitted
        calc.lookupValue("T1", dephasing_time_intra );
        calc.lookupValue("T2", dephasing_time_inter );
        calc.lookupValue("shotNumber", shotNumber);
        calc.lookupValue("diagnostic", diagnostic);
        calc.lookupValue("RKorder", RKorder);

        calc.lookupValue("w1rx", rx1);
        calc.lookupValue("w1ry", ry1);
        calc.lookupValue("w1rz", rz1);

        calc.lookupValue("w2rx", rx2);
        calc.lookupValue("w2ry", ry2);
        calc.lookupValue("w2rz", rz2);

        calc.lookupValue("w3rx", rx3);
        calc.lookupValue("w3ry", ry3);
        calc.lookupValue("w3rz", rz3);

        calc.lookupValue("w4rx", rx4);
        calc.lookupValue("w4ry", ry4);
        calc.lookupValue("w4rz", rz4);

        calc.lookupValue("rmax", rmax);

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
            cerr << "Not-implemented RK order: " << RKorder << "\n";
            exit(EXIT_FAILURE);
            break;
    }
}
//checked
void SBEs::InitializeLaser(const libconfig::Setting *_laser)
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

SBEs::~SBEs()
{
    Delete2D<complex>(dmatrix, this->kmesh->Ntotal, Nband*Nband);
    Delete2D<complex>(initMatrix, this->kmesh->Ntotal, Nband*Nband);

    delete[] hamiltonian;
    delete[] newdMatrix;
    delete[] temp1Matrix, temp2Matrix;
    delete[] dephasingMatrix;
    delete[] uMatrix, ctransuMatrix;
    Delete2D<complex>(rkMatrix, NumRK, Nband*Nband);
    Delete2D<complex>(jMatrix, Ndim, Nband*Nband);
    Delete2D<complex>(dipoleMatrix, Ndim, Nband*Nband);



    Delete3D<complex>(sdMatrix, this->kmesh->Ntotal, Ndim, Nband*Nband);
    Delete3D<complex>(pMatrix, this->kmesh->Ntotal, Ndim, Nband*Nband);
    Delete2D<double>(edispersion, this->kmesh->Ntotal, Nband);

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

    delete[] tempIsuppz, tempEval, tempHmat;
}
////////////////////////////////////////////////////////// time propagation ////////
void SBEs::RunSBEs(int _kindex, double  _time)
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
////////////////////////////////////////////////////////// time propagation ////////
// array<double, Ndim> SBEs::GenCurrent(int _kindex, double _time)
tupleSpinCharge SBEs::GenCurrent(int _kindex, double _time)
{
    array<double, Ndim> outCurrent;
    array<complex, 4> spinDMatrix;
    // array<double, Ndim, Ndim> spinCurrent;
    std::array<std::array<double, Ndim>, Ndim> spinCurrent;
    std::array<std::array<double, Ndim>, 4> outCurrent_bands;
    // string test='test';
    fill(outCurrent.begin(), outCurrent.end(), 0.);
    //fill(outCurrent_bands.begin(), outCurrent_bands.end(), 0.);

    for (int ntmp = 0; ntmp<Ndim;++ntmp)
    {
      for (int mtmp = 0; mtmp<Ndim;++mtmp)
      {
        spinCurrent[ntmp][mtmp] = 0.;
      }
       for (int mtmp = 0; mtmp<4;++mtmp)
      {
        outCurrent_bands[ntmp][mtmp] = 0.;
      }

    }
    if (gauge == GaugeType::LengthWannier)
    {
        auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
        material->GenJMatrix(jMatrix, _tkp);

        // J = Tr{densitymatrix j}
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            for (int i = 0; i < Nband; ++i)
            {
                for (int ii = 0; ii < Nband; ++ii)
                {
                outCurrent_bands[iaxis][i] += real( dmatrix[_kindex][Nband*i + ii] * jMatrix[iaxis][Nband*ii + i] );
                }
                for (int j = 0; j < Nband; ++j)
                {
                    outCurrent[iaxis] += real( dmatrix[_kindex][Nband*i + j] * jMatrix[iaxis][Nband*j + i] );
                    //   //   //   //   //   //   //   //  spinCurrent
                    for (int saxis =0; saxis<Ndim;++saxis )
                    {
                      // int All = 4;
                      std::fill(spinDMatrix.begin(), spinDMatrix.end() , 0.);
                      for (int i = 0; i < 2; ++i)
                    	{
                    		for (int j = 0; j < 2; ++j)
                    		{
                    			for (int l = 0; l < 2; ++l)
                    			{
                    				spinDMatrix[i*2 + j] += dmatrix[_kindex][i*2 + l] * jMatrix[iaxis][l*2 + j];
                    			}
                    		}
                    	}

                      // MatrixMult(spinDMatrix,dmatrix[_kindex],sigma[saxis],2,2,2)
                      spinCurrent[iaxis][saxis] += real( dmatrix[_kindex][Nband*i + j] * jMatrix[iaxis][Nband*j + i] );
                    }
                      //   //   //   //   //   //   //   //   //   //   //   //   //
                     // cout << real(dmatrix[_kindex][Nband*i + j])<<' '<<real(jMatrix[iaxis][Nband*j + i]) <<'\n';
                }
            }
        }
    }
    else
    {
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            for (int i = 0; i < Nband; ++i)
            {
                for (int j = 0; j < Nband; ++j)
                {
                    outCurrent[iaxis] += real( dmatrix[_kindex][Nband*i + j] * pMatrix[_kindex][iaxis][Nband*j + i] );
                    // cout << dmatrix[_kindex][Nband*i + j]  <<'\n';
                }
            }
        }
    }
    // cout << outCurrent[0]<<'\n';
    // cout << 'check this' << material->Nval<<'\n';
    if (gauge == GaugeType::VelocityHamiltonian)
    {
        auto avector = fpulses->avlaser(_time);
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            outCurrent[iaxis] += material->Nval * avector[iaxis];
        }
    }

    return tupleSpinCharge (outCurrent,outCurrent_bands,spinCurrent);
}

void SBEs::GenDifferentialDM(complex *out, complex *input, int _kindex, double _time)
{
    array<double, Ndim> evector = fpulses->elaser(_time);
    auto avector = fpulses->avlaser(_time);
    auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);

    if (gauge == GaugeType::LengthWannier)
    {
        material->GenHamiltonian(hamiltonian, _tkp);
        MatrixMult(out, hamiltonian, input, Nband);
    }
    else if (gauge == GaugeType::VelocityHamiltonian)
    {
        fill(temp1Matrix, temp1Matrix+Nband*Nband, 0.);
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int iaxis = 0; iaxis < Ndim; ++iaxis)
                {
                    temp1Matrix[m*Nband + n] += avector[iaxis] * pMatrix[_kindex][iaxis][m*Nband + n];
                }
            }
        }
        for (int m = 0; m < Nband; ++m)
        {
            temp1Matrix[m*Nband + m] += edispersion[_kindex][m];
        }
        MatrixMult(out, temp1Matrix, input, Nband);
    }
    else if (gauge == GaugeType::LengthHamiltonian)
    {
        cout << "Not implemeted yet";
    }


    /*if (!isDipoleZero)
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
    }*/
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
    if (gauge == GaugeType::LengthWannier)
    {
        GenUMatrix(uMatrix, _tkp);
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
    else
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                out[m*Nband + n] -= (input[m*Nband + n] - initMatrix[_kindex][m*Nband + n]) * dephasingMatrix[m*Nband + n];
            }
        }
    }
}

void SBEs::GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint)
{
    material->GenHamiltonian(tempHmat, _kpoint);
    int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, tempHmat, Nband, 0., 0., 0., 0.,
                                eps, &num_of_eig, tempEval, _ustore, Nband, tempIsuppz );
    // cout << "inside GenUMatrix\n";
    //////////// extra
    // ofstream myfile;
    // myfile.open ("./test.dat");
    // myfile << "Writing this to a file.\n";
    // myfile.close();
  // for (int m = 0; m < Nband; ++m)
  //   {
  //     for (int n = 0; n < Nband; ++ n)
  //     {
  //     // cout << tempHmat[m*Nband + n] << "     ";
  //     }
  //     // cout << endl;
  //     }
    //////////// extra
    if (info != 0)
    {
        std::cerr << "Problem in lapack, info = " << info <<  std::endl;
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++ n)
            {
                std::cout << tempHmat[m*Nband + n] << "     ";
            }
            std::cout << endl;
        }
        for (int m = 0; m < Nband; ++m)
        {
            std::cout << tempEval[m] << "      ";
        }
        std::cout << endl;
        exit(EXIT_FAILURE);
    }
}

array<double, Ndim> SBEs::GenKpulsA(array<double, Ndim> _kpoint, double time)
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

void SBEs::WannierToHamiltonian(complex *out, complex *input, int _kindex, double _time)
{
    auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    GenUMatrix(uMatrix, _tkp);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            ctransuMatrix[m*Nband + n] = conj( uMatrix[n*Nband + m] );
        }
    }

    MatrixMult(temp1Matrix, input, uMatrix, Nband);
    MatrixMult(out, ctransuMatrix, temp1Matrix, Nband);
}

void SBEs::PrintInfo()
{
    cout << "===============================================\n";
    cout << "Target: " <<targetMaterial << endl;
    cout << "===============================================\n";

    material->PrintMaterialInformation();
    fpulses->Print_LaserInfo();

    cout << "===============================================\n";
    cout << "Current gauge: ";
    switch (gauge)
    {
    case GaugeType::LengthWannier:
        cout << "LengthWannier\n";
        break;
    case GaugeType::LengthHamiltonian:
        cout << "LengthHamiltonian\n";
        break;
    case GaugeType::VelocityHamiltonian:
        cout << "VelocityHamiltonian\n";
        break;

    default:
        cerr << "Undefined gauge??\n";
        exit(EXIT_FAILURE);
        break;
    }

    cout << "Inital value type: ";
    switch (initType)
    {
    case InitialValueType::UniformValence:
        cout << "UniformValence, Nval = " << material->Nval << std::endl;
        break;
    case InitialValueType::FermiDirac:
        cout << "FermiDirac, Ef = " << Ef << std::endl;
        break;
    case InitialValueType::Custom:
        cout << "Custom, make sure that you have implemented GenInitialValue " << std::endl;
        break;

    default:
        break;
    }
}
