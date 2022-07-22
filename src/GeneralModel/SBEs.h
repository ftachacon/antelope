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

#include "mpi.h"

#include <libconfig.h++>

// To prevnet include <complex.h> --> which undef complex
// As you see re-define keyword is very dangerous (not exactly keword...)
#undef complex
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#define complex std::complex<double>

#include "momaxis.h"
#include "laser.h"
#include "constant.h"
#include "utility.h"

#include "system_list.h"

class SBEs
{
public:
    const GaugeType gauge;
    InitialValueType initType;

    complex *dmatrix;
    complex **initMatrix; 

    BaseMaterial *material;  ///< dipole&Hamiltonian generator

    complex **rkmatrix; ///< buffers for RK methods
    MPI_Win *win_dm; ///< window for dmatrix
    MPI_Win *win_rk; ///< window for rkmatrix

    // allocated memory with no value - to reduce costs for making temporoary variables
    complex *hamiltonian;
    complex **dipoleMatrix, **jMatrix;
    complex *newdMatrix;
    complex *temp1Matrix, *temp2Matrix;

    complex *uMatrix, *ctransuMatrix;
    double *dephasingMatrix;

    complex ***pMatrix;
    double **edispersion;

    laser *fpulses;

    int Nband;
    std::array<int, Ndim> Nk;

    int kstart, kend;   ///< start and end index of k-space

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
    std::string targetMaterial;
    std::array<double, Ndim> ksfactor;
    int shotNumber;
    int diagnostic;
private:
    double dt;

    // used for lapack routines
    int num_of_eig;
    complex *tempLapackHmat; // separate vraible to prevent unexpected overwriting of hamiltonian matrix
    double *tempEval;
    int *tempIsuppz;
    double eps;


public:
    double dephasing_time_inter, dephasing_time_intra;          // dephasing time, intrabnad dephasing has to be changed later
    double dephasing_factor_inter, dephasing_factor_intra;      // dephasing factor(1/dephasing time), to avoid expensive float division operation and remove if statement

    bool isDipoleZero;                                  ///< default is true. Wannier dipole is non-zero in very rare cases.
    bool isInterIntra;                                  ///< default is true. Separate inter/intra when this is true. Total current is in inter when false.

    SBEs(const libconfig::Setting * _cfg, GaugeType _gauge, complex *_dmatrix, complex **_rkmatrix,
        MPI_Win *_win_dm, MPI_Win *_win_rk, int _kstart, int _kend);
    ~SBEs();
    
    void RunSBEs(double _time);

    std::array<double, Ndim> GenCurrent(int _kindex, double _time);
    std::tuple<std::array<double, Ndim>, std::array<double, Ndim> > GenInterIntraCurrent(int _kindex, double _time);

    void GenDifferentialDM(complex *out, complex *input, int _kindex, double _time);

    void InitializeGeneral(const libconfig::Setting *_calc);     ///< Initialize general stuff inside constructor
    void InitializeLaser(const libconfig::Setting *_laser);     ///< Initialize laser inside constructor

    void GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint);    ///< Generate U matrix

    std::array<double, Ndim> GenKpulsA(std::array<double, Ndim> _kpoint, double time);

    void WannierToHamiltonian(complex *out, complex *input, int _kindex, double _time); ///< Convert Wannier representation to Hamiltonian representation

    void PrintInfo();   ///< Print calculation information

};

SBEs::SBEs(const libconfig::Setting * _cfg, GaugeType _gauge, complex *_dmatrix, complex **_rkmatrix,
    MPI_Win *_win_dm, MPI_Win *_win_rk, int _kstart, int _kend) : gauge(_gauge)
{
    kstart = _kstart;   kend = _kend;
    dmatrix = _dmatrix; rkmatrix = _rkmatrix; // do not allocate, just get shared memory window pointer
    win_dm = _win_dm; win_rk = _win_rk;

    const libconfig::Setting &cfg = (*_cfg);
    InitializeGeneral( &cfg["calc"]);
    InitializeLaser( &cfg["laser"]);
    // Check target, get Nband
    try
    {
        //targetMaterial = cfg.lookup("target").c_str();
        if (!cfg.lookupValue("target", targetMaterial))
        {
            std::cerr << "No 'target' in input file\n";
            std::exit(EXIT_FAILURE);
        }
        
        //Nband = cfg[targetMaterial.c_str()].lookup("Nband");
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        std::cerr << "No 'target' in inputParam.cfg or no 'Nband' parameter" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    //cout << "Target :  " << targetMaterial << std::endl;

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
    material = InitializeMaterial(targetMaterial, _cfg);
    
    Nband = material->Nband;

    auto [bzaxes, bzori] = material->GenBrillouinzone();

    // apply ksfactor
    bzaxes[0] *= ksfactor[0];   bzaxes[1] *= ksfactor[0];   bzaxes[2] *= ksfactor[0];
    bzaxes[3] *= ksfactor[1];   bzaxes[4] *= ksfactor[1];   bzaxes[5] *= ksfactor[1];
    bzaxes[6] *= ksfactor[2];   bzaxes[7] *= ksfactor[2];   bzaxes[8] *= ksfactor[2];

    if (cfg["calc"].exists("BZaxes"))
    {
        const libconfig::Setting &bzaxes_setting = cfg["calc"]["BZaxes"];
        for (int i = 0; i < bzaxes_setting.getLength(); ++i)
        {
            bzaxes[i] = bzaxes_setting[i];
        }
    }
    if (cfg["calc"].exists("BZorigin"))
    {
        const libconfig::Setting &bzorigin_setting = cfg["calc"]["BZorigin"];
        for (int i = 0; i < bzorigin_setting.getLength(); ++i)
        {
            bzori[i] = bzorigin_setting[i];
        }
    }

    // generate k-grid
    kmesh = new momaxis( Nk, bzaxes, bzori );

    /*if (isWannier90)
    {
        dynamic_cast<Wannier90*>(material)->CalculateKMesh(kmesh);
    }*/
    
    initMatrix = Create2D<complex>(kend-kstart, Nband*Nband);

    hamiltonian = new complex[Nband * Nband];
    temp1Matrix = new complex[Nband * Nband];
    temp2Matrix = new complex[Nband * Nband];
    newdMatrix = new complex[Nband * Nband];

    jMatrix = Create2D<complex>(Ndim, Nband*Nband);
    dipoleMatrix = Create2D<complex>(Ndim, Nband*Nband);
    uMatrix = new complex[Nband * Nband];
    ctransuMatrix = new complex[Nband * Nband];
    dephasingMatrix = new double[Nband * Nband];

    pMatrix = Create3D<complex>(kend-kstart, Ndim, Nband*Nband);
    edispersion = Create2D<double>(kend-kstart, Nband);

    tempIsuppz = new int[2*Nband];
    tempEval = new double[Nband];
    tempLapackHmat = new complex[Nband*Nband];

    thermalE = 300.;
    thermalE /= 3.15775024804e5;     // hartree energy (4.3597447222071×10−18) / Boltzman constant(1.380649×10−23)

    eps = 1.0e-18;

    // Temporary variables 
    double *tempInitVal = new double[Nband];
    double tempInitSum = 0;

    std::fill(&dmatrix[kstart*Nband*Nband], &dmatrix[kstart*Nband*Nband] + (kend-kstart)*Nband*Nband, 0.);
    
    for (int k = 0; k < kend-kstart; ++k)
    {
        int kglobal = k + kstart; // global index = kglobal, local index = k
        // edispersion generation & temporary umatrix generation
        int num_of_eig;
        material->GenHamiltonian(hamiltonian, kmesh->kgrid[kglobal]);
        int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, hamiltonian, Nband, 0., 0., 0, 0, 
                                eps, &num_of_eig, edispersion[k], uMatrix, Nband, tempIsuppz );

        if (info != 0)
        {
            std::cerr << "Problem in lapack, info = " << info <<  std::endl;
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++ n)
                {
                    std::cout << hamiltonian[m*Nband + n] << "     ";
                }
                std::cout << std::endl;
            }
            for (int m = 0; m < Nband; ++m)
            {
                std::cout << edispersion[k][m] << "      ";
            }
            std::cout << std::endl;
            std::exit(EXIT_FAILURE);
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
                    dmatrix[kglobal*Nband*Nband + m*Nband + m] = 1.0;
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
                dmatrix[kglobal*Nband*Nband + m*Nband + m] = tempInitVal[m] / tempInitSum;
            }
            break;

        case InitialValueType::Custom:
            material->GenInitialValue(&dmatrix[kglobal*Nband*Nband], kmesh->kgrid[kglobal]);
            break;
        
        default:
            std::cerr<<"Undefined Initial type\n";
            std::exit(EXIT_FAILURE);
            break;
        }

        if (gauge == GaugeType::WannierMoving || gauge == GaugeType::WannierStatic)
        {
            MatrixMult(temp1Matrix, &dmatrix[kglobal*Nband*Nband], ctransuMatrix, Nband);
            MatrixMult(&dmatrix[kglobal*Nband*Nband], uMatrix, temp1Matrix, Nband);
        }

        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                initMatrix[k][m*Nband + n] = dmatrix[kglobal*Nband*Nband + m*Nband + n];
            }
        }

        // pmatrix generation
        material->GenJMatrix(jMatrix, kmesh->kgrid[kglobal]);
        
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            MatrixMult(temp1Matrix, jMatrix[iaxis], uMatrix, Nband);
            MatrixMult(pMatrix[k][iaxis], ctransuMatrix, temp1Matrix, Nband);
        }
    }

    delete[] tempInitVal;

    // allocation for rk realted variables
    tRK = new double[NumRK-1];
    bRK = new double[NumRK];
    aRK = new double*[NumRK-1];


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
            throw std::out_of_range("Undefined Runge-Kutta order");
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

void SBEs::InitializeGeneral(const libconfig::Setting *_calc)
{
    // default values of parameters
    isDipoleZero = true;
    isInterIntra = true;
    std::fill(ksfactor.begin(), ksfactor.end(), 1.0);
    RKorder = 5;
    dephasing_time_inter = -1.;
    dephasing_time_intra = -1.;
    shotNumber = 0;
    isInterIntra = false;

    const libconfig::Setting &calc = (*_calc);
    std::fill(Nk.begin(), Nk.end(), 1);
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
        //calc.lookupValue("RKorder", RKorder); // disable temporary because NumRK is hard-coded.
        calc.lookupValue("InterIntra", isInterIntra );
        // ksfactor - value & array
        const libconfig::Setting &ksfactor_setting = calc["ksfactor"];
        for (int i = 0; i < ksfactor_setting.getLength(); ++i)
        {
            ksfactor[i] = ksfactor_setting[i];
        }
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        std::cerr << "Setting not found while search calc group \n";
        std::exit(EXIT_FAILURE);
    }
    catch(const libconfig::SettingTypeException &tex )
    {
        std::cerr << "Type exception in calc group\n";
        std::exit(EXIT_FAILURE);
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
            std::cerr << "Not-implemented RK order: " << RKorder << "\n";
            std::exit(EXIT_FAILURE);
            break;
    }
}

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
        std::string t_env_name;
        for (int i = 0; i < pulseNum; ++i)
        {
            if ( pulses[i].lookupValue("E0", t_E0) && pulses[i].lookupValue("ellip", t_ellip) && pulses[i].lookupValue("w0", t_w0) 
                && pulses[i].lookupValue("ncycles", t_ncycles) && pulses[i].lookupValue("cep", t_cep) && pulses[i].lookupValue("t0", t_t0) 
                && pulses[i].lookupValue("phix", t_phix) && pulses[i].lookupValue("env_name", t_env_name) )
            {
                t_cep *= pi/180.0;
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
                std::cerr << "Error while initilizng laser pulses. Check config file\n";
                std::exit(EXIT_FAILURE);
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
        std::cerr << "No laser, check config file \n";
        std::exit(EXIT_FAILURE);
    }
    fpulses->Initialize(dt, offset_before, offset_after);
}

SBEs::~SBEs()
{
    Delete2D<complex>(initMatrix, kend-kstart, Nband*Nband);

    delete[] hamiltonian;
    delete[] newdMatrix;
    delete[] temp1Matrix, temp2Matrix;
    delete[] dephasingMatrix;
    delete[] uMatrix, ctransuMatrix;
    Delete2D<complex>(jMatrix, Ndim, Nband*Nband);
    Delete2D<complex>(dipoleMatrix, Ndim, Nband*Nband);

    Delete3D<complex>(pMatrix, kend-kstart, Ndim, Nband*Nband);
    Delete2D<double>(edispersion, kend-kstart, Nband);

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

    delete[] tempIsuppz, tempEval, tempLapackHmat;
}

void SBEs::RunSBEs(double  _time)
{
    // always fence before differential equation, for input array.
    MPI_Win_fence(0, *win_dm);
    // k0 = f(tn, yn)
    for (int k = kstart; k < kend; ++k)
        GenDifferentialDM( &rkmatrix[0][k*Nband*Nband], &dmatrix[k*Nband*Nband], k, _time);

    for (int irk = 1; irk < NumRK; ++irk)
    {
        // NumRK --> buffer
        std::copy( &dmatrix[kstart*Nband*Nband], 
            &dmatrix[kstart*Nband*Nband] + (kend-kstart)*Nband*Nband, 
            &rkmatrix[NumRK][kstart*Nband*Nband]);
        for (int k = kstart; k < kend; ++k)
        {
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    for (int jrk = 0; jrk < irk; ++jrk)
                    {
                        // yn + dt *\sum_{j=0}^{i-1} aRk[i-1][j] * k_j
                        rkmatrix[NumRK][k*Nband*Nband + m*Nband + n] 
                            += fpulses->dt*aRK[irk-1][jrk] * rkmatrix[jrk][k*Nband*Nband + m*Nband + n];
                    }
                }
            }
        }
        MPI_Win_fence(0, win_rk[NumRK]);
        // k_i = f(tn + dt*t[i-1], yn + dt*\sum_{j=0}^{i-1} aRk[i-1][j] * k_j)
        for (int k = kstart; k < kend; ++k)
        {
            GenDifferentialDM(&rkmatrix[irk][kstart*Nband*Nband], &rkmatrix[NumRK][kstart*Nband*Nband], 
                k, _time + tRK[irk-1]*fpulses->dt);
        }
    }
    for (int k = kstart; k < kend; ++k)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int irk = 0; irk < NumRK; ++irk)
                {
                    // y_{n+1} = y_n + dt*\sum_{i=0}^{N-1} b[i] * k_i
                    dmatrix[k*Nband*Nband + m*Nband + n] += fpulses->dt*bRK[irk] * rkmatrix[irk][k*Nband*Nband + m*Nband + n];
                }
            }
        }
    }
}

std::array<double, Ndim> SBEs::GenCurrent(int _kindex, double _time)
{
    std::array<double, Ndim> outCurrent;
    std::fill(outCurrent.begin(), outCurrent.end(), 0.);
    std::array<double, Ndim> _tkp;

    // K+A(t) for moving frame and k for static frame
    if (gauge == GaugeType::WannierMoving)
    {
        _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    }
    else
    {
        _tkp = kmesh->kgrid[_kindex];
    }
    material->GenJMatrix(jMatrix, _tkp);

    // Wannier ues jMatrix
    if (gauge == GaugeType::WannierMoving || gauge == GaugeType::WannierStatic)
    {
        // J = Tr{densitymatrix j}
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            for (int i = 0; i < Nband; ++i)
            {
                for (int j = 0; j < Nband; ++j)
                {
                    outCurrent[iaxis] += real( dmatrix[_kindex*Nband*Nband + Nband*i + j] * jMatrix[iaxis][Nband*j + i] );
                }
            }
        }
    }
    else // both Hamiltonian and Velocity uses momentum matrix element
    {
        GenUMatrix(uMatrix, _tkp);
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                ctransuMatrix[m*Nband + n] = conj(uMatrix[n*Nband + m]);
            }
        }
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            // jMatrix --> -pMatrix (generate momentum matrix element from jMaxtirx, -1 from electric charge)
            MatrixMult(temp1Matrix, jMatrix[iaxis], uMatrix, Nband);
            MatrixMult(temp2Matrix, ctransuMatrix, temp1Matrix, Nband);
            for (int i = 0; i < Nband; ++i)
            {
                for (int j = 0; j < Nband; ++j)
                {
                    outCurrent[iaxis] += real( dmatrix[_kindex*Nband*Nband + Nband*i + j] * temp2Matrix[Nband*j + i] );
                }
            }
        }
    }
    // for VG, j = -(p + A(t)), additional A(t) component added
    if (gauge == GaugeType::Velocity)
    {
        auto avector = fpulses->avlaser(_time);
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            outCurrent[iaxis] += material->Nval * avector[iaxis];
        }
    }
    
    return outCurrent;
}
std::tuple<std::array<double, Ndim>, std::array<double, Ndim> > SBEs::GenInterIntraCurrent(int _kindex, double _time)
{
    std::array<double, Ndim> interCurrent, intraCurrent;
    std::fill(interCurrent.begin(), interCurrent.end(), 0.);
    std::fill(intraCurrent.begin(), intraCurrent.end(), 0.);
    std::array<double, Ndim> _tkp;

    complex *dmatrix_at_k = &dmatrix[_kindex*Nband*Nband];

    // K+A(t) for moving frame and k for static frame
    if (gauge == GaugeType::WannierMoving)
    {
        _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    }
    else
    {
        _tkp = kmesh->kgrid[_kindex];
    }
    material->GenJMatrix(jMatrix, _tkp);

    // U and U^dagger
    GenUMatrix(uMatrix, _tkp);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            ctransuMatrix[m*Nband + n] = conj(uMatrix[n*Nband + m]);
        }
    }

    // Wannier-->Hamiltonian: U^dagger * rho * U = rho^{H}
    if (gauge == GaugeType::WannierMoving || gauge == GaugeType::WannierStatic)
    {
        MatrixMult(temp1Matrix, dmatrix_at_k, uMatrix, Nband);
        MatrixMult(newdMatrix, ctransuMatrix, temp1Matrix, Nband);
    }
    else
    {
        copy(dmatrix_at_k, dmatrix_at_k+Nband*Nband, newdMatrix);
    }
        
    for (int iaxis = 0; iaxis < Ndim; ++iaxis)
    {
        // jMatrix --> -pMatrix (generate momentum matrix element from jMaxtirx, -1 from electric charge)
        MatrixMult(temp1Matrix, jMatrix[iaxis], uMatrix, Nband);
        MatrixMult(temp2Matrix, ctransuMatrix, temp1Matrix, Nband);
        for (int i = 0; i < Nband; ++i)
        {
            for (int j = 0; j < Nband; ++j)
            {
                if (i == j)
                    intraCurrent[iaxis] += real( newdMatrix[Nband*i + j] * temp2Matrix[Nband*j + i] );
                else
                    interCurrent[iaxis] += real( newdMatrix[Nband*i + j] * temp2Matrix[Nband*j + i] );
            }
        }
    }
    // for VG, j = -(p + A(t)), additional A(t) component added
    if (gauge == GaugeType::Velocity)
    {
        auto avector = fpulses->avlaser(_time);
        for (int iaxis = 0; iaxis < Ndim; ++iaxis)
        {
            intraCurrent[iaxis] += material->Nval * avector[iaxis];
        }
    }
    
    return make_tuple(interCurrent, intraCurrent);
}

void SBEs::GenDifferentialDM(complex *out, complex *input, int _kindex, double _time)
{
    std::array<double, Ndim> evector = fpulses->elaser(_time);
    auto avector = fpulses->avlaser(_time);
    auto _tkp = GenKpulsA(kmesh->kgrid[_kindex], _time);
    int local_kindex = _kindex - kstart;

    // d\rho/dt = -i [H0 + E(t)D, \rho] + ...
    // hamiltonian = H0 + E(t)D
    if (gauge == GaugeType::WannierMoving)
    {
        material->GenHamiltonian(hamiltonian, _tkp);
        if (!material->isDipoleZero)
        {
            material->GenDipole(dipoleMatrix, _tkp);
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    hamiltonian[m*Nband + n] += (evector[0]*dipoleMatrix[0][m*Nband + n] 
                                                + evector[1]*dipoleMatrix[1][m*Nband + n]
                                                + evector[2]*dipoleMatrix[2][m*Nband + n]);
                }
            }
        }
        // out = (H0 + E(t)D)\rho
        MatrixMult(out, hamiltonian, input, Nband);
    }
    // d\rho/dt = -i [H0 + A(t)P, \rho] + ...
    // hamiltonian = H0 + A(t)P
    else if (gauge == GaugeType::Velocity)
    {
        std::fill(hamiltonian, hamiltonian+Nband*Nband, 0.);
        for (int m = 0; m < Nband; ++m)
        {
            hamiltonian[m*Nband + m] += edispersion[local_kindex][m];
            for (int n = 0; n < Nband; ++n)
            {
                for (int iaxis = 0; iaxis < Ndim; ++iaxis)
                {
                    hamiltonian[m*Nband + n] += avector[iaxis] * pMatrix[local_kindex][iaxis][m*Nband + n];
                }
            }
        }
        // out = (H0 + A(t)P)\rho
        MatrixMult(out, hamiltonian, input, Nband);
    }
    else if (gauge == GaugeType::WannierStatic)
    {
        // material->GenHamiltonian(hamiltonian, kmesh->kgrid[_kindex]);
        // if (!material->isDipoleZero)
        // {
        //     material->GenDipole(dipoleMatrix, kmesh->kgrid[_kindex]);
        //     for (int m = 0; m < Nband; ++m)
        //     {
        //         for (int n = 0; n < Nband; ++n)
        //         {
        //             hamiltonian[m*Nband + n] += (evector[0]*dipoleMatrix[0][m*Nband + n] 
        //                                         + evector[1]*dipoleMatrix[1][m*Nband + n]
        //                                         + evector[2]*dipoleMatrix[2][m*Nband + n]);
        //         }
        //     }
        // }
        // // out = (H0 + E(t)D)\rho
        // MatrixMult(out, hamiltonian, input, Nband);
    }
    else
    {
        std::cout << "Not implemeted yet";
        std::exit(EXIT_FAILURE);
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
            //-i[H, \rho]_mn = -i(H\rho - \rho*H)_mn = -i ((H\rho) - (H\rho)^\dagger)_mn
            out[m*Nband + n] -= conj(hamiltonian[n*Nband + m]);
            out[m*Nband + n] *= -I;
        }
    }
    // Add dephainsg process here!
    if (gauge == GaugeType::WannierMoving || gauge == GaugeType::WannierStatic)
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
                out[m*Nband + n] -= (input[m*Nband + n] - initMatrix[local_kindex][m*Nband + n]) * dephasingMatrix[m*Nband + n];
            }
        }
    }
}

void SBEs::GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint)
{
    material->GenHamiltonian(tempLapackHmat, _kpoint); 
    int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, tempLapackHmat, Nband, 0., 0., 0, 0, 
                                eps, &num_of_eig, tempEval, _ustore, Nband, tempIsuppz );
    if (info != 0)
    {
        std::cerr << "Problem in lapack, info = " << info <<  std::endl;
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++ n)
            {
                std::cout << tempLapackHmat[m*Nband + n] << "     ";
            }
            std::cout << std::endl;
        }
        for (int m = 0; m < Nband; ++m)
        {
            std::cout << tempEval[m] << "      ";
        }
        std::cout << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

std::array<double, Ndim> SBEs::GenKpulsA(std::array<double, Ndim> _kpoint, double time)
{
    auto avector = fpulses->avlaser(time);
    std::array<double, Ndim> tkp = {_kpoint[0]+avector[0], _kpoint[1]+avector[1], _kpoint[2]+avector[2]};
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
    std::cout << "===============================================\n";
    std::cout << "Target: " <<targetMaterial << std::endl;
    std::cout << "===============================================\n";

    material->PrintMaterialInformation();
    fpulses->Print_LaserInfo();

    std::cout << "===============================================\n";
    std::cout << "Current gauge: ";
    switch (gauge)
    {
    case GaugeType::WannierMoving:
        std::cout << "WannierMoving\n";
        break;
    case GaugeType::WannierStatic:
        std::cout << "WannierStatic\n";
        break;
    case GaugeType::Velocity:
        std::cout << "Velocity\n";
        break;
    case GaugeType::HamiltonianStatic:
        std::cout << "HamiltonianStatic\n";
        break;
    
    default:
        std::cerr << "Undefined gauge??\n";
        std::exit(EXIT_FAILURE);
        break;
    }

    std::cout << "Inital value type: ";
    switch (initType)
    {
    case InitialValueType::UniformValence:
        std::cout << "UniformValence, Nval = " << material->Nval << std::endl;
        break;
    case InitialValueType::FermiDirac:
        std::cout << "FermiDirac, Ef = " << Ef << std::endl;
        break;
    case InitialValueType::Custom:
        std::cout << "Custom, make sure that you have implemented GenInitialValue " << std::endl;
        break;

    default:
        break;
    }
}