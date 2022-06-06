// Class with perform SBE 
#include <iostream>

#ifndef sbesVG_h
#define sbesVG_h

#include <string>
#include <iostream>
#include <fstream>

#include "momaxis.h"
#include "constant.h"
#include "utility.h"

using namespace std;

class sbesVG
{
public:
    // edispersion - (kindex, bandindex)
    // pmatrix - (axis, kindex, bandindex*bandindex)
    complex ***pmatrix;
    double **edispersion;

    complex **dmatrix;

    // allocated memory with no value - to reduce costs for making temporoary variables
    complex *hamiltonian;
    complex *tempMatrix;
    complex **rkMatrix; 

    laser *fpulses;

    int Ngrid[Ndim];
    double Kaxis[Ndim][Ndim];
    double StartK[Ndim];
    int Nband, Nvb;

    // variables for Runge-Kutta method
    int RKorder, NumRK; // order of accumulated error in RK, and number of stages in RK algorithm
    double *tRK;
    double *bRK;
    double **aRK;

    momaxis *kgrid;

    double dephasing_time_inter, dephasing_time_intra;          // dephasing time, intrabnad dephasing has to be changed later
    double dephasing_factor_inter, dephasing_factor_intra;      // dephasing factor(1/dephasing time), to avoid expensive float division operation and remove if statement

    bool isAllocated;   // check whether memory is allocated at destructor

    sbesVG();
    ~sbesVG();

    // Allocating and initializing array, inside ReadInputParam
    void InitializeArray();

    // Reading files
    // return 1 for success, 0 for failure
    int ReadInputParam(string paramfilepath);
    int ReadPMatrix(string pmatfilepath, int _axis);
    int ReadEDispersion(string engfilepath);

    void RunSBEs(int _kindex, int _ktime);

    // functions return current. also real is x' and imag is y' convention
    // For later 3d case, we have to make conversion between laerser (x', y') and solid structyre (x, y, z)
    complex GenInterCurrent(int _kindex);
    complex GenIntraCurrent(int _kindex);

    void GenDifferentialWF(complex *out, int _kindex, double _time);
    void GenDifferentialDM(complex *out, complex *input, int _kindex, double _time);

    // I need to change below routines to more appropriate interface or name
    // All matrix = Nband * Nband
    // M - matrix, s - scalar, p - plus
    void M_MMmMM(complex *out, complex *_A, complex *_X, complex *_B, complex *_Y);
    void M_sMpsM(complex *out, complex _x, complex *_A, complex _y, complex *_B);

};

sbesVG::sbesVG()
{
    isAllocated = false;
    
    RKorder = 4;

    Kaxis[0][0] = 1.;   Kaxis[0][1] = 0.;   Kaxis[0][2] = 0.;
    Kaxis[1][0] = 0.;   Kaxis[1][1] = 1.;   Kaxis[1][2] = 0.;
    Kaxis[2][0] = 0.;   Kaxis[2][1] = 0.;   Kaxis[2][2] = 1.;

    StartK[0] = 0.;   StartK[0] = 0.;   StartK[0] = 0.;

    // for input error checking
    Nband = 0;  Nvb = 0;
    for (int i = 0; i < Ndim; ++i)  Ngrid[i] = 1;
}

sbesVG::~sbesVG()
{
    if (isAllocated)
    {
        delete kgrid;
        Delete3D<complex>(pmatrix, Ndim, this->kgrid->Ntotal, Nband*Nband);
        Delete2D<double>(edispersion, this->kgrid->Ntotal, Nband);
        Delete2D<complex>(dmatrix, this->kgrid->Ntotal, Nband*Nband);

        delete[] hamiltonian;
        delete[] tempMatrix;
        Delete2D<complex>(rkMatrix, NumRK, Nband*Nband);

        delete[] tRK;
        delete[] bRK;
        for (int i = 0; i < NumRK-1; ++i)
        {
            delete[] aRK[i];
        }
        delete[] aRK;
    }
}

void sbesVG::InitializeArray()
{
    kgrid = new momaxis(Ngrid, &Kaxis[0][0], StartK);
    kgrid->integral_method("Trapz");
    pmatrix = Create3D<complex>(Ndim, this->kgrid->Ntotal, Nband*Nband);
    edispersion = Create2D<double>(this->kgrid->Ntotal, Nband);
    dmatrix = Create2D<complex>(this->kgrid->Ntotal, Nband*Nband);

    hamiltonian = new complex[Nband * Nband];
    tempMatrix = new complex[Nband * Nband];

    for (int k = 0; k < this->kgrid->Ntotal; ++k)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                if (m == n && m < Nvb)
                    dmatrix[k][m*Nband + n] = 1.0;
                else
                    dmatrix[k][m*Nband + n] = 0.0;
            }
        }
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

    isAllocated = true;
}

void sbesVG::RunSBEs(int _kindex, int  _ktime)
{
    double time = fpulses->atmin + fpulses->dt * _ktime;
    GenDifferentialDM( rkMatrix[0], dmatrix[_kindex], _kindex, time);
    for (int irk = 0; irk < NumRK-1; ++irk)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                tempMatrix[m*Nband + n] = dmatrix[_kindex][m*Nband + n];
                for (int jrk = 0; jrk < irk+1; ++jrk)
                {
                    tempMatrix[m*Nband + n] += fpulses->dt*aRK[irk][jrk] * rkMatrix[jrk][m*Nband + n];
                }
            }
        }
        GenDifferentialDM(rkMatrix[irk+1], tempMatrix, _kindex, time + tRK[irk]*fpulses->dt);
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
    /*double tempDiffOcc = 0;
    for (int m = 1; m < Nband; ++m)
    {
        tempDiffOcc += real(dmatrix[_kindex][m*Nband + m]);
    }
    dmatrix[_kindex][0] = 1.0 - tempDiffOcc;*/
}

// check: 2D -> 3D
complex sbesVG::GenIntraCurrent(int _kindex)
{
    double jx = 0;
    double jy = 0;
    for (int m = 0; m < Nband; ++m)
    {
        jx -= real(pmatrix[0][_kindex][m*Nband + m] * dmatrix[_kindex][m*Nband + m]);
    }
    for (int m = 0; m < Nband; ++m)
    {
        jy -= real(pmatrix[1][_kindex][m*Nband + m] * dmatrix[_kindex][m*Nband + m]);
    }
    return complex(jx, jy);
}
complex sbesVG::GenInterCurrent(int _kindex)
{
    double jx = 0;
    double jy = 0;
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            //if (m > n) continue;
            if (m == n) continue;
            //jx += real(pmatrix[0][_kindex][n*Nband + m] * dmatrix[_kindex][m*Nband + n]);
            jx -= real(pmatrix[0][_kindex][m*Nband + n] * dmatrix[_kindex][m*Nband + n]);
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            //if (m > n) continue;
            if (m == n) continue;
            //jy += real(pmatrix[1][_kindex][n*Nband + m] * dmatrix[_kindex][m*Nband + n]);
            jy -= real(pmatrix[1][_kindex][m*Nband + n] * dmatrix[_kindex][m*Nband + n]);
        }
    }
    return complex(jx, jy);
}

// Functions for reading data from files

int sbesVG::ReadInputParam(string paramfilepath)
{
    ifstream inputParam(paramfilepath.data(), ios::in);
    if (inputParam.fail())
    {
        cout << "Error in ReadInputParam: Can't open file " << paramfilepath << "\n";
        inputParam.close();
        return 0; 
    }
    string inputline, paramName;
    bool isNgridSet = false;
    bool isStartKSet = false;
    int npulses;
    int itemp, idim = 0;
    double offset_after = 0, offset_before = 0;
    double E0, wfreq, ncycles, ellip, cep, t0, theta0;
    double dt;
    string env_name;
    while(getline(inputParam, inputline))
    {
        if (inputline.empty())  continue;
        trim_left(inputline);
        if (inputline.at(0)=='#')   continue;
        istringstream lineFstream(inputline);
        lineFstream >> paramName;

        if (paramName == "npulses")
        {
            lineFstream >> npulses;
            this->fpulses = new laser(npulses);
            for (int i = 0; i < npulses; ++i)
            {
                getline(inputParam, inputline);
                if (inputline.empty())
                {
                    --i;
                    continue;
                }
                trim_left(inputline);
                if (inputline.at(0)=='#')
                {
                    --i;
                    continue;
                }
                istringstream linePstream(inputline);
                linePstream >> E0 >> wfreq >> ncycles >> ellip >> cep >> t0 >> theta0 >> env_name;
                fpulses->PulseParam[i].Initialize(E0, ellip, wfreq, ncycles, cep, t0, theta0, env_name);
            }
        }
        else if (paramName == "Npoints")
        {
            while(lineFstream >> itemp)
            {
                Ngrid[idim] = itemp;
                idim++;
            }
            isNgridSet = true;
        }
        else if (paramName == "Nband") lineFstream >> Nband;     // Nx, Ny points
        else if (paramName == "Nvb") lineFstream >> Nvb;     // Nx, Ny points
        else if (paramName == "dt") lineFstream >> dt;                          // Time step
        else if (paramName == "T1") lineFstream >> dephasing_time_intra;                          // Dephasing time T1
        else if (paramName == "T2") lineFstream >> dephasing_time_inter;                          // Dephasing time T2
        //else if (paramName == "diagnostic") lineFstream >> diagnostic;          // diagonostic
        //else if (paramName == "ksfactor") lineFstream >> ksfactor;              // extend BZ with ksfactor ratio
        //else if (paramName == "shotNumber") lineFstream >> shotNumber;          // Number of snapshot, snapshot is disabled with shotNumber <=0
        else if (paramName == "offset_before") lineFstream >> offset_before;
        else if (paramName == "offset_after") lineFstream >> offset_after;
        else if (paramName == "RKorder") lineFstream >> RKorder;
        else if (paramName == "Vectors")
        {
            if (!isNgridSet)
            {
                cout << "Error in ReadInputParam: Try to read Vectors before reading N\n";
                return 0;
            }
            for (int i = 0; i < idim; ++i)
            {
                getline(inputParam, inputline);
                if (inputline.empty())
                {
                    --i;
                    continue;
                }
                trim_left(inputline);
                if (inputline.at(0)=='#')
                {
                    --i;
                    continue;
                }
                istringstream linePstream(inputline);
                for (int j = 0; j < idim; ++j)
                {
                    linePstream >> Kaxis[i][j];
                }
            }
        }
        else if (paramName == "StartK")
        {
            if (!isNgridSet)
            {
                cout << "Error in ReadInputParam: Try to read StartK before reading N\n";
                return 0;
            }
            for (int i = 0; i < idim; ++i)
            {
                lineFstream >> StartK[i];
            }
            isStartKSet = true;
        }
        else
        {
            cout << "Warning in ReadInputParam: undefined parameter name: " << paramName << "\n";
            //return 0;
        }
    }
    inputParam.close();
    if (!isStartKSet)
    {
        for (int i = 0; i < idim; ++i) for (int j = 0; j < idim; ++j)   StartK[i] -= Kaxis[j][i]/2;
    }

    if (!isNgridSet || Nband == 0 || Nvb == 0)
    {
        cout << "Error in ReadKGrid: Number paramters are not defined or improperly defined\n";
        cout << "Ngrid : " << Ngrid[0] << ", " << Ngrid[1] << ", " << Ngrid[2] << ", Nband : " << Nband << ", Nvb: " << Nvb << "\n";
        return 0;
    }
    // some post processing
    // dephainsg factor
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
    fpulses->laser_pulses(dt, offset_before, offset_after);

    InitializeArray();

    return 1;
}

int sbesVG::ReadEDispersion(string engfilepath)
{
    ifstream engfile(engfilepath.data(), ios::in);
    if (engfile.fail())
    {
        cout << "Problems in ReadEDispersion\n";
        engfile.close();
        exit(1); 
    }
    for (int i = 0; i < this->kgrid->N[0]; ++i)
    {
        for (int j = 0; j < this->kgrid->N[1]; ++j)
        {
            for (int k = 0; k < this->kgrid->N[2]; ++k)
            {                
                for (int m = 0; m < Nband; ++m)
                {
                    engfile >> edispersion[this->kgrid->index(&i, &j, &k)][m];
                }
            }
        }
    }
    engfile.close();
}
int sbesVG::ReadPMatrix(string pmatfilepath, int _axis)
{
    ifstream pmatfile(pmatfilepath.data(), ios::in|ios::binary);
    int kindex = 0;
    if (pmatfile.fail())
    {
        cout << "Problems in ReadPMatrix\n";
        pmatfile.close();
        exit(1);
    }
    for (int i = 0; i < this->kgrid->N[0]; ++i)
    {
        for (int j = 0; j < this->kgrid->N[1]; ++j)
        {
            for (int k = 0; k < this->kgrid->N[2]; ++k)
            {
                kindex = this->kgrid->index(&i, &j, &k);
                for (int m = 0; m < Nband; ++m)
                {
                    for (int n = 0; n < Nband; ++n)
                    {
                        if (m > n) continue;
                        pmatfile.read(reinterpret_cast<char*>(&pmatrix[_axis][kindex][m*Nband + n]), sizeof(complex));
                        pmatrix[_axis][kindex][n*Nband + m] = conj(pmatrix[_axis][kindex][m*Nband + n]);
                    }
                }
            }
        }
    }
    pmatfile.close();
}

// Temporary functions for matrix operation - switch to BLAS (ex: cblas_dgemm) later
// C = AX + BY

void sbesVG::GenDifferentialDM(complex *out, complex *input, int _kindex, double _time)
{
   /* GenDifferentialWF(hamiltonian, _kindex, _time);
    M_MMmMM(out, hamiltonian, input, input, hamiltonian );
    
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            //if (m > n) continue;
            if (m == n) out[m*Nband + n] = real(out[m*Nband + n]);
            //if (m == n)
            //    out[m*Nband + n] -= dephasing_factor_intra * input[m*Nband + n];
            //else
            else    out[m*Nband + n] -= dephasing_factor_inter * input[m*Nband + n];
        }
    }*/

    complex avector = fpulses->avlaser(&_time);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            if (m == n)
            {
                out[m*Nband + n] = 0;
                for (int l = 0; l < Nband; ++ l)
                {
                    if (m == l) continue;
                    out[m*Nband + n] += 2.* imag((real(avector) * pmatrix[0][_kindex][m*Nband + l] + imag(avector) * pmatrix[1][_kindex][m*Nband + l]) * input[m*Nband + l]);
                }
            }
            else
            {
                if (m > n) continue;
                out[m*Nband + n] = I * (edispersion[_kindex][m] - edispersion[_kindex][n]) * input[m*Nband + n];
                for (int l = 0; l < Nband; ++l)
                {
                    out[m*Nband + n] += I * ((real(avector) * pmatrix[0][_kindex][l*Nband + m] + imag(avector) * pmatrix[1][_kindex][l*Nband + m]) * input[l*Nband + n] 
                        - (real(avector) * pmatrix[0][_kindex][n*Nband + l] + imag(avector) * pmatrix[1][_kindex][n*Nband + l]) * input[m*Nband + l]);
                }
                out[m*Nband + n] -= dephasing_factor_inter * input[m*Nband + n];
            }
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            if (m > n)
                out[m*Nband + n] = conj(out[n*Nband + m]);
        }
    }
}
void sbesVG::GenDifferentialWF(complex *out, int _kindex, double _time)
{
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            out[m*Nband + n] = 0;
        }
    }
    complex avector = fpulses->avlaser(&_time);
    for (int m = 0; m < Nband; ++m)
    {
        out[m*Nband + m] -= I * edispersion[_kindex][m];
        for (int n = 0; n < Nband; ++n)
        {
            out[m*Nband + n] -= I * (real(avector) * pmatrix[0][_kindex][m*Nband + n] + imag(avector) * pmatrix[1][_kindex][m*Nband + n]); 
        }
    }
}
void sbesVG::M_MMmMM(complex *out, complex *_A, complex *_X, complex *_B, complex *_Y)
{
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            if (m > n) continue;
            out[m*Nband + n] = 0;
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            if (m > n) continue;
            for (int l = 0; l < Nband; ++l)
            {
                out[m*Nband+n] += _A[m*Nband+l]*_X[l*Nband+n] - _B[m*Nband+l]*_Y[l*Nband+n];
            }
        }
    }
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            if (m > n)
            {
                out[m*Nband + n] = conj(out[n*Nband + m]);
            }
        }
    }
}

void sbesVG::M_sMpsM(complex *out, complex _x, complex *_A, complex _y, complex *_B)
{
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            out[m*Nband+n] = _x*_A[m*Nband+n] + _y*_B[m*Nband+n];
        }
    }
}


#endif