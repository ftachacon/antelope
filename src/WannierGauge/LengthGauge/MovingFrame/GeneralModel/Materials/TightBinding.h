/// General tight binding class
/**
 * This class will handle every material except some special analytical models.
 * @author Dasol Kim
*/
#pragma once

#include "../wannier_system.h"
#include "../utility.h"
#include "../momaxis.h"

#include <string>
#include <iostream>
#include <fstream>
#include <numeric>

// To prevnet include <complex.h> --> which undef complex
// As you see re-define keyword is very dangerous
#undef complex
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#define complex complex<double>

class TightBinding : public WannierMaterial
{
public:
    /// Bravais lattice vector and Reciprocal lattice vectors
    std::array<std::array<double, Ndim>, Ndim> vec_lattice, vec_reciprocal;

    int Nrpts;                              ///< Number of Bravais lattice vector used in system
    double *weight;                         ///< Weight of unit cell 
    int **indexRvec;                        ///< index of lattice vectors  (R = n1 a1 + n2 a2 + n3 a3 ), indexRvec[irpt][0] = n1
    int indexRorigin;                       ///< index which R = 0

    bool isDipoleZero;                      ///< if true choose e^ik(R+t) basis and make Wannier dipole zero

    // Note that interal index order of Wannier90 is ham_w[j, i, irpt] = <0j|H|Ri>
    // Remeber that c++ memory order is first-slow, last-fast and the fortran order is first-fast, last-slow  (row vs column order)
    complex ***ham_w;          ///< Hamiltonian components in Wannier basis, ham_w[irpt][i][j] = <0i|H|Rj> (Wannier90 order)
    complex ****pos_w;         ///< Position operator in Wannier basis, pos_w[irpt][i][j][0] = <0i|x|Rj> (1 - y, 2 - z)
    double ****rvec;                        ///< value of r vector ( r = R + tn - tm ), rvec[irpt][m][n][0] = (R+tn-tm).x
    double **Rvec;                          ///< value of R vector, Rvec[irpt][0] = Rx

    double Volume;                          ///< volume of unit cell

    momaxis *kmesh;                         ///< momentum axis used in calculation
    bool isMeshAllocated;                   ///< check pre-calculated matrices are allocated or not

    // Pre-calculated matrices for each k points
    complex **pre_ham;
    complex ***pre_pos;
    complex ***pre_jmat;

    // calculated differential of dipole and Hamiltonian
    complex ***pre_d1ham;           ///< 1st order talyor for Hamiltonian, pre_d1ham[kindex][m*Nband+n][i] = \sum_R Hmn(R)*R_i
    complex ****pre_d1pos;          ///< 1st order talyor for dipole, pre_d1pos[kindex][i][m*Nband+n][j] = \sum_R <0m|r_i|Rn>*R_j
    complex ****pre_d1jmat;         ///< 1st order talyor for derivative of Hamiltonian, pre_d1jmat[kindex][i][m*Nband+n][j] = -im *sum_R R_i*R_j*Hmn(R)
    //complex ****pre_d2ham;
    //complex ****pre_d2pos;

    int Nval;                               ///< Number of valence band 
    double FermiE;                          ///< Fermi energy in atomic unit
    bool isFermiEUsed;                      ///< if true Fermi energy is used to generate initial value. if false valence band number is used.
    double thermalE;                        ///< temperature in atomic unit (300K is defalut)

    double eps;                             ///< error epsilon used in this class

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    // temporary matrix used in internal storage
    // size = Nband * Nband
    complex *tempHamiltonian, *tempUmat, *tempctransUmat;
    complex *tempMatrix;

    double *tempEigval;     // size = Nband
    int *isuppz;                         // size = Nband*2
    

    TightBinding( const libconfig::Setting *params );
    ~TightBinding();
    void Allocate();
    void SetReciprocalBasis();
    void GenInitialValue(complex *_dmstore, std::array<double, Ndim> _kpoint) override;
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override;
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    void GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void GenHamiltonianPrimitive(complex *_hstore, std::array<double, Ndim> _kpoint);
    void GenDipolePrimitive(complex **_dstore, std::array<double, Ndim> _kpoint);
    void GenJMatrixPrimitive(complex **_jstore, std::array<double, Ndim> _kpoint);

    void PrintMaterialInformation() override;

    void CalculateKMesh(momaxis *_kmesh);
    std::tuple<int, std::array<double, Ndim> > GetNearestK(std::array<double, Ndim> _kpoint);
};

TightBinding::TightBinding( const libconfig::Setting *params )
{
    isMeshAllocated = false;
    try
    {
        string w90_file = params->lookup("dataName");
        // If FermiE exists, valence band number is ignored
        isFermiEUsed =  params->lookupValue("FermiE", FermiE);
        if (isFermiEUsed)
            FermiE = FermiE / au_eV;
        else
            Nval = params->lookup("Nval");
        thermalE = 300.;
        params->lookupValue("temperature", thermalE);
        thermalE /= 3.15775024804e5;     // hartree energy (4.3597447222071×10−18) / Boltzman constant(1.380649×10−23)

        isDipoleZero = false;
        //isDipoleZero = true;
        //params->lookupValue("isDipoleZero", isDipoleZero);
        ifstream w90_data(w90_file.c_str());

        int itemp, jtemp;
        double ftemp1, ftemp2;
        if (w90_data.is_open())
        {
            string line;
            getline(w90_data, line);    // First line is header
            // lattice_vector
            for (int i = 0; i < Ndim; ++i)
            {
                w90_data >> vec_lattice[i][0] >> vec_lattice[i][1] >> vec_lattice[i][2];
            }
            for (int i = 0; i < Ndim; ++i) for (int j = 0; j < Ndim; ++j) vec_lattice[i][j] /= au_angstrom;
            w90_data >> Nband >> Nrpts; // num_wann, nrpts in Wannier90
            Allocate(); // allocate
            for (int irpt = 0; irpt < Nrpts; ++irpt)
            {
                w90_data >> weight[irpt];
                weight[irpt] = 1.0/weight[irpt];
            }
            // Read Hamiltonian components
            for (int irpt = 0; irpt < Nrpts; ++irpt)
            {
                // Read R vector
                w90_data >> indexRvec[irpt][0] >> indexRvec[irpt][1] >> indexRvec[irpt][2];
                // read Hamitonian components corresponding to <0m|H|Rn>
                for (int m = 0; m < Nband*Nband; ++m)
                {
                    w90_data >> jtemp >> itemp >> ftemp1 >> ftemp2;
                    ham_w[irpt][jtemp-1][itemp-1] = (ftemp1 + I*ftemp2)/au_eV;
                }
            }
            // Read position matrix elements
            for (int irpt = 0; irpt < Nrpts; ++irpt)
            {
                // Read R vector
                // Below is repeated behavior - maybe I can add error-checking routine here
                w90_data >> indexRvec[irpt][0] >> indexRvec[irpt][1] >> indexRvec[irpt][2];
                if (indexRvec[irpt][0] == 0 && indexRvec[irpt][1] == 0 && indexRvec[irpt][2] == 0)
                    indexRorigin = irpt;
                for (int m = 0; m < Nband*Nband; ++m)
                {
                    w90_data >> jtemp >> itemp; 
                    for (int i = 0; i < Ndim; ++i)
                    {
                        w90_data >> ftemp1 >> ftemp2;
                        pos_w[irpt][jtemp-1][itemp-1][i] = (ftemp1 * I*ftemp2)/au_angstrom;
                    }
                }
            }
        }

        // construct r vectors
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int irpt = 0; irpt < Nrpts; ++irpt)
                {
                    for (int iaxis = 0; iaxis < Ndim; ++iaxis)
                        rvec[irpt][m][n][iaxis] = indexRvec[irpt][0]*vec_lattice[0][iaxis] 
                                                + indexRvec[irpt][1]*vec_lattice[1][iaxis] 
                                                + indexRvec[irpt][2]*vec_lattice[2][iaxis]
                                                + real(pos_w[indexRorigin][n][n][iaxis] - pos_w[indexRorigin][m][m][iaxis]);
                }
            }
        }

        // construct R vectors
        for (int irpt = 0; irpt < Nrpts; ++irpt)
        {
            for (int iaxis = 0; iaxis < Ndim; ++iaxis)
            {
                Rvec[irpt][iaxis] = indexRvec[irpt][0]*vec_lattice[0][iaxis] 
                                    + indexRvec[irpt][1]*vec_lattice[1][iaxis] 
                                    + indexRvec[irpt][2]*vec_lattice[2][iaxis];
            }
        }

        // 'S' for safe minimum (1/sfmin != nan)
        // 'E' for machine epsilon
        eps = LAPACKE_dlamch('S');
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    SetReciprocalBasis();
}

void TightBinding::Allocate()
{
    weight = new double[Nrpts];
    indexRvec = Create2D<int>(Nrpts, Ndim);
    rvec = Create4D<double>(Nrpts, Nband, Nband, Ndim);
    Rvec = Create2D<double>(Nrpts, Ndim);
    ham_w = Create3D<complex>(Nrpts, Nband, Nband);
    pos_w = Create4D<complex>(Nrpts, Nband, Nband, Ndim);

    tempHamiltonian = new complex[Nband*Nband];
    tempUmat = new complex[Nband*Nband];
    tempctransUmat = new complex[Nband*Nband];
    tempMatrix = new complex[Nband*Nband];
    tempEigval = new double[Nband];
    isuppz = new int[2*Nband];
}
TightBinding::~TightBinding()
{
    delete[] weight;
    Delete2D(indexRvec, Nrpts, Ndim);
    Delete4D(rvec, Nrpts, Nband, Nband, Ndim);
    Delete2D(Rvec, Nrpts, Ndim);
    Delete3D(ham_w, Nrpts, Nband, Nband);
    Delete4D(pos_w, Nrpts, Nband, Nband, Ndim);

    delete[] tempHamiltonian, tempUmat, tempctransUmat, tempMatrix, tempEigval, isuppz;

    if (isMeshAllocated)
    {
        Delete2D<complex>(pre_ham, kmesh->Ntotal, Nband*Nband);
        Delete3D<complex>(pre_pos, kmesh->Ntotal, Ndim, Nband*Nband);
        Delete3D<complex>(pre_jmat, kmesh->Ntotal, Ndim, Nband*Nband);
        Delete3D<complex>(pre_d1ham, kmesh->Ntotal, Nband*Nband, Ndim);
        Delete4D<complex>(pre_d1pos, kmesh->Ntotal, Ndim, Nband*Nband, Ndim);
        Delete4D<complex>(pre_d1jmat, kmesh->Ntotal, Ndim, Nband*Nband, Ndim);
    }
}

void TightBinding::GenInitialValue(complex *_dmstore, std::array<double, Ndim> _kpoint)
{
    fill(_dmstore, _dmstore + Nband*Nband, 0.);
    GenUMatrix(tempUmat, _kpoint);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            tempctransUmat[m*Nband + n] = conj( tempUmat[n*Nband + m] );
        }
    }

    // set valence = 1., conduction = 0. initial condition
    // Fermi-Dirac distribution is used when FermiE is applied
    // tempEigval is calculated inside the GenUMatrix, be careful about order or sideeffects.
    for (int m = 0; m < Nband; ++m)
    {
        if (isFermiEUsed)
        {
            _dmstore[m*Nband + m] = 1.0 / ( exp( (tempEigval[m] - FermiE) / thermalE ) + 1.0 );
        }
        else
        {
            if (m < Nval)
                _dmstore[m*Nband + m] = 1.;
        }
    }
    MatrixMult(tempMatrix, _dmstore, tempctransUmat, Nband);
    MatrixMult(_dmstore, tempUmat, tempMatrix, Nband);
}

void TightBinding::GenHamiltonianPrimitive(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    fill(_hstore, _hstore + Nband*Nband, 0.);
    if (isDipoleZero)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int irpt = 0; irpt < Nrpts; ++irpt)
                {
                    _hstore[m*Nband + n] += weight[irpt]*ham_w[irpt][m][n] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), rvec[irpt][m][n], 0.) );
                }
            }
        }
    }
    else
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int irpt = 0; irpt < Nrpts; ++irpt)
                {
                    _hstore[m*Nband + n] += weight[irpt]*ham_w[irpt][m][n] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), Rvec[irpt], 0.) );
                }
            }
        }
    }
}
void TightBinding::GenDipolePrimitive(complex **_dstore, std::array<double, Ndim> _kpoint)
{
    for (int i = 0; i < Ndim; ++i)
        fill(_dstore[i], _dstore[i] + Nband*Nband, 0.);
    if (!isDipoleZero)
    {
        for (int i = 0; i < Ndim; ++i)
        {
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    for (int irpt = 0; irpt < Nrpts; ++irpt)
                    {
                        _dstore[i][m*Nband + n] += weight[irpt]*pos_w[irpt][m][n][i] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), Rvec[irpt], 0.) );
                    }
                }
            }
        }
    }
}

void TightBinding::GenJMatrixPrimitive(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double temp;
    for (int i = 0; i < Ndim; ++i)
        fill(_jstore[i], _jstore[i] + Nband*Nband, 0.);
    if (isDipoleZero)
    {
        for (int i = 0; i < Ndim; ++i)
        {
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    for (int irpt = 0; irpt < Nrpts; ++irpt)
                    {
                        _jstore[i][m*Nband + n] -= I*rvec[irpt][m][n][i]* ham_w[irpt][m][n] 
                            * weight[irpt]*exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), rvec[irpt][m][n], 0.) );
                    }
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < Ndim; ++i)
        {
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    for (int irpt = 0; irpt < Nrpts; ++irpt)
                    {
                        _jstore[i][m*Nband + n] -= I*Rvec[irpt][i]* ham_w[irpt][m][n] 
                            * weight[irpt]*exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), Rvec[irpt], 0.) );
                    }
                }
            }
        }
    }
}

void TightBinding::GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint)
{
    GenHamiltonian(tempHamiltonian, _kpoint);
    //lapack_int LAPACKE_zheevr( int matrix_layout, char jobz, char range, char uplo, 
    //    lapack_int n, lapack_complex_double* a, lapack_int lda, double vl, double vu, lapack_int il, lapack_int iu, 
    //    double abstol, lapack_int* m, double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* isuppz );
    int num_of_eig; 
    int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, tempHamiltonian, Nband, 0., 0., 0., 0., 
                                eps, &num_of_eig, tempEigval, _ustore, Nband, isuppz );
    if (info != 0)
    {
        std::cerr << "Problem in lapack, info = " << info << std::endl;
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++ n)
            {
                std::cout << tempHamiltonian[m*Nband + n] << "     ";
            }
            std::cout << endl;
        }
        for (int m = 0; m < Nband; ++m)
        {
            std::cout << tempEigval[m] << "      ";
        }
        std::cout << endl;
        exit(1);
    }
}

void TightBinding::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    auto [kindex, _dkp] = GetNearestK(_kpoint);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            _hstore[m*Nband + n] = pre_ham[kindex][m*Nband + n];
        }
    }

    // 1st - order
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            for (int i = 0; i < Ndim; ++i)
            {
                _hstore[m*Nband + n] += pre_d1ham[kindex][m*Nband + n][i] * _dkp[i];
            }
        }
    }
}

void TightBinding::GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint)
{
    auto [kindex, _dkp] = GetNearestK(_kpoint);
    for (int i = 0; i < Ndim; ++i)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                _dstore[i][m*Nband + n] = pre_pos[kindex][i][m*Nband + n];
            }
        }
    }

    // 1st - order
    for (int i = 0; i < Ndim; ++i)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int j = 0; j < Ndim; ++j)
                {
                    _dstore[i][m*Nband + n] += pre_d1pos[kindex][i][m*Nband + n][j] * _dkp[j];
                }
            }
        }
    }
}

void TightBinding::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    auto [kindex, _dkp] = GetNearestK(_kpoint);
    for (int i = 0; i < Ndim; ++i)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                _jstore[i][m*Nband + n] = pre_jmat[kindex][i][m*Nband + n];
            }
        }
    }

    // 1st - order
    for (int i = 0; i < Ndim; ++i)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int j = 0; j < Ndim; ++j)
                {
                    _jstore[i][m*Nband + n] += pre_d1jmat[kindex][i][m*Nband + n][j] * _dkp[j];
                }
            }
        }
    }
}
std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > TightBinding::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void TightBinding::SetReciprocalBasis()
{
    Volume = 0.;
    Volume += vec_lattice[0][0] * (vec_lattice[1][1]*vec_lattice[2][2] - vec_lattice[1][2]*vec_lattice[2][1])
            + vec_lattice[0][1] * (vec_lattice[1][2]*vec_lattice[2][0] - vec_lattice[1][0]*vec_lattice[2][2])
            + vec_lattice[0][2] * (vec_lattice[1][0]*vec_lattice[2][1] - vec_lattice[1][1]*vec_lattice[2][0]);
    vec_reciprocal[0] = CrossProduct(vec_lattice[1], vec_lattice[2]);
    vec_reciprocal[1] = CrossProduct(vec_lattice[2], vec_lattice[0]);
    vec_reciprocal[2] = CrossProduct(vec_lattice[0], vec_lattice[1]);
    for (int i = 0; i < Ndim; ++i) for (int j = 0; j < Ndim; ++j) vec_reciprocal[i][j] /= Volume;

    BZorigin[0] = 0.;   BZorigin[1] = 0.;   BZorigin[2] = 0.;
    for (int i = 0; i < Ndim; ++i) for (int j = 0; j < Ndim; ++j) BZaxis[i*Ndim + j] = vec_reciprocal[i][j];
}

void TightBinding::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << "Wannier90 data \n";
    cout << "============================================\n";
}

void TightBinding::CalculateKMesh(momaxis *_kmesh)
{
    kmesh = _kmesh;

    // pre-calculated Hamiltonian and dipole 
    pre_ham = Create2D<complex>(kmesh->Ntotal, Nband*Nband);
    pre_pos = Create3D<complex>(kmesh->Ntotal, Ndim, Nband*Nband);
    pre_jmat = Create3D<complex>(kmesh->Ntotal, Ndim, Nband*Nband);

    for (int kindex = 0; kindex < kmesh->Ntotal; ++kindex)
    {
        GenHamiltonianPrimitive(pre_ham[kindex], kmesh->kgrid[kindex]);
        GenDipolePrimitive(pre_pos[kindex], kmesh->kgrid[kindex]);
        GenJMatrixPrimitive(pre_jmat[kindex], kmesh->kgrid[kindex]);
    }
    
    // taylor of Hamiltonian, dipole and current - first order
    pre_d1ham = Create3D<complex>(kmesh->Ntotal, Nband*Nband, Ndim);
    pre_d1pos = Create4D<complex>(kmesh->Ntotal, Ndim, Nband*Nband, Ndim);
    pre_d1jmat = Create4D<complex>(kmesh->Ntotal, Ndim, Nband*Nband, Ndim);

    for (int kindex = 0; kindex < kmesh->Ntotal; ++kindex)
    {
        for (int m = 0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                for (int irpt = 0; irpt < Nrpts; ++irpt)
                {
                    for (int iaxis = 0; iaxis < Ndim; ++iaxis)
                    {
                        //pre_d1ham[kindex][m*Nband + n][iaxis] += ham_w[irpt][m][n] * Rvec[irpt][iaxis];
                        pre_d1ham[kindex][m*Nband + n][iaxis] += weight[irpt]* ham_w[irpt][m][n]  * (indexRvec[irpt][iaxis]*2*pi);
                        for (int jaxis = 0; jaxis < Ndim; ++jaxis)
                        {
                            //pre_d1pos[kindex][iaxis][m*Nband + n][jaxis] += pos_w[irpt][m][n][iaxis] * Rvec[irpt][jaxis];
                            pre_d1pos[kindex][iaxis][m*Nband + n][jaxis] += weight[irpt]* pos_w[irpt][m][n][iaxis] * (indexRvec[irpt][jaxis]*2*pi);
                            //pre_d1jmat[kindex][iaxis][m*Nband + n][jaxis] -= I* Rvec[irpt][iaxis] * Rvec[irpt][jaxis] * ham_w[irpt][m][n];
                            pre_d1jmat[kindex][iaxis][m*Nband + n][jaxis] -= weight[irpt]* I* Rvec[irpt][iaxis] * (indexRvec[irpt][jaxis]*2*pi) * ham_w[irpt][m][n];
                        }
                    }
                }
            }
        }
    }
    isMeshAllocated = true; 
}

std::tuple<int, std::array<double, Ndim> > TightBinding::GetNearestK(std::array<double, Ndim> _kpoint)
{
    array<double, Ndim> bcoeff;
    fill(bcoeff.begin(), bcoeff.end(), 0.0);
    for (int i = 0; i < Ndim; ++i)
    {
        for (int j =0; j < Ndim; ++j)
        {
            bcoeff[i] += vec_lattice[i][j] * _kpoint[j];
        }
        bcoeff[i] /= (2*pi);
    }
    for (int i = 0; i < Ndim; ++i)
    {
        // range 0<= < 1
        if (bcoeff[i] < 0 || bcoeff[i] >= 1)
            bcoeff[i] -= floor(bcoeff[i]);
        bcoeff[i] *= kmesh->N[i];
    }
    array<int, Ndim> indexArr = { int(round(bcoeff[0])),  int(round(bcoeff[1])), int(round(bcoeff[2]))};
    for (int i = 0; i < Ndim; ++i)
    {
        bcoeff[i] -= indexArr[i];
        bcoeff[i] /= kmesh->N[i];
    }
    for (int i = 0; i < Ndim; ++i) indexArr[i] = indexArr[i]%kmesh->N[i];
    return std::make_tuple(kmesh->index(indexArr),  bcoeff);
}