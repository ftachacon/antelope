/// General tight binding class
/**
 * This class will handle every material except some special analytical models.
 * @author Dasol Kim
*/
#pragma once

#include "../wannier_system.h"
#include "../utility.h"

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

    int Nband;                              ///< Number of Wannier function per unit cell = Number of bands
    int Nrpts;                              ///< Number of Bravais lattice vector used in system
    double *weight;                         ///< Weight of unit cell 
    int **indexRvec;                        ///< index of lattice vectors  (R = n1 a1 + n2 a2 + n3 a3 ), indexRvec[irpt, 0] = n1
    int indexRorigin;                       ///< index which R = 0

    bool isDipoleZero;                      ///< if true choose e^ik(R+t) basis and make Wannier dipole zero

    // Note that interal index order of Wannier90 is ham_w[j, i, irpt] = <0j|H|Ri>
    // Remeber that c++ memory order is first-slow, last-fast and the fortran order is first-fast, last-slow  (row vs column order)
    complex ***ham_w;          ///< Hamiltonian components in Wannier basis, ham_w[irpt][i][j] = <0i|H|Rj> (Wannier90 order)
    complex ****pos_w;         ///< Position operator in Wannier basis, pos_w[irpt][i][j][0] = <0i|x|Rj> (1 - y, 2 - z)
    double ****rvec;                        ///< value of r vector ( r = R + tn - tm ), rvec[irpt][m][n][0] = (R+tn-tm).x
    double **Rvec;                          ///< value of R vector, Rvec[irpt][0] = Rx

    double Volume;                          ///< volume of unit cell

    int Nval;                               ///< Number of valence band (changing to fermi energy expression later?)

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

    void PrintMaterialInformation() override;
};

TightBinding::TightBinding( const libconfig::Setting *params )
{
    try
    {
        string w90_file = params->lookup("dataName");
        Nval = params->lookup("Nval");
        isDipoleZero = true;
        params->lookupValue("isDipoleZero", isDipoleZero);
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
                    ham_w[irpt][jtemp][itemp] = ftemp1 * I*ftemp2;
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
                        pos_w[irpt][jtemp][itemp][i] = ftemp1 * I*ftemp2;
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
    // degeneracy between valence and conduction bands are not considered yet
    for (int m = 0; m < Nband; ++m)
    {
        if (m < Nval)
            _dmstore[m*Nband + m] = 1.;
    }
    MatrixMult(tempMatrix, _dmstore, tempctransUmat, Nband);
    MatrixMult(_dmstore, tempUmat, tempMatrix, Nband);
}

void TightBinding::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
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
                    _hstore[m*Nband + n] += ham_w[irpt][m][n] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), rvec[irpt][m][n], 0.) );
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
                    _hstore[m*Nband + n] += ham_w[irpt][m][n] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), Rvec[irpt], 0.) );
                }
            }
        }
    }
}
void TightBinding::GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint)
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
                        _dstore[i][m*Nband + n] += pos_w[irpt][m][n][i] * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), rvec[irpt][m][n], 0.) );
                    }
                }
            }
        }
    }
}

void TightBinding::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
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
                            * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), rvec[irpt][m][n], 0.) );
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
                            * exp( I* std::inner_product(_kpoint.begin(), _kpoint.end(), Rvec[irpt], 0.) );
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
                                Nband*Nband, tempHamiltonian, Nband, 0., 0., 0., 0., 
                                eps, &num_of_eig, tempEigval, _ustore, Nband, isuppz );
    if (info != 0)
    {
        std::cerr << "Problem in lapack, info = " << info << std::endl;
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
