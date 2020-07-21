/// Hamiltonian and dipole generator for Bi2Se3 surface states
/**
 * Equation for this class comes from PhysRevB.84.115413 and Paper of Denitsa (add refer after publihsed)
 * This class assumes system is pure 2d - no slant plane
 * @author Dasol Kim
 * @author Alexis Agustín  Chacón Salazar
*/
#pragma once

#include "../wannier_system.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

// To prevnet include <complex.h> --> which undef complex
// As you see re-define keyword is very dangerous
#undef complex
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#define complex complex<double>

class BieSe3surf : public WannierMaterial
{
public:
    // Haldnae parameters
    double t1, t2, M0, phi0, a0; 
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[4];    ///< coefficent of pauli matrix

    double hcomp[6];
    double acomp[4];    ///< A0, A11, A12, A14
    double bcomp[4];    ///< B0, B11, B12, B14

    double m11;

    double angle_a0[3], angle_b0[3];

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    BieSe3surf( const libconfig::Setting *params );
    void SetBasis();
    void GenInitialValue(complex *_dmstore, std::array<double, Ndim> _kpoint) override;
    void GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint) override;
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

BieSe3surf::BieSe3surf( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;
    /*if (params->lookupValue("t1", t1)
        && params->lookupValue("t2", t2)
        && params->lookupValue("Mt2", M0)
        && params->lookupValue("phi0", phi0)
        && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;
        M0 *= t2; 
    }
    else
    {
        cerr << "Some Haldane paramters are missing" << endl;
        exit(EXIT_FAILURE);
    }*/
    // hard-coded paramter
    acomp[0] = -0.0255 / au_eV;
    acomp[1] = 0.1937 / au_eV;
    acomp[2] = 0.2240 / au_eV;
    acomp[3] = 0.0551 / au_eV;

    bcomp[0] = 0.0164 / au_eV;
    bcomp[1] = 0.1203 / au_eV;
    bcomp[2] = 0.3264 / au_eV;
    bcomp[3] = 0.0 / au_eV;

    m11 = -1.6978 / au_eV;

    a0 = 4.14 / au_angstrom;
    SetBasis();

    // lattice vector: (a0, 0), (-a0/2, sqrt(3.)*a0/2)
    // --> reciprocal : (2*pi/a0, 2*pi/a0/sqrt(3.) ), (0., 4*pi/a0/sqrt(3.))
    // same with Haldane notation (sqrt(3.)*a0 --> a0)
    //double kxMax = pi/sqrt(3.)/a0;  double kyMax = 2.*pi/3./a0;
    double kxMax = pi/a0;  double kyMax = 2.*pi/sqrt(3.)/a0;
    BZaxis = { 4*kxMax,       0,         0,
                  0,       4*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void BieSe3surf::GenInitialValue(complex *_dmstore, std::array<double, Ndim> _kpoint)
{
    fill(_dmstore, _dmstore + Nband*Nband, 0.);
    GenBcomp(_kpoint);

    // double Bnorm = sqrt( Bcomp[1]*Bcomp[1] + Bcomp[2]*Bcomp[2] + Bcomp[3]*Bcomp[3] );
    // if (Bnorm > eps)
    // {
    //     _dmstore[0] = 0.5 - Bcomp[3]/(2.*Bnorm);        _dmstore[1] = -(Bcomp[1] - I*Bcomp[2])/(2.*Bnorm);
    //     _dmstore[2] = conj(_dmstore[1]);                _dmstore[3] = 1.0 - _dmstore[0];
    // }
    // else
    // {
    //     // B1, B2 --> 0, B3 < 0 case
    //     // _dmstore[0] = 0.5 - Bcomp[3]/(2.*Bnorm);        _dmstore[1] = -(Bcomp[1] - I*Bcomp[2])/(2.*Bnorm);
    //     // _dmstore[2] = conj(_dmstore[1]);                _dmstore[3] = 1.0 - _dmstore[0];
    //     _dmstore[0] = 1.0;        _dmstore[1] = 0.0;
    //     _dmstore[2] = 0.0;        _dmstore[3] = 1.0;
    // }
    complex *tempUmat = new complex[Nband*Nband];
    complex *tempctransUmat = new complex[Nband*Nband];
    complex *tempHmat = new complex[Nband*Nband];
    double *tempEval = new double[Nband];
    int *tempIsuppz = new int[2*Nband];
    GenHamiltonian(tempHmat, _kpoint);
    //lapack_int LAPACKE_zheevr( int matrix_layout, char jobz, char range, char uplo, 
    //    lapack_int n, lapack_complex_double* a, lapack_int lda, double vl, double vu, lapack_int il, lapack_int iu, 
    //    double abstol, lapack_int* m, double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* isuppz );
    int num_of_eig; 
    int info = LAPACKE_zheevr( LAPACK_ROW_MAJOR, 'V', 'A', 'U',
                                Nband, tempHmat, Nband, 0., 0., 0., 0., 
                                eps, &num_of_eig, tempEval, tempUmat, Nband, tempIsuppz );
    if (info != 0)
    {
        std::cerr << "Problem in lapack, info = " << info << std::endl;
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
        exit(1);
    }
    GenUMatrix(tempUmat, _kpoint);
    for (int m = 0; m < Nband; ++m)
    {
        for (int n = 0; n < Nband; ++n)
        {
            tempctransUmat[m*Nband + n] = conj( tempUmat[n*Nband + m] );
        }
    }

    
    double FermiE = 0.0 / au_eV;
    double thermalE = 300.;
    thermalE /= 3.15775024804e5;     // hartree energy (4.3597447222071×10−18) / Boltzman constant(1.380649×10−23)

    // set valence = 1., conduction = 0. initial condition
    // Fermi-Dirac distribution is used when FermiE is applied
    // tempEigval is calculated inside the GenUMatrix, be careful about order or sideeffects.
    for (int m = 0; m < Nband; ++m)
    {
        //if (isFermiEUsed)
        {
            _dmstore[m*Nband + m] = 1.0 / ( exp( (tempEval[m] - FermiE) / thermalE ) + 1.0 );
        }
        // else
        // {
        //     if (m < Nval)
        //         _dmstore[m*Nband + m] = 1.;
        // }
    }
    // tempHmat is used as temporary matrix for below
    MatrixMult(tempHmat, _dmstore, tempctransUmat, Nband);
    MatrixMult(_dmstore, tempUmat, tempHmat, Nband);
    
    delete[] tempUmat, tempctransUmat, tempHmat, tempEval, tempIsuppz;
    
}

void BieSe3surf::GenUMatrix(complex *_ustore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);
    double Bnorm = sqrt( Bcomp[1]*Bcomp[1] + Bcomp[2]*Bcomp[2] + Bcomp[3]*Bcomp[3] );
    double waveNorm = sqrt( 2.*Bnorm* (Bnorm + Bcomp[3]) );
    if (waveNorm > eps)
    {
        _ustore[0] = (Bcomp[3] + Bnorm)/waveNorm;            _ustore[1] = -(Bcomp[1] - I*Bcomp[2])/waveNorm;
        _ustore[2] = (Bcomp[1] + I*Bcomp[2])/waveNorm;       _ustore[3] = (Bcomp[3] + Bnorm)/waveNorm;
    }
    else
    {
        // B1, B2 --> 0, B3 < 0 case
        _ustore[0] = 0.;            _ustore[1] = 1.;
        _ustore[2] = 1.;            _ustore[3] = 0.;
    }
}

void BieSe3surf::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    _hstore[0] = Bcomp[0] + Bcomp[3];     _hstore[1] = Bcomp[1] - I*Bcomp[2];
    _hstore[2] = Bcomp[1] + I*Bcomp[2];   _hstore[3] = Bcomp[0] - Bcomp[3];
}

void BieSe3surf::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderHcomp[5];    double yderHcomp[5];
    double xderBcomp[4];    double yderBcomp[4];
    fill(xderBcomp, xderBcomp+4, 0.);
    fill(yderBcomp, yderBcomp+4, 0.);
    fill(xderHcomp, xderHcomp+5, 0.);
    fill(yderHcomp, yderHcomp+5, 0.);
    double cosAvecSum = 0.;
    double sinAvecSum = 0.;
    double cosBvecSum = 0.;
    double sinBvecSum = 0.;

    double derxcosAvecSum = 0.;
    double derxsinAvecSum = 0.;
    double derxcosBvecSum = 0.;
    double derxsinBvecSum = 0.;
    double derycosAvecSum = 0.;
    double derysinAvecSum = 0.;
    double derycosBvecSum = 0.;
    double derysinBvecSum = 0.;
    for (int i = 0; i < 3; ++i)
    {
        angle_a0[i] = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0[i] = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];

        cosAvecSum += cos(angle_a0[i]);
        sinAvecSum += sin(angle_a0[i]);
        cosBvecSum += cos(angle_b0[i]);
        sinBvecSum += sin(angle_b0[i]);

        derxcosAvecSum -= vec_a[i][0] * sin(angle_a0[i]);
        derxsinAvecSum += vec_a[i][0] * cos(angle_a0[i]);
        derxcosBvecSum -= vec_b[i][0] * sin(angle_b0[i]);
        derxsinBvecSum += vec_b[i][0] * cos(angle_b0[i]);

        derycosAvecSum -= vec_a[i][1] * sin(angle_a0[i]);
        derysinAvecSum += vec_a[i][1] * cos(angle_a0[i]);
        derycosBvecSum -= vec_b[i][1] * sin(angle_b0[i]);
        derysinBvecSum += vec_b[i][1] * cos(angle_b0[i]);
    }

    double w = -2*pi/3;

    xderHcomp[0] = 2*acomp[0] * derxcosAvecSum + 2*bcomp[0] * derxcosBvecSum;
    xderHcomp[1] = -2*acomp[3] * sin(w)*(vec_a[1][0]*cos(angle_a0[1]) - vec_a[2][0]*cos(angle_a0[2]))
                    + 2*bcomp[3] * (vec_b[0][0]*cos(angle_b0[0]) 
                    + cos(w)*(vec_b[1][0]*cos(angle_b0[1]) - vec_b[2][0]*cos(angle_b0[2])) );
    xderHcomp[2] = -2*bcomp[3] * sin(w)*(vec_b[1][0]*cos(angle_b0[1]) - vec_b[2][0]*cos(angle_b0[2]))
                    - 2*acomp[3] * (vec_a[0][0]*cos(angle_a0[0]) 
                    + cos(w)*(vec_a[1][0]*cos(angle_a0[1]) - vec_a[2][0]*cos(angle_a0[2])) );
    xderHcomp[3] = 2*acomp[2] * derxsinAvecSum;
    xderHcomp[4] = -2*bcomp[2] * derxsinBvecSum;
    xderHcomp[5] = 2*acomp[1] * derxcosAvecSum + 2*bcomp[1] * derxcosBvecSum;

    yderHcomp[0] = 2*acomp[0] * derycosAvecSum + 2*bcomp[0] * derycosBvecSum;
    yderHcomp[1] = -2*acomp[3] * sin(w)*(vec_a[1][1]*cos(angle_a0[1]) - vec_a[2][1]*cos(angle_a0[2]))
                    + 2*bcomp[3] * (vec_b[0][1]*cos(angle_b0[0]) 
                    + cos(w)*(vec_b[1][1]*cos(angle_b0[1]) - vec_b[2][1]*cos(angle_b0[2])) );
    yderHcomp[2] = -2*bcomp[3] * sin(w)*(vec_b[1][1]*cos(angle_b0[1]) - vec_b[2][1]*cos(angle_b0[2]))
                    - 2*acomp[3] * (vec_a[0][1]*cos(angle_a0[0]) 
                    - cos(w)*(vec_a[1][1]*cos(angle_a0[1]) - vec_a[2][1]*cos(angle_a0[2])) );
    yderHcomp[3] = 2*acomp[2] * derysinAvecSum;
    yderHcomp[4] = -2*bcomp[2] * derysinBvecSum;
    yderHcomp[5] = 2*acomp[1] * derycosAvecSum + 2*bcomp[1] * derycosBvecSum;

    // hcomp[0] = 2*acomp[0] * cosAvecSum + 2*bcomp[0] * cosBvecSum;
    // hcomp[1] = -2*acomp[3] * sin( w*(sin(angle_a0[1]) - sin(angle_a0[2])) ) 
    //             + 2*bcomp[3] * ( sin(angle_b0[0]) + cos( w*(sin(angle_b0[1]) - sin(angle_b0[2])) ) );
    // hcomp[2] = -2*bcomp[3] * sin( w*(sin(angle_b0[1]) - sin(angle_b0[2])) )
    //             -2*acomp[3] * (sin(angle_a0[0]) + cos( w*(sin(angle_a0[1]) + sin(angle_a0[2])) ) );
    // hcomp[3] = 2*acomp[2] * sinAvecSum;
    // hcomp[4] = -2*bcomp[2] * sinBvecSum;
    // hcomp[5] = 2*acomp[1] * cosAvecSum + 2*bcomp[1] * cosBvecSum + m11;

    double hcomp5gamma = 6*acomp[1] + 6*bcomp[1] + m11;
    double bcoeff = sqrt(1. - bcomp[0]*bcomp[0]/(bcomp[1]*bcomp[1]));
    // Bcomp[0] = hcomp[0] + bcomp[0]/bcomp[1] * (-hcomp[5] + hcomp5gamma);
    // Bcomp[1] = bcoeff * hcomp[1];
    // Bcomp[2] = bcoeff * hcomp[2];
    // Bcomp[3] = bcoeff * hcomp[3];

    xderBcomp[0] = xderHcomp[0] - bcomp[0]/bcomp[1] * xderHcomp[5];
    xderBcomp[1] = bcoeff * xderHcomp[1];
    xderBcomp[2] = bcoeff * xderHcomp[2];
    xderBcomp[3] = bcoeff * xderHcomp[3];

    yderBcomp[0] = yderHcomp[0] - bcomp[0]/bcomp[1] * yderHcomp[5];
    yderBcomp[1] = bcoeff * yderHcomp[1];
    yderBcomp[2] = bcoeff * yderHcomp[2];
    yderBcomp[3] = bcoeff * yderHcomp[3];

    _jstore[0][0] = -(xderBcomp[0] + xderBcomp[3]);      _jstore[1][0] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[0][1] = -(xderBcomp[1] - I*xderBcomp[2]);    _jstore[1][1] = -(yderBcomp[1] - I*yderBcomp[2]);
    _jstore[0][2] = -(xderBcomp[1] + I*xderBcomp[2]);    _jstore[1][2] = -(yderBcomp[1] + I*yderBcomp[2]);
    _jstore[0][3] = -(xderBcomp[0] - xderBcomp[3]);      _jstore[1][3] = -(yderBcomp[0] - yderBcomp[3]);
    

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > BieSe3surf::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void BieSe3surf::SetBasis()
{
    // z-components are ignored
    vec_a[0][0] = a0;                   vec_a[0][1] = 0.;
    vec_a[1][0] = -a0/2;                vec_a[1][1] = sqrt(3.)/2*a0;
    vec_a[2][0] = -a0/2;                vec_a[2][1] = -sqrt(3.)/2*a0;

    vec_b[0][0] = 0.;                   vec_b[0][1] = sqrt(3.)*a0/3;
    vec_b[1][0] = -a0/2;                vec_b[1][1] = -sqrt(3.)/6*a0;
    vec_b[2][0] = a0/2;                 vec_b[2][1] = -sqrt(3.)/6*a0;
}

void BieSe3surf::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << "Bi2Se3 surface - paramters are hard-coded\n";
    cout << "============================================\n";
}

void BieSe3surf::GenBcomp(std::array<double, Ndim> _kpoint)
{
    fill(Bcomp, Bcomp+4, 0.);
    double cosAvecSum = 0.;
    double sinAvecSum = 0.;
    double cosBvecSum = 0.;
    double sinBvecSum = 0.;
    for (int i = 0; i < 3; ++i)
    {
        angle_a0[i] = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0[i] = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];

        cosAvecSum += cos(angle_a0[i]);
        sinAvecSum += sin(angle_a0[i]);
        cosBvecSum += cos(angle_b0[i]);
        sinBvecSum += sin(angle_b0[i]);
    }

    double w = -2*pi/3;
    
    hcomp[0] = 2*acomp[0] * cosAvecSum + 2*bcomp[0] * cosBvecSum;
    hcomp[1] = -2*acomp[3] * sin(w)*(sin(angle_a0[1]) - sin(angle_a0[2]))  
                + 2*bcomp[3] * (sin(angle_b0[0]) + cos(w)*(sin(angle_b0[1]) - sin(angle_b0[2])) );
    hcomp[2] = -2*bcomp[3] * sin(w)*(sin(angle_b0[1]) - sin(angle_b0[2])) 
                - 2*acomp[3] * (sin(angle_a0[0]) + cos(w)*(sin(angle_a0[1]) + sin(angle_a0[2])) );
    hcomp[3] = 2*acomp[2] * sinAvecSum;
    hcomp[4] = -2*bcomp[2] * sinBvecSum;
    hcomp[5] = 2*acomp[1] * cosAvecSum + 2*bcomp[1] * cosBvecSum + m11;

    double hcomp5gamma = 6*acomp[1] + 6*bcomp[1] + m11;
    double bcoeff = sqrt(1. - bcomp[0]*bcomp[0]/(bcomp[1]*bcomp[1]));
    Bcomp[0] = hcomp[0] + bcomp[0]/bcomp[1] * (-hcomp[5] + hcomp5gamma);
    Bcomp[1] = bcoeff * hcomp[1];
    Bcomp[2] = bcoeff * hcomp[2];
    Bcomp[3] = bcoeff * hcomp[3];
}
