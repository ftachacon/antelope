/**
 * @brief Hamiltonian and dipole generator for Bi2Se3 surface states
 * @details Equation for this class comes from 10.1103/PhysRevA.103.023101
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class BieSe3surf : public BaseMaterial
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
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

BieSe3surf::BieSe3surf( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;   isDipoleZero = true;
    
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

    double klength = 4*pi/(3*a0);
    BZaxis = { klength*3./2,  klength*sqrt(3.)/2,     0,
                     0,       klength*sqrt(3.),       0,
                     0,                0,             0};


    // double kxMax = 2*pi/a0;  double kyMax = sqrt(3.)*pi/a0;
    // BZaxis = { 2*kxMax,       0,         0,
    //              kxMax,    2*kyMax,      0,
    //               0,          0,         0};
    // BZorigin = {-1.5*kxMax, -kyMax, 0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
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
    std::fill(xderBcomp, xderBcomp+4, 0.);
    std::fill(yderBcomp, yderBcomp+4, 0.);
    std::fill(xderHcomp, xderHcomp+5, 0.);
    std::fill(yderHcomp, yderHcomp+5, 0.);
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

        derxcosAvecSum += -vec_a[i][0] * sin(angle_a0[i]);
        derxsinAvecSum += vec_a[i][0] * cos(angle_a0[i]);
        derxcosBvecSum += -vec_b[i][0] * sin(angle_b0[i]);
        derxsinBvecSum += vec_b[i][0] * cos(angle_b0[i]);

        derycosAvecSum += -vec_a[i][1] * sin(angle_a0[i]);
        derysinAvecSum += vec_a[i][1] * cos(angle_a0[i]);
        derycosBvecSum += -vec_b[i][1] * sin(angle_b0[i]);
        derysinBvecSum += vec_b[i][1] * cos(angle_b0[i]);
    }

    double w = -2.*pi/3.;

    xderHcomp[0] = 2.*acomp[0] * derxcosAvecSum + 2.*bcomp[0] * derxcosBvecSum;
    
    xderHcomp[1] = -2.*acomp[3] * sin(w)*( vec_a[1][0]*cos(angle_a0[1]) - vec_a[2][0]*cos(angle_a0[2]) ) +
                   + 2.*bcomp[3] * (  vec_b[0][0]*cos(angle_b0[0])
                   + cos(w)*( vec_b[1][0]*cos(angle_b0[1]) + vec_b[2][0]*cos(angle_b0[2]) ) ) ;
    
    
    xderHcomp[2] = -2.*bcomp[3] * sin( w )*( vec_b[1][0]*cos(angle_b0[1]) - vec_b[2][0]*cos(angle_b0[2]) )
                    - 2.*acomp[3] * ( vec_a[0][0]*cos(angle_a0[0])
                     + cos( w )*( vec_a[1][0]*cos(angle_a0[1])
                                + vec_a[2][0]*cos(angle_a0[2]) )
                                     );
    
    
    xderHcomp[3] = 2.*acomp[2] * derxsinAvecSum;
    
    xderHcomp[4] = -2.*bcomp[2] * derxsinBvecSum;
    
    xderHcomp[5] = 2.*acomp[1] * derxcosAvecSum + 2.*bcomp[1] * derxcosBvecSum;

    
    
    yderHcomp[0] = 2.*acomp[0] * derycosAvecSum + 2.*bcomp[0] * derycosBvecSum;
    

    
    yderHcomp[1] = -2.*acomp[3] * sin(w)*( vec_a[1][1]*cos(angle_a0[1]) - vec_a[2][1]*cos(angle_a0[2]) ) +
    + 2.*bcomp[3] * (  vec_b[0][1]*cos(angle_b0[0])
                     + cos(w)*(  vec_b[1][1]*cos(angle_b0[1])
                               + vec_b[2][1]*cos(angle_b0[2]) ) ) ;
    
    
    yderHcomp[2] = -2.*bcomp[3] * sin( w )*( vec_b[1][1]*cos(angle_b0[1]) - vec_b[2][1]*cos(angle_b0[2]) )
    - 2.*acomp[3] * ( vec_a[0][1]*cos( angle_a0[0] )
                     + cos( w )*(  vec_a[1][1]*cos(angle_a0[1])
                                 + vec_a[2][1]*cos(angle_a0[2]) ) );
    
    
    yderHcomp[3] =  2.*acomp[2] * derysinAvecSum;
    yderHcomp[4] = -2.*bcomp[2] * derysinBvecSum;
    yderHcomp[5] = 2.*acomp[1] * derycosAvecSum + 2.*bcomp[1] * derycosBvecSum;


    double hcomp5gamma = 6.*(acomp[1] + bcomp[1]) + m11;
    double bcoeff = sqrt(1. - bcomp[0]*bcomp[0]/(bcomp[1]*bcomp[1]));


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
    std::cout << "============================================\n";
    std::cout << "Bi2Se3 surface - paramters are hard-coded\n";
    std::cout << "============================================\n";
}

void BieSe3surf::GenBcomp(std::array<double, Ndim> _kpoint)
{
    std::fill(Bcomp, Bcomp+4, 0.);
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
    
    double w = -2.*pi/3.;
    
    hcomp[0] = 2.*acomp[0] * cosAvecSum + 2.*bcomp[0] * cosBvecSum;
    
    hcomp[1] = -2.*acomp[3] * sin( w)*( sin(angle_a0[1]) - sin(angle_a0[2]) )
    + 2.*bcomp[3] * ( sin(angle_b0[0]) + cos(w)*(sin(angle_b0[1]) - sin(angle_b0[2]) ) ) ;
    
    hcomp[2] = -2.*bcomp[3] * sin( w )*( sin(angle_b0[1]) - sin(angle_b0[2]) )
    -2.*acomp[3] * ( sin(angle_a0[0]) + cos( w )*(sin(angle_a0[1]) + sin(angle_a0[2]) )  );
    
    hcomp[3] = 2.*acomp[2] * sinAvecSum;
    
    hcomp[4] = -2.*bcomp[2] * sinBvecSum;
    
    hcomp[5] = 2.*acomp[1] * cosAvecSum + 2.*bcomp[1] * cosBvecSum + m11;
    
    double hcomp5gamma = 6.*acomp[1] + 6.*bcomp[1] + m11;
    double bcoeff = sqrt(1. - bcomp[0]*bcomp[0]/(bcomp[1]*bcomp[1]));
    
    Bcomp[0] = hcomp[0] + bcomp[0]/bcomp[1] * (-hcomp[5] + hcomp5gamma);
    
    Bcomp[1] = bcoeff * hcomp[1];
    
    Bcomp[2] = bcoeff * hcomp[2];
    
    Bcomp[3] = bcoeff * hcomp[3];
}
