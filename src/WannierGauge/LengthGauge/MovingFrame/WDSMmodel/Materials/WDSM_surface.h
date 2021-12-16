/// Hamiltonian and dipole generator for WDSM model
/**
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class WDSM_surface : public BaseMaterial
{
public:
    // Haldnae parameters
    double t1, t2, M0, phi0, a0;
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[6];
    double angle_a0, angle_b0;

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    WDSM_surface( const libconfig::Setting *params );
    // void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

WDSM_surface::WDSM_surface( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;
    if (params->lookupValue("t1", t1)
        && params->lookupValue("t2", t2)
        && params->lookupValue("Mt2", M0)
        && params->lookupValue("phi0", phi0)
        && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;
        M0 = 0;

    }
    else
    {
        cerr << "Some WDSM paramters are missing" << endl;
        exit(EXIT_FAILURE);
    }
    // int M0=0;
    // SetBasis();

    double kxMax = pi/a0;  double kyMax = pi/a0; double kzMax = pi/a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void WDSM_surface::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    _hstore[0] = Bcomp[0] ;   _hstore[1] = Bcomp[1]+I*Bcomp[4];
    _hstore[2] = Bcomp[2]+I*Bcomp[5] ;   _hstore[3] = Bcomp[3];
}

void WDSM_surface::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[6];    double yderBcomp[6]; double zderBcomp[6];
    fill(Bcomp, Bcomp+6, 0.);
    fill(xderBcomp, xderBcomp+6, 0.);
    fill(yderBcomp, yderBcomp+6, 0.);
    fill(zderBcomp, yderBcomp+6, 0.);

    for (int i = 0; i < 2; ++i)
    {
        Bcomp[0] += -cos( _kpoint[i]*a0 );
        Bcomp[3] += + cos( _kpoint[i]*a0 );
    }
    Bcomp[1] += sin( _kpoint[0]*a0 );
    Bcomp[2] += sin( _kpoint[0]*a0 );
    Bcomp[0] += 2-M0;
    Bcomp[3] += -2+ M0;
    Bcomp[4] += - sin( _kpoint[1]*a0 );
    Bcomp[5] += sin( _kpoint[1]*a0 );


    xderBcomp[0] += -cos( _kpoint[0]*a0);
    xderBcomp[1] += sin( _kpoint[0]*a0);
    xderBcomp[2] += sin( _kpoint[0]*a0);
    xderBcomp[3] += cos( _kpoint[0]*a0);

    yderBcomp[0] += - cos( _kpoint[1]*a0);
    yderBcomp[1] += - sin( _kpoint[1]*a0);
    yderBcomp[2] += sin( _kpoint[1]*a0);
    yderBcomp[3] += - cos( _kpoint[1]*a0);


    _jstore[0][0] = - xderBcomp[0];
   _jstore[1][0] = - yderBcomp[0];

    _jstore[0][1] = - xderBcomp[1];
   _jstore[1][1] = - I*yderBcomp[1];
 
   _jstore[0][2] = - xderBcomp[2];
   _jstore[1][2] = -  I*yderBcomp[2];

    _jstore[0][3] = - xderBcomp[3];
   _jstore[1][3] = - yderBcomp[3];


}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > WDSM_surface::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

// void WDSM::SetBasis()
// {
//     vec_a[0][0] = 1.;                   vec_a[0][1] = a0;
//     vec_a[1][0] = 0.;                   vec_a[1][1] = -0.5*a0;
//     vec_a[2][0] = 0. ;                   vec_a[2][1] = -0.5*a0;
//
//     vec_b[0][0] = 0.0;                   vec_b[0][1] = 0.;
//     vec_b[1][0] = 1.0;       vec_b[1][1] = 3./2.*a0;
//     vec_b[2][0] = 0.0;       vec_b[2][1] = -3./2.*a0;
// }

void WDSM_surface::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << "WDSM Model paramters\n";
    cout << "a0         = " << a0 << " a.u. \n";
    cout << "t1         = " << t1 << " a.u. \n";
    cout << "t2         = " << t2 << " a.u. \n";
    cout << "M0         = " << M0 << " a.u. \n";
    cout << "phi0         = " << phi0 << " rad \n";
    cout << "============================================\n";
}

void WDSM_surface::GenBcomp(std::array<double, Ndim> _kpoint)
{
    fill(Bcomp, Bcomp+6, 0.);
    for (int i = 0; i < 2; ++i)
    {
        Bcomp[0] += -cos( _kpoint[i]*a0 );
        Bcomp[3] += + cos( _kpoint[i]*a0 );
    }
    Bcomp[1] += sin( _kpoint[0]*a0 );
    Bcomp[2] += sin( _kpoint[0]*a0 ) ;
    Bcomp[0] += 2-M0-1;
    Bcomp[3] += -2+ M0+1;
    Bcomp[4] += - sin( _kpoint[1]*a0 );
    Bcomp[5] += sin( _kpoint[1]*a0 );

}
