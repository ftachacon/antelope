/// Hamiltonian and dipole generator for MultiWSM model
/**
 * Some contents of this file comes from 'solidstructure.h' written by Alexis Chacon
 * This class assumes system is pure 2d - no slant plane
 * @author Dasol Kim
 * @author Alexis Agustín  Chacón Salazar
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>


class MultiWSM : public BaseMaterial
{
public:
    // Haldnae parameters
    double n, v, alphan, phi0, a0; 
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[4];
    double angle_a0, angle_b0;

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    MultiWSM( const libconfig::Setting *params );
    void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

MultiWSM::MultiWSM( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;   isDipoleZero = true;
    if (params->lookupValue("n", n)
        && params->lookupValue("v", v)
        && params->lookupValue("alphan", alphan)
        && params->lookupValue("phi0", phi0)
        && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;
    }
    else
    {
        std::cerr << "Some MultiWSM paramters are missing" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    SetBasis();

    double kxMax = pi/a0;//pi/sqrt(3.)/a0;  
    double kyMax = pi/a0;//2.*pi/3./a0; 
    double kzMax = pi/a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         2*kzMax};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void MultiWSM::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    _hstore[0] = Bcomp[0] + Bcomp[3];     _hstore[1] = Bcomp[1] - I*Bcomp[2];
    _hstore[2] = Bcomp[1] + I*Bcomp[2];   _hstore[3] = Bcomp[0] - Bcomp[3];
}

void MultiWSM::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[4];   
    double yderBcomp[4];
    double zderBcomp[4];
    double phip=atan2(_kpoint[1],_kpoint[0]);
    double psquare=_kpoint[0]*_kpoint[0]+_kpoint[1]*_kpoint[1];

    std::fill(Bcomp, Bcomp+4, 0.);
    std::fill(xderBcomp, xderBcomp+4, 0.);
    std::fill(yderBcomp, yderBcomp+4, 0.);
    std::fill(zderBcomp, zderBcomp+4, 0.);    
/*
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );

        xderBcomp[0] += vec_b[i][0] * sin( angle_b0 );      yderBcomp[0] += vec_b[i][1] * sin( angle_b0 );
        xderBcomp[1] += vec_a[i][0] * sin( angle_a0 );      yderBcomp[1] += vec_a[i][1] * sin( angle_a0 );
        xderBcomp[2] += vec_a[i][0] * cos( angle_a0 );      yderBcomp[2] += vec_a[i][1] * cos( angle_a0 );
        xderBcomp[3] += vec_b[i][0] * cos( angle_b0 );      yderBcomp[3] += vec_b[i][1] * cos( angle_b0 );
    }

    Bcomp[0] *= 2.* v * cos( phi0 );       Bcomp[1] *= n;
    Bcomp[2] *= +n;                        Bcomp[3] = alphan - 2.*v*sin( phi0 )*Bcomp[3]; */

    xderBcomp[0] = 0.;     
    yderBcomp[0] = 0.;
    zderBcomp[0] = 0.;

    xderBcomp[1] = n * alphan * _kpoint[0] * cos( n * phip )*pow( psquare , n/2.-1 )
                 + n * alphan * sin( n * phip )*pow( psquare , n/2. )*_kpoint[1]/psquare;

    yderBcomp[1] =  n * alphan * _kpoint[1] * cos( n * phip )*pow( psquare , n/2.-1 )
                  - n * alphan * sin( n * phip ) * pow( psquare , n/2. )*_kpoint[0]/psquare;
    zderBcomp[1] = 0.;

    xderBcomp[2] = +n * alphan * _kpoint[0] * sin( n * phip )*pow( psquare , n/2.-1 ) 
                  - n * alphan * _kpoint[1] * cos( n * phip )*pow( psquare , n/2. )/psquare;
                  
    yderBcomp[2] = +n * alphan * _kpoint[1] * sin( n * phip )*pow( psquare , n/2.-1 ) 
                   +n * alphan * _kpoint[0] * cos( n * phip )*pow( psquare , n/2. )/psquare;
    zderBcomp[2] = 0.;

    xderBcomp[3] = 0;     
    yderBcomp[3] = 0;
    zderBcomp[3] = v;

    _jstore[0][0] = -(xderBcomp[0] + xderBcomp[3]);      
    _jstore[1][0] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[2][0] = -(zderBcomp[0] + zderBcomp[3]);    

    _jstore[0][1] = -(xderBcomp[1] - I*xderBcomp[2]);    
    _jstore[1][1] = -(yderBcomp[1] - I*yderBcomp[2]);
    _jstore[2][1] = -(zderBcomp[1] - I*zderBcomp[2]);

    _jstore[0][2] = -(xderBcomp[1] + I*xderBcomp[2]);    
    _jstore[1][2] = -(yderBcomp[1] + I*yderBcomp[2]);
    _jstore[2][2] = -(zderBcomp[1] + I*zderBcomp[2]);

    _jstore[0][3] = -(xderBcomp[0] - xderBcomp[3]);     
    _jstore[1][3] = -(yderBcomp[0] - yderBcomp[3]);
    _jstore[2][3] = -(zderBcomp[0] - zderBcomp[3]);    

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > MultiWSM::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void MultiWSM::SetBasis()
{
    vec_a[0][0] = 0.;                   vec_a[0][1] = a0;
    vec_a[1][0] = -sqrt(3.)/2.*a0;       vec_a[1][1] = -0.5*a0;
    vec_a[2][0] = sqrt(3.)/2.*a0;        vec_a[2][1] = -0.5*a0;

    vec_b[0][0] = sqrt(3.)*a0;          vec_b[0][1] = 0.;
    vec_b[1][0] = -sqrt(3.)/2.*a0;       vec_b[1][1] = 3./2.*a0;
    vec_b[2][0] = -sqrt(3.)/2.*a0;       vec_b[2][1] = -3./2.*a0;
}

void MultiWSM::PrintMaterialInformation()
{
    std::cout << "============================================\n";
    std::cout << "MultiWSM Model paramters\n";
    std::cout << "a0         = " << a0 << " a.u. \n";
    std::cout << "n         = " << n << " a.u. \n";
    std::cout << "v         = " << v << " a.u. \n";
    std::cout << "alphan         = " << alphan << " a.u. \n";
    std::cout << "phi0         = " << phi0 << " rad \n";
    std::cout << "============================================\n";
}

void MultiWSM::GenBcomp(std::array<double, Ndim> _kpoint)
{
    std::fill(Bcomp, Bcomp+4, 0.);
   /* for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );
    }

    Bcomp[0] *= 2.* t2 * cos( phi0 );       Bcomp[1] *= t1;
    Bcomp[2] *= t1;                         Bcomp[3] = alphan - 2.*t2*sin( phi0 )*Bcomp[3];
    */
   Bcomp[0] = 0.;
   Bcomp[1] = alphan* pow( _kpoint[0]*_kpoint[0] + _kpoint[1]*_kpoint[1]  ,n/2.)*cos(n* atan2( _kpoint[1],_kpoint[0] ) ) ;
   Bcomp[2] = alphan* pow( _kpoint[0]*_kpoint[0] + _kpoint[1]*_kpoint[1]  ,n/2.)*sin(n* atan2( _kpoint[1],_kpoint[0] ) ) ;
   Bcomp[3] = v * _kpoint[2];
}
