/// Hamiltonian and dipole generator for 2-band square lattice including Wilson mass term
/**
 * This class assumes system is pure 2d - no slant plane
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class WilsonMass : public BaseMaterial
{
public:
    // Wilson mass parameters
    double t, delta, mu, a0, a1;

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[4];

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    WilsonMass( const libconfig::Setting *params );
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

WilsonMass::WilsonMass( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;
    if (params->lookupValue("t", t)
        && params->lookupValue("delta", delta)
        && params->lookupValue("mu", mu) 
        && params->lookupValue("a0", a0) )
    {
        if (!params->lookupValue("a1", a1)) a1 = a0;
    }
    else
    {
        std::cerr << "Some Haldane paramters are missing" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double kxMax = pi/a0;  double kyMax = pi/a1;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void WilsonMass::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    _hstore[0] = Bcomp[0] + Bcomp[3];     _hstore[1] = Bcomp[1] - I*Bcomp[2];
    _hstore[2] = Bcomp[1] + I*Bcomp[2];   _hstore[3] = Bcomp[0] - Bcomp[3];
}

void WilsonMass::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[4];    double yderBcomp[4];

    Bcomp[0] = 0.;          
    Bcomp[1] = delta*sin(a0*_kpoint[0]);       Bcomp[2] = delta*sin(a1*_kpoint[1]);
    Bcomp[3] = -2.*t*(cos(a0*_kpoint[0]) + cos(a1*_kpoint[1])) - mu;

    xderBcomp[0] = 0.;                             yderBcomp[0] = 0.;
    xderBcomp[1] = a0*delta*cos(a0*_kpoint[0]);    yderBcomp[1] = 0.;
    xderBcomp[2] = 0.;                             yderBcomp[2] = a1*delta*cos(a1*_kpoint[1]);
    xderBcomp[3] = 2.*t *a0* sin(a0*_kpoint[0]);   yderBcomp[3] = 2.*t *a1* sin(a1*_kpoint[1]);

    _jstore[0][0] = -(xderBcomp[0] + xderBcomp[3]);      _jstore[1][0] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[0][1] = -(xderBcomp[1] - I*xderBcomp[2]);    _jstore[1][1] = -(yderBcomp[1] - I*yderBcomp[2]);
    _jstore[0][2] = -(xderBcomp[1] + I*xderBcomp[2]);    _jstore[1][2] = -(yderBcomp[1] + I*yderBcomp[2]);
    _jstore[0][3] = -(xderBcomp[0] - xderBcomp[3]);      _jstore[1][3] = -(yderBcomp[0] - yderBcomp[3]);
    

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > WilsonMass::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void WilsonMass::PrintMaterialInformation()
{
    std::cout << "============================================\n";
    std::cout << "WilsonMass Model paramters\n";
    std::cout << "t          = " << t << " a.u. \n";
    std::cout << "mu         = " << mu << " a.u. \n";
    std::cout << "delta      = " << delta << " a.u. \n";
    std::cout << "a0         = " << a0 << " a.u. \n";
    std::cout << "a1         = " << a1 << " a.u. \n";
    std::cout << "============================================\n";
}

void WilsonMass::GenBcomp(std::array<double, Ndim> _kpoint)
{
    Bcomp[0] = 0.;          
    Bcomp[1] = delta*sin(a0*_kpoint[0]);       Bcomp[2] = delta*sin(a1*_kpoint[1]);
    Bcomp[3] = -2.*t*(cos(a0*_kpoint[0]) + cos(a1*_kpoint[1])) - mu;
}
