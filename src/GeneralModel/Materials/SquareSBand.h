/// Hamiltonian and dipole generator for simple 2D square lattice model
/**
 * 2D square lattice with s orbital only. 
 * Inversion symmetric, time-reversal symmetric, band gap at the corner of the edge.
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class SquareSBand : public BaseMaterial
{
public:
    // parameters
    double t1, M0, a0; 

    // internal variables for storing temporary values
    double Bcomp[4];

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    SquareSBand( const libconfig::Setting *params );
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

SquareSBand::SquareSBand( const libconfig::Setting *params )
{
    Nband = 2;  Nval = 1;   isDipoleZero = true;
    if (params->lookupValue("t1", t1)
        && params->lookupValue("M0", M0)
        && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;
    }
    else
    {
        std::cerr << "Some SquareSBand paramters are missing" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double kxMax = pi/a0;  double kyMax = pi/a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};
}

void SquareSBand::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    _hstore[0] = Bcomp[0] + Bcomp[3];     _hstore[1] = Bcomp[1] - I*Bcomp[2];
    _hstore[2] = Bcomp[1] + I*Bcomp[2];   _hstore[3] = Bcomp[0] - Bcomp[3];
}

void SquareSBand::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    std::fill(_jstore[0], _jstore[0]+4, 0.);
    std::fill(_jstore[1], _jstore[1]+4, 0.);

    // Note J = -P (minus from electron charge sign)
    _jstore[0][1] = 2*t1*a0*sin(a0*_kpoint[0]);
    _jstore[1][1] = 2*t1*a0*sin(a0*_kpoint[1]); 

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > SquareSBand::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void SquareSBand::PrintMaterialInformation()
{
    std::cout << "============================================\n";
    std::cout << "SquareSBand Model paramters\n";
    std::cout << "a0         = " << a0 << " a.u. \n";
    std::cout << "t1         = " << t1 << " a.u. \n";
    std::cout << "M0         = " << M0 << " a.u. \n";
    std::cout << "============================================\n";
}

void SquareSBand::GenBcomp(std::array<double, Ndim> _kpoint)
{
    std::fill(Bcomp, Bcomp+4, 0.);

    Bcomp[0] = 0.;       Bcomp[1] = 2*t1*(cos(a0*_kpoint[0]) + cos(a0*_kpoint[1]));
    Bcomp[2] = 0;        Bcomp[3] = M0;
}
