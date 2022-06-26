/**
 * @brief  Hamiltonian and dipole generator for XYH lattice
 * @details Model comes from Phys. Rev. B 105, 075421 (2022).
 *          XYHlattice = px,py orbitals in Honeycomb lattice.
 *          This model only includes NN interaction.
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class XYHlattice : public BaseMaterial
{
public:
    // a0 is lattice constant (do not follow Haldane's specific notation a0=a0/sqrt(3))
    // t_parallel --> tl, t_perpendicular --> tr
    double a0, tl, tr, lambda; 
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    //double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    // double Bcomp[4];
    complex tempphase[3];
    complex alpha_plus, alpha_minus, beta;
    double angle_a0, angle_b0;

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    XYHlattice( const libconfig::Setting *params );
    void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpointm);
};

XYHlattice::XYHlattice( const libconfig::Setting *params )
{
    Nband = 4;  Nval = 2;
    if (params->lookupValue("tl", tl)
        && params->lookupValue("tr", tr)
        && params->lookupValue("lambda", lambda)
        && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;
    }
    else
    {
        cerr << "Some XYHlattice paramters are missing" << endl;
        exit(EXIT_FAILURE);
    }

    SetBasis();

    double kxMax = pi/a0;  double kyMax = 2.*pi/sqrt(3.)/a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void XYHlattice::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    fill(_hstore, _hstore+Nband*Nband, 0.);

    GenBcomp(_kpoint);

    // diagonal
    _hstore[0] = lambda;    _hstore[5] = lambda;    _hstore[10] = -lambda;  _hstore[15] = -lambda;

    // block of off-diagonal
    _hstore[2] = (3./4 * tl + 1./4 * tr) * alpha_plus + tr * beta;      _hstore[3] = sqrt(3.)/4 * (tl - tr) * alpha_minus;
    _hstore[6] = _hstore[3];                                            _hstore[7] = (1./4 * tl + 3./4 * tr) * alpha_plus + tl * beta;

    _hstore[8] = conj(_hstore[2]);      _hstore[9] = conj(_hstore[6]);
    _hstore[12] = conj(_hstore[3]);     _hstore[13] = conj(_hstore[7]);
}

void XYHlattice::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    // Bcomp - exp(-ik a1), exp(-ik a2), exp(-ik a3)
    // alpha_pm = Bcomp[0] \pm Bcomp[1], beta = Bcomp[2]
    complex dx_Bcomp[3];    complex dy_Bcomp[3];
    // fill(dx_Bcomp, dx_Bcomp+3, 0.);
    // fill(dy_Bcomp, dy_Bcomp+3, 0.);
    for (int i = 0; i < 3; ++i)
    {
        tempphase[i] = exp(-I* (_kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1]) );
        dx_Bcomp[i] = -I*vec_a[i][0] * tempphase[i];
        dy_Bcomp[i] = -I*vec_a[i][1] * tempphase[i];
    }

    //Bcomp[0] *= 2.* t2 * cos( phi0 );       Bcomp[1] *= t1;
    //Bcomp[2] *= t1;                         Bcomp[3] = M0 - 2.*t2*sin( phi0 )*Bcomp[3];

    fill(_jstore[0], _jstore[0] + 16, 0.);
    fill(_jstore[1], _jstore[0] + 16, 0.);

    // Acomp --> alpha_plus, alpha_minus, beta
    complex dx_Acomp[3];    complex dy_Acomp[3];

    dx_Acomp[0] = dx_Bcomp[0] + dx_Bcomp[1];
    dx_Acomp[1] = dx_Bcomp[0] - dx_Bcomp[1];
    dx_Acomp[2] = dx_Bcomp[2];

    dy_Acomp[0] = dy_Bcomp[0] + dy_Bcomp[1];
    dy_Acomp[1] = dy_Bcomp[0] - dy_Bcomp[1];
    dy_Acomp[2] = dy_Bcomp[2];
    // block of off-diagonal
    _jstore[0][2] = (3./4 * tl + 1./4 * tr) * dx_Acomp[0] + tr * dx_Acomp[2];      
    _jstore[0][3] = sqrt(3.)/4 * (tl - tr) * dx_Acomp[1];
    _jstore[0][6] = _jstore[0][3];                                            
    _jstore[0][7] = (1./4 * tl + 3./4 * tr) * dx_Acomp[0] + tl * dx_Acomp[2];

    _jstore[0][8] = conj(_jstore[0][2]);
    _jstore[0][9] = conj(_jstore[0][6]);
    _jstore[0][12] = conj(_jstore[0][3]);     
    _jstore[0][13] = conj(_jstore[0][7]);

    _jstore[1][2] = (3./4 * tl + 1./4 * tr) * dy_Acomp[0] + tr * dy_Acomp[2];      
    _jstore[1][3] = sqrt(3.)/4 * (tl - tr) * dy_Acomp[1];
    _jstore[1][6] = _jstore[1][3];                                            
    _jstore[1][7] = (1./4 * tl + 3./4 * tr) * dy_Acomp[0] + tl * dy_Acomp[2];

    _jstore[1][8] = conj(_jstore[1][2]);
    _jstore[1][9] = conj(_jstore[1][6]);
    _jstore[1][12] = conj(_jstore[1][3]);     
    _jstore[1][13] = conj(_jstore[1][7]);

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > XYHlattice::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void XYHlattice::SetBasis()
{
    vec_a[0][0] = 0.;                   vec_a[0][1] = a0/sqrt(3.);
    vec_a[1][0] = -1./2*a0;       vec_a[1][1] = -0.5*a0/sqrt(3.);
    vec_a[2][0] = 1./2*a0;        vec_a[2][1] = -0.5*a0/sqrt(3.);

    // vec_b[0][0] = sqrt(3.)*a0;          vec_b[0][1] = 0.;
    // vec_b[1][0] = -sqrt(3.)/2*a0;       vec_b[1][1] = 3./2.*a0;
    // vec_b[2][0] = -sqrt(3.)/2*a0;       vec_b[2][1] = -3./2.*a0;
}

void XYHlattice::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << "Kane-Mele Model paramters\n";
    cout << "a0         = " << a0 << " a.u. \n";
    cout << "tl         = " << tl << " a.u. \n";
    cout << "tr         = " << tr << " a.u. \n";
    cout << "lambda     = " << lambda << " a.u. \n";
    cout << "============================================\n";
}

void XYHlattice::GenBcomp(std::array<double, Ndim> _kpoint)
{
    for (int i = 0; i < 3; ++i)
    {
        tempphase[i] = exp(-I* (_kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1]) );
    }

    alpha_plus = tempphase[0] + tempphase[1];
    alpha_minus = tempphase[0] - tempphase[1];
    beta = tempphase[2];
}
