/**
 * @brief Hamiltonian and dipole generator for Kane-Mele model
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class KaneMele : public BaseMaterial
{
public:
    // Haldnae parameters
    double t1, t2, M0, phi0, a0; 
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[4];
    double angle_a0, angle_b0;

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    KaneMele( const libconfig::Setting *params );
    void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpointm, double _phi0);
};

KaneMele::KaneMele( const libconfig::Setting *params )
{
    Nband = 4;  Nval = 2;
    if (params->lookupValue("t1", t1)
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
        cerr << "Some KaneMele paramters are missing" << endl;
        exit(EXIT_FAILURE);
    }

    SetBasis();

    double kxMax = pi/sqrt(3.)/a0;  double kyMax = 2.*pi/3./a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void KaneMele::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    fill(_hstore, _hstore+Nband*Nband, 0.);

    GenBcomp(_kpoint, phi0);

    _hstore[0] = Bcomp[0] + Bcomp[3];     _hstore[1] = Bcomp[1] - I*Bcomp[2];
    _hstore[4] = Bcomp[1] + I*Bcomp[2];   _hstore[5] = Bcomp[0] - Bcomp[3];

    GenBcomp(_kpoint, -phi0);

    _hstore[10] = Bcomp[0] + Bcomp[3];     _hstore[11] = Bcomp[1] - I*Bcomp[2];
    _hstore[14] = Bcomp[1] + I*Bcomp[2];   _hstore[15] = Bcomp[0] - Bcomp[3];
}

void KaneMele::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[4];    double yderBcomp[4];
    fill(Bcomp, Bcomp+4, 0.);
    fill(xderBcomp, xderBcomp+4, 0.);
    fill(yderBcomp, yderBcomp+4, 0.);
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        //Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        //Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );

        xderBcomp[0] += vec_b[i][0] * sin( angle_b0 );      yderBcomp[0] += vec_b[i][1] * sin( angle_b0 );
        xderBcomp[1] += vec_a[i][0] * sin( angle_a0 );      yderBcomp[1] += vec_a[i][1] * sin( angle_a0 );
        xderBcomp[2] += vec_a[i][0] * cos( angle_a0 );      yderBcomp[2] += vec_a[i][1] * cos( angle_a0 );
        xderBcomp[3] += vec_b[i][0] * cos( angle_b0 );      yderBcomp[3] += vec_b[i][1] * cos( angle_b0 );
    }

    //Bcomp[0] *= 2.* t2 * cos( phi0 );       Bcomp[1] *= t1;
    //Bcomp[2] *= t1;                         Bcomp[3] = M0 - 2.*t2*sin( phi0 )*Bcomp[3];

    xderBcomp[0] *= -2.*t2 * cos(phi0);     yderBcomp[0] *= -2.*t2 * cos(phi0);
    xderBcomp[1] *= -t1;                    yderBcomp[1] *= -t1;
    xderBcomp[2] *= t1;                     yderBcomp[2] *= t1;
    xderBcomp[3] *= -2.*t2 * sin(phi0);     yderBcomp[3] *= -2.*t2 * sin(phi0);

    _jstore[0][0] = -(xderBcomp[0] + xderBcomp[3]);      _jstore[1][0] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[0][1] = -(xderBcomp[1] - I*xderBcomp[2]);    _jstore[1][1] = -(yderBcomp[1] - I*yderBcomp[2]);
    _jstore[0][4] = -(xderBcomp[1] + I*xderBcomp[2]);    _jstore[1][4] = -(yderBcomp[1] + I*yderBcomp[2]);
    _jstore[0][5] = -(xderBcomp[0] - xderBcomp[3]);      _jstore[1][5] = -(yderBcomp[0] - yderBcomp[3]);

    // phi0 --> -phi0, sin(phi0) --> -sin(phi0)

    _jstore[0][10] = -(xderBcomp[0] - xderBcomp[3]);      _jstore[1][10] = -(yderBcomp[0] - yderBcomp[3]);
    _jstore[0][11] = -(xderBcomp[1] - I*xderBcomp[2]);    _jstore[1][11] = -(yderBcomp[1] - I*yderBcomp[2]);
    _jstore[0][14] = -(xderBcomp[1] + I*xderBcomp[2]);    _jstore[1][14] = -(yderBcomp[1] + I*yderBcomp[2]);
    _jstore[0][15] = -(xderBcomp[0] + xderBcomp[3]);      _jstore[1][15] = -(yderBcomp[0] + yderBcomp[3]);
}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > KaneMele::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void KaneMele::SetBasis()
{
    vec_a[0][0] = 0.;                   vec_a[0][1] = a0;
    vec_a[1][0] = -sqrt(3.)/2*a0;       vec_a[1][1] = -0.5*a0;
    vec_a[2][0] = sqrt(3.)/2*a0;        vec_a[2][1] = -0.5*a0;

    vec_b[0][0] = sqrt(3.)*a0;          vec_b[0][1] = 0.;
    vec_b[1][0] = -sqrt(3.)/2*a0;       vec_b[1][1] = 3./2.*a0;
    vec_b[2][0] = -sqrt(3.)/2*a0;       vec_b[2][1] = -3./2.*a0;
}

void KaneMele::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << "Kane-Mele Model paramters\n";
    cout << "a0         = " << a0 << " a.u. \n";
    cout << "t1         = " << t1 << " a.u. \n";
    cout << "t2         = " << t2 << " a.u. \n";
    cout << "M0         = " << M0 << " a.u. \n";
    cout << "phi0         = " << phi0 << " rad \n";
    cout << "============================================\n";
}

void KaneMele::GenBcomp(std::array<double, Ndim> _kpoint, double _phi0)
{
    fill(Bcomp, Bcomp+4, 0.);
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );
    }

    Bcomp[0] *= 2.* t2 * cos( phi0 );       Bcomp[1] *= t1;
    Bcomp[2] *= t1;                         Bcomp[3] = M0 - 2.*t2*sin( _phi0 )*Bcomp[3];
}
