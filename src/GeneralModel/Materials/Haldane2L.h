/**
 * @brief Hamiltonian and dipole generator for Haldane2L model (Bilayer Haldane) 
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class Haldane2L : public BaseMaterial
{
public:
    // Haldnae parameters
    // layer top (0) and bottom (1) have seprate parameters. Except a0 (it should be implemented as a seprate model.)
    // t11 - interlayer orthogonal param (hopping between directly above and below atoms)
    double t1[2], t2[2], M0[2], phi0[2], a0, t11; 
    complex h_inter[4];
    double B2sign[2];     // B2 sign for each layer (-1 for complex conjugate)
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    // internal variables for storing temporary values
    double Bcomp[4];
    double Bcomp_bottom[4];
    double angle_a0, angle_b0;

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    // polytype - 2H or 3R
    std::string polytype;
    // check whether top and bottom have same paramters.
    bool is_top_bottom_same;

    Haldane2L( const libconfig::Setting *params );
    void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

Haldane2L::Haldane2L( const libconfig::Setting *params )
{
    Nband = 4;  Nval = 2;   isDipoleZero = true;
    if (params->lookupValue("t1", t1[0])
        && params->lookupValue("t2", t2[0])
        //&& params->lookupValue("Mt2", M0)
        && params->lookupValue("phi0", phi0[0])
        && params->lookupValue("a0", a0)
        && params->lookupValue("t11", t11)
        && params->lookupValue("polytype", polytype) )
    {
        a0 /= au_angstrom;
        if (!params->lookupValue("M0", M0[0]))
        {
            if (params->lookupValue("Mt2",M0[0]))
            {
                M0[0] *= t2[0];
            }
            else
            {
                std::cerr << "M0 param missing" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }

        // bottom layer paramters default
        is_top_bottom_same = true;
        t1[1] = t1[0];  t2[1] = t2[0];
        phi0[1] = phi0[0]; M0[1] = M0[0];

        // bottom layer paramters
        if (params->lookupValue("t1_bottom", t1[1])
            || params->lookupValue("t2_bottom", t2[1])
            || params->lookupValue("phi0_bottom", phi0[1]))
        {
            is_top_bottom_same = false;
        }
        if (params->lookupValue("M0_bottom", M0[1]))
        {
            is_top_bottom_same = false;
        }
        else if (params->lookupValue("Mt2_bottom",M0[1]))
        {
            M0[1] *= t2[1];
            is_top_bottom_same = false;
        }
    }
    else
    {
        std::cerr << "Some Haldane2L paramters are missing" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // defulat B2 sign
    B2sign[0] = 1.;  B2sign[1] = 1.;

    // Notation: A, B position in matrix is unchanged - Atop, Btop, Abottom, Bbottom
    // Stacking follows PhysRevB.83.235312 (ribeiro2011)
    if (polytype == "AA")
    {
        h_inter[0] = t11;   h_inter[1] = 0.;
        h_inter[2] = 0.;    h_inter[3] = t11;
    }
    else if (polytype == "AB" || polytype == "AB1")
    {
        h_inter[0] = 0.;    h_inter[1] = t11;
        h_inter[2] = 0.;    h_inter[3] = 0.;
    }
    else if (polytype == "AB2" || polytype == "3R")
    {
        h_inter[0] = 0.;     h_inter[1] = 0.;
        h_inter[2] = t11;    h_inter[3] = 0.;
    }
    else if (polytype == "AA'" || polytype == "2H")
    {
        h_inter[0] = 0.;     h_inter[1] = t11;
        h_inter[2] = t11;    h_inter[3] = 0.;
        B2sign[1] = -1.;
    }
    else if (polytype == "A'B")
    {
        h_inter[0] = 0.;    h_inter[1] = 0.;
        h_inter[2] = 0.;    h_inter[3] = t11;
        B2sign[0] = -1.;
    }
    else if (polytype == "AB'")
    {
        h_inter[0] = t11;    h_inter[1] = 0.;
        h_inter[2] = 0.;     h_inter[3] = 0.;
        B2sign[1] = -1.;
    }
    else
    {
        std::cerr << "Undefined polytype" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    SetBasis();

    //double kxMax = pi/sqrt(3.)/a0;  double kyMax = 2.*pi/3./a0;
    double klength = 4*pi/(3*sqrt(3.)*a0);
    BZaxis = { klength*3./2,  klength*sqrt(3.)/2,     0,
                     0,       klength*sqrt(3.),       0,
                     0,                0,             0};

    // using reciprocal vector
    // corresponding to: (1, 0) *sqrt(3)a0,   (-1/2, sqrt(3)/2) *sqrt(3)a0
    // BZaxis = { 2*pi/sqrt(3.)/a0,     2*pi/3./a0,        0,
    //                    0,            4*pi/3./a0,        0,
    //                    0,               0,              0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void Haldane2L::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);

    // block diagonal part (intra-layer part)
    _hstore[0] = Bcomp[0] + Bcomp[3];               _hstore[1] = Bcomp[1] - B2sign[0]*I*Bcomp[2];
    _hstore[4] = Bcomp[1] + B2sign[0]*I*Bcomp[2];   _hstore[5] = Bcomp[0] - Bcomp[3];

    _hstore[10] = Bcomp_bottom[0] + Bcomp_bottom[3];              _hstore[11] = Bcomp_bottom[1] - B2sign[1]*I*Bcomp_bottom[2];
    _hstore[14] = Bcomp_bottom[1] + B2sign[1]*I*Bcomp_bottom[2];  _hstore[15] = Bcomp_bottom[0] - Bcomp_bottom[3];

    // inter-layer interaction part
    _hstore[2] = h_inter[0];    _hstore[3] = h_inter[1];
    _hstore[6] = h_inter[2];    _hstore[7] = h_inter[3];

    _hstore[8] = conj(h_inter[0]);    _hstore[9] = conj(h_inter[2]);
    _hstore[12] = conj(h_inter[1]);   _hstore[13] = conj(h_inter[3]);
}

void Haldane2L::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[4] = {0};    double yderBcomp[4] = {0};
    //std::fill(Bcomp, Bcomp+4, 0.);
    // std::fill(xderBcomp, xderBcomp+4, 0.);
    // std::fill(yderBcomp, yderBcomp+4, 0.);
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        // Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        // Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );

        xderBcomp[0] += vec_b[i][0] * sin( angle_b0 );      yderBcomp[0] += vec_b[i][1] * sin( angle_b0 );
        xderBcomp[1] += vec_a[i][0] * sin( angle_a0 );      yderBcomp[1] += vec_a[i][1] * sin( angle_a0 );
        xderBcomp[2] += vec_a[i][0] * cos( angle_a0 );      yderBcomp[2] += vec_a[i][1] * cos( angle_a0 );
        xderBcomp[3] += vec_b[i][0] * cos( angle_b0 );      yderBcomp[3] += vec_b[i][1] * cos( angle_b0 );
    }

    // Bcomp[0] *= 2.* t2[0] * cos( phi0[0] );    Bcomp[1] *= t1[0];
    // Bcomp[2] *= t1[0];                         Bcomp[3] = M0[0] - 2.*t2[0]*sin( phi0[0] )*Bcomp[3];

    xderBcomp[0] *= -2.*t2[0] * cos(phi0[0]);  yderBcomp[0] *= -2.*t2[0] * cos(phi0[0]);
    xderBcomp[1] *= -t1[0];                    yderBcomp[1] *= -t1[0];
    xderBcomp[2] *= t1[0];                     yderBcomp[2] *= t1[0];
    xderBcomp[3] *= -2.*t2[0] * sin(phi0[0]);  yderBcomp[3] *= -2.*t2[0] * sin(phi0[0]);

    _jstore[0][0] = -(xderBcomp[0] + xderBcomp[3]);               _jstore[1][0] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[0][1] = -(xderBcomp[1] - B2sign[0]*I*xderBcomp[2]);   _jstore[1][1] = -(yderBcomp[1] - B2sign[0]*I*yderBcomp[2]);
    _jstore[0][4] = -(xderBcomp[1] + B2sign[0]*I*xderBcomp[2]);   _jstore[1][4] = -(yderBcomp[1] + B2sign[0]*I*yderBcomp[2]);
    _jstore[0][5] = -(xderBcomp[0] - xderBcomp[3]);               _jstore[1][5] = -(yderBcomp[0] - yderBcomp[3]);


    // bottom layer
    //std::fill(Bcomp, Bcomp+4, 0.);
    std::fill(xderBcomp, xderBcomp+4, 0.);
    std::fill(yderBcomp, yderBcomp+4, 0.);
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];

        xderBcomp[0] += vec_b[i][0] * sin( angle_b0 );      yderBcomp[0] += vec_b[i][1] * sin( angle_b0 );
        xderBcomp[1] += vec_a[i][0] * sin( angle_a0 );      yderBcomp[1] += vec_a[i][1] * sin( angle_a0 );
        xderBcomp[2] += vec_a[i][0] * cos( angle_a0 );      yderBcomp[2] += vec_a[i][1] * cos( angle_a0 );
        xderBcomp[3] += vec_b[i][0] * cos( angle_b0 );      yderBcomp[3] += vec_b[i][1] * cos( angle_b0 );
    }

    xderBcomp[0] *= -2.*t2[1] * cos(phi0[1]);  yderBcomp[0] *= -2.*t2[1] * cos(phi0[1]);
    xderBcomp[1] *= -t1[1];                    yderBcomp[1] *= -t1[1];
    xderBcomp[2] *= t1[1];                     yderBcomp[2] *= t1[1];
    xderBcomp[3] *= -2.*t2[1] * sin(phi0[1]);  yderBcomp[3] *= -2.*t2[1] * sin(phi0[1]);

    _jstore[0][10] = -(xderBcomp[0] + xderBcomp[3]);               _jstore[1][10] = -(yderBcomp[0] + yderBcomp[3]);
    _jstore[0][11] = -(xderBcomp[1] - B2sign[1]*I*xderBcomp[2]);   _jstore[1][11] = -(yderBcomp[1] - B2sign[1]*I*yderBcomp[2]);
    _jstore[0][14] = -(xderBcomp[1] + B2sign[1]*I*xderBcomp[2]);   _jstore[1][14] = -(yderBcomp[1] + B2sign[1]*I*yderBcomp[2]);
    _jstore[0][15] = -(xderBcomp[0] - xderBcomp[3]);               _jstore[1][15] = -(yderBcomp[0] - yderBcomp[3]);
    

}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > Haldane2L::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

void Haldane2L::SetBasis()
{
    vec_a[0][0] = 0.;                   vec_a[0][1] = a0;
    vec_a[1][0] = -sqrt(3.)/2*a0;       vec_a[1][1] = -0.5*a0;
    vec_a[2][0] = sqrt(3.)/2*a0;        vec_a[2][1] = -0.5*a0;

    vec_b[0][0] = sqrt(3.)*a0;          vec_b[0][1] = 0.;
    vec_b[1][0] = -sqrt(3.)/2*a0;       vec_b[1][1] = 3./2.*a0;
    vec_b[2][0] = -sqrt(3.)/2*a0;       vec_b[2][1] = -3./2.*a0;
}

void Haldane2L::PrintMaterialInformation()
{
    std::cout << "============================================\n";
    std::cout << "Haldane2L Model paramters\n";
    std::cout << "a0         = " << a0 << " a.u. \n";
    std::cout << "t1         = " << t1[0] << ", " << t1[1] << " a.u. \n";
    std::cout << "t2         = " << t2[0] << ", " << t2[1] << " a.u. \n";
    std::cout << "M0         = " << M0[0] << ", " << M0[1] << " a.u. \n";
    std::cout << "B2 sign    = " << B2sign[0] << ", " << B2sign[1] << " \n";
    std::cout << "phi0       = " << phi0[0] << ", " << phi0[1] << " rad \n";
    std::cout << "t11        = " << t11 << " a.u. \n";
    std::cout << "polytype   = " << polytype << "\n";
    std::cout << "============================================\n";
}

void Haldane2L::GenBcomp(std::array<double, Ndim> _kpoint)
{
    std::fill(Bcomp, Bcomp+4, 0.);
    for (int i = 0; i < 3; ++i)
    {
        angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
        angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
        Bcomp[0] += cos( angle_b0 );    Bcomp[1] += cos( angle_a0 );
        Bcomp[2] += sin( angle_a0 );    Bcomp[3] += sin( angle_b0 );
    }

    Bcomp[0] *= 2.* t2[0] * cos( phi0[0] );       Bcomp[1] *= t1[0];
    Bcomp[2] *= t1[0];                            Bcomp[3] = M0[0] - 2.*t2[0]*sin( phi0[0] )*Bcomp[3];

    if (is_top_bottom_same)
    {
        Bcomp_bottom[0] = Bcomp[0];    Bcomp_bottom[1] = Bcomp[1];
        Bcomp_bottom[2] = Bcomp[2];    Bcomp_bottom[3] = Bcomp[3];
    }
    else
    {
        std::fill(Bcomp_bottom, Bcomp_bottom+4, 0.);
        for (int i = 0; i < 3; ++i)
        {
            angle_a0 = _kpoint[0]*vec_a[i][0] + _kpoint[1]*vec_a[i][1];
            angle_b0 = _kpoint[0]*vec_b[i][0] + _kpoint[1]*vec_b[i][1];
            Bcomp_bottom[0] += cos( angle_b0 );    Bcomp_bottom[1] += cos( angle_a0 );
            Bcomp_bottom[2] += sin( angle_a0 );    Bcomp_bottom[3] += sin( angle_b0 );
        }

        Bcomp_bottom[0] *= 2.* t2[1] * cos( phi0[1] );       Bcomp_bottom[1] *= t1[1];
        Bcomp_bottom[2] *= t1[1];                            Bcomp_bottom[3] = M0[1] - 2.*t2[1]*sin( phi0[1] )*Bcomp[3];
    }
}
