/// Material property generator for TMDCs - 3 band approximation
/**
 * @brief Material property generator for TMDCs - 3 band approximation. 
 * @details Equation for this class comes from Phys. Rev. B 88, 085433. Spinless - 3band, soc - 6 band.
 * @author Dasol Kim
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>

class TMDCs : public BaseMaterial
{
public:
    // TB parameters
    double e1, e2, t0, t1, t2, t11, t12, t22, r0, r2, r11, r12, u0, u1, u2, u11, u12, u22, r1;
    double a0, socparam;
    double vec_a[3][2];          ///< 3  NN vectors (from A to B )
    double vec_b[3][2];          ///< 3 NNN vectors (from A to A, or B to B)

    double eps;                  ///< to avoid singular norm

    bool isSOC;
    std::string material;

    double angle_a0[3], angle_b0[3];

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    complex *tempUmat;
    complex *tempctransUmat;
    complex *tempHmat;
    double *tempEval;
    int *tempIsuppz;

    TMDCs( const libconfig::Setting *params );
    ~TMDCs();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;
};

TMDCs::TMDCs( const libconfig::Setting *params )
{
    isSOC = false;
    params->lookupValue("isSOCEnabled", isSOC);

    if (isSOC)
    {
        Nband = 6;  Nval = 2;
    }
    else
    {
        Nband = 3;  Nval = 1;
    }
    isDipoleZero = true;
    
    //string tempmaterial = params->lookup("Material");
    //material = tempmaterial;
    params->lookupValue("Material", material);
    // hard-coded paramter - GGA
    double gga_params[6][19] = {0.683, 1.707, -0.146, -0.114, 0.506, 0.085, 0.162, 0.073, 0.060, -0.236,
                                0.067, 0.016, 0.087, -0.038, 0.046, 0.001, 0.266, -0.176, -0.150,
                                0.717, 1.916, -0.152, -0.097, 0.590, 0.047, 0.178, 0.016, 0.069, -0.261,
                                0.107, -0.003, 0.109, -0.054, 0.045, 0.002, 0.325, -0.206, -0.163,
                                0.684, 1.546, -0.146, -0.130, 0.432, 0.144, 0.117, 0.075, 0.039, -0.209,
                                0.069, 0.052, 0.060, -0.042, 0.036, 0.008, 0.272, -0.172, -0.150,
                                0.728, 1.655, -0.146, -0.124, 0.507, 0.117, 0.127, 0.015, 0.036, -0.234,
                                0.107, 0.044, 0.075, -0.061, 0.032, 0.007, 0.329, -0.202, -0.164,
                                0.588, 1.303, -0.226, -0.234, 0.036, 0.400, 0.098, 0.017, 0.003, -0.025,
                                -0.169, 0.082, 0.051, 0.057, 0.103, 0.187, -0.045, -0.141, 0.087,
                                0.697, 1.380, -0.109, -0.164, 0.368, 0.204, 0.093, 0.038, -0.015, -0.209,
                                0.107, 0.115, 0.009, -0.066, 0.011, -0.013, 0.312, -0.177, -0.132};
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 19; ++j)
        {
            gga_params[i][j] /= au_eV;
        }
    }
    int materialIndex = -1;
    // MoS2
    if (material == "MoS2")
    {
        a0 = 3.190/au_angstrom; materialIndex = 0;  socparam = 0.073/au_eV;
    }
    else if (material == "WS2")
    {
        a0 = 3.191/au_angstrom; materialIndex = 1;  socparam = 0.211/au_eV;
    }
    else if (material == "MoSe2")
    {
        a0 = 3.326/au_angstrom; materialIndex = 2;  socparam = 0.091/au_eV;
    }
    else if (material == "WSe2")
    {
        a0 = 3.325/au_angstrom; materialIndex = 3;  socparam = 0.228/au_eV;
    }
    else if (material == "MoTe2")
    {
        a0 = 3.557/au_angstrom; materialIndex = 4;  socparam = 0.107/au_eV;
    }
    else if (material == "WTe2")
    {
        a0 = 3.560/au_angstrom; materialIndex = 5;  socparam = 0.237/au_eV;
    }
    else
    {
        std::cout << "Unknowm material: " << material << std::endl;
        std::exit(EXIT_FAILURE);
    }
    

    // assign paramters
    e1 = gga_params[materialIndex][0];
    e2 = gga_params[materialIndex][1];
    t0 = gga_params[materialIndex][2];
    t1 = gga_params[materialIndex][3];
    t2 = gga_params[materialIndex][4];
    t11 = gga_params[materialIndex][5];
    t12 = gga_params[materialIndex][6];
    t22 = gga_params[materialIndex][7];
    r0 = gga_params[materialIndex][8];
    r1 = gga_params[materialIndex][9];
    r2 = gga_params[materialIndex][10];
    r11 = gga_params[materialIndex][11];
    r12 = gga_params[materialIndex][12];
    u0 = gga_params[materialIndex][13];
    u1 = gga_params[materialIndex][14];
    u2 = gga_params[materialIndex][15];
    u11 = gga_params[materialIndex][16];
    u12 = gga_params[materialIndex][17];
    u22 = gga_params[materialIndex][18];

    //SetBasis();

    tempUmat = new complex[Nband*Nband];
    tempctransUmat = new complex[Nband*Nband];
    tempHmat = new complex[Nband*Nband];
    tempEval = new double[Nband];
    tempIsuppz = new int[2*Nband];

    // lattice vector: (a0, 0), (-a0/2, sqrt(3.)*a0/2)
    // --> reciprocal : (2*pi/a0, 2*pi/a0/sqrt(3.) ), (0., 4*pi/a0/sqrt(3.))
    // same with Haldane notation (sqrt(3.)*a0 --> a0)
    //double kxMax = pi/sqrt(3.)/a0;  double kyMax = 2.*pi/3./a0;
    double kxMax = pi/a0;  double kyMax = 2.*pi/sqrt(3.)/a0;
    BZaxis = { 2*kxMax,       0,         0,
                  0,       2*kyMax,      0,
                  0,          0,         0};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

TMDCs::~TMDCs()
{
    delete[] tempUmat, tempHmat, tempEval, tempIsuppz, tempctransUmat;
}

void TMDCs::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    double r3 = sqrt(3.);
    double alpha = _kpoint[0]*a0/2; double beta = _kpoint[1]*a0*r3/2;

    double ca[4];   double sa[4];   double cb[4];   double sb[4];
    for (int i = 0; i < 4; ++i)
    {
        ca[i] = cos((i+1)*alpha);  cb[i] = cos((i+1)*beta);
        sa[i] = sin((i+1)*alpha);  sb[i] = sin((i+1)*beta);
    }

    double V0 = e1 + 2*t0*(2*ca[0]*cb[0] + ca[1]) + 2*r0*(2*ca[2]*cb[0] + cb[1]) + 2*u0*(2*ca[1]*cb[1] + ca[3]);

    complex V1 = (-2*r3*t2*sa[0]*sb[0] + 2*(r1+r2)*sa[2]*sb[0] -2*r3*u2*sa[1]*sb[1]) 
                + I*(2*t1*sa[0]*(2*ca[0]+cb[0]) + 2*(r1-r2)*sa[2]*cb[0] + 2*u1*sa[1]*(2*ca[1]+cb[1]));

    complex V2 = (2*t2*(ca[1] - ca[0]*cb[0]) - 2/r3*(r1+r2)*(ca[2]*cb[0] - cb[1]) + 2*u2*(ca[3] - ca[1]*cb[1])) 
                + I*(2*r3*t1*ca[0]*sb[0] + 2/r3*sb[0]*(r1-r2)*(ca[2]+2*cb[0]) + 2*r3*u1*ca[1]*sb[1]);

    double V11 = e2 + (t11+3*t22)*ca[0]*cb[0] + 2*t11*ca[1] + 4*r11*ca[2]*cb[0] + 2*(r11+r3*r12)*cb[1] + (u11+3*u22)*ca[1]*cb[1] + 2*u11*ca[3];

    complex V12 = (r3*(t22-t11)*sa[0]*sb[0] + 4*r12*sa[2]*sb[0] + r3*(u22-u11)*sa[1]*sb[1]) 
                + I*(4*t12*sa[0]*(ca[0]-cb[0]) + 4*u12*sa[1]*(ca[1]-cb[1]));

    double V22 = e2 + (3*t11+t22)*ca[0]*cb[0] + 2*t22*ca[1] + 2*r11*(2*ca[2]*cb[0]+cb[1]) 
                + 2/r3*r12*(4*ca[2]*cb[0]-cb[1]) + (3*u11+u22)*ca[1]*cb[1] + 2*u22*ca[3];

    if (isSOC)
    {
        std::fill(_hstore, _hstore + Nband*Nband, 0.);
        _hstore[0] = V0;                _hstore[1] = V1;                        _hstore[2] = V2;
        _hstore[6] = conj(V1);          _hstore[7] = V11;                       _hstore[8] = V12 + I*socparam;
        _hstore[12] = conj(V2);         _hstore[13] = conj(V12) - I*socparam;   _hstore[14] = V22;

        _hstore[21] = V0;                _hstore[22] = V1;                        _hstore[23] = V2;
        _hstore[27] = conj(V1);          _hstore[28] = V11;                       _hstore[29] = V12 - I*socparam;
        _hstore[33] = conj(V2);          _hstore[34] = conj(V12) + I*socparam;    _hstore[35] = V22;
    }
    else
    {
        _hstore[0] = V0;        _hstore[1] = V1;        _hstore[2] = V2;
        _hstore[3] = conj(V1);  _hstore[4] = V11;       _hstore[5] = V12;
        _hstore[6] = conj(V2);  _hstore[7] = conj(V12); _hstore[8] = V22;
    }
}

void TMDCs::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    // diff components

    double r3 = sqrt(3.);
    double alpha = _kpoint[0]*a0/2; double beta = _kpoint[1]*a0*r3/2;

    double ca[4];   double sa[4];   double cb[4];   double sb[4];
    for (int i = 0; i < 4; ++i)
    {
        ca[i] = cos((i+1)*alpha);  cb[i] = cos((i+1)*beta);
        sa[i] = sin((i+1)*alpha);  sb[i] = sin((i+1)*beta);
    }

    double dca[4];   double dsa[4];   double dcb[4];   double dsb[4];
    for (int i = 0; i < 4; ++i)
    {
        dca[i] = -(i+1)*a0/2*sin((i+1)*alpha);  dcb[i] = -(i+1)*a0*r3/2*sin((i+1)*beta);
        dsa[i] = (i+1)*a0/2*cos((i+1)*alpha);   dsb[i] = (i+1)*a0*r3/2*cos((i+1)*beta);
    }

    double dxV0 = 2*t0*(2*dca[0]*cb[0] + dca[1]) + 2*r0*(2*dca[2]*cb[0]) + 2*u0*(2*dca[1]*cb[1] + dca[3]);

    complex dxV1 = (-2*r3*t2*dsa[0]*sb[0] + 2*(r1+r2)*dsa[2]*sb[0] -2*r3*u2*dsa[1]*sb[1]) 
                + I*(2*t1*dsa[0]*(2*ca[0]+cb[0]) + 2*t1*sa[0]*(2*dca[0]) + 2*(r1-r2)*dsa[2]*cb[0] + 2*u1*dsa[1]*(2*ca[1]+cb[1]) + 2*u1*sa[1]*(2*dca[1]));

    complex dxV2 = (2*t2*(dca[1] - dca[0]*cb[0]) - 2/r3*(r1+r2)*(dca[2]*cb[0]) + 2*u2*(dca[3] - dca[1]*cb[1])) 
                + I*(2*r3*t1*dca[0]*sb[0] + 2/r3*sb[0]*(r1-r2)*(dca[2]) + 2*r3*u1*dca[1]*sb[1]);

    double dxV11 = (t11+3*t22)*dca[0]*cb[0] + 2*t11*dca[1] + 4*r11*dca[2]*cb[0] + (u11+3*u22)*dca[1]*cb[1] + 2*u11*dca[3];

    complex dxV12 = (r3*(t22-t11)*dsa[0]*sb[0] + 4*r12*dsa[2]*sb[0] + r3*(u22-u11)*dsa[1]*sb[1]) 
                + I*(4*t12*dsa[0]*(ca[0]-cb[0]) + 4*t12*sa[0]*(dca[0]) + 4*u12*dsa[1]*(ca[1]-cb[1]) + 4*u12*sa[1]*(dca[1]));

    double dxV22 = (3*t11+t22)*dca[0]*cb[0] + 2*t22*dca[1] + 2*r11*(2*dca[2]*cb[0]) 
                + 2/r3*r12*(4*dca[2]*cb[0]) + (3*u11+u22)*dca[1]*cb[1] + 2*u22*dca[3];
    

    double dyV0 = 2*t0*(2*ca[0]*dcb[0]) + 2*r0*(2*ca[2]*dcb[0] + dcb[1]) + 2*u0*(2*ca[1]*dcb[1]);

    complex dyV1 = (-2*r3*t2*sa[0]*dsb[0] + 2*(r1+r2)*sa[2]*dsb[0] -2*r3*u2*sa[1]*dsb[1]) 
                + I*(2*t1*sa[0]*(dcb[0]) + 2*(r1-r2)*sa[2]*dcb[0] + 2*u1*sa[1]*(dcb[1]));

    complex dyV2 = (2*t2*(-ca[0]*dcb[0]) - 2/r3*(r1+r2)*(ca[2]*dcb[0] - dcb[1]) + 2*u2*(-ca[1]*dcb[1])) 
                + I*(2*r3*t1*ca[0]*dsb[0] + 2/r3*dsb[0]*(r1-r2)*(ca[2]+2*cb[0]) + 2/r3*sb[0]*(r1-r2)*(2*dcb[0]) + 2*r3*u1*ca[1]*dsb[1]);

    double dyV11 = (t11+3*t22)*ca[0]*dcb[0] + 4*r11*ca[2]*dcb[0] + 2*(r11+r3*r12)*dcb[1] + (u11+3*u22)*ca[1]*dcb[1];

    complex dyV12 = (r3*(t22-t11)*sa[0]*dsb[0] + 4*r12*sa[2]*dsb[0] + r3*(u22-u11)*sa[1]*dsb[1]) 
                + I*(4*t12*sa[0]*(-dcb[0]) + 4*u12*sa[1]*(-dcb[1]));

    double dyV22 = (3*t11+t22)*ca[0]*dcb[0] + 2*r11*(2*ca[2]*dcb[0]+dcb[1]) 
                + 2/r3*r12*(4*ca[2]*dcb[0]-dcb[1]) + (3*u11+u22)*ca[1]*dcb[1];

    
    if (isSOC)
    {
        _jstore[0][0] = dxV0;            _jstore[0][1] = dxV1;            _jstore[0][2] = dxV2;
        _jstore[0][6] = conj(dxV1);      _jstore[0][7] = dxV11;           _jstore[0][8] = dxV12;
        _jstore[0][12] = conj(dxV2);     _jstore[0][13] = conj(dxV12);    _jstore[0][14] = dxV22;

        _jstore[1][0] = dyV0;            _jstore[1][1] = dyV1;            _jstore[1][2] = dyV2;
        _jstore[1][6] = conj(dyV1);      _jstore[1][7] = dyV11;           _jstore[1][8] = dyV12;
        _jstore[1][12] = conj(dyV2);     _jstore[1][13] = conj(dyV12);    _jstore[1][14] = dyV22;

        _jstore[0][21] = dxV0;            _jstore[0][22] = dxV1;            _jstore[0][23] = dxV2;
        _jstore[0][27] = conj(dxV1);      _jstore[0][28] = dxV11;           _jstore[0][29] = dxV12;
        _jstore[0][33] = conj(dxV2);      _jstore[0][34] = conj(dxV12);     _jstore[0][35] = dxV22;

        _jstore[1][21] = dyV0;            _jstore[1][22] = dyV1;            _jstore[1][23] = dyV2;
        _jstore[1][27] = conj(dyV1);      _jstore[1][28] = dyV11;           _jstore[1][29] = dyV12;
        _jstore[1][33] = conj(dyV2);      _jstore[1][34] = conj(dyV12);     _jstore[1][35] = dyV22;
    }
    else
    {
        _jstore[0][0] = dxV0;           _jstore[0][1] = dxV1;           _jstore[0][2] = dxV2;
        _jstore[0][3] = conj(dxV1);     _jstore[0][4] = dxV11;          _jstore[0][5] = dxV12;
        _jstore[0][6] = conj(dxV2);     _jstore[0][7] = conj(dxV12);    _jstore[0][8] = dxV22;

        _jstore[1][0] = dyV0;           _jstore[1][1] = dyV1;           _jstore[1][2] = dyV2;
        _jstore[1][3] = conj(dyV1);     _jstore[1][4] = dyV11;          _jstore[1][5] = dyV12;
        _jstore[1][6] = conj(dyV2);     _jstore[1][7] = conj(dyV12);    _jstore[1][8] = dyV22;
    }
}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > TMDCs::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

// void BieSe3surf::SetBasis()
// {
//     // z-components are ignored
//     vec_a[0][0] = a0;                   vec_a[0][1] = 0.;
//     vec_a[1][0] = -a0/2;                vec_a[1][1] = sqrt(3.)/2*a0;
//     vec_a[2][0] = -a0/2;                vec_a[2][1] = -sqrt(3.)/2*a0;

//     vec_b[0][0] = 0.;                   vec_b[0][1] = sqrt(3.)*a0/3;
//     vec_b[1][0] = -a0/2;                vec_b[1][1] = -sqrt(3.)/6*a0;
//     vec_b[2][0] = a0/2;                 vec_b[2][1] = -sqrt(3.)/6*a0;
// }

void TMDCs::PrintMaterialInformation()
{
    std::cout << "============================================\n";
    std::cout << "TMDC - target: " << material << std::endl;
    std::cout << "============================================\n";
}