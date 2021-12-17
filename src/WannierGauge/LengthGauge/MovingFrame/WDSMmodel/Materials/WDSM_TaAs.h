/// Hamiltonian and dipole generator for  WDSM_TaAs model
/**
 * 
*/
#pragma once

#include "../BaseMaterial.h"
#include "../utility.h"

#include <string>
#include <iostream>
#include <fstream>

class WDSM_TaAs : public BaseMaterial
{
public:
    // Haldnae parameters
    double t, delta,my,mz, a0;
      ///< 3 NNN vectors (from A to A, or B to B)

    double eps,angle_a0, angle_b0;

    double Bcomp[16];
    double Bcomp_c[16];

    std::array<double, 9> BZaxis;
    std::array<double, 3> BZorigin;

    WDSM_TaAs( const libconfig::Setting *params );
    // void SetBasis();
    void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) override;
    void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) override {};    // Do nothing!
    void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) override;
    std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone(  ) override;

    void PrintMaterialInformation() override;

    void GenBcomp(std::array<double, Ndim> _kpoint);
};

WDSM_TaAs::WDSM_TaAs( const libconfig::Setting *params )
{
    Nband = 4;  Nval = 1;
    if (params->lookupValue("my", my)
        && params->lookupValue("mz", mz)
      && params->lookupValue("t", t)
      && params->lookupValue("delta", delta)
    && params->lookupValue("a0", a0) )
    {
        a0 /= au_angstrom;

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
                  0,          0,         2*kzMax};
    BZorigin = {0, 0, 0};

    eps = 1.0e-18;
}

void WDSM_TaAs::GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint)
{
    GenBcomp(_kpoint);
    _hstore[0]=Bcomp[0]+I*Bcomp_c[0]; _hstore[1]=Bcomp[1]+I*Bcomp_c[1]; _hstore[2]=Bcomp[2]+I*Bcomp_c[2]; _hstore[3]=Bcomp[3]+I*Bcomp_c[3];
    _hstore[4]=Bcomp[4]+I*Bcomp_c[4]; _hstore[5]=Bcomp[5]+I*Bcomp_c[5]; _hstore[6]=Bcomp[6]+I*Bcomp_c[6]; _hstore[7]=Bcomp[7]+I*Bcomp_c[7];
    _hstore[8]=Bcomp[8]+I*Bcomp_c[8]; _hstore[9]=Bcomp[9]+I*Bcomp_c[9]; _hstore[10]=Bcomp[10]+I*Bcomp_c[10]; _hstore[11]=Bcomp[11]+I*Bcomp_c[11];
    _hstore[12]=Bcomp[12]+I*Bcomp_c[12]; _hstore[13]=Bcomp[13]+I*Bcomp_c[13]; _hstore[14]=Bcomp[14]+I*Bcomp_c[14];_hstore[15]=Bcomp[15]+I*Bcomp_c[15];
    // for (int i=0;i<16;i++){
    //   _hstore[i]=Bcomp[i]+I*Bcomp_c[i];
    // }
}

void  WDSM_TaAs::GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint)
{
    double xderBcomp[16];  double xderBcomp_c[16];  double yderBcomp[16]; double yderBcomp_c[16]; double zderBcomp[16]; double zderBcomp_c[16];
    // cout<<"here wdsm\n";
    // fill(Bcomp, Bcomp+16, 0.);
    fill(xderBcomp, xderBcomp+16, 0.);
    fill(yderBcomp, yderBcomp+16, 0.);
    fill(zderBcomp, zderBcomp+16, 0.);

    fill(xderBcomp_c, xderBcomp_c+16, 0.);
    fill(yderBcomp_c, yderBcomp_c+16, 0.);
    fill(zderBcomp_c, yderBcomp_c+16, 0.);

    ////// version 0
    xderBcomp[2] += -t*sin( _kpoint[0]*a0);
    xderBcomp[7] += -t*sin( _kpoint[0]*a0);
    xderBcomp[8] += -t*sin( _kpoint[0]*a0);
    xderBcomp[13] +=-t*sin( _kpoint[0]*a0);


    yderBcomp[2] += t*my*sin( _kpoint[1]*a0);
    yderBcomp[7] += t*my*sin( _kpoint[1]*a0);
    yderBcomp[8] += t*my*sin( _kpoint[1]*a0);
    yderBcomp[13] += t*my*sin( _kpoint[1]*a0);

    yderBcomp_c[2] += -t*cos( _kpoint[1]*a0);
    yderBcomp_c[7] += -t*cos( _kpoint[1]*a0);
    yderBcomp_c[8] += t*cos( _kpoint[1]*a0);
    yderBcomp_c[13] += t*cos( _kpoint[1]*a0);

    yderBcomp_c[3] += delta*t*sin( _kpoint[1]*a0);
    yderBcomp_c[6] += delta*t*sin( _kpoint[1]*a0);
    yderBcomp_c[9] += -delta*t*sin( _kpoint[1]*a0);
    yderBcomp_c[12] += -delta*t*sin( _kpoint[1]*a0);


    zderBcomp[2] += t*mz*sin( _kpoint[2]*a0);
    zderBcomp[7] += t*mz*sin( _kpoint[2]*a0);
    zderBcomp[8] += t*mz*sin( _kpoint[2]*a0);
    zderBcomp[13] += t*mz*sin( _kpoint[2]*a0);

    zderBcomp[1] += t*cos( _kpoint[2]*a0);
    zderBcomp[4] += t*cos( _kpoint[2]*a0);
    zderBcomp[11] += -t*cos( _kpoint[2]*a0);
    zderBcomp[14] += -t*cos( _kpoint[2]*a0);
    // cout<<"here wdsm\n";

    _jstore[0][0]=-xderBcomp[0];
    _jstore[0][1]=-xderBcomp[1];
    _jstore[0][2]=-xderBcomp[2];
    _jstore[0][3]=-xderBcomp[3];
    _jstore[0][4]=-xderBcomp[4];
    _jstore[0][5]=-xderBcomp[5];
    _jstore[0][6]=-xderBcomp[6];
    _jstore[0][7]=-xderBcomp[7];
    _jstore[0][8]=-xderBcomp[8];
    _jstore[0][9]=-xderBcomp[9];
    _jstore[0][10]=-xderBcomp[10];
    _jstore[0][11]=-xderBcomp[11];
    _jstore[0][12]=-xderBcomp[12];
    _jstore[0][13]=-xderBcomp[13];
    _jstore[0][14]=-xderBcomp[14];
    _jstore[0][15]=-xderBcomp[15];

    _jstore[1][0]=-yderBcomp[0]-I*yderBcomp_c[0];
    _jstore[1][1]=-yderBcomp[1]-I*yderBcomp_c[1];
    _jstore[1][2]=-yderBcomp[2]-I*yderBcomp_c[2];
    _jstore[1][3]=-yderBcomp[3]-I*yderBcomp_c[3];
    _jstore[1][4]=-yderBcomp[4]-I*yderBcomp_c[4];
    _jstore[1][5]=-yderBcomp[5]-I*yderBcomp_c[5];
    _jstore[1][6]=-yderBcomp[6]-I*yderBcomp_c[6];
    _jstore[1][7]=-yderBcomp[7]-I*yderBcomp_c[7];
    _jstore[1][8]=-yderBcomp[8]-I*yderBcomp_c[8];
    _jstore[1][9]=-yderBcomp[9]-I*yderBcomp_c[9];
    _jstore[1][10]=-yderBcomp[10]-I*yderBcomp_c[10];
    _jstore[1][11]=-yderBcomp[11]-I*yderBcomp_c[11];
    _jstore[1][12]=-yderBcomp[12]-I*yderBcomp_c[12];
    _jstore[1][13]=-yderBcomp[13]-I*yderBcomp_c[13];
    _jstore[1][14]=-yderBcomp[14]-I*yderBcomp_c[14];
    _jstore[1][15]=-yderBcomp[15]-I*yderBcomp_c[15];

    _jstore[1][0]=-zderBcomp[0];
    _jstore[1][1]=-zderBcomp[1];
    _jstore[1][2]=-zderBcomp[2];
    _jstore[1][3]=-zderBcomp[3];
    _jstore[1][4]=-zderBcomp[4];
    _jstore[1][5]=-zderBcomp[5];
    _jstore[1][6]=-zderBcomp[6];
    _jstore[1][7]=-zderBcomp[7];
    _jstore[1][8]=-zderBcomp[8];
    _jstore[1][9]=-zderBcomp[9];
    _jstore[1][10]=-zderBcomp[10];
    _jstore[1][11]=-zderBcomp[11];
    _jstore[1][12]=-zderBcomp[12];
    _jstore[1][13]=-zderBcomp[13];
    _jstore[1][14]=-zderBcomp[14];
    _jstore[1][15]=-zderBcomp[15];


    // for (int itmp = 0; itmp<16; itmp++){
    //     // cout<<-xderBcomp[i]-I*xderBcomp_c[i]<<endl;
    //     _jstore[0][itmp]=-xderBcomp[itmp]-I*xderBcomp_c[itmp];
    //     _jstore[1][itmp]=-yderBcomp[itmp]-I*yderBcomp_c[itmp];
    //     _jstore[2][itmp]=-zderBcomp[itmp]-I*zderBcomp_c[itmp];
    // }
    // cout<<"here wdsm 1\n";


}

std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> >  WDSM_TaAs::GenBrillouinzone(  )
{
    return std::make_tuple(BZaxis, BZorigin );
}

// void  WDSM_TaAs::SetBasis()
// {
//     vec_a[0][0] = 1.;                   vec_a[0][1] = a0;
//     vec_a[1][0] = 0.;                   vec_a[1][1] = -0.5*a0;example
//     vec_a[2][0] = 0. ;                   vec_a[2][1] = -0.5*a0;
//
//     vec_b[0][0] = 0.0;                   vec_b[0][1] = 0.;
//     vec_b[1][0] = 1.0;       vec_b[1][1] = 3./2.*a0;
//     vec_b[2][0] = 0.0;       vec_b[2][1] = -3./2.*a0;
// }

void  WDSM_TaAs::PrintMaterialInformation()
{
    cout << "============================================\n";
    cout << " WDSM_TaAs Model paramters\n";
    cout << "a0         = " << a0 << " a.u. \n";
    cout << "my         = " << my << " a.u. \n";
    cout << "mz         = " << mz << " a.u. \n";
    cout << "t          = " << t << " a.u. \n";
    cout << "delta         = " << delta << " a.u. \n";
    cout << "============================================\n";
}

void  WDSM_TaAs::GenBcomp(std::array<double, Ndim> _kpoint)
{

    fill(Bcomp, Bcomp+16, 0.);
    fill(Bcomp_c, Bcomp_c+16, 0.);


    Bcomp[1]+=t*sin(_kpoint[2]*a0);
    Bcomp[4]+=t*sin(_kpoint[2]*a0);

    Bcomp[11]+=-t*sin(_kpoint[2]*a0);
    Bcomp[14]+=-t*sin(_kpoint[2]*a0);




    Bcomp[2]+=t*cos(_kpoint[0]*a0)+t*my*(1-cos(_kpoint[1]*a0))+ t*mz*(1- cos(_kpoint[2]*a0));
    Bcomp[7]+= t*cos(_kpoint[0]*a0)+t*my*(1-cos(_kpoint[1]*a0))+ t*mz*(1- cos(_kpoint[2]*a0));

    Bcomp_c[2]-= t*sin(_kpoint[1]*a0);
    Bcomp_c[3]-= t*cos(_kpoint[1]*a0)*delta;
    Bcomp_c[6]-= t*cos(_kpoint[1]*a0)*delta;
    Bcomp_c[7]-= t*sin(_kpoint[1]*a0);



    Bcomp[8]+= t*cos(_kpoint[0]*a0)+t*my*(1- cos(_kpoint[1]*a0))+t*mz*(1- cos(_kpoint[2]*a0));
    Bcomp[13]+= t*cos(_kpoint[0]*a0)+t*my*(1- cos(_kpoint[1]*a0))+t*mz*(1- cos(_kpoint[2]*a0));

    Bcomp_c[8]+= t*sin(_kpoint[1]*a0);
    Bcomp_c[9]+= t*cos(_kpoint[1]*a0)*delta;
    Bcomp_c[12]+= t*cos(_kpoint[1]*a0)*delta;
    Bcomp_c[13]+= t*sin(_kpoint[1]*a0);
    // cout<<Bcomp_c[9]<<endl;

}
