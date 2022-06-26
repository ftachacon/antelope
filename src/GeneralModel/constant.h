/// collection of constants used in projects
/**
 * @author Alexis Agustín  Chacón Salazar
 * @author Dasol Kim
*/

#ifndef CONSTANT_H
#define CONSTANT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <new>
#include <vector>
#include <string>
#include <complex>
//#define MKL_Complex16 std::complex<double>
#define complex std::complex<double>

using namespace std;


const complex I=complex(0.,1.);	                //Imaginary basic number i=sqrt(-1)
const double lightC_au = 137.036;				//Speed of the light in atomic units 
const double one_by_lightC_au = 1./lightC_au;	//Inverse of the light speed

const double pi		= 4.0*atan(1.0);				//Pi definition
const double dospi	=	2.*pi;                  // 6.2831853071795862;
const double charge_electron_au = -1.;			//Electronic charge in atomic unit
const double mass_proton_au=1836.;				//Proton mass in atomic units 
const double BohrRadius = 5.29177211e-11;		//Bohr's radius in a.u.
const double time_SI = 2.418884326505e-17;		//Equivalencia de tiempo en segundos de una unidad atómica
const double lightC_SI = 2.99792458e8;			//Light speed international sistem units

const double au_angstrom = 0.529177210903;      ///< 1 a.u. in \f$\AA\f$
const double au_eV =  27.211386245988;          ///< 1 a.u. in eV

constexpr int Ndim      = 3;
const int Ngrad     = 2;

enum class GaugeType
{
    LengthWannier,
    LengthHamiltonian,
    VelocityHamiltonian
};


#endif  /* CONSTANT_H */


