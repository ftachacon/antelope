/**
 * @file laser.h
 * Laser-related things are implemented in this head file.
 */

#pragma once

#include <stdlib.h>
#include <string>
#include "constant.h"
#include <vector>
#include <array>
#include <algorithm>
using namespace std;

const double gaussian_width_multiply = 2.0;
const double gaussian_factor = gaussian_width_multiply*gaussian_width_multiply* 4*log(2);

enum Envelope
{
    ENVELOPE_GAUSSIAN,
    ENVELOPE_COS4
};

/// representation of single pulse
class Pulse
{
public:
    double E0;                              ///< Electric Field peak strength (a.u.)
    double w0;			                    ///< Centrl frequency of pulse (a.u.)
    double ellip;   	                    ///< Ellipticity of pulse (-1.0 <= ellip <= 1.0). ellip > 0 is right-handed polarization.
	double ncycles;		                    ///< If envelope is gaussian, number of cycles per FWHM. If envelope is cos-n squares, number of cycles per pulse duration.
    double cep;		                        ///< Carrier envelope phase of pulse (rad)
    double t0;                              ///< Mid-time of pulse (a.u.)
    double thetaz;                          ///< Polar angle of propagation direction (rad)
    double phiz;                            ///< Azimuthal angle of propagation direction (rad)
    double phix;                            ///< Azimuthal angle of semi-major axis in transformed coordinate (rad)
	string envname;	                        ///< Name of envelope type used in the pulse ("cos4", "gauss")

// below variables are only for internal data processing or printing
    double I0;			                    // Intensity per pulse (W/cm^2)
    double period0;		                    //Period per pulse (a.u.)
    double E0x, E0y;	                    //Electric Field Component in the x, y-direction
	double A0x, A0y;	                    //Vector potential component in the x, y-direction
    double twidth;		                    // Time bandwidth per pulse (-twidth <= <= twidth is active region)

    Envelope envtype;               

    double envfactor;                       // This variable is denfied to reduce the repeated calculation in f_envelope. ex) g(envfactor*t)
    double inverse_w0;                      // = 1/w0. Remember multiplication is faster than division

    array<double, Ndim> xdir, ydir;         // direction of laser x-axis and y-axis in real sapce
public:
    /// Constructor for pulse, 
    /**
     * Propagation direction is dentermined by two paramter thetaz and phiz. 
     * They are just polar and azimuthal angle in spherical coordinatte, so propagation direction is
     * \f$ (\sin\theta_{z}\cos\phi_{z}, \sin\theta_{z}\sin\phi_{z}, \cos\theta_{z})\f$.
     * For the definition of LCP, and RCP I choose view from source. Therefore, for left or right hand rule,
     * direction of thumb is same as propagtion direction (out of source).
     * To doing standard 2d calculation (x-y plane), set thetaz = 0 and phiz = 0, or just omit them from the constructor. 
     */
    Pulse(double _E0, double _w0, double _ellip, double _ncycle, double _cep,  double _t0, double _phix, string _env_name, double _thetaz = 0.0, double _phiz = 0.0 )
        : E0(_E0), w0(_w0), ellip(_ellip), ncycles(_ncycle), cep(_cep), t0(_t0), phix(_phix), envname(_env_name), thetaz(_thetaz), phiz(_phiz)
    {
        I0 = E0*E0 * 3.5e16;

        period0    =  dospi/w0; // laser period or cycle
        
        if (envname == "cos4")  envtype = ENVELOPE_COS4;
        else if (envname == "gauss") envtype = ENVELOPE_GAUSSIAN;
        else
        {
            cerr << "invalid envelope type\n";
            exit(EXIT_FAILURE);
        }
        
        // laser time duration
        switch (envtype)
        {
        case ENVELOPE_GAUSSIAN:
            twidth = gaussian_width_multiply * ncycles*period0;
            envfactor = -gaussian_factor/(twidth*twidth);
            break;
        
        // same for every cos n-square envelope
        default:
            twidth = ncycles*period0/2;
            envfactor = w0/(2*ncycles);
            break;
        }
        

        double factor = sqrt( 1.0 / (1.0 + ellip*ellip) );
        E0x = factor * E0;
        E0y = ellip * factor * E0;

        if (w0 > 1.0e-30)
        {
            A0x = E0x / w0; A0y = E0y / w0;
        }
        else
        {
            A0x = 0.;   A0y = 0.;
        }

        // Coordinate transformation
        // x'_0 = R_z(phiz) R_y(thetaz) x, y'_0 = R_z(phiz) R_y(thetaz) y, z' = R_z(phiz) R_y(thetaz) z 
        // x' = R_z'(phix)  x'_0,   y' = R_z'(phix)  y'_0,
        // --> x' = cos(phix) x'_0 + sin(phix) y'_0, y' = -sin(phix) x'_0 + cos(phix) y'_0,
        array<double, Ndim> trans_x0 = {cos(thetaz)*cos(phiz), cos(thetaz)*sin(phiz), -sin(thetaz)};
        array<double, Ndim> trans_y0 = {-sin(phiz), cos(phiz), 0.};

        for (int i = 0; i < Ndim; ++i)
        {
            xdir[i] =  cos(phix) * trans_x0[i] + sin(phix) * trans_y0[i];
            ydir[i] = -sin(phix) * trans_x0[i] + cos(phix) * trans_y0[i];
        }
    };

    /// generate envelope of pulse
    /** 
     * A(t) = A0 f(t) cos(w0t - cep) = E0/w0 f(t) cos(w0t - cep),
     * E(t) = -dA(t)/dt = - E0/w0 f'(t) cos(w0t - cep) - E0 f(t) sin(w0t - cep).
     * @param _t    time. envelope has maximum value when _t=0. Therefore, t0 has to be subtracted before the function.
     * @param _isDt if true return time derivative of envelope. if false return envelope
     */
    double f_envelope(double _t, bool _isDt)
    {
        switch (envtype)
        {
        case ENVELOPE_COS4:
            if (abs(_t) > twidth) return 0;
            if (_isDt) return -4.* envfactor *sin(_t* envfactor) * pow( cos(_t* envfactor ), 3. );
            else return pow( cos(_t *envfactor ), 4. );
        
        // default is gausssian
        default:
            //if (_isDt) return -2*gaussian_factor* _t/(twidth*twidth)*exp( -gaussian_factor* _t*_t/(twidth*twidth) );
            //else return exp( -gaussian_factor* _t*_t/(twidth*twidth) );
            if (_isDt) return 2* _t * envfactor *exp( envfactor * _t*_t );
            else return exp( envfactor * _t*_t );
        }
    };
    /// generate vector potential in laser 2d-plane
    /// Ax(t=t0) = 0 when cep = 0
    complex avlaser_2d( double _t)
    {
        double kt = _t - t0;
        return complex(-A0x * sin(kt * w0 - cep ) * f_envelope(kt, false),
                A0y * cos(kt * w0 - cep ) * f_envelope(kt, false));
    };
    /// generate electric field in laser 2d-plane
    /// E0x = Ex(t=t0) when cep = 0
    complex elaser_2d( double _t )
    {
        double kt = _t - t0;
        return complex(E0x*cos(kt* w0 - cep ) * f_envelope(kt, false)
                - A0x * sin(kt * w0 - cep ) * f_envelope(kt, true),
                E0y * sin(kt * w0 - cep ) * f_envelope(kt, false)
                + A0y * cos(kt * w0 - cep ) * f_envelope(kt, true));
    }
    array<double, Ndim> avlaser(double _t)
    {
        complex out2d = avlaser_2d(_t);
        return {real(out2d)*xdir[0] + imag(out2d)*ydir[0],
            real(out2d)*xdir[1] + imag(out2d)*ydir[1],
            real(out2d)*xdir[2] + imag(out2d)*ydir[2] };
    }
    array<double, Ndim> elaser(double _t)
    {
        complex out2d = elaser_2d(_t);
        return {real(out2d)*xdir[0] + imag(out2d)*ydir[0],
            real(out2d)*xdir[1] + imag(out2d)*ydir[1],
            real(out2d)*xdir[2] + imag(out2d)*ydir[2] };
    }

    void Print_PulseInfo()
    {
        cout << "-------------------------------------------------------------------\n";
        cout << "E0 = " << E0 << ", I0 = " << I0 << "(W/cm^2) , w0 = " << w0 << ", ellip = " << ellip << ", ncycles = " << ncycles << ", cep = " << cep << endl;
        cout << "t0 = " << t0 << ", phix = " << phix << ", (theta, phi)z = (" << thetaz << ", " << phiz << "), envelope = " << envname << endl;
        cout << "-------------------------------------------------------------------\n";
    }
};

class laser 
{

public:

	vector<Pulse> pulses;       ///< collection of pulses
    
    int Nt;                     ///< number of time steps


	double dt;					///< size of time step
	double blaser;				///< Offset time before of the pulses
	double alaser;				///< Offset time after of the pulses
    
	double minus0, major0;		///< minimum and maximum time in pulse active domains
    double atmin, atmax;        ///< minimum and maximum time in total domain (include offset)
    array<double, Ndim> initial_avec, initial_evec;     // store the value at t = atmin
    
    
    void Initialize(double _dt, double _blaser, double _alaser);
    
    array<double, Ndim> avlaser( double _t );
    array<double, Ndim> elaser( double _t );

    void Print_LaserInfo();
    
    
};

void laser::Initialize(double _dt = 0.1, double _blaser = 0., double _alaser = 0.)
{
    dt = _dt;   blaser = _blaser;   alaser = _alaser;

    // set time domain
    vector<double> startTime, endTime;
    for (const auto &el : pulses)
    {
        startTime.push_back(el.t0 - el.twidth );
        endTime.push_back(el.t0 + el.twidth );
    }
    minus0 = *min_element(startTime.begin(), startTime.end());
    major0 = *max_element(endTime.begin(), endTime.end());
    atmin = minus0 - blaser;
    atmax = major0 + alaser;

    initial_avec = avlaser( atmin );
    initial_evec = elaser( atmin );
    
    Nt   = floor( (atmax - atmin)/dt );
}


array<double, Ndim> laser::elaser( double _t )
{
    array<double, Ndim> outVal = {0., 0., 0.};
    for (auto &el : pulses)
    {
        array<double, Ndim> pulseVal = el.elaser(_t);
        for (int i = 0; i < Ndim; ++i) outVal[i] += pulseVal[i];
    }
    for (int i = 0; i < Ndim; ++i) outVal[i] -= initial_evec[i];
    return outVal;
}
array<double, Ndim> laser::avlaser( double _t )
{
    array<double, Ndim> outVal = {0., 0., 0.};
    for (auto &el : pulses)
    {
        array<double, Ndim> pulseVal = el.avlaser(_t);
        for (int i = 0; i < Ndim; ++i) outVal[i] += pulseVal[i];
    }
    for (int i = 0; i < Ndim; ++i) outVal[i] -= initial_avec[i];
    return outVal;
}

void laser::Print_LaserInfo()
{
    cout << "\n\n========================================";
    cout << "\nNtime = " << Nt << endl;
    for (int i = 0; i < pulses.size(); ++i)
    {
        cout << "Pulse " << i+1 << endl;
        pulses[i].Print_PulseInfo();
    }
    cout << "=========================================== \n\n";
}
