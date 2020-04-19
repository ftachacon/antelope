/*
 *  laser.cpp laser pulses train differents characteristics by pulse
 *
 *  Created by Alexis Chacon
 *
 */

#ifndef LASER_H
#define LASER_H
#include <stdlib.h>
#include <string>
#include "timegrid.h"
#include "timeobject.h"
#include "constant.h"
#include <vector>
using namespace std;

const int gaussian_width_multiply = 2.12;
const double gaussian_factor = gaussian_width_multiply*gaussian_width_multiply* 4*log(2);
enum Envelope
{
    ENVELOPE_GAUSSIAN,
    ENVELOPE_SIN2,
    ENVELOPE_RSIN2,
    ENVELOPE_REC,
    ENVELOPE_CONST
};

struct LaserParam
{
public:
    double I0;			                    //Intensity per pulse
	double e;			                    //Ellipticity per pulse
    double E0;                              //Electric Field Magnitude


    double E0x;			                    //Electric Field Component in the x-direction
	double E0y;			                    //Electric Field Component in the y-direction

    double w0;			                    //Frequency per pulse
	double period0;		                    //Period per pulse
	double cycles0;		                    //Cycles number per pulse


	double cep0;		                    //Carrier envelope phase per pulse
	double phi_x;		                    //Additional phase in the x-direction given by elliptical polarization 
    double phi_y;                           //Additional phase in the y-direction given by elliptical polarization 

    double twidth;		                    //Time bandwidth per pulse
    double t0;                              // Mid-time of pulses
    double theta0;
    
	Envelope envelope;	//Envelope type used

    // Essential initial values - I0, e, w0, cycles0, cep0, phi_rel, t0, theta0
    void Initialize(double _E0, double _e, double _w0, int _ncycle, double _cep,  double _t0, double _theta0, string _env_name)
    {
        E0 = _E0;   e = _e; w0 = _w0;   cycles0 = _ncycle;  cep0 = _cep;  t0 = _t0;   theta0 = _theta0;
        if (_env_name == "sin2")    envelope = ENVELOPE_SIN2;
        else    envelope = ENVELOPE_GAUSSIAN;
    }
    void PreProcessing()
    {
        
        I0 = E0*E0 * 3.5e16;

        period0    =  dospi/w0; // laser period or cycle
        
        // laser time duration
        switch (envelope)
        {
        case ENVELOPE_SIN2:
            twidth = cycles0*period0;
            break;
        
        default:
            twidth = gaussian_width_multiply * cycles0*period0;
            break;
        }
        
        // The ellipticity is defined by e = b/a. and -1 <= e <= 1
        // theta0 is degree between major axis and x-axis
        // If there is no additional phase, E(t=t0) = {E0 cos(theta0), E0 sin(theta0)}
        // i.e. additional phase is adjusted to electric field is at major axis at t=t0
        // Also note that theta0 and pi-theta0 represents same major axis, but phase of electric field is changed
        // i.e. E(theta0 = pi-at0) = -E(theta0 = at0)
        // If e = 0, the laser is circularly polarized independent of value of theta0
        // however, theta0 give additional CEP. additional CEP = theta0 in this case
        // e < 0 --> right helicity, e > 0 --> lefet helicity, e = 0 --> linear polarization

        // Note that here we assume that laser filed is propagating in +z direction
        // i.e. E is proportional to cos(kz - wt + phase), so cos(wt - phase) is used in code.

        // For e= 0 and phi_rel = 0, the laser is linearly polarized along x-direction
        // For e= 0 and phi_rel = 0, and theta != 0 (different to zero), the laser is linearly
        //polarized along a stringht line with theta angle with respect to positive x-direction
        // For e= 1 and phi_rel = pi/2., the laser is circularly polarized
        // Between e = (0, 1) and phi_rel = pi/2., the laser is ellipticaly polarized with major axis along x
        
        // not fully tested yet!

        double temp_a = 1.0 / sqrt(1.0 + e*e) * E0;
        double temp_b = temp_a * e;

        E0x = sqrt(pow(temp_a*cos(theta0), 2) + pow(temp_b*sin(theta0), 2));
        E0y = sqrt(pow(temp_a*sin(theta0), 2) + pow(temp_b*cos(theta0), 2));

        phi_x = atan2(-temp_b*sin(theta0), temp_a*cos(theta0));
        phi_y = atan2( temp_b*cos(theta0), temp_a*sin(theta0));
    };
};

class laser 
{

public:

	
    
    int Npulses;				//Total pulses number
	int Nmaxt;					//Total iteration number
    int NewNt;

    vector<LaserParam> PulseParam; // paramters for pulses
    
	double t01;					//Initial time to first pulse
	double dt;					//Time step
	double blaser;				//Time before of the pulses
	double alaser;				//Time after of the pulses

    double ncycles2;
    
	double major0;				//Major value
	double minus0;				//Minus value
    double width;
    complex a_ef0,a_ef1;
    complex a_af0,a_af1;
    double atmin, atmax;
    
    
	vector<double> clock0;      //Initial time per pulse
	vector<double> clock1;      //The time until half per pulse
	vector<double> clock2;		//The time until final per pulse
	
	
    
    
	vector<int> kclock0;		//k-th value to the initial time per pulse 
	vector<int> kclock1;		//k-th value to the time until half per pulse
	vector<int> kclock2;		//k-th value to the time until final per pulse
    

    timegrid g;					//Object of the time grid
	timeobject efield;			//Object to the electric field for total pulses
	timeobject avector;			//Object to the vector potential for total pulses

    
    
	timeobject *ef;				//Object to the electric field per pulse
	timeobject *env;			//Object to the envelope per pulse
	timeobject *av;				//Object to the vector potential per pulse
	
    
    
	timeobject *av_int;			//Object to the integral of the vector potential per pulse
	timeobject *avsq_int;		//Object to the square integral of the vector potential per pulse	
	
    
    
    
	/*==========================*/
	/*      MAIN FUNCTIONS
	/*==========================*/
    laser(){};
    laser(int _Npulses);
	void init_laser(int _Npulses);      							                                // Creator Object
	//~laser();														//Destructor
	//void laser_pulses(double _dt, double _t01, double _blaser, double _alaser);	        // LASER PULSES TRAIN
    void laser_pulses(double _dt, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
    
    
	/*=========================*/
	/*   SECUNDARY FUNCTIONS
	/*=========================*/
	inline void Initialize_Amplitude_Period();			// Initialize max amplitude and period
	inline void Set_StartTime_EndTime();		        // "Start" and "end" time per pulse
	
    
    void Laser_Grid();									// Object time axis and field_set
	void Evaluation_Laser_Pulse();						// Evaluation pulses train
	
    
    void Sum_Pulses();									// Pulses sum
	void Set_Vector_Potential();						// Vector potential per pulse
	
    
    void Vector_Potential();							// Total vector potential
	void put_on_envelope();								// Generator of Envelope of the pulses

    double f_envelope(double _t, LaserParam const *_param1, bool _isDt);    // Generate enevelop function
    
    
    complex avlaser( double const *t );
    complex elaser( double const *t );
    inline complex avlaser( double const *t , LaserParam const *_param1);
    inline complex elaser( double const *t , LaserParam const *_param1);
    
    
    void Set_Av_Integral();
	inline void Set_index_pulses();	
	
	
	double qmajor(vector<double>& v);					//Finding the maximum value of a vector
	double qminus(vector<double>& v);					//Finding the minimum value ot a vector
    
    void laser_outs_EA( FILE *laserout, int skiper );
    void copy_laser( laser *lcopy );

    void Print_LaserInfo();                              // Printing information of laser
    
    
};


//==============================================================================//
			/*=== MAIN FUNCTIONS ===*/

/*===========================================================          		
		/*=== OBJECT'S LASER CONSTRUCTOR  ===*/
laser::laser(int _Npulses)
{
    init_laser(_Npulses);
}
void laser::init_laser(int _Npulses)
{
    Npulses = _Npulses;

    clock0.resize( Npulses, 0.0 );
    clock1.resize( Npulses, 0.0 );
    clock2.resize( Npulses, 0.0 );
    
    kclock0.resize( Npulses, 0.0 );
    kclock1.resize( Npulses, 0.0 );
    kclock2.resize( Npulses, 0.0 );
    
    PulseParam.resize(Npulses);
    
    //ef  = new timeobject[Npulses];			// reserve memory
    //env = new timeobject[Npulses];			// reserve memory
    //av  = new timeobject[Npulses];
    
    
    
    //av_int    = new timeobject[Npulses];		// reserve memory
    //avsq_int  = new timeobject[Npulses];		// reserve memory	  	  
    
    a_ef0 = 0.;
    a_af0 = 0.;
    a_ef1 = 0.;
    a_af1 = 0.;
}//End initialize variable  




/*void laser::copy_laser( laser *lcopy )
{
    
    if (lcopy->Npulses==NULL)
    {
        printf ("\n*************************\nError with the number of pulses\n*************\n");
        exit (EXIT_FAILURE);
    }
    
    for ( int i = 0; i < lcopy->Npulses; i++ )
    {
        
        
        I0[i]            = lcopy->I0[i];        // Intensity W/cm^2
        w0[i]            = lcopy->w0[i];        // Central frequency
        
        
        cycles0[i]       = lcopy->cycles0[i];  // Cycles number
        cep0[i]          = lcopy->cep0[i];     // Carrier Envelop Phase
        
        
        e[i]             = lcopy->e[i];        // Elliptical of the pulse
        phi_rel[i]       = lcopy->phi_rel[i];  // Relative phase between the polarization Ex and Ey
        
        
        theta0[i]        = lcopy->theta0[i] ;   // Angle between of the laser with respect to x direction in radians
        
        
        envelope[i]      = lcopy->envelope[i];    // Envelop name
        
        
        
        if ( i > 1 )
            delay0[i-1]     = lcopy->delay0[i-1];
        
        
    }
    

    // Making the linear polarization pulse Ex
    laser_pulses( lcopy->dt, lcopy->t01,  lcopy->blaser, lcopy->alaser );  //Routine that computes the laser
    
    
}*/



/*===========================================================          		
 /*===  LASER PULSES FUNCTION  ===*/
void laser::laser_pulses(double _dt, double _blaser, double _alaser)
{
    
    
   	/*======= Pulse's Parameter ======*/
   	//t01     = _t01;									//	Start time first pulse
   	dt      = abs( _dt );			                    //	Time step
    
    
    
    blaser  = abs(_blaser);	                        //	Time before the laser or the pulse train
   	alaser  = abs(_alaser);	                        //	Time after the laser or the pulse train
	
    
    
    
   	//================================//
	Initialize_Amplitude_Period();					//	Initialize max amplitude and period        
   	Set_StartTime_EndTime();	 	     	        //	"Start" and "end" time by each pulse
	
    
    
    
   	major0	= qmajor(clock2);	     	                //	Maximum time to all pulse
   	minus0	= qminus(clock0);          	            //	Minimun time to all pulse     
	
    
    
    
	Set_index_pulses();
	
    
    
    
    	Laser_Grid();       	     	                //	Objet time axis
   //	Evaluation_Laser_Pulse();          	            //	Evaluation pulses
	
    
    
    
    //	Sum_Pulses();	     	     	                //	Pulses sum (build pulses train)
   	//Set_Vector_Potential();		     	            //	Vector potential by each pulse
    
    
    
    
   	//Vector_Potential();								//	Vector potential to pulses train
    
    

}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/***********************************************************/
		/*========= SECUNDARY FUNCTIONS ===========*/
/***********************************************************/

//== Function initialize max amplitude and period ==// 
inline void laser::Initialize_Amplitude_Period()
{
	for( int kpulse = 0; kpulse < Npulses; kpulse++ )
	{
        PulseParam[kpulse].PreProcessing();
	}

}


/***********************************************************/
  //== Function "start" and "end" time by each pulse ==//
/***********************************************************/
inline void laser::Set_StartTime_EndTime()
{
    
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
        
		/*if(kpulse==0)
		{
            
			clock0[kpulse]   = t01;
            
			clock1[kpulse]   = clock0[kpulse] + twidth[kpulse]/2.;
            
			clock2[kpulse]   = clock0[kpulse] + twidth[kpulse];
            
		}
		else
		{
            
			clock0[kpulse]  = clock2[kpulse-1]
            
								+ delay0[kpulse-1]
            
			                    - (twidth[kpulse-1] + twidth[kpulse])/2.;
            
            
			clock1[kpulse]  = clock0[kpulse] + twidth[kpulse]/2.;
			
            clock2[kpulse]  = clock0[kpulse] + twidth[kpulse];
			
		}*/
        clock0[kpulse] = PulseParam[kpulse].t0 - PulseParam[kpulse].twidth;
        clock1[kpulse] = PulseParam[kpulse].t0;
        clock2[kpulse] = PulseParam[kpulse].t0 + PulseParam[kpulse].twidth;
        
	}//End loop
 
    
}//End function



/***********************************************************/
    /* Setting indexes to maxima,
                    intital
                        and final time for each pulse*/
/***********************************************************/
inline void laser::Set_index_pulses()
{
	
	
    for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
        
		kclock0[kpulse] = floor( abs( clock0[kpulse] - (minus0-blaser) )/dt )  + 1;
        
		kclock1[kpulse] = floor( abs( clock1[kpulse] - (minus0-blaser) )/dt )  + 1;
        
		kclock2[kpulse] = floor( abs( clock2[kpulse] - (minus0-blaser) )/dt )  + 1;
        
	}//End loop
    
	

    
}//end of setting pulses indexes


/***********************************************************/
                //== Time axis and field set  ==//
/***********************************************************/
void laser::Laser_Grid()   
{
    
    double axisinitialtime  = minus0-blaser;
    double timelength       = major0 + alaser - axisinitialtime;
    
    
	/*Nmaxt = ceil(  timelength /dt );
    
	if ( Nmaxt%2!=0 )
		Nmaxt = Nmaxt+1;
	
    
    
    //setting time grid
	g.set_grid( Nmaxt
               , dt
               , axisinitialtime );
    
    
    
    //setting array for time axis
	efield.put_on_grid( g );
	avector.put_on_grid( g );
	
    
    
	//Setting fields and evelopes memories
    for (int kfield=0; kfield <Npulses;kfield++)
	{
        
		env[kfield].put_on_grid(g);
		
        ef[kfield].put_on_grid(g);
		
        av[kfield].put_on_grid(g);
        
	}*/

    a_ef1=0.;
    a_af1=0.;
    
    /*width   = twidth[0];
    atmin   = -1.35*width;
    atmax   = +1.35*width;
    double tempMin=-width;*/
    double tempMin = minus0;
    atmin = axisinitialtime;
    atmax = major0 + alaser;
    //*/
    /*
    width   = twidth[0]/( 2.*sqrt( 2.*log(2.) ) );
    atmin   = -5.*width;
    atmax   = +5.*width;
    
    double tempMin=atmin;//*/

    a_ef1 = elaser(  &(atmin) );
    a_af1 = avlaser( &(atmin) );
    
    NewNt   = floor( (atmax - atmin)/dt );

    if ( NewNt%2 !=0 ) NewNt+=1; 
    
}
//== End of Time axis and field set  ==//
/***********************************************************/


/***********************************/
//==Building pulses train ==//
/***********************************/
// void laser::Evaluation_Laser_Pulse()
// {
    
    
// 	double arg1;  // Phase x
// 	double arg2;  // Phase y
	
    
    
// 	put_on_envelope( );
	
    
// 	//Evaluating train of pulses
//     for( int kpulse=0; kpulse<Npulses; kpulse++ )
// 	{
        
        
// 		for( int ktime=0; ktime<Nmaxt; ktime++ )
// 		{
            
// 			arg1 = w0[kpulse]*( g.t[ktime] - clock1[kpulse] ) - cep0[kpulse];           // Phase x
			
//             arg2 = arg1 + phi_rel[kpulse];                                              // Phase y
				
            
//             ef[kpulse].f[ktime]  = complex(
                                           
//                                            real( env[kpulse].f[ktime] )*sin( arg1 ), //x-laser field component
                                           
//                                            imag( env[kpulse].f[ktime] )*sin( arg2 )  //y-laser field component
                                           
//                                            );
            
			
// 			av[kpulse].f[ktime] = ef[kpulse].f[ktime] ;

//         }
        
//     }
    
    
// }//End of pulses
/***********************************/
//E(t) = (E0x,E0y)*exp( aarg1 )*sin( arg1 )
/***********************************/





// /*************************************************
//                 Building envelope
//         *************************************************/
// void laser::put_on_envelope()
// {
// 	double arg0;	
// 	double gaussian_factor=64.*log10(2.);
//     double fgauss;
//     double ts;

    
// 	vector<string> envelopes (5);
	
// 	envelopes[0]="rect";
// 	envelopes[1]="sin2";		
// 	envelopes[2]="rsin2";
// 	envelopes[3]="gauss";	
// 	envelopes[4]="konst";	
	
    
// 	for(int kpulse=0;kpulse<Npulses;kpulse++)
// 	{
        
        
        
//         /****************************/
//         //Rectangle Envelope
//         /****************************/
// 		if( envelope[kpulse] == ENVELOPE_REC )
// 		{
            
            
// 			for ( int ktime = 0; ktime < Nmaxt; ktime++ )
// 			{
				
// 				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse])
//                 {
                
//                     env[kpulse].f[ktime] = 0.;
                    
//                 }
// 				else
//                 {
                    
//                     env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] );
                
//                 }
				
				
// 			}
            
            
// 		}
//         /****************************/
//         //End Rectangle Envelope
//         /****************************/

        
        
        
        
//         /****************************/
//         //sin2 Envelope
//         /****************************/
// 		if ( envelope[kpulse] == ENVELOPE_SIN2 )
            
            
// 			for( int ktime = 0; ktime < Nmaxt; ktime++ )
// 			{
            
                
//                 ts      =  g.t[ktime] - clock0[kpulse] ;
// 				arg0    = w0[kpulse] * ts/cycles0[kpulse]/2.;

                
//                 env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] )*sin(arg0)*sin(arg0);

				
// 			}
//         /****************************/
// 		    //End sin2 Envelope//
//         /****************************/
        
        
        
        
//         /*************************************************/
// 		//Rectangle Multiplied by a sin2 Envelope (rsin2)
//         /*************************************************/
// 		if (envelope[kpulse]==ENVELOPE_RSIN2) 
// 		{
            
            
            
// 			for(int ktime=0;ktime<Nmaxt;ktime++)
// 			{
				
                
                
// 				if ( g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse] )
//                 {
                    
                    
//                     env[kpulse].f[ktime] = 0.0;
                    
                    
// 				} 
// 				else
// 				{
                
//                     ts      = g.t[ktime] - clock0[kpulse] ;
                    
// 					arg0    = w0[kpulse]*ts/cycles0[kpulse]/2.;
					
//                     env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] )*sin(arg0)*sin(arg0);

                    
// 				}
                
                
                
// 			}
            
            
            
// 		}/********************************************************/
//              //End Rectangle  Multiplied sin2 Envelop (rsin2)
//         /********************************************************/
        
        
        
//         /*************************************************/
//                     //Gaussian Envelope
//         /*************************************************/
// 		if (envelope[kpulse]==ENVELOPE_GAUSSIAN)
			
		  
		  
// 			for(int ktime=0;ktime<Nmaxt;ktime++)
// 			{
				
//                 ts      = g.t[ktime] - clock1[kpulse] ;
                
                
//                 arg0    = gaussian_factor*ts*ts/twidth[kpulse]/twidth[kpulse];
                
                
// 				fgauss = exp( -arg0 );
                
                
//                 env[kpulse].f[ktime] = complex( E0x[kpulse] ,E0y[kpulse] )*fgauss;
                
                
// 			}
//         /*************************************************/
//                     //End Gaussian Envelope
//         /*************************************************/
        
        
        
//         /*************************************************/
//                     //Constant Envelope
//         /*************************************************/
// 		if (envelope[kpulse]==ENVELOPE_CONST) 
// 		{
            
// 			for(int ktime=0;ktime<Nmaxt;ktime++)
// 			{
                
//                 env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] );

// 			}
            
// 		}//End Constant Envelope
        
// 	}//End number of pulses
    
    
//     /*
//      a_af0 =  complex( -E0x[0]*sin( (*t)* w0[0] - cep0[0] )
//      , E0y[0]*cos( (*t)* w0[0] - cep0[0] )
//      )*(1.-pow( sin(  (*t)* w0[0]/4./cycles0[0]/1.5 ),2.))/w0[0];// - a_af1;
     
     
//      }//*/

    
    
    
// }
/*************************************************/
                //End of laser or pulse train envelopes dk = Length/(Nx-1)
        /*************************************************/

double laser::f_envelope(double _t, LaserParam const * p1, bool _isDt)
{
   
    switch (p1->envelope)
    {
    case ENVELOPE_SIN2:
        if (abs(_t) > p1->twidth) return 0;
        if (_isDt) return -p1->w0*sin(p1->w0*_t/2./p1->cycles0)/(4*p1->cycles0);
        else return pow( cos(_t* p1->w0/(4*p1->cycles0) ), 2. );
        break;
    
    default:
        if (_isDt) return -2*gaussian_factor* _t*exp( -gaussian_factor* _t*_t/p1->twidth/p1->twidth )/(p1->twidth*p1->twidth);
        else return exp( -gaussian_factor* _t*_t/p1->twidth/p1->twidth );
        break;
    }
}

complex laser::elaser( double const *t )
{
    //if (*t < minus0 || *t > major0) return 0;
    a_ef0=0.;
    double kt;
    LaserParam *p1;
    for (int kpulse = 0; kpulse < Npulses; ++kpulse)
    {
        p1= &(PulseParam[kpulse]);
        kt = *t - p1->t0;
        a_ef0 += elaser(&kt, p1);
    }
    return a_ef0-  a_ef1;
}
inline complex laser::elaser( double const *t, LaserParam const *p1)
{
    return complex(p1->E0x*cos((*t)* p1->w0 - p1->cep0 - p1->phi_x) * f_envelope((*t), p1, false)
            + p1->E0x/p1->w0 * sin((*t) * p1->w0 - p1->cep0 - p1->phi_x) * f_envelope((*t), p1, true),
            p1->E0y * cos((*t) * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope((*t), p1, false)
            + p1->E0y/p1->w0 * sin((*t) * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope((*t), p1, true));
}
// Value of phase is adjust to E0x = E0(t=t0, cep = 0)
// So E(cos, sin) --> A(sin, -cos)
complex laser::avlaser( double const *t )
{
    //if (*t < minus0 || *t > major0) return 0;
    a_af0=0.;
    double kt;
    LaserParam *p1;
    for (int kpulse = 0; kpulse < Npulses; ++kpulse)
    {
        p1= &(PulseParam[kpulse]);
        kt = *t - p1->t0;
        a_af0 += avlaser(&kt, p1);
    }
    return a_af0-  a_af1;
}
inline complex laser::avlaser( double const *t, LaserParam const *p1)
{
    return complex(-p1->E0x/p1->w0 * sin((*t) * p1->w0 - p1->cep0 - p1->phi_x) * f_envelope((*t), p1, false),
            -p1->E0y/p1->w0 * sin((*t) * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope((*t), p1, false));
}


//== Function pulses sum (build pulses train) ==//
void laser::Sum_Pulses()
{
    
	//Start loop sum pulses
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
        
        efield.f[ktime] = 0.;

		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
            
            efield.f[ktime]+= ef[kpulse].f[ktime];

		}
        
	}//End sum set pulses
    
}//Building train of attosecond pulses




//== Function vector potential by each pulse ==//
void laser::Set_Vector_Potential()
{
   
    
	for( int kpulse = 0; kpulse < Npulses; kpulse++ )
	{
        
		av[kpulse].integrateRK4();			//av[kpulse].integrate();
        
		complex aa= av[kpulse].f[0];

		
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
            
            av[kpulse].f[ktime]-= aa;
            
            av[kpulse].f[ktime]*= -1.;
			
		}
        
	}
    
    
}
/*************************************************/
//== End of vector potential for each pulse ==//
/*************************************************/



//== Function Integral of the vector potential by each pulse ==//
void laser::Set_Av_Integral()
{
    
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
        
		av_int[kpulse].put_on_grid(g);
		avsq_int[kpulse].put_on_grid(g);
		
		
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
            av_int[kpulse].f[ktime] = av[kpulse].f[ktime];

			
            avsq_int[kpulse].f[ktime] = complex( real(av[kpulse].f[ktime])*real(av[kpulse].f[ktime]),
                                                 imag(av[kpulse].f[ktime])*imag(av[kpulse].f[ktime]) );

			
		}
		
		av_int[kpulse].integrateRK4();//av[kpulse].integrate(); 
		avsq_int[kpulse].integrateRK4();//av[kpulse].integrate(); 				
		//*/
	}
	
} //End integral of the vector potential





//== Function: vector potential to pulses train ==//
void laser::Vector_Potential()
{
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
        avector.f[ktime] = 0. ;

        
		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
            
            avector.f[ktime]+= av[kpulse].f[ktime];

            
		}	
	}//End sum set pulses       
}



//Start function major
double laser::qmajor(vector<double>& v)
{
	int N=v.size();
	double may = v[0];
	
	for (int k=1;k<N;k++)   
		if (v[k]>may) may=v[k];
	
	return may;
}//End major

//Start function minus
double laser::qminus(vector<double>& v)
{
	int N=v.size();
	double min=v[0];
	
	for (int k=1;k<N;k++)
		if (v[k]<min) min=v[k];	
	
	return min;
}//End minus


//Output for the total electric field and vector potential, x and y components
void laser::laser_outs_EA( FILE *laserout, int skiper=1 )
{
    for ( int ktime=0; ktime < Nmaxt/skiper; ktime++)
        fprintf( laserout,"%.16e %.16e %.16e %.16e %.16e\n",
                g.t[ ktime*skiper ],
                real( efield.f[ ktime*skiper ] ),
                imag( efield.f[ ktime*skiper ] ),
                real( avector.f[ ktime*skiper ] ),
                imag( avector.f[ ktime*skiper ] )
                );
    
    fflush( laserout );
}

void laser::Print_LaserInfo()
{
    cout << "\n\n*************************";
	cout << "\nNtime = " << Nmaxt << endl;
	cout << "Memory = "<< ( ( 2 + 3*Npulses )*24 + 8 )*Nmaxt*1e-6  <<  "  Mb" << endl;
    cout << "*************************\n\n";
}

/*========================== END  SECUNDARY FUNCTIONS =============================*/

#endif
//END
