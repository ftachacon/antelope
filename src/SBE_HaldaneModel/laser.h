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

class laser 
{

public:

	
    
    int Npulses;				//Total pulses number
	int Nmaxt;					//Total iteration number
    int NewNt;

    
    
	double t01;					//Initial time to first pulse
	double dt;					//Time step
	double blaser;				//Time before of the pulses
	double alaser;				//Time after of the pulses

    
    
	double major0;				//Major value
	double minus0;				//Minus value
    double width;
    double a_ef0,a_ef1;
    double a_af0,a_af1;
    double atmin, atmax;
    
    
	vector<double> clock0;      //Initial time per pulse
	vector<double> clock1;      //The time until half per pulse
	vector<double> clock2;		//The time until final per pulse
	vector<double> twidth;		//Time bandwidth per pulse
	
    
    
	vector<int> kclock0;		//k-th value to the initial time per pulse 
	vector<int> kclock1;		//k-th value to the time until half per pulse
	vector<int> kclock2;		//k-th value to the time until final per pulse
	
	
    
	vector<double> I0;			//Intensity per pulse
	vector<double> e;			//Ellipticity per pulse
    vector<double> E0;          //Electric Field Magnitude
	
    
    
    vector<double> E0x;			//Electric Field Component in the x-direction
	vector<double> E0y;			//Electric Field Component in the y-direction
	
	
    
    vector<double> w0;			//Frequency per pulse
	vector<double> period0;		//Period per pulse
	vector<double> cycles0;		//Cycles number per pulse
	
    
    
    
	vector<double> cep0;		//Carrier envelope phase per pulse
	vector<double> phi_rel;		//Relative phase between Ex and Ey components of the electric field of the laser pulse
	
    
    
    vector<double> delay0;		//Delay between consecutive pulses
    vector<double> theta0;
    
    
    
	vector<string> envelope;	//Envelope name used
	
    

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
	laser(int _Npulses);      							                                // Creator Object
	//~laser();														//Destructor
	void laser_pulses(double _dt, double _t01, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
    
    
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
    
    
    void avlaser( double *t );
    void elaser( double *t );
    
    
    void Set_Av_Integral();
	inline void Set_index_pulses();	
	
	
	double qmajor(vector<double>& v);					//Finding the maximum value of a vector
	double qminus(vector<double>& v);					//Finding the minimum value ot a vector
    
    void laser_outs_EA( FILE *laserout, int skiper );
    void copy_laser( laser *lcopy );
    
    
};


//==============================================================================//
			/*=== MAIN FUNCTIONS ===*/

/*===========================================================          		
		/*=== OBJECT'S LASER CONSTRUCTOR  ===*/
laser::laser(int _Npulses)
{
	  Npulses = _Npulses;
    
      clock0.resize( Npulses, 0.0 );
	  clock1.resize( Npulses, 0.0 );
	  clock2.resize( Npulses, 0.0 );
	  
	  kclock0.resize( Npulses, 0.0 );
	  kclock1.resize( Npulses, 0.0 );
	  kclock2.resize( Npulses, 0.0 );
	  
	  twidth.resize(Npulses,0.0);
    
      I0.resize(Npulses,0.0);
	  e.resize(Npulses,0.0);	  
    
      envelope.resize(Npulses);
	  
      for (int kpulse=0; kpulse<Npulses; kpulse++)
		 envelope[kpulse] ="gauss";
    
      E0.resize(Npulses,0.0);
	  E0x.resize(Npulses,0.0);
	  E0y.resize(Npulses,0.0);	  
    
      theta0.resize(Npulses,0.);
	  w0.resize(Npulses,0.0);
	  period0.resize(Npulses,0.0);
      
	  
      
	  cycles0.resize(Npulses,0.0); 
	  cep0.resize(Npulses,0.0);	  
	  phi_rel.resize(Npulses,0.0);
	  
	  
	  
	  ef  = new timeobject[Npulses];			// reserve memory
	  env = new timeobject[Npulses];			// reserve memory
	  av  = new timeobject[Npulses];
     
	  
	  
	  av_int    = new timeobject[Npulses];		// reserve memory
	  avsq_int  = new timeobject[Npulses];		// reserve memory	  	  
	  
      a_ef0 = 0.;
      a_af0 = 0.;
      a_ef1 = 0.;
      a_af1 = 0.;
    
    
	  if (Npulses > 1 )
		  delay0.resize(Npulses-1,0.0);
      
    
    
  }//End initialize variable  




void laser::copy_laser( laser *lcopy )
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
    
    
}



/*===========================================================          		
 /*===  LASER PULSES FUNCTION  ===*/
void laser::laser_pulses(double _dt, double _t01, double _blaser, double _alaser)
{
    
    
   	/*======= Pulse's Parameter ======*/
   	t01     = _t01;									//	Start time first pulse
   	dt      = abs(_dt);			                    //	Time step
    
    
    
    blaser  = abs(_blaser);	                        //	Time before the laser or the pulse train
   	alaser  = abs(_alaser);	                        //	Time after the laser or the pulse train
	
    
    
    
   	//================================//
	Initialize_Amplitude_Period();					//	Initialize max amplitude and period        
   	Set_StartTime_EndTime();	 	     	        //	"Start" and "end" time by each pulse
	
    
    
    
   	major0	= qmajor(clock2);	     	                //	Maximum time to all pulse
   	minus0	= qminus(clock0);          	            //	Minimun time to all pulse     
	
    
    
    
	Set_index_pulses();
	
    
    
    
    Laser_Grid();       	     	                //	Objet time axis
   	Evaluation_Laser_Pulse();          	            //	Evaluation pulses
	
    
    
    
    Sum_Pulses();	     	     	                //	Pulses sum (build pulses train)
   	Set_Vector_Potential();		     	            //	Vector potential by each pulse
    
    
    
    
   	Vector_Potential();								//	Vector potential to pulses train
    
    

}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/***********************************************************/
		/*========= SECUNDARY FUNCTIONS ===========*/
/***********************************************************/

//== Function initialize max amplitude and period ==// 
inline void laser::Initialize_Amplitude_Period()
{
    
    //The ellipticity is defined by e = E0y/E0x.
    //That ellipticity param can vary between -1 and +1
    // For e= 0 and phi_rel = 0, the laser is linearly polarized along x-direction
    // For e= 0 and phi_rel = 0, and theta != 0 (different to zero), the laser is linearly
    //polarized along a stringht line with theta angle with respect to positive x-direction
    // For e= 1 and phi_rel = pi/2., the laser is circularly polarized
    // Between e = (0, 1) and phi_rel = pi/2., the laser is ellipticaly polarized with major axis along x
    
    double factor;
    
    
	for( int kpulse = 0; kpulse < Npulses; kpulse++ )
	{
        
        
        //Electric field amplitude
        E0[kpulse]         =  sqrt( I0[kpulse]/3.5e16 );
        
        period0[kpulse]    =  dospi/w0[kpulse]; // laser period or cycle
        
        twidth[kpulse]       =  cycles0[kpulse]*period0[kpulse]; // laser time duration
        
		factor              =  1.0/( 1.0 + e[kpulse]*e[kpulse] ); // ellipticity factor
        
        
        
        //Condition for linear or elliptical polarized laser field
        if (e[kpulse] == 0 )
        {
            
            
            phi_rel[kpulse] = 0.;

            E0x[ kpulse ]        =  E0[ kpulse ]*cos( theta0[ kpulse ] );
            
            E0y[ kpulse ]        =  E0[ kpulse ]*sin( theta0[ kpulse ] );
            
            
        }
        else
        {
            
           
            E0x[kpulse]        =  sqrt( factor )*E0[kpulse];
            
            E0y[kpulse]        =  e[kpulse]*sqrt( I0[kpulse]*factor/3.5e16 );
            
            
        }
        
	}


    
}


/***********************************************************/
  //== Function "start" and "end" time by each pulse ==//
/***********************************************************/
inline void laser::Set_StartTime_EndTime()
{
    
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
        
		if(kpulse==0)
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
			
		}
        
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
    
    
	Nmaxt = ceil(  timelength /dt );
    
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
        
	}

    
    width   = twidth[0]/(2.*sqrt( 2.*log(2.) ));
    atmin   = -5.*width;
    atmax   = +5.*width;
    NewNt   = floor( (atmax - atmin)/dt );
    
    avlaser( &atmin );
    elaser(  &atmin );
    a_ef1 = a_ef0;
    a_af1 = a_af0;
    
    
    cout << "\n\n*************************";
	cout << "\nNtime = " << Nmaxt << endl;
	cout << "Memory = "<< ( ( 2 + 3*Npulses )*24 + 8 )*Nmaxt*1e-6  <<  "  Mb" << endl;
    cout << "*************************\n\n";
    
    
}
//== End of Time axis and field set  ==//
/***********************************************************/



/***********************************/
//==Building pulses train ==//
/***********************************/
void laser::Evaluation_Laser_Pulse()
{
    
    
	double arg1;  // Phase x
	double arg2;  // Phase y
	
    
    
	put_on_envelope( );
	
    
	//Evaluating train of pulses
    for( int kpulse=0; kpulse<Npulses; kpulse++ )
	{
        
        
		for( int ktime=0; ktime<Nmaxt; ktime++ )
		{
            
			arg1 = w0[kpulse]*( g.t[ktime] - clock1[kpulse] ) - cep0[kpulse];           // Phase x
			
            arg2 = arg1 + phi_rel[kpulse];                                              // Phase y
				
            
            ef[kpulse].f[ktime]  = complex(
                                           
                                           real( env[kpulse].f[ktime] )*sin( arg1 ), //x-laser field component
                                           
                                           imag( env[kpulse].f[ktime] )*sin( arg2 )  //y-laser field component
                                           
                                           );
            
			
			av[kpulse].f[ktime] = ef[kpulse].f[ktime] ;

        }
        
    }
    
    
}//End of pulses
/***********************************/
//E(t) = (E0x,E0y)*exp( aarg1 )*sin( arg1 )
/***********************************/





/*************************************************
                Building envelope
        *************************************************/
void laser::put_on_envelope()
{
	double arg0;	
	double gaussian_factor=64.*log10(2.);
    double fgauss;
    double ts;

    
	vector<string> envelopes (5);
	
	envelopes[0]="rect";
	envelopes[1]="sin2";		
	envelopes[2]="rsin2";
	envelopes[3]="gauss";	
	envelopes[4]="konst";	
	
    
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
        
        
        
        /****************************/
        //Rectangle Envelope
        /****************************/
		if( envelope[kpulse] == envelopes[0] )
		{
            
            
			for ( int ktime = 0; ktime < Nmaxt; ktime++ )
			{
				
				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse])
                {
                
                    env[kpulse].f[ktime] = 0.;
                    
                }
				else
                {
                    
                    env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] );
                
                }
				
				
			}
            
            
		}
        /****************************/
        //End Rectangle Envelope
        /****************************/

        
        
        
        
        /****************************/
        //sin2 Envelope
        /****************************/
		if ( envelope[kpulse] == envelopes[1] )
            
            
			for( int ktime = 0; ktime < Nmaxt; ktime++ )
			{
            
                
                ts      =  g.t[ktime] - clock0[kpulse] ;
				arg0    = w0[kpulse] * ts/cycles0[kpulse]/2.;

                
                env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] )*sin(arg0)*sin(arg0);

				
			}
        /****************************/
		    //End sin2 Envelope//
        /****************************/
        
        
        
        
        /*************************************************/
		//Rectangle Multiplied by a sin2 Envelope (rsin2)
        /*************************************************/
		if (envelope[kpulse]==envelopes[2]) 
		{
            
            
            
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
                
                
				if ( g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse] )
                {
                    
                    
                    env[kpulse].f[ktime] = 0.0;
                    
                    
				} 
				else
				{
                
                    ts      = g.t[ktime] - clock0[kpulse] ;
                    
					arg0    = w0[kpulse]*ts/cycles0[kpulse]/2.;
					
                    env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] )*sin(arg0)*sin(arg0);

                    
				}
                
                
                
			}
            
            
            
		}/********************************************************/
             //End Rectangle  Multiplied sin2 Envelop (rsin2)
        /********************************************************/
        
        
        
        /*************************************************/
                    //Gaussian Envelope
        /*************************************************/
		if (envelope[kpulse]==envelopes[3])
			
		  
		  
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
                ts      = g.t[ktime] - clock1[kpulse] ;
                
                
                arg0    = gaussian_factor*ts*ts/twidth[kpulse]/twidth[kpulse];
                
                
				fgauss = exp( -arg0 );
                
                
                env[kpulse].f[ktime] = complex( E0x[kpulse] ,E0y[kpulse] )*fgauss;
                
                
			}
        /*************************************************/
                    //End Gaussian Envelope
        /*************************************************/
        
        
        
        /*************************************************/
                    //Constant Envelope
        /*************************************************/
		if (envelope[kpulse]==envelopes[4]) 
		{
            
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
                
                env[kpulse].f[ktime] = complex( E0x[kpulse], E0y[kpulse] );

			}
            
		}//End Constant Envelope
        
	}//End number of pulses
    
    
    
    
}
/*************************************************/
                //End of laser or pulse train envelopes
        /*************************************************/

void laser::elaser( double *t )
{
    
    a_ef0 = E0[0]*cos( (*t)* w0[0] - cep0[0] )*exp( -.5* ( *t )*( *t )/width/width ) -a_ef1;
    
}

void laser::avlaser( double *t )
{
    
    a_af0 = -E0[0]/w0[0]*sin( (*t)* w0[0] - cep0[0] )*exp( -.5* ( *t )*( *t )/width/width )-a_af1;
    
    /*-(.5*I*self.alpha0*self.E0*np.sqrt(pi/2.)*(-(special.erfi( (I*t + (self.alpha0**2)*self.w0)/(np.sqrt(2)*self.alpha0))*(np.cos(self.cep0) - I*np.sin(self.cep0))) +
                                                 special.erfi((-I*t + (self.alpha0**2)*self.w0)/(np.sqrt(2.)*self.alpha0))*(np.cos(self.cep0) + I*np.sin(self.cep0))))/np.exp( (self.alpha0**2)*(self.w0**2)/ 2.);
*/
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

/*========================== END  SECUNDARY FUNCTIONS =============================*/

#endif
//END
