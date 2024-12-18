//
//  testing.cpp
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/20/19.
//

//#include "testing.hpp"


// Standar C/C++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <iomanip>
//#include "mpi.h"
//#include "mkl.h"


// My Own Headers
#include "constant.h"
#include "timegrid.h"
#include "timeobject.h"
#include "laser.h"
#include "momaxis.h"
//#include "operators.h"
#include "solidstructure.h"
//#include "observables.h"
//#include "sbeqs.h"

#define MASTER 0    /* task id of master task or node */
#define energy_factor 27.2

void initialize_laser_parameters( double eparam[],  string ename, laser &lobject );
void setting_complex_memory(complex **_TempPointer, int _Ntime);
void der_pi( double *eg, double *xig_ef, double *T2, complex *ROmega, complex *cohPi0, complex *cpi);
void der_nc( complex *ROmega, complex *cpi, complex *nc );

double iparam   = 1;
double jparam   = 0; 
double kparam   = 0; 
int ncyparam    = 3;
int nxparam     = 1;
int nyparam     = 1;
double dtparam  = 1.;
int main( int argc, char *argv[] )
{
    
    
    iparam        = atof( argv[1] );         //Magnetic flux
    jparam        = atof( argv[2] );         //On Side potential ratio to t2, i.e M/t2
    nxparam       = atoi( argv[3] );         //Nx points
    nyparam       = atoi( argv[4] );         //Nx points
    kparam        = atof( argv[5] );         //Laser field strength 
    ncyparam      = atoi( argv[6] );         //Number of cycles at FWHM 
    dtparam       = atof( argv[7] );         //Time step, param dt
     
    
    cout << "\n\n\n#*============================================================*#\n";
    cout << "       High order harmonic generation in the \n       Topological Haldane Model (HM), Chern Insulator " << endl;
    cout << "       via the numerical integration of the \n       Semiconductor Bloch Equation (SBEs) \n          by Alexis Chacon June 27, 2019....\n";
    cout << "#*=========================================================*#\n\n\n";
    
    
    
    


    //###############################
    //Variables
    int N[Ndim], n;
    int ktime, i, j, l;
    int Nrk4 = 4;
    
    


    
    //Wavefunction Bloch Gauge choice
    int gauge   = 1;
    int gauge2  = 2 ;
    
    
    
    
    /*****************************************s
      HM. Parameters for Honney Comb lattice
     ****/
    double la0            = 1.;         //Lattice constant in Angstroms
    double a0             = la0/0.529;  //Lattice constant in atomic units (a.u.)
    double t1             = 0.075;      //Hopping parameters in a.u.
    double t2             = t1/3.;
    double phi0           = iparam;        //magnetic flux or phase in radians
    double M0             = jparam*t2;
    double chernN0        = 0.;
    double gauge_ky_down  = -0.16;
    double gauge_ky_up    = +0.90;
    
    
    
    /*************************************s
     
     Momentum box via momaxis class
     
     *****/ //
    N[0]        = nxparam;  //Nkx -- number of points on x-direction
    N[1]        = nyparam;  //Nky -- number of points on y-dir
    N[2]        = 1 ;
    



    double kmaxs0[Ndim];
    kmaxs0[0]   = 2.*pi/sqrt( 3. )/a0; // Momentum grid max on x-dir a.u.
    kmaxs0[1]   = 2.*pi/3./a0;         // Momentum grid max on y-dir a.u.
    kmaxs0[2]   = 0. ;
    

    
    /*************************************s
    
    //Laser Parameters
     
     */

    int npulses         = 1;             // Number of pulses
    double wfreq        = 0.014;         // Frequency in atomic units (a.u.)
    double period0      = dospi/wfreq;   // Optical Period in a.u.
    double E0           = kparam;       // Laser field strenth
    double I0           = 0.;           //Intensity
    double ncycles      = ncyparam;     // Number of Optical cycles at FWHM or at 1/e
    double cep          = 0.;           // Carrier envelope phase
    
    //Time step dt
    double dt           = dtparam;           // Time step
    
    double tstart       = 0.;
    double offset       = 8.*period0;
    double outset       = 8.*period0;
    double ellipticity  = 0.;
    double relativephase= 0.;
    double ltheta0      = 0.;           //Inclination angle
    double T2           = period0;      //Dephasing time in the crystal
    
    
    
    
    
    double kx = 1., ky = 1.;
    double dxdt=0., dydt=0.;
    
    
    
    
    
    complex *inter_rad[2];
    complex *intra_rad[2];
    complex *ConductionOccup;
    complex *Coherence_Mom_Int;
    
    
    
    

    double t[Nrk4], at, tmin;
    double nv;
    double tkx[Nrk4], tky[Nrk4], tgv_v[2]={0.,0.}, tgv_c[2]={0.,0.};
    double teg[Nrk4], txig[Nrk4], tef[Nrk4], taf[Nrk4];
    complex tdip_cv[2], tROmega[Nrk4], k_pi[Nrk4], k_nc[Nrk4];
    complex temp_pi0,complex_ef,_half_pi;
    
    
    
    for (i=0; i<Nrk4; i++)
    {
        k_nc[i]=0.;
        k_pi[i]=0.;
    }
    
    
    
    
    complex *rkutta_nc, *rkutta_pi; //Occupation and Coherence variables
    
    
    
    
    
    //Opening output files!
    FILE *laserout, *mout, *fout, *engout;
    FILE *dipout, *conout, *curvaout, *gvelout;
    FILE *inh_out, *occup_out, *sparamout;
    FILE *simulation_out;
    
    
    
    
    laserout        = fopen( "outlaserdata.dat","w"); // Laser field characteristics,time axis,
    mout            = fopen( "mgrid.dat","w");
    fout            = fopen( "derout.dat","w");
    engout          = fopen( "edispersion.dat","w" );
    dipout          = fopen( "dipoles.dat","w");
    conout          = fopen( "connection.dat","w");
    curvaout        = fopen( "curvature.dat","w");
    gvelout         = fopen( "gvelocities.dat","w");
    inh_out         = fopen( "coherence.dat", "w" );
    occup_out       = fopen( "occupationCB.dat", "w" );
    sparamout       = fopen( "setOfparameters.dat", "w" );
    
    simulation_out  = fopen("full_integrated_currents.dat","w");
    
    
    


    
    /********************************************************
     
     Creating LASER PULSE and defining other parameters
     
     *****/
    


    laser fpulse(   npulses );             // Constructor of Laser Pulses

    I0      = E0*E0*3.5e16;             // Intensity W/cm^2

    string env_name = "gauss";          // Name envel: "rect", "sin2", "gauss" or "rsin2"
    


    
    
    //Collecting parameters
    double laserparam[]  = { I0, wfreq, ncycles, cep, ellipticity, relativephase, ltheta0, dt, tstart, offset, outset };

    

    
    //Initializing laser parameters ...
    initialize_laser_parameters( laserparam, env_name, fpulse );
    tmin            = fpulse.atmin;
    

    
    
    
    
    //Output of the laser field
    //fpulse.laser_outs_EA( laserout );
    cout << "\n*=============================*\nLaser Parameters:\n";
    cout << "E0           = "     << E0    << " a.u.  I0 = "  << I0   << endl;
    cout << "w0           = "     << wfreq << " a.u."     << endl;
    cout << "A0           = "     << E0/wfreq << " a.u."  << endl;
    cout << "NCycles      = "     << ncycles << " No. Opt. Cycles"<<endl;
    cout << "Up           = "     << E0*E0/wfreq/wfreq/4. << " a.u." << endl;
    cout << "Time-step dt = "     << dt << " a.u." <<  endl;
    cout << "New-No. Of Total time Steps = " << fpulse.NewNt << "\n---\n";
    
    


    fprintf(sparamout,"\n%s","#Laser-Parameters,  E0,          I0,              w0,             Ncycles,           CEP,     No.-time-steps,  Time-step-dt "); 
    fprintf(sparamout,"\n          %e      %e      %e      %e      %e      %d      %e      \n\n", E0, I0, wfreq,  ncycles, cep, fpulse.NewNt, dt );
    
    
    
    
    //##########################
    //Coherence and occupations
    setting_complex_memory( &rkutta_pi, fpulse.NewNt );
    setting_complex_memory( &rkutta_nc, fpulse.NewNt);
    
    for (int itemp=0; itemp<2; itemp++ )
    {
        setting_complex_memory( &inter_rad[itemp], fpulse.NewNt );
        setting_complex_memory( &intra_rad[itemp], fpulse.NewNt );
    }
    

    setting_complex_memory( &ConductionOccup, fpulse.NewNt );
    setting_complex_memory( &Coherence_Mom_Int, fpulse.NewNt);

    
    
    
    
    
    

    
    
    /*************************************
     
        Momentum axis construction
     
     ***/
    
    momaxis g( N, kmaxs0 ); //momentum object
    
    
    
    
    /*****************************************************
     
      Crystal structure object for the Haldane model
     
     *********/
    solidstructure cs( &g, &a0, &T2 );
    

    

    
    //Output momentum
    g.mom_outputs( mout, 1, 1, 1 );
    
   
    
    /*************************************
     
       Parameters of Haldane Model
     
     *********/
    
    cs.haldane_model_params( &t1, &t2, &phi0, &M0, gauge );
    
    
    
    
    //###############################
    //Computing Berry Features, Energy dispersion, and dipoles...
    cs.lattice_structure( );
    
    
    
    //###############################
    //Energy dispersion output
    cs.energy_vc_output(    engout,   1,  1,  1 );
    cs.dipole_cv_output(    dipout,   1,  1,  1 );
    cs.connection_c_output( conout,   1,  1,  1 );
    cs.curvature_c_output(  curvaout, 1,  1,  1 );
    cs.group_vel_vc_output( gvelout,  1,  1,  1 );
    
    
    chernN0 = cs.ChernNumber( );
    
    //###############################
    cs.Set_Of_B_BGrad( &kx, &ky );
    cs.zBerryCurvaCV( );
    cs.group_velCV( );
    
    
    
    cout << "\n\n\n*=============================*\n";
    cout << "Haldane model structure";
    cout << "\na0           = " << a0 << " a.u. (Lattice Const.)";
    cout << "\nt1           = " << t1 << " a.u.";
    cout << "\nt2           = " << t2;
    cout << "\nM0           = " << M0;
    cout << "\nphi0         = " << phi0 << endl;
    cout << "\ngauge        = " << cs.gauge;
    cout << "\n-----\nTesting Structure quantities at k =(1,1)\neg(1,1)     = " << cs.eg0 << " a.u.";
    cout << "\ngroup_c0     = " << cs.group_c0[0];
    cout << "\ngroup_c1     = " << cs.group_c0[1];
    cout << "\ngroup_v0     = " << cs.group_v0[0];
    cout << "\ngroup_v1     = " << cs.group_v0[1];
    cout << "\nBCon         = ("   << cs.chig0[0]/2. << "," << cs.chig0[1]/2. << ")";
    cout << "\nxdcv         = "    << cs.dip_cv0[0];
    cout << "\nydcv         = "    << cs.dip_cv0[1] << "";
    cout << "\nBCurv        = "    << cs.zBcurva_c0 ;
    cout << "\nChernNo.     = "    << chernN0<< "\n";
    cout << "MinGap  	   = " << cs.MinEGap << " a.u."; 
    cout << "\nMaxGap  	   = " << cs.MaxEGap << " a.u. \n" ; 
    
    
    
    
    //###############################
    //Fresh-dephasing
    cs.T2   = period0;
    cout << "Dephasing, T2  = " << cs.T2 << " a.u.\n\n\n";
    
    
    fprintf(sparamout,"\n\n%s","#Momentum Grid and Haldane Parameters, dkx (a.u.),  dky (a.u.),  Nx,  Ny, a0 (a.u.),  t1(a.u.),  t2 "); 
    fprintf(sparamout,"\n          %e      %e      %d      %d      %e      %e      %e      \n\n", g.dk[0], g.dk[1], g.N[0], g.N[1], a0, t1, t2 );



    fprintf(sparamout,"\n\n%s","#Haldane Parameters, phi0 (rad.),  M0 (a.u.),   Min-eg(a.u.),        Max-eg,          Chern No.,        T2 (a.u.),    gauge "); 


    fprintf(sparamout,"\n          %e      %e      %e      %e      %e      %e      %d      \n\n", phi0, M0, cs.MinEGap, cs.MaxEGap, chernN0, cs.T2, cs.gauge );
    fflush(sparamout );
    


    //###############################
    //terminal Output of grid
    cout <<"\n\n*===================================*\nMomentum Grid Features";
    cout << "\ndkx          = "              << g.dk[0];
    cout << "\ndky          = "              << g.dk[1];
    cout << "\ndkz          = "              << g.dk[2];

    
    cout <<"\n\nNo. of points";
    cout << "\nNx           = "               << g.N[0];
    cout << "\nNy           = "               << g.N[1];
    cout << "\nNz           = "               << g.N[2];
    
    
    cout <<"\n\nMinima momentum points";
    cout << "\nkxmin        = "            << g.k[0][0];
    cout << "\nkymin        = "            << g.k[1][0];
    cout << "\nkzmin        = "            << g.k[2][0];
    cout << "\n\nky         = "             << g.k[1][0];
    cout << "\n/*************************/\n";


    cout << "\n\n\n\n#*********==============**********#\nNtime     = " << fpulse.NewNt << endl;
    cout << "dt        =   " << dt << " a.u."<< endl << endl;
    
    
    
    
    //###############################
   //Saving or Output of laser
    for ( n = 0; n < fpulse.NewNt; n++)
    {
        
        at = fpulse.atmin + fpulse.dt*n;
        
        
        fpulse.elaser( &at );
        fpulse.avlaser( &at );
        
        
        fprintf( laserout,"%e    %e    %e\n", at, fpulse.a_ef0, fpulse.a_af0 );

        
    }

    
    cout << "\n\n/************************************************/\n";
    cout << "Starting Momentum and Time integration loops\n (Nx,Ny,Nt) = ("<< g.N[0] << " , "<< g.N[1] <<" , " <<fpulse.NewNt << " ) \n\n";
    
    cout.precision(5);
    cout << "\nktime     " << "       ( kx  ,  ky )      " << "         nc      " <<  "          pi(a.u.)    " << "  gauge  \n" << scientific;
    //###############################
    //Momentum and time integration
    for( j = 0; j < g.N[1]; j++ )
    {
    
        
        ky          = g.k[1][j];
        
        
        cs.gauge    = gauge;
        if ( ky >=gauge_ky_down && ky <=gauge_ky_up )
            cs.gauge = gauge2;
        
        
        
        
        //Momentum kx direction Loop
        for( i = 0; i < g.N[0] ; i++ )
        {
            
            
            
            
            //#####################################
            //Initial conditions for val.
            // and cond. bands occupations, and coherence
            rkutta_nc[0]          = complex( 0., 0. );
            rkutta_pi[0]          = complex( 0., 0. );
            
            
            
            //#####################################
            //Time integration loop
            for( ktime = 1; ktime<fpulse.NewNt; ktime++ )
            {
                
                
                
                t[0]  = tmin + dt*double(ktime-1);
                t[1] = t[0]  + dt/2.;
                t[2] = t[1];
                t[3] = t[0]  + dt;
                
                
                //#####################################
                //Evaluation of SBEs variables
                    //at tn-1, tn-1 + dt/2., and, tn-1 + dt
                for( l = 0; l < Nrk4; l++ )
                {
                    
                    
                    
                    fpulse.elaser(  &t[l] );
                    fpulse.avlaser( &t[l] );
                        
                        
                    tef[l]      = fpulse.a_ef0;
                    taf[l]      = fpulse.a_af0;
                    tkx[l]      = g.k[0][i] + taf[l];
                        
                        
                        
                    cs.Set_Of_B_BGrad( &tkx[l], &ky );
                    

 
                    
                    if ( l == 0 )
                        {
                        
                            complex_ef  = complex( tef[l], 0.);
                        
                            cs.group_velCV( );
                            cs.zBerryCurvaCV( );
                        
                            tdip_cv[0]      = cs.dip_cv0[0];
                            tdip_cv[1]      = cs.dip_cv0[1];
                        
                            tgv_c[0]        = cs.group_c0[0];
                            tgv_c[1]        = cs.group_c0[1];
                        
                            tgv_v[0]        = cs.group_v0[0];
                            tgv_v[1]        = cs.group_v0[1];
                        
                        }
                    
                    
                    teg[l]      = cs.eg0;
                    txig[l]     = cs.chig0[0]*tef[l];
                    tROmega[l]  = cs.dip_cv0[0]*tef[l];
                        
                        
                    
                    
                }//End loop of evaluation of KR4 vars.
                //#####################################

                
                
                
                
                
                //###############################################################
                //Coherence of the occupations via Runge-Kutta 4th order (RK4)

                der_pi( &teg[0], &txig[0], &cs.T2, &tROmega[0], &rkutta_pi[ktime-1], &k_pi[0] );
                
                
                temp_pi0 = rkutta_pi[ktime-1] + k_pi[0]*dt/2.;
                der_pi( &teg[1], &txig[1], &cs.T2, &tROmega[1], &temp_pi0, &k_pi[1] );
                
                
                temp_pi0    = rkutta_pi[ktime-1] + k_pi[1]*dt/2.;
                der_pi( &teg[2], &txig[2], &cs.T2, &tROmega[2], &temp_pi0, &k_pi[2] );
                
                
                temp_pi0    = rkutta_pi[ktime-1] + k_pi[2]*dt;
                der_pi( &teg[3], &txig[3], &cs.T2, &tROmega[3], &temp_pi0, &k_pi[3] );
                
                
                
                //RK4 Formula
                rkutta_pi[ktime] = rkutta_pi[ktime-1]+dt/6.*(k_pi[0]+(k_pi[1]+k_pi[2])*2.+k_pi[3]);
                
                
                
                
                
                
                //####################################################
                //Occupation via RK4
                der_nc( &tROmega[0], &rkutta_pi[ktime-1], &k_nc[0] );
                _half_pi = (rkutta_pi[ktime-1] + rkutta_pi[ktime])/2.;
                der_nc( &tROmega[1], &_half_pi, &k_nc[1] );
                der_nc( &tROmega[2], &_half_pi, &k_nc[2] );
                der_nc( &tROmega[3], &rkutta_pi[ktime], &k_nc[3] );
                
                
                
                rkutta_nc[ktime] = rkutta_nc[ktime-1]+dt/6.*(k_nc[0]+(k_nc[1]+k_nc[2])*2.+k_nc[3]);
                
                
                
                
                
                
                
                
                
                //##############################
                //Momentum integrals --Inter-intra-band-contribution dipole radiation --
                inter_rad[0][ktime-1]+= 2.*real( conj( tdip_cv[0] )*rkutta_pi[ktime-1] );
                inter_rad[1][ktime-1]+= 2.*real( conj( tdip_cv[1] )*rkutta_pi[ktime-1] );
                
                

                
                
                //##########################################
                //Evaluation of the anomalous velocity and
                //classical velocities
                cs.anomalous_velCV( complex_ef );
                nv = 1. -  real( rkutta_nc[ktime-1] );
                
                
                
                
                
                
                intra_rad[0][ktime-1]+= tgv_v[0] * nv + tgv_c[0] * real( rkutta_nc[ktime-1] );
                intra_rad[1][ktime-1]+= (tgv_v[1] - cs.yanomalous_v0) * nv + (tgv_c[1] - cs.yanomalous_c0) * real(rkutta_nc[ktime-1]);
                

                
                
                
                //##############################################
                //Evaluation of occupations and coherence
                ConductionOccup[ktime-1]+=   rkutta_nc[ktime-1];
                Coherence_Mom_Int[ktime-1]+= rkutta_pi[ktime-1];
               
                
                
                
            }//End of Time integration loop
            
            
            
            //OutPut for diagnostic
            if( i%((g.N[0]-1)/4)==0 )
            {
                cout << "\n" << ktime-2 << "   " << g.k[0][i] << "   " << ky << "       " << real(ConductionOccup[ktime-2]*g.dV) << "       " << real(Coherence_Mom_Int[ktime-2])*g.dV << "           "<<cs.gauge;
            
            }
        
            
        }//End of Loop on x-direction
        
        
        
        cout << "\n---\n";
        
        
        
    }//End of Loop on y-direction
    
    
    
    
    
    //#######################################
    //Output inter, intra band currents,
    //occupation and coherence and vs time

    fprintf(simulation_out,"%s    %s                   %s                    %s                 %s                   %s                    %s                   %s\n ", "#","#time(a.u.)","#Jx_er(a.u.)","#Jy_er(a.u.)","#Jx_ra(a.u.)","#Jy_ra(a.u.)", "#nc(a.u.)","#Coherent(a.u.)"  );

    ktime=0;
    fprintf(simulation_out,"%.16e       %.16e         %.16e         %.16e         %.16e         %.16e         %.16e\n ", fpulse.atmin, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV );
    
    for( ktime = 1; ktime<fpulse.NewNt-1; ktime++ )
    {
        
        
        at = fpulse.atmin + dt*ktime;
        
        
        
        dxdt=( real(inter_rad[0][ktime+1]) - real(inter_rad[0][ktime-1]) )*g.dV/dt/2.;
        
        
        
        
        dydt=( real(inter_rad[1][ktime+1]) - real(inter_rad[1][ktime-1]) )*g.dV/dt/2.;
        
        
        
        
        fprintf(simulation_out,"%.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e\n ", at, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV );
        
        
        
    }
    
    ktime=fpulse.NewNt-2;
    at = fpulse.atmin + dt*ktime;
    fprintf(simulation_out,"%.16e       %.16e       %.16e      %.16e       %.16e      %.16e    %.16e\n", at, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV );
    
    
    cout <<"\nNo. of total time steps = " << fpulse.NewNt << endl;
    
    
    
    
    fclose(  fout     );
    fclose(  conout   );
    fclose(  dipout   );
    fclose(  curvaout );
    fclose(  engout   );
    
    
    for (int itemp=0; itemp < 2; itemp++)
    {
        free( inter_rad[itemp] );
        free( intra_rad[itemp] );
        
    }
    
    free( ConductionOccup );
    free( Coherence_Mom_Int );
    free(rkutta_nc);
    free(rkutta_pi);
    cout << "\n***********************\nEND OF THE PROGRAM\n";
    cout << "=====================\n\n";
    
    
};


//initializing laser parameters
void initialize_laser_parameters( double eparam[],  string ename, laser &lobject )
{
    
    
    lobject.I0[0]            = eparam[0];        // Intensity W/cm^2
    lobject.w0[0]            = eparam[1];        // Central frequency
    
    
    
    lobject.cycles0[0]       = eparam[2];        // Cycles number
    lobject.cep0[0]          = eparam[3];        // Carrier Envelop Phase
    
    
    
    lobject.e[0]             = eparam[4];        // Elliptical of the pulse
    lobject.phi_rel[0]       = eparam[5];        // Relative phase between the polarization Ex and Ey
    lobject.theta0[0]        = eparam[6];       // Angle between of the laser with respect to x direction in radians
    
    
    lobject.envelope[0]      = ename;           // Envelop name
    
    
    // Making the linear polarization pulse Ex
    lobject.laser_pulses( eparam[7], eparam[8],  eparam[9], eparam[10] );  //Routine that computes the laser fields,
   
}


void der_pi( double *eg, double *xig_ef, double *T2, complex *ROmega, complex *cohPi0, complex *cpi)
{
    
    *cpi = -I*( ( *eg + *xig_ef - I/(*T2) ) *( *cohPi0 ) + (*ROmega) );
    
}


void der_nc( complex *ROmega, complex *cpi, complex *nc )
{
    
    *nc = -2.*imag( conj( *ROmega )*( *cpi ) );
    
}

void setting_complex_memory(complex **_TempPointer, int _Ntime)
{
    *_TempPointer      = (complex*)malloc( _Ntime*sizeof(complex) );
    memset( *_TempPointer,     0,    sizeof(complex)*_Ntime );
}
