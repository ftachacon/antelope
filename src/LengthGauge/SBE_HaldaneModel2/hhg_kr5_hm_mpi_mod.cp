//---------------
//  testing.cpp
//...
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
#include "mpi.h"
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

int const Nrk4 = 4; //Runge-Kutta 4th coheff.
int const Nrk5 = 6; //Runge-Kutta 5th coheff.

void initialize_laser_parameters( double eparam[],  string ename, laser &lobject );
void setting_complex_memory(complex **_TempPointer, long long int _Ntime);
void setting_double_memory(double **_TempPointer, long long int _Ntime );
void der_pi( double const *eg, double const *xig_ef, double const *T2, complex const *ROmega, complex const *cohPi0, complex *cpi);
void der_nc( complex const *ROmega, complex const *cpi, double *nc );


complex der_pi_RK5( double const *eg, double const *xig_ef, double const *T2, complex const *ROmega, complex const *cohPi0);
complex jaco_pi_RK5( double const *eg, double const *xig_ef, double const *T2 );

double der_nc_RK5( complex const *ROmega, complex const *cpi );

void My_MPI_CLX_SUM( complex *in, complex *inout, int *len, MPI_Datatype *dptr );




//void My_2D_Domain_Reconstruction(int rank, int block, int ist, int iend, long long int *xshift, int *rmin_int_rank, int *imin_rank, int *jmin_rank, int *imax_rank, int *jmax_rank, int *xi, int *yj, long long int *iglobal, double *rmin_rank, double *rmax_rank, double *rmin_decimal_rank);
//void My_2D_Domain_Decomposition(int rank, int Nprocessors, int *block, int *ist, int *iend );




int ncyparam    = 3;
int nxparam     = 1;
int nyparam     = 1;



double iparam   = 1.;
double jparam   = 0.; 
double kparam   = 0.; 
double dtparam  = 1.;
double dephasing= 1000.;
double ellip    = 0.;
double treg0    = 0.;
int tgauge2     = 0.;
int  trflag = 0;




//################################
int main( int argc, char *argv[] )
{
    
    
    
    
    //################################
    //Numerical parameter and set
    iparam        = atof( argv[1] );         //Magnetic flux
    jparam        = atof( argv[2] );         //On Side potential ratio to t2, i.e M/t2
    nxparam       = atoi( argv[3] );         //Nx points
    nyparam       = atoi( argv[4] );         //Nx points
    kparam        = atof( argv[5] );         //Laser field strength 
    ncyparam      = atoi( argv[6] );         //Number of cycles at FWHM 
    dtparam       = atof( argv[7] );         //Time step, param dt
    dephasing     = atof( argv[8] );         //Dephasing time
    ellip         = atof( argv[9] );         //Laser field Ellipticity
    treg0         = atof( argv[10] );        //Regularization parameter
    tgauge2       = atoi( argv[11] );        //gauge2s
    trflag        = atoi( argv[12] );        //controlling regul. type, 0, non, 1 local, 2 taylor
    
    
    //String Type of variables, envelope and integration rule
    string env_name     = "gauss";          //Name envel: "rect", "sin2", "gauss" or "rsin2"
    string IntionMethod = "Trapz";       // "Trapz" or "Simpson" or "NoneTrapz"/"NoneSimpson"
    
    
	MPI_Op MPI_CLX_SUM;
    
    
	//Fundamental initialization of MPI
	MPI_Init(&argc, &argv);	
	
    
    
	//Creating our MPI-operation to sum up complex arrays
	MPI_Op_create((MPI_User_function *) My_MPI_CLX_SUM, true, &MPI_CLX_SUM);	
	
    
    //#############################
	//Getting rank and size variables
	int rank=0, size=1;
	int Nprocessors=1,NumberOfDoubles=1,NumberOfComplex=1;
    


    cout << "\n\n\n#*============================================================*#\n";
    cout << "       High order harmonic generation in the \n       Topological Haldane Model (HM), Chern Insulator " << endl;
    cout << "       via the numerical integration of the \n       Semiconductor Bloch Equation (SBEs) \n          by Alexis Chacon June 27, 2019....\n";
    cout << "#*=========================================================*#\n\n\n";
    
    
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );	//Getting the local id or rank process
	MPI_Comm_size( MPI_COMM_WORLD, &size );	//Getting the size or number of proccesses or nodes
	MPI_Barrier( MPI_COMM_WORLD );	
	
    
    
    //#############################
	//Number of processes
	Nprocessors		= size;
	NumberOfDoubles = 1;
	NumberOfComplex = 1;
	

    
    //#############################
	MPI_Barrier( MPI_COMM_WORLD );		
	// -- Printing on terminal the PARAMETERs -- //	
	if (rank==MASTER)
	{	cout << "\n-------------\nPARAMETERS\n";
		cout << "\niparam = "		  	<< iparam ;
		cout << "\njparam = "         	<< jparam ;
		cout << "\nNxparam = "         	<< nxparam ;
		cout << "\nNyparam = "         	<< nyparam ;        
		cout << "\n---\nTotal of Processors= " << Nprocessors << "\n";
	}
	//End of Part 1	
	
    


    //###############################
    //Variables
    int N[Ndim]={1,1,1};
    int n=0;
    int ktime=0, itemp=0, jtemp=0, ltemp=0,ktemp=0;
    int jstart=0, jend=0, jstep=0;
    long long int NTotal=1;


    
    //Wavefunction Bloch Gauge choice
    int gauge1   = 1;
    int gauge2  = tgauge2 ;
    
    //if (treg0==0.)
    //    gauge2=1;
    
    
    
    //#############################
    /*****************************************
     
      HM. Parameters for Honney Comb lattice
     
     ****/
    double la0            = 1.;         //Lattice constant in Angstroms
    double a0             = la0/0.529;  //Lattice constant in atomic units (a.u.)
    double t1             = 0.075;      //NN Hopping parameters in a.u.
    double t2             = t1/3.;      //Amplitude of the Next Neirest Neibour of HM
    double phi0           = iparam;     //Magnetic flux or phase of the NNN hoping parameter
    double M0             = jparam*t2;  //Local or on-site potential
    double chernN0        = 0.;
    
    
     //Controlling momentum gauge variation window
    double gauge_ky_down  =  0.;//-0.16;
    double gauge_ky_up    =  0.;//0.90
    double kmaxs0[Ndim]   = {0.,0.,0.};
    double kyShift[Ngrad] = {0.,0.}; //{-0.38,gauge_ky_down};

    
    
    
    
    //#############################
    /*************************************
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
    double tstart       = 0.;                // Starting time of the laser
    double offset       = 8.*period0;
    double outset       = 8.*period0;
    double ellipticity  = ellip;
    double relativephase= 0.;
    double ltheta0      = 0.;               //Inclination angle
    
    
    double T2           = dephasing;        //Dephasing time in the crystal
    
    double kx = 1., ky = 1.; 
    double dxdt=0., dydt=0.;                //Time-derivative for Jx and Jy inter-band contributions
    
    double  zberrycurvature0=0;
    double  classical_velocity_v[Ngrad]={0.,0.}, classical_velocity_c[Ngrad]={0.,0.};
    double  anomalous_vel_v[Ngrad]={0.,0.}, anomalous_vel_c[Ngrad]={0.,0.};
    double  group_velocity_v[Ngrad]={0.,0.}, group_velocity_c[Ngrad]={0.,0.};
    
    
    double  *group_rad[Ngrad];
    double  *anomalous_rad[Ngrad];
    
    complex *inter_rad[Ngrad];
    complex *intra_rad[Ngrad];
    complex *ConductionOccup;
    complex *Coherence_Mom_Int;

    double t[Nrk4], at=0, nv=0.;
    double tkx[Nrk4], tky[Nrk4];
    double teg[Nrk4], txig[Nrk4];
    double k_nc[Nrk4];
    
    double t5[Nrk5];
    double tkx5[Nrk5], tky5[Nrk5];
    double teg5[Nrk5], txig5[Nrk5];
    //double k_nc[Nrk5];
    
    
    complex tgv_v[Ngrad] = {0.,0.};
    complex tgv_c[Ngrad] = {0.,0.};
    complex tdip_cv[Ngrad]={0.,0.};
    
    complex tef[Nrk4], taf[Nrk4];
    complex k_pi[Nrk4], tROmega[Nrk4];
    complex up_dip[Nrk4], down_dip[Nrk4];
    complex temp_pi0[Nrk4],_half_pi=0.;
    
    complex tef5[Nrk5], taf5[Nrk5];
    complex k_pi5[Nrk5], tROmega5[Nrk5];
    complex temp_pi5[Nrk5], _half_pi5=0.;
    
    
    //Occupation and Coherence variables
    complex *rkutta_pi;
    complex rkutta_pi_ktimeplus=0.;

    double *rkutta_nc;
    double rkutta_nc_ktimeplus=0.;

    
    //Initializing RK4 variables
    for (n=0; n<Nrk5; n++)
    {
        if (n<Nrk4)
        {
            t[n]    = 0.;
            tkx[n]  = 0.;
            tky[n]  = 0.;
            teg[n]  = 0.;
            txig[n] = 0.;
            tef[n]  = 0.;
            taf[n]  = 0.;
            tROmega[n]=0.;
            k_nc[n] = 0.;
            k_pi[n] = 0.;
            
        };
        
        t5[n]    = 0.;
        tkx5[n]  = 0.;
        tky5[n]  = 0.;
        teg5[n]  = 0.;
        txig5[n] = 0.;
        tef5[n]  = 0.;
        taf5[n]  = 0.;
        tROmega5[n]=0.;
        k_pi5[n] = 0.;
        
        
    }
    
    
    /*************************************s
     
     Momentum box via momaxis class
     
     *****/ //
    N[0]        = nxparam;  //Nkx -- number of points on x-direction
    N[1]        = nyparam;  //Nky -- number of points on y-dir
    N[2]        = 1 ;
    
    //Total number k-space points
    NTotal      = N[0]*N[1]*N[2];
    
    kmaxs0[0]       = 2.*pi/sqrt( 3. )/a0; // Momentum grid max on x-dir a.u.
    kmaxs0[1]       = 2.*pi/3./a0;         // Momentum grid max on y-dir a.u.
    kmaxs0[2]       = 0. ;
    
    
    //Controlling BZ gauge windows and limits
    gauge_ky_down  = 0.01;
    gauge_ky_up    = +0.70;
    
/*
    if (ellip==0.)
    {
        
        gauge_ky_down  = -0.16;//0.01;//
        gauge_ky_up    = 0.85;//+0.70;//
        
    }//*/

    
    kyShift[0]     = -0.35;
    kyShift[1]     = gauge_ky_down;
    
    
    
    //######################################
    //Opening output files!
    FILE *laserout, *mout, *fout, *engout;
    FILE *dipout, *conout, *curvaout, *gvelout;
    FILE *inh_out, *occup_out, *sparamout;
    FILE *simulation_out;
    
    
    
    
    if ( rank == MASTER )
    {
        
        
        laserout        = fopen( "outlaserdata.dat","w" ); // Laser field characteristics,time axis,
        mout            = fopen( "mgrid.dat","w" );
        
        engout          = fopen( "edispersion.dat","w" );
        dipout          = fopen( "dipoles.dat","w" );
        
        conout          = fopen( "connection.dat","w" );
        curvaout        = fopen( "curvature.dat","w" );
        
        gvelout         = fopen( "gvelocities.dat","w" );
        sparamout       = fopen( "setOfparameters.dat", "w" );
        
        
    }
    
    
    
    
    
    char full_rad[200] = "full_integrated_currents_rank";
    char dat0[]        = ".dat";
    char ctemp0[50];
    
    sprintf( ctemp0 , "_%.6d" , rank );

    strcat( full_rad, ctemp0 );
    strcat( full_rad, dat0 );
    
    
    
    simulation_out  = fopen( full_rad, "w" );
    
    
    

    
    /********************************************************
     
     Creating LASER PULSE and defining other parameters
     
     *****/
    


    laser fpulse(   npulses );             // Constructor of Laser Pulses

    I0      = E0*E0*3.5e16;             // Intensity W/cm^2


    


    
    
    //Collecting parameters
    double laserparam[]  = { I0, wfreq, ncycles, cep, ellipticity, relativephase, ltheta0, dt, tstart, offset, outset };

    

    
    //Initializing laser parameters ...
    initialize_laser_parameters( laserparam, env_name, fpulse );
    

    
    
    
    if ( rank==MASTER )
        {
    //Output of the laser field
    //fpulse.laser_outs_EA( laserout );
            
            
        cout << "\n*=============================*\nLaser Parameters:\n";
        cout << "E0           = "     << E0    << " a.u.  I0 = "  << I0   << endl;
        cout << "w0           = "     << wfreq << " a.u."     << endl;
        cout << "A0           = "     << E0/wfreq << " a.u."  << endl;
        cout << "NCycles      = "     << ncycles << " No. Opt. Cycles"<<endl;
        cout << "Up           = "     << E0*E0/wfreq/wfreq/4. << " a.u." << endl;
        cout << "ellipticity  = "     << ellip << endl;
        cout << "\nTime-step dt = "     << dt << " a.u." <<  endl;
        cout << "New-No. Of Total time Steps = " << fpulse.NewNt << "\n---\n";

        fprintf(sparamout,"\n%s","#Laser-Parameters,  E0,          I0,              w0,             Ncycles,           CEP,     No.-time-steps,  Time-step-dt "); 
        fprintf(sparamout,"\n          %e      %e      %e      %e      %e      %d      %e      \n\n", E0, I0, wfreq,  ncycles, cep, fpulse.NewNt, dt );
            
        }
    
    
    
    
    //##########################
    //Coherence and occupations
    setting_complex_memory( &rkutta_pi, NTotal );
    setting_double_memory(  &rkutta_nc, NTotal );
    
    cout << endl;
    for (itemp =0 ; itemp<NTotal; itemp++)
    {
        
        rkutta_pi[itemp]=0.;
        rkutta_nc[itemp]=0.;

    }
    
    for ( itemp=0; itemp<Ngrad; itemp++ )
    {        
        
        setting_complex_memory( &inter_rad[itemp], fpulse.NewNt+1 );
        setting_complex_memory( &intra_rad[itemp], fpulse.NewNt+1 );
        
        
        setting_double_memory( &group_rad[itemp], fpulse.NewNt+1 );
        setting_double_memory( &anomalous_rad[itemp], fpulse.NewNt+1 );
        
        for (jtemp=0; jtemp<fpulse.NewNt+1; jtemp++)
        {
            
            inter_rad[itemp][jtemp]=0.;
            intra_rad[itemp][jtemp]=0.;
            
            group_rad[itemp][jtemp]=0.;
            anomalous_rad[itemp][jtemp]=0.;
            
        }
        
    }
    

    setting_complex_memory( &ConductionOccup, fpulse.NewNt+1 );
    setting_complex_memory( &Coherence_Mom_Int, fpulse.NewNt+1 );

    
    for (jtemp=0; jtemp<fpulse.NewNt+1; jtemp++)
    {
        
        ConductionOccup[jtemp]=0.;
        Coherence_Mom_Int[jtemp]=0.;
    
    }
    
    

    
    
    /*************************************
     
        Momentum axis construction
     
     ***/
    
    momaxis g( N, kmaxs0, &a0 ); //momentum object
    //g.integral_method("Simpson");
    g.integral_method( IntionMethod );
    
    
    
    
    if (rank==MASTER)
    {
        
        //Output momentum
        g.mom_outputs( mout, 1, 1, 1 );
        
    
        cout << "\n\n==============================\n" ;
        cout << "Momentum integral method for:" << endl ;
        cout << IntionMethod << endl;
        /*cout <<"\nRectangular method or " ;
        cout <<"\nTrapezoidal method or " ;
        cout <<"\nSimpson Method\n===================\n" ;*///
        
        
        cout << "First Raw of weight matrix: ";
        ktemp=0;
        for ( int j =0; j < 1; j++ )
            for ( int i=0; i<N[0]; i++ )
                cout << "\nweight[ " <<i << ", " << j << "] = " << g.weight[ g.index( &i, &j, &ktemp ) ];

        cout << "\n\nLast Raw of weight matrix: ";
        ktemp=0;
        for ( int j =N[1]-1; j < N[1]; j++ )
            for ( int i=0; i<N[0]; i++ )
            cout << "\nweight[ " <<i << ", " << j << "] = " << g.weight[ g.index( &i, &j, &ktemp ) ];
        //*/
        
        
    }

    
    
    
    /*****************************************************
     
      Crystal structure object for the Haldane model
     
     *********/
    solidstructure cs( &g, &a0, &T2 );
    


    
    /*************************************
     
       Parameters of Haldane Model
     
     *********/
    
    cs.haldane_model_params( &t1, &t2, &phi0, &M0, gauge1 );
    
    
    cs.reg0     = treg0;
    cs.re_width = 0.2;
    cs.rflag    = trflag;
    
    
    //###############################
    //Computing Berry Features, Energy dispersion, and dipoles...
    cs.lattice_structure( );
    
    
    
    
    if ( rank == MASTER )
    {
        
    //###############################
    //Energy dispersion output
        cs.energy_vc_output(    engout,   1,  1,  1 );
        cs.dipole_cv_output(    dipout,   1,  1,  1 );
        cs.connection_c_output( conout,   1,  1,  1 );
        cs.curvature_c_output(  curvaout, 1,  1,  1 );
        cs.group_vel_vc_output( gvelout,  1,  1,  1 );
        
        fclose( engout );
        fclose( dipout );
        fclose( conout );
        fclose( curvaout );
        fclose( gvelout );
        
    }//*/
    
    
    
    if (rank == MASTER)
    {

        chernN0 = cs.ChernNumber( );
    
    //###############################
        cs.Set_Of_B_BGrad( &kx, &ky );
        cs.zBerryCurvaCV( );
        cs.group_velCV( );
    
    
    
        cout << "\n\n\n*=============================*\n";
        cout << "Haldane model structure";
        cout << "\na0           = "     << a0 << " a.u. (Lattice Const.)";
        cout << "\nt1           = "     << t1 << " a.u.";
        cout << "\nt2           = "     << t2;
        cout << "\nM0           = "     << M0;
        cout << "\nphi0         = "     << phi0 ;
        cout << "\ngauge1        = "     << cs.gauge;
        cout << "\ngauge2        = "     << gauge2;
        cout << "\n\n\n\n-----\nTesting Structure quantities at k =(1,1)\neg(1,1)     = " << cs.eg0 << " a.u.";
        cout << "\ngroup_c0     = "     << cs.group_c0[0];
        cout << "\ngroup_c1     = "     << cs.group_c0[1];
        cout << "\ngroup_v0     = "     << cs.group_v0[0];
        cout << "\ngroup_v1     = "     << cs.group_v0[1];
        cout << "\nBCon         = ("    << cs.chig0[0]/2. << "," << cs.chig0[1]/2. << ")";
        cout << "\nxdcv         = "     << cs.dip_cv0[0];
        cout << "\nydcv         = "     << cs.dip_cv0[1] << " ";
        cout << "\nBCurv        = "     << cs.zBcurva_c0 ;
        cout << "\nChernNo.     = "     << chernN0<< "\n";
        cout << "MinGap  	   = "      << cs.MinEGap << " a.u."; 
        cout << "\nMaxGap  	   = "      << cs.MaxEGap << " a.u. \n" ;
        
        
        
        cout << "\n---\nRegularization parameters\n";
        cout << "type of flag = "   << cs.rflag;
        cout << "\neps0 = "         << cs.reg0;
        cout << "\nwith-eps0 = "    << cs.re_width << endl << endl;
    
    
        
    //###############################
    //Fresh-dephasing
    //cs.T2   = period0;
        cout << "Dephasing, T2  = "     << cs.T2 << " a.u.\n\n\n";
    
    
        fprintf(sparamout,"\n\n%s","#Momentum Grid and Haldane Parameters, dkx (a.u.),  dky (a.u.),  Nx,  Ny, a0 (a.u.),  t1(a.u.),  t2 "); 
        fprintf(sparamout,"\n          %e      %e      %d      %d      %e      %e      %e      \n\n", g.dk[0], g.dk[1], g.N[0], g.N[1], a0, t1, t2 );



        fprintf(sparamout,"\n\n%s","#Haldane Parameters, phi0 (rad.),  M0 (a.u.),   Min-eg(a.u.),        Max-eg,          Chern No.,        T2 (a.u.),    gauge1 ");


        fprintf(sparamout,"\n          %e      %e      %e      %e      %e      %e      %d      \n\n", phi0, M0, cs.MinEGap, cs.MaxEGap, chernN0, cs.T2, cs.gauge );


        
        fprintf(sparamout,"\n\n%s","#Parameters,   No.cores,  ellipticity,        regularization,      reg.-width,     gauge_ky_down,     gauge_ky_up,     gauge2 ");
        
        fprintf(sparamout,"\n              %d          %e         %e        %e      %e      %e      %d      \n\n", Nprocessors, ellip, cs.reg0, cs.re_width, gauge_ky_down, gauge_ky_up, gauge2 );
        
        
        fflush(sparamout );
        
        
        fclose(sparamout);
        
        
        

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
        cout << "\nNtotal=Nx*Ny*Nz= "             << NTotal;
        
        
        cout <<"\n\nMinima momentum points";
        cout << "\nkxmin        = "            << g.k[0][0];
        cout << "\nkymin        = "            << g.k[1][0];
        cout << "\nkzmin        = "            << g.k[2][0];
        cout << "\n\nky         = "            << g.k[1][0];
        cout << "\n/*************************/\n";
        
        
        
        //###############################
        //Saving or Output of laser
        for ( n = 0; n < fpulse.NewNt; n++)
           {
        
               
               
            at = fpulse.atmin + dt*n;
        
        
            fpulse.a_ef0 = fpulse.elaser( &at );
            fpulse.a_af0 = fpulse.avlaser( &at );
        
               
        
            fprintf( laserout,"%e    %e    %e    %e    %e\n", at, real(fpulse.a_ef0), imag(fpulse.a_ef0 ), real(fpulse.a_af0), imag(fpulse.a_af0) );


               
               
            }
        
        
        
        fflush( laserout );
        fclose( laserout );

        
        
    }
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    
    //################################################################################
    //Splitting Ny direction integration on different cores by MPI 1D scheme
    jstep   =  ceil ( double(g.N[1])/double(Nprocessors) );
    jstart  =  rank * jstep;
    jend    =  ( rank + 1 ) * jstep;
        
    
    if ( rank+1==Nprocessors )
        jend = g.N[1];
    //########################################
    
    
    
    
    
    cout << "\nrank = " << rank << "  uses j = " << jstart << ";  to  j = " << jend << "  with Nprocesses = " << Nprocessors <<endl << endl;
    
    
    if( Nprocessors > g.N[1] )
    {
    
         cout << "\n************\nNprocesses have to be less than the numbers of ky points, please, change it accordling! \n\n" ;
         exit(0);
        
         
    }
    
    
    if(rank==MASTER)
    {
        

        
        cout.precision( 5 );
        
        
        
        
        cout << "\n\n\n//############################################//" ;
        cout <<"\n#*********==============**********#\n ";
        cout << "Beginning of Momentum/Time integration loops\n (Nx,Ny,Nt)";
        cout << " = ( "<< g.N[0] << " , " << g.N[1] <<" , " <<fpulse.NewNt << " ) " << scientific;
        cout << "\nMomentum Steps and Time Step";
        cout << "\n ( dkx, dky, dt )     =   (" << g.dk[0] << ",  " << scientific;
        cout << g.dk[1] << ",  " << dt << "  )  a.u." << scientific << endl;
        cout <<"############################################"<<endl;
        
        
        cout << "\n\nktime     " << "       ";
        cout << "         nc      " <<  "        pi(a.u.)    ";
        cout << "   gauge1      id-rank \n" << scientific;
        
    }
    
    

    
    //#####################################
    //#####################################
    //Time integration loop
    for( ktime = 0; ktime<fpulse.NewNt; ktime++ )
    {


        
        
        
        
        //#############################################
        // Defining set of times for Runge-Kutta-4th
        // at tn, tn + dt/2., and, tn + dt, where tn =  t0 + n*dt, n=0,1,...,Nx/y
        t[0]    = fpulse.atmin + dt*double(ktime);
        t[1]    = t[0]  + dt/2.;
        t[2]    = t[1];
        t[3]    = t[0]  + dt;
        
        //Runge-Kutta 5th, RK5, time-steps
        t5[0]    = t[0] ;
        t5[1]    = t[0] + dt/4. ;
        t5[2]    = t[0] + 3.*dt/8. ;
        t5[3]    = t[0] + 12.*dt/13. ;
        t5[4]    = t[3] ;
        t5[5]    = t[1] ;
        
        
        
        //#####################################
        // Evaluation of SBEs variables
        // At tn, tn + dt/2., and, tn + dt, where tn = t0 + n*dt
        // here, t0 is the initial time-grid point
        // n is a time-index, corresponing to n -> ktime in the code
        // dt, time-step
        for( ltemp = 0; ltemp < Nrk5; ltemp++ )
        {
            
            if (ltemp<Nrk4)
            {
                
                tef[ ltemp ] = fpulse.elaser(  &t[ ltemp ] );
                taf[ ltemp ] = fpulse.avlaser( &t[ ltemp ] );
                
            }

            tef5[ ltemp ] = fpulse.elaser(  &t5[ ltemp ] );
            taf5[ ltemp ] = fpulse.avlaser( &t5[ ltemp ] );
            
            
        }//End loop of evaluation of KR4 vars.
        //#####################################
        
        
        
        
        
        
    //###############################
    //Momentum and time integration
    for( jtemp = jstart; jtemp < jend; jtemp++ )
    {            
        

        //Checking wavefunction gauge region for ky axis
        cs.gauge    = gauge1;
        if ( g.k[1][jtemp] > gauge_ky_down & g.k[1][jtemp] <=gauge_ky_up )
            cs.gauge = gauge2;
        

        
        //EVALUATION OF RK5 and RK4 VARIABLES FOR SBEs along ky
        for ( ltemp = 0; ltemp < Nrk5; ltemp++ )
        {
            
            
            //Evaluating RK4 for occupations
            if (ltemp < Nrk4)
            {
                
                
                //Evaluationg RK4 variables for ky at time, t[ltemp]
                tky[ ltemp ]      = g.k[1][jtemp] + imag( taf[ ltemp ] );
                
                
                //Shifting BZ, if is needed along  ky
                if ( g.k[1][jtemp] >= kyShift[0] & g.k[1][jtemp] <= kyShift[1] )//  & ellip !=0. )
                    tky[ ltemp ]+= g.DeltaKp[1]  ;
                
                
                
                
            }//End evaluation RK4
            
            
            
            //Evaluationg RK5 variables for ky at time, t[ltemp]
            tky5[ ltemp ]      = g.k[1][jtemp] + imag( taf5[ ltemp ] );
            
            //Shifting BZ if it is needed along ky
            if ( g.k[1][jtemp] >= kyShift[0] & g.k[1][jtemp] <= kyShift[1] )//& ellip !=0. )
                tky5[ ltemp ]+= g.DeltaKp[1]  ;
            
            
        }//END OF THE EVALUATION OF RK5 and RK4 VARIABLES FOR SBEs along ky
            
        
        
        
        //###############################        
        //Momentum kx direction Loop
        for( itemp = 0; itemp < g.N[0] ; itemp++ )
        {
            


            //###############################
            //Evaluating variables for kx at time tn
            tkx[ 0 ]      = g.k[0][ itemp ] + real( taf[ 0 ] );


            //###############################        
            //Shifting BZ, if is needed along kx  at tn
            if ( g.k[1][jtemp] >= kyShift[0] & g.k[1][jtemp] <= kyShift[1] )// & ellip !=0. )
                tkx[ 0 ]+= g.DeltaKp[0];
                
            
            

            //#################################################
            // Computing set of dipoles, Berry connections,
            // curvatures and energy dispersions
            cs.Set_Of_B_BGrad( &tkx[ 0 ], &tky[ 0 ] );
            
            
            //Group velocity evaluations and Rabbi Freq., E(t)\cdot d_cv(t)
            cs.group_velCV( );
            tROmega[0]    = cs.dip_cv0[0] * real( tef[ 0 ] ) + cs.dip_cv0[1] * imag( tef[ 0 ] ) ;
            
            
            
            
            
            
            
            
            //##########################################################
            //MOMENTUM INTEGRALS FOR INTER AND INTRA BAND POL. AND CUR.
            //###############################
            
            
            //##########################################################
             //CALCULATION OF INTER-BAND POLARIZATION
            //###########################################
            inter_rad[0][ktime]+= 2.*real( conj( cs.dip_cv0[0] )*rkutta_pi[ g.index( &itemp, &jtemp, &ktemp) ] )*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ] ;
            
            
            
            inter_rad[1][ktime]+= 2.*real( conj( cs.dip_cv0[1] )*rkutta_pi[ g.index( &itemp, &jtemp, &ktemp) ] )*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            //END OF INTER-BAND CALCULATIONS
            //############################################
            
            
            
            
            
            
            
            
            //######################################
            //CALCULAITON OF INTRABAND CURRENTS
            //######################################

            //#####################################
            // Evaluation of the anomalous velocity
            //  at time tn
            cs.anomalous_velCV( tef[ 0 ], cs.zBcurva_c0 );
            

            
            
            
            //################################################
            //Valence band occupations  at time tn
            nv              = 1. -  rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ] ;
            //################################
            
            
            
            
            
            //*************************************************/
            //Preparing variables for intra-band current
            
            //Group velocities * occupations at time tn
            group_velocity_v[0] =  cs.group_v0[0]  * nv * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            group_velocity_v[1] =  cs.group_v0[1]  * nv * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            

            group_velocity_c[0] =  cs.group_c0[0]  * rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ]  * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            group_velocity_c[1] =  cs.group_c0[1]  * rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ]  * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            
            
            //Anomalous velocities * occupations at time tn
            anomalous_vel_v[0] = cs.xanomalous_v0 * nv * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            anomalous_vel_v[1] = cs.yanomalous_v0 * nv * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            
            
            anomalous_vel_c[0] = cs.xanomalous_c0 * rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ] * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            anomalous_vel_c[1] = cs.yanomalous_c0 * rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ] * g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            
            
            
            //Classical velocities  at time tn
            classical_velocity_v[0] = group_velocity_v[0] + anomalous_vel_v[0];
            classical_velocity_v[1] = group_velocity_v[1] + anomalous_vel_v[1];
            
            classical_velocity_c[0] = group_velocity_c[0] + anomalous_vel_c[0];
            classical_velocity_c[1] = group_velocity_c[1] + anomalous_vel_c[1];
            
            
            //Radiation by group velocity
            group_rad[0][ktime]+=  group_velocity_v[0] + group_velocity_c[0];
            group_rad[1][ktime]+=  group_velocity_v[1] + group_velocity_c[1];
            
            //Radiation by anomalous velocity
            anomalous_rad[0][ktime]+= anomalous_vel_v[0] + anomalous_vel_c[0];
            anomalous_rad[1][ktime]+= anomalous_vel_v[1] + anomalous_vel_c[1];
            //********************************************************************/
            
            
            
            
            
            //################################
            //INTRABAND CURRENTs
            intra_rad[0][ktime]+= classical_velocity_v[0] + classical_velocity_c[0] ;
            intra_rad[1][ktime]+= classical_velocity_v[1] + classical_velocity_c[1] ;
            //################################
            
            
            
            
            //##############################################
            //Evaluation of occupations and coherence
            ConductionOccup[ktime]+=   rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ]*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            
            Coherence_Mom_Int[ktime]+= rkutta_pi[ g.index( &itemp, &jtemp, &ktemp ) ]*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            //##############################################
            
            
            
            
            //###################VERY-IMPORANT-CODE-PART#################\\
            //###################################################
            // Time & k-space PROPAGATION A.C. SCHEME OF SBEs
            //.... Runge-Kutta propagator, RK4 and RK5 ....
            // FINDING:
            // (1) coherence  pi(kx,ky) at tn + dt and,
            // (2) occupation nc(kx,ky) at tn + dt
            //##################################################
            
            
            //EVALUATION OF RK5 and RK4 VARIABLES FOR SBEs
            for ( ltemp = 0; ltemp < Nrk5; ltemp++ )
            {
                
                
                //Evaluating RK4 for occupations
                if (ltemp > 0 & ltemp < Nrk4)
                {
                    
                    //Evaluationg RK4 variables for kx at the time t[ltemp]
                    tkx[ ltemp ]      = g.k[0][itemp] + real( taf[ ltemp ] );
                    
                    
                    //Shifting BZ, if is needed along kx
                    if ( g.k[1][jtemp] >= kyShift[0] & g.k[1][jtemp] <= kyShift[1]  )//& ellip !=0. )
                        tkx[ ltemp ]+= g.DeltaKp[0]  ;

                    
                    
                    //Computing set of dipoles, Berry connections, curvatures and energy dispersions
                    cs.Set_Of_B_BGrad( &tkx[ ltemp ], &tky[ ltemp ] );
                    
                    
                    tROmega[ltemp]    = cs.dip_cv0[0] * real( tef[ ltemp ] ) + cs.dip_cv0[1] * imag( tef[ ltemp ] ) ;
                    
                    
                }//End evaluation RK4
                
                
                
                //Evaluationg RK5 variables for kx at time, t[ltemp]
                tkx5[ ltemp ]      = g.k[0][itemp] + real( taf5[ ltemp ] );
                
                
                //Shifting BZ if is needed along kx
                if ( g.k[1][jtemp] >= kyShift[0] & g.k[1][jtemp] <= kyShift[1]  )//& ellip !=0. )
                     tkx5[ ltemp ]+= g.DeltaKp[0]  ;
                
                
                //Computing set of dipoles, Berry connections, curvatures and energy dispersions
                cs.Set_Of_B_BGrad( &tkx5[ ltemp ], &tky5[ ltemp ] );
                
                
                //Energy dispersion difference between cond. and val. bands
                teg5[ ltemp ]      = cs.eg0;
                txig5[ ltemp ]     = cs.chig0[0]*real( tef5[ltemp] )+cs.chig0[1]*imag( tef5[ltemp] );
                tROmega5[ltemp]    = cs.dip_cv0[0]*real( tef5[ltemp] )+cs.dip_cv0[1]*imag(tef5[ltemp] );
                
                
            
            }//END OF EVALUATION OF RK5 and RK4 VARIABLES
                
        
            
            
            
            //###############################################################
            //Coherence pi(kx,ky,t) of the occupations nc(kx,ky,t) via RK5
            //###########################################
            
            temp_pi5[0] = rkutta_pi[g.index( &itemp, &jtemp, &ktemp ) ];
            der_pi( &teg5[0], &txig5[0], &cs.T2, &tROmega5[0], &temp_pi5[0], &k_pi5[0] );
            
            
            
            temp_pi5[1] = temp_pi5[0] + dt*k_pi5[0]/4.;
            der_pi( &teg5[1], &txig5[1], &cs.T2, &tROmega5[1], &temp_pi5[1], &k_pi5[1] );
            
            
            
            temp_pi5[2] = temp_pi5[0] + dt*(k_pi5[0]*3./32.+ k_pi5[1]*9./32.);
            der_pi( &teg5[2], &txig5[2], &cs.T2, &tROmega5[2], &temp_pi5[2], &k_pi5[2] );
            
            
            
            temp_pi5[3]    = temp_pi5[0] + (k_pi5[0]*1932./2197. - k_pi5[1]*7200./2197.
                            + k_pi5[2]*7296./2197.)*dt;
            der_pi( &teg5[3], &txig5[3], &cs.T2, &tROmega5[3], &temp_pi5[3], &k_pi5[3] );
            
            
            
            temp_pi5[4]    = temp_pi5[0] +  (k_pi5[0]*439./216.  -  k_pi5[1]*8.
                            + k_pi5[2]*3680./513. - k_pi5[3]*845./4104.)*dt;
            der_pi( &teg5[4], &txig5[4], &cs.T2, &tROmega5[4], &temp_pi5[4], &k_pi5[4] );
            

            
            temp_pi5[5]    = temp_pi5[0] +(- k_pi5[0]*8./27. + k_pi5[1]*2.
                            - k_pi5[2]*3544./2565.  +  k_pi5[3]*1859./4104. - k_pi5[4]*11./40.)*dt;
            der_pi( &teg5[5], &txig5[5], &cs.T2, &tROmega5[5], &temp_pi5[5], &k_pi5[5] );
            
            
            
            
            
            //########################################################
            //RK5 solution for the coherence pi(kx,ky) at time:
            //tn + n*dt
            rkutta_pi_ktimeplus = rkutta_pi[ g.index( &itemp, &jtemp, &ktemp )  ]
                                +  (16./135.*k_pi5[0]      +  6656./12825.*k_pi5[2]
                                +  k_pi5[3]*28561./56430. -  k_pi5[4]*9./50.   + k_pi5[5]*2./55.)*dt;

            //#################################
            
            
            
            
            
            
            //####################################################
            //Occupation via RK4
            //###############################
            
            der_nc( &tROmega[0], &rkutta_pi[ g.index( &itemp, &jtemp, &ktemp ) ], &k_nc[0] );
            
            _half_pi = ( rkutta_pi_ktimeplus + rkutta_pi[ g.index( &itemp, &jtemp, &ktemp ) ] )/2.;
            der_nc( &tROmega[1], &_half_pi, &k_nc[1] );
            
            der_nc( &tROmega[2], &_half_pi, &k_nc[2] );
            
            der_nc( &tROmega[3], &rkutta_pi_ktimeplus, &k_nc[3] );
            
            
            
            
            //########################################################
            //RK4 solution for the conduction band occupation;
            // nc(kx,ky), at the time: tn + n*dt
            rkutta_nc_ktimeplus = rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ] + dt/6.*( k_nc[0] + ( k_nc[1] + k_nc[2] )*2. + k_nc[3] );
            //###################################
            
            
            
            
            //saving tn + dt pi(kx,ky) and nc(kx,ky)
            rkutta_pi[ g.index( &itemp, &jtemp, &ktemp ) ] = rkutta_pi_ktimeplus;
            rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ] = rkutta_nc_ktimeplus;
        
        
            
            
            
        }//End of Loop on x-direction
        
        
        

        
    }//End of Loop on y-direction
    
        
        
        
        
    
    //OutPut for diagnostic
    if( ktime%((fpulse.NewNt-1)/200)==0 )
    {
            
        if ( g.k[1][jtemp-1] >= kyShift[0] & g.k[1][jtemp-1] <= kyShift[1] )
        {
                
            printf( "\n%.6d / %.7d        %.3e        %.3e        %.2d        %.6d     %s  %.3f %s %.3f %s", ktime, fpulse.NewNt-1, real(ConductionOccup[ktime]*g.dV  ), real(Coherence_Mom_Int[ktime])*g.dV, cs.gauge, rank, "Shifted BZ (",tkx[ 0 ],",",tky[ 0 ],")" );
            
        }
        else
        {
                
            printf( "\n%.6d / %.7d        %.3e        %.3e        %.2d        %.6d     %s  %.3f %s %.3f %s", ktime, fpulse.NewNt-1, real(ConductionOccup[ktime]*g.dV  ), real(Coherence_Mom_Int[ktime])*g.dV, cs.gauge, rank, "Non-Shift BZ (",tkx[ 0 ],",",tky[ 0 ],")" );
                
        }
            
    }//END OF DIAGNOSTIC

        
        
 
        
}//End of Time integration loop
    
    
    
    //#######################################
    //Output inter, intra band currents,
    //occupation and coherence and vs time
    
    if (rank == -1) 
    {
        
        fprintf(simulation_out,"%s    %s                   %s                    %s                 %s                   %s                    %s                       %s                    %s                 %s                   %s                    %s\n ", "#","#time(a.u.)","#Jx_er(a.u.)","#Jy_er(a.u.)","#Jx_ra(a.u.)","#Jy_ra(a.u.)", "#nc(a.u.)","#Coherent(a.u.)", "#vgx" , "#vgy", "#vax", "#vay" );
        
    }
    
    
    ktime=0;
    fprintf(simulation_out,"%.16e       %.16e         %.16e         %.16e         %.16e         %.16e         %.16e         %.16e         %.16e         %.16e         %.16e\n ", fpulse.atmin, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV, group_rad[0][ktime]*g.dV, group_rad[1][ktime]*g.dV, anomalous_rad[0][ktime]*g.dV, anomalous_rad[1][ktime]*g.dV  );
    
    
    
    for( ktime = 1; ktime<fpulse.NewNt-1; ktime++ )
    {
        
        
        at  = fpulse.atmin + dt*ktime;
        
        
        
        dxdt = real( inter_rad[0][ktime+1] - inter_rad[0][ktime-1] )*g.dV/dt/2.;
        
        
        
        
        dydt = real( inter_rad[1][ktime+1] - inter_rad[1][ktime-1] )*g.dV/dt/2.;
        
        
        
        
        fprintf(simulation_out,"%.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e       %.16e\n ", at, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV, group_rad[0][ktime]*g.dV, group_rad[1][ktime]*g.dV, anomalous_rad[0][ktime]*g.dV, anomalous_rad[1][ktime]*g.dV );
        
        
        
    }
    
    
    
    ktime   = fpulse.NewNt-2;
    at      = fpulse.atmin + dt*(ktime+1);
    fprintf(simulation_out,"%.16e       %.16e        %.16e       %.16e       %.16e       %.16e      %.16e       %.16e       %.16e       %.16e       %.16e\n", at, dxdt , dydt  , real(intra_rad[0][ktime])*g.dV, real(intra_rad[1][ktime])*g.dV, real(ConductionOccup[ktime])*g.dV, real(Coherence_Mom_Int[ktime])*g.dV, group_rad[0][ktime]*g.dV, group_rad[1][ktime]*g.dV, anomalous_rad[0][ktime]*g.dV, anomalous_rad[1][ktime]*g.dV );
    
    
    //cout <<"\nNo. of total time steps = " << fpulse.NewNt << endl;
    
    fflush(simulation_out);
    
    
    //fclose(  fout     );
    if (rank==MASTER)    
    {       
        fclose(  conout   );
        fclose(  dipout   );
        fclose(  curvaout );
        fclose(  engout   );
    }   
    
    
    for ( itemp=0; itemp < Ngrad; itemp++)
    {
        free( inter_rad[itemp] );
        free( intra_rad[itemp] );
        
        free( group_rad[itemp] );
        free( anomalous_rad[itemp] );

        
    }
    
    free( ConductionOccup );
    free( Coherence_Mom_Int );
    free( rkutta_nc );
    free( rkutta_pi );
    
    
    fclose(simulation_out);
    
    
    cout << "\n***********************\nEND OF THE PROGRAM with No. of cores = " << Nprocessors << "  uses j = " << jstart << ";  to  j = " << jend << "  with rank = "  << rank << "\n\n" <<  endl;
    cout << "=====================\n\n";
    
    
    MPI_Finalize( )	;	
	//exit(0); 


    
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


void der_pi( double const *eg, double const *xig_ef, double const *T2, complex const *ROmega, complex const *cohPi0, complex *cpi)
{
    
    *cpi = -I*( ( *eg + *xig_ef - I/(*T2) ) *( *cohPi0 ) + (*ROmega) );
    
}

void der_nc( complex const *ROmega, complex const *cpi, double *nc )
{
    
    *nc = -2.*imag( conj( *ROmega )*( *cpi ) );
    
}


complex der_pi_RK5( double const *eg, double const *xig_ef, double const *T2, complex const *ROmega, complex const *cohPi0)
{
    
    return -I*( ( *eg + *xig_ef - I/(*T2) ) *( *cohPi0 ) + (*ROmega) );
    
}


complex jaco_pi_RK5( double const *eg, double const *xig_ef, double const *T2 )
{
    
    return -I * ( *eg + *xig_ef - I/(*T2) )  ;
    
}
//*/

double der_nc_RK5( complex const *ROmega, complex const *cpi )
{
    
    return -2.*imag( conj( *ROmega )*( *cpi ) );
    
}


void setting_complex_memory(complex **_TempPointer, long long int _Ntime)
{
    *_TempPointer      = (complex*)malloc( _Ntime*sizeof(complex) );
    memset( *_TempPointer,     0,    sizeof(complex)*_Ntime );
    /*for (int i=0; i<_Ntime; i++)
    {
        _TempPointer[i]=complex(0.,0.);
        //_TempPointer++;
    }//*/
}


void setting_double_memory(double **_TempPointer, long long int _Ntime)
{
    *_TempPointer      = (double*)malloc( _Ntime*sizeof(double) );
    memset( *_TempPointer,     0,    sizeof(double)*_Ntime );
    
   /* for (int i=0; i<_Ntime; i++)
    {
        _TempPointer[i]=0.;
        //_TempPointer++;
    }//*/
    
}



//MPI Operation complex sum of arrays
void My_MPI_CLX_SUM( complex *in, complex *inout, int *len, MPI_Datatype *dptr )
{
	for (int itemp = 0; itemp < *len; itemp++ )
	{
		*inout+= *in;
		in++;
		inout++;
	}
}


