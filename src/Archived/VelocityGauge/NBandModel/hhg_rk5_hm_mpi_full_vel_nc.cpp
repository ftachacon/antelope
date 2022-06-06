/**
 * @file hhg_rk5_hm_mpi_full_vel_nc.cpp
 * @author Dasol Kim
 * @author Alexis Agustín  Chacón Salazar
 * @brief 
 * @version 0.1
 * @date 2020-12-06
 * 
 * @copyright Copyright (c) 2020
 * 
 */


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

#include <unistd.h>

// My Own Headers
#include "constant.h"
#include "timegrid.h"
#include "timeobject.h"
#include "laser.h"
#include "momaxis.h"
//#include "operators.h"
//#include "solidstructure.h"
//#include "observables.h"
//#include "sbeqs.h"
#include "sbesVG.h"
#include "utility.h"

#define MASTER 0    /* task id of master task or node */
#define energy_factor 27.2

//int const Nrk4 = 4; //Runge-Kutta 4th coheff.
//int const Nrk5 = 6; //Runge-Kutta 5th coheff.

void initialize_laser_parameters( double eparam[],  string ename, laser &lobject );
void setting_complex_memory(complex **_TempPointer, long long int _Ntime);
void setting_double_memory(double **_TempPointer, long long int _Ntime );
//void der_pi( const double *eg, const double *xig_ef, const double *T2, const  complex *ROmega, const  complex *cohPi0, const double *nc, complex *cpi);
void der_pi( const double  *eg, const double *avg, const double *T2, const complex *apcv, const complex *cohPi0, const double *_nc, complex *cpi);
void der_nc( const complex *ROmega, const complex *cpi, double *nc );
//void der_nc2( const complex *ROmega, const complex *cpi, const double *T1, const double *tnc,  double *nc );
void der_nc2( const complex *apcv, const complex *cpi, const double *T1, const double *tnc,  double *nc );

complex der_pi_RK5( double const *eg, double const *xig_ef, double const *T2, complex const *ROmega, complex const *cohPi0);
complex jaco_pi_RK5( double const *eg, double const *xig_ef, double const *T2 );

double der_nc_RK5( complex const *ROmega, complex const *cpi );

void My_MPI_CLX_SUM( complex *in, complex *inout, int *len, MPI_Datatype *dptr );



//void My_2D_Domain_Reconstruction(int rank, int block, int ist, int iend, long long int *xshift, int *rmin_int_rank, int *imin_rank, int *jmin_rank, int *imax_rank, int *jmax_rank, int *xi, int *yj, long long int *iglobal, double *rmin_rank, double *rmax_rank, double *rmin_decimal_rank);
//void My_2D_Domain_Decomposition(int rank, int Nprocessors, int *block, int *ist, int *iend );


//int Nband     = 2;
//int Nvb    = 1;

int ncyparam    = 3;
int nxparam     = 1;
int nyparam     = 1;
int tgauge2     = 0.;
int  trflag     = 0;
int  fbz_shift  = 0;


double iparam   = 1.;
double jparam   = 0.; 
double kparam   = 0.; 
double dtparam  = 1.;
double dephasing= 1000.;
double ellip    = 0.;
double treg0    = 0.;

double ky_shift_down  = 0.;
double gbox_ky_down = 0.;
double gbox_ky_up   = 0.;
double ksfactor = 1;
double flag_occ_dephasing = 0.;
double flag_t2 = 0.;
int diagnostic = 0;
int shotNumber   = 0;

//################################
int main( int argc, char *argv[] )
{
    MPI_Op MPI_CLX_SUM;
    
    
	//Fundamental initialization of MPI
	MPI_Init(&argc, &argv);	
	
    
	//Creating our MPI-operation to sum up complex arrays
	MPI_Op_create((MPI_User_function *) My_MPI_CLX_SUM, true, &MPI_CLX_SUM);	
	
    //#############################
	//Getting rank and size variables
	int rank=0, size=1;
	int Nprocessors=1,NumberOfDoubles=1,NumberOfComplex=1;


	MPI_Comm_rank( MPI_COMM_WORLD, &rank );	//Getting the local id or rank process
	MPI_Comm_size( MPI_COMM_WORLD, &size );	//Getting the size or number of proccesses or nodes
	MPI_Barrier( MPI_COMM_WORLD );	


    //String Type of variables, envelope and integration rule
    string env_name     = "gauss";          //Name envel: "rect", "sin2", "gauss" or "rsin2"
    string IntionMethod = "Trapz";//       // "Trapz" or "Simpson" or "NoneTrapz"/"NoneSimpson"


	if (rank==MASTER)
	{
		cout << "\n\n\n#*============================================================*#\n";
		cout << "       High order harmonic generation in the \n       Topological Haldane Model (HM), Chern Insulator " << endl;
		cout << "       via the numerical integration of the \n       Semiconductor Bloch Equation (SBEs) \n          by Alexis Chacon June 27, 2019....\n";
		cout << "#*=========================================================*#\n\n\n";
	}

    sbesVG *sbe;
    sbe = new sbesVG();

    // txt  - input file with text file, energy relation dispersion is text file
    // bin - input file with binary file, mometum matrix should be binary file

    // Reading files, allocating memory
    if (!sbe->ReadInputParam("input.txt") 
    || !sbe->ReadPMatrix("pmatx.bin" , 0) || !sbe->ReadPMatrix("pmaty.bin", 1)
    || !sbe->ReadEDispersion("edispersion.txt")) // add the binary reader
    {
        cout << "Errors in reading data files\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    double dV = sbe->kgrid->dV;
    int Nband = sbe->Nband;
    int Nvb = sbe->Nvb;
    
    //#############################
	//Number of processes
	Nprocessors		= size;
	NumberOfDoubles = 1;
	NumberOfComplex = 1;
	

    //#############################
	MPI_Barrier( MPI_COMM_WORLD );		
	// -- Printing on terminal the PARAMETERs -- //	
	/*if (rank==MASTER)
	{	cout << "\n-------------\nPARAMETERS\n";
		cout << "\niparam = "		  	<< iparam ;
		cout << "\njparam = "         	<< jparam ;
		cout << "\nNxparam = "         	<< nxparam ;
		cout << "\nNyparam = "         	<< nyparam ;        
		cout << "\n---\nTotal of Processors= " << Nprocessors << "\n";
	}*/
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
    int gauge2   = tgauge2 ;
    
    
     //Controlling momentum gauge variation window
    double gauge_ky_down  =  0.;//-0.16;
    double gauge_ky_up    =  0.;//0.90
    double kmaxs0[Ndim]   = {0.,0.,0.};
    double kyShift[Ngrad] = {0.,0.}; //{-0.38,gauge_ky_down};

    
    
    
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
    //complex *ConductionOccup;
    //complex *Coherence_Mom_Int;
    complex **density_matrix_integrated;



    double at=0, nv=0.;

    complex tgv_v[Ngrad] = {0.,0.};
    complex tgv_c[Ngrad] = {0.,0.};
    complex tdip_cv[Ngrad]={0.,0.};
    
    
    //Occupation and Coherence variables
    complex *rkutta_pi;
    complex rkutta_pi_ktimeplus=0.;

    double *rkutta_nc;
    double rkutta_nc_ktimeplus=0.;

    // Number of k-grid for each process
    // and their displacement in global buffer
    int *blockcounts;
    int *displs;

    int kindex;

    complex tempCurrent;
    
    //######################################
    //Opening output files!
    FILE *laserout, *mout, *fout, *engout;
    FILE *dipout, *conout, *curvaout, *gvelout;
    FILE *inh_out, *occup_out, *sparamout, *lparamout;
    FILE *simulation_out;
    FILE *density_out, *densitytime_out, *polarization_out;
    FILE *intraj_out, *interj_out;
    
    
    if ( rank == MASTER )
    {
        
        
        laserout         = fopen( "outlaserdata.dat","w" ); // Laser field characteristics,time axis,
        //mout             = fopen( "mgrid.dat","w" );
        
        //sparamout        = fopen( "setOfparameters.dat", "w" );
        lparamout        = fopen( "laserParameters.dat", "w" );

        interj_out       = fopen( "interband_dipole_full_evol.dat", "w");
        intraj_out       = fopen( "intraband_current_full_evol.dat", "w");
        occup_out        = fopen( "occupation__full__evol.dat", "w");
        
    }
    
    

    if (rank == MASTER)
    {
        //Output of the laser field
        //fpulse.laser_outs_EA( laserout );

        /*cout << "\n*=============================*\nLaser Parameters:\n";
        cout << "E0           = " << E0 << " a.u.  I0 = " << I0 << endl;
        cout << "w0           = " << wfreq << " a.u." << endl;
        cout << "A0           = " << E0 / wfreq << " a.u." << endl;
        cout << "NCycles      = " << ncycles << " No. Opt. Cycles" << endl;
        cout << "Up           = " << E0 * E0 / wfreq / wfreq / 4. << " a.u." << endl;
        cout << "ellipticity  = " << ellip << endl;*/
        cout << "\nTime-step dt = " << sbe->fpulses->dt << " a.u." << endl;
        cout << "New-No. Of Total time Steps = " << sbe->fpulses->Nt << "\n---\n";

        sbe->fpulses->Print_LaserInfo();
    }

    //##########################
    //Coherence and occupations
    
    for ( itemp=0; itemp<Ngrad; itemp++ )
    {        
        
        setting_complex_memory( &inter_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &intra_rad[itemp], sbe->fpulses->Nt );
        
        
        
        for (jtemp=0; jtemp<sbe->fpulses->Nt; jtemp++)
        {
            
            inter_rad[itemp][jtemp]=0.;
            intra_rad[itemp][jtemp]=0.;
            
            
        }
        
    }

    density_matrix_integrated = Create2D<complex>(sbe->fpulses->Nt, Nband*Nband);
    
    
    /*************************************
     
        Momentum axis construction
     
     ***/

    
    if (rank == MASTER)
    {
        
    //###############################
    //Fresh-dephasing
    //cs.T2   = period0;
       /* cout << "Dephasing, T2  = "     << cs.T2 << " a.u.\n\n\n";
    
        fprintf(sparamout,"\n\n%s","#Time grid, No.-time-steps,  Time-step-dt,     No.-pulses     ---,     ---,     ---,     --- ");
        fprintf(sparamout,"\n          %d      %e      %d      %e      %e      %e      %e      \n\n", fpulse.NewNt, dt, fpulse.Npulses, 1.0, 1.0, 1.0, 1.0 );
    
        fprintf(sparamout,"\n\n%s","#Momentum mesh & H.M. params., dkx (a.u.),  dky (a.u.),  Nx,  Ny, a0 (a.u.),  t1(a.u.),  t2 ");
        fprintf(sparamout,"\n          %e      %e      %d      %d      %e      %e      %e      \n\n", g.dk[0], g.dk[1], g.N[0], g.N[1], a0, t1, t2 );



        fprintf(sparamout,"\n\n%s","#H.M. params.:, phi0 (rad.),  M0 (a.u.),   Min-eg(a.u.),        Max-eg,          Chern No.,        T2 (a.u.),    gauge1 ");


        fprintf(sparamout,"\n          %e      %e      %e      %e      %e      %e      %d      \n\n", phi0, M0, cs.MinEGap, cs.MaxEGap, chernN0, cs.T2, cs.gauge );


        
        fprintf(sparamout,"\n\n%s","#Params.,   No.cores,  ellipticity,        regularization,      reg.-width,     gauge_ky_down,     gauge_ky_up,     gauge2 ");
        
        fprintf(sparamout,"\n              %d          %e         %e        %e      %e      %e      %d      \n\n", Nprocessors, ellip, cs.reg0, cs.re_width, gauge_ky_down, gauge_ky_up, gauge2 );
        

        fprintf(sparamout,"\n\n%s","#Params.,   BZ-shift-Switch,     ky-shift-0,       ky-shift-f,      BZ.-kxmax,      BZ.-kymax,      snapshot rate,     --- ");
        
        fprintf(sparamout,"\n              %d          %e         %e        %e      %e      %d      %e      \n\n", fbz_shift, kyShift[0], kyShift[1], kmaxs0[0], kmaxs0[1], shotNumber, 0. );

        
        fflush(sparamout );
        
        
        fclose(sparamout);*/

        fprintf(lparamout, "\n%s", "#Laser-Parameters,  E0,          I0,              w0,             Ncycles,           ellipticity,           CEP,           t0,           theta0,     envelope (0 - gauss, 1 - sin2) ");
        LaserParam tempLaserParam;
        for (int i = 0; i < sbe->fpulses->Npulses; ++i)
        {
            tempLaserParam = sbe->fpulses->PulseParam[i];
            fprintf(lparamout, "\n          %e      %e      %e      %e      %e      %e      %e      %e      %d   \n\n", 
                tempLaserParam.E0, tempLaserParam.I0, tempLaserParam.w0, tempLaserParam.cycles0, 
                tempLaserParam.e, tempLaserParam.cep0, tempLaserParam.t0, tempLaserParam.theta0, tempLaserParam.envelope);
        }
        fflush(lparamout );
        fclose(lparamout);
        
        
        

        //###############################
        //terminal Output of grid
        /*cout <<"\n\n*===================================*\nMomentum Grid Features";
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
        cout << "\n\nky         = "            << g.k[1][0];*/
        cout << "\n/*************************/\n";
        
        
        
        //###############################
        //Saving or Output of laser
        complex Etemp, Atemp;
        for ( n = 0; n < sbe->fpulses->Nt; n++)
        {
            
            at = sbe->fpulses->atmin + sbe->fpulses->dt*n;


            Etemp = sbe->fpulses->elaser(  &at );
            Atemp = sbe->fpulses->avlaser( &at );
                

            fprintf( laserout,"%e    %e    %e    %e    %e\n", at, real(Etemp), imag(Etemp), real(Atemp), imag(Atemp) );
            
            
        }
        
        fflush( laserout );
        fclose( laserout );   
        
    }
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    
    blockcounts = new int[Nprocessors];
    displs = new int[Nprocessors];

    //int Npoints = sbe->kgrid->N[0];
    int Npoints = sbe->kgrid->Ntotal;
    for (itemp = 0; itemp < Nprocessors; ++itemp)
    {
        //ParaRange(g.N[1], 0, Nprocessors, itemp, &jstart, &jend);
        ParaRange(Npoints, 0, Nprocessors, itemp, &jstart, &jend);

        //blockcounts[itemp] = (jend - jstart) * g.N[0];
        blockcounts[itemp] = (jend - jstart);
        //displs[itemp] = jstart * g.N[0];
        displs[itemp] = jstart;
    }
    //ParaRange(g.N[1], 0, Nprocessors, rank, &jstart, &jend);
    ParaRange(Npoints, 0, Nprocessors, rank, &jstart, &jend);
    // if (rank == MASTER)
    // {
    //     for (int i = 0; i < Nprocessors; ++i) cout << blockcounts[i] << " ";
    //     cout << endl;
    //     for (int i = 0; i < Nprocessors; ++i) cout << displs[i] << " ";
    //     cout << endl;

    // }
    //########################################
    
    
    
    
    
    cout << "\nrank = " << rank << "  uses j = " << jstart << ";  to  j = " << jend << "  with Nprocesses = " << Nprocessors <<endl << endl;
    
    
    if(rank==MASTER)
    {
        

        
        cout.precision( 5 );
        
        
        
        
        cout << "\n\n\n//############################################//" ;
        cout <<"\n#*********==============**********#\n ";
        cout << "Beginning of Momentum/Time integration loops\n (Nx,Ny,Nt)";
        //cout << " = ( "<< g.N[0] << " , " << g.N[1] <<" , " <<fpulse.NewNt << " ) " << scientific;
        cout << "\nMomentum Steps and Time Step";
        //cout << "\n ( dkx, dky, dt )     =   (" << g.dk[0] << ",  " << scientific;
        //cout << g.dk[1] << ",  " << dt << "  )  a.u." << scientific << endl;
        cout <<"############################################"<<endl;
        
        
        cout << "\n\nktime     " << "       ";
        cout << "         nc      " <<  "        pi(a.u.)    ";
        cout << "   gauge1      id-rank \n" << scientific;
        
        cout << "\n\n valence band number / total band number: " << sbe->Nvb << "   /   " << sbe->Nband;
    }

    ktemp = 0;
    int tempy = (sbe->kgrid->N[1]-1)/2;
    int tempz = 0;
    //int tempSK = 0; jstart = tempSK;    jend = tempSK+1;
    /*if (rank == MASTER)
    {
        for (int m = 0; m < Nband; ++ m)
        {
            cout << "m, E(m) : " << m << ", " << sbe->edispersion[tempSK][m] << "\n";
        }
        for (int m =0; m < Nband; ++m)
        {
            for (int n = 0; n < Nband; ++n)
            {
                cout << "m, n, px(m,n), py(m,n) : " << m << ", " << n << ", " << sbe->pmatrix[0][tempSK][m*Nband + n] << sbe->pmatrix[1][tempSK][m*Nband + n] << "\n";
            }
        }
    }*/
    //#####################################
    //#####################################
    //Time integration loop
    for( ktime = 0; ktime<sbe->fpulses->Nt; ktime++ )
    { 
    
        //###############################
        //Momentum and time integration
        for( jtemp = jstart; jtemp < jend; jtemp++ )
        {          
            //kindex = g.index( &itemp, &jtemp, &ktemp);
            kindex = jtemp;
            //kindex = sbe->kgrid->index(&jtemp, &tempy, &tempz);


            
            //##########################################################
            //MOMENTUM INTEGRALS FOR INTER AND INTRA BAND POL. AND CUR.
            //##########################################################
            //##########################################
            //CALCULATION OF INTER-BAND POLARIZATION
            //##########################################
            //inter_rad[0][ktime]+= 2.*real( conj( pmatrix[0][kindex] )*rkutta_pi[ g.index( &itemp, &jtemp, &ktemp) ] )*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ] ;
            
            
            //inter_rad[1][ktime]+= 2.*real( conj( pmatrix[1][kindex] )*rkutta_pi[ g.index( &itemp, &jtemp, &ktemp) ] )*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];
            tempCurrent = sbe->GenInterCurrent(kindex) * sbe->kgrid->weight[ kindex ] ;
            
            inter_rad[0][ktime] += real(tempCurrent);
            inter_rad[1][ktime] += imag(tempCurrent);
            //END OF INTER-BAND CALCULATIONS
            //############################################

            
            
            //################################
            //INTRABAND CURRENTs
            //intra_rad[0][ktime]+= classical_velocity_v[0] + classical_velocity_c[0] ;
            //intra_rad[1][ktime]+= classical_velocity_v[1] + classical_velocity_c[1] ;

            tempCurrent = sbe->GenIntraCurrent(kindex) * sbe->kgrid->weight[ kindex ] ;
            
            intra_rad[0][ktime] += real(tempCurrent);
            intra_rad[1][ktime] += imag(tempCurrent);
            //END OF INTER-BAND CALCULATIONS
            //################################
            
            
            
            
            //##############################################
            //Evaluation of occupations and coherence
            // Only keep 
            //ConductionOccup[ktime]+=   rkutta_nc[ g.index( &itemp, &jtemp, &ktemp ) ]*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];            
            //Coherence_Mom_Int[ktime]+= rkutta_pi[ g.index( &itemp, &jtemp, &ktemp ) ]*g.weight[ g.index(  &itemp, &jtemp, &ktemp ) ];

            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    //if (m > n) continue;
                    density_matrix_integrated[ktime][m*Nband + n] += sbe->dmatrix[kindex][m*Nband + n] * sbe->kgrid->weight[ kindex ];
                }
            }

            sbe->RunSBEs(kindex, ktime);

            //##############################################
    
        }//End of Loop on x-direction
            
        //###################################
        //OutPut for diagnostic
        if( ktime%( (sbe->fpulses->Nt-1)/200 )==0 )
        {
            /*printf( "\n%.6d / %.7d        %.3e        %.3e        %.3e        %.3e        %.3e        %.3e        %.3e ", ktime, sbe->fpulses->Nt-1, 
                real(density_matrix_integrated[ktime][(Nvb-1)*Nband + (Nvb-1)]  ),
                real(density_matrix_integrated[ktime][Nvb*Nband + Nvb]   ), 
                real(density_matrix_integrated[ktime][(Nvb+1)*Nband + (Nvb+1)]  ),
                real(density_matrix_integrated[ktime][(Nvb-1)*Nband + Nvb] ), imag(density_matrix_integrated[ktime][(Nvb-1)*Nband + Nvb]  ),
                real(density_matrix_integrated[ktime][(Nvb-1)*Nband + (Nvb+1)] ), imag(density_matrix_integrated[ktime][(Nvb-1)*Nband + (Nvb+1)]  )  );*/
            printf( "\n%.6d / %.7d        %.3e        %.3e ", ktime, sbe->fpulses->Nt-1, 
                real(density_matrix_integrated[ktime][Nvb*Nband + Nvb]   ), real(density_matrix_integrated[ktime][(Nvb-1)*Nband + Nvb] ));
            
            /*if ( g.k[1][jtemp-1] >= kyShift[0] & g.k[1][jtemp-1] <= kyShift[1] & fbz_shift == 1)
            {
                    
                printf( "\n%.6d / %.7d        %.3e        %.3e        %.2d        %.6d     %s  %.3f %s %.3f %s", ktime, fpulse.NewNt-1, real(ConductionOccup[ktime]*g.dV  ), real(Coherence_Mom_Int[ktime])*g.dV, cs.gauge, rank, "Shifted BZ (",tkx5[ 0 ],",",tky5[ 0 ],")" );
                
            }
            else
            {
                    
                printf( "\n%.6d / %.7d        %.3e        %.3e        %.2d        %.6d     %s  %.3f %s %.3f %s", ktime, fpulse.NewNt-1, real(ConductionOccup[ktime]*g.dV  ), real(Coherence_Mom_Int[ktime])*g.dV, cs.gauge, rank, "Non-Shift BZ (",tkx5[ 0 ],",",tky5[ 0 ],")" );
                    
            }*/
                
        }//END OF DIAGNOSTIC
        //###################################

        if (shotNumber > 0)    
            MPI_Barrier( MPI_COMM_WORLD );
        
    }//End of Time integration loop

    //##############################################
    // inter, intrabnad currents, occupation, cohereneces are integrated through every sub-grids (all-processors)
    // inter, intra - currents and radiation by group and anomalous velocity
    MPI_Barrier( MPI_COMM_WORLD );
    for (int i = 0; i < 2; ++i)
    {
        if (rank == MASTER)
        {
            MPI_Reduce(MPI_IN_PLACE, inter_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Reduce(inter_rad[i], inter_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
        }
        MPI_Barrier( MPI_COMM_WORLD );
        if (rank == MASTER)
        {
            MPI_Reduce(MPI_IN_PLACE, intra_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Reduce(intra_rad[i], intra_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
        }
        MPI_Barrier( MPI_COMM_WORLD );
    }
    if (rank == MASTER)
    {
        MPI_Reduce(MPI_IN_PLACE, &density_matrix_integrated[0][0], sbe->fpulses->Nt * Nband*Nband, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Reduce(&density_matrix_integrated[0][0], &density_matrix_integrated[0][0], sbe->fpulses->Nt * Nband*Nband, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    
    double occup_temp = 0;
    if (rank == MASTER)
    {
        cout << "\n\nWriting output-data\n\n";

        for (ktime = 0; ktime < sbe->fpulses->Nt; ++ktime)
        {
            fprintf(interj_out, "%.16e    %.16e\n", real(inter_rad[0][ktime]) * dV, real(inter_rad[1][ktime]) * dV);
            fprintf(intraj_out, "%.16e    %.16e\n", real(intra_rad[0][ktime]) * dV, real(intra_rad[1][ktime]) * dV);
            occup_temp = 0;
            for (int m = 0; m < Nband; ++m)
            {
                fprintf(occup_out, "%.16e    ", real(density_matrix_integrated[ktime][m*Nband + m]) / (Npoints-1));
                occup_temp += real(density_matrix_integrated[ktime][m*Nband + m]);
            }
            fprintf(occup_out, "%.16e\n", occup_temp / (Npoints-1));
        }
        fflush(interj_out);
        fflush(intraj_out);
        fflush(occup_out);
    }
    
    //fclose(  fout     );
    if (rank==MASTER)    
    {       

        fclose( interj_out );
        fclose( intraj_out );
        fclose( occup_out);
        
    }   
    
    
    for ( itemp=0; itemp < Ngrad; itemp++)
    {
        free( inter_rad[itemp] );
        free( intra_rad[itemp] );
    
        
    }

    delete[] blockcounts;
    delete[] displs;

    Delete2D<complex>(density_matrix_integrated, Nband*Nband, sbe->fpulses->Nt);
    
    
    //fclose(simulation_out);
    
    
    cout << "\n***********************\nEND OF THE PROGRAM with No. of cores = " << Nprocessors << "  uses j = " << jstart << ";  to  j = " << jend << "  with rank = "  << rank << "\n\n" <<  endl;
    cout << "=====================\n\n";
    
    
    MPI_Finalize( )	;	
	exit(0);


    
};


void der_pi( const double  *eg, const double *avg, const double *T2, const complex *apcv, const complex *cohPi0, const double *_nc, complex *cpi)
{
    
    *cpi = -I*( ( *eg + *avg - I/(*T2) * flag_t2 ) *( *cohPi0 ) + (*apcv)*( 1. - 2.*(*_nc ) )   );
    
}

void der_nc( const complex *ROmega, const complex *cpi,  double *nc )
{
    
    *nc = -2.*imag( conj( *ROmega )*( *cpi ) );
    
}


void der_nc2( const complex *apcv, const complex *cpi, const double *T1, const double *tnc,  double *nc )
{
    
    *nc = -2.*imag( conj( *apcv )*( *cpi ) ) - (*tnc)/(*T1) * flag_occ_dephasing;
    
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
