##########################################
#By Alexis A. Chac\'on S.
#July 12, 2019
################################


The C/C++ version still does not recieve any input file. This is in process of development.






##########################################
DESCRIPTION OF INPUTS in the main program
###################


The input-parameters are given and writen in the main file, hhg_kr4_hm_mpi.cpp. This program generate an execuable that recieve 11 parameters via terminal, those are described by:



//################################
int main( int argc, char *argv[] )
{

    //################################
    //Numerical parameter and set
    iparam        = atof( argv[1] );         //Magnetic flux, or phase of the second coupling parameter of Haldane Model,
    		    	  	  	     //ie. t2*e^{phi0}, iparam is phi0, in can be a double or integer data type
					     
    jparam        = atof( argv[2] );         //On-site potential ratio to t2, i.e M/t2
    
    nxparam       = atoi( argv[3] );         //Nx points along kx-direction

    nyparam       = atoi( argv[4] );         //Nx points along ky-direction 

    kparam        = atof( argv[5] );         //Laser field strength 

    ncyparam      = atoi( argv[6] );         //Number of optical cycles at FWHM 

    dtparam       = atof( argv[7] );         //Time step, param dt

    dephasing     = atof( argv[8] );         //Phenomelogical Dephasing time

    ellip         = atof( argv[9] );         //Laser field Ellipticity

    treg0         = atof( argv[10] );        //Regularization parameter for dipole/connection discontinuity/singularity

    tgauge2       = atoi( argv[11] );        //gauge2 modification of wavefunction gauge from 1 to 2
...


Other Haldane Model (HM) parameters are defined and fixed about the lines 165 -- 185, i.e.

    //#############################
    /*****************************************
     
      HM. Parameters for Honney Comb lattice
     
     ****/
    double la0            = 1.;         //Lattice constant in Angstrom 
    double a0             = la0/0.529;  //Lattice constant in atomic units (a.u.)
    double t1             = 0.075;      //Neirest Neighbour hopping parameters of HM in a.u.
    double t2             = t1/3.;      //Amplitude of the Next Neirest Neighbour of HM
    double phi0           = iparam;     //magnetic flux or phase in radians 
    double M0             = jparam*t2;  //Local or on-site potential





######################
The Laser parameters and time-steps are defined in the lines: 204 -- 235

    //#############################
    /*************************************
    
    //Laser Parameters
     
     */

    int npulses         = 1;             // Number of pulses
    double wfreq        = 0.014;         // Frequency in atomic units (a.u.)
    double period0      = dospi/wfreq;   // Optical Period in a.u.
    double E0           = kparam;       // Laser field strenth
    double I0           = 0.;           // Intensity
    double ncycles      = ncyparam;     // Number of Optical cycles at FWHM or at 1/e
    double cep          = 0.;           // Carrier envelope phase
    
    
    
    
    //Time step dt and other laser params
    double dt           = dtparam;      // Time step
    
    double tstart       = 0.;           //initial-time
    double offset       = 8.*period0;   //time-lenght before laser 
    double outset       = 8.*period0;   //time-length after laser
    double ellipticity  = ellip;        //laser field Ellipticity
    double relativephase= 0.;           //laser field relative components
    double ltheta0      = 0.;           //Inclination angle
    
    string env_name = "gauss";          //Name envel: "rect", "sin2", "gauss" or "rsin2"
    
    
    double T2           = dephasing;     //Dephasing time in the crystal





##########################################

OUTPUTs
##############################


1) mgrid.dat, save the momentum axes in a single column, where the first Nx+3 elements are, skiper, No. Of points for kx direction and momentum step along kx, the next points up to Nx+3 are the kx axis, similar structure for ky, from the position or index point Nx+4 up to Nx+Ny+6, Ny is the number of points along ky direction,...s



2) edispersion.dat, contains the energy dispersion as a function of momentum (kx,ky) in the gnuplot format for the  splot command 
  # kx in the 1rst column in a.u.
  # ky in the 2nd column in a.u.
  # valence band VB energies in a.u. in the 3rd column 
  # conduction band CB energies in a.u. in the 4th column 
  # energy diff between VB and CB in a.u.  in the 5th column 



3) dipoles.dat, contains the dipoles as a function of momentum (kx,ky) in the gnuplot format
  # kx in the 1rst column in a.u.
  # ky in the 2nd column in a.u.
  # real part of x-comp. dipole matrix element VB->CB in a.u. in the 3rd column 
  # imag part of x-comp. dipole matrix element VB->CB in a.u. in the 4th column 
  # real part of y-comp. dipole matrix element VB->CB in a.u. in the 5th column 
  # imag part of y-comp. dipole matrix element VB->CB in a.u. in the 6th column 
  


4) connection.dat, contains the Berry Connection of CB as a function of momentum (kx,ky) in the gnuplot format
  # kx in the 1rst column in a.u.
  # ky in the 2nd column in a.u.
  #  x-comp. Two times Berry Connection of the CB in a.u. in the 3rd column (Berry connection difference between CB and VB)
  #  y-comp. Two times Berry Connection of the CB in a.u. in the 4th column 



5) curvature.dat, contains the energy dispersion as a function of momentum (kx,ky) in the gnuplot format
  # kx in the 1rst column in a.u.
  # ky in the 2nd column in a.u.
  #  z-comp. Berry Curvature of the CB in a.u. in the 3rd column 




6) gvelocities.dat, contains the group velocity VB and CB as a function of momentum (kx,ky) in the gnuplot format
  # kx in the 1rst column in a.u.
  # ky in the 2nd column in a.u.
  # x-comp. Group val. VB in a.u. in the 3rd column 
  # y-comp. Group val. VB in a.u. in the 4th column 
  # x-comp. Group val. CB in a.u. in the 5th column 
  # y-comp. Group val. CB in a.u. in the 6th column




7) setOfparameters.dat, contains the set of input parameters for the laser and Haldane model, or Hexagonal lattice
   Additional information of the min and max direct band gap, Chern number of the CB, number of cores used in the MPI
   calculation is also found in this file. 


 
8) full_integrated_currents_rank*, this is the most important output and contains 
   # 1rst col. "#","#time(a.u.)",
   # 2nd col.  "#Jx_er(a.u.)", 	x-comp. Of the inter-band current  
   # 3rd col.  "#Jy_er(a.u.)",  y-comp. Of the inter-band current 
   # 4th col.  "#Jx_ra(a.u.)",  x-comp. Of the intra-band current
   # 5th col.  "#Jy_ra(a.u.)",  y-comp. Of the intra-band current
   # 6th col.  "#nc",  	        CB occupation Momentum integral 
   # 7th col.  "#Coherent(a.u.)", coherent 
   # 8th col.  "#vgx" , 	x-comp. Of the group velocity contribution to the intra-current
   # 9th col.  "#vgy",  	y-comp. Of the group velocity contribution to the intra-current 
   # 10th col. "#vax", 		x-comp. Of the anomalous velocity contribution to the intra-current 
   # 11th col. "#vay", 		y-comp. Of the anomalous velocity contribution to the intra-current 



9) outlaserdata.dat, this output contains information of the electric field of the laser and its vector potential
   # 1rst col. "#","#time(a.u.)",
   # 2nd col.  "#Ex(a.u.)",  x-comp. of the electric field   
   # 3rd col.  "#Ey(a.u.)",  y-comp. of the electric field 
   # 4th col.  "#Ax(a.u.)",  x-comp. of the vector potential  
   # 5th col.  "#Ay(a.u.)",  y-comp. of the vector potential

 
