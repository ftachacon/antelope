/// main cpp file
/**
 * @author Dasol Kim
 * @author Alexis Agustín  Chacón Salazar
*/

// Standard C/C++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

// 3rd-party headeres
#include "mpi.h"
//#include "mkl.h"
#include <libconfig.h++>
#include <unistd.h>

// My Own Headers
#include "constant.h"
#include "laser.h"
#include "momaxis.h"
#include "SBEs.h"
#include "utility.h"

#define MASTER 0    /* task id of master task or node */


void setting_complex_memory(complex **_TempPointer, long long int _Ntime);
void setting_double_memory(double **_TempPointer, long long int _Ntime );

void My_MPI_CLX_SUM( complex *in, complex *inout, int *len, MPI_Datatype *dptr );

string targetMaterial;

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
	int rank=0, Nprocessors=1;

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );	//Getting the local id or rank process
	MPI_Comm_size( MPI_COMM_WORLD, &Nprocessors );	//Getting the size or number of proccesses or nodes
	MPI_Barrier( MPI_COMM_WORLD );	


    // Test libconfigc++
    // try-catch code are copied from libconfig example
    libconfig::Config cfg;
    try
    {
        cfg.readFile("inputParam.cfg");
    }
    catch(const libconfig::FileIOException &fioex)
    {
        cerr << "I/O error while reading file." << endl;
        return(EXIT_FAILURE);
    }
    catch(const libconfig::ParseException &pex)
    {
        cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << endl;
        return(EXIT_FAILURE);
    }

    // Read gauge from config
    string gaugeString; GaugeType tempGauge;
    cfg.lookupValue("gauge", gaugeString);
    if (gaugeString == "Wannier" || gaugeString == "LengthWannier")
    {
        tempGauge = GaugeType::LengthWannier;
    }
    else if (gaugeString == "Hamiltonian" || gaugeString == "LengthHamiltonian")
    {
        tempGauge = GaugeType::LengthHamiltonian;
    }
    else if (gaugeString == "Velocity" || gaugeString == "VelocityHamiltonian")
    {
        tempGauge = GaugeType::VelocityHamiltonian;
    }
    else
    {
        if (rank == MASTER) cout << "Warning: No gauge type is specified. Choose LengthWannier for default gauge.\n";
        tempGauge = GaugeType::LengthWannier;
    }

    SBEs *sbe = new SBEs(&cfg.getRoot(), tempGauge);

    double dV = sbe->kmesh->dV;
    int Nband = sbe->Nband;

	MPI_Barrier( MPI_COMM_WORLD );	

    // if (rank == MASTER)
    // {
    //     ofstream tempInit, currentInitx, currentInity;
    //     tempInit.open("InitWannier.dat");
    //     currentInitx.open("InitCurrentx.dat");
    //     currentInity.open("InitCurrenty.dat");
    //     array<double, Ndim> tempCurrent;
    //     for (int i = 0; i < sbew->kmesh->N[0]; ++i)
    //     {
    //         for (int j = 0; j < sbew->kmesh->N[1]; ++j)
    //         {
    //             tempInit << sbew->dmatrix[sbew->kmesh->index(i, j, 0)][0] << "   ";
    //             tempCurrent = sbew->GenCurrent(sbew->kmesh->index(i, j, 0), sbew->fpulses->atmin);
    //             currentInitx << tempCurrent[0] << "   " ;
    //             currentInity << tempCurrent[1] << "   " ;
    //         }
    //         tempInit << endl;
    //         currentInitx << endl;
    //         currentInity << endl;
    //     }
    // }	


    //###############################
    //Variables
    int jstart, jend, ktemp;
    complex *inter_rad[Ngrad];
    complex *intra_rad[Ngrad];
    complex **density_matrix_integrated;

    double at=0, nv=0.;

    complex tgv_v[Ngrad] = {0.,0.};
    complex tgv_c[Ngrad] = {0.,0.};
    complex tdip_cv[Ngrad]={0.,0.};
    
    // Number of k-grid for each process
    // and their displacement in global buffer
    int *blockcounts;
    int *displs;

    int kindex;

    array<double, Ndim> tempCurrent;
    
    //######################################
    //Opening output files!
    FILE *laserout, *mout, *fout, *engout;
    FILE *inh_out, *occup_out, *sparamout, *lparamout;
    FILE *simulation_out;
    FILE *density_out, *densitytime_out, *polarization_out;
    FILE *intraj_out, *interj_out;
    
    
    if ( rank == MASTER )
    {
        laserout         = fopen( "outlaserdata.dat","w" ); // Laser field characteristics,time axis,
        lparamout        = fopen( "laserParameters.dat", "w" );
        interj_out       = fopen( "interband_dipole_full_evol.dat", "w");
        intraj_out       = fopen( "intraband_current_full_evol.dat", "w");
        occup_out        = fopen( "occupation__full__evol.dat", "w");
    }
    
    

    if (rank == MASTER)
    {
        sbe->PrintInfo();
    }

    //##########################
    //Coherence and occupations
    
    for (int itemp=0; itemp<Ngrad; itemp++ )
    {        
        
        setting_complex_memory( &inter_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &intra_rad[itemp], sbe->fpulses->Nt );
        
        for (int jtemp=0; jtemp<sbe->fpulses->Nt; jtemp++)
        {
            inter_rad[itemp][jtemp]=0.;
            intra_rad[itemp][jtemp]=0.;
        }
        
    }

    density_matrix_integrated = Create2D<complex>(sbe->fpulses->Nt, Nband*Nband);

    
    if (rank == MASTER)
    {
        //###############################
        //Saving or Output of laser
        array<double, Ndim> Etemp, Atemp;
        for (int n = 0; n < sbe->fpulses->Nt; n++)
        {
            at = sbe->fpulses->atmin + sbe->fpulses->dt*n;

            Etemp = sbe->fpulses->elaser(  at );
            Atemp = sbe->fpulses->avlaser( at );
                
            fprintf( laserout,"%e    %e    %e    %e    %e    %e    %e\n", at, Etemp[0], Etemp[1], Etemp[2], Atemp[0], Atemp[1], Atemp[2] );

        } 
        fflush( laserout );
        fclose( laserout );   
    }
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    
    blockcounts = new int[Nprocessors];
    displs = new int[Nprocessors];

    int Npoints = sbe->kmesh->Ntotal;
    for (int itemp = 0; itemp < Nprocessors; ++itemp)
    {
        ParaRange(Npoints, 0, Nprocessors, itemp, &jstart, &jend);

        blockcounts[itemp] = (jend - jstart);
        displs[itemp] = jstart;
    }
    ParaRange(Npoints, 0, Nprocessors, rank, &jstart, &jend);
 
    
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
    }

    ktemp = 0;
    double currenttime, currentcycle;
    complex *temp_density_matrix_integrated = new complex[Nband*Nband];
    double occup_temp = 0;
    double maxcycles = 2*sbe->fpulses->pulses[0].twidth / sbe->fpulses->pulses[0].period0;
    //#####################################
    //#####################################
    //Time integration loop
    for(int ktime = 0; ktime<sbe->fpulses->Nt; ktime++ )
    { 
        currenttime = sbe->fpulses->atmin + sbe->fpulses->dt * ktime;
        //###############################
        //Momentum and time integration
        for(int jtemp = jstart; jtemp < jend; jtemp++ )
        {
            kindex = jtemp;

            tempCurrent = sbe->GenCurrent(kindex, currenttime);

            inter_rad[0][ktime] += tempCurrent[0] * sbe->kmesh->weight[ kindex ] ;
            inter_rad[1][ktime] += tempCurrent[1] * sbe->kmesh->weight[ kindex ] ;

            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    //if (m > n) continue;
                    density_matrix_integrated[ktime][m*Nband + n] += sbe->dmatrix[kindex][m*Nband + n] * sbe->kmesh->weight[ kindex ];
                }
            }
            sbe->RunSBEs(kindex, currenttime);

            //##############################################

        }
        // Snapshots - including diagnostics and density plots 
        // if (shotNumber > 0 && ktime % shotNumber == 0)
        if (ktime % 200 == 0)
        {
            for (int ij = 0; ij < Nband * Nband; ++ij)
            {
                temp_density_matrix_integrated[ij] = density_matrix_integrated[ktime][ij];
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == MASTER)
            {
                MPI_Reduce(MPI_IN_PLACE, &temp_density_matrix_integrated[0], Nband*Nband, 
                    MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Reduce(&temp_density_matrix_integrated[0], &temp_density_matrix_integrated[0], Nband*Nband, 
                    MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }
            if (rank == MASTER)
            {
                currentcycle = currenttime - (sbe->fpulses->pulses[0].t0 - sbe->fpulses->pulses[0].twidth);
                if (currentcycle < 0) 
                {
                    currentcycle = 0.;
                }
                else
                {
                    currentcycle /= sbe->fpulses->pulses[0].period0;
                }
                
                if (currentcycle > maxcycles)
                {
                    currentcycle = maxcycles;
                }
                occup_temp = 0;
                for (int m = 0; m < Nband; ++m)
                {
                    occup_temp += real(temp_density_matrix_integrated[m*Nband + m]);
                }
                occup_temp /= Npoints;
                cout << ktime << " / " << sbe->fpulses->Nt << ", " << currentcycle << " / " << maxcycles;
                cout << ", nc = " << temp_density_matrix_integrated[sbe->material->Nval*Nband + sbe->material->Nval] / static_cast<double>(Npoints);
                cout << ", ntotal = " << occup_temp << endl;
            }
        }

    } //End of Time integration loop
    delete[] temp_density_matrix_integrated;

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
    
    if (rank == MASTER)
    {
        cout << "\n\nWriting output-data\n\n";

        for (int ktime = 0; ktime < sbe->fpulses->Nt; ++ktime)
        {
            fprintf(interj_out, "%.16e    %.16e\n", real(inter_rad[0][ktime]) * dV, real(inter_rad[1][ktime]) * dV);
            fprintf(intraj_out, "%.16e    %.16e\n", real(intra_rad[0][ktime]) * dV, real(intra_rad[1][ktime]) * dV);
            occup_temp = 0;
            for (int m = 0; m < Nband; ++m)
            {
                fprintf(occup_out, "%.16e    ", real(density_matrix_integrated[ktime][m*Nband + m]) / Npoints);
                occup_temp += real(density_matrix_integrated[ktime][m*Nband + m]);
            }
            fprintf(occup_out, "%.16e\n", occup_temp / Npoints);
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
    
    
    for (int itemp=0; itemp < Ngrad; itemp++)
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

// Plot 2d
// void Plot2D( string _filename, int skip0 = 1, int skip1 = 1, int skip2 = 1 )