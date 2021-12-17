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
#include <time.h>
#include <tuple>
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

typedef std::tuple<array<double, Ndim>,  array<array<double, Ndim>, Nband> , array<array<double, Ndim>, Ndim> > tupleSpinCharge;

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
 clock_t tStart = clock();
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
        cout << "debug.\n";
    }
    else
    {
        if (rank == MASTER) cout << "Warning: No gauge type is specified. Choose LengthWannier for default gauge.\n";
        tempGauge = GaugeType::LengthWannier;
    }
    // cout << "debug.before.sbe\n";
    SBEs *sbe = new SBEs(&cfg.getRoot(), tempGauge);
    // cout << "debug.after.sbe\n";
    double dV = sbe->kmesh->dV;
    int Nband = sbe->Nband;

	MPI_Barrier( MPI_COMM_WORLD );

    //###############################
    //Variables
    int jstart, jend, ktemp;
    complex *inter_rad[Ngrad];
    complex *intra_rad[Ngrad];
    complex *w1_rad[Ngrad];
    complex *w2_rad[Ngrad];
    complex *w3_rad[Ngrad];
    complex *w4_rad[Ngrad];
    complex *inter_rad_spin[Ngrad][Ndim];
    complex *inter_bands[Ngrad][4];
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

    //######################################
    //Opening output files!
    FILE *laserout, *mout, *fout, *engout;
    FILE *inh_out, *occup_out, *sparamout, *lparamout;
    FILE *simulation_out;
    FILE *density_out, *densitytime_out, *polarization_out;
    FILE *intraj_out, *interj_out, *interj_out_spin,*interj_out_bands ;
    FILE *w1j_out, *w2j_out,*w3j_out, *w4j_out ;
    array<double, Ndim> tempCurrent;
    array<array<double, Ndim>, 4> tempCurrent_bands;
    array<array<double, Ndim>, Ndim> tempSpinCurrent;
    if ( rank == MASTER )
    {
        laserout         = fopen( "outlaserdata.dat","w" ); // Laser field characteristics,time axis,
        lparamout        = fopen( "laserParameters.dat", "w" );
        interj_out       = fopen( "interband_dipole_full_evol.dat", "w");
        intraj_out       = fopen( "intraband_current_full_evol.dat", "w");
        interj_out_spin  = fopen( "spin_current.dat", "w");
        occup_out        = fopen( "occupation__full__evol.dat", "w");
        interj_out_bands= fopen( "current_bands_resolved.dat", "w");
        w1j_out       = fopen( "partial_current_w1j.dat", "w");
        w2j_out       = fopen( "partial_current_w2j.dat", "w");
        w3j_out       = fopen( "partial_current_w3j.dat", "w");
        w4j_out       = fopen( "partial_current_w4j.dat", "w");
}


    if (rank == MASTER)
    {
        sbe->PrintInfo();
    }
    //##########################
    //Coherence and occupations
    for (int itemp=0; itemp<Ngrad; itemp++ )
    {
        setting_complex_memory( &w1_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &w2_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &w3_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &w4_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &inter_rad[itemp], sbe->fpulses->Nt );
        setting_complex_memory( &intra_rad[itemp], sbe->fpulses->Nt );

        for (int saxis =0; saxis<Ndim;++saxis)
        {
        setting_complex_memory( &inter_rad_spin[itemp][saxis], sbe->fpulses->Nt );
        }
        for (int saxis =0; saxis<Nband;++saxis)
        {
        setting_complex_memory( &inter_bands[itemp][saxis], sbe->fpulses->Nt );
        }


        for (int jtemp=0; jtemp<sbe->fpulses->Nt; jtemp++)
        {
            inter_rad[itemp][jtemp]=0.;
            intra_rad[itemp][jtemp]=0.;
            for (int saxis =0; saxis<Ndim;++saxis)
            {
            inter_rad_spin[itemp][saxis][jtemp]=0.;
            }
           for (int saxis =0; saxis<Nband;++saxis)
        {
         inter_bands[itemp][saxis][jtemp]=0.;
          }
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
        cout << "Beginning of Momentum/Time integration loops\n (Nx,Ny,Nz?,Nt)";
        //cout << " = ( "<< g.N[0] << " , " << g.N[1] <<" , " <<fpulse.NewNt << " ) " << scientific;
        cout << "\nMomentum Steps and Time Step";
        //cout << "\n ( dkx, dky, dt )     =   (" << g.dk[0] << ",  " << scientific;
        //cout << g.dk[1] << ",  " << dt << "  )  a.u." << scientific << endl;
        cout <<"############################################"<<endl;
    }

    ktemp = 0;
    double currenttime, currentcycle;
    complex *temp_density_matrix_integrated = new complex[Nband*Nband];
    complex *temp_1d_integrated_occupation = new complex[sbe->Nk[0]];
   // complex *temp_2d_occupation = new complex[sbe->Nk[0],Nband*Nband];
   FILE *dm1dout;
    if (rank == MASTER)
    {
        dm1dout = fopen("occupation_1d_integration.dat", "w");
    }
    double occup_temp = 0;
    double maxcycles = 2*sbe->fpulses->pulses[0].twidth / sbe->fpulses->pulses[0].period0;

    int shotFrequency = int(sbe->fpulses->pulses[0].period0/10/sbe->fpulses->dt);
    //#####################################
    //#####################################
    //Time integration loop
    // cout << "debug.before.loop\n";
    for(int ktime = 0; ktime<sbe->fpulses->Nt; ktime++ )
    {
        currenttime = sbe->fpulses->atmin + sbe->fpulses->dt * ktime;
        //###############################
        //Momentum and time integration
        // cout << "debug.loop.time\n";
        for(int jtemp = jstart; jtemp < jend; jtemp++ )
        {
            kindex = jtemp;

            tupleSpinCharge temp = sbe->GenCurrent(kindex, currenttime);

            tempCurrent = std::get<0>(temp);
            tempSpinCurrent = std::get<2>(temp);
            tempCurrent_bands = std::get<1>(temp);


	//        cout << tempCurrent[2]<<"debug.loop.current\n";


	          inter_rad[0][ktime] += tempCurrent[0] * sbe->kmesh->weight[ kindex ] ;
            inter_rad[1][ktime] += tempCurrent[1] * sbe->kmesh->weight[ kindex ] ;
            inter_rad[2][ktime] += tempCurrent[2] * sbe->kmesh->weight[ kindex ] ;

            auto kxi=sbe->kmesh->kgrid[kindex][0];
            auto kyi=sbe->kmesh->kgrid[kindex][1];
            auto kzi=sbe->kmesh->kgrid[kindex][2];

            auto rx1 = sbe->rx1;
            auto ry1 = sbe->ry1;
            auto rz1 = sbe->rz1;

            auto rx2 = sbe->rx2;
            auto ry2 = sbe->ry2;
            auto rz2 = sbe->rz2;

            auto rx3 = sbe->rx3;
            auto ry3 = sbe->ry3;
            auto rz3 = sbe->rz3;

            auto rx4 = sbe->rx4;
            auto ry4 = sbe->ry4;
            auto rz4 = sbe->rz4;

            double rmax2;

            double distance_w1, distance_w2,distance_w3, distance_w4;

	          distance_w1 = (kzi-rz1)*(kzi-rz1)+(kyi-ry1)*(kyi-ry1)+(kxi-rx1)*(kxi-rx1);
	          distance_w2 = (kzi-rz2)*(kzi-rz2)+(kyi-ry2)*(kyi-ry2)+(kxi-rx2)*(kxi-rx2);
            distance_w3 = (kzi-rz3)*(kzi-rz3)+(kyi-ry3)*(kyi-ry3)+(kxi-rx3)*(kxi-rx3);
            distance_w4 = (kzi-rz4)*(kzi-rz4)+(kyi-ry4)*(kyi-ry4)+(kxi-rx4)*(kxi-rx4);

            rmax2=rmax*rmax;

            for (int iband =0; iband < Nband; ++iband)
            {
              inter_bands[0][iband][ktime] += tempCurrent_bands[0][iband] * sbe->kmesh->weight[ kindex ] ;
              inter_bands[1][iband][ktime] += tempCurrent_bands[1][iband] * sbe->kmesh->weight[ kindex ] ;
              inter_bands[2][iband][ktime] += tempCurrent_bands[2][iband] * sbe->kmesh->weight[ kindex ] ;

              if (iband == 1 || iband == 2)
              {
                for ( int iaxis=0; iaxis<Ndim; iaxis++)
                {
                  if (distance_w1 < rmax2 )
                  {
                   w1_rad[iaxis][ktime] += tempCurrent_bands[iaxis][iband] * sbe->kmesh->weight[ kindex ] ;
                   }

                 else if (distance_w2 < rmax2 )
                  {
                   w2_rad[iaxis][ktime] += tempCurrent_bands[iaxis][iband] * sbe->kmesh->weight[ kindex ] ;
                   }

                 else if (distance_w3 < rmax2 )
                 {
                   w3_rad[iaxis][ktime] += tempCurrent_bands[iaxis][iband] * sbe->kmesh->weight[ kindex ] ;
                   }
                 else if (distance_w4 < rmax2 )
                 {
                   w4_rad[iaxis][ktime] += tempCurrent_bands[iaxis][iband] * sbe->kmesh->weight[ kindex ] ;
                   }
                }
            }
		}


            //////// spinCurrent
            for (int iaxis =0; iaxis < Ndim; ++iaxis)
            {
              for (int saxis =0; saxis < Ndim; ++saxis)
              {
                inter_rad_spin[saxis][iaxis][ktime] +=  tempSpinCurrent[saxis][iaxis] * sbe->kmesh->weight[ kindex ] ;
              }
            }
            //////// spinCurrent

           // cout << "debug.loop.current2\n";
            for (int m = 0; m < Nband; ++m)
            {
                for (int n = 0; n < Nband; ++n)
                {
                    //if (m > n) continue;
                    density_matrix_integrated[ktime][m*Nband + n] += sbe->dmatrix[kindex][m*Nband + n] * sbe->kmesh->weight[ kindex ];
      //          cout <<ktime<<" "<<kindex<<" "<<n<<" "<<m<<" "<<density_matrix_integrated[ktime][m*Nband + n]<<"\n";
   		         }
            }
            sbe->RunSBEs(kindex, currenttime);

            //##############################################

        }
        // Snapshots - including diagnostics and density plots
        // if (shotNumber > 0 && ktime % shotNumber == 0)
        if (ktime % shotFrequency == 0)
        {
            // temporary routine - integrating occupation
            for (int m = 0; m < sbe->Nk[0]; ++m)
            {
                temp_1d_integrated_occupation[m] = 0.;
            }
            for(int jtemp = jstart; jtemp < jend; jtemp++ )
            {
                auto ki = sbe->kmesh->index(jtemp);
                if (sbe->gauge == GaugeType::LengthWannier)
                {
                    sbe->WannierToHamiltonian(sbe->newdMatrix, sbe->dmatrix[jtemp], jtemp, currenttime);
                   // for (int m = 0; m < Nband; ++m)
                   // {
       		         //    for (int n = 0; n < Nband; ++n)
        	         //     {
                   //       cout<<ktime<<" "<<ki[0]<<" "<<ki[1]<<" "<<ki[2]<<" "<<m<<" "<<n<<" "<< real(sbe->newdMatrix[Nband*m + n])<<" "<< imag(sbe->newdMatrix[Nband*m + n])<<"\n";
                   //     }
   	               //   }

                    temp_1d_integrated_occupation[ki[0]] += sbe->newdMatrix[sbe->material->Nval*Nband + sbe->material->Nval] * sbe->kmesh->weight[ jtemp ];
                }
                else
                {
                    temp_1d_integrated_occupation[ki[0]] += sbe->dmatrix[jtemp][sbe->material->Nval*Nband + sbe->material->Nval] * sbe->kmesh->weight[ jtemp ];
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == MASTER)
            {
                MPI_Reduce(MPI_IN_PLACE, &temp_1d_integrated_occupation[0], sbe->Nk[0],
                    MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
          }
            else
            {
                MPI_Reduce(&temp_1d_integrated_occupation[0], &temp_1d_integrated_occupation[0], sbe->Nk[0],
                    MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }
            if (rank == MASTER)
            {
                for (int m = 0; m < sbe->Nk[0]; ++m)
                {
                    fprintf(dm1dout, "%16e    ", real(temp_1d_integrated_occupation[m]) / static_cast<double>(sbe->Nk[1]) );
                }
                fprintf(dm1dout, "\n");
            }
            // -----------

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
    delete[] temp_1d_integrated_occupation;

    //##############################################
    // inter, intrabnad currents, occupation, cohereneces are integrated through every sub-grids (all-processors)
    // inter, intra - currents and radiation by group and anomalous velocity
    MPI_Barrier( MPI_COMM_WORLD );
    for (int i = 0; i < 3; ++i)
    {
        if (rank == MASTER)
        {
            MPI_Reduce(MPI_IN_PLACE, inter_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, intra_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w1_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w2_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w3_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w4_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            for (int j=0; j<3 ;j++)
            {
              MPI_Reduce(MPI_IN_PLACE, inter_rad_spin[i][j], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }

          for (int j=0; j<Nband ;j++)
            {
              MPI_Reduce(MPI_IN_PLACE, inter_bands[i][j], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Reduce(MPI_IN_PLACE, inter_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, intra_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w1_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w2_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w3_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, w4_rad[i], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
           for (int j=0; j<3 ;j++)
            {
              MPI_Reduce(MPI_IN_PLACE, inter_rad_spin[i][j], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
            }
                   for (int j=0; j<Nband ;j++)
            {

              MPI_Reduce(MPI_IN_PLACE, inter_bands[i][j], sbe->fpulses->Nt, MPI_DOUBLE_COMPLEX, MPI_CLX_SUM, MASTER, MPI_COMM_WORLD);
              }

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
            fprintf(interj_out, "%.16e    %.16e   %.16e\n", real(inter_rad[0][ktime]) * dV, real(inter_rad[1][ktime]) * dV,real(inter_rad[2][ktime]) * dV);

            fprintf(interj_out_spin, "%.16e    %.16e   %.16e  %.16e    %.16e   %.16e %.16e    %.16e   %.16e\n",  real(inter_rad_spin[0][0][ktime]) * dV, real(inter_rad_spin[1][1][ktime]) * dV,real(inter_rad_spin[2][2][ktime]) * dV,real(inter_rad_spin[0][1][ktime]) * dV,real(inter_rad_spin[0][2][ktime]) * dV,real(inter_rad_spin[1][0][ktime]) * dV,real(inter_rad_spin[1][2][ktime]) * dV,real(inter_rad_spin[2][0][ktime]) * dV ,real(inter_rad_spin[2][1][ktime]) * dV);


            fprintf(intraj_out, "%.16e  %.16e  %.16e\n", real(intra_rad[0][ktime]) * dV, real(intra_rad[1][ktime]) * dV, real(intra_rad[2][ktime]) * dV);
            fprintf(w1j_out, "%.16e  %.16e  %.16e\n", real(w1_rad[0][ktime]) * dV, real(w1_rad[1][ktime]) * dV, real(w1_rad[2][ktime]) * dV);
            fprintf(w2j_out, "%.16e  %.16e  %.16e\n", real(w2_rad[0][ktime]) * dV, real(w2_rad[1][ktime]) * dV, real(w2_rad[2][ktime]) * dV);
            fprintf(w3j_out, "%.16e  %.16e  %.16e\n", real(w3_rad[0][ktime]) * dV, real(w3_rad[1][ktime]) * dV, real(w3_rad[2][ktime]) * dV);
            fprintf(w4j_out, "%.16e  %.16e  %.16e\n", real(w4_rad[0][ktime]) * dV, real(w4_rad[1][ktime]) * dV, real(w4_rad[2][ktime]) * dV);

            fprintf(interj_out_bands, "%.16e  %.16e  %.16e  %.16e %.16e  %.16e  %.16e  %.16e %.16e  %.16e  %.16e  %.16e \n",real(inter_bands[0][0][ktime]) * dV,real(inter_bands[0][1][ktime]) * dV,real(inter_bands[0][2][ktime]) * dV,real(inter_bands[0][3][ktime]) * dV,real(inter_bands[1][0][ktime]) * dV, real(inter_bands[1][1][ktime]) * dV,real(inter_bands[1][2][ktime]) * dV,real(inter_bands[1][3][ktime]) * dV,real(inter_bands[2][0][ktime]) * dV,real(inter_bands[2][1][ktime]) * dV,real(inter_bands[2][2][ktime]) * dV, real(inter_bands[2][3][ktime]) * dV);
            occup_temp = 0;

            for (int m = 0; m < Nband; ++m)
            {
                fprintf(occup_out, "%.16e    ", real(density_matrix_integrated[ktime][m*Nband + m]) / Npoints);
                occup_temp += real(density_matrix_integrated[ktime][m*Nband + m]);
            }
            fprintf(occup_out, "%.16e\n", occup_temp / Npoints);
        }
        fflush(interj_out);
        fflush(interj_out_spin);
        fflush(intraj_out);
        fflush(occup_out);
    }
    //fclose(  fout     );
    if (rank==MASTER)
    {

        fclose( interj_out );
        fclose( intraj_out );
        fclose( occup_out);
        fclose( dm1dout );

    }


    for (int itemp=0; itemp < Ngrad; itemp++)
    {
        free( inter_rad[itemp] );
        free( intra_rad[itemp] );
        for (int jtemp=0; jtemp < Ndim; jtemp++)
        {
          free( inter_rad_spin[itemp][jtemp] );
        }


    }
    delete[] blockcounts;
    delete[] displs;

    Delete2D<complex>(density_matrix_integrated, Nband*Nband, sbe->fpulses->Nt);


    //fclose(simulation_out);


    cout << "\n***********************\nEND OF THE PROGRAM with No. of cores = " << Nprocessors << "  uses j = " << jstart << ";  to  j = " << jend << "  with rank = "  << rank << "\n\n" <<  endl;
    cout << "=====================\n\n";

   cout << "Time taken: "<< (double)(clock() - tStart)/CLOCKS_PER_SEC ;
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
};


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
