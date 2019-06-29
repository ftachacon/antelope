//
//  momaxis.h
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
//

#ifndef momaxis_h
#define momaxis_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <complex>
#include "constant.h"

#define complex complex<double>
using namespace std;


class momaxis
{
    
public:
    
    double **k, tkx, tky;                             //kx, ky, and kz mesh
    
    double *dk;                             //grid or mesh steps
    
    double *kmaxs;                           //Number of points
    
    int *N;
    
    long long int Ntotal;
    
    double dV;
    
    
    momaxis( int *_N
            ,double *_kmaxs );                // Creator/Constructor Object ...
    
    ~momaxis( );                              // Destructor ...
    
    
    void basic_mem();                        // Setting memory ...
    
    void set_grid_mem( );                    // Setting memory ...
    
    void steps_sizes( int *_dx
                     ,double *_kmaxs );            // Momentum calculation
    
    
    void box( int dir_index );                // Basic axis or box ...
    void set_brillouin_zone_grid( );         // Creating axis ...
    
    
    void mom_outputs( FILE *output
                     , int skip0
                     , int skip1
                     , int skip2 );          //Momentum output
    
    
    long int index( int *i
                    ,int *j
                    ,int *l );
    
    
    void checker( int *_N );
    
};




//Constructor
momaxis::momaxis( int *_N, double *_kmaxs )
{
    
    steps_sizes( _N, _kmaxs );
    set_brillouin_zone_grid( );
    
}




//Destructor
momaxis::~momaxis()
{
    
    free(   dk   );
    free( kmaxs  );
    free(   N    );
    
    
    for (int i = 0; i < Ndim; i++)
    {
        
        free( k[i] );
        
    }
    
    free( k );
}




//Momentum steps and sizes
void  momaxis::steps_sizes(int *_N, double *_kmax )
{
    
    /*
     
     Conditional on kmaxs
     
     N[1]=N[2] has to be one to reduce the dimensionality
     
     of the problem to 1D,
     
     N[2]= 1 to reduce the dimentionality of the problem to 2D
     
     then, it follows:
     
     */
    
    checker( _N );
    
    basic_mem( );
    
    
    for ( int i=0; i<Ndim; i++ )
    {
        
        kmaxs[i]    = _kmax[i];
        N[i]        = _N[i];
        
    }
    
    
    //Building condition for symmetric grid
    for ( int i=0; i<Ndim; i++ )
        dk[i] = 2.*kmaxs[i] / double( N[i] -1 );
    
    
    
    //1D condition
    if( N[1]==1 && N[2]==1 )
    {
    
        dk[1] = 1.;
        dk[2] = 1.;
    
    }
    
    
    //2D condition
    if( N[2]==1 )
        dk[2] = 1.;
    
    
    Ntotal  = N[0]*N[1]*N[2];
    dV      = dk[0]*dk[1]*dk[2];
    
    set_grid_mem( );

    
}



void momaxis::checker( int *_N )
{
    if (_N[0] <= 0)
    {
        cout << "\n\n/****************OJO***************/\nNx should be larger or eq. to one, Nx >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
    
    if (_N[1] <= 0)
    {
        cout << "\n\n/***************OJO******************/\nNy should be larger or eq. to one, Ny >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
    if (_N[2] <= 0)
    {
        cout << "\n\n/*****************OJO***************/\nNz should be larger or eq. to one, Nz >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
}


//Creating basic memory
void momaxis::basic_mem()
{
    
    
    dk      = ( double* ) malloc( Ndim*sizeof( double ) );
    kmaxs   = ( double* ) malloc( Ndim*sizeof( double ) );
    N       = (    int* ) malloc( Ndim*sizeof(    int ) );
    
    
    memset( dk, 1, Ndim*sizeof( double ) );
    memset( kmaxs, 0, Ndim*sizeof( double ) );
    memset( N, 1, Ndim*sizeof( int ) );
    
    k   = ( double** ) malloc( Ndim*sizeof( double*) );
    
    
}




//Setting axis or box memory
void momaxis::set_grid_mem( )
{
    
    
    for (int i=0; i<Ndim; i++)
    {
        
        k[i] = ( double* )malloc( N[i] *sizeof( double ) );
        memset( k[i], 0, N[i]*sizeof( double ) );
        
    }
    

}




//Box or grid  along dir-index
void momaxis::box(int dir_index )
{
    
    
    if (N[dir_index]>1)
    {
        
    
        double kmin = -kmaxs[dir_index];
        
        
        for ( int i = 0; i < N[dir_index]; i++ )
            k[dir_index][i] = kmin + dk[dir_index]*i;
        
        
    }
    else
    {
        
        
        k[dir_index][0] = 0.;
        
        
    }
    
    //k[1][0] = -0.1   ;
    //k[1][0]=0.2;
}


//index calculation
long int momaxis::index( int *i, int *j, int *l )
{
    
    return N[1]*N[0]*(*l) + N[0]*(*j) + (*i);
    
}



//First Brillouine Zone
void momaxis::set_brillouin_zone_grid()
{

    
    //x-direction
    box( 0 );
    
    
    //y-direction
    box( 1 );
    
    
    //z-direction
    box( 2 );
    
    
}



//Outputs of grid
void momaxis::mom_outputs( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
    fprintf( output,"%d \n", skip0  );
    fprintf( output,"%d \n", N[0]  );
    fprintf( output,"%e \n\n", dk[0]  );
    
    
    
    for ( int i=0; i<N[0]/skip0; i++ )
        fprintf( output, "%e \n", k[0][ i*skip0 ] );
    
    
    fprintf( output,"\n%d \n", skip1  );
    fprintf( output,"%d \n", N[1]  );
    fprintf( output,"%e \n", dk[1]  );
    fprintf(output,"\n");
    
    
    for ( int i=0; i<N[1]/skip1; i++ )
        fprintf( output, "%e \n", k[1][ i*skip1 ] );
    
    
    fprintf( output,"\n%d \n", skip2  );
    fprintf( output,"%d \n", N[2]  );
    fprintf( output,"%e \n", dk[2]  );
    fprintf(output,"\n");
    
    
    for ( int i=0; i<N[2]/skip2; i++ )
        fprintf( output, "%e \n", k[2][ i*skip2 ] );
    
    
    fflush( output );
    
    
}

#endif /* momaxis_h */
