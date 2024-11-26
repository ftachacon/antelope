/// class for momentum grid 
/**
 * This class get axes and start point of Brillouin zone and generate 
 * Axes can be non-orthogonal. K-points are equally spaces.
 * Now first and last point are not same (step is divided by (N), not (N-1) )
 * @todo Fix integration method - simpson, write new momoutput 
 * @author Alexis Agustín  Chacón Salazar
 * @author Dasol Kim
 */


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
#include <vector>
#include <numeric>
#include "constant.h"
#include "utility.h"


class momaxis
{
    
public:
    
    std::array<double, Ndim> *kgrid; 
    std::array<double, Ndim*Ndim> BzAxes;
    std::array<double, Ndim> BzOrigin;
    
    //double *dk;                             //grid or mesh steps
    double *weight;
    double total_weight;
    
    std::array<int, Ndim> N;
    
    // Now internal (index 0 <= < Npartial) corresponds to global index (idx_start <= < idx_end)
    int Ntotal; // Total number of grid points
    int Npartial; // Number of grid points belonging to this process
    int idx_start, idx_end; // Starting and ending index of grid points belonging to this process
    
    double Volume, dV;
    std::string imethod;
    
    momaxis( std::array<int, Ndim> _N ,std::array<double, Ndim*Ndim> _axisvec, std::array<double, Ndim> _origink, int _idx_start=-1, int _idx_end=-1 )
        : N(_N), BzAxes(_axisvec), BzOrigin(_origink)
    {
        idx_start = _idx_start; idx_end = _idx_end;
        Npartial = idx_end - idx_start;

        steps_sizes( _N, _axisvec, _origink );
        set_brillouin_zone_grid( );
        //imethod = "Trapz";
    };
    
    ~momaxis( );                              // Destructor ...

    void basic_mem();                        // Setting memory ...

    void steps_sizes( std::array<int, Ndim> _N ,std::array<double, Ndim*Ndim> _axisvec, std::array<double, Ndim> _origink );            // Momentum calculation

    void set_brillouin_zone_grid( );         // Creating axis ...

    void mom_outputs( FILE *output
                     , int skip0
                     , int skip1
                     , int skip2 );          //Momentum output
    
    
    int index( const int i ,const int j ,const int k )
    {
        return N[1]*N[2]*i + N[2]*j + k;
    };
    int index( const std::array<int, Ndim> _indexArr )
    {
        return N[1]*N[2]*_indexArr[0] + N[2]*_indexArr[1] + _indexArr[2];
    };
    std::array<int, Ndim> index( int _kindex)
    {
        return { _kindex/(N[1]*N[2]), (_kindex/N[2])%N[1], _kindex%N[2] };
    };
 
    //void Simpson();
    //void integral_method( const string NMethod );
    void checker( );

    void print_info();                       // Printing information of momentum grid
    
};


/*
void momaxis::integral_method( const string NMethod )
{
    
    imethod = NMethod;
    if (imethod == "Simpson")
    {
        Simpson( );
    }
    else
    {
        fill(weight, weight + Ntotal, 1.);
    }  
}*/


//Destructor
momaxis::~momaxis()
{
    delete[] weight;
    delete[] kgrid;
}







void momaxis::checker(  )
{
    if ( N[0] <= 0 || N[1] <= 0 || N[2] <= 0)
    {
        std::cerr << "Number of grid is negative or zero. check N = (" << N[0] << ", " << N[1] << ", " << N[2] << ")\n";
        std::exit(EXIT_FAILURE);
    }
    if (idx_start < 0 || idx_end < 0 || idx_start >= Ntotal || idx_end > Ntotal || idx_start > idx_end)
    {
        std::cerr << "Index is out of range. check idx_start = " << idx_start << ", idx_end = " << idx_end << "\n";
        std::exit(EXIT_FAILURE);
    }
}


//Creating basic memory
void momaxis::basic_mem()
{
    kgrid = new std::array<double, Ndim>[Npartial];

    weight = new double[Npartial];
    std::fill(weight, weight + Npartial, 1.); // This default value acts as trapzoid

    double temp_weight;
    int kindex = 0;
    for (int i = 0; i < N[0]; ++i)
    {
        for (int j = 0; j < N[1]; ++j)
        {
            for (int k = 0; k < N[2]; ++k)
            {
                temp_weight = 1.0;
                if (N[0] != 1 && (i ==0 || i == N[0]-1))
                {
                    temp_weight /= 2;
                }
                if (N[1] != 1 && (j == 0 || j == N[1]-1))
                {
                    temp_weight /= 2;
                }
                if (N[2] != 1 && (k == 0 || k == N[2]-1))
                {
                    temp_weight /= 2;
                }
                kindex = index(i, j, k) - idx_start;
                if (0 <= kindex && kindex < Npartial)
                    weight[ kindex ] = temp_weight;
            }
        }
    }
    total_weight = 0.;
    for (int i = 0; i < Npartial; ++i)
    {
        total_weight += weight[i];
    }
}




//Momentum steps and sizes
void  momaxis::steps_sizes( std::array<int, Ndim> _N ,std::array<double, Ndim*Ndim> _axisvec, std::array<double, Ndim> _origink )
{
    Ntotal = N[0] * N[1] * N[2];

    // if partial is not set, set Npartial to Ntotal
    if (idx_start == -1)
    {
        Npartial = Ntotal;
        idx_start = 0;
        idx_end = Ntotal;
    }
    
    checker( );
    
    basic_mem( );
    
    std::vector<std::array<double, Ndim> > axesList;
    for (int i = 0; i < Ndim; ++i)
    {
        if (N[i] != 1)
        {
            axesList.push_back({BzAxes[Ndim*i + 0], BzAxes[Ndim*i + 1], BzAxes[Ndim*i + 2]});
        }
    }
    std::array<double, Ndim> tempArray;
    switch (axesList.size())
    {
    case 1: // length of first vector
        Volume = sqrt( std::inner_product(axesList[0].begin(), axesList[0].end(), axesList[0].begin(), 0.0) );
        break;
    case 2: // length of cross product between two vectors
        tempArray = CrossProduct( axesList[0], axesList[1] );
        Volume = sqrt( std::inner_product(tempArray.begin(), tempArray.end(), tempArray.begin(), 0.0) );
        break;
    case 3: // triple product
        tempArray = CrossProduct( axesList[0], axesList[1] );
        Volume = sqrt( std::inner_product(tempArray.begin(), tempArray.end(), axesList[2].begin(), 0.0) );
        break;
    default: // single point
        // cerr << "Something is wrong while calculating Brillouin zone volume\n";
        // exit(EXIT_FAILURE);
        std::cerr << "Warning: You are using single point condition.\n";
        Volume = 1.;
        break;
    }
    
    dV = abs(Volume) / Ntotal;
    
}

/*
void momaxis::Simpson()
{
    
    int i_aux, j_aux;
    int ia0 = 0;
    int ie0 = N[0]-1;
    
    int ja0 = 0;
    int je0 = N[1]-1;
    
    int le=0;
    
    weight[ index( &ia0, &ja0, &le) ] = 1./9.;
    weight[ index( &ie0, &ja0, &le) ] = 1./9.;
    
    weight[ index( &ia0, &je0, &le) ] = 1./9.;
    weight[ index( &ie0, &je0, &le) ] = 1./9.;
    
    
    for ( int i = 1; i < int( (N[0]+1)/2.); i++ )
    {
        

        i_aux=2*i-1;
        
        weight[ index( &i_aux, &ja0, &le) ] = 4./9.;
        weight[ index( &i_aux, &je0, &le) ] = 4./9.;
        
    }
    
    
    for ( int i=1; i < int( N[0]/2.); i++ )
    {
        i_aux=2*i;
        weight[ index( &i_aux, &ja0, &le) ] = 2./9.;
        weight[ index( &i_aux, &je0, &le) ] = 2./9.;
    }
    
    
    for ( int j = 1; j < int( (N[1]+1) /2. ); j++ )
    {
        
        j_aux = 2*j-1;
        weight[ index( &ia0, &j_aux, &le) ] = 4./9.;
        weight[ index( &ie0, &j_aux, &le) ] = 4./9.;
        
    }
    
    for ( int j = 1; j < int( N[1]/2. ); j++ )
    {
        
        j_aux = 2*j;
        weight[ index( &ia0, &j_aux, &le) ] = 2./9.;
        weight[ index( &ie0, &j_aux, &le) ] = 2./9.;
        
    }

    
    
    
    for ( int j=1; j < int( (N[1]+1)/2. ); j++ )
        for ( int i=1; i < int( (N[0]+1)/2. ); i++ )
        {
            
            i_aux=2*i-1;
            j_aux=2*j-1;
            
            weight[ index( &i_aux, &j_aux, &le) ] = 16./9.;
            
        }
    
    for ( int j=1; j < int( N[1]/2. ); j++ )
        for ( int i=1; i < int( (N[0]+1)/2. ); i++ )
        {
            j_aux=2*j;
            i_aux=2*i-1;
            weight[ index( &i_aux, &j_aux, &le) ] = 8./9.;
        }
    
    for ( int j=1; j < int( (N[1]+1)/2. ); j++ )
        for ( int i=1; i < int( N[0]/2. ); i++ )
        {

            i_aux=2*i;
            j_aux=2*j-1;
            weight[ index( &i_aux, &j_aux, &le) ] = 8./9.;
        }
    

    
    
    for ( int j=1; j < int( N[1]/2. ); j++ )
        for ( int i=1; i < int( N[0]/2. ); i++ )
        {
            i_aux=2*i;
            j_aux=2*j;
            weight[ index( &i_aux, &j_aux, &le) ] = 4./9.;
        }
    
    
}*/




//First Brillouine Zone
void momaxis::set_brillouin_zone_grid()
{
    double ifrac, jfrac, kfrac;
    int kindex;
    for (int i = 0; i < N[0]; ++i)
    {
        if (N[0] == 1)
            ifrac = 0.;
        else
            ifrac = static_cast<double>(i)/(N[0]-1) - 0.5;
        for (int j = 0; j < N[1]; ++j)
        {
            if (N[1] == 1)
                jfrac = 0.;
            else
                jfrac = static_cast<double>(j)/(N[1]-1) - 0.5;
            for (int k = 0; k < N[2]; ++k)
            {
                if (N[2] == 1)
                    kfrac = 0.;
                else
                    kfrac = static_cast<double>(k)/(N[2]-1) - 0.5;
                // Now internal (index 0 <= < Npartial) corresponds to global index (idx_start <= < idx_end)
                kindex = index(i, j, k) - idx_start;
                if (0 <= kindex && kindex < Npartial)
                    kgrid[kindex] = {BzOrigin[0] + ifrac*BzAxes[0*Ndim + 0] + jfrac*BzAxes[1*Ndim + 0] + kfrac*BzAxes[2*Ndim + 0], 
                                        BzOrigin[1] + ifrac*BzAxes[0*Ndim + 1] + jfrac*BzAxes[1*Ndim + 1] + kfrac*BzAxes[2*Ndim + 1], 
                                        BzOrigin[2] + ifrac*BzAxes[0*Ndim + 2] + jfrac*BzAxes[1*Ndim + 2] + kfrac*BzAxes[2*Ndim + 2] };
            }
        }
    }
    
}



//Outputs of grid
void momaxis::mom_outputs( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
}

void momaxis::print_info()
{
    std::cout << "\n==========================================\n";
    std::cout << "Brillouin zone volume = " << Volume <<  " a.u.\n";
    std::cout << "\n==========================================\n";
    //cout << "Second K' point = ("<< Kpoint2[0] << " , " << Kpoint2[1] << " ) a.u.\n";
    //cout << "Delta-K' point = (" << DeltaKp[0] << " , " << DeltaKp[1] << " ) a.u.\n========================\n\n";
}

#endif /* momaxis_h */
