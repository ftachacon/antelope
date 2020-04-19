// utility.h
// Consider using boost library instead of using this header file.

#ifndef utility_h
#define utility_h

#include <string>

const std::string whiteSpaces( " \f\n\r\t\v" );

// ========================================
// Utility realted to string 
// Might be replaced by using boost/algorithm/string/trim.hpp

void trim_left( std::string& str,
      const std::string& trimChars = whiteSpaces )
{
   std::string::size_type pos = str.find_first_not_of( trimChars );
   str.erase( 0, pos );
}

void trim_right( std::string& str,
      const std::string& trimChars = whiteSpaces )
{
   std::string::size_type pos = str.find_last_not_of( trimChars );
   str.erase( pos + 1 );    
}

void trim( std::string& str, const std::string& trimChars = whiteSpaces )
{
   trim_left( str, trimChars );
   trim_right( str, trimChars );
}

// Distribute jobs through processors
// divide domain offset <= x < totalLength + offset
// Note that if mpi_size > totalLength && mpi_rank >= totalLength, 
// then istart = totalLength, iend = totalLength
void ParaRange(int totalLength, int offset, int mpi_size, int mpi_rank, int *istart, int *iend)
{
    int iwork1, iwork2;

    iwork1 = totalLength / mpi_size;
    iwork2 = totalLength % mpi_size;

    *istart = iwork1 * mpi_rank + min(iwork2, mpi_rank) + offset;
    *iend = *istart + iwork1;
    if (iwork2 > mpi_rank)
        *iend = *iend + 1;
}


// Replace below functions with boost multiarray feature later
template <class T> T **Create2D(int N1, int N2)
{
	T *arr1d = new T[N1*N2];
	for (int i = 0; i < N1*N2; ++i)
		arr1d[i] = 0;
	T **arr = new T*[N1];
	for (int i0 = 0; i0 < N1; ++i0)
	{
		arr[i0] =  arr1d
				+ i0 * N2;
	}
	return arr;
};
template <class T> void Delete2D(T **arr, int N1, int N2)
{
	delete[] &(arr[0][0]);
	delete[] arr;
};
template <class T> T ***Create3D(int N1, int N2, int N3)
{
	T *arr1d = new T[N1*N2*N3];
	for (int i = 0; i < N1*N2*N3; ++i)
		arr1d[i] = 0;
	T ***arr = new T**[N1];
	for (int i0 = 0; i0 < N1; ++i0)
	{
		arr[i0] = new T*[N2];
		for (int i1 = 0; i1 < N2; ++i1)
		{
			arr[i0][i1] = arr1d
				+ i1 * N3
				+ i0 * N3 * N2;
		}
	}
	return arr;
};
template <class T> void Delete3D(T ***arr, int N1, int N2, int N3)
{
	delete[] &(arr[0][0][0]);
	for (int i = 0; i < N1; ++i)
	{
		delete[] arr[i];
	}
	delete[] arr;
};

template <class T> T ****Create4D(int N1, int N2, int N3, int N4)
{
	T *arr1d = new T[N1*N2*N3*N4];
	for (int i = 0; i < N1*N2*N3*N4; ++i)
		arr1d[i] = 0;
	T ****arr = new T***[N1];
	for (int i = 0; i < N1; ++i)
	{
		arr[i] = new T**[N2];
		for (int j = 0; j < N2; ++j)
		{
			arr[i][j] = new T*[N3];
			for (int k = 0; k < N3; ++k)
			{
				arr[i][j][k]  = arr1d
					+ k * N4
					+ j * N4 * N3
					+ i * N4 * N3 * N2;
			}
		}
	}
	return arr;
};
template <class T> void Delete4D(T ****arr, int N1, int N2, int N3, int N4)
{
	delete[] &(arr[0][0][0][0]);
	for (int i = 0; i < N1; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			delete[] arr[i][j];
		}
		delete[] arr[i];
	}
	delete[] arr;
};
template <class T> T *****Create5D(int N1, int N2, int N3, int N4, int N5)
{
	T *arr1d = new T[N1*N2*N3*N4*N5];
	for (int i = 0; i < N1*N2*N3*N4*N5; ++i)
		arr1d[i] = 0;
	T *****arr = new T****[N1];
	for (int i = 0; i < N1; ++i)
	{
		arr[i] = new T***[N2];
		for (int j = 0; j < N2; ++j)
		{
			arr[i][j] = new T**[N3];
			for (int k = 0; k < N3; ++k)
			{
				arr[i][j][k] = new T*[N4];
				for (int l = 0; l < N4; ++l)
				{
					arr[i][j][k][l] = arr1d
						+ l * N5
						+ k * N5 * N4
						+ j * N5 * N4 * N3
						+ i * N5 * N4 * N3 * N2;
				}
			}
		}
	}
	return arr;
};
template <class T> void Delete5D(T *****arr, int N1, int N2, int N3, int N4, int N5)
{
	delete[] &(arr[0][0][0][0][0]);
	for (int i = 0; i < N1; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			for (int k = 0; k < N3; ++k)
			{
				delete[] arr[i][j][k];
			}
			delete[] arr[i][j];
		}
		delete[] arr[i];
	}
	delete[] arr;
};
template <class T> T ******Create6D(int N1, int N2, int N3, int N4, int N5, int N6)
{
	T *arr1d = new T[N1*N2*N3*N4*N5*N6];
	for (int i = 0; i < N1*N2*N3*N4*N5*N6; ++i)
		arr1d[i] = 0;
	T ******arr = new T*****[N1];
	for (int i = 0; i < N1; ++i)
	{
		arr[i] = new T****[N2];
		for (int j = 0; j < N2; ++j)
		{
			arr[i][j] = new T***[N3];
			for (int k = 0; k < N3; ++k)
			{
				arr[i][j][k] = new T**[N4];
				for (int l = 0; l < N4; ++l)
				{
					arr[i][j][k][l] = new T*[N5];
					for (int m = 0; m < N5; ++m)
					{
						arr[i][j][k][l][m] = arr1d
							+ m * N6
							+ l * N6* N5
							+ k * N6* N5 * N4
							+ j * N6* N5 * N4 * N3
							+ i * N6* N5 * N4 * N3 * N2;
					}
				}
			}
		}
	}
	return arr;
};
template <class T> void Delete6D(T ******arr, int N1, int N2, int N3, int N4, int N5, int N6)
{
	delete[] &(arr[0][0][0][0][0][0]);
	for (int i = 0; i < N1; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			for (int k = 0; k < N3; ++k)
			{
				for (int l = 0; l < N4; ++l)
				{
					delete[] arr[i][j][k][l];
				}
				delete[] arr[i][j][k];
			}
			delete[] arr[i][j];
		}
		delete[] arr[i];
	}
	delete[] arr;
};

#endif