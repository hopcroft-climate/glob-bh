// Tridiagonal solver Hirsch p506
// Solves Ax = B; with A = AA[k]X[k-1] + BB[k]x[k] + CC[k]x[k+1] = FF[k];
// and x is Temperature (in this case)
// Uses Thomas algorithm for forward and backward substitutions.  
#include <iostream>
#include <sstream>
//#include "tnt_array2d.h"
//#include "tnt_array2d_utils.h"
//#include "tnt_array1d.h"
//#include "tnt_array1d_utils.h"
using namespace std;
using namespace TNT;
TNT::Array1D<double> tsolver_bc(Array1D<double>BB, Array1D<double> plus, Array1D<double> minus, Array1D<double>f, int n, double Tsurface)
{    
     f[0] = Tsurface;
     f[1] = f[1] - Tsurface * plus[0]; 
        
     // For known surface Temperature (ie implementing b.c's) must replace the leading edges with 0's and element [0][0] with  1
        
     BB[0]= 1.0; 
     minus[0] = 0;
     plus[0] = 0;


     
    TNT::Array1D<double> FF;
    FF = TNT::Array1D<double> (n, 0.0);
    FF =f;
    
    TNT::Array1D<double>AA = Array1D<double> (n, 0.0);
    TNT::Array1D<double>CC = Array1D<double> (n, 0.0);
    
    for (int i =0; i<n-1; i++)
    {
        AA[i+1] = minus[i];
    }
   
    for (int i = 0; i<n-1; i++)
    {
        CC[i] = plus[i];
    }
   
    //Start substitutions:
    BB[0] = 1/BB[0];
    AA[0] = FF[0]*BB[0];         
    for ( int k=1; k<n; k++)
    {   
         CC[k-1] = CC[k-1] * BB[k-1];
         BB[k] = BB[k] -AA[k]*CC[k-1];
         BB[k] = 1/BB[k];
         AA[k] = (FF[k] -AA[k]*AA[k-1])* BB[k];
    }  
    FF[n-1] = AA[n-1];
    for (int i = n-1; i>0; i--)
    {
        FF[i] = AA[i] - CC[i]*FF[i+1];  
        
    }
    
    return FF;
}



