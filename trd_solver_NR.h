// Num Recipes in C++ pg 54.
// 12/05/05
#include <iostream>
#include <sstream>
//#include "tnt_array2d.h"
//#include "tnt_array2d_utils.h"
//#include "tnt_array1d.h"
//#include "tnt_array1d_utils.h"
using namespace std;
using namespace TNT;
TNT::Array1D<double> tsolver_NR(Array1D<double>b, Array1D<double> A, Array1D<double> C, Array1D<double>r, int n)
{   
    TNT::Array1D<double> gam = Array1D<double> (n, 0.0);
    TNT::Array1D<double> u = Array1D<double> (n, 0.0);
    TNT::Array1D<double> c = Array1D<double> (n-1, 0.0);
    TNT::Array1D<double> a = Array1D<double> (n, 0.0);
    for (int i = 0; i<n-1; i++)
    {
        c[i] = C[i];
    }
    for (int i = 1; i<n; i++)
    {
        a[i] = A[i-1];
    }
    
    double bet = 0.0;
    if (b[0] == 0.0) cerr << " Error1 in tridiag" << endl;
    
    u[0] = r[0]/(bet = b[0]);
    
    for (int j =1; j<n; j++)
    {
        gam[j] = c[j-1]/bet;
        bet = b[j] - a[j]*gam[j];
        
        if (bet == 0.0) cerr << " Error 2 in tridiag " << endl;
        u[j] = (r[j] - a[j] * u[j-1])/bet;
    }
    for (int j= (n-2); j>=0; j--)
    {
        u[j] -= gam[j+1] * u[j+1];           // Backsubstitution
    }
    
    return u ;
}



