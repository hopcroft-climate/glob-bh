// local K -matrix calculator for assembly to Global K matrix
//#include <iostream>
//
//#include "tnt_array2d.h"
//#include "tnt_array2d_utils.h"
//using namespace std;

//______________________________________________________________________________
// local K -matrix calculator for assembly to Global K matrix
TNT::Array2D<double> local_K ( double A, bool surfaceTerm, double a)
{ 
      TNT::Array2D<double> surfMatrix = TNT::Array2D<double> (2, 2, 0.0);
      TNT::Array2D<double> k = TNT::Array2D<double> (2, 2, 0.0);
      
     // stringstream sout;
      //sout.str("  2   2    1    -1     -1     1 ") ;   // generic local element calc
      //sout>> k;
      //sout.clear();
       k[0][0] = 1;
      k[1][0] = -1;
      k[0][1] = -1;
      k[1][1] = 1;
      switch (surfaceTerm)
      {
      case false:    
            //cerr << " k is "  << k << endl;
            for (int i =0; i < 2; i++)
            {  
                 for(int j = 0; j<2; j++)
                 {
                       k[i][j]  = k[i][j] * A ;
                 }
            }      
            //cout << " In case 0 " << endl;
            break;       
      case true:   
           // sout.str(" 2    2     1     0     0    0 ") ;                    // for including the surface heat exchange
           // sout >> surfMatrix;
           // sout.clear();
            surfMatrix[0][0] = 1;
		surfMatrix[0][1] = 0;
		surfMatrix[1][0] = 0;
		surfMatrix[1][1] = 0;
            //cerr << " SurfaceMatrix is " << surfMatrix << endl;
            for (int i =0; i < 2; i++)
            {  
                for(int j = 0; j<2; j++)
                {
                    k[i][j]  = a*A*surfMatrix[i][j]; 
                }
            }   
            //cerr << " k is "  << k << endl;            
            // cout << " In case 1 " << endl;
            
            break;
      }     
    
      return k;
}
//______________________________________________________________________________
// local M -matrix calculator for assembly to Global M matrix

TNT::Array2D<double> local_M (double p,double c,  double A)
{ 
    TNT::Array2D<double> m;
    m = TNT::Array2D<double> (2, 2, 0.0);
    //stringstream sout;
    //sout.str("  2   2    2    1     1     2 ") ;                      // generic local element calc
    //sout>> m;
    //sout.clear();
    m[0][0] = 2;
    m[1][0] = 1;
    m[0][1] = 1;
    m[1][1] = 2;   

    for (int i =0; i < 2; i++)
    {  
         for(int j = 0; j<2; j++)
         {
               m[i][j]  = m[i][j]*p*c/(6);                          // (* dx/dx)
         }
    }      
    return m;
}


