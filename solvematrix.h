#include <iostream>
#include <cmath>
//#include "nr.h"


    void ludcmp(Array2D<double> &a, Array1D<double> &indx, double &d, int n)
    {
          const double TINY =1E-20;
          int i= 0, imax= 0, j =0, k=0;
          double big = 0, dum =0, sum= 0, temp= 0;
//          int n= a.nrows();
          Array1D<double> vv = Array1D<double> (n, 0.0);                   // vv stores the implicit scaling of each row
          d=1.0;                            // no row interchange yet
    
    for(i=0; i<n; i++)                    // loop over rows to get implicit scaling information
    {
    	
        big = 0.0;
        for (j=0;j<n; j++){
                if ((temp=fabs(a[i][j])) > big) big = temp;
            }    
        if (big == 0.0) cerr << "Singular matrix in routine LUdcmp";
        
        vv[i] = 1.0/big;
    }
    for (j=0; j<n; j++)                              //Loop over Crout's method
    {
    	
        for (i=0; i<j; i++)                            // Eqn 2.3.12 except i=j
        {
            sum= a[i][j];
            for (k=0; k<i; k++) sum -=a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;                                    // Initialize search for largest pivot element
        for (i=j; i<n; i++) {
            sum = a[i][j];                        // i=j of eqn 2.3.12 and i=j+1...N-1 of 2.3.13
            for(k=0;k<j;k++) sum -=a[i][j]*a[k][j];
            a[i][j] = sum;
            if ((dum=vv[i]*fabs(sum))>=big){
                big=dum;
                imax = i;
            }
        }
        if (j !=imax)                             // do we need to interchange rows?
        {
            for (k=0; k<n; k++)                        // Yes, do so...
            {                
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] =dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        // If the pivot element is zero the matrix is singular ( at least to the precision of the algorithm) for some
        // applications on singular matrices it is better to set TINY = 0.0;
        if (j !=n-1)
        {
            dum = 1.0/(a[j][j]);
            for (i=j+1; i<n; i++) a[i][j] *=dum;
        }
    }    
                 // Go back for the next column in the reduction
}

//==========================================================================================
// Crout's algorithm
int LU(Array2D<double> A, int n, Array2D<double>& alpha, Array2D<double>& beta)
{
	      //cerr << " A " << A << endl;
       double sum = 0.0;
       int iter= 0;
       
       for (int i=0; i<n; i++)
       {
           alpha[i][i] = 1.0;
       }
       
       for (int j=0; j<n;j++)
       {
           for (int i=0; i<=j; i++)
           {
               sum = 0;
               iter  = i-1;
               if (i==0) {iter  =0;}
              
               for (int k=0; k<=(iter); k++)
               {
                   sum += alpha[i][k] * beta[k][j];
                   
               }
              
                beta[i][j] = A[i][j] -sum; 
           }
           
           for (int ii=j+1; ii<n; ii++)
           {
              
               sum = 0.0;
               for (int k=0; k<=j-1; k++)
               {
                   sum += alpha[ii][k]*beta[k][j];
                   
               }
               alpha[ii][j] = (1/beta[j][j]) * (A[ii][j]  - sum);

              
           }
           
           
       }
           
       return 0;
}
//=================================================================================================================



// back substitution for LU decmp



   void lubksb(Array2D<double> &a, Array1D<int>  &indx, Array1D<double> &b, int iterations)
   {
      int i , ii=0, ip, j;
      double sum;
      
      int n = iterations;
      for (i=0; i<n;i++)                  // When ii set positive value, it becomes index of 1st non vanishing element of b
      {                                    // We no do forward substitution using eqn 2.3.6 mush unscramble the
        ip = indx[i];                    // permutation as we go.
        sum = b[ip];
        b[ip]= b[i];
        if(ii !=0)
                for (j=ii-1;j<i; j++)  sum -= a[i][j]*b[j];
        else if (sum !=0.0)                 // A nonzero element was encountered, so from now on we will have to do the sums
                ii = i+1;                   // in the loop above
        b[i] = sum;
      }
      for (i=n-1;i>=0;i--)                  // Now we do backsubstitution as in eqn 2.3.7
      {
        sum = b[i];
        for (j=i+1; j<n; j++) sum -=a[i][j]*b[j];
        b[i] = sum/a[i][j];                // store a solution of the vector X
      }
   }

