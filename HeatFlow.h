int prop_hf(	Array1D<double> & HtFnew, Array1D<double>  HtF,  Array1D<double> & TEQnew, 
									 Array1D<double>  Teq, double alpha_H, double alpha_Teq, int nBH, int hu )
									 
		    
{		
		 GaussianRNG up;
		 
		 Array1D<double> alpha = Array1D<double> (2, 0.0);
		 alpha[0] = alpha_H;
		 alpha[1] = alpha_Teq;
		 Array2D<double> chol = Array2D<double> (2,2, 0.0);
		 Array1D<double> update = Array1D<double> (2, 0.0);
		 chol[0][0]  = 1.0;
		 chol[0][1] = 0.0;
		 //chol[1][1]=  0.2500;
		 //chol[1][0] = -0.9000 ;
		 chol[1][1]=  0.2800;
		 chol[1][0] = -0.9600 ;
	     
	     for (int i=0; i<2; i++)
	     {
	     	for (int j=0; j<2; j++)
	     	{
	     		update[i] += chol[i][j] * up.generate() * alpha[i];
	     	}
	     	
	     }
	     
	     HtFnew[hu] += update[0] ;
	     //cerr << "update[0] " << update[0] << endl;
	     TEQnew[hu]	  += update[1] ;
		 //cerr << HtF << "\t" << HtFnew << "\t" << Teq << "\t" << TEQnew << endl;	
		
		 if (TEQnew[hu] > 15.0 || TEQnew[hu] < 0.0)
		 {
		 	HtFnew[hu] -= update[0] ;
		 	TEQnew[hu] -=update[1];
		 }
		 if (HtFnew[hu]> 0.15 || HtFnew[hu] < 0.0150)
		 {
		 	HtFnew[hu] -= update[0] ;
		 	TEQnew[hu] -=update[1];
		 }
		 cerr << "HtF[hu]" <<HtF[hu] << "\t TEQ[hu] " << Teq[hu] << endl;
		 cerr << "HtFnew[hu]" <<HtFnew[hu] << "\t TEQnew[hu] " << TEQnew[hu] << endl;
		 
 return 0;		 
 
}


