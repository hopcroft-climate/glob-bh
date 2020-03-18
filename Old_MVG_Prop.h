	
   		
//==============================================================	
   			
   	/* for (int j=0; j<sample; j++)
    {
    	Tvec.assign(nP[j], 0.0);
    	tvec.assign(nP[j], 0.0);
    	for (int i=0; i<nP[j] ; i++)
	    {
	    		inPrT >> Tvec[i];	
	    		inPrt >> tvec[i];
	    }
    }
    for (int j=sample; j<(sample+1) ; j++)
    {
    	Tvec.assign(nP[j], 0.0);
    	tvec.assign(nP[j],0.0);
    	for (int i=0; i<nP[j] ; i++)
	    {
	    		inPrT >> Tvec[i];
	    		inPrt >> tvec[i];
	    }
	    inLLss >> ll;
	    inPRss >> prior;     
	    UKnew[j-sample].set_history(Tvec,nP[j]);
	    UKtnew[j-sample].set_history(tvec,nP[j]);
    }
    for (int j=sample + 1; j<MCMCiterations; j++)
    {
    	Tvec.assign(nP[j], 0.0);
    	tvec.assign(nP[j], 0.0);
    	for (int i=0; i<nP[j] ; i++)
	    {
	    		inPrT >> Tvec[i];	
	    		inPrt >> tvec[i];
	    }
    }
    */

     		

        	
   /* prop = false;
			if (nP[0] > 0)
		     {
			
				
				for (int i=1; i<(nP[0]-1); i++)
				{
				//	cerr << tvec[i] << "\t";
					if (tvec[i] == 0.0)
					{    
						prop = false;
						beta[number] = 1.0;
						i++;
					}
					else{
						prop = true;
						i++;
					}
				}
				//cerr << endl;
		    }
		    else{
				prop = false;
				beta[number] = 0.0;
			}
				
			cerr << "PROP "<< prop << endl;
			
			*/
			
			//=================================================================================================

	 	
/*double Evidence( double beta, int N, double ll, double prior)
{
	beta = beta/(N);
    beta = 1/beta;
    beta =log(beta;
    double Ev =  - beta;
    Ev +=  ll + prior;
    //now add the 
   
	return Ev;
}
*/

                	
                	
                	
                	
                	
                	
                	
                	
                	
	/*
                	
                	
                	
                	
                	
//=========================================================================================================
	
Array2D<double> MCMC_SS::correlations(Array2D<double> Tarray, int iterations, int MCMCiterations)
{
	Array2D<double> C = Array2D<double> (iterations,iterations,0.0);
	// compute mean and correlation matrix
	Array1D<double> Tm = Array1D<double> (iterations, 0.0);
	
	//ofstream outC("CorrelationMatrix.txt");
	
	for (int i=0; i<iterations; i++)
	{
		Tm[i] = 0.0;
		for (int j=0; j<MCMCiterations; j++)
		{
			Tm[i] += Tarray[i][j];
		}
		Tm[i] = Tm[i]/MCMCiterations;
	}
   
	for (int i=0; i<iterations; i++)
	{
		//ofstream outT("tarrya.txt");
		 for (int j=0; j<iterations; j++)
		 {
		 	for (int k=0; k<MCMCiterations; k++)
		 	{
		 		C[i][j] += (1/double(MCMCiterations-1))*(Tarray[j][k]-Tm[j])*(Tarray[j][k]-Tm[j]);
		 		//cerr <<(Tarray[j][k]-Tm[j])*(Tarray[j][k]-Tm[j]) <<"\t";
		 		//outT << Tarray[j][k] << "\t";
		 	}
		 	
		 }
		 
	}
	
    return C;	
}
//=====================================================================
double MCMC_SS::propProb( Array2D<double> C, int sampleNumber, 
      vector<double> Wmcmc, vector<double> tmcmc, int points)
{
	
	
	Array2D<double> inv = Array2D<double> (iterations, iterations, 0.0);
    //determinant of covariance matrix
    double det = determinant(C,iterations);
    cerr << "det " << det << endl;
    // covariance matrix inverse
    inv = invert(C , iterations);
    ofstream outIn("Inverse.txt");
    outIn << inv;
    outIn.close();
	
	
	double PP = 0.0;
	// calculate the posterior probability from the MVG distribution.
	Array1D<double> Tsurf = Array1D<double> (iterations, 0.0);
    RJ_SS cov;
	Tsurf = cov.Interpolate_SS( Wmcmc,  tmcmc, points );
    
    for (int i=0; i<iterations; i++)     
    {
    	for (int j=0; j<iterations; j++)
    	{
    			inv[i][j] = inv[i][j]/det;
    			PP += exp(-0.5* Tsurf[i] *inv[i][j]*  Tsurf[j]);
    	}
    }
    
	double mul = (2*3.1471);
	mul = 1/(det* pow(mul,iterations/2));
	PP = mul * PP;
	cerr << "here in pp, PP= " << PP << endl;
	return PP;
}

//=====================================================================
double determinant(Array2D<double> A, int a)
{
	JAMA::Cholesky<double> chol(A);

	if (chol.is_spd())
		A= chol.getL();
		
  	else
		cout << "factorization was not complete.\n";
		double det=0.0;
		for (int i=0; i<a; i++)
		{
			 det = A[i][i] * A[i][i];
		}
		return det;
}
//=====================================================================
Array2D<double> invert(Array2D<double> A,  int a)
{
	 Array2D<double> inv = Array2D<double> (a,a,0.0);
	 Array2D<double> l = Array2D<double> (a,a,0.0);
	 Array2D<double> u = Array2D<double> (a,a,0.0);
	 
	 Array1D<int> indx = Array1D<int> (a,0.0);
	 Array1D<double> col = Array1D<double> (a,0.0);
	 
	 LU( A, a,  l,  u);
	 
	 for (int j=0; j<a; j++)
	 {
	 	for (int i=0; i<a; i++) col[i] = 0.0;
	 	col[j] = 1.0;
	 	lubksb(A,indx,col,a);
	 	for (int i=0; i<a; i++) inv[i][j] = col[i];
	 }
	 return inv;
}
//=======================================================================
		*/
