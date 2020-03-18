 void sample_dimensions(Array1D<int> k_sample,int k_sample_length, int& k, double& k_prob);

void times_histomatrix(vector<double> times, int points, Array1D<double>&  time_histomatrix, 
				vector<double> & times_posterior, int & tot);

  void Temps_histomatrix(vector<double> times, vector<double> Temp, int points, 
            Array2D<double> & Temp_hist,Array2D<double> & Temp_posterior, int& Tot, Array1D<double> Tsurf, int j);

void sample_Temps(vector<double>& Temp, Array2D<double> Temp_hist, int k,Array2D<double>   Temp_posterior, 
				double& Temp_prob,vector<int> locats, int MCMCiterations );
 
 int sample_times(vector<double> times, vector<double> times_posterior, vector<double> TIMES_OLD);

double abs(double a,double b);
 
 
 void  GST_Proposal (vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew,
                vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta)
{
       int MCMCiterations = 100;      
       vector<double> Tvec;
       vector<double> tvec;
       Array2D<double> Tarray = Array2D<double> (iterations, 100, 0.0);
       char filename[nBH+11];
	 
	 
       strcpy(filename, "Samples/");
       for (int i=8; i<nBH+8; i++)
       {
       	  filename[i] = (char) (GST_id[i-8]+48);
       }
       filename[nBH+8] =(char) 49;
       filename[nBH+9] = (char) GST_id[number]+48;
       filename[nBH+10] = (char) 0;
       // cerr << filename << endl;
       ifstream inPrT(filename);
       filename[nBH+8] =(char) 48;
       //  cerr << filename << endl;
       ifstream inPrt(filename);
       filename[nBH+8] =(char) 83;
       //  cerr << filename << endl;
       ifstream inPrS(filename);
       //filename[nBH+8] = (char) 76;
       //ifstream ink(filename);
	   
	  
       Array1D<double> nPd =Array1D<double> (MCMCiterations,0.0);
       Array1D<int> nP =Array1D<int> (MCMCiterations,0);
    
       for (int j1=0; j1<MCMCiterations; j1++)
       {
    	  inPrS >> nPd[j1];
          nP[j1] = int(nPd[j1]);
       }
       inPrS.close();      
                
       Array1D<int> k_sample = Array1D<int> (MCMCiterations, 0);
       for (int i=0; i<MCMCiterations; i++)
       {
    	   k_sample[i] = nP[i];
    	   	   
       }
       if (k_sample[MCMCiterations-1] > 0) 
       {
    	    prop = true;     
	        Array2D<double> Temp_hist = Array2D<double> (iterations, 100, 0.0);
     	    Array2D<double>  Temp_posterior = Array2D<double> (iterations, MCMCiterations, 0.0);
		    
		    // sample dimensions:----------------------------------
		    //-----------------------------------------------------
		    int k=0;
		    double k_prob = 0.0;
		    sample_dimensions(k_sample,MCMCiterations,  k,  k_prob);
		    
		    cerr  <<" k is set to "<< k << endl;
		    //sample times:----------------------------------------
		    //-----------------------------------------------------
		    Array1D<double> times_hist = Array1D<double> (iterations, 0.0);
		    vector<double> times_posterior;
		    vector<double> times;
		    vector<int> INDEX_OLD;
		    times.assign(k, 0.0);
		    int tot = 0;
		    int Tot =0;
		    Array1D<double> Tsurf = Array1D<double> (iterations, 0.0);
		    for (int jj=0; jj<MCMCiterations; jj++)
		    {
		    	Tvec.assign(k_sample[jj], 0.0);
		    	tvec.assign(k_sample[jj], 0.0);
		    	for (int i=0; i<k_sample[jj] ; i++)
			    {
			    		inPrT >> Tvec[i];	
			    		inPrt >> tvec[i];
			    }
			    times_histomatrix(tvec, k_sample[jj], times_hist,times_posterior, tot);
			    // Interpolate the Tvec first
			    RJ_SS Tpost;
			    Tsurf = Tpost.Interpolate_SS(Tvec, tvec, k_sample[jj] );
			    
			    Temps_histomatrix(tvec, Tvec,   nP[jj],  Temp_hist, Temp_posterior, Tot, Tsurf, jj);
		    }
		    inPrT.close();
		    inPrt.close();
		    for (int i=0; i<iterations; i++) 
		    {
		        //normalise the time values pdf
		    	times_hist[i] = double( times_hist[i])/double(tot);
		        for(int j=0; j<100; j++)
		        {
		                // normalise the Temperature/time pdf
		        	    Temp_hist[i][j] =double( Temp_hist[i][j])/double(Tot);
		        }
		    }
		    cerr << "total is "  << Tot <<endl;
		    
		    // select (k-2) times:
		    times[0] = 756.0;
		    times[k-1] = 0.0;
		    double time_prob = 1.0;
		    vector<int> locats;
			int INDEX =0;
			vector<double> times_old;
			time_prob =1.0;
		    if (k>2)
		    {
		    	 
		    	 for (int i=0; i<(k-2); i++)
		    	 {
		    	    INDEX = sample_times(times, times_posterior, times_old);
		    	 	//INDEX_OLD.push_back(INDEX);		    	 	
		    	 	times[i+1] = times_posterior[INDEX];
		    	 	times_old.push_back(times[i+1]); 	
		    	 	 	
		    	 	for (int z = 0; z<iterations; z++)
		    	    {
		    	 	   if (times[i+1] >= (z*12) && times[i+1] < ((z+1)*12))  // 12 is (dt/(sec/ year))
		    	 	   {
		   	  			  time_prob = time_prob* times_hist[z];
		   	  			  locats.push_back(z);
		   	  			  //locats_check(locats,z, times, times_posterior, INDEX_OLD);
		    	 	    }
		    	  	}
		    	 }
		    }
		    cerr << "size locats here " << locats.size() << endl;
		    // need to sort the time-vector in descending values 
		    sort(times.begin(), times.end());
		    reverse(times.begin(), times.end());
		    // need to also sort the locations vector (int) descending values
		    if (locats.size() > 1)
		    {
		       sort(locats.begin(),locats.end());
		       reverse(locats.begin(), locats.end());
		    }
		    
		    for (int i=0; i<k; i++)
		    {
		    	cerr <<"times" <<   times[i] << "\t ";
		    }
		    cerr << endl;
		    for (int i=0; i<(k-2); i++)
		    {
		       cerr <<"locats" <<   locats[i] << "\t ";
		    }
		    cerr << endl;
		    // set the proposal to the sample time-vector of length k.
		    UKtnew[0].set_history(times,k);
			
		    // sample temperatures---------------------------------
		    //-----------------------------------------------------
		    double Temp_prob = 1.0;
		    
		    sample_Temps(Tvec,Temp_hist, k,  Temp_posterior,Temp_prob, locats, MCMCiterations);
		    
		    //cerr << "k " << k << endl;
		    UKnew[0].set_history(Tvec,k);
		    for (int i=0; i<k; i++)
		    {
		    	cerr << times[i] << "\t";
		    	cerr << Tvec[i] << "\t";
		    }
		    cerr << endl;
		    beta[number] = time_prob * k_prob* Temp_prob;
		    cerr <<"p(time) " << time_prob << "\t p(k) " << k_prob << "\t p(T) " << Temp_prob << endl;
		    cerr << "beta(number) " << beta[number] << endl;
		    //===================
		    
		    
		
	}
	else
	{
		prop = false;
		beta[number] = 1.0;
		
	}
	

 }
 
//=====================================================================
//=====================================================================


void sample_dimensions(Array1D<int> k_sample,int k_sample_length, int& k, double& k_prob)
 {
 	Array1D<double> k_hist = Array1D<double> (20 , 0.0);
 		
 	for (int i=0 ; i<k_sample_length; i++)
 	{
 		//cerr << k_sample[i] << "\t";
 		for (int j=0; j< 20; j++)
 		{
 			if (k_sample[i] == j) k_hist[j] ++;   //increment pdf for dimensions
 		}
 	}
 	
 	UniformRNG choose;
 	//choose a dimension from the posterior sample
 	int index =  int(choose.generate()*double(k_sample_length));   
 	//------------------------------------------------------------
 	//cerr << "k_sample_length" << k_sample_length << endl;
 	//cerr << k_sample ;
 	//cerr << index << endl;
 	//------------------------------------------------------------
 	
 	k = int(k_sample[index]);
 	
 	//calculated the (normalised) probability of selecting this k value.
 	k_prob = double(k_hist[k])/double(k_sample_length); // normalise the pdf for dimensions
 	
 	    
 	    
 }
 
//=====================================================================
   
   void times_histomatrix(vector<double> times, int points, Array1D<double>&  time_histomatrix, 
                vector<double> & times_posterior, int & tot)
   {
   	
   	  
     for (int i=0; i<iterations;i++) 
     {
     	for (int j=1; j<points-1; j++)// don't sample the end points.
     	{
     		if (times[j] >= i*12 && times[j] < (i+1)*12)
     		{
     		  // create a time values pdf	
     		  time_histomatrix[i] ++;
     		  tot ++;
     		}
     		// create a the posterior time values sample
     		times_posterior.push_back(times[j]);
     	}
     }
   	  
   	 
   }

//=====================================================================
void Temps_histomatrix(vector<double> times, vector<double> Temp, int points, 
            Array2D<double> & Temp_hist,Array2D<double> & Temp_posterior, int& Tot, Array1D<double> Tsurf, int j)
{
	for (int i=0; i<iterations; i++)
	{
		for (int k=0; k<100; k++)   //given the 100 Temperature discretisations
		{
			if (Tsurf[i] >=(double(k)/10 -5)&& Tsurf[i] <(double(k+1)/10 -5))			//gives range [-5,5]
			{
				Temp_hist[i][k] ++;  //calculate the Temp/time pdf
				Tot++;
				Temp_posterior[i][j] += Tsurf[i];
			}
         }	
	}
}
//==================================================================
void sample_Temps(vector<double>& Temp, Array2D<double> Temp_hist, int k,Array2D<double> Temp_posterior, 
                  double& Temp_prob,vector<int> locats, int MCMCiterations )
{
	Temp_prob = 1.0;
	Temp.resize(k, 0.0);
	cerr << "size locats in sample Temps" << locats.size() << endl;
	vector<double> temporary;
    reverse(locats.begin(), locats.end());
	int l =0;
	for (int i=0; i< k; i++)	
	{
		l = 0;
		if (i>0 && i < (k-1))
		{
			l = locats[i-1];
		}
		else{
			if (i==0) l = (iterations-1);			//time positions not time values
			if (i==k-1) l = 0;
		}
		 cerr << "L " << l << endl;
		UniformRNG choose;
        
                double index = choose.generate()*double(MCMCiterations);
                int INDEX = int(index);
		        cerr << "INDEX " << INDEX << "\t index " << index << endl;
                cerr << "Temp post " <<Temp_posterior[l][INDEX] << endl;
                Temp[i] =  Temp_posterior[l][INDEX] ;
    	        for (int j=0; j<100; j++) 
    	        {
    	    	    if (Temp[i]>=(double(j)/10 -5) && Temp[i] <(double(j+1)/10 -5))
    	            {
    	    		
    	              Temp_prob =Temp_prob* Temp_hist[l][j];
    	            }
    	        }
		    
	}	
}
//===========================================================================
int sample_times(vector<double> times, vector<double> times_posterior, vector<double> times_old)
{
            UniformRNG choose_times;
            double index = double(choose_times.generate()*times_posterior.size());
    	 	
    	 	int INDEX  = int(index);
    	 	for (int j=0; j<int(times_old.size()); j++)
    	 	{
    	 	    if (abs(times_posterior[INDEX],times_old[j]) < 12.0)
     	 	    {
    	 	         //INDEX = sample_times(times, times_posterior,  times_old);
    	 	         index = double(choose_times.generate()*choose_times.generate()*times_posterior.size());
    	 	         INDEX  = int(index);
    	 	         cerr <<"in this loop " << endl;
    	 	         cerr << "INDEX " << INDEX << endl;
    	 	         cerr <<"times_posterior[INDEX] "<< times_posterior[INDEX] << "\t" << times_old[j] << endl;
    	 	     }
    	 	}
    	 	return INDEX;
}
//==============================================================================
double abs(double a,double b)
{
	double ans = (a-b)*(a-b);
	ans = sqrt(ans);
	return ans;
}

