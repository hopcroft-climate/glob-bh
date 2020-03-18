 //void sample_dimensions(Array1D<int> k_sample,int k_sample_length, int& k, double& k_prob);

void times_histomatrix(vector<double> times, int points, Array1D<double>&  time_histomatrix, 
				vector<double> & times_posterior, int & tot, double dt_new, int iterations_new);

  void Temps_histomatrix(vector<double> times, vector<double> Temp, int points, 
            Array2D<double> & Temp_hist,Array2D<double> & Temp_posterior,  
            Array1D<double> Tsurf, int j,int discretisations, double dt_new );

//void sample_Temps(vector<double>& Temp, Array2D<double> Temp_hist, int k,Array2D<double>   Temp_posterior, 
				//double& Temp_prob,vector<int> locats, int MCMCiterations );
 
 
 
 
 void  GST_Proposal (vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew,
                vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta,
                        Array2D<char> Labels)
{
       
       int MCMCiterations = 10000;    
       int discretisations = 100;  // for temperatures
       
       double dt_new =2.0;
       int iterations_new = 301;
      
       vector<double>  Tvec;
       vector<double>  tvec;
       Array2D<double> Tarray = Array2D<double> (iterations_new, discretisations, 0.0);
       
	    string FileName;
       string firstpart ="Samples/";
       string finalpart =".txt";
       
	
        int count_id = 0;
        for (int st=0; st<nBH; st++)
            {
            	if ((GST_id[st]-1) == number)
            	{	
            		
                        firstpart ="Samples/";
                        finalpart =".txt";
                        ostringstream tmpStr;
                        tmpStr << (st); 
                        string boreholel= tmpStr.str();
                        FileName.append(firstpart);
                        FileName.append(boreholel);
                       //FileName.append(finalpart);
            		count_id ++;
            	}
            }
        
       
       //filename[count_id+9] = (char) GST_id[number]+48;
       string FileNameGST=FileName;
       string FileNamet=FileName;
       string FileNameS=FileName;
       string FileNameHTF=FileName;
       
       FileNameGST.append("G.txt");
       ifstream inPrT(FileNameGST.c_str());
        FileNamet.append("T.txt");
       ifstream inPrt(FileNamet.c_str());
       FileNameS.append("S.txt");
       ifstream inPrS(FileNameS.c_str());
       FileNameHTF.append("H.txt");
       ifstream inHp(FileNameHTF.c_str());
        
	   //==================
	  
	  
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
       vector<double> GSTvec;
	   vector<double> GSTtvec ;
       if (k_sample[MCMCiterations-1] > 0 && k_sample[MCMCiterations-1] < 20) 
       {
    	    prop = true;     
	        Array2D<double> Temp_hist = Array2D<double> (iterations_new, discretisations, 0.0);
     	    Array2D<double>  Temp_posterior = Array2D<double> (iterations_new, MCMCiterations, 0.0);
		   
		    
		    // sample dimensions:----------------------------------
		    //-----------------------------------------------------
		    int k=0;
		    double k_prob = 0.0;
		   // sample_dimensions(k_sample,MCMCiterations,  k,  k_prob);
		    
		   //cerr  <<" k is set to "<< k << endl;
		    //sample times:----------------------------------------
		    //-----------------------------------------------------
		    Array1D<double> times_hist = Array1D<double> (iterations_new, 0.0);
		    vector<double> times_posterior;
		    vector<double> times;
		    times.resize(k, 0.0);
		    int tot = 0;
		   
		    Array1D<double> Tsurf = Array1D<double> (iterations_new, 0.0);
		    
			UniformRNG chooseGSThistory;
		    double Sampled = chooseGSThistory.generate()*(MCMCiterations-1);
		    int sampled = int(Sampled);
		    Array1D<double> GST = Array1D<double> (iterations_new, 0.0);
		    
			for (int jj=0; jj<MCMCiterations; jj++)
		    {
		    	Tvec.assign(k_sample[jj], 0.0);
		    	tvec.assign(k_sample[jj], 0.0);
		    	for (int i=0; i<k_sample[jj] ; i++)
			    {
			    		inPrT >> Tvec[i];	
			    		inPrt >> tvec[i];
			    		
			    }
			    if (tvec[0] == 0.0)
			    {
			    	prop = false;
			    	jj = MCMCiterations-1;
			    }
			    
			    if (k_sample[jj] >2)
			    {
			    	times_histomatrix(tvec, k_sample[jj], times_hist,times_posterior, tot,dt_new,iterations_new);
			    }
			    
			    
			    // Interpolate the Tvec first
			    RJ_SS Tpost;
			    Tsurf = Tpost.Interpolate_SS(Tvec, tvec, k_sample[jj] );
			    Temps_histomatrix(tvec, Tvec,   nP[jj],  Temp_hist, Temp_posterior, 
			    			Tsurf, jj, discretisations, dt_new );
			    if (jj == sampled)
			    {
			    	//for (int i=0; i<iterations_new; i++)
			    	//{
			    	//	GST[i]  = Tsurf[i] ;
			    		
			    	//}
			    	GSTvec = Tvec;
			    	GSTtvec = tvec;
			    	k = GSTvec.size();
			    }
		    }
		//	for (int i=0; i<k; i++)
		//	{
		//		cerr << "Chosen GSTory " ;
		//		cerr << GSTtvec[i] << "\t" << GSTvec[i] << endl;
		//	}

		    Array1D<double> Tot = Array1D<double> (iterations_new, 0.0);
			ofstream outTPost("Tpost.txt");
			for (int i=0; i<iterations_new; i++)
			{
				for (int k2=0; k2<discretisations; k2++)   //given the discretisations Temperature discretisations
				{
						Tot[i] += Temp_hist[i][k2];
						
				}

				for (int j=0; j<discretisations; j++)
				{
					outTPost << Temp_hist[i][j] <<"\t";
					if (Temp_hist[i][j] >0)
					{
						Temp_hist[i][j] = double(Temp_hist[i][j])/Tot[i];
					}
					else Temp_hist[i][j] = 0.0;
				
					
				}
				outTPost << endl;
			}
			outTPost.close();
		    double pptimes = 0.0;
		    for (int i=0; i<iterations_new; i++) 
		    {
		        //normalise the time values pdf
				//cerr << "times-hist " << times_hist[i] << endl;
		    	
		    	times_hist[i] = double( times_hist[i])/double(tot);
		         pptimes += times_hist[i];
				 
		    }
		  //  cerr <<"Sum_TIMES " << pptimes << endl;
		    
			//=================================================================================
		    
		    inPrT.close();
		    inPrt.close();
		    // select (k-2) times:
		    times.assign(k, 0.0);
		    times = GSTtvec;
		    double time_prob = 1.0;
		    vector<int> locats;
			
		    if (k>2)
		    {
				for (int i=0; i<(iterations_new) ;i++) 
				{
     				for (int j=1; j<k-1; j++)// don't sample the end points.
     				{
     					if (times[j] >= double(i)*dt_new && times[j] < double(i+1)*dt_new)
     					{
		    	 			time_prob = time_prob* times_hist[i];
		   	  				locats.push_back(i);
		    		 	}
		    		}
				}
		    }
		    //cerr << "size locats " << locats.size() << endl;
		    // need to sort the time-vector in descending values 
		    sort(times.begin(), times.end());
		    reverse(times.begin(), times.end());
		    // need to also sort the locations vector (int) descending values
		    sort(locats.begin(),locats.end());
		    reverse(locats.begin(), locats.end());
		    // set the proposal to the sample time-vector of length k.
		    
		    double Temp_prob = 1.0;
		    for (int i=0; i< k; i++)	
			{
				int l = 0;
				if (i>0 && i < (k-1))
				{
					l = roundnew(times[i]/dt_new);
				}
				else{
					if (i==0) l = iterations_new-1;
					if (i==k-1) l = 0;
				}
		
		    	for (int j=0; j<discretisations; j++) 
    	    	{
    	    		   if (GSTvec[i]>=(double(j)/double(discretisations/10.0) -5) && GSTvec[i] <(double(j+1)/double(discretisations/10.0) -5))
    	     	      {
    	    				Temp_prob =Temp_prob* Temp_hist[l][j];
    	  //  				cerr << "Temp hist "<< Temp_hist[l][j] << endl;
    	     	      }
    	    	}
		    }
		    k = GSTvec.size();
		    k_prob = 1.0;
		    UKtnew[0].set_history(times,k);
		    UKnew[0].set_history(GSTvec,k);
		    
		 //   for (int i=0; i<k; i++)
		  //  {
		  //  	cerr << times[i] << "\t";
		  //  	cerr << Tvec[i] << "\t";
		  //  }
		   // cerr << endl;
		    beta[number] = time_prob * k_prob* Temp_prob;
			
			//cerr << "k "<< k << endl;
		   // cerr <<"p(time) " << time_prob << "\t p(k) " << k_prob << "\t p(T) " << Temp_prob << endl;
		   // cerr << "beta(number) " << beta[number] << endl;
		   // cerr <<" number " << number << endl;
		    //===================
		    
		    if (beta[number] >0  &&  beta[number] <= 1.0)
		    {
		    ;	
		    }
		    
		    else
		    {
		    ofstream outerror("Tpost_error.txt");
		    
		    	for (int i=0; i<iterations_new; i++)
		    	{
		    		for (int j=0; j<discretisations; j++)
					{
					   outerror << Temp_hist[i][j] <<"\t";
					}
					outerror << endl;
		    	}
		    	outerror.close();
		    	beta[number] = -10.0;
		    	 
		    	
		    }
		    
		
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
 	Array1D<double> k_hist = Array1D<double> (20, 0.0);
 		
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
 	int index =  int(choose.generate()*k_sample_length);   
 	k = int(k_sample[index]);
 	//calculated the (normalised) probability of selecting this k value.
 	
 	double ksum = 0.0;
 	for (int i=0; i<20; i++)    
 	{
 		k_hist[i]= double(k_hist[i])/double(k_sample_length);
 		ksum += k_hist[i];
 	}
 	cerr <<"KSUM " << ksum << endl;
    
    k_prob = double(k_hist[k]); // normalise the pdf for dimensions
 	  
 }
 
//=====================================================================
   
   void times_histomatrix(vector<double> times, int points, Array1D<double>&  time_histomatrix, 
                vector<double> & times_posterior, int & tot, double dt_new, int iterations_new)
   {
   	
   	  
     for (int i=0; i<(iterations_new) ;i++) 
     {
     	for (int j=1; j<points-1; j++)// don't sample the end points.
     	{
     		if (times[j] >= double(i)*dt_new && times[j] < double(i+1)*dt_new)
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
  void Temps_histomatrix(vector<double> times, vector<double> Tvec, int points, 
            Array2D<double> & Temp_hist,Array2D<double> & Temp_posterior,
            Array1D<double> Tsurf, int j, int discretisations , double dt_new)
  {
	 
	 for (int i=0; i<int(Tvec.size()); i++)
	 {
		int l = roundnew(times[i]/dt_new);
		//cerr <<" l " <<l << endl;
		for (int k=0; k<discretisations; k++)   //given the discretisations Temperature discretisations
		{
			if (Tvec[i] >=(double(k)/double(double(discretisations)/10.0) -5.0)&& Tvec[i] <(double(k+1)/double(double(discretisations)/10.0) -5.0))
			{
				Temp_hist[l][k] ++;  //calculate the Temp/time pdf
				
				Temp_posterior[i][j] += Tvec[i];
			
			}

         }	
		
	}

	
}
//==================================================================
/*void sample_Temps(vector<double>& Temp, Array2D<double> Temp_hist, int k,Array2D<double> Temp_posterior, 
                  double& Temp_prob,vector<int> locats, int MCMCiterations )
{
	Temp_prob = 1.0;
	Temp.resize(k, 0.0);
	//cerr << "size locats " << locats.size() << endl;
	vector<double> temporary;
	int l =0;
	for (int i=0; i< k; i++)	
	{
		l = 0;
		if (i>0 && i < (k-1))
		{
			l = locats[i-2];
		}
		else{
			if (i==0) l = iterations_new-1;
			if (i==k-1) l = 0;
		}
		
		UniformRNG choose;
        
                double index = choose.generate()*MCMCiterations;
                int INDEX = int(index);
		
                Temp[i] =  Temp_posterior[l][INDEX] ;
    	        for (int j=0; j<discretisations; j++) 
    	        {
    	    	    if (Temp[i]>=(double(j)/5 -5) && Temp[i] <(double(j+1)/5 -5))
    	            {
    	    		
    	              Temp_prob =Temp_prob* Temp_hist[l][j];
    	            }
    	        }
		    
	}	
	//check
	double pp = 0.0;
	for (int i=0; i<discretisations; i++)
	{
		pp += Temp_hist[15][i];
	}
	cerr << "PP " << pp << endl;
}
//==================================================================
		      	*/	
		      	
