
double Calc_prior(vector<double> Wmcmc, vector<double> tmcmc, int points,double L);
//double FE(TNT::Array1D<double> Tsurface, double TEq,double HtF, double noise, bool Synth, int BHnumber);
double Analytical(Array1D<double> Tsurf, double Teq, double HtF, double noise, 
            bool Synth, int BHnumber, int points, Borehole_info EuropeSet);

double Likelihood_prior( vector<GST> UKnew, vector<GST> UKtnew, vector<int> GST_id, 
                  Array1D<double> HtF,Array1D<double> Teq, double noise, bool Synth, int nC,int nBH, 
                  double& Prior, vector<Borehole_info> EuropeSet )


{
	
	
	Array2D<double> Tsurface = Array2D<double> (iterations,nC, 0.0);
	Array1D<double> Tsurf = Array1D<double> (iterations, 0.0);
	int m = 0;
	double LL = 0.0;
	double LLP = 0.0;
	vector<double> Wmcmc;
	vector<double> tmcmc; 
	Prior = 0.0;
	//cerr << "here in FE" << endl;
	
	// Calculate temperatures in each borehole, using the nC GST histories (Tsurface).
	double L = iterations ;
	LL = 0.0;
	double r=0.0;
	for (int i=0; i<nC; i++)
	{
		//cerr << "here in loop " << endl;
		Wmcmc.resize(UKnew[i].get_size());
	  	tmcmc.resize(UKtnew[i].get_size());
		
		Wmcmc = UKnew[i].get_history();
    	tmcmc = UKtnew[i].get_history();
		
		Interpolate(Wmcmc,  tmcmc, UKnew[i].get_size(), Tsurface, i);
	     
	     }
	     
	     
	     
		LL = 0.0;
	
	for (int i=0; i<nBH; i++)
	{
		// must subtract 1 from GST_id, since all vectors etc are 0 reference
		m = GST_id[i];
		//cerr <<"GST history to use " << m << endl;
		for (int j=0; j<iterations; j++)
		{
			Tsurf[j] = Tsurface[j][m-1];
			
		}
		
		 //add probabilities because they are LOG likelihoods.
		
		r= Analytical(Tsurf,Teq[i],HtF[i],  noise,  
		               Synth, i, UKnew[m-1].get_size(), EuropeSet[i]);     
        LL +=r;
	}
	
	LL = (LL);
			
		LL= -1.0;	
			
			
			
	// Calculate the prior probability for each borehole 
	Prior = 0.0;
	double p = 0.0;
	for (int i=0; i<nC; i++)
	{
		
		Wmcmc.resize(UKnew[i].get_size());
		Wmcmc = UKnew[i].get_history();
		
		tmcmc.resize(UKtnew[i].get_size());
		tmcmc = UKtnew[i].get_history();
	
		p= Calc_prior(Wmcmc,tmcmc, Wmcmc.size(), L);
		Prior +=p;
	}
	
	
	// LL * Prior ( add when logarithmic)
	cerr << "LL " << LL << endl;
	LLP=0.0;
	double LgPr = Prior;

	LLP = LL + LgPr;
	
	return LLP;
}
//========================================================================================================


//=========================================================================================================

double Calc_prior(vector<double> Wmcmc, vector<double> tmcmc, int points,double L)
{
   double Prior = 0.0;
   
       double P = 0.0;
       Array1D<double> prior = Array1D<double> (points, 0.0);
       
       for (int j=0; j<points; j++)
       {
          P += 0.5*(Wmcmc[j] - prior[j])*(Wmcmc[j] - prior[j])*(1/(PstdDev*PstdDev));
        //  cerr <<"Wmcmc[j] " << Wmcmc[j] << endl;
       //   cerr << "priorj " << prior[j] << endl;
       } 
       
       double mul =0.0; 
       mul =  (pow((2*3.14159265),(points/2)))*(pow((PstdDev*PstdDev),(points/2)));
       mul = (1/(mul));
       mul = log(mul);
       Prior = mul+ (-P) ;
       
       // Calculate priors on the times.
       double t_prior = 1.0;
       
       for (int j=1; j<points ; j++)
       {
       	  t_prior  = t_prior*(-tmcmc[j]+tmcmc[j-1]);
          //cerr <<tmcmc[j] <<"\t" << tmcmc[j-1] << endl;
       }
       
       t_prior = (t_prior*factorial(points) )/(pow(L,points));
       t_prior = log(t_prior);
       
       //-------------------------------------------------------
      // cerr <<"tprior" << t_prior << endl;
    //   cerr <<"Tprior" << Prior << endl;
       Prior = t_prior + Prior;
       
       
       //Prior on number of points in GST history
	   // Uniform 
	   double Pgst = 1.0/18.0;
	   Pgst = log(Pgst);
       
       Prior += Pgst;
       return Prior;   // returning log of the prior
}

//===============================================================================================================

double Analytical(Array1D<double> Tsurf, double Teq, double HtF, double noise, 
            bool Synth, int BHnumber, int points, Borehole_info EuropeSet)

{  
     for (int i=0; i<iterations ;i++)
     {
     	Tsurf[i] += Teq;
     	//cerr << "TMcmc " << Tmcmc[i] << endl;
     }
    
     
     Array1D<double> T = Array1D<double> (elements, 0.0);
     T = TFe(Tsurf, T, points, elements, HtF,Teq,dx, EuropeSet.upperLayer, EuropeSet.length_C, BHnumber, 
                                     EuropeSet.length_T, EuropeSet.years,EuropeSet.DataZ,EuropeSet.cond,EuropeSet.years);
    
     //==================================
   
   
    
     double weight = 0.0;
     Array1D<double> misfit = Array1D<double> (EuropeSet.length_T, 0.0);
     double valfun = 0.0;
     
     for (int ijk =(EuropeSet.offset); ijk<(EuropeSet.length_T +EuropeSet.offset); ijk++)
       {
       
       	 weight  = 0.0;
       	 weight  =  -0.5*(( EuropeSet.Data[ijk-EuropeSet.offset]-T[ijk] )*(1/(noise*noise))*( EuropeSet.Data[ijk-EuropeSet.offset]-T[ijk]));
       //	 cerr << ( T[ijk] ) << "\t EuropeSet.Data " << EuropeSet.Data[ijk-EuropeSet.offset] << endl;
         valfun += weight ;
         misfit[ijk-EuropeSet.offset] = EuropeSet.Data[ijk-EuropeSet.offset] - T[ijk];
         
         
       }
       //cerr << endl;
       //cerr <<" valfun FE " << valfun << endl;
      // double mul =0.0;
       //mul =  (pow((2*3.14159265),(BH2l/2)))*(pow((noise*noise),(BH2l/2)));
      // mul = 1/mul;
       //mul = 7.30281e-91;
      // mul = log(mul);
       // valfun += mul;
       
       
      char OutName[10];
      
      ofstream outm(OutName);
      //if (BHnumber < 23)
      //{
        
		OutName[0] = (char) (79);		//o
		OutName[1] = (char) (85);		//u
		OutName[2] = (char) (84);		//t
		OutName[3] = (char) (47);		// 'forward slash'
		OutName[4] = (char) (77);		//M
		OutName[5] = (char) (70);		//f
		OutName[6] = (char) (47);		// 'forward slash'
		OutName[7] = (char) (65+BHnumber);		// label
		OutName[8] = (char) (65+BHnumber);
		OutName[9] = (char) 0;
       
       for (int i=0; i<EuropeSet.length_T; i++)
       {
       	 outm << misfit[i] << endl;
       }
       outm.close();  
     //}
     //else
     //{
     //	;
     // }
      
               
      
      /*ofstream outPr("Out/Tprof.txt");
       	for (int is=0; is<elements; is++)
       	{
       		
       		outPr << T[is] << endl;
       	}
       	outPr.close();*/
      
         
    for (int i=0; i<iterations; i++)
    {
    	Tsurf[i] -= Teq;
    } 
   
    for (int i=0; i<iterations; i++)
    {
    	if (Tsurf[i] <-10.0)
    	{
    		cerr <<"Tsurf < -10.0 " << endl;
    		valfun = -999.999;
    	}
    }
  
        return (valfun) ;
}
