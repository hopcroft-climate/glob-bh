
//double Determinant(Array2D<double> a,int n);
//void CoFactor(Array2D<double> a,int n,Array2D<double>& b);
//void Transpose(Array2D<double> a,int n);
//double determinant(Array2D<double> A, int a);
//Array2D<double> invert(Array2D<double> A,  int a);
//double Evidence( double beta, int N, double Ll, double Prior);
// Added sampling routines 9.8.07

class MCMC_SS
{
	public:
		int Update_SS(vector<GST>& UKnew,  vector<GST>& UKtnew,
		        vector<double> Cx_new, vector<double> Cy_new, vector<int> GST_id,
		        double noise, int nBH, int nC, Array1D<double> HeatFlow, Array1D<double> Tequil,
		         int number, int itmax, vector<double>& beta, 
		                     vector<Borehole_info> EuropeSet, Array2D<char> Labels);
	 	 
		 int Proposal(vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew, 
		 vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta,  Array2D<char> Labels);
		 
		
};

int MCMC_SS::Proposal(vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew,
                vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta,  Array2D<char> Labels)
{
	
	
	beta[number]=0.0;
	UKnew.resize(1);
	UKtnew.resize(1);
	//---------------------
	Label con2Code;
	
	//con2Code.set_label(GST_id,GST_id.size());
	//con2Code.get_label();
	//con2Code.get_size();
	//GST_id = con2Code.get_code(GST_id.size());
	
	//---------------------
	
    cerr <<"here  in  Proposal not MCMC";
    GST_Proposal ( GST_id,  nBH, nC, UKnew, UKtnew,  prop,  number,  beta, Labels);
    
  
	return 0;
}

//==============================================================================================
//==============================================================================================
//==============================================================================================
//==============================================================================================
//==============================================================================================
//==============================================================================================
//==============================================================================================
//==============================================================================================

int MCMC_SS::Update_SS(vector<GST>& UKnew,  vector<GST>& UKtnew,
		        vector<double> Cx_new, vector<double> Cy_new, vector<int> GST_id,
		        double noise, int nBH, int nC, Array1D<double> HeatFlow, Array1D<double> Tequil,
		         int number, int itmax, vector<double>& beta,  
		         vector<Borehole_info> EuropeSet, Array2D<char> Labels)
{
	  //Array2D<char> Labels = Array2D<char> (nBH,2,'a');
	  
	  vector<double> W;
	  vector<double> told;
      
    // hyperparameter std dev of normal distribution propr on gsth
      double  pp=1.0000;
      double  ppnew=pp;
          
    
      beta[number] = 0.0;
      cerr << "number " << number << endl;
      cerr << "burn-in " << burn << endl;
      int points =2;
      double tmax =  iterations * (dt/(3600*365.25*24));
     
	  vector<double> times;
	  times.assign(points, 0.0);
	  times.resize(points, 0.0);
	  W.resize(points, 0.0);
	  W.assign(points, 0.0);
	
		
		//Label con2Code;
		
		//con2Code.set_label(GST_id,GST_id.size());
		//con2Code.get_label();
		//con2Code.get_size();
		//GST_id = con2Code.get_code(GST_id.size());
			/*for (int i=0; i<int(GST_id.size()); i++)
			{
			cerr << GST_id[i] << " ";
			}
			cerr << endl;*/
		//cerr << number <<  Cx_new[number] << Cy_new[number] << endl;
		
	    //--------------------------------------  
        
                string FileName;
       string firstpart ="Samples/";
       string finalpart =".txt";
       ostringstream tmpStr;
       
	
      
        int count_id = 0;
        int last_id = 0;
        for (int st=0; st<nBH; st++)
            {
            	if ((GST_id[st]-1) == number)
            	{	
            		
                        cerr << "st " << st << endl;
            		count_id ++;
			last_id = st;
            	}
            }
	    // -----------------------------------------
    
     Array2D<double> bh_locat =Array2D<double> (nBH,2);
     ifstream inpos;
     inpos.open("DATA/BH_locats.txt");
     for(int i=0; i<nBH; i++)
     {
        inpos >> bh_locat[i][0] >> bh_locat[i][1];
         if( bh_locat[i][0] < 0.0)
        { 
          bh_locat[i][0]+=360.0;
        }
       // cerr << bh_locat[i][0]<< "\t" << bh_locat[i][1] << endl;
     }
     inpos.close();
    // Set up the grid of points (converted to a 1x1 square on output)
    // -----------------------------------------
     int mg=73;   //change to 5 deg   longitude
     int jgg=37;   //change to 5 deg    latitude
    
      double dlon = 5.0;
      double dlat = 5.0;
    // -----------------------------------------
    
      Array1D<double> longitude =Array1D<double> (mg,0.0);
     Array1D<double> latitude =Array1D<double> (jgg,0.0);
     for (int i=0; i<mg; i++)
     {
       longitude[i]=(i)*dlon;
       //cerr << longitude[i] << endl;
     }
     for (int j=0; j<jgg; j++)
     {
       latitude[j]=90.0-(j)*dlat;
      // cerr << latitude[j] << endl;
     }
       
      for (int i=0; i<mg-1; i++)
       {
        for (int j=0; j<jgg-1; j++)
        {
          if(bh_locat[last_id][0]>=longitude[i] && bh_locat[last_id][0]<longitude[i+1])
          { if(bh_locat[last_id][1]<latitude[j] && bh_locat[last_id][1]>=latitude[j+1])
           { 
              cerr << "lon-lat "<< longitude[i] <<" " << latitude[j] << endl;
            }
          }
         }
        }
       
	    
     firstpart ="Samples/";
                        finalpart =".txt";
                        
                        tmpStr << number; 
                        string boreholel= tmpStr.str();
                        FileName.append(firstpart);
                        FileName.append(boreholel);
                       //FileName.append(finalpart);
       string FileNameGST=FileName;
       string FileNamet=FileName;
       string FileNameS=FileName;
       string FileNamePP=FileName;
       string FileNameNA=FileName;
       string FileNameHTF=FileName;
       
       //cerr <<FileNameGST.c_str() << endl;
       //beta[number] =1;
       //return points;
       
       FileNameGST.append("G.txt");
       ofstream outPrT(FileNameGST.c_str());
        FileNamet.append("T.txt");
       ofstream outPrt(FileNamet.c_str());
       FileNameS.append("S.txt");
       ofstream outPrS(FileNameS.c_str());
       FileNamePP.append("PP.txt");
       ofstream outPrPP(FileNamePP.c_str());
       FileNameNA.append("NA.txt");
       ofstream outPrNA(FileNameNA.c_str());
       FileNameHTF.append("H.txt");
       ofstream outHp(FileNameHTF.c_str());
       
        string FileName_llfinal;
	FileName_llfinal.append(firstpart);
 	FileName_llfinal.append(boreholel);
	FileName_llfinal.append(finalpart);
	ofstream outll_final;
        outll_final.open(FileName_llfinal.c_str());
        double tot_length=0.0;

        times[0] = tmax;   
        vector<double>::iterator  q ;
        
        
        times[0] = tmax;
        times[1] =  0.0;
        told.assign(points, 0.0);
        told = times;
        
      RJ_SS Synth1;
      
      vector<double> Wnew; 
      Wnew.assign(points, 0.0); 
      Wnew = W;
      int newpoints = points;
      
      vector<double> tnew;
      tnew.assign(points, 0.0);
      tnew = told;
      
      bool Synth = false;
      
     
      vector<double> Wbest; 
      Wbest.assign(points, 0.0);
      
      vector<double> tbest; 
      tbest.assign(points, 0.0);
      
      vector<double>::iterator ff;
      Array1D<double> Tsurface = Array1D<double> (iterations, 0.0);
    
      Array1D<double> Reduced = Array1D<double> ( BH2l, 0.0);
      Array1D<double> Data_R = Array1D<double> ( BH2l, 0.0);
      
      double u1=0;
      double u2 = 0;
      double JR = 1.0;
      double alpha = 5e-5 ; controls the size of the Birth and Death temperature steps.
      double ll = 0;
      double wminus =0.0, wplus = 0.0;int births =0;
      int deaths = 0;
      int perturbs =0;
      int Tchanges =0;
      int birthsa =0;
      int deathsa = 0;
      int perturbsa =0;
      int Tchangea = 0;
      int Hchanges = 0;
      int Hchangea = 0;
      int PPchanges=0;
      int PPchangea=0;
      int NAchanges=0;
      int NAchangea=0;
    
      bool BIRTH = false, DEATH = false, PERTURB = false, PERTURB_T = false, PERTURB_H= false;
      bool PERTURB_PP = false;
      bool PERTURB_NA = false;
      int MCMCiterations= itmax;
      Array2D<double> Tarray = Array2D<double> (iterations,MCMCiterations, 0.0);
      double llold = 1.0;
      double priorold = 1.0;    
      double prior = 1.0; 
      double old_func = 0.0;
      double old_acceptance = 0.0;
      
      Array1D<double>  TEQnew =Array1D<double> (nBH, 0.0);
      Array1D<double> HtFnew =Array1D<double> (nBH, 0.0);


      Array1D<double>  noise_array =Array1D<double> (nBH, 0.05);
      Array1D<double> noise_array_new =Array1D<double> (nBH, 0.05);

      Array2D<double> misfit_array =Array2D<double> (nBH,elements, 0.0);
      Array2D<double> misfit_array_new =Array2D<double> (nBH,elements, 0.0);

      double alpha_H = 1.0e-003;
      double alpha_Teq = 0.35;
      
      old_func = 0.0;
      prior =0.0;

 	

      for (int st=0; st<nBH; st++)
      {
        	//cerr << "GST id " << GST_id[st]  << endl;
		  //cerr <<"st " << st << endl;
        	
        	if ((GST_id[st]-1) == number)
        	{
        		
                        old_func +=   FE_SS(W,times, Tequil[st],HeatFlow[st], noise, Synth, 
        		                              st,points, EuropeSet[st],0,1,misfit_array_new,noise_array[st]); // Functional 
        		prior += Prior_SS(Wnew,tnew, newpoints,pp);
        	}
			
       }
       ll=old_func ;
       old_acceptance = Likelihood(llold, priorold, prior,  noise, Wnew, newpoints, tnew, told, 0, BIRTH, DEATH,PERTURB, 
                          ll, wminus, wplus, u1, u2, alpha, JR,0.0, pp,pp, noise_array,noise_array,nBH);
      
      priorold = prior;
      llold = ll;
      double newT;
      int x =0;
      int x1 = 0;
      int v =0;
      vector<double>::iterator fiter;
      vector<double>::iterator  rr;
      
      UniformRNG b;
      int Accept=0;
      int counter = 0;
      
      ofstream outpoints("Subsampler/Dimensions.txt");
      ofstream outll("Subsampler/Ll.txt");
      ofstream outT("Subsampler/SurfTemp.txt");
      ofstream outH("Subsampler/HeatFlow.txt");

      times = told;
      bool half = false;

       
      
    
      for (int MCMCit=1; MCMCit<=itmax; MCMCit++)
      { 
      	
             if ((MCMCit % 1000)==0)
             {
             	cerr << 100.00*double(MCMCit)/double(itmax)<<"%: " ;
                cerr << 100.0*double(birthsa)/double(births) << "\t" << 100.0*double(deathsa)/double(deaths)  << "\t" << 
			100.0*double(perturbsa)/double	(perturbs)  << "\t" << 100.0*double(Tchangea)/double(Tchanges)  
			<< "\t" << 100.0*double(Hchangea)/double(Hchanges) 
                    << "\t" << 100.0*double(PPchangea)/double(PPchanges)  
                    << "\t" << 100.0*double(NAchangea)/double(NAchanges) 
                    << endl;
             	
             }
      	     BIRTH = false; 
      	     DEATH =false;
      	     PERTURB= false;
      	     PERTURB_T= false;
             PERTURB_H = false;
             newpoints = points;
             Wnew.resize(points);
             Wnew =W;
          	 tnew.resize(points);
          	 tnew = told;
             UniformRNG ccd;
             double Yd = ccd.generate();
             
             x = 0;
             x = Synth1.Prob(points, x, JR, told);
             HtFnew = HeatFlow;
             TEQnew = Tequil;
             
	     //cerr << "x " << x << endl;
	     
           switch (x)
           {
               case 1:{
                  //cerr <<"BIRTH" << endl;
                  BIRTH = true;
                  newpoints = points +1;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Birth(W, Wnew, tnew, told, tmax,newT, newpoints,v, wminus, wplus, u1,u2, alpha , JR);
                  break;
               }
                  
              case 2:{
                  //cerr <<"DEATH" << endl;
                  DEATH = true;
                  newpoints = points -1;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Death(W, Wnew, tnew, told, tmax,newT, points, newpoints,v, wminus, wplus,u1,u2,  alpha );
                  break;
              }
               case 3:{
                  PERTURB = true;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  newpoints = points;
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Peturb(W, Wnew, tnew,told,tmax,newT, points, newpoints ,v);
                  break;
               }
               case 4:
                   {
                   PERTURB_T = true;
                   newpoints = points;
                   x1 = Synth1.Element(points, x, v, JR,told);
                   tnew.resize(newpoints);
                   Wnew.resize(newpoints);
                   Synth1.Peturb_T(W, Wnew, tnew,told,tmax,newT, points, newpoints ,v, half);
                   break;
                   }
               case 5:
                   {
                   PERTURB_H = true;
                   newpoints = points;
                   vector<int> elements_partition;
                   elements_partition.resize(0);
                   for (int st=0; st<nBH; st++)
     		       {
        	         if ((GST_id[st]-1) == number)
        	          {                 
			             elements_partition.push_back(st);
                      }
		           }
		           UniformRNG HF;
                   double Hv = HF.generate()*(int(elements_partition.size()));
                   int hu = int(Hv);
                   //cerr << "hu " << hu << endl;
                   hu=elements_partition[hu];
		           int hv = hu + 1; 
                   //cerr<< "hu+++++++++++++++++++++++++++++++++++++++++ " << hu << endl;
		           HtFnew  = HeatFlow.copy();
                   TEQnew = Tequil.copy();
                   //x1 = Synth1.Element(newpoints, x, v, JR,told);
                   tnew.resize(newpoints);
                   Wnew.resize(newpoints);
                   tnew = told;
                   Wnew = W;
                   Synth1.Peturb_H(HtFnew[hu],HeatFlow[hu],TEQnew[hu],Tequil[hu],alpha_H,alpha_Teq, hu,hv);
                   break;  
                   }
               case 6:
                   {
                   // Sample hyper prior parameter:
                   PERTURB_PP = true;
                   newpoints = points;
                    tnew.resize(newpoints);
                   Wnew.resize(newpoints);
                   tnew = told;
                   Wnew = W;
       HtFnew  = HeatFlow.copy();
                   TEQnew = Tequil.copy();
                   GaussianRNG php;
                   u2 = php.generate()  * 0.1;
                   ppnew = pp + u2;
                       if(ppnew <= 0){
                           ppnew=pp;
                       }
                   break;
                   }
               case 7:
                   {
                   PERTURB_NA = true;
                   newpoints = points;
       HtFnew  = HeatFlow.copy();
                   TEQnew = Tequil.copy();
                   vector<int> elements_partition;
                   elements_partition.resize(0);
                   for (int st=0; st<nBH; st++)
     		       {
        	         if ((GST_id[st]-1) == number)
        	          {                 
			             elements_partition.push_back(st);
                      }
		           }
		   UniformRNG HF;
                   double Hv = HF.generate()*(int(elements_partition.size()));
                   int hu = int(Hv);
                   //cerr << "hu " << hu << endl;
                   hu=elements_partition[hu];
		           int hv = hu + 1; 
                   tnew.resize(newpoints);
                   Wnew.resize(newpoints);
                   tnew = told;
                   Wnew = W;
		   GaussianRNG php;
                   u2 = php.generate()  * 0.01;
                   noise_array_new=noise_array.copy();
                   noise_array_new[hu]=noise_array[hu] +  u2;
                   if(noise_array_new[hu] <= 0.000){
                           noise_array_new[hu]=noise_array[hu];
                       }
		   break;  
                   }
           }
             
           for (int ij=0; ij<newpoints; ij++)
           {
             	if (Wnew[ij] < -10.0 )
             	{
             	        cerr << "newpoints " << newpoints << endl;
             			cerr << "Wnew[ij] " << Wnew[ij] << endl;
             			cerr << "x" << x << endl;
             	        return -1 ;
             	}
            }
            
            double newFunc= 0.0;;
            double newAcceptance;
            
            newFunc = 0.0;
            prior = 0.0;
            for (int st=0; st<nBH; st++)
            {
            	
            	
            	//cerr << "GST id " << GST_id[st] -1 << "\t ss " << ss << endl;
            	if ((GST_id[st]-1) == number)
            	{
            		
				   newFunc += FE_SS(Wnew,tnew, TEQnew[st],HtFnew[st], noise, Synth, st,newpoints, EuropeSet[st],MCMCit,burn,misfit_array,noise_array_new[st]);     // Functional 
            	   prior += Prior_SS(Wnew,tnew, newpoints,ppnew);

            	}
            }
            ll = newFunc;
            // ll = -10.0;
	     newAcceptance = Likelihood(llold, priorold, prior, noise, Wnew, newpoints, tnew, 
                                told, v, BIRTH, DEATH,PERTURB, 
                              ll, wminus, wplus, u1, u2, alpha, JR,Yd,pp,ppnew,noise_array,noise_array_new,nBH);
              
             
             if (newFunc == -999.999)
             {
              cerr << "Tsurface error " << -999.99 << endl;
                
             	return -1;
             }
             double p = 0.0;
             p = newAcceptance;
             double U= b.generate();
             counter ++;
              if (BIRTH == true)
              	 {  births ++;}
              if (DEATH == true)
              	 {  deaths ++;}
              if (PERTURB ==true)
              	  {  perturbs ++;}
              if (PERTURB_T ==true)
              	  {  Tchanges ++;}
              if (PERTURB_H ==true)
              	  {  Hchanges ++;}  
             if (PERTURB_PP ==true)
              	  {  PPchanges ++;}  
             
             if (PERTURB_NA ==true)
              	  {  NAchanges ++;}  
             //Accept?
             U = log(U);
              if (U<p)
              {
              	  if (BIRTH == true)
              	  {  birthsa ++;}
              	  if (DEATH == true)
              	  {  deathsa ++;}
              	  if (PERTURB ==true)
              	  {  perturbsa ++;}
              	  if (PERTURB_T ==true)
              	  {  Tchangea ++;}
              	  if (PERTURB_H ==true)
              	  {  Hchangea ++;}
                  if (PERTURB_PP ==true)
              	  {  PPchangea ++;}  
		  if (PERTURB_NA ==true)
              	  {  NAchangea ++;}
                  llold = ll;
              	  priorold = prior;
          	      points = newpoints;
                  W.resize(points);
          	      W = Wnew;
          	      told.resize(points);
                  told = tnew;
          	      HeatFlow = HtFnew.copy();
                  Tequil = TEQnew.copy();
                  llold = ll; 
                  pp=ppnew;
		  noise_array=noise_array_new.copy();
                  //Acceptance = newAcceptance;
                  Accept++;
                 
                 
              }
              else 
              { 
              	  ll= llold;
              	  prior = priorold;
              	  newpoints = points;
              	  Wnew.resize(points);
              	  Wnew = W;
              	  tnew.resize(points);
              	  tnew = told;
                  HtFnew = HeatFlow.copy();
                  TEQnew = Tequil.copy();
              	  Tsurface  = Synth1.Interpolate_SS(W,  told,newpoints );
                  ppnew=pp;
                  noise_array_new=noise_array.copy();
             }
        
             
             
	         
             
              //cerr << pp << endl;
             // Output post-burn in Samples to file named by GSTid.
             if (MCMCit >=(itmax-burn)  )
	     //if (MCMCit >=(0)  )
             {
             	 
             	 if (MCMCit % 1000 ==0 )
                {
		          outll << ll << endl;	        
		          outpoints << points << endl;
		          for (int ii=0; ii<iterations; ii++)
                  {
                      outT <<  Tsurface[ii] << "\t" ;
                  }
                     outT << endl;
                    for (int st=0; st<nBH; st++)
                    {
            	   if ((GST_id[st]-1) == number)
                   {
                  outHp << HeatFlow[st] << "\t" << Tequil[st] << "\t";
                   }
	               }
		       outHp << endl;
                 }
                 
                 
                 outPrS  << points <<endl;
             	 for (int i=0; i<points ; i++)
             	 {
             	 	outPrT << W[i] <<"\t";
             	 }
                 outPrT << endl;
                 
                 
             	 for (int i=0; i<points ; i++)
             	 {
             	 	outPrt << told[i] <<"\t";
             	 }
                 outPrt << endl;
                 
                outPrPP  << pp <<endl;
                for (int st=0; st<nBH ; st++)
             	 {
                   if ((GST_id[st]-1) == number)
                   {		   
                  outPrNA << noise_array[st] << "\t" ;
                 }
		}
 		 outPrNA << endl;		
                
             }
             
            if(MCMCit%1000==0.0){
              tot_length=0.0;
              for (int st=0; st<nBH; st++)
              {
                if ((GST_id[st]-1) == number)
            	{
                   tot_length+=EuropeSet[st].length_T;
		          }
	          }	
              outll_final << ll << "\t" <<  tot_length << endl;
 	         } 
      	 
      } // end i<itmax
      // -----------------------------------------------------------------------------------------------------------	 
	 
      
     outT.close();
     outpoints.close();
     outH.close();
     outPrT.close();
     outPrt.close();
     outPrS.close();
     outPrPP.close();
     outPrNA.close();
     outll_final.close();
	 
    
    //-----------------------------------------------------------
      //-----------------------------------------------------------
      // Now calculate the posterior mean GST history, the 95% credible intervals
	  //  and the posterior distribution (histomatrix) of the GST history 
      //int MCMCiterations =itmax; 
	/*  int tsteps =iterations;  
	
      int PosteriorSize = double(MCMCiterations-burn)/100;
      int PosteriorBurn = burn/100;
      int smpl = PosteriorSize-PosteriorBurn;
     
      Array2D<double> Temp = Array2D<double> ( PosteriorSize,tsteps, 0.0);
     
    cerr << "MCMCiterations " << MCMCiterations << endl;  
    cerr << "PosteriorSize " << PosteriorSize << endl;
      
	
	ifstream in("Subsampler/SurfTemp.txt");
	for (int i=0; i< PosteriorSize; i++)
	{
		for (int j =0; j<tsteps; j++)
		{
			in >> Temp[i][j]; 
			//cerr  << Temp[i][j] << "\t";
		} 
		//cerr << endl;
	}
	
	//cerr << "End Tsurf read";
	double Tmin = -5.0;
	
	double Trange = 10.0;
	int Tlevels =100;
	Array1D<double> Tsum = Array1D<double> (tsteps, 0.0);
	Array1D<double> l = Array1D<double> (Tlevels+1, 0.0);
	Array2D<int> counterP = Array2D<int> ( tsteps, Tlevels, 0);
	
	firstpart ="Posteriors/Post_";
                        finalpart =".txt";
                        
                        //tmpStr << number; 
                        // boreholel= tmpStr.str();
    string FileNameTm;                    
    FileNameTm.append(firstpart);
                        FileNameTm.append(boreholel);
                       //FileName.append(finalpart);
        FileNameGST=FileNameTm;
       
       
       FileNameGST.append("_HistoMatrix.txt");
       ofstream outHist(FileNameGST.c_str());
    
    //ofstream outHist("Posteriors/HistoMatrix.txt");
	for (int ii=1; ii<=Tlevels; ii++)
	{
		l[ii] = Tmin + ii* Trange/(Tlevels);
		//cerr << l[ii] << endl;
		
	}
	vector<double> min;
	min.assign(tsteps,-5.0);
	vector<double> max;
	max.assign(tsteps,5.0);
	
	Array2D<double> credible_int = Array2D<double> (tsteps, 2, 0.0);
	
	for (int i=PosteriorBurn; i<PosteriorSize; i++)
	{ 
		for (int j=0; j<tsteps; j++)
		{  
			for (int k=0; k<Tlevels+1; k++)
			{
				if ((Temp[i][j] >= l[k]) && (Temp[i][j] < l[k+1]))
				{
					counterP[j][k] ++;
					//cerr << Temp[i][j] << "\t";
					Tsum[j] += Temp[i][j];
					
					if (Temp[i][j] < min[j])
						min[j] = Temp[i][j];
					
					if (Temp[i][j] > max[j])
					{ max[j] = Temp[i][j];
					}
				}	
			}
		}
	} 
	
	/*
	//Output to file the posterior mean of the GST history:	
       firstpart ="Posteriors/Post_";
                        finalpart =".txt";
                        
                        //tmpStr << number; 
                         //boreholel= tmpStr.str();
    string FileNameH;                    
    FileNameH.append(firstpart);
                        FileNameH.append(boreholel);
                       //FileName.append(finalpart);
        FileNameGST=FileNameH;
       
       
       FileNameGST.append("_Tmean.txt");
       ofstream outTmean(FileNameGST.c_str());
       //ofstream outTmean("Posteriors/Tmean.txt");
        cerr << "------------------------------" << endl;
	cerr << "Writing Tmean.txt ..." << endl;
        for (int j=0; j<tsteps; j++)
	{
		Tsum[j] = double(Tsum[j])/ double(smpl);
		outTmean << Tsum[j]<< endl;
	}
	outTmean.close();
	
	// -----------------  OUT histomatrix-----------------------------------
	cerr << "Writing HistoMatrix.txt ..." << endl;
	for (int k=Tlevels-1; k>=0; k--)
    {
	    for (int j=0; j<tsteps; j++)
	    {
			outHist << counterP[j][k] << "\t";
		}
		outHist << endl;
	}
	outHist.close();
	//---------------------------- Credible intervals --------------------
    int smpl5 = int ( 0.05*smpl);
    int smpl95 =int( 0.95*smpl);
    vector<double> sorted;
    sorted.assign(smpl, 0.0);
    Array2D<double> Tlim = Array2D<double> (2,tsteps, 0.0);
    for (int i=0; i<tsteps; i++)
    {
    	sorted.assign(smpl,0.0);
    	 int k = 0;
    	for (int j=PosteriorBurn; j<PosteriorSize; j++)
    	{
    		
    		sorted[k] = Temp[j][i];	
    		k++;
    	}
    	sort(sorted.begin(),sorted.end());
    	Tlim[0][i] = sorted[smpl5];
    	Tlim[1][i] = sorted[smpl95];
    	
    }
    cerr << "Writing Crediblelimits.txt ..." << endl;
   firstpart ="Posteriors/Post_";
                        finalpart =".txt";
                        
                        //tmpStr << number; 
                        // boreholel= tmpStr.str();
    string FileNameC;                    
    FileNameC.append(firstpart);
                        FileNameC.append(boreholel);
                       //FileName.append(finalpart);
    cerr <<"FileNameC " << FileNameC << endl;   
    FileNameGST=FileNameC;
        
       
       FileNameGST.append("_Crediblelimits.txt");
    cerr <<"Finename" << FileNameGST << endl;   
    ofstream out95(FileNameGST.c_str());
    //ofstream out95("Posteriors/Crediblelimits.txt");
    for (int i=0; i<tsteps; i++)
    {
    	out95 << Tlim[0][i] <<"\t" << Tlim[1][i] << endl;
    }
    out95.close();
    
    cerr << "RJ-MCMC programme finished." << endl;
    cerr << "------------------------------" << endl;
    */
    
     bool prop ;
     //calculate proposal by sampling from the sub-RJ-MCMC output
     GST_Proposal ( GST_id,  nBH, nC, UKnew, UKtnew,  prop,  number,  beta, Labels);
 
   
     return points;
}


