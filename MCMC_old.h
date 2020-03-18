
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
		int MCMC_SS::Update_SS(vector<GST>& UKnew,  vector<GST>& UKtnew,
		        vector<double> Cx_new, vector<double> Cy_new, vector<int> GST_id,
		        double noise, int nBH, int nC, Array1D<double> HeatFlow, Array1D<double> Tequil,
		         int number, int itmax, vector<double>& beta,  Array1D<int>& BoreholeLength,Array1D<int>& BoreholeOffset,Array1D<double>& UpperLayer,
        Array1D<double>& El_off, Array1D<int> CondLength, Array1D<double> Years);
	 	 
		 int MCMC_SS::Proposal(vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew, 
		 vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta);
		 
		
};

int MCMC_SS::Proposal(vector<int> GST_id, int nBH,int nC, vector<GST>& UKnew,
                vector<GST>& UKtnew, bool & prop, int number, vector<double>& beta)
{
	
	
	beta[number]=0.0;
	UKnew.resize(1);
	UKtnew.resize(1);
	//---------------------
	Label con2Code;
	
	con2Code.set_label(GST_id,GST_id.size());
	con2Code.get_label();
	con2Code.get_size();
	GST_id = con2Code.get_code(GST_id.size());
	
	//---------------------
	
    cerr <<"here  in  Proposal not MCMC";
    GST_Proposal ( GST_id,  nBH, nC, UKnew, UKtnew,  prop,  number,  beta);
    
  
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
		         int number, int itmax, vector<double>& beta,  Array1D<int>& BoreholeLength, 
		         Array1D<int>& BoreholeOffset, Array1D<double>& UpperLayer, Array1D<double>& El_off,
		          Array1D<int> CondLength, Array1D<double> Years)
{
	  
	  UKnew.resize(1);
	  UKtnew.resize(1);
	  vector<double> W;
	  vector<double> told;
      
      //cerr <<"number = " << number << endl;
      
      beta[number] = 0.0;
     
      int points =6;
      //cerr << "points " << points << endl;
      
      double tmax =  iterations * (dt/(3600*365.25*24));
      points = 6;
 
	  vector<double> times;
	  times.assign(points, 0.0);
	  times.resize(points, 0.0);
	  W.resize(points, 0.0);
	  W.assign(points, 0.0);
	
	//---------------------
	Label con2Code;
	
	
	con2Code.set_label(GST_id,GST_id.size());
	con2Code.get_label();
	con2Code.get_size();
	GST_id = con2Code.get_code(GST_id.size());
	for (int i=0; i<int(GST_id.size()); i++)
	{
		cerr << GST_id[i] << "\t";
	}
	cerr << endl;
	
	
	//---------------------   
        char filename[nBH+11];
      
        strcpy(filename, "Samples/");
        for (int i=8; i<nBH+8; i++)
        {
       	    filename[i] = (char) (GST_id[i-8]+48);
        }
        filename[nBH+8] =(char) 49;
        filename[nBH+9] = (char) GST_id[number]+48;
        filename[nBH+10] = (char) 0;
        cerr << filename << endl;
        ofstream outPrT(filename);
        filename[nBH+8] =(char) 48;
        cerr << filename << endl;
	    ofstream outPrt(filename);
	    filename[nBH+8] =(char) 83;
	    cerr << filename << endl;
        ofstream outPrS(filename);
		
	    filename[nBH+8] = (char) 76;
        ofstream outk(filename);
        
        times[0] = tmax;   
        vector<double>::iterator  q ;
        double t =0.0;
        for (q = times.begin(); q!=times.end(); q++)
        {
          *q = tmax - t* tmax/(points-1);
          t +=1;
             cerr <<"time " << *q <<"\t";
        }
        cerr<<endl;
        told = times;
        told[0] = 756;
        told[1] = 700;
        told[2] = 600;
        told[3] = 500;
        told[4] = 400; 
        told[5] = 0;  
      
       //double Func;
      //------------
      RJ_SS Synth1;
      //------------
      
      vector<double> Wnew; 
      Wnew.assign(points, 0.0); 
      Wnew = W;
      int newpoints = points;
      
      vector<double> tnew;
      tnew.assign(points, 0.0);
      tnew = told;
      
      bool Synth = false;
     // double Acceptance;
     
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
      double alpha = 0.0;
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
      bool BIRTH = false, DEATH = false, PERTURB = false, PERTURB_T = false, PERTURB_H= false;
      int MCMCiterations= 50;
      Array2D<double> Tarray = Array2D<double> (iterations,MCMCiterations, 0.0);
      double llold = 1.0;
      double priorold = 1.0;    
      double prior = 1.0; 
      double old_func = 0.0;
      double old_acceptance = 0.0;
      
      Array1D<double>  TEQnew =Array1D<double> (nBH, 0.0);
      Array1D<double> HtFnew =Array1D<double> (nBH, 0.0);
      double alpha_H = 8.0e-004;
      double alpha_Teq = 0.35;
      
          old_func = 0.0;
          prior =0.0;
            for (int st=0; st<nBH; st++)
            {
            	//cerr << "GST id " << GST_id[st]  << endl;
            	if ((GST_id[st]-1) == number)
            	{
            		old_func +=   FE_SS(W,times, Tequil[st],HeatFlow[st], noise, Synth, st,points,  BoreholeLength[st], 
            		BoreholeOffset[st], UpperLayer[st], El_off[st], CondLength[st], Years[st] ); // Functional 
            		
            	    prior += Prior_SS(Wnew,tnew, newpoints);
            	    
            	}
            }
           
             ll=old_func ;
             old_acceptance = Likelihood(llold, priorold, prior,  noise, Wnew, newpoints, tnew, told, 0, BIRTH, DEATH,PERTURB, 
                              ll, wminus, wplus, u1, u2, alpha, JR,0.0);
  
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
      
      int count =1;
      for (int MCMCit=1; MCMCit<itmax; MCMCit++)
      { 
      	
             if ((MCMCit % 100) ==0)
             {
             	cerr << count ;
             	count++;
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
             
           switch (x)
           {
               case 1:
                  //cerr <<"BIRTH" << endl;
                  BIRTH = true;
                  newpoints = points +1;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Birth(W, Wnew, tnew, told, tmax,newT, newpoints,v, wminus, wplus, u1,u2, alpha , JR);
                  break;
                  
              case 2:
                  //cerr <<"DEATH" << endl;
                  DEATH = true;
                  newpoints = points -1;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Death(W, Wnew, tnew, told, tmax,newT, points, newpoints,v, wminus, wplus,u1,u2,  alpha );
                  break;
               case 3:
                  PERTURB = true;
                  x1 = Synth1.Element(points, x, v, JR, told);
                  newpoints = points;
                  tnew.resize(newpoints);
                  Wnew.resize(newpoints);
                  Synth1.Peturb(W, Wnew, tnew,told,tmax,newT, points, newpoints ,v);
                  break;
               case 4:
                   PERTURB_T = true;
                   newpoints = points;
                   x1 = Synth1.Element(newpoints, x, v, JR,told);
                   tnew.resize(newpoints);
                   Wnew.resize(newpoints);
                   Synth1.Peturb_T(W, Wnew, tnew,told,tmax,newT, points, newpoints ,v, half);
                   break;
               case 5:
                   PERTURB_H = true;
                   newpoints = points;
                   UniformRNG HF;
                   double Hv = HF.generate()*(nBH);
                   int hu = int(Hv);
                   int hv = hu + 1;
                                    
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
            		
				   newFunc += FE_SS(Wnew,tnew, TEQnew[st],HtFnew[st], noise, Synth, st,newpoints, BoreholeLength[st], 
				                             BoreholeOffset[st], UpperLayer[st],El_off[st], CondLength[st], Years[st]);     // Functional 
            	   prior += Prior_SS(Wnew,tnew, newpoints);
            	}
            }
            ll = newFunc;
             newAcceptance = Likelihood(llold, priorold, prior, noise, Wnew, newpoints, tnew, 
                                told, v, BIRTH, DEATH,PERTURB, 
                              ll, wminus, wplus, u1, u2, alpha, JR,Yd);
              
           
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
                 
             }
        
             outpoints << points << endl;
             for (int ii=0; ii<iterations; ii++)
             {
               outT <<  Tsurface[ii] << "\t" ;
             }
             outT << endl;
             
             outH << HeatFlow << "\t" << Tequil << endl;  
             
             
             // Output post-burn in Samples to file named by GSTid.
             if (MCMCit >1000 && MCMCit %5 ==0)
             {
             	 
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
                 
                 outk << points << endl;
                
             }
             outll << ll << endl;
             
      	 
      } // end i<itmax
      // -----------------------------------------------------------------------------------------------------------	 
	 
      
     outT.close();
     outpoints.close();
     outH.close();
     outPrT.close();
     outPrt.close();
     outPrS.close();
     outk.close();
	 
     bool prop ;
     //calculate proposal by sampling from the sub-RJ-MCMC output
     GST_Proposal ( GST_id,  nBH, nC, UKnew, UKtnew,  prop,  number,  beta);
 
   
     return points;
}


