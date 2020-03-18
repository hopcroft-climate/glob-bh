//Bayesian Partition Modelling
//14.03.07
//3-12.8.07  added probabilities of sampling a particular thermal history model

// 27.2.2018: Completely reconfigured from above to fixed spatial partitions.
//             Now includes sampling the measurement noise assumed for each borehole site
//              and in each partition, the assumed width of the prior on the temperature history.
//         With a noise of 0.05K and prior of 1.0K, the NH mean GST is ~0.6 K warming since 1750
//         Doubling the prior std dev. to 2.0K increases the GST NH warming to 1.1 K
//         As does using a uniform prior over a wider interval.
//         Using adaptive prior, will hopefully more accurately account for uncertainty and misift in data.
// 13.3.2018: Bug fix: get rid of all calls to con2code, as this mixes up the ordering of the spatial partitions, by re-writing 
//             over GST_ID!        

  
//#include "includesMCMC.h"
#include <iostream> 
#include <vector>
#include <fstream>  
#include <cmath>     
#include <vector>       
#include <algorithm>              
#include <string>              

#include "TNT/tnt_array2d.h"
#include "TNT/tnt_array2d_utils.h"
#include "TNT/tnt_array1d.h"
#include "TNT/tnt_array1d_utils.h"
//#include "C:\workspace\TNT\jama_cholesky.h" 
       
  
#include <time.h>                      // define time()
  using namespace TNT; 
       
using namespace std;
#include "BoreholeInfo.h"

#include "baseRNG.h"
#include "Data.h" 
#include "parameters.h"
 
  
#include "GSTclass.h"
#include "Voronoiclass.h"
#include "LabelClass.h"

#include "HeatFlow.h"
#include "RJ.h"
#include "local_matrices.h"
#include "elementclass.h" 
#include "trd_solver_bc.h" 
#include "trd_solver_NR.h" 
#include "Tlinear.h"  
#include "FE.h" 
// RJMCMC sub sampler for GST history values.
#include "Subsampler/SS_RJ.h"
#include "Subsampler/func.h"
#include "Subsampler/likeli.h"
#include "GST_Proposal_Mark_III.h"     
#include "MCMC_Mark_IV.h"

 
int main () 
{  
	bool Synth = false;    //Synth=true when running forward model to produce synthetic data.
	
	int nBH =1012;
    
    
	
	int itmax = 1;  //main RJMCMC sampler
	//---------------------------------------------------------
	Array2D<char> Labels = Array2D<char> (nBH, 2,'a');
	Array1D<double> HtF = Array1D<double> (nBH, 0.060);
        Array1D<double> HtFnew = Array1D<double> (nBH, 0.060);
        Array1D<double> Teq = Array1D<double> (nBH, 9.0);
        Array1D<double> Teqnew = Array1D<double> (nBH, 9.0);
	
	 
	 // vector of Borehole info structures for borehole data:
	 vector<Borehole_info> EuropeSet;
	 EuropeSet.resize(nBH);
	 //---------------------------------------------------------------------
	 
	nBH = getdata(EuropeSet, Labels,Teq,HtF,nBH);
	if(nBH==1){
          cerr<< "problem in data " << endl;
	}	
	// Borehole locations
	
	
	//test---------------------------------------------------------------------
	/*vector<double> tester;
	tester.assign(EuropeSet[20].length_T, 1.0);
	ifstream intest;
	
	
	for (int i=0; i<nBH; i++)
	{
		intest.open(EuropeSet[i].name_T.c_str());
		for (int ja=0; ja<EuropeSet[i].length_T; ja++)
		{
		//	cerr <<  ja << "\t";
			intest >> tester[ja];
		//	cerr  << tester[ja] << endl;
		}
		intest.close();
		//cerr << nBH << endl;
	}
	//--------------------------------------------------------------------------
	*/
		
	Array2D<double> BH_LOC = Array2D<double> (nBH, 2);
	
	//data files:
	Array2D<double> data = Array2D<double> (nBH,elements, 0.0);
	
	// Initial number of tessellation points
	int nC = 173;
	int nCnew = nC;
	
	vector<double> Cx; 
        Cx.assign(nC, 0.5);
        vector<double> Cy; 
        Cy.assign(nC, 0.5);
        vector<double> Cx_new; 
        Cx_new.assign(nC, 0.5);
        vector<double> Cy_new; 
        Cy_new.assign(nC, 0.5);
	// Input labels for grouping in Partitions:
	ifstream inBHl("DATA/Labels.txt");
	for (int l=0; l<nBH; l++)
	{
		inBHl >> Labels[l][0] >> Labels[l][1];
		cerr << Labels[l][0] << Labels[l][1]  << endl;
	}
        inBHl.close();
	ifstream ingrid("BH_get_grid_squares/bhlocats.txt");
        for (int i=0; i<nC; i++)
	{
	    ingrid>> Cx[i] >> Cy[i] ;
            cerr << Cx[i] << "\t" << Cy[i] << endl;
	}
        cerr << endl;
	ingrid.close();
	
	//---------------------------------------------------------
	//GST histories for each tessellation point
	vector<int> GST_id; 
	GST_id.assign(nBH,0);
	vector<int> GST_id_new; 
	GST_id_new.assign(nBH,0);
	vector<bool> SubSample;
	SubSample.assign(nC,false);
	
	vector<int> GST_id_number;
	GST_id_number.assign(nC, 0);
	
	
  
    Teqnew = Teq;
    HtFnew  = HtF;
    
    bool BIRTH = false; 
    bool DEATH =false;
    bool PERTURB= false;
    bool MOVE= false;
    bool PERTURB_GST = false;
  
    bool LLeval = true;
  
    ofstream outacc("Out/Acceptance.txt");
    
    int Accept =0;
    int birthsa=0; 
    int deathsa =0; 
    int perturbsa =0; 
    int movesa=0; 
    int GSTmovesa=0; 

    int births = 0;
    int deaths =0; 
    int perturbs =0; 
    int moves=0; 
    int GSTmoves=0; 
    double d= 0.0;
  
    //--------------------------------------------------------
   /* ofstream outGST;
    outGST.open("Out/GST.txt");
    ofstream outGSTt;
    outGSTt.open("Out/GSTt.txt");
    ofstream outHtF;
    outHtF.open("Out/HtF.txt");
    ofstream outnC;
    outnC.open("Out/nVCentres.txt");
    
    ofstream outLL;
    outLL.open("Out/LLprior.txt");
    ofstream outCx;
    outCx.open("Out/CentresX.xls");
    ofstream outCy;
    outCy.open("Out/CentresY.xls");
    ofstream outID;
    outID.open("Out/GST_ID.txt");
    ofstream outHP;
    outHP.open("Out/HyperPrior.txt");
    ofstream outBeta;
    outBeta.open("Out/Beta.txt");*/
    //---------------------------------------------------------
    vector<double> Wmcmc;
    vector<double> tmcmc;
    
    double LLP = 0.0;
    double LLP_old = 0.0;
	double Prior = 0.0;
	double Prold = 0.0;
	double JRatio = 1.0;
	double JR =1.0;
	double length = 0.0; 
	//---------------------------------------------------------
	vector<double> beta;
	beta.resize(nC, 1.0);
	vector<double> beta_new;
	beta_new.resize(nCnew, 1.0);
	
	//double L = (iterations * dt)/(24.0*3600.0*365.25);
	
	double mu_new = 1;
	vector<double> init;
	init.assign(2,0.0);

	// Define BH locations over a square grid 1.0 X 1.0---
	ifstream inbh_pos;
        inbh_pos.open("BH_get_grid_squares/bh_pos.txt");
	for (int b=0; b<nBH; b++)
	{
		//BH_LOC[b][0] = EuropeSet[b].bx;
		//BH_LOC[b][1] = EuropeSet[b].by;
                inbh_pos >> BH_LOC[b][0] >> BH_LOC[b][1] ;
		EuropeSet[b].bx=BH_LOC[b][0];
		EuropeSet[b].by=BH_LOC[b][1];
                cerr << "Position: " << EuropeSet[b].bx <<"\t" <<   EuropeSet[b].by<< endl;
	}
	 
	Cx_new = Cx;
	Cy_new = Cy;
	
    //Vector of Objects for variable dimensionality
	
	vector<GST> UK(nC);
	vector<GST> UKt(nC);
	vector<GST> UKnew(nC);
	vector<GST> UKtnew(nC);
    
    int number = 0;
    vector<GST> UKtest;
    vector<GST>UKttest; 
    
    bool prop;
    
	vector<double> flat;
    flat.assign(6,0.0);
	
	vector<double> init_t;
	init_t.resize(2,0.0);
	init_t[0] = 756.0;
	init_t[1] = 0.0;
	
	//==================================================================
	// Test initial model
	Voronoi test;
	for (int i=0; i<nC; i++)
	{
	     d = test.get_dist(Cx_new[i],Cy_new[i],BH_LOC[3][0],BH_LOC[3][1],nC);
	}
	
	// Which region is borehole i in?
	Voronoi ident;
	ifstream ingstid;
        ingstid.open("BH_get_grid_squares/gstid.txt");
        for (int i=0; i<nBH; i++)
	{
		//GST_id[i]  = ident.which_region(Cx,Cy,BH_LOC[i][0],BH_LOC[i][1],nC);
		ingstid >> GST_id[i];
                cerr << GST_id[i] << " ";
	}
ingstid.close();
	cerr << endl;
	Label con2Code;
	
		
		//con2Code.set_label(GST_id,GST_id.size());
		//con2Code.get_label();
		//con2Code.get_size();
		//GST_id = con2Code.get_code(GST_id.size());
			for (int i=0; i<int(GST_id.size()); i++)
			{
			cerr << GST_id[i] << "\n";
			}
			cerr << endl;	

		GST_id_new= GST_id;
	
	for (int ic=0; ic<nC; ic++)
    {	
         cerr << Cx[ic] << Cy[ic] << endl;
    }
   
	
	int minsofar;
	double minsofar2;	
	double maxzsofar;
	int hole ;
		
		for (int ic=0; ic<nC; ic++)
	{
	
	  cerr << "partition i: " << ic << endl;				
		minsofar=0;
		minsofar2=1000.000;
                maxzsofar=0.000;
   	    for (int st=0; st<nBH; st++)
            {
            	if ((GST_id[st]-1) == ic)
            	{	
            		
                        cerr << "st " << st << endl;
			cerr << "GST_id" << GST_id[st] << endl;
 			
			cerr << EuropeSet[st].years << endl;
            		cerr << "Europe offset " << EuropeSet[st].offset << endl;
            		if(EuropeSet[st].years>minsofar){
				hole=st;
                                minsofar=EuropeSet[st].years;
				}
			
			if(EuropeSet[st].DataZ[0]<minsofar2){
				hole=st;
                                minsofar2=EuropeSet[st].DataZ[0];
				}
                               if(EuropeSet[st].DataZ[EuropeSet[st].length_T-1]>maxzsofar){
				hole=st;
                                maxzsofar=EuropeSet[st].DataZ[EuropeSet[st].length_T-1];
				}		
			
                      
            	}
            }
		
		cerr << "i:nC , minsofar " <<  "\t" << minsofar << endl;		                
                cerr << "i:nC , z(0) " <<  "\t" << minsofar2 << endl;		
                cerr << "max z " << "\t" << maxzsofar << endl;
		}
		//return 0;
		
	// ------------------------------------------
	// Main partition loop for sub-sampler RJ-MCMC
	// ------------------------------------------
	
	for (int ic=0; ic<113; ic++)
    {	
    	  UKtest.resize(1);
	  UKttest.resize(1);
         // MCMC_SS Subchain1;
          
          prop = false;
          
         // Subchain1.Proposal(GST_id,nBH,nC,UKtest,UKttest, prop, ic, beta, Labels);
          
          if (prop ==false)
          {
  	          UK.resize(nCnew);
              UKt.resize(nCnew);
  	          MCMC_SS Subchain;
  	         // cerr << "here nBH " << nBH  << endl;
              int k =6;
              //cerr<< "EuropeSet size " << EuropeSet.size() << endl;
              k = Subchain.Update_SS(UKtest,UKttest, Cx, Cy, GST_id,noise,nBH,nC,  HtFnew,  Teqnew, ic, SSit, 
              								beta, EuropeSet, Labels);
              //cerr << " K " << k << endl;
              
              if (beta[ic] == 0.0)
              {
              	cerr <<" beta = 0 !" << endl;
              	return 0;
              }
              if (k<1)
              {
                  cerr <<"ERROR in subchain " << endl;
              	  return 0;
              }
          } 
            vector<double> Tvec;
         	vector<double> tvec;
         	
         	
         	Tvec.assign(UKtest[0].get_size(),0.0);
         	tvec.assign(UKttest[0].get_size(),0.0);
         	Tvec = UKtest[0].get_history();
         	tvec = UKttest[0].get_history();
         	UK[ic].set_history(Tvec,UKtest[0].get_size());
         	UKt[ic].set_history(tvec,UKttest[0].get_size());
         	
        }
        
    	UKnew = UK;
    	UKtnew = UKt;	
	    
	    //cerr << "Beta " ;
	    
	    beta_new.assign(beta.size(), 0.0);
	    for (int i=0; i<nC; i++)
	    {
	    	cerr << beta[i] <<"\t";
	    	beta_new[i] = beta[i];
	    }
	    cerr << endl;
	    
        LLP = Likelihood_prior(UK,UKt, GST_id, HtF,Teq, noise, Synth, nC,nBH,Prior,EuropeSet);
        cerr << "initial LLP " << LLP << endl;
        
        // stop if running forward model only
        if (Synth == true)
        {
        	return 0;
        }
	    
	    
	    return 0;
	
	//===================================================================
	//===================================================================
	
	for (int MCMCit =0; MCMCit<itmax; MCMCit++) // start MCMC iterations
	{
             BIRTH = false; 
      	     DEATH =false;
      	     PERTURB= false;
      	     MOVE= false;
             PERTURB_GST = false;
			 
    	 //make old:
    	 LLP_old = LLP;
    	 Prold = Prior;
    	 			
    	 Cx_new = Cx;
    	 Cy_new = Cy;
    	 
    	 beta_new  = beta;
    	 
    	 // cerr <<" size Cx, Cxnew " << Cx.size() << "\t" << Cx_new.size() << endl;
    	 // cerr <<" size Cy, Cynew " << Cy.size() << "\t" << Cy_new.size() << endl;
    		
    	 GST_id= GST_id_new;
    	 SubSample.resize(nC,false);
    	 
    	 for (int i=0; i<nC; i++)
             {	
           	    UKnew[i].set_history(UK[i].get_history(), UK[i].get_size());
           	    UKtnew[i].set_history(UKt[i].get_history(),UKt[i].get_size());
             }
    	 for (int ih=0; ih<nBH; ih++)
    	 {
    	 	HtFnew[ih] = HtF[ih];
    	 	Teqnew[ih] = Teq[ih];
    	 }
    	 
    	 Voronoi Move;
    	 UniformRNG u;
    	 RJ xis ;
    	 
    	 int x;
    	 xis.choose(nC, x, JR, nBH);
    			 
    			  
     //=============================================================================
     //=============================================================================
     cerr << "X " << x <<endl;
     UniformRNG H;
     double HU  ;
     int hu  ;
     double V ;
	 int v ;
	 
           switch (x)
           {
              
              case 1:
              //Voronoi Birth
             
              BIRTH =true;
              nCnew = nC +1;
              cerr <<"I BIRTH ..........nC, nCnew " << nC << "\t" << nCnew << endl;
              
              beta_new.resize(nCnew);
              
			  //cerr << "nC, nCnew " << nC << "\t" << nCnew << endl;
					
			    Cx_new.resize(nCnew);
		        Cy_new.resize(nCnew);
				//cerr <<"now call move.birth "<< endl;
				
				Move.birth(Cx_new,Cy_new,Cx,Cy,nCnew, length);
				
				UKnew.resize(nCnew);
				UKtnew.resize(nCnew);
				for (int i=0; i<(nCnew-1); i++)
				{
					UKnew[i] = UK[i];
					UKtnew[i] = UKt[i];
					beta_new[i] = beta[i];
				}
				
				UKnew[nCnew-1].set_history(init,2);
				UKtnew[nCnew-1].set_history(init_t,2);
				
        break;
              
              //=============================================================================
              case 2:
              //Voronoi Death
              DEATH =true;
              nCnew = nC -1;
              cerr <<"II DEATH ...............nC, nCnew " << nC << "\t" << nCnew << endl;
              
			
				V = u.generate()*(nC);
			    v = int(V);
				//cerr << " v " << v << endl;
				Move.death(Cx_new,Cy_new,Cx,Cy,nCnew,v,length);
				
				beta_new.resize(nCnew);
				UKnew.resize(nCnew);
				UKtnew.resize(nCnew);
				for (int i=0; i<v; i++)
				{
					UKnew[i] = UK[i];
					UKtnew[i] = UKt[i];
					beta_new[i] = beta[i];
				}
				for (int i=v+1; i<(nCnew+1); i++)
				{
					UKnew[i-1] = UK[i];
					UKtnew[i-1] = UKt[i];
					beta_new[i-1] = beta[i];			//index for beta_new is i-1.  11.9.07
				}
				break;
              
              //=============================================================================
              case 3:
              //Voronoi Move
              // cerr <<"III MOVE ..............."<<endl;
              MOVE = true;
              nCnew = nC;
              V = u.generate()*(nC);
		        v = int(V);
			if (v >= nC)
			{
				cerr << "v greater than nC " << endl;
				return 0;
			}
			Cx_new.resize(nCnew);
		        Cy_new.resize(nCnew);
		        beta_new.resize(nCnew);
		        beta_new = beta;
			Move.move(Cx_new,Cy_new,Cx,Cy,nC,v);
			UKnew.resize(nC);
		              UKtnew.resize(nC);
			UKnew = UK;
		        UKtnew = UKt;
              break;
              
              //=============================================================================
              case 4:
              //Voronoi Perturb
             // cerr <<"IV PETURB ..............."<<endl;
              PERTURB = true;
              nCnew = nC;
              Cx_new.resize(nCnew);
	          Cy_new.resize(nCnew);
	           beta_new.resize(nCnew);
	           beta_new = beta;
	           Move.perturb(Cx_new,Cy_new,Cx,Cy,nC,length);
	           UKnew.resize(nC);
	           UKtnew.resize(nC);
	           UKnew = UK;
	           UKtnew = UKt;
	        break;
              
              //=============================================================================
            case 5:
             
             PERTURB_GST =true;
             //update the heatflow and teq for one borehole:
             HtFnew = HtF.copy();
             Teqnew = Teq.copy();
              HU = H.generate()*nBH;
              hu = int(HU);
             prop_hf(HtFnew,  HtF,   Teqnew,Teq, alpha_H,  alpha_Teq,  nBH,  hu );
             
            
              
             break;
              //=============================================================================
              case 6:
              //Update GST 
              //cerr <<"V GST UPDATE..............."<<endl;
              
              PERTURB_GST =true;
              nCnew = nC;
              Cx_new.resize(nCnew);
              Cy_new.resize(nCnew);
              UKnew.resize(nC);
              UKtnew.resize(nC);
              UKnew = UK;
              UKtnew = UKt;
              
              HtFnew = HtF.copy();
              Teqnew = Teq.copy();
              beta_new.resize(nCnew);
              beta_new = beta;
              
              UniformRNG wchR;
              double IC = wchR.generate()*nC;
              int ic = int(IC);
              ic= GST_id[ic]-1;
              
              for (int i=0; i<nBH; i++)
	          {
	              GST_id[i]  = ident.which_region(Cx_new,Cy_new,BH_LOC[i][0],BH_LOC[i][1],nCnew);
	    				
	          }
			//---------------------------------------------------------  
			MCMC_SS Subchain1;
			
			double Resample = wchR.generate();
			
			if (Resample < 0.5)
			{
				Subchain1.Proposal(GST_id,nBH,nCnew,UKtest,UKttest, prop, ic, beta_new, Labels);
			}
			else{
				prop = false;
			
			}
	        
       		//cerr << "after proposal prop is "<< prop << endl;
       		if (prop ==false)
             {
             	
                  UKnew = UK;
                  UKtnew = UKt;
                 	
      		      //UKtest.resize(1);
      		      //UKttest.resize(1);
      		      // cerr <<"go to subchain" << endl;
      		    
                  int k = Subchain1.Update_SS(UKtest,UKttest,Cx_new, Cy_new, GST_id,noise,nBH,  nC, HtFnew,  
                                                                         Teqnew, ic,SSit, beta_new, EuropeSet,Labels);
                  //cerr << " k " << k << endl;
                  if (k<1)
                  {
                      cerr <<"ERROR in subchain " << endl;
                  	return 0;
                  }
              }
            
              vector<double> Tvec;
              vector<double> tvec;
              Tvec.assign(UKtest[0].get_size(),0.0);
              tvec.assign(UKttest[0].get_size(),0.0);
              Tvec = UKtest[0].get_history();
              tvec = UKttest[0].get_history();
              UKnew[ic].set_history(Tvec,UKtest[0].get_size());
              UKtnew[ic].set_history(tvec,UKttest[0].get_size());
              
              for (int i=0; i<nC; i++)
              {
              	UKnew[i].get_history();
              }
              
              //??  UKnew,UK,UKtnew,UKt
              break;
             
              
              //=============================================================================
              
              
              
            }
              
              //=================================================================================
           
	Label con2Code;
	
	//con2Code.set_label(GST_id,GST_id.size());
	//con2Code.get_label();
	//con2Code.get_size();
	//fGST_id = con2Code.get_code(GST_id.size());
			
            //================================================================================
                
	SubSample.assign(nCnew,false);
	for (int i=0; i<nCnew; i++)
	{
		d = test.get_dist(Cx_new[i],Cy_new[i],BH_LOC[0][0],BH_LOC[0][1],nC);
	}
		
  		Voronoi ident;
  	
  		for (int i=0; i<nBH; i++)
  		{
  			GST_id_new[i]  = ident.which_region(Cx_new,Cy_new,BH_LOC[i][0],BH_LOC[i][1],nCnew);
  			
  		}
        ident.any_data(GST_id_new, GST_id, UKnew,  UKtnew,  UK, UKt, 
                                 nCnew,nC,  nBH, Cx_new, Cy_new, Cx, Cy, beta , beta_new);
  		
  		
  		//---------------------
	Label con2Code1;
	
	
	//con2Code1.set_label(GST_id_new,GST_id_new.size());
	//con2Code1.get_label();
	//con2Code1.get_size();
	//GST_id_new = con2Code1.get_code(GST_id_new.size());
			   		 
  		
  		for (int i1=0; i1<nBH; i1++)
  		{
  			if (GST_id_new[i1] != GST_id[i1])
  			{
  				SubSample[GST_id_new[i1]-1] = true;
  				SubSample[GST_id[i1]-1]=true;
  				   					
  			}
  		}
  		
  		//---------------------------------------------------------
  		LLeval = false;
  		if (PERTURB_GST == true)
		{
			LLeval = true;
		}
  		
  		
  		for (int ic=0; ic<nCnew; ic++)
  		{
    		if (SubSample[ic] ==true)
    		{
    			LLeval = true;
    			 bool prop =true;
    			 UKtest.resize(1);
    			 UKttest.resize(1);
    			 MCMC_SS Subchain1;
  		                 
    			 Subchain1.Proposal(GST_id_new,nBH,nCnew,UKtest,UKttest, prop, ic, beta_new, Labels);
  		                 
    			 if (prop ==false)
                         {
	    		    UKnew.resize(nCnew);
                    UKtnew.resize(nCnew);
	    		    number = ic;
	    		    // cerr <<"New data config: subsample: number=" << number  << endl;
	    		    MCMC_SS Subchain;
	    		    int k = Subchain.Update_SS(UKtest,UKttest, Cx_new, Cy_new, GST_id_new,noise,nBH,  
	                 								nCnew,  HtFnew,  Teqnew, ic, SSit, beta_new,EuropeSet,Labels);
	                //cerr << " K " << k << endl;
	    		    if (k<1)
	    		    {
	                       cerr <<"ERROR in subchain " << endl;
	                 	return 0;
	    		    }
                         } 
               
               	vector<double> Tvec;
               	vector<double> tvec;
               	Tvec.assign(UKtest[0].get_size(),0.0);
               	tvec.assign(UKttest[0].get_size(),0.0);
               	Tvec = UKtest[0].get_history();
               	tvec = UKttest[0].get_history();
               	UKnew[ic].set_history(Tvec,UKtest[0].get_size());
               	UKtnew[ic].set_history(tvec,UKttest[0].get_size());
                
    	    }
  		}
  		
  		//---------------------------------------------------------------------------
  		double accept=0.0;
  		if (LLeval == true)
		{
  			LLP = Likelihood_prior(UKnew,UKtnew, GST_id_new, HtFnew,Teqnew, noise, Synth, nCnew,nBH,Prior,EuropeSet);
  		
        	//cerr << "LL X p =" << LLP << endl;
        
       
        	//acceptance term
        	RJ BPM;
        	accept = BPM.accept(LLP,LLP_old, Prold, JR,BIRTH, DEATH,PERTURB,MOVE, 
                                        JRatio, nCnew,nC, length,mu_new, beta_new, beta);
                                        }
		else{
			accept = log(1.0);
		}
        //---------------------------------------------------------
        UniformRNG acceptanceU;
        double U= acceptanceU.generate();
        //acceptance criterion
        U = log(U);
        
        if (BIRTH == true)
              {  births ++;}
              if (DEATH == true)
              {  deaths ++;}
              if (PERTURB ==true)
              {  perturbs ++;}
              if (MOVE ==true)
              {  moves ++;}
              if (PERTURB_GST ==true)
              {  GSTmoves ++;}
        //---------------------------------------------------------
        if (U < accept)
        {
          Accept ++;
          if (BIRTH == true)
            	  {  birthsa ++;}
            	  if (DEATH == true)
            	  {  deathsa ++;}
            	  if (PERTURB ==true)
            	  {  perturbsa ++;}
            	  if (MOVE ==true)
            	  {  movesa ++;}
            	  if (PERTURB_GST ==true)
            	  {  GSTmovesa ++;}
        	
        	//accept new model
        	//cerr <<"ACCEPT Model**************************************" << endl;
        //	cerr <<"**************************************************" << endl;
        	nC = nCnew;
        	UK.resize(nCnew);
        	UKt.resize(nCnew);
        	for (int i=0; i<nCnew; i++)
        	{
        		UK[i].set_history(UKnew[i].get_history(),UKnew[i].get_size());
        		UKt[i].set_history(UKtnew[i].get_history(), UKtnew[i].get_size());
        	}
        	for (int ih=0; ih<nBH; ih++)
        	{
        		HtF[ih] = HtFnew[ih];
        		Teq[ih] = Teqnew[ih];
        	}
                Prold = Prior;
            
                Cx.resize(nCnew);
                Cy.resize(nCnew);
        	
        	Cx = Cx_new;
        	Cy = Cy_new;
        	Cx_new.resize(nCnew);
        	Cy_new.resize(nCnew);
        	LLP_old = LLP;
        	
        	GST_id = GST_id_new;
            
                beta.resize(nCnew,0.0);
                beta = beta_new;
                
            
        }
        //---------------------------------------------------------
        else
        {
			//reject new model, reinstate current model
        	//cerr <<"REJECT Model:::::::::::::::::::::::::::::" << endl;
        	//cerr <<":::::::::::::::::::::::::::::::::::::::::" << endl;
        	nCnew = nC;
        	UKnew.resize(nC);
        	UKtnew.resize(nC);
        	for (int i=0; i<nC; i++)
        	{
        		UKnew[i].set_history(UK[i].get_history(),UK[i].get_size());
        		UKtnew[i].set_history(UKt[i].get_history(),UKt[i].get_size());
        	}
        	for (int ih=0; ih<nBH; ih++)
        	{
        		HtFnew[ih] = HtF[ih];
        		Teqnew[ih] = Teq[ih];
        	}
        	
        	Prior = Prold;
        	Cx_new.resize(nC);
        	Cy_new.resize(nC);
        	//cerr <<" Reject: size Cx, Cxnew " << Cx.size() << "\t" << Cx_new.size() << endl;
		    //cerr <<" Reject: size Cy, Cynew " << Cy.size() << "\t" << Cy_new.size() << endl;
			
        	Cx_new = Cx;
        	Cy_new = Cy;
        	Cx.resize(nC);
        	Cy.resize(nC);
        	LLP = LLP_old;
        	
        	GST_id_new = GST_id;
        	beta_new.resize(nC, 0.0);
        	beta_new = beta;
        	
        	
        }
        //---------------------------------------------------------
        for (int i=0; i<nC; i++)
        {
        	if (beta[i] == 0.0)
        	{
        		cerr << "beta not defined " << endl;
        		return 0;
        	}
        }
        //---------------------------------------------------------
        // output UK, UKt,C, LL, nC, UK.size()?
        for (int io=0; io<nC; io ++)
  		{
  			Wmcmc = UK[io].get_history();
  			tmcmc = UKt[io].get_history();
  			
  			/*for (int jo=0; jo<(UK[io].get_size()); jo++) 
  			{
  				outGST << Wmcmc[jo] << "\t";
  				outGSTt << tmcmc[jo] <<"\t";
  			}
  			outGST << endl;
  			outGSTt << endl;*/
  		}
  		//---------------------------------------------------------
  		
     /*         for (int ho =0; ho<nBH; ho++)
	{
		outHtF << HtF[ho] <<"\t";
		
	}   
	for (int ho =0; ho<nBH; ho++)
	{
		outHtF << Teq[ho] <<"\t";
		
	} 
	outHtF << endl;
     
	outLL << LLP << endl;
       
	for (int i=0; i<nC ; i++)
	{
	    outBeta << beta[i] <<"\t";
	}
	outBeta << endl;*/
	//---------------------------------------------------------
      /*  if (MCMCit % 10 == 0)
        {
        	 outnC << nC << endl;
	        for (int io=0; io<nC; io++)
	        {
	             outCx << Cx[io] << "\t" ;
	       //      cerr <<"Cx " <<  Cx[io] << "\t" ;
	        }
	       // cerr << endl;
	         for (int io=0; io<nC; io++)
	        {
	             outCy << Cy[io] << "\t" ;
	      //       cerr <<"Cy " <<  Cy[io] << "\t" ;
	        }
	       // cerr << endl;
	        outCx << endl;
	        outCy << endl;
	        
	        for (int i=0; i< nBH; i++)
	        {
	        	outID << GST_id[i] << "\t" ;
	        }
	        outID << endl;
       
          	cerr<<"Acceptance rate  " <<100.0*double(Accept)/double(MCMCit+1)<<"%" << endl;
          	outacc << Accept << "/" << (MCMCit+1) << "\t";
          	outacc <<  birthsa << "/" << births << "\t";
 		    outacc << deathsa<< "/" << deaths << "\t";
  		    outacc << perturbsa<< "/" << perturbs << "\t";
            outacc << movesa<< "/" << moves << "\t";
            outacc << GSTmovesa << "/" << GSTmoves<< endl;
	}*/
	        	
	//---------------------------------------------------------   
	}
      /*  outGST.close();
        outGSTt.close();
        outHtF.close();
        outnC.close();
        outLL.close();
        outCx.close();
        outCy.close();
        outID.close();
        outBeta.close();*/
    //-----------------------------------------------------------
    
    cerr << "RJ-MCMC programme finished." << endl;
    cerr << "------------------------------" << endl;
       return 0;
}



