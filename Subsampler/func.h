
// Calculates forward model or the gradient by the adjoint method.
//#include <math.h>
//#include "quad.h"
// uses regularization with integral of gradient of m squared
     
     
    
double FE_SS(vector<double> Tmcmc, vector<double> times, double Teq, double HtF, double noise, bool Synth, 
					 int BHnumber, int points, Borehole_info EuropeSet,int MCMCit,int burn,
					 Array2D<double> misfit_array,double noiseA)
{   
     
    for (int i=0; i<points; i++)
    {
    	 Tmcmc[i] += Teq;
    	// cerr << "Tsurf " << Tmcmc[i] << endl;
    }
    //cerr << "NOise in likelihood :" << noiseA << endl;   
    int elements_local = EuropeSet.length_T+ 5 + 20;
	//cerr << "elements_local " << elements_local << endl;
    Array1D<double> T = Array1D<double> (elements_local, 0.0);
    RJ_SS newone;
    Array1D<double> Tsurf = Array1D<double> (iterations, 0.0);
    Tsurf = newone.Interpolate_SS(Tmcmc, times, points );
       if(Tsurf[0]==-999.99){  return -999.999;}
       
	T = TFe(Tsurf, T, points, elements_local, HtF,Teq,dx, EuropeSet.upperLayer, EuropeSet.length_C, BHnumber, 
	              EuropeSet.length_T, EuropeSet.years,EuropeSet.DataZ,EuropeSet.cond,EuropeSet.years );
    	
     // cerr <<"EuropeSet.Data " << FileName << EuropeSet.Data ;
     double weight = 0.0;
     Array1D<double> misfit = Array1D<double> (EuropeSet.length_T, 0.0);
     double valfun = 0.0;
    // misfit_array = Array1D<double> (elements, 0.0);
     for (int ijk =(EuropeSet.offset); ijk<(EuropeSet.length_T +EuropeSet.offset); ijk++)
       {
       	 weight  = 0.0;
       	 weight  =  -0.5*(( EuropeSet.Data[ijk-EuropeSet.offset]-T[ijk] )*(1/(noiseA*noiseA))*( EuropeSet.Data[ijk-EuropeSet.offset]-T[ijk]));
	 valfun += weight ;
         misfit[ijk-EuropeSet.offset] = T[ijk] -EuropeSet.Data[ijk-EuropeSet.offset];
         //cerr << EuropeSet.Data[ijk-EuropeSet.offset] << "\t " << T[ijk] << endl;
	 
	 misfit_array[BHnumber][ijk-EuropeSet.offset]=misfit[ijk-EuropeSet.offset];
	 
       }
       
       
       
       //double mul =0.0; 
       //mul =  (pow((2*3.14159265),(EuropeSet.length_T/2)))*(pow((noiseA*noiseA),(EuropeSet.length_T/2)));
       //mul = (1.000/(mul));
       //mul = log(mul);

        long double mul =0.0; 
        long double exponent = double(EuropeSet.length_T/2.000)*log(noiseA*noiseA);
        mul = -double(EuropeSet.length_T/2)*log(noiseA*noiseA)- exponent;

       valfun = mul+ (valfun) ;
       
       /*
       char OutName[10];
      
        
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
      ofstream outm(OutName);
       
       for (int i=0; i<EuropeSet.length_T; i++)
       {
       	 outm << misfit[i] << endl;
       }
       outm.close();  
     */
     
     	
       
       /*ofstream outPr("Out/Tprof.txt");
       	for (int is=0; is<elements; is++)
       	{
       		
       		outPr << T[is] << endl;
       	}
       	outPr.close();*/
      
       if (Synth== true)
       {
         	
	        GaussianRNG NOISE;
	        double nn = 0.0;
       	    ofstream outN("Data/noise.txt");
        	ofstream outd("Data/M.txt");
         	for (int is=0; is<elements; is++)
        	{
       		  nn = 0.0;
       		  nn = NOISE.generate()*0.1;
       		  T[is] += nn;
       		  outd << T[is] << endl;
       		  outN << nn << endl;
        	}
        	outd.close();
        	outN.close();
       }
    
   for (int i=0; i<points; i++)
    {
    	 Tmcmc[i]-= Teq;
    }
    //cerr << "iterations " << iterations << endl;
    //cerr << "Tsurface " << Tsurface << endl;
    for (int i=0; i<points; i++)
    {
    	if (Tmcmc[i] < -30.0)
    	{   
    		cerr << Tmcmc[i] << endl;
    		cerr <<"Tsurf < -10.0 " << endl;
    		valfun = -999.999;
    	}
    }
  
       /* string FileName;
	string firstpart ="Samples/";
                       string  finalpart =".txt";
                    ostringstream tmpStr;     
                        tmpStr << BHnumber; 
                        string boreholel= tmpStr.str();
                        FileName.append(firstpart);
                        FileName.append(boreholel);
                       //FileName.append(finalpart);
          string FileNameHTF=FileName;
          FileNameHTF.append("M.txt");
	  ofstream outT(FileNameHTF.c_str(), ios::out | ios::app);
	 
         if (MCMCit >=(burn)  )
	     //if (MCMCit >=(0)  )
             {
             	 
             	 if (MCMCit % 10 ==0 )
                {
		          
		      for (int ijk =0; ijk<EuropeSet.length_T; ijk++)
                   {
                
                      outT <<  misfit[ijk ] << "\t" ;
                   }
                     outT << endl;
                 }
	      }
               outT.close();
  */
  
        //cerr <<" valfun " << valfun << endl;
        return (valfun) ;
}

//================================================================
double Prior_SS(vector<double> Wmcmc, vector<double> tmcmc, int points, double pp)
{
	double Prior = 0.0;
   
       double P = 0.0;
       Array1D<double> prior = Array1D<double> (points, 0.0);
       
        // This is an optional prior
        // that approximates the Huang et al 2000 average warming
        // of 1.1C since AD 1500:
        //   for (int j=0; j<points; j++)
        //   {
        //       prior[j] = 1.1 * (  tmcmc[j])/600.0 ;
        //       
        //   }
        // end of optional prior
	// This is an optional prior
        // that approximates the Pollack&Smerdon 2004 average warming
        // of 1.05C since AD 1500:
         //  for (int j=0; j<points; j++)
          // {
           //    prior[j] = 5e-6 * pow(tmcmc[j],2)  -0.0014*tmcmc[j] +0.0583;
            //   cerr << tmcmc[j] << " " << prior[j] << "\t";
          // }
	 //  cerr << endl;
         //end of optional prior
	
       
       
       double PstdDev=pp;
       
       for (int j=0; j<points; j++)
       {
	PstdDev=pp*(1.000+ ((2015.00 - tmcmc[j])/2015.00));     
        //cerr << "tmcmc[j], PstdDev "<< tmcmc[j] << "\t" << PstdDev << endl;     
	P += 0.500*(Wmcmc[j] - prior[j])*(Wmcmc[j] - prior[j])*(1.000/(PstdDev*PstdDev));
          //cerr <<"tmcmc[j] " << tmcmc[j] << endl;
         // cerr << "priorj " << prior[j] << endl;
       } 
       
       double mul =0.0; 
       mul =  (pow((2*3.14159265),(points/2)))*(pow((PstdDev*PstdDev),(points/2)));
       mul = (1.000/(mul));
       mul = log(mul);
       Prior = mul+ (-P) ;
       
       // Calculate priors on the times.
      
       //cerr <<"factorial(points)" << factorial(points);
      // cerr << "pow(L,points) " << pow(L,points) << endl;
      // cerr <<" t_prior " << t_prior << endl;
       //-------------------------------------------------------
       
         //-------------------------------------------------------
    // This is an optional uniform prior with fixed +-3 C bounds:
   /* P=0.0;  
    double factorial =1.0;
    for (int j=0; j<points; j++)
    {
          factorial *= (j+1);
	  }
    for (int j=0; j<points; j++)
    {
	  
	  // Equation 3.16 of book Denison et al 2001.
	  P += log (factorial * pow(6.0,double(-points))* (1./(points+1)));
	  
	  // effectively rule out temperatures outside of this range.
	  if(abs(Wmcmc[j]) > 3.0){
	     P+= log(10000000.0);
	     
	     }
	     
    }
    
    Prior = -P;*/
    //-------------------------------------------------------
    
      
   return Prior;
}

