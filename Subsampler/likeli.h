// Evaluate likelihood

double Likelihood(double llold, double prold,double& Prior, double noise, vector<double> Wmcmc, int Points, 
          vector<double> tmcmc,vector<double> told, int j, bool Birth, bool Death, bool Perturb, double& ll
          ,double& wminus, double& wplus, double u1, double u2, double beta, double JRatio, double Yd,
 	   double pp, double ppnew, Array1D<double> noise_array,Array1D<double> noise_array_new,int nBH)
{
       
       double PriorRatio = 1.00;
       double J= 1.00; 
       //double L = tmcmc[0] -tmcmc[Points-1];
      
       //cerr << "Points " << Points << endl;      
       if (Birth == true) 
       {
       	  
       	   PriorRatio = (((tmcmc[j-1]-tmcmc[j])*(tmcmc[j] - tmcmc[j+1]))/(told[j-1] -told[j]));
       	   // PropRatio = (double(Points+1)/L);
       	   J =   wminus * beta;
       } 
       if (Death == true)
       {
       	   
       	   PriorRatio = (((told[j-1]-told[j])*(told[j] - told[j+1]))/(tmcmc[j-1] -tmcmc[j]));
       
       	    PriorRatio = 1.000/PriorRatio;
       	    //  cerr <<"beta " << beta << endl;
       	   J = 1.00/(wminus*beta);
       }
       if (Perturb ==true)
       { 
       	   PriorRatio = ((tmcmc[j-1] - tmcmc[j])*(tmcmc[j] - tmcmc[j+1]))/((told[j-1] - told[j])*(told[j] - told[j+1]));
       	   J =1;
       }
       
       
       double alpha = 0.0;
       alpha =(ll-llold)+ (Prior)- (prold)+ log(PriorRatio) + log(JRatio) + log(J)  ;
       

double hyperprior_ratio=0.0;
double hyperlikeli_ratio=0.0;
// Log normal with mu =0.0, sigma =0.25
	double mu=0.15;
	double stdhp=0.15;
   	hyperprior_ratio = -1.000*(pow(log(ppnew)- mu,2))/(2*pow(stdhp,2)) -( -1.000*(pow(log(pp)- mu,2))/(2.000*pow(stdhp,2)));
// Log normal with mu =0.0, sigma =0.13
	mu=-1.3;	
	stdhp=0.5;
   	for (int st=0; st<nBH; st++){
         hyperlikeli_ratio += -1.000*(pow(log(noise_array_new[st])- mu,2))/(2*pow(stdhp,2)) -( -1.000*(pow(log(noise_array[st])- mu,2))/(2.000*pow(stdhp,2)) );
        }
      //  cerr <<"hyper prior ratio "<< hyperprior_ratio << endl;

      //  cerr <<"hyper likeli ratio "<< hyperlikeli_ratio << endl;
alpha+=hyperprior_ratio;

alpha+=hyperlikeli_ratio;
       return   alpha;
}

