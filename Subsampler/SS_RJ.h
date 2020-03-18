// generic class for a RJMCMC update object
// 4 move types, Birth/Death/Move time point/Perturb temp. point
int roundnew(double a);
class RJ_SS
    { 
          public:
                 int SetPoints(int points);      
                
                 int Birth(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double>& told,
				double tmax,double newT, int newpoints , int& v,double& wminus, double& wplus, 
				double& u1, double& u2, double& alpha, double& JR  );
				
     		 vector<double>Death(vector<double> W,vector<double>& Wnew,vector<double>& tnew,vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int v ,double& wminus, double& wplus,
				double& u1, double& u2 , double& alpha );
				
		 int Peturb(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& v  );
		    
		 int Peturb_T(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& element, bool& half  );
		    
		 int Peturb_H(double& HtFnew, double& HtF, double& TEQnew, double& Teq, double alph_H, double alph_Teq, int a, int b );
		    
		 int Prob(int points, int& x, double& JR,vector<double> tmcmc);
		    
		 Array1D<double>Interpolate_SS(vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints );
           
                 int Element(int points, int x, int& v, double& JR, vector<double> tmcmc);
            
          
          private:
			  int k; // no of points in Tsurf
			 
 
                   
    };

//--------------------------------------------------------------------------------------------------------
    int RJ_SS::SetPoints(int points)
    {      k =points; 
    	return 0;  
    }

    
//--------------------------------------------------------------------------------------------------------
	int RJ_SS::Birth(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double>& told,
				double tmax,double newT, int newpoints , int& v, double& wminus, double& wplus, 
				double& u1, double& u2, double& alpha, double& JR )
	{
       
          
           
        
 //          cerr << "told " << endl;
		   
//		   for (int i=0; i<newpoints-1; i++)
//		   { cerr << told[i] << "\t" ;}
//		   cerr << endl;
		   
//		    cerr << "W " << endl;
//		   for (int i=0; i<newpoints-1; i++)
//		   { cerr << W[i] << "\t" ; }
//		   cerr << endl;
       
     
     
     if (v >0)  
     {for (int i=0;i<v; i++)
         {
         	tnew[i] = told[i] ;
          }
     }
     else
     {
      	 tnew[0] = told[0];
     }
     for (int i=v+1; i<newpoints; i++)
     {
     	tnew[i] = told[i-1];
     	
     }
//   tnewError:
     UniformRNG ccU;
     
     u1 = ccU.generate();
     tnew[v] = told[v]+((dt/(24.0*3600.0*365.25)))+ u1* (told[v-1]- told[v] -(2*(dt/(24.0*3600.0*365.25))));      
   
       
     wminus = told[v-1]- told[v];   
     for (int i=0; i<v; i++)
     {
     	Wnew[i] = W[i];
     }
     for (int i=v+1; i<newpoints; i++)
     {
     	Wnew[i] = W[i-1];
     }
     GaussianRNG ccG;
     u2 = ccG.generate();
  
     int iter= -1;
     Array1D<double>Tsurface = Array1D<double> (iterations+1, 0.0);
     Array1D<int> pos = Array1D<int> (2, 0);
    
     pos[0] = int(told[v-1]/((dt/(24.0*3600.0*365.25))));
     pos[1] = int(told[v]/((dt/(24.0*3600.0*365.25))));
     
     for (int i=pos[0]; i>=pos[1]; i--)
     {
         	Tsurface[i] = W[v-1] + ( i-pos[0] )*(W[v-1] - W[v])/(pos[0] - pos[1]);
         	if (int(tnew[v]/((dt/(24.0*3600.0*365.25)))) == i)
         	{
          		iter = i;
          	}
          	
     }
     //alpha =5E-5;
     Wnew[v] = Tsurface[iter] + (u2)*alpha;
       
     	
 	if (Wnew[v] < -10.0 || Wnew[v] > 10.0)
 	{
 		cerr << "outside limit " ;
 		Wnew[v] = W[v];
 	}
 	if (Wnew[v+1] < -10.0 || Wnew[v+1] > 10.0)
 	{
 		cerr << "outside limit "  << endl;
 		Wnew[v+1] = W[v];
 	}
	 
	 
	
//	 		   cerr << "tnew " << endl;
//		   for (int i=0; i<newpoints; i++)
//		   {cerr << tnew[i] << "\t" ;}
//		   cerr << endl;
		   
//		     cerr << "Wnew " << endl;
//		   for (int i=0; i<newpoints; i++)
//		   { cerr << Wnew[i] << "\t" ;}
//		   cerr << endl;
          
	
 	return 0;
        
         
	}


//--------------------------------------------------------------------------------------------------------
vector<double>RJ_SS::Death(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> told,
				double tmax,double newT, int points, int& newpoints, int v , double& wminus, double& wplus,
				   double& u1, double& u2 , double& alpha )
	{ 
	
	       
	
//          cerr << "told " << endl;
//		   for (int i=0; i<points; i++)
//		   { cerr << told[i] << "\t" ;}
//		   cerr << endl;
		   
//		    cerr << "W " << endl;
//		   for (int i=0; i<points; i++)
//		   { cerr << W[i] << "\t" ; }
//		   cerr << endl;
		  
		    
		     for (int i=0; i<v; i++)
		     {
		         tnew[i] = told[i];   
		         
		     }
		     
		     for (int i=v; i<newpoints; i++)  
		     {
		     	tnew[i] = told[i+1];
		     	
		     }
		    
		     for (int i=0; i<v; i++)
		     {
		      	Wnew[i] = W[i] ;
		      	
		     }
		     for (int i=v; i<newpoints; i++)  
		     {
		     	Wnew[i] = W[i+1];
		        
		     }
		     wminus = told[v]  - told[v+2];
//		   cerr << "tnew " << endl;
//		   for (int i=0; i<newpoints; i++)
//		   {cerr << tnew[i] << "\t" ;}
//		   cerr << endl;
		   
//		     cerr << "Wnew " << endl;
//		   for (int i=0; i<newpoints; i++)
//		   { cerr << Wnew[i] << "\t" ;}
//		   cerr << endl;
		 
		   //alpha =5E-5;
		   return tnew;
	}
   
//--------------------------------------------------------------------------------------------------------
int RJ_SS::Peturb(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& v  )
{
	         vector<double>::iterator q; 
              
             GaussianRNG c;      
             tnew.assign(points,0.0);
             Wnew.assign(points,0.0);
             Wnew = W;
             tnew = times;
             double u1 = c.generate();
             Array1D<double> alpha = Array1D<double> (newpoints,0.0);
	            
             for (int ik=0; ik<newpoints ; ik++)
	            {
	            	
	            	alpha[ik] = 0.2;
	            	//alpha[ik]  =  0.4*exp(double(ik+1)/newpoints);			//0.4 9.8.2009
	            //	cerr << alpha[ik] << "\t"; 
	            }
	            //cerr << endl; 
             
             double add = (tnew[v-1] - tnew[v+1])*u1*alpha[v];
          
             if (  ( (tnew[v] + add) > (tnew[v-1]-(dt/(24.0*3600.0*365.25))) || ((tnew[v] + add) < (tnew[v+1]+((dt/(24.0*3600.0*365.25)))))) )
             { 
             //	cerr <<"PETURB time TOO CLOSE" << endl;
             //	cerr <<tnew[v-1] << "\t" << tnew[v] + add << "\t" << tnew[v+1] << endl;
             	add = 0.0;
             }
             
             tnew[v] = tnew[v] +   add;
             
             return 0;
}
		
		
//--------------------------------------------------------------------------------------------------------
int RJ_SS::Peturb_T(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& element, bool& half  )
{
	       
	        bool MVG =false;
	        int u;
		    int v;
		    
	        if (MVG == false)
	        {
	        	
		         vector<double>::iterator q; 
	             GaussianRNG c;     
	             tnew.assign(points,0.0);
	             Wnew.assign(points,0.0);
	             Wnew = W;
	            // for (int i=0; i<points; i++)
	            // {
	           //  	cerr << Wnew[i] << "\t" ;
	           //  }
	           //  cerr << endl;
	             tnew = times;
	             double u = c.generate();
	             
	             Array2D<double> mm = Array2D<double> (iterations+1, 2, 0.0);
	         
		   		Array1D<double> alpha = Array1D<double> (newpoints,0.3500);
	            
	            for (int ik=0; ik<newpoints ; ik++)
	            {
	                alpha[ik]  = 0.05;
		    	//alpha[ik]  =  0.1*exp(double(ik+1)/newpoints);			 //change to 0.1 9.8.2009
	             //	cerr << alpha[ik] << "\t"; 
	            }
	            //cerr << endl;
	            Wnew[element] =Wnew[element] + ((u)* alpha[element]); 
		   		 
		   		 if (Wnew[element] <-10.0 || Wnew[element] > 10.0)
		   		 {
		   		 	cerr <<Wnew[element] << " is outside range at v = " << element<< endl;
		   		 	Wnew[element] -= (u*alpha[element]);
		   		 }
	        }
	        //==============================================================
	        else
	        {
	        	
	        	UniformRNG HALF;
	        	double dice = HALF.generate();
	        	if (dice < 0.5)
	        	{
	        		half = true;
	        	}
	        	else
	        	    half = false;
	        	
		        if (half ==false)
		        {
		        	u = 0; 
		        	v = int(double(newpoints)/2.0);
		        	
		        }
		        if (half == true)
		        { 
		        	u = int(double(newpoints)/2.0);
		        	v= newpoints;
		        	
		        }
		       
		        // cerr <<"u "<< u << "\t" << "v " << v << endl;
		        //Multi-variate-Gaussian update:
		        UniformRNG aa;
		        vector<double>::iterator q; 
	            
	            tnew.assign(points,0.0);
	            Wnew.assign(points,0.0);
	            Wnew = W;
	            tnew = times;
	            reverse(Wnew.begin(), Wnew.end());
	            reverse(tnew.begin(), tnew.end());
	            reverse(W.begin(), W.end());
	            reverse(times.begin(), times.end());
	            
	            Array1D<int> TIME = Array1D<int> (newpoints, 0);
	            Array2D<double> Chol = Array2D<double> (150, 150, 0.0);
	            Array2D<double> chol = Array2D<double> (newpoints, newpoints, 0.0);
	            ifstream inC("chol.txt");
	            for (int i=0; i<150; i++)
	            {
	            	for (int j=0; j<150; j++)
	            	{
	            		inC >> Chol[i][j];
	            	}
	            }
	            //cerr <<"tnew[0] "<< tnew[0] << endl;
	            for (int i=0; i<newpoints; i++)
	            {
	            	TIME[i] =int(tnew[i]/50.0);	
	            }
	            TIME[newpoints-1] -=1;
	            
	            for (int i=0; i<newpoints; i++)
	            {
	            	for (int j=0; j<newpoints; j++)
	            	{
	            		chol[i][j] = Chol[TIME[i]][TIME[j]];
	                }	               
	            }	
	            
	        	//Proposal
	            Array1D<double> y = Array1D<double> (newpoints, 0.0);
	            Array1D<double> alpha = Array1D<double> (newpoints,0.300);
	            for (int i=0; i<newpoints ; i++)
	            {
	            	alpha[i]  =  0.05*exp(0.25*double(i)/newpoints);
	            }
	          
	             
	            GaussianRNG a; 
	            for (int ui=0; ui<points; ui++)
			    {    
			         y[ui] = a.generate()* alpha[ui];
			    }
			   
			    for (int i_ch=u; i_ch<v; i_ch++)
			    {  
			   	    for (int j_ch=0; j_ch<newpoints; j_ch++)
			        { 
			            Wnew[i_ch] += chol[i_ch][j_ch]*y[j_ch] ;
			            
			        }
			        
			        if (Wnew[i_ch] >15.0|| Wnew[i_ch] < 5.0) 
			           {    cerr <<"bad proposal" << endl;
			         	    Wnew[i_ch] = W[i_ch];} 
			    }	
			    
	            reverse(Wnew.begin(), Wnew.end());
	            reverse(tnew.begin(), tnew.end());
	            reverse(W.begin(), W.end());
	            reverse(times.begin(), times.end());
	        } 
            
		    return 0;
}
//--------------------------------------------------------------------------------------------------------
int RJ_SS::Peturb_H(double& HtFnew, double& HtF, double& TEQnew, double& Teq, 
                       double alpha_H, double alpha_Teq, int a, int b )
{		
		 
	GaussianRNG up;
	//	 for (int i=hu; i<hv; i++)
	//	 {	
	//	 	HtFnew[i] = HtF[i] + Update.generate()*alpha_H;
	//	 	if (HtFnew[i] > 0.125    || HtFnew[i] < 0.005)
	//	 	{
	//	 		HtFnew[i] = HtF[i];
	//	 	}
	//	 }
		 
		 
	//	 for (int i=eu; i<ev; i++)
	//	 {	
	//	 	TEQnew[i] = Teq[i] + Update.generate()*alpha_Teq;
	//	    if (TEQnew[i] > 15.0    || TEQnew[i] < 0.0)
	//	 	{
	//	 		TEQnew[i] = Teq[i];
	//	 	}
	//	 }
		 Array1D<double> alpha = Array1D<double> (2, 0.0);
		 alpha[0] = alpha_H;
		 alpha[1] = alpha_Teq;
		 Array2D<double> chol = Array2D<double> (2,2, 0.0);
		 Array1D<double> update = Array1D<double> (2, 0.0);
		 chol[0][0]  = 1.0;
		 chol[0][1] = 0.0;
		 chol[1][1]=  0.6614;
		 chol[1][0] = -0.7500;
	     	     
	     for (int i=0; i<2; i++)
	     {
	     	for (int j=0; j<2; j++)
	     	{
	     		update[i] += chol[i][j] * up.generate() * alpha[i];
	     	}
	     	
	     }
	     
	     HtFnew = update[0] + HtF;
	     TEQnew	   = update[1] + Teq;  	
	     //cerr << "HtFnew " << HtFnew <<"\t old:" << HtF << endl;
	     // cerr << "TEQnew " << TEQnew <<"\t old:" << Teq << endl;	
return 0;	 
}
		 
		 
		 
//--------------------------------------------------------------------------------------------------------

int RJ_SS::Prob(int points, int& x, double& JR, vector<double> tmcmc)
{           
	 int maxPoints =40;
	 double total_options;
         double min_total_options;
	 double max_total_options;
	 UniformRNG c;
     UniformRNG d;               
     //which element to perturb 
     
 
    
     
     // what RJ?
     total_options=7.0;
 	 double X = (total_options)* d.generate();
    
     JR = 1.0;
     if (X <1.00000) 
       {x =1;}					// Birth
     if (X >=1.00000 && X<2.00000) 
       {x =2; }					//Death
     if (X >= 2.00000 && X <3.00000) 	
       { x = 3;}				//perturb the time point
     if (X >= 3.00000 && X <=4.00000) 
       { x = 4;}				// Perturb the temperature value
     if (X > 4.00000 && X <=5.00000) 
       { x = 5;}  				// perturb both the heatflow and the Teq
     if (X > 5.00000 && X <=6.00000) 
       { x = 6;}  				// perturb prior width on temperature history
      if (X > 6.00000) 
       { x = 7;}   				// perturb width on likelihood function
       //cerr <<"x " << x << endl;
       
 	 if (points ==2)  // if only 2 points, force Birth or perturb Temperature or perturb Heat flow:
       {   
      	    min_total_options=2.0;
	    X = min_total_options* d.generate();
      	   
      	  if (X <1.00000) 
          {x =1;		//Birth step
            JR = (1.0/total_options)/(1.0/min_total_options);
	    }
        if (X >=1.00000 && X<2.00000) 
           {x =4;          // Perturb temperature
           JR = 1.0;}
          if (X >=2.00000 && X<3.00000) 
           {x =5; 		// Perturb heatflow/Teq
           JR = 1.0;}
          
    	}
     // if max points, force death or pertur T or t or perturb H:
     if (points>maxPoints){ 
    	 max_total_options=3.0;
	 X = max_total_options* d.generate();
 	 if (X <1.00000) 
       {x =2;
       JR = (1.0/total_options)/(1.0/max_total_options);
       }
     if (X >=1.00000 && X<2.00000) 
       {x =3;
       	JR = 1.0; }
     if (X >= 2.00000 && X <3.00000) 
       { x = 4;
       JR = 1.0;}	
 	 if (X >= 3.00000 && X <4.00000) 
       { x = 5;
       JR = 1.0;}	
 		 
 		 
 	 }
       
     return x;
}
//---------------------------------------------------------------------------------------------------------
int RJ_SS::Element(int points, int x, int& v, double& JR, vector<double> tmcmc)
{  
	 
	 
	 BirthError:
	 UniformRNG c;
	 double V = 0.0;
	 
	 if (x ==1)
	 {
	 	V =  c.generate();   
     	V = V*(points-1);
     	v = int(V);
     	v = v+1;
	 }
	 if (x ==2 || x==3)
	 {
	 	V =  c.generate();   
     	V = V*(points-2);
     	v = int(V);
     	v = v+1;
	 }
	 if (x ==4)
	 {
	 	V =  c.generate();   
     	V = V*(points);
     	v = int(V);
     	v = v;
	 }
     //cerr << "v " << v << endl; 
    // cerr << "x " << x << endl;
     if ( x == 1 )
     {	if (abs(tmcmc[v] - tmcmc[v-1]) < (dt/(24.0*3600.0*365.25)))
     	{
     		cerr << "STEPS TOO CLOSE for BIRTH " << endl;
     		cerr << tmcmc[v] << "\t" << tmcmc[v-1] << endl;
     		
     		goto BirthError;
     	}
     }
     	
     if ( x == 3 )
     {	if (abs(tmcmc[v+1] - tmcmc[v-1]) <(dt/(24.0*3600.0*365.25)))
     	{
     		cerr << "STEPS TOO CLOSE for PETURB time " << endl;
     		cerr << tmcmc[v+1] << "\t" << tmcmc[v-1] << endl;
     		cerr << "dt "<< dt << endl;
     		cerr << "dt/dt'" << ((dt/(24.0*3600.0*365.25)))<<endl;
     		goto BirthError;
     	}
     }
	return v;
	
}


//--------------------------------------------------------------------------------------------------------
Array1D<double> RJ_SS::Interpolate_SS(vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints )
{
	vector<double>::iterator  q ;
		
        //for (q=Wmcmc.begin(); q!=Wmcmc.end(); q++) 
     //        {
     //        	cerr << *q  <<"\t";
      //       	
     //        }
    //         cerr << endl;
	//for (q=tmcmc.begin(); q!=tmcmc.end(); q++) 
   //          {
    //         	cerr << *q  <<"\t";
  //           }
    //         cerr << endl;
	vector<double>::iterator qq;

	
	Array1D<double> time = Array1D<double>(iterations, 0.0);
	vector<double> Tsurf;
	Array1D<double> Tsurface = Array1D<double>(iterations, 0.0);
    Tsurf.assign(iterations, 0.0);
   
    Tsurface = Array1D<double> (iterations, 0.0);
    Array1D<int> pos = Array1D<int> (NEWpoints, 0);
    
    for (int i=0; i<NEWpoints; i++)
    {
    	pos[i]= roundnew(tmcmc[i]/(dt/(24.0*3600.0*365.25)));
    	//cerr <<"pos "<< pos[i] << endl;
    }
    
    for (int j=0; j<NEWpoints-1; j++)
    {
       for (int i=pos[j]; i>pos[j+1]; i--)
       {
          	//cerr << "i "<< i << endl;
          	Tsurface[i-1] = Wmcmc[j] + ( i-pos[j])*(Wmcmc[j] - Wmcmc[j+1])/(pos[j] - pos[j+1]);
       }
    }   
    Tsurface[iterations - 1] = Wmcmc[0];
    
    for (int iter=0; iter<iterations; iter++)
    {
    	Tsurf[iter] = Tsurface[iter];
        
    }
    
    vector<double>::iterator jac;
    /*ofstream outTs("Out/Tsurf.txt");
    for (jac=Tsurf.begin(); jac!=Tsurf.end(); jac++)
    {
    	outTs << *jac << endl;    	
    }
    outTs.close();*/
   
    int a = Tsurf.size();
    Array1D<double> error =Array1D<double> (1, -999.99);
    if (a !=iterations)
    {
    	cerr << "Interpolation error "<< endl;
    	return error;
    }
    
    for (int i=0; i<NEWpoints; i++)
    {
    	if ((Wmcmc[i] <-15.0)  || (Wmcmc[i]>45.0))
    	{cerr << "W error "<< Wmcmc[i] << endl;
    		// Comment out for prior test:
		//return error
		;}
    	
    }
    return Tsurface;
}
//-------------------------------------------------------------------------------------------------------
int roundnew(double a) 
{
return int(a + 0.5);
}


