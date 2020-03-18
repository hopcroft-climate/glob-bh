// generic class for a RJMCMC update object

class RJ
    { 
          public:
              int SetPoints(int points);      
                
             int Birth(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double>& told,
				double tmax,double newT, int newpoints , int& v,double& wminus, double& wplus );
				
     		vector<double>Death(vector<double> W,vector<double>& Wnew,vector<double>& tnew,vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int v ,double& wminus, double& wplus );
				
			 int Peturb(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& v  );
		    int Prob(int points, int x, int& v);
		    
		    Array1D<double>Interpolate(vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints );
           
            
          
          private:
			  int k; // no of points in Tsurf
			 

                   
    };

//--------------------------------------------------------------------------------------------------------
    int RJ::SetPoints(int points)
    {      k =points; 
    	return 0;  
    }

    
//--------------------------------------------------------------------------------------------------------
	int RJ::Birth(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double>& told,
				double tmax,double newT, int newpoints , int& v, double& wminus, double& wplus )
	{
     
    // for (int i=0; i<newpoints-1; i++)
     //	{
     //		cerr << "Wold " << W[i] << "\t" ;
     		
     //	}
     	//cerr << endl;
     	//for (int i=0; i<=newpoints-1; i++)
     	//{
     	//	cerr << "told " << told[i] << "\t" ;
     		
     	//}
     
     if (v >0)  
     {for (int i=0;i<=v; i++)
         {
         	tnew[i] = told[i] ;
          }
     }
     else
     {
      	 tnew[0] = told[0];
     }
     for (int i=v+2; i<=newpoints; i++)
     {
     	tnew[i] = told[i-1];
     	
     }
     UniformRNG ccU;
     double u1 = ccU.generate();
     tnew[v+1] = tnew[v+2] + u1* (told[v]- told[v+1]);      
     wminus = tnew[v+1] - tnew[v+2];
     wplus = tnew[v] - tnew[v+1];
     for (int i=0; i<=v; i++)
     {
     	Wnew[i] = W[i];
     }
     for (int i=v+1; i<newpoints; i++)
     {
     	Wnew[i] = W[i-1];
     }
     GaussianRNG ccG;
     double u2 = ccG.generate();
     Wnew[v] = W[v] + u2*0.5;
     Wnew[v+1] = W[v] - u2*0.5;
   
     	if (Wnew[v] < 5.0 || Wnew[v] > 15.0)
     	{
     		cerr << "outside limit " ;
     		Wnew[v] = W[v];
     	}
     	if (Wnew[v+1] < 5.0 || Wnew[v+1] > 15.0)
     	{
     		cerr << "outside limit "  << endl;
     		Wnew[v+1] = W[v];
     	}
     	
    // 	for (int i=0; i<newpoints; i++)
    // 	{
     //		cerr << "Wnew " << Wnew[i] << "\t" ;
     		
    // 	}
     //	cerr << endl;
    /// 	for (int i=0; i<=newpoints; i++)
    // 	{
     //		cerr << "tnew " << tnew[i] << "\t" ;
    // 		
    // 	}
     	return 0;
        
         
	}


//--------------------------------------------------------------------------------------------------------
vector<double>RJ::Death(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times,
				double tmax,double newT, int points, int& newpoints, int v , double& wminus, double& wplus )
	{ 
		     
		     if (v==0)
		     {
		     	tnew[v] = times[v];
		     }
		     else{
		        for (int i=0; i<v; i++)
		        {
		           tnew[i] = times[i];   
		        }
		     }
		     for (int i=v+1; i<=points; i++)  
		     {
		     	tnew[i] = times[i+1];
		     }
		    
		     for (int i=0; i<=v; i++)
		     {
		      	Wnew[i] = W[i] ;
		      	
		     }
		     for (int i=v+1; i<newpoints; i++)  
		     {
		     	Wnew[i] = W[i];
		        
		     }
		     Wnew[v] = ((times[v] - times[v+1])*W[v] + (times[v+1] - times[v+2])*W[v+1])/(times[v] - times[v+2]);
		     wplus = (times[v] - times[v+1]);
		     wminus =(times[v+1] - times[v+2]); 
		     return tnew;
	}
   
//--------------------------------------------------------------------------------------------------------
int RJ::Peturb(vector<double> W,vector<double>& Wnew,vector<double>& tnew, vector<double> times, 
				double tmax,double newT, int points, int& newpoints, int& v  )
{
	         
	        
	         vector<double>::iterator q; 
             
	           
	         UniformRNG cc;
             UniformRNG d;               
             
             UniformRNG c;     
             tnew.assign(points+1,0.0);
             Wnew.assign(points,0.0);
             Wnew = W;
             tnew = times;
            
             tnew[v+1] =tnew[v+2]+ (tnew[v+1]-tnew[v+2])*(c.generate());
            
		     return 0;
}
		
		
//--------------------------------------------------------------------------------------------------------
int RJ::Prob(int points, int x, int& v )
{           
	
	 UniformRNG c;
     UniformRNG d;               
     //which element to use
     double V =  c.generate();   
     V = V*(points-1);
     v = int(V);
     
     // what RJ?
 	 double X = 3* d.generate();
    
    
     if (X <1.00000) 
       {x =1;}
     if (X >=1.00000 && X<2.00000) 
       {x =2; }
     if (X >= 2.00000) 
       { x = 3;}	
 	 if (points ==1)  // if only 2 points, force Birth
       {   
      	  x= 2;
    	}
    	
     if (points>14){ 
    	if(d.generate() <0.5) {	x = 3;}
    	else x = 1;
 		 }
       
     return x;
}

//--------------------------------------------------------------------------------------------------------
Array1D<double> RJ::Interpolate(vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints )
{
	
	vector<double>::iterator qq;
    
	
	Array1D<double> time = Array1D<double>(iterations, 0.0);
	vector<double> Tsurf;
	Array1D<double> Tsurface = Array1D<double>(iterations, 0.0);
    Tsurf.assign(iterations, 0.0);
   
    Tsurface = Array1D<double> (iterations+1, 0.0);
    Array1D<int> pos = Array1D<int> (NEWpoints+1, 0);
    
    for (int i=0; i<=NEWpoints; i++)
    {
    	pos[i]= int(tmcmc[i]/10);
    }
    
    for (int j=0; j<NEWpoints; j++)
    {
       for (int i=pos[j]; i>pos[j+1]; i--)
       {
          	Tsurface[i] = Wmcmc[j];   
       }
    }   
    
    for (int iter=1; iter<=iterations; iter++)
    {
    	Tsurf[iter-1] = Tsurface[iter];
    	
    }
    reverse(Tsurf.begin(),Tsurf.end());
    vector<double>::iterator jac;
    ofstream outTs("Out/Tsurf.txt");
    for (jac=Tsurf.begin(); jac!=Tsurf.end(); jac++)
    {
    	outTs << *jac << endl;
    	
    }
    outTs.close();
   
    int a = Tsurf.size();
    Array1D<double> error =Array1D<double> (1, -999.99);
    if (a !=iterations)
    {
    	
    	
    	return error;
    }
    Tsurface[0] = Tsurface[1];
    for (int i=0; i<iterations; i++)
    {
    	if (Tsurface[i] <1.0)
    	{
    		return error;
    	}
    }
          	   
    return Tsurface;
}



