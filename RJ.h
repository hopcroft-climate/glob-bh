int factorial(int number) ;

class RJ
{
		public:
			double accept (double LL,double LLold,double Prold, double JR,bool BIRTH,bool DEATH,bool PERTURB,bool MOVE, 
                  double JRatio, int nC, int nCold, double length, double mu, vector<double> beta_new,vector<double> beta);
			
			void interpolate (vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints, 
                   Array2D<double>& Tsurface, int nn);
		    
		    void choose(int points, int& x, double& JR, int nBH );
		
};

double RJ::accept(double LL,double LLold,double Prold, double JR,bool Birth,bool Death,bool Perturb,bool Move, 
                    double JRatio, int nC, int nCold, double length, double mu, vector<double> beta_new,vector<double> beta)
{

	
	//double PriorRatio = 1.00;		//Prior ratio
  //  double J= 1.00; 				//det(Jacobian)
  //  double PropRatio = 1.0;			//Proposal ratio
    
   // double CpriorR = 1.0;
      
       
       
       
    
    //double L = tmcmc[0] -tmcmc[Points-1];
    
     if (Birth == true) 
       {
       	  //PropRatio = (nC)/(1.0*1.0);
       	//  PropRatio =  1.0/(double(nC)/double(length));
       	 
       } 
       if (Death == true)
       {
       //	   PropRatio = (double(nCold))/(double(length));
       	   //PropRatio = 1/PropRatio;
	//	   cerr << "length " << length << endl;
       }
       if (Perturb ==true)
       { 
       	//  PropRatio = length;
	//	  cerr << "Perturb PropRatio " << PropRatio << endl;
       }
       
       if (Move ==true)
       { ;
       }
      
       double alpha = 0.0;
       //alpha =(LL-LLold)+ log(Prior)+ log(PriorRatio) + log(JRatio) + log(J)  - log(prold);
      
       // Account for the probability of choosing a certian GST
        double betaProp = 1.0;
        for (int i=0; i<nCold; i++)
        { 
        	betaProp = betaProp*beta[i];
        }
        double betaProp_new = 1.0;
        for (int i=0; i<nC; i++)
        { 
        	betaProp_new = betaProp_new*beta_new[i];
        }
        
        betaProp= log(betaProp);
        betaProp_new= log(betaProp_new);
        cerr << "Marg Post.new " << betaProp_new <<"\t old " << betaProp << endl;
        
        alpha =(LL-LLold) - (betaProp_new - betaProp)  +  log(JRatio);  ;
        
       // cerr <<"J " << log(J) << endl;
        cerr <<"JumpRatio " << log(JRatio) << endl;
        //cerr <<"ProprRatio " << log(PropRatio) << endl;
        cerr <<"LL " << LL << endl;
        cerr <<"LLold " << LLold << endl;
        
	
       return alpha;

}

//==================================================================================================


 void Interpolate(vector<double>& Wmcmc, vector<double> tmcmc, int NEWpoints, 
                   Array2D<double>& Tsurface, int nn )
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

	
	//Array1D<double> time = Array1D<double>(iterations, 0.0);
	vector<double> Tsurf;
	//Array1D<double> Tsurface = Array1D<double>(iterations, 0.0);
    Tsurf.assign(iterations, 0.0);
   
    //Tsurface = Array1D<double> (iterations, 0.0);
    Array1D<int> pos = Array1D<int> (NEWpoints, 0);
    
    for (int i=0; i<NEWpoints; i++)
    {
    	pos[i]= int(tmcmc[i]/(dt/(24.0*3600.0*365.25)));
    	//cerr <<"pos "<< pos[i] << endl;
    }
    
    for (int j=0; j<NEWpoints-1; j++)
    {
       for (int i=pos[j]; i>pos[j+1]; i--)
       {
          	//cerr << "i "<< i << endl;
          	Tsurface[i-1][nn] = Wmcmc[j] + ( i-pos[j])*(Wmcmc[j] - Wmcmc[j+1])/(pos[j] - pos[j+1]);
       }
    }   
    Tsurface[iterations - 1][nn] = Wmcmc[0];
    
    for (int iter=0; iter<iterations; iter++)
    {
    	Tsurf[iter] = Tsurface[iter][nn];
        
    }
    
    vector<double>::iterator jac;
    ofstream outTs("Out/Tsurf.txt");
    for (jac=Tsurf.begin(); jac!=Tsurf.end(); jac++)
    {
    	outTs << *jac << endl;
    	
    }
    outTs.close();
   
    int a = Tsurf.size();
    //Array1D<double> error =Array1D<double> (1, -999.99);
   // int error = 0;
    if (a !=iterations)
    {
    	cerr << "Interpolation error "<< endl;
    	//return error;
    }
    
    for (int i=0; i<NEWpoints; i++)
    {
    	if ((Wmcmc[i] >-10.0)  && (Wmcmc[i]< 25.0))
    	{
    		;
    	}
    	else
    	{cerr << "W error but do nothing "<< endl;
    		//return error;
    		}
    	
    }
          	   
   
    return ;
}


//===========================================================================================================

void RJ::choose(int points, int& x, double& JR, int nBH )
{
     
	 UniformRNG c;
     UniformRNG d;               
      
     // what RJ?
 	 double X = (6)* d.generate();			
    
     JR = 1.0;
     if (X <1.00000) 
       {x =1;}				//Birth
     if (X >=1.00000 && X<2.00000) 
       {x =2; }				//Death
     if (X >= 2.00000 && X <3.00000) 
       { x = 3;}		    //Perturb			
     if (X >=3.00000 && X <4.00000)		
       { x = 4;}			//Move
     if (X >=4.00000 && X <5.00000)		
       { x = 5;}			// UPdate heat flow pair
     if (X >= 5.00000)	
       { x = 6;}			//GST update
     
       
 	 if (points ==1)  // if only 1 tessellation, force Birth or perturb/move/subSample
       {   
      	   X = 4* d.generate();
      	   
      	  if (X <1.00000) 
          {	x =1;
          	JR = 0.25/0.2;
          }
          if (X >=1.00000 && X<2.00000) 
           {x =3; 
           JR = 1.0;}
           if (X >=2.00000 && X<3.00000) 
           {x =4; 
           JR = 1.0;}
           if (X >=3.00000 && X<4.00000) 
           {x =6; 
           JR = 1.0;}
          
    	}
     // if max points, force death or perturb
     if (points==(nBH))
     { 
     	 
         X = 4* d.generate();
 	     if (X <1.00000) 
         {   x =2;
             JR = 0.25/0.2;
         }
        if (X >=1.00000 && X<2.00000) 
            {x =3;
       	     JR = 1.0; }
        if (X >= 2.00000 && X <3.00000) 
           { x = 4;
           JR = 1.0;}	
 		 if (X >= 3.00000 && X <4.00000) 
           { x = 6;
           JR = 1.0;}
 		 
 	 }
       

}

int factorial(int number) 
{
	int temp;
	if(number <= 1) return 1;
    temp = number * factorial(number - 1);
	return temp;
}
		
