class Voronoi
    { 
      public:
         double get_dist(double Cx,double Cy, double Bx,double By, int nC);
         int birth(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC, double& length);
         int death(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC,int v, double& length);
         int perturb(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC,double& length);
         int move(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC, int v);
              
         int which_region(vector<double> Cx, vector<double> Cy,double BHx, double BHy, int nC);
         int any_data(vector<int>& GST_id_new,vector<int>& GST_id, vector<GST>& Uknew, vector<GST>& UKtnew, vector<GST>& UK, 
                                            vector<GST>&  UKt, int& nCnew, int nC, int nBH,
                         vector<double>& Cx_new, vector<double>& Cy_new,vector<double>& Cx,vector<double>& Cy, vector<double>&  beta, vector<double> & beta_new );
          private:
			  int k; // no of points in Tsurf
			 
 
                   
    };
    
double Voronoi::get_dist(double Cx,double Cy, double Bx,double By, int nC)
{
	double xd = (Cx - Bx)*(Cx- Bx);
	double yd = (Cy - By)*(Cy- By);
	
	double d = xd + yd;
	d = sqrt(d);
	
	
	return d;
}
//==============================================================================================================

int Voronoi::birth(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, 
												int nCnew, double& length)
{
    
    
   // if ((nCnew -1) > 1)
  //  {
	//	    UniformRNG where;
	//	    double Where1 = where.generate()*(nCnew-1);
	//	   // double Where2 = where.generate()*(nCnew-1);
		    
	//	    int where1 = int(Where1);
	//	   // int where2 = int(Where2);
		    
	//	   double Wherex = Cx[Where1];
	//	   double Wherey = Cy[Where1];
		   
	//	   Cx_new[nCnew-1] = Wherex + where.generate()*0.05;
	//	   Cy_new[nCnew-1] = Wherey + where.generate()*0.05;
	//	   length = 0.05*0.05;
 //   }
 //   else
 //   {
     	UniformRNG Cgen;
    	double x1 = Cgen.generate();
    	double y1 = Cgen.generate();
    	Cx_new[nCnew-1] = x1;
    	Cy_new[nCnew-1] = y1;
    
    	length = 1;
    	
  //  }
    //cerr <<"Cxnew, CYnew sizes: " << Cx_new.size() << "\t" << Cy_new.size() << endl;
    
    for (int i=0; i<(nCnew-1); i++)
	{
		Cx_new[i] = Cx[i];
		Cy_new[i] = Cy[i];
	}
	
	
	
	
	
	return 0;
	
}
//==============================================================================================================

int Voronoi::death(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nCnew, int v, double& length)
{
    //cerr << "Centres original" ;
   // for (int i=0; i<nCnew+1; i++)
   // {
   // 	cerr << Cx[i] << "\t"<< Cy[i] << endl;
   // }
   // cerr <<"end" << endl;
    	
    
    
    Cx_new.resize(nCnew);
    Cy_new.resize(nCnew);
    
    for (int i=0; i<v; i++)
    { 
    	Cx_new[i] = Cx[i] ;
    	Cy_new[i] = Cy[i];
    }
    for (int i=v+1; i<(nCnew+1); i++)
    { 
    	Cx_new[i-1] = Cx[i] ;
    	Cy_new[i-1] = Cy[i];
    }	
    length = 1.0;
    
    
  //  cerr <<"new centres " ;
   // for (int i=0; i<nCnew; i++)
  //  {
  //  	cerr << Cx_new[i] << "\t"<< Cy_new[i] << endl;
  //  }
	//cerr <<"end" << endl;
	
	return 0;
	
}
//==============================================================================================================

int Voronoi::perturb(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC, double& length)
{
    
        Cx_new = Cx;
    Cy_new = Cy;
    double st_val = 0.50;
    UniformRNG Cperturb;
    double V = Cperturb.generate()*(nC*2);
    int v = int(V);

	double sigma_x = st_val;
	double sigma_y = st_val;

	GaussianRNG XYpert;

    if (v >= nC)
	{
		 v-= nC;
		// Propose new positions x
		double u1 = XYpert.generate()*sigma_x;
		Cx_new[v] += u1;
	    // is the proposal outside the boundaries ? prpb = 0.02
		if (Cx_new[v] > 1.0 || Cx_new[v]< 0.0 )
		{
			cerr <<"PERTURB:outside space x" << endl;
			Cx_new[v]  -=u1;
			
		}
	}
	else
	{
		// Propose new positions y
		double u2 = XYpert.generate()*sigma_y;
	    Cy_new[v] += u2;
	    // is the proposal outside the boundaries ? prpb = 0.02
		if (Cy_new[v] > 1.0 || Cy_new[v]< 0.0 )
		{
			cerr <<"PERTURB:outside space y" << endl;
			Cy_new[v]  -=u2;
			
		}
	}
	length = 1.0;
   
   /* Cx_new = Cx;
    Cy_new = Cy;
    double st_val = 0.15;
    UniformRNG Cperturb;
    double V = Cperturb.generate()*(nC);
    int v = int(V);
    double sigma_x = st_val;
    double sigma_y = st_val;
	double s_x_old = st_val;
	double s_y_old = st_val;
    //double r2 = (sigma_x*sigma_x + sigma_y*sigma_y);
    // If proposed move is too large for the area
    
    
	
	
	// is the current position close to the boundary?
	if (Cx[v] + 3*sigma_x >= 1.0) 
	{
		s_x_old = (1.0- Cx_new[v])/3.0;
		
	}
	if (Cx[v] - 3*sigma_x <=0.0) 
	{
		s_x_old = ( Cx_new[v])/3.0;
	}
    
    if (Cy[v] + 3*sigma_y >= 1.0)
	{
		s_y_old = (1.0- Cy_new[v])/3.0;
	}
	if (Cy[v] - 3*sigma_y <= 0.0) 
	{
		s_y_old = ( Cy_new[v])/3.0;
	}
    
	// Propose new positions
	GaussianRNG XYpert;
    double u1 = XYpert.generate()*s_x_old;
    double u2 = XYpert.generate()*s_y_old;
    Cx_new[v] += u1;
    Cy_new[v] += u2;
	

	// is the new position close to a boundary
     if (Cx_new[v] + 3*sigma_x >= 1.0) 
	{
		sigma_x = (1.0- Cx_new[v])/3.0;
		
	}
	if (Cx_new[v] - 3*sigma_x <=0.0) 
	{
		sigma_x = ( Cx_new[v])/3.0;
	}
    
    if (Cy_new[v] + 3*sigma_y >= 1.0)
	{
		sigma_y = (1.0- Cy_new[v])/3.0;
	}
	if (Cy_new[v] - 3*sigma_y <= 0.0) 
	{
		sigma_y = ( Cy_new[v])/3.0;
	}

	// is the proposal outside the boundaries ? prpb = 0.02
	if (Cx_new[v] > 1.0 || Cx_new[v]< 0.0)
	{
		cerr <<"PERTURB:outside space x" << endl;
		Cx_new[v]  -=u1;
		sigma_x = st_val;
	}
	if (Cy_new[v] > 1.0 || Cy_new[v]< 0.0)
	{
		cerr <<"PERTURB:outside space y" << endl;
		Cy_new[v]  -=u2;
		sigma_y = st_val;
	}	
	
	

	// calculate the proposal probabilities
	// q(m  | m')	prob of coming back to initial position
	double pm = (1/sigma_x)*(1/sigma_y)* exp(-(Cx_new[v] - Cx[v])*(Cx_new[v] - Cx[v])*(1/(sigma_x*sigma_x)) + (Cy_new[v] - Cy[v])*(Cy_new[v] - Cy[v])*(1/(sigma_x*sigma_x)));
    // q(m' | m) probability of proposing new position
	double pm_old = 	(1/s_x_old)*(1/s_y_old)* exp(-(Cx_new[v] - Cx[v])*(Cx_new[v] - Cx[v])*(1/(s_x_old*s_x_old)) + (Cy_new[v] - Cy[v])*(Cy_new[v] - Cy[v])*(1/(s_y_old*s_y_old)));
		
	 length = pm/pm_old;
	 cerr << "PERTURB Proposal RATIO: " << length << endl;
    */
	return 0;
	
}
//==============================================================================================================

int Voronoi::move(vector<double>& Cx_new, vector<double>& Cy_new,vector<double> Cx, vector<double> Cy, int nC, int v)
{
    
    Cx_new = Cx;
    Cy_new = Cy;
    
    
    
    UniformRNG XYpert;
    double u1 = XYpert.generate();
    double u2 = XYpert.generate();
    
    Cx_new[v] = u1;
    Cy_new[v] = u2;
	
	
	return 0;
	
}

//==============================================================================================================
int Voronoi::which_region(vector<double> Cx, vector<double> Cy,double BHx, double BHy, int nC)
{
	// which Region is a borehole in?
	double small = 100.0;
	double d;
	int id =0;
	
	Voronoi test;
	
	for (int i=0; i<nC; i++)
	{
		// Get distance to Borehole 1
		d = test.get_dist(Cx[i],Cy[i],BHx,BHy,nC);
		//cerr << "Cx[i] " << Cx[i] << "\t d " << d << endl;
		if (d < small) 
		{
			small = d;
			id = i;
		}
	
	}
	id ++;
	//cerr <<"This bh in C: " <<  id << endl;
    return id;	
}

//==============================================================================================================
int Voronoi::any_data(vector<int>& GST_id_new,vector<int>& GST_id, vector<GST>& UKnew, vector<GST>& UKtnew, vector<GST>& UK, 
                       vector<GST>&  UKt, int& nCnew, int nC, int nBH,
                         vector<double>& Cx_new, vector<double>& Cy_new,vector<double>& Cx,vector<double>& Cy, vector<double>&  beta, vector<double> & beta_new )
{

	        vector<int> GSTsort;
    		GSTsort = GST_id_new;
    		bool logi2;
    		vector<bool> logi;
    		
    		sort(GSTsort.begin(), GSTsort.end());
    		//cerr << "GSTsort" << endl;
    		//for (int i=0; i<nBH; i++)
    		//{
    		//	cerr << GSTsort[i] << "\t" ;
    		//}
    		
    		logi.assign(nCnew, false);
			for (int j=1; j<=nCnew; j++)
			{
				logi[j-1] = false;
				for (int i=0; i<nBH; i++)
				{
				    if (j==GSTsort[i])
					{
						logi[j-1] = true;
						j++;
					}
					else{
						logi[j-1] = false;
					}
				//	cerr << "j " << j << "logi" << logi[j-1] << endl;
				}
			}
			logi2 = true;
    		//cerr << "nCnew " << nCnew << endl;
    		for (int j=0; j<nCnew; j++)
    		{
    		//	cerr <<"logi "<< logi[j] << endl;
    			if (logi[j] == false)
    			    {
    			        logi2 = false;
    			    }
    		}
    		//	cerr << "logi2 " << logi2 << endl;
    			if (logi2 ==false)
    			{
	    				cerr << "->Change to region with no data "<< endl;
	    				
	    				nCnew = nC;
			        	UKnew.resize(nC);
			        	UKtnew.resize(nC);
			        	for (int i=0; i<nC; i++)
			        	{
			        		UKnew[i].set_history(UK[i].get_history(),UK[i].get_size());
			        		UKtnew[i].set_history(UKt[i].get_history(),UKt[i].get_size());
			        	}
			        	
			        	
			        	Cx_new.resize(nC);
			        	Cy_new.resize(nC);
			        	
			        	Cx_new = Cx;
			        	Cy_new = Cy;
			        	
			        	GST_id_new = GST_id;
						beta_new = beta;

    			}
    		
	return 0;
}
