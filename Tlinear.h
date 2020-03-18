
//====================================================================================================================
//====================================================================================================================
//====================================================================================================================
//====================================================================================================================
double round_new(double a);
Array1D<double> TFe (Array1D<double> Tsurface, Array1D<double>T, int points, int elements, double HtF, double Teq,double dx,
  double  UpperLayer, int CondLength, int  BHnumber, int BHdataLength, double years, vector<double> z, double conductivity_value,
   double years_offset)
{
	
	
	//=========================
	
	

	 
	 //cerr << "BHnumber " << BHnumber << endl;
     //     cerr << "BHdataLength " << BHdataLength << endl;
	vector<double> zn;
	  Array1D<double> kappa = Array1D<double> (elements, 0.0);

	  string EName;
	  //-------------------------------------------------------
	 
	  //cerr << "zName" << zName << endl;
	 // cerr << "BHnumber "<< BHnumber << endl;
	  //cerr << "BHdataLength "<<BHdataLength<<  endl;
	  //cerr << "z.size() " << z.size()  << endl;
	  Array1D<double> dz =Array1D<double> (elements, 0.0);

	 
		  zn.resize(elements);
	          //cerr <<"UpperLayer " <<UpperLayer << endl;
		  for (int i=0; i<5; i++)
		  {
			  zn[i] = UpperLayer*double(i);

		  }
		  for (int i=5; i<(BHdataLength+5); i++)
		  {
			  zn[i] = z[i-5];
		  }
                  for (int i=(BHdataLength+5); i<elements; i++)
		  {
			  zn[i] = zn[i-1]+dx;
		  }
		
		  //--------------------------------------------------------------
                 
	   Array1D<double> cond = Array1D<double> (elements, 0.0);
	cond=Array1D<double>(elements,conductivity_value);
	
     Element oneD;
    double a = 1.0;     
    oneD.SetWidth(dx);
    oneD.GetWidth();
    oneD.SetPhys(2000, 1000, 1.0,  a);    // density, ht_cap, area,  Alpha
   
    oneD.Get_k_matrix(1);
    oneD.Get_m_matrix();
    
    
    //****************Physical Parameters ************************** 
    TNT::Array1D<double> layer= TNT::Array1D<double> (elements, dx);        //dx = element spacing ( 0 is at surface)
    //TNT::Array1D<double> z= TNT::Array1D<double> (elements, 0.0);  
    //TNT::Array1D<double> kappa= TNT::Array1D<double> (elements, 1.0);  
    TNT::Array1D<double> depths= TNT::Array1D<double> (elements, 0.0);  
    
    
    kappa = cond.copy();
   
   // for (int i=0; i<elements; i++)
   // {
   //      cerr << 	kappa[i] << endl;
   // }
    
    for (int i=1; i<(elements); i++)
    {
       
	 layer[i-1] =  zn[i] - zn[i-1];
    }
    
    /*ofstream outLayer("Out/Layer.txt");
    outLayer << layer << endl;
    outLayer.close();*/
    
  //   cerr << "HtF, Teq" << HtF << "\t" << Teq << endl;

    //*********************** K *********************************
    Array1D<double> diag_k      = Array1D<double> (elements, 0.0);
    Array1D<double> diag_kplus  = Array1D<double> (elements-1, 0.0);
    Array1D<double> diag_kminus = Array1D<double> (elements-1, 0.0);
    
   
    for (int i=0; i<elements-1; i++)
    {
        diag_k[i]     += oneD.Get_k_matrix(0)[0][0]*kappa[i]/(layer[i]);
        diag_k[i+1]   += oneD.Get_k_matrix(0)[1][1]*kappa[i]/(layer[i]);
        diag_kplus[i]  = oneD.Get_k_matrix(0)[0][1]*kappa[i]/(layer[i]);
        
    }
    diag_k[0] += oneD.Get_k_matrix(1)[0][0];                             // no * kappa[i] !
    
 
    //***************** M **************************************
    Array1D<double> diag_m      = Array1D<double> (elements, 0.0);
    Array1D<double> diag_mplus  = Array1D<double> (elements-1, 0.0);
    Array1D<double> diag_mminus = Array1D<double> (elements-1, 0.0);
    
    for (int i=0; i<elements-1; i++)
    {
        oneD.SetPhys(2000,1000,1.0,a);
        diag_m[i]     += oneD.Get_m_matrix()[0][0]*layer[i];                       // missed multiplication by layer[i] or L: 22.03.05
        diag_m[i+1]   += oneD.Get_m_matrix()[1][1]*layer[i];
        diag_mplus[i]  = oneD.Get_m_matrix()[0][1]*layer[i];
       // diag_mminus[i] = oneD.Get_m_matrix()[1][0]*layer[i];
       depths[i+1] += layer[i] + depths[i] ;
       
    }
    
    
    //*************** f *******************************************
    TNT::Array1D<double> f = Array1D<double> (elements, 0.0);
     
     //%%%%%%%%%% Assign Heatflow value %%%%%%%%%%%%%%%
     
     f[elements -1] = HtF;
     //%%%%%%%%%% Assign Heatflow value %%%%%%%%%%%%%%%
    
     
    //********************* Triadiagonal Solver ******************************
    
    //Solve system for equillibrium distribution using K, f
    T = Array1D<double> (elements, 0.0);
    // initial steady state condition
                        
    T[0] = Teq;
 
    Array1D<double> solution = Array1D<double> (elements, 0.0);
    
    Array1D<double> diag_kpass  = Array1D<double> (elements, 0.0);
    Array1D<double> kpluspass   = Array1D<double> (elements, 0.0);
    Array1D<double> kminuspass  = Array1D<double> (elements, 0.0); 
    
    diag_kpass=diag_k.copy();
    kpluspass=diag_kplus.copy();

    solution = tsolver_bc(diag_kpass, kpluspass, kpluspass, f, elements, T[0]); 
  
  //  ofstream outis("Out/initialsolution.txt");
  //  for (int i=0; i<elements; i++)
  //  {
  //      outis << solution[i] << endl;
        
  //  }
 //   outis.close();
    T = solution;
    //**************** f (heating) -vector *****************
    f = Array1D<double> (elements, 0.0);  // reset f
    f[elements-1] = HtF;
    
    //___________________________________ Crank Nicholson forward and reverse _______________________________________
    
    Array1D<double> diag_A  = Array1D<double> (elements,   0.0);
    Array1D<double> A_plus  = Array1D<double> (elements-1, 0.0);         
    Array1D<double> B       = Array1D< double>(elements,   0.0);
    Array1D<double> diag_b  = Array1D<double> (elements,   0.0); 
    Array1D<double> b_plus  = Array1D<double> (elements-1, 0.0);
    
    
    
          
    //-------------------------------------------------------------------- 
    double offset_decades=0;
    offset_decades=double(year_today-years_offset)/(dt/(24.0*3600.0*365.25));     // 'decades' are 6 years (not 10) since forward model timestep =6yrs. 19.8.2009
    offset_decades=round_new(offset_decades);
   // cerr << "offset_decades " << years_offset << "\t" << offset_decades << endl;
    int iterations_calc=int(iterations -offset_decades);
    // time - step the equillibrium distribution forward in time.               
    for (int steps = 0; steps <(iterations_calc); steps++)        
    {      
              
              diag_A = Array1D<double> (elements,   0.0);
              A_plus = Array1D<double> (elements-1, 0.0);
              diag_b = Array1D<double> (elements,   0.0); 
              b_plus = Array1D<double> (elements-1, 0.0);
              B =      Array1D< double>(elements,   0.0);
              // set up A and B matrices (see Lewis page p101,103,104,95.            
              for (int i =0; i<elements; i++)
              {  
                   diag_A[i] = diag_m[i] + e*dt*(diag_k[i]); 
                   diag_b[i] = diag_m[i] -(1-e)*dt*(diag_k[i]);
              } 
              for ( int i = 0; i<(elements-1); i++)
              {
                   A_plus[i] = diag_mplus[i] + e*dt*(diag_kplus[i]);
                   b_plus[i] = diag_mplus[i] -(1-e)*dt*diag_kplus[i];
                  
              }   
              
              if (steps == 0) f[0] = a*Tsurface[0];
              else 
              {
                 f[0] = e*a*(Tsurface[steps]) + (1-e)*a* (Tsurface[steps-1]);
              }
              
              for (int i=1; i<(elements-1); i++)
              {
                  B[i] = b_plus[i-1] *T[i-1] +diag_b[i]*T[i] + b_plus[i] * T[i+1];
              }
              
              B[0] = diag_b[0]*T[0] + b_plus[0]*T[1];
              B[elements-1] = b_plus[elements-2]*T[elements-2] + diag_b[elements-1]*T[elements-1];
                 
              for (int i=0; i<elements; i++)
              {
                  B[i] += (dt*f[i]);  
              }
             
              T = tsolver_NR(diag_A, A_plus, A_plus, B, elements);
              
              //outTT << T[elements-1] << endl;
        } 
        return T;
}
double round_new(double a) 
{
return int(a + 0.500);
}


