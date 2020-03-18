//--------------------------------------- Element class -------------------------------------------
    class Element
    {
          public:
              double Elementlabel();      
              double SetWidth(double dx);     
              double GetWidth();  
              TNT::Array2D<double> Get_m_matrix();
              TNT::Array2D<double> Get_k_matrix(bool surfaceterm);
              double SetPhys(double density, double area, double ht_cap, double Alpha);
              double GetPhys();      
          private:
              double elementWidth;
            
              TNT::Array2D<double> local_m;
              TNT::Array2D<double> local_k;
              double p;
              double alpha;
              double k; 
              double c;   
              double A;        
    };
    double Element::SetWidth(double dx)
    {      elementWidth = dx;   
    		return 0;}

    double Element::GetWidth()
    { //cerr << " Element width: " << elementWidth << "m" << endl; 
    	return 0; 
    }       
    double Element::SetPhys(double density, double ht_cap, double area, double Alpha)
    {      p  = density;
           c = ht_cap;
           A = area;
           alpha = Alpha;
           return 0;
    }
    double Element::GetPhys()
    {      cerr << "density : " << p << endl;
           cerr << "area: " << A << endl;
           cerr << "heat capacity: " << c << endl;
           cerr << "alpha: " << alpha << endl;
           return 0;
    }
    TNT::Array2D<double> Element::Get_m_matrix()
    {      local_m = local_M(p, c,  A);
           //cerr << " local m_matrix " << local_m << endl;
           return local_m;       }      
    TNT::Array2D<double> Element::Get_k_matrix(bool surfaceterm)
    {       local_k = local_K( A,  surfaceterm, alpha);
            //cerr << " local k_matrix " << local_k << endl;
            return local_k;      }
   //--------------------------------------------------------------------------------------------------------------------
   
   
   
