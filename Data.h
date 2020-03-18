istream& read_data(istream& in, vector<double>& d);
double round (double a);
double zhere ;
double anoisenew=0.0;
GaussianRNG ccG;
int getdata(vector<Borehole_info>& EuropeSet, Array2D<char>& Labels,Array1D<double>& Teq,Array1D<double>& HtF,int nBH)
{
	
	
	

	/*ifstream inBH("DATA/Data.txt");
	inBH >> nBH; 
	
	
	    for (int i=0; i<nBH; i++)
	    {
	    	inBH >> EuropeSet[i].length_T ;
	    }*/
	    for (int i=0; i<nBH; i++)
	    {
	    	EuropeSet[i].offset = 5;
	    }
	   
	   /*
	    for (int i=0; i<nBH; i++)
	    {
	    	//inBH >> EuropeSet[i].el_off ;
                 EuropeSet[i].el_off=0.0;
	    }
	    for (int i=0; i<nBH; i++)
	    {
	    	//inBH >> EuropeSet[i].length_C ;
                EuropeSet[i].length_C=0.0;
	    }
	    inBH.close(); */
      
       
      
      //========Get the filenames for the European boreholes==============
	 ifstream inBHn("DATA/BHN.txt");
	 for ( int ii=0; ii<nBH; ii++)
	 {
	 	getline(inBHn, EuropeSet[ii].name_T); 
	 	int end = int(EuropeSet[ii].name_T.length());
	    string temp  = EuropeSet[ii].name_T;
	    
	    EuropeSet[ii].name_T.resize(end);
	    for (int j=0; j<end-1; j++)
	    { 
	    	EuropeSet[ii].name_T[j] = temp[j];
	    }
	 	cerr << EuropeSet[ii].name_T  << endl;
	 }
	 inBHn.close();
	
	//=========Get Data============================================
     
       
    char FileName[11];
    string EName;
    string ext;
    FileName[0] = (char) (68);		//D
	FileName[1] = (char) (65);		//a
	FileName[2] = (char) (84);		//t
	FileName[3] = (char) (65);		//a
	FileName[4] = (char) (47);		// 'forward slash'
	FileName[5] = (char) (84);		//T
	FileName[6] = (char) (90);		//Z
	FileName[7] = (char) (47);		// 'forward slash'
	
	
	// Read the data vectors.
	
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    {
    	
	        EName = "DATA/TZ/";
	    	ext = ".txt";
	    	EName += EuropeSet[BHnumber].name_T ;
	    	//EName += EuropeSet[nBH-1].name_T ;
	    	cerr <<   EuropeSet[BHnumber].name_T   << endl;
                EName +=ext;
	    	//EuropeSet[BHnumber].name_T = EName;
	        
	         
    	
    	        ifstream indata;
	     	indata.open(EName.c_str());
	     	vector<double> x;
	     	read_data(indata,  x);
	     	
	     	EuropeSet[BHnumber].DataZ.clear();
	     	EuropeSet[BHnumber].Data.clear();
		
		
	     	
	     	for (int jj=0; jj<int(x.size()); jj++)
	     	{
	     		anoisenew =ccG.generate()*0.05;
			//cerr << "anoisenew " << anoisenew << endl;
			
			if (jj % 2 == 0)
	     		{ 
	     			EuropeSet[BHnumber].DataZ.push_back( x[jj]);
				zhere = x[jj];
				//cerr << "z " << zhere << endl;
	     		}
	     		else 
			EuropeSet[BHnumber].Data.push_back( x[jj]);
			// Set measurements to linear function of depth (i.e. heat flow but no climate change)
			//cerr << "T "<<  x[jj] << endl;
			//cerr << "synth T "<<  zhere*30.0/1000.0 << endl;
			//EuropeSet[BHnumber].Data.push_back(zhere*30.0/1000.0 + anoisenew);
			
	     	}
			
			EuropeSet[BHnumber].upperLayer = EuropeSet[BHnumber].DataZ[0];
//			for (int ins=0; ins<5; ins++)
//			{
				//EuropeSet[BHnumber].DataZ.insert(EuropeSet[BHnumber].DataZ.begin()+1,                 
                                  //       (0.2*EuropeSetBHnumber].upperLayer*ins));
//			}		
	     	
	     	
	     	EuropeSet[BHnumber].length_T = int(EuropeSet[BHnumber].Data.size());
	     	//cerr << EuropeSet[BHnumber].length_T  << endl;
	     	//for (int ie=0; ie<EuropeSet[BHnumber].length_T; ie++)
	     //	{cerr <<    EuropeSet[BHnumber].DataZ[ie] << "\t" << EuropeSet[BHnumber].Data[ie] << endl;	    	}
	     	
	     
     
    }
   
    cerr << "Upperlayers " ;
    for (int i=0; i<nBH; i++)
	    {
                EuropeSet[i].upperLayer=EuropeSet[i].DataZ[0]/5.0;
                if(EuropeSet[i].upperLayer==0||EuropeSet[i].upperLayer>100.0){
                cerr << i << "Data.h: upperLayer=0" << endl;
                 return 1;
	   }
//                cerr <<  EuropeSet[i].upperLayer << "\t";
	    }
 //   cerr << endl;

    cerr << "Cond " ;
    
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    {
    	        EName = "DATA/COND/";
	    	ext = ".conductivity.txt";
	    	EName += EuropeSet[BHnumber].name_T ;
	    	EName +=ext;
               ifstream incond;                             
               incond.open(EName.c_str());            
               incond >> EuropeSet[BHnumber].cond;
               //EuropeSet[BHnumber].cond=2.0;
	       if(EuropeSet[BHnumber].cond==0||EuropeSet[BHnumber].cond>10.0){
                 cerr << EName << " Data.h: cond=0" << endl;
                 //return 1;
     	       }   
               incond.close();
    }
    
    cerr << "HtF:" ;
    
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    {
    	        EName = "DATA/DTDZ/";
	    	ext = ".dTdz.txt";
	    	EName += EuropeSet[BHnumber].name_T ;
	    	EName +=ext;
	    	//EuropeSet[BHnumber].name_T = EName;                
//                cerr << "EName  "<<EName << endl;                
                ifstream inhtf;           
                inhtf.open(EName.c_str());
                //cerr << "EName "<< EName << endl;
                inhtf >> HtF[BHnumber];
		//HtF[BHnumber]=30.0;
               
                if(HtF[BHnumber]==0||HtF[BHnumber]>100.0){
                // cerr << EName << " Data.h: HtF=0" << endl;
                 //return 1;
     	       } 
        HtF[BHnumber]=EuropeSet[BHnumber].cond*HtF[BHnumber]/1000.0;
        if(HtF[BHnumber]<=0||HtF[BHnumber]>1.0){
                 cerr << EName << " Data.h: HtF=0"<< HtF[BHnumber] << endl;
                 return 1;
     	       } 
                //cerr <<  HtF[BHnumber] << "\t";   
                inhtf.close();  
    }
   
   
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    {
    	        EName = "DATA/TEQ/";
	    	ext = ".teq.txt";
	    	EName += EuropeSet[BHnumber].name_T ;
	    	EName +=ext;
	    	//EuropeSet[BHnumber].name_T = EName;
                //cerr << "EuropeSet[BHnumber].name_T  "<<EName << endl;
                ifstream inteq;
                inteq.open(EName.c_str());
                inteq >> Teq[BHnumber];
		//Teq[BHnumber] =0.0;
		//Teq[BHnumber] = Teq[BHnumber] + 1.0;
                //cerr <<  Teq[BHnumber] << "\t";
                inteq.close();
       if(Teq[BHnumber]<-10||Teq[BHnumber]>35){
          cerr << EName << " Data.h: Teq=0" << Teq[BHnumber] << endl;
                 //return 1;
    }
   
// Convert dTdz to heat flux:
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    { 

  
     	       } 
       //cerr << "HtF " <<HtF[BHnumber] << "\t" ;
}
    //ofstream outyears;
    //outyears.open("DATA/years_profiles.txt");
    for (int BHnumber=0; BHnumber<nBH; BHnumber++)
    {
    	        EName = "DATA/YEAR/";
	    	ext = ".year.txt";
	    	EName += EuropeSet[BHnumber].name_T ;
	    	EName +=ext;
	    	EuropeSet[BHnumber].name_T = EName;
                ifstream inyear;
                inyear.open(EuropeSet[BHnumber].name_T.c_str());
                inyear >> EuropeSet[BHnumber].years;
                if(EuropeSet[BHnumber].years<=0){
                    cerr << EName << " Data.h: year=0" << endl;
                 return 1;
     	       } 
                //EuropeSet[BHnumber].years = 2000;
		cerr <<  EuropeSet[BHnumber].years << "\t";
		
               // outyears <<  EuropeSet[BHnumber].years <<endl;
     }
     //outyears.close();
// cerr << endl;
    // Define BH locations over a square grid 1.0 X 1.0---
	
	ifstream inBHloc("DATA/BH_locats.txt");
	for (int b=0; b<nBH; b++)
	{
		inBHloc >> EuropeSet[b].bx;
                EuropeSet[b].bx=EuropeSet[b].bx/5.0;
                inBHloc >> EuropeSet[b].by;
                EuropeSet[b].by=EuropeSet[b].by/5.0;
                EuropeSet[b].bx=round(EuropeSet[b].bx)+1;
                EuropeSet[b].by=round(EuropeSet[b].by)+1;
                // Spherical globe:
                if(EuropeSet[b].bx>73){EuropeSet[b].bx=1;}
                if(EuropeSet[b].by>37){EuropeSet[b].by=1;}
                 cerr << EuropeSet[b].bx << "\t"  <<EuropeSet[b].by << endl;
		
		// Now convert from long-lat to hadCm3 grid points: 1:96 and 1:73:
                
	}

        cerr << "finish Data.h" << endl;

	
	
      
    return nBH;
}
//-------------------------------------------------------------

istream& read_data(istream& in, vector<double>& d)
{
	if (in){
		d.clear();
		
		double x;
		while( in>> x)
		  d.push_back(x);
		  
		  in.clear();
	}
	return in;
}

double round(double a)
{ 
    if((a-int(a))>=0.5){ a=int(a)+1;}
   else {a=int(a);}
	 return a;
} 
