

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>         

#include "../TNT/tnt_array2d.h"
#include "../TNT/tnt_array2d_utils.h"
#include "../TNT/tnt_array1d.h"
#include "../TNT/tnt_array1d_utils.h"

using namespace TNT;   
using namespace std;

int main () 
{ 
     int nBH =1012;
     Array2D<double> bh_locat =Array2D<double> (nBH,2);
     ifstream inpos;
     inpos.open("../DATA/BH_locats.txt");
     for(int i=0; i<nBH; i++)
     {
        inpos >> bh_locat[i][0] >> bh_locat[i][1];
         if( bh_locat[i][0] < 0.0)
        { 
          bh_locat[i][0]+=360.0;
        }
        cerr << bh_locat[i][0]<< "\t" << bh_locat[i][1] << endl;
     }
     inpos.close();
    
    
    // -----------------------------------------
    // Set up the grid of points (converted to a 1x1 square on output)
    // -----------------------------------------
     int mg=73;   //change to 5 deg   longitude
     int jgg=37;   //change to 5 deg    latitude
    
      double dlon = 5.0;
      double dlat = 5.0;
    // -----------------------------------------
    
    
     Array1D<double> longitude =Array1D<double> (mg,0.0);
     Array1D<double> latitude =Array1D<double> (jgg,0.0);
     for (int i=0; i<mg; i++)
     {
       longitude[i]=(i)*dlon;
       cerr << longitude[i] << endl;
     }
     for (int j=0; j<jgg; j++)
     {
       latitude[j]=90.0-(j)*dlat;
       cerr << latitude[j] << endl;
     }
     int ic=0;
     Array2D<bool> grid = Array2D<bool> (mg,jgg,false);
     for (int b=0; b<nBH; b++)
     {
      for (int i=0; i<mg-1; i++)
       {
        for (int j=0; j<jgg-1; j++)
        {
          if(bh_locat[b][0]>=longitude[i] && bh_locat[b][0]<longitude[i+1])
          { if(bh_locat[b][1]<latitude[j] && bh_locat[b][1]>=latitude[j+1])
           { 
              grid[i][j] =true;
              
            }
          }
         }
        }
       }
	ofstream outgstid;
      
       ofstream outgr;
       ofstream outbhl;
	 Array2D<int> centres =Array2D<int> (173,2,0.0);
       ic=0;
       outgr.open("colocbh.txt");
       outbhl.open("bhlocats.txt");
       for (int i=0; i<mg; i++)
       {
        for (int j=0; j<jgg; j++)
        {
         if(grid[i][j]==true)
         {
            if(i+1<10 && j+1<10){
            outgr <<"0"<< i+1 << "_0" << j+1 << endl;}
            if(i+1<10 && j+1>=10){
            outgr <<"0"<< i+1 << "_" << j+1 << endl;}
            if(i+1>=10 && j+1<10){
            outgr << i+1 << "_0" << j+1 << endl;}
            if(i+1>=10 && j+1>=10){
            outgr <<i+1 << "_" << j+1 << endl;}
            outbhl << double(i+1)/double(mg) <<"\t" << double(j+1)/double(jgg) <<endl;
            centres[ic][0] = i+1;
            centres[ic][1] = j+1;
	    ic=ic+1;
          } 
         }
        }
        outgr.close();
        outbhl.close();
        
       ofstream outbh_pos;
       outbh_pos.open("bh_pos.txt");
      
	 for (int b=0; b<nBH; b++)
     {
      for (int i=0; i<mg-1; i++)
       {
        for (int j=0; j<jgg-1; j++)
        {
          
          if(bh_locat[b][0]>=longitude[i] && bh_locat[b][0]<longitude[i+1])
          { if(bh_locat[b][1]<latitude[j] && bh_locat[b][1]>=latitude[j+1])
           { 
             outbh_pos << double(i+1)/double(mg) <<"\t" << double(j+1)/double(jgg) <<endl;
	      
             
            }
            
          }
         }
        }
       }
        outbh_pos.close();

      
 Array1D<int> gst_id = Array1D<int> (1012,-1);
       int il;
	int jl; 
        for (int b=0; b<nBH; b++)
     {

        
for (int ic=0; ic<173; ic++)
{       
	il=centres[ic][0]-1;
	jl=centres[ic][1]-1;
	
		if(bh_locat[b][0]>=longitude[il] && bh_locat[b][0]<longitude[il+1])
          { if(bh_locat[b][1]<latitude[jl] && bh_locat[b][1]>=latitude[jl+1])
{
	gst_id[b] = ic+1;
	cerr << "bh site " << b << endl;
	cerr << "gst_id[b] "<< gst_id[b] << endl;
	cerr << bh_locat[b][0] << "\t" << bh_locat[b][1] << endl;

}
}
}
      }
outgstid.open("gstid.txt");
      
	for (int b=0; b<nBH; b++)
{
cerr << gst_id[b] << endl;
outgstid << gst_id[b] << endl;
}
outgstid.close();

     return 0;
}

