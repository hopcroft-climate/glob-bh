#ifndef BOREHOLEINFO_H_
#define BOREHOLEINFO_H_

#endif /*BOREHOLEINFO_H_*/
 
 
 struct Borehole_info{
	 	int length_T;				//length temperature vector
	 	int offset;			//number offset elements (always 5)
	 	double upperLayer;		//offset to surface from 1st useable measurement (divided by 5)
	 	double el_off;			// elevation above sea level * lapse rate
	 	int length_C;			// length of cond vector
	 	int years;				// years offset from most recent log
	 	vector<double> Data;	// Temperature measurement
	 	vector<double> DataZ;	// depths
	 	string name_T;			//name of data file
	 	string name_z;			//name of z file
	 	string name_C;			//name of cond file
	 	double bx;				// x co-ordinate of borehole
	 	double by;				// y co-ord borehole
	 	double cond;				// borehole mean conductivity for use in inversion  29.7.2009
                vector<double> air_transfer;            // coefficient of heat transfer air to ground
	 };
	 
