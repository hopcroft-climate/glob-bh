// parameterisation of forward problems
const int iterations =123;                 // number of time steps

const int elements =480;                    // number of nodes in finite element model

const double a  = 10.00000;                    // surface heat transfer coefficient

const double dt =5.0*24.0*3600.0*365.25;        // time step length

const double dx  =10.0;                       // nodal spacing metres

const double e  =2.00/3.00;                       // adjust for different time-stepping schemes 
                                                   //eg Galerkin for 2/3 and Crank-Nicholson for 1/2.

//Lengths of input vectors
const int BH2l = 150;
//const int BH3l = 1111;
//const int BH4l = 821;
//const int BH5l = 807;
//const int BH7l = 647;

// Number of elements above first Temperature measurement.
const int d_f_meas_2 = 10;
//const int d_f_meas_3 = 12;
//const int d_f_meas_4 = 10;
//const int d_f_meas_5 = 15;
//const int d_f_meas_7 = 15;

// parameters of the inverse problem


const double PstdDev=1.0;
const double noise =0.05;     //Std noise on Temperature data
const double alpha_H = 2e-4;			//for updating the HeatFlows (option 6);
const double alpha_Teq = 0.15;
const int SSit = 250000;   //sub-sampler iterations
const int burn =  100000;
const double alpha = 5.0E-5;   // For updating the birth/death moves
      
const double year_today=2015;
