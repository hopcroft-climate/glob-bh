#include "baseRNG.h"
//#include "math_functions.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Constructors
BaseRNG::BaseRNG(void){}
BaseRNG::~BaseRNG(void){}
int BaseRNG::_SEED;
void BaseRNG::set_seed(int i) {_SEED =i;}



double BaseRNG::uniform()
//Uniform random generator taken from RAN3 (NR)
{
	static int inext,inextp;
	static int iff=0;
	const int MBIG=1000000000,MSEED=161803398,MZ=0;
	const double FAC=(1.0/MBIG);
	static int ma[56];
	int i,ii,k,mj,mk;

	if (_SEED < 0 || iff == 0) {
		iff=1;
		mj=labs(MSEED-labs(_SEED));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < int(MZ)) mk += MBIG;
			mj=ma[ii];
		}
		for (k=0;k<4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < int(MZ)) ma[i] += MBIG;
			}
			inext=0;
			inextp=31;
			_SEED=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < int(MZ)) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}




void BaseRNG::set_params(double *_params) 
{
	RNGparams = _params;
}



int BaseRNG::get_seed() {return _SEED;}


/*------------------------------------
Gaussian RNG Class 
--------------------------------------*/
GaussianRNG::~GaussianRNG(void)
{
  delete [] RNGparams;
}

GaussianRNG::GaussianRNG(void)
{
	RNGparams = new double[2];
	RNGparams[0] = 0;
	RNGparams[1] =1.0;
}
void GaussianRNG::set_params(double _params[2])
{
    RNGparams[0] = _params[0];
    RNGparams[1] = _params[1];
}
void GaussianRNG::get_params(void)
{
	std::cerr << "Mean: " << RNGparams[0] <<'\t'
		<< "Std Dev: " << RNGparams[1] << '\n';
}
double GaussianRNG::generate(void)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (_SEED < 0) iset=0;
	if (iset == 0) 
	{
		do 
		{
			v1=2.0*uniform()-1.0;
			v2=2.0*uniform()-1.0;
			rsq=v1*v1+v2*v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac*RNGparams[1] + RNGparams[0];
	} 
	else 
	{
		iset=0;
		return gset*RNGparams[1] + RNGparams[0];
	}

}



//Truncated Gaussian definition


TruncatedGaussian::TruncatedGaussian(void){}
TruncatedGaussian::~TruncatedGaussian(void){}

double TruncatedGaussian::generate(void)
{
  double output = -4;
  while (output < -3 | output > 3)
  {
    output = GaussianRNG::generate();
  }
  return output;
}




//------------------------------------
//   Uniform RNG Class 
//--------------------------------------
UniformRNG::~UniformRNG(void){}
UniformRNG::UniformRNG(void){}
double UniformRNG::generate(void)
{
	return uniform();
}
