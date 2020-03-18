

#ifndef _BASERNG_H_INCLUDED_
#define _BASERNG_H_INCLUDED_

class BaseRNG
{
protected:
	static int _SEED;
	double *RNGparams;
public:
	BaseRNG(void);
	virtual ~BaseRNG(void);
	void set_seed(int i);
	double uniform();
	virtual double generate(void)=0;
	virtual void set_params(double *params);
	int get_seed();
};


class GaussianRNG: public BaseRNG
{
public:
	GaussianRNG(void);
	virtual ~GaussianRNG(void);
	double generate(void);
	void get_params(void);
	virtual void set_params(double params[2]);
};


class UniformRNG: public BaseRNG
{
public:
	UniformRNG(void);
	virtual ~UniformRNG(void);
	double generate(void);
};


class TruncatedGaussian : public GaussianRNG
{
  public:
    TruncatedGaussian(void);
    virtual ~TruncatedGaussian(void);
    double generate(void);
};


#endif
