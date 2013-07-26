#include "ComplexArithmetic.h"
#include <math.h>
	
// c = a + b
int CAdd(struct Complex a, struct Complex b, struct Complex *c)
{
	c->re=a.re+b.re;
	c->im=a.im+b.im;
	return 0;
}

// c = a * b
int CMult(struct Complex a, struct Complex b, struct Complex *c)
{
	struct Complex t; // use temporary so a/b and c may be the same
	t.re=a.re*b.re-a.im*b.im;
	t.im=a.re*b.im+b.re*a.im;
	*c=t;
	return 0;
}

// e = exp{i * expArg}
int CExpi(double expArg, struct Complex *e)
{
	e->re=cos(expArg);
	e->im=sin(expArg);
	return 0;
}

// e = exp{expArg}
// Reduces to CExpi(-i*expArg) if expArg is pure imaginary
int CExp(struct Complex expArg, struct Complex *e)
{
	double t=exp(expArg.re);
	e->re=t*cos(expArg.im);
	e->im=t*sin(expArg.im);
	return 0;
}
