#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>


/* For boundary condition handling */
#define PERIODIC        1
#define SYMMETRIC       2

/* For the type of wavelet decomposition */
#define WAVELET     1   /* The standard decomposition */
#define STATION     2   /* The stationary decomposition */

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

/*
 * ACCESSC handles negative accesses, as well as those that exceed the number
 * of elements
 */

#define ACCESS(image, size, i, j)       *(image + (i)*(size) + (j))
#define ACCESSC(c, firstC, lengthC, ix, bc) *(c+reflect(((ix)-(firstC)),(lengthC),(bc)))
#define ACCESSD(l, i)   *(Data + (*LengthData*(l)) + (i))
#define POINTD(l,i) (Data + (*LengthData*(l)) + (i))
#define POINTC(l,i) (Carray +(*LengthData*(l)) + (i))




/*
 * comwr:  Do 1D complex wavelet reconstruction
 */

void comwr(double *CR, double *CI, int *LengthC,
	   double *DR, double *DI, int *LengthD,
	double *HR, double *HI, double *GR, double *GI, int *LengthH,
	int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
double *CR::             Input data, and the subsequent smoothed data 
double *CI::             Input data, and the subsequent smoothed data 
int *LengthC::           Length of C array                            
double *DR::             The wavelet coefficients                     
double *DI::             The wavelet coefficients                     
int *LengthD::           Length of D array                            
double *HR::             The smoothing filter H                       
double *HI::             The smoothing filter H                       
double *GR::             The highpass filter H                       
double *GI::             The highpass filter H                       
int *LengthH::           Length of smoothing filter                   
int *levels::            The number of levels in this decomposition   
int *firstC::            The first possible C coef at a given level   
int *lastC::             The last possible C coef at a given level    
int *offsetC::           Offset from C[0] for certain level's coeffs  
int *firstD::            The first possible D coef at a given level   
int *lastD::             The last possible D coef at a given level    
int *offsetD::           Offset from D[0] for certain level's coeffs  
int *type::       	 The type of wavelet decomposition        
int *bc::        	 Method of boundary correction        
int *error::             Error code                                   
 *---------------------*/

{
register int next_level, at_level;
register int verbose;   /* Printing messages, passed in error       */

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc);

if (*error == 1)
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
        break;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
        break;
    }

if (verbose) Rprintf("Building level: ");

*error = 0;

for(next_level = 1; next_level <= *levels; ++next_level)    {

    
    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level - 1; 

    comcbr( (CR+*(offsetC+at_level)),
        (CI+*(offsetC+at_level)),
        (int)(*(lastC+at_level) - *(firstC+at_level) + 1),
        (int)(*(firstC+at_level)),
        (int)(*(lastC+at_level)),
        (DR+*(offsetD+at_level)),
        (DI+*(offsetD+at_level)),
        (int)(*(lastD+at_level) - *(firstD+at_level) + 1),
        (int)(*(firstD+at_level)),
        (int)(*(lastD+at_level)),
        HR, HI, GR, GI,
        (int)*LengthH,
        (CR+*(offsetC+next_level)),
        (CI+*(offsetC+next_level)),
        (int)(*(lastC+next_level) - *(firstC+next_level)+1),
                (int)(*(firstC+next_level)),
                (int)(*(lastC+next_level)),
        (int)(*type),
        (int)(*bc) );
    }
if (verbose)
    Rprintf("\n");

return;

}
/*
 * COMCBR: Does the reconstruction convolution
 */

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc)
{
register int n,k;
register int cfactor;
double sumCR, sumCI, sumDR, sumDI;
double a,b,c,d,e,f;
void commul(double a, double b, double c, double d, double *e, double *f);
int reflect(int n, int lengthC, int bc);


switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumCR = 0.0;
    sumCI = 0.0;
    sumDR = 0.0;
    sumDI = 0.0;

    while( cfactor*k <= n ) {

        a = *(HR + n - cfactor*k);
        b = *(HI + n - cfactor*k);

        c = ACCESSC(c_inR, firstCin, LengthCin, k, bc);
        d = ACCESSC(c_inI, firstCin, LengthCin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumCR += e;
        sumCI += f;

        /* Now D part */

        a = *(GR + n - cfactor*k);
        b = *(GI + n - cfactor*k);

        c = ACCESSC(d_inR, firstDin, LengthDin, k, bc);
        d = ACCESSC(d_inI, firstDin, LengthDin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumDR += e;
        sumDI += f;

        ++k;
        }

    sumCR += sumDR;
    sumCI += sumDI;

    ACCESSC(c_outR, firstCout, LengthCout, n, bc) = sumCR;
    ACCESSC(c_outI, firstCout, LengthCout, n, bc) = sumCI;
    }

}
/*
 * Multiplication: (a+ib)(c+id) = ac-bd +i(bc+ad) = e + i f
 */

void commul(double a, double b, double c, double d, double *e, double *f)
{
*e = (a*c - b*d);
*f = (b*c + a*d);
}
