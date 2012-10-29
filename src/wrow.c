#include <R.h>
#include <stdio.h>
#include <stdlib.h>

/* Error condition      */
#define ERROR           (-1)
#define OK              (0)
 
/* For boundary condition handling */
#define PERIODIC        1
#define SYMMETRIC       2
 
/* For the type of wavelet decomposition */
#define WAVELET         1       /* The standard decomposition */
#define STATION         2       /* The stationary decomposition */
/* Threshold types */
#define HARD    1
#define SOFT    2
 
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
 * waverecons:  Do 1D wavelet reconstruction
 */

/* MAN: changed type of variables from long to int to remove segmentation fault
 */
void wavereconsow(C, D, H, LengthH, levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD, type, bc, error)
double *C;              /* Input data, and the subsequent smoothed data */
double *D;              /* The wavelet coefficients                     */
double *H;              /* The smoothing filter H                       */
int *LengthH;          /* Length of smoothing filter                   */
int *levels;           /* The number of levels in this decomposition   */
int *firstC;           /* The first possible C coef at a given level   */
int *lastC;            /* The last possible C coef at a given level    */
int *offsetC;          /* Offset from C[0] for certain level's coeffs  */
int *firstD;           /* The first possible D coef at a given level   */
int *lastD;            /* The last possible D coef at a given level    */
int *offsetD;          /* Offset from D[0] for certain level's coeffs  */
int *type;     /* The type of wavelet decomposition        */
int *bc;       /* Which boundary handling are we doing     */
int *error;            /* Error code                                   */
{
register int next_level, at_level;
register int verbose;   /* Printing messages, passed in error       */

void conbarow();

if (*error == 1l)
    verbose = 1;
else
    verbose = 0;

if (verbose)
   Rprintf("wavereconsow\n");

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose)Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose)Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
       Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose)Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose)Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose)Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
    }

if (verbose)Rprintf("Building level: ");

*error = 0l;

for(next_level = 1; next_level <= *levels; ++next_level)    {

    
    if (verbose)
       Rprintf("%d ", next_level);

    at_level = next_level - 1; 

    conbarow( (C+*(offsetC+at_level)),
        (int)(*(lastC+at_level) - *(firstC+at_level) + 1),
        (int)(*(firstC+at_level)),
        (D+*(offsetD+at_level)),
        (int)(*(lastD+at_level) - *(firstD+at_level) + 1),
        (int)(*(firstD+at_level)),
        H,
        (int)*LengthH,
        (C+*(offsetC+next_level)),
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
 * CONBAR: Does the reconstruction convolution
 */

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

void conbarow(c_in, LengthCin, firstCin,
       d_in, LengthDin, firstDin,
       H, LengthH,
       c_out, LengthCout, firstCout, lastCout, type, bc)
double *c_in;
int LengthCin;
int firstCin;
double *d_in;
int LengthDin;
int firstDin;
double *H;
int LengthH;
double *c_out;
int LengthCout;
int firstCout;      /* This determines summation over n     */
int lastCout;       /* and this does too                */
int type;       /* The type of wavelet reconstruction       */
int bc;
{
register int n,k;
register int cfactor;
double sumC, sumD;

int reflect();

switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;      /* MAN: added to cover all possibilities, but shouldn't happen */ 
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumC = 0.0;

    while( cfactor*k <= n ) {

        sumC += *(H + n - cfactor*k)*ACCESSC(c_in, firstCin, LengthCin,
                    k, bc);

        ++k;
        }

    /* Now do D part */

    k = CEIL(n-1);

    sumD = 0.0;

    while( cfactor*k <= (LengthH +n -2) )   {

        sumD += *(H+1+cfactor*k-n) * ACCESSC(d_in, firstDin, LengthDin,
                    k, bc);

        ++k;

        }

    if (n & 1)      /* n odd */
        sumC -= sumD;
    else
        sumC += sumD;

    ACCESSC(c_out, firstCout, LengthCout, n, bc) += sumC;
    }

}
/* Works out reflection, as REFLECT, but reports access errors */
int reflect(n, lengthC, bc)
int n;
int lengthC;
int bc;
{

if ((n >= 0) && (n < lengthC))
    return(n);
else if (n<0)   {
    if (bc==PERIODIC)   {
        /*
        n = lengthC+n;
        */
        n = n%lengthC + lengthC*((n%lengthC)!=0);
        if (n < 0)      {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: left info from right\n");
            error("Error occured. Stopping.\n");
            }
        else
            return(n);
        }

    else if (bc==SYMMETRIC) {
        n = -1-n;
        if (n >= lengthC)       {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            error("Error occured. Stopping.\n");
            }
        else
            return(n);
        }

    else    {
        REprintf("reflect: Unknown boundary correction");
        REprintf(" value of %d\n", bc);
            error("Error occured. Stopping.\n");
        }

    }
else    {
    if (bc==PERIODIC)   {
        /*
       Rprintf("periodic extension, was %d (%d) now ",n,lengthC);
        n = n - lengthC; 
        */
        n %= lengthC;
        /*
       Rprintf("%d\n", n);
        */
        if (n >= lengthC)   {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: right info from left\n");
            error("Error occured. Stopping.\n");
            }
        else
            return(n);
        }
    else if (bc==SYMMETRIC) {
        n = 2*lengthC - n - 1;
        if (n<0)        {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            error("Error occured. Stopping.\n");
            }
        else
            return(n);
        }
    else    {
        REprintf("reflect: Unknown boundary correction\n");
            error("Error occured. Stopping.\n");
        }


    }
/* Safety */
REprintf("reflect: SHOULD NOT HAVE REACHED THIS POINT\n");
error("Error occured. Stopping.\n");
return(0); /* for lint only */
}
