/*
Sanjay R. Kharche.
Part of PM3 codes.
October 29, 2019.

Program to generate input for the UQ calculation.

https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

As soon as you have a few (say more than 10k instances), the LHS is ensured/I assume. LHS is not done explicitly.

November 2020.
I modified the program to construct output for the Frontiers 0D model.
This model takes inputs numbered from 0 to 53. First 4 are definite, next 4 to 53 are random numbers. Control has mean = 1.
*/


/* my standard sundials headers, constants, and macros that will fit Petsc. Petsc cannot handle cvodes, at least not installations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <cvode/cvode.h>                  	/* prototypes for CVODE fcts., consts.  	*/
#include <nvector/nvector_serial.h>       	/* serial N_Vector types, fcts., macros 	*/
#include <cvode/cvode_dense.h>            	/* prototype for CVDense                		*/
#include <sundials/sundials_dense.h>      /* definitions DlsMat DENSE_ELEM        	*/
#include <sundials/sundials_types.h>      /* definition of type realtype          		*/
#include <sundials/sundials_math.h>

#include <byteswap.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */
#define RTOL        RCONST(1.0e-6)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-3)        /* scalar absolute tolerance components */


/**************************************************Random number generator.****************************************************************************/
/************************************************************************************************************************************************************************/
/*
Notes on the random number generator:
generates a random number on [0,1]-real-interval
genrand_real1() gives you a double between 0 and 1. This is probably what you will use.
*/
/* For Random Numbers */
/* Period parameters */
#define N 				624
#define M 				397
#define MATRIX_A 		0x9908b0dfUL   	/* constant vector a 		*/
#define UPPER_MASK 	0x80000000UL 	/* most significant w-r bits 	*/
#define LOWER_MASK 	0x7fffffffUL 		/* least significant r bits 	*/

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

int solution_broke = 0, somerandomnumber = 2333788678;

/*
The random number generator.
*/

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/**************************** end of random number generator.***********************/

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main (int argc, char *argv[])
{
	  FILE		*output, *random;

	  int fileNum, parNum, numIt;

		// seed the random number generator with system time.
		/* Initialising the random number generator is important. */
		random = fopen("/dev/urandom","rb"); // you need this location.
		fread(&somerandomnumber, sizeof(unsigned int), 1, random); // somerandomnumber has to be global for some reason.
		fclose(random);
		init_genrand(somerandomnumber);

 double U1, U2;
  double X1, X2, Z1, temperature;
  int dialysisoffOn;
  // mu is the mean, sigma is the standard deviation. sigma is got from coeffVariation.
  double mu[55], sigma[55], coeffVariation[55]; // mu is the mean, sigma is the standard deviation.

parNum = 55;

	for(numIt = 0; numIt < parNum; numIt++){
		mu[numIt]	 		= 1.0; // for randoms centered around 1.0.
		coeffVariation[numIt]	= 0.25; // see Richard Clayton 2019/2018 papers, he also did something 2020 August with Peter Coveney.
		sigma[numIt] 			= mu[numIt] * coeffVariation[numIt];
	}

/*************************************************************************************/
// manipulate any special mu's and coeffVariation's here.
double alpha = 0.7; // 70% stenosis equation
int renalFailureOnOff;
if(atoi(argv[4]) ==0){
  renalFailureOnOff == atoi(argv[4]); // temp placeholder.
}
else if(atoi(argv[4]) ==1){
  mu[16] = mu[17] = mu[18] = mu[19] = mu[20] = mu[21] = 1.5; // Roughly multiplying by 4
  mu[33] = mu[34] = mu[35] = mu[36] = mu[37] = mu[38] = 0.5; // Roughly multiplying by 0.35
}else if (atoi(argv[4]) ==2){
  mu[16] = mu[18] = mu[20]  = 1.5; // Roughly multiplying by 4
  mu[33] = mu[35] = mu[37]  = 0.5; // Roughly multiplying by 0.35
}else if (atoi(argv[4]) ==3){
  mu[16] = mu[18] = mu[20] = mu[22] = mu[24] = mu[26] = 1.5; // Roughly multiplying by 4
  mu[33] = mu[35] = mu[37] = mu[39] = mu[41] = mu[43] = 0.5; // Roughly multiplying by 0.35
}else{
  printf("Homogenous vs heterogenous input not properly provided. Needs to be 1 or 2 or 3. Got %d\n", atoi(argv[4]));
}

/*************************************************************************************/

output = fopen("job_file","w");

temperature = atof(argv[2]);
dialysisoffOn = atof(argv[3]);

for(fileNum=0; fileNum < 5000; fileNum++){
// for(fileNum=0; fileNum < 12; fileNum++){
	// name of binary.
	fprintf(output, "./hdt %d\t%f\t%d\t", fileNum, temperature, dialysisoffOn);

	for(numIt = 4; numIt < parNum-1; numIt++){

	if(fileNum==0) Z1 = 1.0; // case 0 is always the control.
	else{
	do{
		U1 	= genrand_real1()	;	U2 	= genrand_real1()	;
		X1 	= sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
		X2 	= sqrt(-2.0 * log(U1)) * sin(2.0 * M_PI * U2);
		if(genrand_real1()>0.5)
		Z1 	= (mu[numIt] + sigma[numIt] * X1); // move to mu mean with sigma standard deviation.
		else
		Z1 	= (mu[numIt] + sigma[numIt] * X2); // making use of both X1 and X2.
	}while(Z1<=0.1); // a small lower bound keeps your population sensible.
	}
	fprintf(output, "%lf ", Z1 ); // this needs to be between 0 and some small number more than mu.
	}

	fprintf(output, "%d ", atoi(argv[1]) ); // this is argv[54] that controls how much output is given in the line diagrams.
  // fprintf(output, "%d ", atoi(argv[4]) ); // this is argv[54] that controls how much output is given in the line diagrams.

	fprintf(output, "\n");
} // end of i loop.

fclose(output);

	return 0;
}
