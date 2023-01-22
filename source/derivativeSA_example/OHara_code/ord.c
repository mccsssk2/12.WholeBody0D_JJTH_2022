/*
THE BEHAVIOUR HAS CHANGED: PAR 85 (ADP) NOW DEPENDS ON PAR 86 (ATP)

This program now has perturbation of:
1. All ical
2. all INaK
3. All electrolytes
4. Volume
5. Capacitance
6. Temperature.

ORD sundials code.
December 2011.
19/12/2013
Prep this program for the ICaL SA study. Prepping and 2 times running done.
So run the FSA as of now.
1. Make sure ALL ICaL, electrolytes, and conductances are included.
2. Run on Volkenstein.

This is the generalised code for everything - dynamical SA (mutiple params,
MI, variability, Sarkar Sobie matrix).

How to run this program:
1. Remove run_number bits, make all pars switched on to run for FSA
2. Reduce pars to 1 and run for LHS

CINC 2014 EV requirements:

*/

// my standard header and #defines, as of 22 Dec. 2013
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */
#define RTOL  	    RCONST(1.0e-12)   /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-6)   /* scalar absolute tolerance components */
#define MAXSTEPS    500000
#define ZERO        RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */

#define num_beats 150 
#define NS        41       // number of sensitivity states: same as NEQ.
#define NEQ       41       /* number of equations ORD basal is 41, one for charge = \Int_0_dt Istim dt  */
#define DELTAT     0.1     // 0.025
#define amp      -80.0
#define duration   0.5
#define NUMPARAMS  109      // doing nao, ko, cao. there are no other extracellular concentrations. conductances, cm
#define NP         109      // Total p: same as NUMPARAMS
#define pcl        1000.0
#define perturbation 0.80 // 40% perturbation. Just get the numbers.

/* For Random Numbers */
/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

int solution_broke = 0, somerandomnumber = 23334556;

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

/************************************************************************/

typedef struct {
realtype p[NUMPARAMS];
realtype IV; // stimulus
realtype lap_val; // the laplacian, when it comes time

double ICaL, IKr, IK1, IKs, INa, Ito; // outputs
double dvdtf;
double vm_alg;
} *UserData;


/* Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int check_flag(void *flagvalue, char *funcname, int opt);


// these are error weights.
static int ewt(N_Vector y, N_Vector w, void *user_data);

int main(int argc, char *argv[]){

  realtype t, tout;
  int iout, i, stimInt, *plist;
  char *str;
  int NOUT, is, j, var_num;
  realtype *pbar, *sdata; // to get 1D arrays from uS which is a p-D array

 FILE *output, *random, *sa_data;

  booleantype sensi, err_con;
  int sensi_meth;
  int if_perturb, run_number; // 0 if no perturbation, 1 if perturbation. This is an input parameter

  if_perturb = (atoi(argv[2]) - atoi(argv[1])); 

  if(if_perturb<2&&atoi(argv[3])==0) // both, GSA and FSA are 0.
   sensi = TRUE; // use this when doing the sampling for variability/MI/sobie.
  else
   sensi = FALSE;

  sensi_meth = CV_STAGGERED; // CV SIMULTANEOUS, CV STAGGERED, or CV STAGGERED1. These are all int values.
  err_con    = TRUE;

  double sstime = 0.0;

double restingPotential[num_beats];
double     maxPotential[num_beats];
double          dvdtmax[num_beats];
double           tstart[num_beats];
double             tend[num_beats];
double           tend30[num_beats];
double            apd30[num_beats];
double            apd90[num_beats];

double     mincai[num_beats];
double     maxcai[num_beats];
double     minnai[num_beats];
double     maxnai[num_beats];
double     minki[num_beats];
double     maxki[num_beats];

double     mincam[num_beats];
double     maxcam[num_beats];

double     mincajsr[num_beats];
double     maxcajsr[num_beats];

double     mincansr[num_beats];
double     maxcansr[num_beats];

double     mincass[num_beats];
double     maxcass[num_beats];

double caimax_time[num_beats]; // time to max cai
double caimin_time[num_beats]; // time to min cai

double ymin[NEQ][num_beats], ymax[NEQ][num_beats];

double old_v = 0.0;

int flag;

for(run_number = atoi(argv[1]); run_number < atoi(argv[2]); run_number++){

sstime = 0.0;
old_v = 0.0;

  N_Vector y, *yS; // yS is the sensitivity solution with length y and breadth pbar
  void *cvode_mem;
  UserData data; // a structure that will contain the modelling parameters. This is essentially the same as STR

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
  pbar      = NULL;
  plist     = NULL;
  yS        = NULL;

 /* Initialize y. The Y1, Y2,... are defined in the header. */
NOUT = (int)(((num_beats)*pcl)/DELTAT);
/* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  data = (UserData) malloc(sizeof *data); // now it is created.
//---------------------------------------------------------------------------
// State variables
//---------------------------------------------------------------------------
// while we are at it, we can also do SA for initial conditions, especially for Nai, Cai, and RMP
   Ith(y,1) = 0.0;   // CaMKt (millimolar) (in CaMK)
   Ith(y,2) = 0.0;   // d (dimensionless) (in ICaL)
   Ith(y,3) = 1.0;   // fcaf (dimensionless) (in ICaL)
   Ith(y,4) = 1.0;   // fcafp (dimensionless) (in ICaL)
   Ith(y,5) = 1.0;   // fcas (dimensionless) (in ICaL)
   Ith(y,6) = 1.0;   // ff (dimensionless) (in ICaL)
   Ith(y,7) = 1.0;   // ffp (dimensionless) (in ICaL)
   Ith(y,8) = 1.0;   // fs (dimensionless) (in ICaL)
   Ith(y,9) = 1.0;   // jca (dimensionless) (in ICaL)
   Ith(y,10) = 0.0;   // nca (dimensionless) (in ICaL)
   Ith(y,11) = 1.0;   // xk1 (dimensionless) (in IK1)
   Ith(y,12) = 0.0;   // xrf (dimensionless) (in IKr)
   Ith(y,13) = 0.0;   // xrs (dimensionless) (in IKr)
   Ith(y,14) = 0.0;   // xs1 (dimensionless) (in IKs)
   Ith(y,15) = 0.0;   // xs2 (dimensionless) (in IKs)
   Ith(y,16) = 1.0;   // hL (dimensionless) (in INaL)
   Ith(y,17) = 1.0;   // hLp (dimensionless) (in INaL)
   Ith(y,18) = 0.0;   // mL (dimensionless) (in INaL)
   Ith(y,19) = 1.0;   // hf (dimensionless) (in INa)
   Ith(y,20) = 1.0;   // hs (dimensionless) (in INa)
   Ith(y,21) = 1.0;   // hsp (dimensionless) (in INa)
   Ith(y,22) = 1.0;   // j (dimensionless) (in INa)
   Ith(y,23) = 1.0;   // jp (dimensionless) (in INa)
   Ith(y,24) = 0.0;   // m (dimensionless) (in INa)
   Ith(y,25) = 0.0;   // a (dimensionless) (in Ito)
   Ith(y,26) = 0.0;   // ap (dimensionless) (in Ito)
   Ith(y,27) = 1.0;   // iF (dimensionless) (in Ito)
   Ith(y,28) = 1.0;   // iFp (dimensionless) (in Ito)
   Ith(y,29) = 1.0;   // iS (dimensionless) (in Ito)
   Ith(y,30) = 1.0;   // iSp (dimensionless) (in Ito)
   Ith(y,31) = 1.0e-4;   // cai (millimolar) (in intracellular_ions)
   Ith(y,32) = 1.2;   // cajsr (millimolar) (in intracellular_ions)
   Ith(y,33) = 1.2;   // cansr (millimolar) (in intracellular_ions)
   Ith(y,34) = 1.0e-4;   // cass (millimolar) (in intracellular_ions)
   Ith(y,35) = 145.0;   // ki (millimolar) (in intracellular_ions)
   Ith(y,36) = 145.0;   // kss (millimolar) (in intracellular_ions)
   Ith(y,37) = 7.0;   // nai (millimolar) (in intracellular_ions)
   Ith(y,38) = 7.0;   // nass (millimolar) (in intracellular_ions)
   Ith(y,39) = -87.0;   // v (millivolt) (in membrane)
   Ith(y,40) = 0.0;   // Jrelnp (dimensionless) (in ryr)
   Ith(y,41) = 0.0;   // Jrelp (dimensionless) (in ryr)

// for non-zero basal values to be used in the SA
/*
double Y0 = 0.00039026566807932699;
Ith(y,1) = 0.00039026566807932699 ; 
double Y1 = 0.00000000230597871971;
Ith(y,2) = 0.00000000230597871971 ; 
double Y2 = 0.99999999103153447422;
Ith(y,3) = 0.99999999103153447422 ; 
double Y3 = 0.99999999103151115953;
Ith(y,4) = 0.99999999103151115953 ; 
double Y4 = 0.99999999103132797273;
Ith(y,5) = 0.99999999103132797273 ; 
double Y5 = 0.99999999103153403013;
Ith(y,6) = 0.99999999103153403013 ; 
double Y6 = 0.99999999103150982727;
Ith(y,7) = 0.99999999103150982727 ; 
double Y7 = 0.99999999109021420196;
Ith(y,8) = 0.99999999109021420196 ; 
double Y8 = 0.99999999103138348389;
Ith(y,9) = 0.99999999103138348389 ; 
double Y9 = 0.00113655347163258672;
Ith(y,10) = 0.00113655347163258672 ; 
double Y10 = 0.99674267737035215919;
Ith(y,11) = 0.99674267737035215919 ; 
double Y11 = 0.00000793884928804177;
Ith(y,12) = 0.00000793884928804177 ; 
double Y12 = 0.00000789914698405012;
Ith(y,13) = 0.00000789914698405012 ; 
double Y13 = 0.00017396454261471326;
Ith(y,14) = 0.00017396454261471326 ; 
double Y14 = 0.00019142613131950859;
Ith(y,15) = 0.00019142613131950859 ; 
double Y15 = 0.51519148171869688646;
Ith(y,16) = 0.51519148171869688646 ; 
double Y16 = 0.31724251978692835641;
Ith(y,17) = 0.31724251978692835641 ; 
double Y17 = 0.00018601095782504458;
Ith(y,18) = 0.00018601095782504458 ; 
double Y18 = 0.70029572916411186867;
Ith(y,19) = 0.70029572916411186867 ; 
double Y19 = 0.70029566766196471583;
Ith(y,20) = 0.70029566766196471583 ; 
double Y20 = 0.45759562800078662503;
Ith(y,21) = 0.45759562800078662503 ; 
double Y21 = 0.70029536619711063583;
Ith(y,22) = 0.70029536619711063583 ; 
double Y22 = 0.70029519837573639407;
Ith(y,23) = 0.70029519837573639407 ; 
double Y23 = 0.00729750019146782527;
Ith(y,24) = 0.00729750019146782527 ; 
double Y24 = 0.00099682861566740633;
Ith(y,25) = 0.00099682861566740633 ; 
double Y25 = 0.00050791002168451484;
Ith(y,26) = 0.00050791002168451484 ; 
double Y26 = 0.99955914652723587555;
Ith(y,27) = 0.99955914652723587555 ; 
double Y27 = 0.99955914655315847295;
Ith(y,28) = 0.99955914655315847295 ; 
double Y28 = 0.99956024149484123953;
Ith(y,29) = 0.99956024149484123953 ; 
double Y29 = 0.99955959145407469180;
Ith(y,30) = 0.99955959145407469180 ; 
double Y30 = 0.00006874781885822879;
Ith(y,31) = 0.00006874781885822879 ; 
double Y31 = 1.21925426330905573025;
Ith(y,32) = 1.21925426330905573025 ; 
double Y32 = 1.21835017545700541319;
Ith(y,33) = 1.21835017545700541319 ; 
double Y33 = 0.00006752597141811216;
Ith(y,34) = 0.00006752597141811216 ; 
double Y34 = 145.05895282940031165708;
Ith(y,35) = 145.05895282940031165708 ; 
double Y35 = 145.05892942850564963919;
Ith(y,36) = 145.05892942850564963919 ; 
double Y36 = 6.94509983176759693180;
Ith(y,37) = 6.94509983176759693180 ; 
double Y37 = 6.94515885960682766154;
Ith(y,38) = 6.94515885960682766154 ; 
double Y38 = -88.06522773172292772870;
Ith(y,39) = -88.06522773172292772870 ; 
double Y39 = 0.00000006867290946662;
Ith(y,40) = 0.00000006867290946662 ; 
double Y40 = 0.00000008584285561025;
Ith(y,41) = 0.00000008584285561025 ; 
*/
realtype Y0 = 0.01233856388121875973;
Ith(y,1) = 0.01233856388121875973;
realtype Y1 = 0.00000000241195801243;
Ith(y,2) = 0.00000000241195801243;
realtype Y2 = 0.99999999055667077030;
Ith(y,3) = 0.99999999055667077030;
realtype Y3 = 0.99999999055381916246;
Ith(y,4) = 0.99999999055381916246;
realtype Y4 = 0.99981769030135769771;
Ith(y,5) = 0.99981769030135769771;
realtype Y5 = 0.99999999055660626635;
Ith(y,6) = 0.99999999055660626635;
realtype Y6 = 0.99999999055364086065;
Ith(y,7) = 0.99999999055364086065;
realtype Y7 = 0.90979096185495345050;
Ith(y,8) = 0.90979096185495345050;
realtype Y8 = 0.99997646936867445877;
Ith(y,9) = 0.99997646936867445877;
realtype Y9 = 0.00265725224026816228;
Ith(y,10) = 0.00265725224026816228;
realtype Y10 = 0.99679287773489122504;
Ith(y,11) = 0.99679287773489122504;
realtype Y11 = 0.00000822700685507778;
Ith(y,12) = 0.00000822700685507778;
realtype Y12 = 0.45315073159181046281;
Ith(y,13) = 0.45315073159181046281;
realtype Y13 = 0.27209707014694450855;
Ith(y,14) = 0.27209707014694450855;
realtype Y14 = 0.00019560211117187397;
Ith(y,15) = 0.00019560211117187397;
realtype Y15 = 0.49712902301408540273;
Ith(y,16) = 0.49712902301408540273;
realtype Y16 = 0.26649838965841943228;
Ith(y,17) = 0.26649838965841943228;
realtype Y17 = 0.00019284666574059426;
Ith(y,18) = 0.00019284666574059426;
realtype Y18 = 0.69370208879468941987;
Ith(y,19) = 0.69370208879468941987;
realtype Y19 = 0.69368507855741490253;
Ith(y,20) = 0.69368507855741490253;
realtype Y20 = 0.44979450093908501795;
Ith(y,21) = 0.44979450093908501795;
realtype Y21 = 0.69359001815020315806;
Ith(y,22) = 0.69359001815020315806;
realtype Y22 = 0.69352623838091509434;
Ith(y,23) = 0.69352623838091509434;
realtype Y23 = 0.00743827985757666756;
Ith(y,24) = 0.00743827985757666756;
realtype Y24 = 0.00100968330199003122;
Ith(y,25) = 0.00100968330199003122;
realtype Y25 = 0.00051446306448339003;
Ith(y,26) = 0.00051446306448339003;
realtype Y26 = 0.99954418140208534105;
Ith(y,27) = 0.99954418140208534105;
realtype Y27 = 0.99954418891528640234;
Ith(y,28) = 0.99954418891528640234;
realtype Y28 = 0.58960284959288744577;
Ith(y,29) = 0.58960284959288744577;
realtype Y29 = 0.64179781286702364262;
Ith(y,30) = 0.64179781286702364262;
realtype Y30 = 0.00008528593866030478;
Ith(y,31) = 0.00008528593866030478;
realtype Y31 = 1.56130905013931298164;
Ith(y,32) = 1.56130905013931298164;
realtype Y32 = 1.60729787049307093483;
Ith(y,33) = 1.60729787049307093483;
realtype Y33 = 0.00008417753717076495;
Ith(y,34) = 0.00008417753717076495;
realtype Y34 = 143.79;
Ith(y,35) = 143.79;
realtype Y35 = 143.79;
Ith(y,36) = 143.79;
realtype Y36 = 7.12718825862747173971;
Ith(y,37) = 7.23;
realtype Y37 = 7.12727112160550202447;
Ith(y,38) = 7.23;
realtype Y38 = -87.87521686743920668050;
Ith(y,39) = -87.87521686743920668050;
realtype Y39 = 0.00000025146367888839;
Ith(y,40) = 0.00000025146367888839;
realtype Y40 = 0.00000031416695655881;
Ith(y,41) = 0.00000031416695655881;

// Ith(y,42) = 0.0; // the charge given by the Istim

/* Initialise these arrays, the values are used in the measurement conditionals */
int apd_counter = 0;
for(apd_counter = 0; apd_counter < num_beats; apd_counter++){
	 restingPotential[apd_counter] = 10000.0;
	 maxPotential[apd_counter] = -10000.0;
         dvdtmax[apd_counter] = -10000.0;
         tstart[apd_counter] = -10000.0;
         tend[apd_counter] = -10000.0;
         tend30[apd_counter] = -10000.0;
	 apd30[apd_counter] = -100000.0;
         apd90[apd_counter] = -10000.0;
	 mincai[apd_counter] = 100000.0;  maxcai[apd_counter] = -10000.0;
	 minnai[apd_counter] = 100000.0;  maxnai[apd_counter] = -10000.0;
	 minki[apd_counter]  = 100000.0;  maxki[apd_counter]  = -100000.0;

      mincam[apd_counter] = 1000.0;
      maxcam[apd_counter] = -1000.0;

      mincajsr[apd_counter] = 1000.0;
      maxcajsr[apd_counter] = -1000.0;

      mincansr[apd_counter] = 1000.0;
      maxcansr[apd_counter] = -1000.0;

      mincass[apd_counter] = 1000.0;
      maxcass[apd_counter] = -1000.0;

	for(var_num=0;var_num<NEQ; var_num++){
		ymin[var_num][num_beats] =  100000.0;
		ymax[var_num][num_beats] = -100000.0;
	}
}

apd_counter = 0;

    str = malloc(64 * sizeof(char));
/* This file will have AP etc. */
    sprintf(str, "output%08d.dat",run_number);
    output = fopen(str, "w+");
    free(str);

stimInt = (int)(pcl/DELTAT);

/* Initialising the random number generator is important. */
random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);

init_genrand(somerandomnumber);

/*********set up your ical parameters *************************************************************/

// d_inf
data->p[0] = 3.94;  // dinf v1/2 2
data->p[1] = 4.23;  // dinf k FSA = 0? 3
// f_inf
data->p[2] = 19.58; // finf v1/2 4
data->p[3] = 3.696; // 5
data->p[4] = 1e-15; // relevant to Grandi, keeping it the same. 6
data->p[5] = 1e-15; // 7
// taud
data->p[6]  = 0.6;   // 8
data->p[7]  = 0.05;  // low FSA 9
data->p[8]  = 6.0;   // 10
data->p[9]  = 0.09;  // 11 
data->p[10] = 14.0;  // 12
/* Tauf:
   tff = 7.0+1.0/(0.0045*exp(-(Y[38]+20.0)/10.0)+0.0045*exp((Y[38]+20.0)/10.0));
   tfs = 1000.0+1.0/(0.000035*exp(-(Y[38]+5.0)/4.0)+0.000035*exp((Y[38]+5.0)/6.0));
*/
data->p[11] = 7.0;      // 13
data->p[12] = 0.0045;   // 14
data->p[13] = 20.0;     // 15
data->p[14] = 10.0;     // 16
data->p[15] = 0.0045;   // 17
data->p[16] = 20.0;     // 18
data->p[17] = 10.0;     // 19
data->p[18] = 1000.0;   // 20
data->p[19] = 0.000035; // 21
data->p[20] = 5.0;      // 22
data->p[21] = 4.0; 
data->p[22] = 0.000035;
data->p[23] = 5.0;
data->p[24] = 6.0;

// af,fast
data->p[25] = 0.6;

// f,ca,fast and f,ca,slow
/* tfcaf = 7.0+1.0/(0.04*exp(-(Y[38]-4.0)/7.0)+0.04*exp((Y[38]-4.0)/7.0));
   tfcas = 100.0+1.0/(0.00012*exp(-Y[38]/3.0)+0.00012*exp(Y[38]/7.0));
   Afcaf = 0.3+0.6/(1.0+exp((Y[38]-10.0)/10.0));
*/

data->p[26] = 7.0;
data->p[27] = 0.04;
data->p[28] = 4.0;
data->p[29] = 7.0;
data->p[30] = 0.04;
data->p[31] = 4.0;
data->p[32] = 7.0;
data->p[33] = 100.0;
data->p[34] = 0.00012;
data->p[35] = 3.0;
data->p[36] = 0.00012;
data->p[37] = 7.0;
data->p[38] = 0.3;
data->p[39] = 0.6;
data->p[40] = 10.0;
data->p[41] = 10.0;

/*
   tffp = 2.5*tff;     // 42
   tfcafp = 2.5*tfcaf; // 43
*/

data->p[42] = 2.5;
data->p[43] = 2.5;
/*
   tffp = data->[42]*tff;     // 42
   tfcafp = data->[43]*tfcaf; // 43
*/

/*
anca = 1.0/(k2n/km2n+pow(1.0+Kmn/Y[33], 4.0));
dY[9] = anca*k2n-Y[9]*km2n;

anca = 1.0/(k2n/km2n+pow(1.0+data->p[44]/Y[33], 4.0));
dY[9] = anca*k2n-Y[9]*data->p[45];
*/

data->p[44] = 0.002; // kmn
data->p[45] = 1000.0; // k2n

//    PhiCaL = 4.0*vffrt*(Y[33]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
//    PhiCaL = 4.0*vffrt*(data->p[46]*Y[33]*exp(2.0*vfrt)-data->p[47]*cao)/(exp(2.0*vfrt)-1.0);

data->p[46] = 1.0;
data->p[47] = 0.341;

// gamma nai and gamma nai are equal in the basal model. perhaps they are associated with some 0.25 elsewhere.
// PhiCaNa = 1.0*vffrt*(0.75*Y[37]*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
// PhiCaNa = 1.0*vffrt*(data->p[48]*Y[37]*exp(1.0*vfrt)-data->p[49]*nao)/(exp(1.0*vfrt)-1.0);

data->p[48] = 0.75;
data->p[49] = 0.75;

// PhiCaK = 1.0*vffrt*(0.75*Y[35]*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
// PhiCaK = 1.0*vffrt*(data->p[50]*Y[35]*exp(1.0*vfrt)-data->p[51]*ko)/(exp(1.0*vfrt)-1.0);

data->p[50] = 0.75;
data->p[51] = 0.75;

//    fICaLp = 1.0/(1.0+KmCaMK/CaMKa);
// CamKa depends on Y0, which in turn has all the Cai parameters assosciated with it. So there you are.
//    fICaLp = 1.0/(1.0+data->p[52]/CaMKa);

data->p[52] = 0.15;

// PCa_b
data->p[53] = 0.0001;

// proportions
/*
   PCap   = 1.1*PCa;
   PCaNa  = 0.00125*PCa;
   PCaK   = 3.574e-4*PCa;
   PCaNap = 0.00125*PCap;
   PCaKp  = 3.574e-4*PCap;
*/

data->p[54] = 1.1;
data->p[55] = 0.00125;
data->p[56] = 3.574e-4;
data->p[57] = 0.00125;
data->p[58] = 3.574e-4;

// electrolytes
   data->p[59]= 1.8; // cao = 1.8;   // millimolar (in extracellular) larger than most.
   data->p[60] = 5.4; // ko = 5.4;   // millimolar (in extracellular) larger than most
   data->p[61] = 140.0; // nao = 140.0;   // millimolar (in extracellular)

// Conductances up for SA.
/*   PCab   = */ data->p[62] = 2.5e-8;   // milliS_per_microF (in ICab)
/*   GK1_b  = */ data->p[63] = 0.1908;   // milliS_per_microF (in IK1)
/*   GKb_b  = */ data->p[64] = 0.003;   // milliS_per_microF (in IKb)

/*   GKr_b  = */ data->p[65] = 0.046;   // milliS_per_microF (in IKr)
/*   GKs_b  = */ data->p[66] = 0.0034;   // milliS_per_microF (in IKs)
/*   Gncx_b = */ data->p[67] = 0.0008;   // milliS_per_microF (in INaCa_i)
/*   Pnak_b = */ data->p[68] = 30.0;   // milliS_per_microF (in INaK)
// Gnal_b multiplier: control Nai at least.
/*   GNaL_b = */ data->p[69] = 0.5*0.0075;   // milliS_per_microF (in INaL) To reduce drift of sodium in the cell.

/*   PNab  = */ data->p[70] = 3.75e-10;   // milliS_per_microF (in INab)
/*   GNa   = */ data->p[71] = 75.0;   // milliS_per_microF (in INa)
/*   GpCa  = */ data->p[72] = 0.0005;   // milliS_per_microF (in IpCa)
/*   Gto_b = */ data->p[73] = 0.02;   // milliS_per_microF (in Ito)

// structural parameters
/*   cm    = */ data->p[74] = 1.0;   // microF_per_centimeter_squared (in intracellular_ions)
	        data->p[75] = 0.01; // cell length

// INaK parameters.
/*   H= */      data->p[76] = 1.0e-7;   // millimolar (in INaK) the pH of cytosol
/*   Khp= */    data->p[77] = 1.698e-7;   // millimolar (in INaK)
/*   Kki= */    data->p[78] = 0.5;   // per_millisecond (in INaK)
/*   Kko= */    data->p[79] = 0.3582;   // per_millisecond (in INaK)
/*   Kmgatp= */ data->p[80] = 1.698e-7;   // millimolar (in INaK)
/*   Knai0= */  data->p[81] = 9.073;   // millimolar (in INaK)
/*   Knao0= */  data->p[82] = 27.78;   // millimolar (in INaK)
/*   Knap= */   data->p[83] = 224.0;   // millimolar (in INaK)
/*   Kxkur= */  data->p[84] = 292.0;   // millimolar (in INaK)
/*   MgATP= */  data->p[86] = 9.8;   // millimolar (in INaK)

// this is not an independent paramter.
/*   MgADP= */  data->p[85] = 10.0 - data->p[86]; // 0.05;   // millimolar (in INaK)


/*   delta= */  data->p[87] = -0.155;   // millivolt (in INaK)
/*   eP= */     data->p[88] = 4.2;   // dimensionless (in INaK)
/*   k1m= */    data->p[89] = 182.4;   // per_millisecond (in INaK)
/*   k1p= */    data->p[90] = 949.5;   // per_millisecond (in INaK)
/*   k2m= */    data->p[91] = 39.4;   // per_millisecond (in INaK)
/*   k2p= */    data->p[92] = 687.2;   // per_millisecond (in INaK)
/*   k3m= */    data->p[93] = 79300.0;   // per_millisecond (in INaK)
/*   k3p= */    data->p[94] = 1899.0;   // per_millisecond (in INaK)
/*   k4m= */    data->p[95] = 40.0;   // per_millisecond (in INaK)
/*   k4p= */    data->p[96] = 639.0;   // per_millisecond (in INaK)
// iPCa Michalis-Menton
/*   KmCap = */ data->p[97] = 0.0005;   // millimolar (in IpCa)
// SERCA
// the CaM parameters are set in the CaM parameter location.
   data->p[98]  = 0.004375;
   data->p[99]  = 0.00092;
   data->p[100] = 2.75*0.004375;
   data->p[101] = 0.00075;
   data->p[102] = 0.0039375/15.0;

// CaM parameters.
/*   CaMKo  = */ data->p[103] = 0.05;   // dimensionless (in CaMK)
/*   KmCaM  = */ data->p[104] = 0.0015;   // millimolar (in CaMK)
/*   KmCaMK = */ data->p[105] = 0.15;   // millimolar (in CaMK)
/*   aCaMK  = */ data->p[106] = 0.05;   // per_millimolar_per_millisecond (in CaMK)
/*   bCaMK  = */ data->p[107] = 0.00068;   // per_millisecond (in CaMK)

// I like temperature also
/*   T = */ data->p[108] = 310.0;   // kelvin (in physical_constants)

/**************************************************************************************/

//if(argc!=3){
//printf("This program takes 2 inputs, 1. ifperturb (0 or 1) and run number (a non-zero int)\n");
//exit(0);
//}
//else
//{
// input 0 is the binary
// if_perturb    = atoi(argv[1]); // 1st input
// run_number    = atoi(argv[2]); // 2nd input
// printf("%d %d \n",if_perturb,run_number);
//}


if(((atoi(argv[2]) - atoi(argv[1]))>1)&&atoi(argv[3])>0){
	/*
	Now make it all a bit random.
	genrand_real1(); random numbers between 0 and 1 inclusive.

	The things that must be done to the parameters:
	1. Take initial conditions as parameters - these systems have infinite stable stationary points - see Jacquamet
	2. To apply constraints (positive, bounds, no change of sign)
	3. To take the parameters from a Latin Hypercube sample
	*/

	for(i=0;i<NUMPARAMS;i++)
	data->p[i] = data->p[i]*( 1.0 + perturbation*(genrand_real1()-0.5) );

} // end of if_perturb


    // solver
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeInit(cvode_mem, f, 0.0, y);
    /* Error weights. This is used during run time. */
    CVodeWFtolerances(cvode_mem, ewt);
    CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data.
    CVodeSStolerances(cvode_mem, RTOL, ATOL);
    CVDense(cvode_mem, NEQ);
    CVodeSetMaxStep(cvode_mem, DELTAT);
    iout = 0;  tout = DELTAT;
  /* Sensitivity-related settings. Typically, if you have copied everything else correctly, this also just needs pasting in the right place. */
if(sensi==TRUE){
    plist = (int *) malloc(NP * sizeof(int));            // so as of now, this is of size 1
    for(is=0; is<NP; is++) plist[is] = is;               // list of ints from 0 to NP-1. This is the iterator for pbar.
    pbar  = (realtype *) malloc(NP * sizeof(realtype));  // list of SA parameters. In the production case, this could be 1 by 1, gcal, then Ko, ...
    for(is=0; is<NP; is++) pbar[is] = data->p[plist[is]];// this is because pbar cannot be a struct or union, it has to be a 1D array.
    yS = N_VCloneVectorArray_Serial(NP, y);              // so uS is a 2D array, each uS[i] has NS bits in it. This function clones u NS times.
    for(is=0;is<NP;is++) N_VConst(ZERO, yS[is]);         // SA initial conditions, unless they depend on parameters.
    CVodeSensInit1(cvode_mem, NP, sensi_meth, NULL, yS); // The NULL would have fS, if you has an analytical Jacobian and gradient.
    CVodeSensEEtolerances(cvode_mem);                    // set up tolerances
    CVodeSetSensErrCon(cvode_mem, err_con);              // always have err_con switched on.
    CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);  // this gives a central difference based esitmate of the RHS. Using centered gives O(dt^2)
    CVodeSetSensParams(cvode_mem, data->p, pbar, plist); // set up the sensitivity paramters in the model memory
    CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);
} // end of sensi == TRUE

  for(iout = 0; iout <= NOUT; iout++) { // start of time loop.

     // stimulation: Conservative form.
     if((iout%stimInt>0)&&(iout%stimInt<=(int)(duration/DELTAT))) 
		data->IV = amp; 
	else 
		data->IV = -amp*duration/(pcl-duration);


     old_v = Ith(y,39);
     flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1)) break;

     if(sensi==TRUE)
     CVodeGetSens(cvode_mem, &t, yS); // call to CVODEGETSENS assigns the values of uS by querying cvode_mem

/* output */
if(iout%10==0&&apd_counter>(num_beats-2)&&sensi==TRUE){ // I need this just for 1 AP and 1 AP only. do for 2 in test, then 1 for production.

// var by var
for(var_num = 0; var_num < NEQ; var_num++){
    str = malloc(64 * sizeof(char));
    sprintf(str, "ord_sa_var%08d.dat",var_num+1); // matched with Ith(y,)
    sa_data = fopen(str, "a+");
    free(str);
	fprintf(sa_data,"%f %f ",sstime,Ith(y,var_num+1));
	for(j=0;j<NP;j++){
		sdata = NV_DATA_S(yS[j]); // FSA of all variables w.r.t. parameter j
//		if(ymax[var_num][apd_counter]>ymin[var_num][apd_counter])
		if(Ith(y,var_num)>0.000001||Ith(y,var_num)<-0.000001)
		fprintf(sa_data,"%10.10f\t",pbar[j]*sdata[var_num] /* /Ith(y,var_num+1) */ ); /* /(ymax[var_num][apd_counter]-ymin[var_num][apd_counter]) */ // FSA of variable y_{num_var} w.r.t. pj.
		else
		fprintf(sa_data,"%10.10f\t",-1000000.0);

//		else
//		fprintf(sa_data,"%10.10f\t",/* pbar[j]*sdata[var_num]* */ (-1000.0) ); // marked to set it out.
	}
		fprintf(sa_data,"\n");
   fclose(sa_data);
} // end of var_num

// par by par. This can be inside the above NP loop, but for clarity is left independent of it.
for(j=0;j<NP;j++){
    str = malloc(64 * sizeof(char));
    sprintf(str, "ord_sa_par%08d.dat",j);
    sa_data = fopen(str, "a+");
    free(str);
	fprintf(sa_data,"%f %f %f %f ",sstime,Ith(y,39),Ith(y,31),Ith(y,37)); // remove this line afterwards
   sdata = NV_DATA_S(yS[j]); // FSA of all variables w.r.t. parameter j
	for(var_num = 0; var_num < NEQ; var_num++){
//		if(ymax[var_num][apd_counter]>ymin[var_num][apd_counter])
		if(Ith(y,var_num)>0.000001||Ith(y,var_num)<-0.000001)
    		fprintf(sa_data,"%10.10f\t",pbar[j]*sdata[var_num] /* /Ith(y,var_num+1) */ ); /* /(ymax[var_num][apd_counter]-ymin[var_num][apd_counter]) */ // the units are not uniform, so the scale is not uniform.
		else
		fprintf(sa_data,"%10.10f\t",-1000000.0);
//		else
//    		fprintf(sa_data,"%10.10f\t",/* pbar[j]*sdata[var_num]* */ (-1000.0));
	} // end of var_num
	fprintf(sa_data,"\n");
   fclose(sa_data);

} // end of j array


} // end of SA output.


    if(iout%100==0&&apd_counter>(num_beats-10)){
	fprintf(output,"%d %10.10f\t",run_number, sstime);
	for(i=0;i<NEQ;i++) fprintf(output,"%10.10f\t",Ith(y,i+1));
	fprintf(output,"%10.10f\t",data->vm_alg);
	fprintf(output,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t\n",data->ICaL,data->IK1,data->Ito,data->IKs,data->IKr,data->INa);
    }

/* Measurements */
// max and min of all variables. Grandi Ki is constant, but its not a variable any more.
	for(var_num=0;var_num<NEQ; var_num++){
		if( ymin[var_num][apd_counter]>Ith(y,var_num+1) ) ymin[var_num][apd_counter] =  Ith(y,var_num+1);
		if( ymin[var_num][apd_counter]<Ith(y,var_num+1) ) ymax[var_num][apd_counter] =  Ith(y,var_num+1);
	}

if(maxPotential[apd_counter]<Ith(y,39)) maxPotential[apd_counter] = Ith(y,39);                  // max pot
if(maxcai[apd_counter]<Ith(y,31)){ maxcai[apd_counter] = Ith(y,31); caimax_time[apd_counter] = sstime; } // max cai
if(maxnai[apd_counter]<Ith(y,37))  maxnai[apd_counter] = Ith(y,37);                             // max nai
if(maxki[apd_counter]<Ith(y,35))   maxki[apd_counter]  = Ith(y,35);                             // max ki

if(maxcam[apd_counter]<Ith(y,1))     maxcam[apd_counter]    = Ith(y,1);                   // max cam
if(maxcajsr[apd_counter]<Ith(y,32))  maxcajsr[apd_counter]  = Ith(y,32);                  // max cajsr
if(maxcansr[apd_counter]<Ith(y,33))  maxcansr[apd_counter]  = Ith(y,33);                  // max cansr
if(maxcass[apd_counter] <Ith(y,34))  maxcass[apd_counter]   = Ith(y,34);                  // max cass

if(mincai[apd_counter]>Ith(y,31)&& (iout%stimInt>(int)(10.0*duration/DELTAT))   ){ mincai[apd_counter] = Ith(y,31); caimin_time[apd_counter] = sstime; } // mincai
if(minnai[apd_counter]>Ith(y,37)&& (iout%stimInt>(int)(10.0*duration/DELTAT))   ) minnai[apd_counter] = Ith(y,37); // minnai
if(minki[apd_counter]>Ith(y,35) && (iout%stimInt>(int)(10.0*duration/DELTAT))   ) minki[apd_counter]  = Ith(y,35); // minki

if(mincam[apd_counter]>Ith(y,1)   && (iout%stimInt>(int)(10.0*duration/DELTAT))   ) mincam[apd_counter]   = Ith(y,1) ; // mincam
if(mincajsr[apd_counter]>Ith(y,32)&& (iout%stimInt>(int)(10.0*duration/DELTAT))   ) mincajsr[apd_counter] = Ith(y,32); // mincajsr
if(mincansr[apd_counter]>Ith(y,33)&& (iout%stimInt>(int)(10.0*duration/DELTAT))   ) mincansr[apd_counter] = Ith(y,33); // mincansr
if(mincass[apd_counter]>Ith(y,34) && (iout%stimInt>(int)(10.0*duration/DELTAT))   ) mincass[apd_counter]  = Ith(y,34); // mincass

if(dvdtmax[apd_counter]<data->dvdtf){ // dvdtmax
 dvdtmax[apd_counter]=data->dvdtf; tstart[apd_counter] = sstime;
}

if((iout%stimInt>0)&&(iout%stimInt<=1)){ // resting
	restingPotential[apd_counter] = Ith(y,39);
//        printf("%d %f %f \n",apd_counter,Ith(y,39),sstime);
} // end of if

if(old_v>0.3*restingPotential[apd_counter] && Ith(y,39)<=0.3*restingPotential[apd_counter]){
	tend30[apd_counter] = sstime;
}
if(old_v>0.9*restingPotential[apd_counter] && Ith(y,39)<=0.9*restingPotential[apd_counter]){
	tend[apd_counter] = sstime;
	if(apd_counter < num_beats) apd_counter++;
}

      sstime = (double)iout*DELTAT;
      tout += DELTAT;
//    printf("%f\n",sstime);
} // end of time loop.
  fclose(output);

str = malloc(64*sizeof(char));
sprintf(str,"cellData%08d.dat",run_number);
output = fopen(str,"a+");
free(str);

/* Get the biomarkers and put output into a file, for a first run, just take 1 output at the end. */
for(apd_counter = num_beats-3; apd_counter < num_beats-2; apd_counter++){
fprintf(output,"%d ",run_number); // 1
/* Parameters  */
for(i = 0; i < NUMPARAMS; i++) fprintf(output,"%40.40f ",data->p[i]); // 2 to 110
/* Outputs */
apd90[apd_counter] = tend[apd_counter]   - tstart[apd_counter];
apd30[apd_counter] = tend30[apd_counter] - tstart[apd_counter];
//						64       65      66                           67                      68                   69
fprintf(output,"%d %40.40f %40.40f %40.40f %40.40f %40.40f %40.40f %40.40f %40.40f ",apd_counter,pcl,restingPotential[apd_counter],dvdtmax[apd_counter],maxPotential[apd_counter],apd90[apd_counter],apd30[apd_counter],maxcai[apd_counter],mincai[apd_counter]);
fprintf(output,"%40.40f %40.40f %40.40f %40.40f ",maxnai[apd_counter],minnai[apd_counter],maxki[apd_counter],minki[apd_counter]);
fprintf(output,"%40.40f %40.40f %40.40f %40.40f ",maxcam[apd_counter],mincam[apd_counter],maxcajsr[apd_counter],mincajsr[apd_counter]);
fprintf(output,"%40.40f %40.40f %40.40f %40.40f ",maxcansr[apd_counter],mincansr[apd_counter],maxcass[apd_counter],mincass[apd_counter]);
fprintf(output,"%40.40f %40.40f ",caimax_time[apd_counter],caimin_time[apd_counter]); // calculate Cai relaxation from here
fprintf(output,"\n");
} // end of loop

fclose(output);

// output initial conditions.
  N_VDestroy_Serial(y);
  CVodeFree(&cvode_mem);
} // end of run number loop

return 0;
}

/* f routine. Compute function f(t,y) */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
 realtype Y[NEQ];
 realtype dY[NEQ];
 int ii; 

 UserData data;
 data = (UserData) user_data;

  for(ii=0;ii<NEQ;ii++) Y[ii] = Ith(y,ii+1);

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

  double CaMKo;   // dimensionless (in CaMK)
  double KmCaM;   // millimolar (in CaMK)
  double KmCaMK;   // millimolar (in CaMK)
  double aCaMK;   // per_millimolar_per_millisecond (in CaMK)
  double bCaMK;   // per_millisecond (in CaMK)
  double Kmn;   // millimolar (in ICaL)
  double PCa_b;   // dimensionless (in ICaL)
  double k2n;   // per_millisecond (in ICaL)
  double PCab;   // milliS_per_microF (in ICab)
  double GK1_b;   // milliS_per_microF (in IK1)
  double GKb_b;   // milliS_per_microF (in IKb)
  double GKr_b;   // milliS_per_microF (in IKr)
  double GKs_b;   // milliS_per_microF (in IKs)
  double Gncx_b;   // milliS_per_microF (in INaCa_i)
  double KmCaAct;   // millimolar (in INaCa_i)
  double kasymm;   // dimensionless (in INaCa_i)
  double kcaoff;   // per_millisecond (in INaCa_i)
  double kcaon;   // per_millisecond (in INaCa_i)
  double kna1;   // per_millisecond (in INaCa_i)
  double kna2;   // per_millisecond (in INaCa_i)
  double kna3;   // per_millisecond (in INaCa_i)
  double qca;   // dimensionless (in INaCa_i)
  double qna;   // dimensionless (in INaCa_i)
  double wca;   // dimensionless (in INaCa_i)
  double wna;   // dimensionless (in INaCa_i)
  double wnaca;   // dimensionless (in INaCa_i)
  double H;   // millimolar (in INaK)
  double Khp;   // millimolar (in INaK)
  double Kki;   // per_millisecond (in INaK)
  double Kko;   // per_millisecond (in INaK)
  double Kmgatp;   // millimolar (in INaK)
  double Knai0;   // millimolar (in INaK)
  double Knao0;   // millimolar (in INaK)
  double Knap;   // millimolar (in INaK)
  double Kxkur;   // millimolar (in INaK)
  double MgADP;   // millimolar (in INaK)
  double MgATP;   // millimolar (in INaK)
  double Pnak_b;   // milliS_per_microF (in INaK)
  double delta;   // millivolt (in INaK)
  double eP;   // dimensionless (in INaK)
  double k1m;   // per_millisecond (in INaK)
  double k1p;   // per_millisecond (in INaK)
  double k2m;   // per_millisecond (in INaK)
  double k2p;   // per_millisecond (in INaK)
  double k3m;   // per_millisecond (in INaK)
  double k3p;   // per_millisecond (in INaK)
  double k4m;   // per_millisecond (in INaK)
  double k4p;   // per_millisecond (in INaK)
  double GNaL_b;   // milliS_per_microF (in INaL)
  double thL;   // millisecond (in INaL)
  double PNab;   // milliS_per_microF (in INab)
  double Ahf;   // dimensionless (in INa)
  double GNa;   // milliS_per_microF (in INa)
  double hssV1;   // millivolt (in INa)
  double hssV2;   // millivolt (in INa)
  double mssV1;   // millivolt (in INa)
  double mssV2;   // millivolt (in INa)
  double mtD1;   // dimensionless (in INa)
  double mtD2;   // dimensionless (in INa)
  double mtV1;   // millivolt (in INa)
  double mtV2;   // millivolt (in INa)
  double mtV3;   // millivolt (in INa)
  double mtV4;   // millivolt (in INa)
  double GpCa;   // milliS_per_microF (in IpCa)
  double KmCap;   // millimolar (in IpCa)
  double Gto_b;   // milliS_per_microF (in Ito)
  double L;   // centimeter (in cell_geometry)
  double rad;   // centimeter (in cell_geometry)
  double celltype;   // dimensionless (in environment)
  double cao;   // millimolar (in extracellular)
  double ko;   // millimolar (in extracellular)
  double nao;   // millimolar (in extracellular)
  double BSLmax;   // millimolar (in intracellular_ions)
  double BSRmax;   // millimolar (in intracellular_ions)
  double KmBSL;   // millimolar (in intracellular_ions)
  double KmBSR;   // millimolar (in intracellular_ions)
  double cm;   // microF_per_centimeter_squared (in intracellular_ions)
  double cmdnmax_b;   // millimolar (in intracellular_ions)
  double csqnmax;   // millimolar (in intracellular_ions)
  double kmcmdn;   // millimolar (in intracellular_ions)
  double kmcsqn;   // millimolar (in intracellular_ions)
  double kmtrpn;   // millimolar (in intracellular_ions)
  double trpnmax;   // millimolar (in intracellular_ions)
//  double amp;   // microA_per_microF (in membrane)
//  double duration;   // millisecond (in membrane)
  double F;   // coulomb_per_mole (in physical_constants)
  double R;   // joule_per_kilomole_kelvin (in physical_constants)
  double T;   // kelvin (in physical_constants)
  double zca;   // dimensionless (in physical_constants)
  double zk;   // dimensionless (in physical_constants)
  double zna;   // dimensionless (in physical_constants)
  double PKNa;   // dimensionless (in reversal_potentials)
  double bt;   // millisecond (in ryr)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

  double CaMKa;   // millimolar (in CaMK)
  double CaMKb;   // millimolar (in CaMK)
  double Afcaf;   // dimensionless (in ICaL)
  double Afcas;   // dimensionless (in ICaL)
  double Aff;   // dimensionless (in ICaL)
  double Afs;   // dimensionless (in ICaL)
  double ICaK;   // microA_per_microF (in ICaL)
  double ICaL;   // microA_per_microF (in ICaL)
  double ICaNa;   // microA_per_microF (in ICaL)
  double PCa;   // dimensionless (in ICaL)
  double PCaK;   // dimensionless (in ICaL)
  double PCaKp;   // dimensionless (in ICaL)
  double PCaNa;   // dimensionless (in ICaL)
  double PCaNap;   // dimensionless (in ICaL)
  double PCap;   // dimensionless (in ICaL)
  double PhiCaK;   // dimensionless (in ICaL)
  double PhiCaL;   // dimensionless (in ICaL)
  double PhiCaNa;   // dimensionless (in ICaL)
  double anca;   // dimensionless (in ICaL)
  double dss;   // dimensionless (in ICaL)
  double f;   // dimensionless (in ICaL)
  double fICaLp;   // dimensionless (in ICaL)
  double fca;   // dimensionless (in ICaL)
  double fcap;   // dimensionless (in ICaL)
  double fcass;   // dimensionless (in ICaL)
  double fp;   // dimensionless (in ICaL)
  double fss;   // dimensionless (in ICaL)
  double km2n;   // per_millisecond (in ICaL)
  double td;   // millisecond (in ICaL)
  double tfcaf;   // millisecond (in ICaL)
  double tfcafp;   // millisecond (in ICaL)
  double tfcas;   // millisecond (in ICaL)
  double tff;   // millisecond (in ICaL)
  double tffp;   // millisecond (in ICaL)
  double tfs;   // millisecond (in ICaL)
  double tjca;   // millisecond (in ICaL)
  double ICab;   // microA_per_microF (in ICab)
  double GK1;   // milliS_per_microF (in IK1)
  double IK1;   // microA_per_microF (in IK1)
  double rk1;   // millisecond (in IK1)
  double txk1;   // millisecond (in IK1)
  double xk1ss;   // dimensionless (in IK1)
  double GKb;   // milliS_per_microF (in IKb)
  double IKb;   // microA_per_microF (in IKb)
  double xkb;   // dimensionless (in IKb)
  double Axrf;   // dimensionless (in IKr)
  double Axrs;   // dimensionless (in IKr)
  double GKr;   // milliS_per_microF (in IKr)
  double IKr;   // microA_per_microF (in IKr)
  double rkr;   // dimensionless (in IKr)
  double txrf;   // millisecond (in IKr)
  double txrs;   // millisecond (in IKr)
  double xr;   // dimensionless (in IKr)
  double xrss;   // dimensionless (in IKr)
  double GKs;   // milliS_per_microF (in IKs)
  double IKs;   // microA_per_microF (in IKs)
  double KsCa;   // dimensionless (in IKs)
  double txs1;   // millisecond (in IKs)
  double txs2;   // millisecond (in IKs)
  double xs1ss;   // dimensionless (in IKs)
  double xs2ss;   // dimensionless (in IKs)
  double E1_i;   // dimensionless (in INaCa_i)
  double E1_ss;   // dimensionless (in INaCa_i)
  double E2_i;   // dimensionless (in INaCa_i)
  double E2_ss;   // dimensionless (in INaCa_i)
  double E3_i;   // dimensionless (in INaCa_i)
  double E3_ss;   // dimensionless (in INaCa_i)
  double E4_i;   // dimensionless (in INaCa_i)
  double E4_ss;   // dimensionless (in INaCa_i)
  double Gncx;   // milliS_per_microF (in INaCa_i)
  double INaCa_i;   // microA_per_microF (in INaCa_i)
  double INaCa_ss;   // microA_per_microF (in INaCa_i)
  double JncxCa_i;   // millimolar_per_millisecond (in INaCa_i)
  double JncxCa_ss;   // millimolar_per_millisecond (in INaCa_i)
  double JncxNa_i;   // millimolar_per_millisecond (in INaCa_i)
  double JncxNa_ss;   // millimolar_per_millisecond (in INaCa_i)
  double allo_i;   // dimensionless (in INaCa_i)
  double allo_ss;   // dimensionless (in INaCa_i)
  double h10_i;   // dimensionless (in INaCa_i)
  double h10_ss;   // dimensionless (in INaCa_i)
  double h11_i;   // dimensionless (in INaCa_i)
  double h11_ss;   // dimensionless (in INaCa_i)
  double h12_i;   // dimensionless (in INaCa_i)
  double h12_ss;   // dimensionless (in INaCa_i)
  double h1_i;   // dimensionless (in INaCa_i)
  double h1_ss;   // dimensionless (in INaCa_i)
  double h2_i;   // dimensionless (in INaCa_i)
  double h2_ss;   // dimensionless (in INaCa_i)
  double h3_i;   // dimensionless (in INaCa_i)
  double h3_ss;   // dimensionless (in INaCa_i)
  double h4_i;   // dimensionless (in INaCa_i)
  double h4_ss;   // dimensionless (in INaCa_i)
  double h5_i;   // dimensionless (in INaCa_i)
  double h5_ss;   // dimensionless (in INaCa_i)
  double h6_i;   // dimensionless (in INaCa_i)
  double h6_ss;   // dimensionless (in INaCa_i)
  double h7_i;   // dimensionless (in INaCa_i)
  double h7_ss;   // dimensionless (in INaCa_i)
  double h8_i;   // dimensionless (in INaCa_i)
  double h8_ss;   // dimensionless (in INaCa_i)
  double h9_i;   // dimensionless (in INaCa_i)
  double h9_ss;   // dimensionless (in INaCa_i)
  double hca;   // dimensionless (in INaCa_i)
  double hna;   // dimensionless (in INaCa_i)
  double k1_i;   // dimensionless (in INaCa_i)
  double k1_ss;   // dimensionless (in INaCa_i)
  double k2_i;   // dimensionless (in INaCa_i)
  double k2_ss;   // dimensionless (in INaCa_i)
  double k3_i;   // dimensionless (in INaCa_i)
  double k3_ss;   // dimensionless (in INaCa_i)
  double k3p_i;   // dimensionless (in INaCa_i)
  double k3p_ss;   // dimensionless (in INaCa_i)
  double k3pp_i;   // dimensionless (in INaCa_i)
  double k3pp_ss;   // dimensionless (in INaCa_i)
  double k4_i;   // dimensionless (in INaCa_i)
  double k4_ss;   // dimensionless (in INaCa_i)
  double k4p_i;   // dimensionless (in INaCa_i)
  double k4p_ss;   // dimensionless (in INaCa_i)
  double k4pp_i;   // dimensionless (in INaCa_i)
  double k4pp_ss;   // dimensionless (in INaCa_i)
  double k5_i;   // dimensionless (in INaCa_i)
  double k5_ss;   // dimensionless (in INaCa_i)
  double k6_i;   // dimensionless (in INaCa_i)
  double k6_ss;   // dimensionless (in INaCa_i)
  double k7_i;   // dimensionless (in INaCa_i)
  double k7_ss;   // dimensionless (in INaCa_i)
  double k8_i;   // dimensionless (in INaCa_i)
  double k8_ss;   // dimensionless (in INaCa_i)
  double x1_i;   // dimensionless (in INaCa_i)
  double x1_ss;   // dimensionless (in INaCa_i)
  double x2_i;   // dimensionless (in INaCa_i)
  double x2_ss;   // dimensionless (in INaCa_i)
  double x3_i;   // dimensionless (in INaCa_i)
  double x3_ss;   // dimensionless (in INaCa_i)
  double x4_i;   // dimensionless (in INaCa_i)
  double x4_ss;   // dimensionless (in INaCa_i)
  double E1;   // dimensionless (in INaK)
  double E2;   // dimensionless (in INaK)
  double E3;   // dimensionless (in INaK)
  double E4;   // dimensionless (in INaK)
  double INaK;   // microA_per_microF (in INaK)
  double JnakK;   // millimolar_per_millisecond (in INaK)
  double JnakNa;   // millimolar_per_millisecond (in INaK)
  double Knai;   // millimolar (in INaK)
  double Knao;   // millimolar (in INaK)
  double P;   // dimensionless (in INaK)
  double Pnak;   // milliS_per_microF (in INaK)
  double a1;   // dimensionless (in INaK)
  double a2;   // dimensionless (in INaK)
  double a3;   // dimensionless (in INaK)
  double a4;   // dimensionless (in INaK)
  double b1;   // dimensionless (in INaK)
  double b2;   // dimensionless (in INaK)
  double b3;   // dimensionless (in INaK)
  double b4;   // dimensionless (in INaK)
  double x1;   // dimensionless (in INaK)
  double x2;   // dimensionless (in INaK)
  double x3;   // dimensionless (in INaK)
  double x4;   // dimensionless (in INaK)
  double GNaL;   // milliS_per_microF (in INaL)
  double INaL;   // microA_per_microF (in INaL)
  double fINaLp;   // dimensionless (in INaL)
  double hLss;   // dimensionless (in INaL)
  double hLssp;   // dimensionless (in INaL)
  double mLss;   // dimensionless (in INaL)
  double thLp;   // millisecond (in INaL)
  double tmL;   // millisecond (in INaL)
  double INab;   // microA_per_microF (in INab)
  double Ahs;   // dimensionless (in INa)
  double INa;   // microA_per_microF (in INa)
  double fINap;   // dimensionless (in INa)
  double h;   // dimensionless (in INa)
  double hp;   // dimensionless (in INa)
  double hss;   // dimensionless (in INa)
  double hssp;   // dimensionless (in INa)
  double jss;   // dimensionless (in INa)
  double mss;   // dimensionless (in INa)
  double thf;   // millisecond (in INa)
  double ths;   // millisecond (in INa)
  double thsp;   // millisecond (in INa)
  double tj;   // millisecond (in INa)
  double tjp;   // millisecond (in INa)
  double tm;   // millisecond (in INa)
  double IpCa;   // microA_per_microF (in IpCa)
  double AiF;   // dimensionless (in Ito)
  double AiS;   // dimensionless (in Ito)
  double Gto;   // milliS_per_microF (in Ito)
  double Ito;   // microA_per_microF (in Ito)
  double ass;   // dimensionless (in Ito)
  double assp;   // dimensionless (in Ito)
  double delta_epi;   // dimensionless (in Ito)
  double dti_develop;   // dimensionless (in Ito)
  double dti_recover;   // dimensionless (in Ito)
  double fItop;   // dimensionless (in Ito)
  double i;   // dimensionless (in Ito)
  double ip;   // dimensionless (in Ito)
  double iss;   // dimensionless (in Ito)
  double ta;   // millisecond (in Ito)
  double tiF;   // millisecond (in Ito)
  double tiF_b;   // millisecond (in Ito)
  double tiFp;   // millisecond (in Ito)
  double tiS;   // millisecond (in Ito)
  double tiS_b;   // millisecond (in Ito)
  double tiSp;   // millisecond (in Ito)
  double Jleak;   // millimolar_per_millisecond (in SERCA)
  double Jup;   // millimolar_per_millisecond (in SERCA)
  double Jupnp;   // millimolar_per_millisecond (in SERCA)
  double Jupp;   // millimolar_per_millisecond (in SERCA)
  double fJupp;   // dimensionless (in SERCA)
  double upScale;   // dimensionless (in SERCA)
  double Acap;   // centimeter_squared (in cell_geometry)
  double Ageo;   // centimeter_squared (in cell_geometry)
  double vcell;   // microliter (in cell_geometry)
  double vjsr;   // microliter (in cell_geometry)
  double vmyo;   // microliter (in cell_geometry)
  double vnsr;   // microliter (in cell_geometry)
  double vss;   // microliter (in cell_geometry)
  double Jdiff;   // millimolar_per_millisecond (in diff)
  double JdiffK;   // millimolar_per_millisecond (in diff)
  double JdiffNa;   // millimolar_per_millisecond (in diff)
  double Bcai;   // dimensionless (in intracellular_ions)
  double Bcajsr;   // dimensionless (in intracellular_ions)
  double Bcass;   // dimensionless (in intracellular_ions)
  double cmdnmax;   // millimolar (in intracellular_ions)
  double Istim;   // microA_per_microF (in membrane)
  double vffrt;   // coulomb_per_mole (in membrane)
  double vfrt;   // dimensionless (in membrane)
  double EK;   // millivolt (in reversal_potentials)
  double EKs;   // millivolt (in reversal_potentials)
  double ENa;   // millivolt (in reversal_potentials)
  double Jrel;   // millimolar_per_millisecond (in ryr)
  double Jrel_inf;   // dimensionless (in ryr)
  double Jrel_inf_temp;   // dimensionless (in ryr)
  double Jrel_infp;   // dimensionless (in ryr)
  double Jrel_temp;   // dimensionless (in ryr)
  double a_rel;   // millisecond (in ryr)
  double a_relp;   // millisecond (in ryr)
  double btp;   // millisecond (in ryr)
  double fJrelp;   // dimensionless (in ryr)
  double tau_rel;   // millisecond (in ryr)
  double tau_rel_temp;   // millisecond (in ryr)
  double tau_relp;   // millisecond (in ryr)
  double tau_relp_temp;   // millisecond (in ryr)
  double Jtr;   // millimolar_per_millisecond (in trans_flux)

  //---------------------------------------------------------------------------
  // Constants
  //---------------------------------------------------------------------------

// this is where all the data from the calling function are passed to f.
Istim = data->IV;

// CaM parameters.
   CaMKo  = data->p[103];  // 0.05;   // dimensionless (in CaMK)
   KmCaM  = data->p[104];  // 0.0015;   // millimolar (in CaMK)
   KmCaMK = data->p[105]; // 0.15;   // millimolar (in CaMK)
   aCaMK  = data->p[106]; // 0.05;   // per_millimolar_per_millisecond (in CaMK)
   bCaMK  = data->p[107]; // 0.00068;   // per_millisecond (in CaMK)

   Kmn = 0.002;   // millimolar (in ICaL) This one is done for SA
//   PCa_b = 0.0001;   // dimensionless (in ICaL) This is a parameter.
   k2n = 1000.0;   // per_millisecond (in ICaL) This one is done for SA

// INaCa parameters except for max INaCa.
   KmCaAct = 150.0e-6;   // millimolar (in INaCa_i)
   kasymm = 12.5;   // dimensionless (in INaCa_i)
   kcaoff = 5.0e3;   // per_millisecond (in INaCa_i)
   kcaon = 1.5e6;   // per_millisecond (in INaCa_i)
   kna1 = 15.0;   // per_millisecond (in INaCa_i)
   kna2 = 5.0;   // per_millisecond (in INaCa_i)
   kna3 = 88.12;   // per_millisecond (in INaCa_i)
   qca = 0.167;   // dimensionless (in INaCa_i)
   qna = 0.5224;   // dimensionless (in INaCa_i)
   wca = 6.0e4;   // dimensionless (in INaCa_i)
   wna = 6.0e4;   // dimensionless (in INaCa_i)
   wnaca = 5.0e3;   // dimensionless (in INaCa_i)

/* INaK parameters up for SA */
/*   H = 1.0e-7;   // millimolar (in INaK) the pH of cytosol
   Khp = 1.698e-7;   // millimolar (in INaK)
   Kki = 0.5;   // per_millisecond (in INaK)
   Kko = 0.3582;   // per_millisecond (in INaK)
   Kmgatp = 1.698e-7;   // millimolar (in INaK)
   Knai0 = 9.073;   // millimolar (in INaK)
   Knao0 = 27.78;   // millimolar (in INaK)
   Knap = 224.0;   // millimolar (in INaK)
   Kxkur = 292.0;   // millimolar (in INaK)
   MgADP = 0.05;   // millimolar (in INaK)
   MgATP = 9.8;   // millimolar (in INaK)

   delta = -0.155;   // millivolt (in INaK)
   eP = 4.2;   // dimensionless (in INaK)
   k1m = 182.4;   // per_millisecond (in INaK)
   k1p = 949.5;   // per_millisecond (in INaK)
   k2m = 39.4;   // per_millisecond (in INaK)
   k2p = 687.2;   // per_millisecond (in INaK)
   k3m = 79300.0;   // per_millisecond (in INaK)
   k3p = 1899.0;   // per_millisecond (in INaK)
   k4m = 40.0;   // per_millisecond (in INaK)
   k4p = 639.0;   // per_millisecond (in INaK)
*/

   H = data->p[76]; // 1.0e-7;   // millimolar (in INaK) the pH of cytosol
   Khp = data->p[77]; // 1.698e-7;   // millimolar (in INaK)
   Kki = data->p[78]; // 0.5;   // per_millisecond (in INaK)
   Kko = data->p[79]; // 0.3582;   // per_millisecond (in INaK)
   Kmgatp = data->p[80]; // 1.698e-7;   // millimolar (in INaK)
   Knai0 = data->p[81]; // 9.073;   // millimolar (in INaK)
   Knao0 = data->p[82]; // 27.78;   // millimolar (in INaK)
   Knap = data->p[83]; // 224.0;   // millimolar (in INaK)
   Kxkur = data->p[84]; // 292.0;   // millimolar (in INaK)
   MgATP = data->p[86]; // 9.8;   // millimolar (in INaK)

// this is not an independent parameter
   MgADP = data->p[85]; // 0.05;   // millimolar (in INaK)


   delta = data->p[87]; // -0.155;   // millivolt (in INaK)
   eP = data->p[88]; // 4.2;   // dimensionless (in INaK)
   k1m = data->p[89]; // 182.4;   // per_millisecond (in INaK)
   k1p = data->p[90]; // 949.5;   // per_millisecond (in INaK)
   k2m = data->p[91]; // 39.4;   // per_millisecond (in INaK)
   k2p = data->p[92]; // 687.2;   // per_millisecond (in INaK)
   k3m = data->p[93]; // 79300.0;   // per_millisecond (in INaK)
   k3p = data->p[94]; // 1899.0;   // per_millisecond (in INaK)
   k4m = data->p[95]; // 40.0;   // per_millisecond (in INaK)
   k4p = data->p[96]; // 639.0;   // per_millisecond (in INaK)

   thL = 200.0;    // millisecond (in INaL)
   Ahf = 0.99;     // dimensionless (in INa)
   hssV1 = 82.9;   // millivolt (in INa)
   hssV2 = 6.086;  // millivolt (in INa)
   mssV1 = 39.57;  // millivolt (in INa)
   mssV2 = 9.871;  // millivolt (in INa)
   mtD1 = 6.765;   // dimensionless (in INa)
   mtD2 = 8.552;   // dimensionless (in INa)
   mtV1 = 11.64;   // millivolt (in INa)
   mtV2 = 34.77;   // millivolt (in INa)
   mtV3 = 77.42;   // millivolt (in INa)
   mtV4 = 5.955;   // millivolt (in INa)


   L = data->p[75];  // 0.01;   // centimeter (in cell_geometry).

   rad = 0.0011;     // centimeter (in cell_geometry)
   celltype = 0.0;   // dimensionless (in environment) 0 is for endo, 1 is for epi, 2 is for mid.

// conductances up for SA.*******************************************************
   PCab   = data->p[62]; // = 2.5e-8;   // milliS_per_microF (in ICab)
   GK1_b  = data->p[63]; // = 0.1908;   // milliS_per_microF (in IK1)
   GKb_b  = data->p[64]; // = 0.003;   // milliS_per_microF (in IKb)

   GKr_b  = data->p[65]; // = 0.046;   // milliS_per_microF (in IKr)
   GKs_b  = data->p[66]; // = 0.0034;   // milliS_per_microF (in IKs)
   Gncx_b = data->p[67]; // = 0.0008;   // milliS_per_microF (in INaCa_i)
   Pnak_b = data->p[68]; // = 30.0;   // milliS_per_microF (in INaK)
// Gnal_b multiplier: control Nai at least.
   GNaL_b = data->p[69]; // = 0.5*0.0075;   // milliS_per_microF (in INaL) To reduce drift of sodium in the cell.

   PNab  = data->p[70]; // = 3.75e-10;   // milliS_per_microF (in INab)
   GNa   = data->p[71]; // = 75.0;   // milliS_per_microF (in INa)
   GpCa  = data->p[72]; // = 0.0005;   // milliS_per_microF (in IpCa)
   Gto_b = data->p[73]; // = 0.02;   // milliS_per_microF (in Ito)
   cm    = data->p[74]; // = 1.0;   // microF_per_centimeter_squared (in intracellular_ions)

   KmCap = data->p[97]; // 0.0005;   // millimolar (in IpCa)


//*******************************************************************************

/*
   cao = 1.8;   // millimolar (in extracellular)
   ko = 5.4;   // millimolar (in extracellular)
   nao = 140.0;   // millimolar (in extracellular)
*/

   cao = data->p[59]; // 1.8;   // millimolar (in extracellular)
   ko = data->p[60]; // 5.4;   // millimolar (in extracellular)
   nao = data->p[61]; // 140.0;   // millimolar (in extracellular)

   BSLmax = 1.124;   // millimolar (in intracellular_ions)
   BSRmax = 0.047;   // millimolar (in intracellular_ions)
   KmBSL = 0.0087;   // millimolar (in intracellular_ions)
   KmBSR = 0.00087;   // millimolar (in intracellular_ions)
   cmdnmax_b = 0.05;   // millimolar (in intracellular_ions)
   csqnmax = 10.0;   // millimolar (in intracellular_ions)
   kmcmdn = 0.00238;   // millimolar (in intracellular_ions)
   kmcsqn = 0.8;   // millimolar (in intracellular_ions)
   kmtrpn = 0.0005;   // millimolar (in intracellular_ions)
   trpnmax = 0.07;   // millimolar (in intracellular_ions)
//   amp = -80.0;   // microA_per_microF (in membrane)
//   duration = 0.5;   // millisecond (in membrane)
   F = 96485.0;   // coulomb_per_mole (in physical_constants)
   R = 8314.0;   // joule_per_kilomole_kelvin (in physical_constants)
   T = data->p[108]; // 310.0;   // kelvin (in physical_constants)
   zca = 2.0;   // dimensionless (in physical_constants)
   zk = 1.0;   // dimensionless (in physical_constants)
   zna = 1.0;   // dimensionless (in physical_constants)
   PKNa = 0.01833;   // dimensionless (in reversal_potentials)
   bt = 4.75;   // millisecond (in ryr)

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

// Aff = 0.6;
   Aff = data->p[25];

   Afs = 1.0-Aff;
   tjca = 75.0;

/*
Parts of all conductances depend of PCa_b. So if I perturb PCa_b, then all these
are dependent parameters. But, the proportions are still a constant.
How do you justify not perturbing the proportions?
*/

// PCa_b
PCa_b = data->p[53]; // = 0.0001;

// proportions
   if (celltype == 1.0)
      PCa = PCa_b*1.2;
   else if (celltype == 2.0)
      PCa = PCa_b*2.5;
   else
      PCa = PCa_b;

   PCap   = data->p[54]*PCa;
   PCaNa  = data->p[55]*PCa;
   PCaK   = data->p[56]*PCa;
   PCaNap = data->p[57]*PCap;
   PCaKp  = data->p[58]*PCap;

   if (celltype == 1.0)
      GK1 = GK1_b*1.2;
   else if (celltype == 2.0)
      GK1 = GK1_b*1.3;
   else
      GK1 = GK1_b;

   if (celltype == 1.0)
      GKb = GKb_b*0.6;
   else
      GKb = GKb_b;

   if (celltype == 1.0)
      GKr = GKr_b*1.3;
   else if (celltype == 2.0)
      GKr = GKr_b*0.8;
   else
      GKr = GKr_b;

   if (celltype == 1.0)
      GKs = GKs_b*1.4;
   else
      GKs = GKs_b;

   Ahs = 1.0-Ahf;
   h10_i = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
   h11_i = nao*nao/(h10_i*kna1*kna2);
   h12_i = 1.0/h10_i;
   k1_i = h12_i*cao*kcaon;
   k2_i = kcaoff;
   k5_i = kcaoff;

   if (celltype == 1.0)
      Gncx = Gncx_b*1.1;
   else if (celltype == 2.0)
      Gncx = Gncx_b*1.4;
   else
      Gncx = Gncx_b;

   h10_ss = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
   h11_ss = nao*nao/(h10_ss*kna1*kna2);
   h12_ss = 1.0/h10_ss;
   k1_ss = h12_ss*cao*kcaon;
   k2_ss = kcaoff;
   k5_ss = kcaoff;
   b1 = k1m*MgADP;
   a2 = k2p;
   a4 = k4p*MgATP/Kmgatp/(1.0+MgATP/Kmgatp);

   if (celltype == 1.0)
      Pnak = Pnak_b*0.9;
   else if (celltype == 2.0)
      Pnak = Pnak_b*0.7;
   else
      Pnak = Pnak_b;

   thLp = 3.0*thL;

   if (celltype == 1.0)
      GNaL = GNaL_b*0.6;
   else
      GNaL = GNaL_b;

   if (celltype == 1.0)
      Gto = Gto_b*4.0;
   else if (celltype == 2.0)
      Gto = Gto_b*4.0;
   else
      Gto = Gto_b;

   if (celltype == 1.0)
      upScale = 1.3;
   else
      upScale = 1.0;

   vcell = 1000.0*3.14*rad*rad*L;
   Ageo = 2.0*3.14*rad*rad+2.0*3.14*rad*L;
   Acap = 2.0*Ageo;
   vmyo = 0.68*vcell;
   vnsr = 0.0552*vcell;
   vjsr = 0.0048*vcell;
   vss = 0.02*vcell;

// printf("%20.20f \n",vss);
   if (celltype == 1.0)
      cmdnmax = cmdnmax_b*1.3;
   else
      cmdnmax = cmdnmax_b;

   a_rel = 0.5*bt;
   btp = 1.25*bt;
   a_relp = 0.5*btp;

   // time: time (millisecond)

/*
// do your algebriac here.
realtype vm_alg;
realtype Q0;
realtype phi1_const;
realtype phi2_const;
realtype phi3_const;

   phi1_const   = Y[30] - cmdnmax*kmcmdn/pow(kmcmdn+Y[30], 1.0)-trpnmax*kmtrpn/pow(kmtrpn+Y[30], 1.0);
   phi2_const   = Y[33] - BSRmax*KmBSR/pow(KmBSR+Y[33], 1.0)-BSLmax*KmBSL/pow(KmBSL+Y[33], 1.0);
   phi3_const   = Y[31] - csqnmax*kmcsqn/pow(kmcsqn+Y[31], 1.0);

// find the Q0 using Istim = 0

//Q0 = 153.4*Y[38] +(F*vmyo*(-Y[36] - Y[34] - 2*phi1_const) \
//	         + F*vss *( Y[37] + Y[35] - 2*phi2_const) \
//	        -2*F*     (Y[32]*vnsr       - phi3_const*vjsr));
//printf("%30.30f\n",Q0);


dY[41] = Istim; // Y41 has units of microA/microF x milli-seconds = nA/microF x s = nC/micro F = mV
// Calculate vm_alg using Q0 and assign to Y[38]. Output vm_alg	   
Q0 = -13841.338579765548274735920131206512;

vm_alg =  (Q0 -(F*vmyo*(-Y[36] - Y[34] - 2*phi1_const) \
	         + F*vss *( Y[37] + Y[35] - 2*phi2_const) \
	        -2*F*     (Y[32]*vnsr       - phi3_const*vjsr)))/153.4 - Y[41];
Y[38] = vm_alg;
data->vm_alg = vm_alg;
// printf("%f %30.30f \n",Istim,Y[41]);
*/


// dependence on CaM like so:
   CaMKb = CaMKo*(1.0-Y[0])/(1.0+KmCaM/Y[33]);
   CaMKa = CaMKb+Y[0];
   dY[0] = aCaMK*CaMKb*(CaMKb+Y[0])-bCaMK*Y[0];


//   dss = 1.0/(1.0+exp(-(Y[38]+3.94)/4.23));
   dss = 1.0/(1.0+exp(-(Y[38]+data->p[0])/data->p[1]));

//   td = 0.6+1.0/(exp(-0.05*(Y[38]+6.0))+exp(0.09*(Y[38]+14.0)));
     td = data->p[6] + 1.0/(exp(-data->p[7]*(Y[38]+data->p[8]))+exp(data->p[9]*(Y[38]+data->p[10])));

   dY[1] = (dss-Y[1])/td;

//   fss = 1.0/(1.0+exp((Y[38]+19.58)/3.696));

    fss = 1.0/(1.0+exp((Y[38]+data->p[2])/data->p[3]));

/* Tauf:
   tff = 7.0+1.0/(0.0045*exp(-(Y[38]+20.0)/10.0)+0.0045*exp((Y[38]+20.0)/10.0));
   tfs = 1000.0+1.0/(0.000035*exp(-(Y[38]+5.0)/4.0)+0.000035*exp((Y[38]+5.0)/6.0));
data->p[11] = 7.0;
data->p[12] = 0.0045;
data->p[13] = 20.0;
data->p[14] = 10.0;
data->p[15] = 0.0045;
data->p[16] = 20.0;
data->p[17] = 10.0;
data->p[18] = 1000.0;
data->p[19] = 0.000035;
data->p[20] = 5.0;
data->p[21] = 4.0;
data->p[22] = 0.000035;
data->p[23] = 5.0;
data->p[24] = 6.0;
*/
   tff = data->p[11] + 1.0/(data->p[12]*exp(-(Y[38]+data->p[13])/data->p[14])+data->p[15]*exp((Y[38]+data->p[16])/data->p[17]));
   tfs = data->p[18]+1.0/(data->p[19]*exp(-(Y[38]+data->p[20])/data->p[21])+data->p[22]*exp((Y[38]+data->p[23])/data->p[24]));



   dY[5] = (fss-Y[5])/tff;
   dY[7] = (fss-Y[7])/tfs;
   f = Aff*Y[5]+Afs*Y[7];
   fcass = fss;

/*
   tfcaf = 7.0+1.0/(0.04*exp(-(Y[38]-4.0)/7.0)+0.04*exp((Y[38]-4.0)/7.0));
   tfcas = 100.0+1.0/(0.00012*exp(-Y[38]/3.0)+0.00012*exp(Y[38]/7.0));
   Afcaf = 0.3+0.6/(1.0+exp((Y[38]-10.0)/10.0));
*/

   tfcaf = data->p[26] + 1.0/( data->p[27] *exp(-(Y[38]- data->p[28] )/ data->p[29] )+ data->p[30] *exp((Y[38]- data->p[31] )/ data->p[32] ));
   tfcas =  data->p[33] +1.0/( data->p[34] *exp(-Y[38]/ data->p[35] )+ data->p[36] *exp(Y[38]/ data->p[37] ));
   Afcaf =  data->p[38] + data->p[39] /(1.0+exp((Y[38]- data->p[40] )/ data->p[41] ));

   Afcas = 1.0-Afcaf;

   dY[2] = (fcass-Y[2])/tfcaf;
   dY[4] = (fcass-Y[4])/tfcas;
   fca = Afcaf*Y[2]+Afcas*Y[4];
   dY[8] = (fcass-Y[8])/tjca;

/*
   tffp = 2.5*tff;
   tfcafp = 2.5*tfcaf;
*/

   tffp = data->p[42]*tff;     // 42
   tfcafp = data->p[43]*tfcaf; // 43

   dY[6] = (fss-Y[6])/tffp;
   fp = Aff*Y[6]+Afs*Y[7];


   dY[3] = (fcass-Y[3])/tfcafp;
   fcap = Afcaf*Y[3]+Afcas*Y[4];
   km2n = Y[8];

//   anca = 1.0/(k2n/km2n+pow(1.0+Kmn/Y[33], 4.0));
anca = 1.0/(data->p[45]/km2n+pow(1.0+data->p[44]/Y[33], 4.0));

//   dY[9] = anca*k2n-Y[9]*km2n;
 dY[9] = anca*data->p[45]-Y[9]*km2n;

   vffrt = Y[38]*F*F/(R*T);
   vfrt  = Y[38]*F  /(R*T);

// gamma cai multiplies the Y[33] (Cass concentration). It is 1 at basal value. I do not know if its evaluated, or a free parameter.
// gamma cai has to be more than 0.341, which is gamma cao.
//  PhiCaL = 4.0*vffrt*(Y[33]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
 PhiCaL = 4.0*vffrt*(data->p[46]*Y[33]*exp(2.0*vfrt)-data->p[47]*cao)/(exp(2.0*vfrt)-1.0);

//  PhiCaNa = 1.0*vffrt*(0.75*Y[37]*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
 PhiCaNa = 1.0*vffrt*(data->p[48]*Y[37]*exp(1.0*vfrt)-data->p[49]*nao)/(exp(1.0*vfrt)-1.0);

//   PhiCaK = 1.0*vffrt*(0.75*Y[35]*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
 PhiCaK = 1.0*vffrt*(data->p[50]*Y[35]*exp(1.0*vfrt)-data->p[51]*ko)/(exp(1.0*vfrt)-1.0);

//   fICaLp = 1.0/(1.0+KmCaMK/CaMKa);
// CamKa depends on Y0, which in turn has all the Cai parameters assosciated with it. So there you are.
    fICaLp = 1.0/(1.0+data->p[52]/CaMKa);

   ICaL  = (1.0-fICaLp)*PCa*PhiCaL*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCap*PhiCaL*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);
   ICaNa = (1.0-fICaLp)*PCaNa*PhiCaNa*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCaNap*PhiCaNa*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);
   ICaK  = (1.0-fICaLp)*PCaK*PhiCaK*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCaKp*PhiCaK*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);

/***********************************end of ical*****************************************************************/


   ICab = PCab*4.0*vffrt*(Y[30]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);


   xk1ss = 1.0/(1.0+exp(-(Y[38]+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
   txk1 = 122.2/(exp(-(Y[38]+127.2)/20.36)+exp((Y[38]+236.8)/69.33));
   dY[10] = (xk1ss-Y[10])/txk1;
   rk1 = 1.0/(1.0+exp((Y[38]+105.8-2.6*ko)/9.493));
   EK = R*T/F*log(ko/Y[34]);
   IK1 = GK1*sqrt(ko)*rk1*Y[10]*(Y[38]-EK);
   xkb = 1.0/(1.0+exp(-(Y[38]-14.48)/18.34));
   IKb = GKb*xkb*(Y[38]-EK);
   xrss = 1.0/(1.0+exp(-(Y[38]+8.337)/6.789));
   txrf = 12.98+1.0/(0.3652*exp((Y[38]-31.66)/3.869)+4.123e-5*exp(-(Y[38]-47.78)/20.38));
   txrs = 1.865+1.0/(0.06629*exp((Y[38]-34.7)/7.355)+1.128e-5*exp(-(Y[38]-29.74)/25.94));
   Axrf = 1.0/(1.0+exp((Y[38]+54.81)/38.21));
   Axrs = 1.0-Axrf;
   dY[11] = (xrss-Y[11])/txrf;
   dY[12] = (xrss-Y[12])/txrs;
   xr = Axrf*Y[11]+Axrs*Y[12];
   rkr = 1.0/(1.0+exp((Y[38]+55.0)/75.0))*1.0/(1.0+exp((Y[38]-10.0)/30.0));
   IKr = GKr*sqrt(ko/5.4)*xr*rkr*(Y[38]-EK);
   xs1ss = 1.0/(1.0+exp(-(Y[38]+11.6)/8.932));
   txs1 = 817.3+1.0/(2.326e-4*exp((Y[38]+48.28)/17.8)+0.001292*exp(-(Y[38]+210.0)/230.0));
   dY[13] = (xs1ss-Y[13])/txs1;
   xs2ss = xs1ss;
   txs2 = 1.0/(0.01*exp((Y[38]-50.0)/20.0)+0.0193*exp(-(Y[38]+66.54)/31.0));
   dY[14] = (xs2ss-Y[14])/txs2;
   KsCa = 1.0+0.6/(1.0+pow(3.8e-5/Y[30], 1.4));
   EKs = R*T/F*log((ko+PKNa*nao)/(Y[34]+PKNa*Y[36]));
   IKs = GKs*KsCa*Y[13]*Y[14]*(Y[38]-EKs);
   mss = 1.0/(1.0+exp(-(Y[38]+mssV1)/mssV2));
   tm = 1.0/(mtD1*exp((Y[38]+mtV1)/mtV2)+mtD2*exp(-(Y[38]+mtV3)/mtV4));
   dY[23] = (mss-Y[23])/tm;
   hss = 1.0/(1.0+exp((Y[38]+hssV1)/hssV2));
   thf = 1.0/(1.432e-5*exp(-(Y[38]+1.196)/6.285)+6.149*exp((Y[38]+0.5096)/20.27));
   ths = 1.0/(0.009794*exp(-(Y[38]+17.95)/28.05)+0.3343*exp((Y[38]+5.73)/56.66));
   dY[18] = (hss-Y[18])/thf;
   dY[19] = (hss-Y[19])/ths;
   h = Ahf*Y[18]+Ahs*Y[19];
   jss = hss;
   tj = 2.038+1.0/(0.02136*exp(-(Y[38]+100.6)/8.281)+0.3052*exp((Y[38]+0.9941)/38.45));
   dY[21] = (jss-Y[21])/tj;
   hssp = 1.0/(1.0+exp((Y[38]+89.1)/6.086));
   thsp = 3.0*ths;
   dY[20] = (hssp-Y[20])/thsp;
   hp = Ahf*Y[18]+Ahs*Y[20];
   tjp = 1.46*tj;
   dY[22] = (jss-Y[22])/tjp;
   fINap = 1.0/(1.0+KmCaMK/CaMKa);
   ENa = R*T/F*log(nao/Y[36]);
   INa = GNa*(Y[38]-ENa)*pow(Y[23], 3.0)*((1.0-fINap)*h*Y[21]+fINap*hp*Y[22]);
   hca = exp(qca*Y[38]*F/(R*T));
   hna = exp(qna*Y[38]*F/(R*T));
   h1_i = 1.0+Y[36]/kna3*(1.0+hna);
   h2_i = Y[36]*hna/(kna3*h1_i);
   h3_i = 1.0/h1_i;
   h4_i = 1.0+Y[36]/kna1*(1.0+Y[36]/kna2);
   h5_i = Y[36]*Y[36]/(h4_i*kna1*kna2);
   h6_i = 1.0/h4_i;
   h7_i = 1.0+nao/kna3*(1.0+1.0/hna);
   h8_i = nao/(kna3*hna*h7_i);
   h9_i = 1.0/h7_i;
   k3p_i = h9_i*wca;
   k3pp_i = h8_i*wnaca;
   k3_i = k3p_i+k3pp_i;
   k4p_i = h3_i*wca/hca;
   k4pp_i = h2_i*wnaca;
   k4_i = k4p_i+k4pp_i;
   k6_i = h6_i*Y[30]*kcaon;
   k7_i = h5_i*h2_i*wna;
   k8_i = h8_i*h11_i*wna;
   x1_i = k2_i*k4_i*(k7_i+k6_i)+k5_i*k7_i*(k2_i+k3_i);
   x2_i = k1_i*k7_i*(k4_i+k5_i)+k4_i*k6_i*(k1_i+k8_i);
   x3_i = k1_i*k3_i*(k7_i+k6_i)+k8_i*k6_i*(k2_i+k3_i);
   x4_i = k2_i*k8_i*(k4_i+k5_i)+k3_i*k5_i*(k1_i+k8_i);
   E1_i = x1_i/(x1_i+x2_i+x3_i+x4_i);
   E2_i = x2_i/(x1_i+x2_i+x3_i+x4_i);
   E3_i = x3_i/(x1_i+x2_i+x3_i+x4_i);
   E4_i = x4_i/(x1_i+x2_i+x3_i+x4_i);
   allo_i = 1.0/(1.0+pow(KmCaAct/Y[30], 2.0));
   JncxNa_i = 3.0*(E4_i*k7_i-E1_i*k8_i)+E3_i*k4pp_i-E2_i*k3pp_i;
   JncxCa_i = E2_i*k2_i-E1_i*k1_i;
   INaCa_i = 0.8*Gncx*allo_i*(zna*JncxNa_i+zca*JncxCa_i);
   h1_ss = 1.0+Y[37]/kna3*(1.0+hna);
   h2_ss = Y[37]*hna/(kna3*h1_ss);
   h3_ss = 1.0/h1_ss;
   h4_ss = 1.0+Y[37]/kna1*(1.0+Y[37]/kna2);
   h5_ss = Y[37]*Y[37]/(h4_ss*kna1*kna2);
   h6_ss = 1.0/h4_ss;
   h7_ss = 1.0+nao/kna3*(1.0+1.0/hna);
   h8_ss = nao/(kna3*hna*h7_ss);
   h9_ss = 1.0/h7_ss;
   k3p_ss = h9_ss*wca;
   k3pp_ss = h8_ss*wnaca;
   k3_ss = k3p_ss+k3pp_ss;
   k4p_ss = h3_ss*wca/hca;
   k4pp_ss = h2_ss*wnaca;
   k4_ss = k4p_ss+k4pp_ss;
   k6_ss = h6_ss*Y[33]*kcaon;
   k7_ss = h5_ss*h2_ss*wna;
   k8_ss = h8_ss*h11_ss*wna;
   x1_ss = k2_ss*k4_ss*(k7_ss+k6_ss)+k5_ss*k7_ss*(k2_ss+k3_ss);
   x2_ss = k1_ss*k7_ss*(k4_ss+k5_ss)+k4_ss*k6_ss*(k1_ss+k8_ss);
   x3_ss = k1_ss*k3_ss*(k7_ss+k6_ss)+k8_ss*k6_ss*(k2_ss+k3_ss);
   x4_ss = k2_ss*k8_ss*(k4_ss+k5_ss)+k3_ss*k5_ss*(k1_ss+k8_ss);
   E1_ss = x1_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E2_ss = x2_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E3_ss = x3_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E4_ss = x4_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   allo_ss = 1.0/(1.0+pow(KmCaAct/Y[33], 2.0));
   JncxNa_ss = 3.0*(E4_ss*k7_ss-E1_ss*k8_ss)+E3_ss*k4pp_ss-E2_ss*k3pp_ss;
   JncxCa_ss = E2_ss*k2_ss-E1_ss*k1_ss;
   INaCa_ss = 0.2*Gncx*allo_ss*(zna*JncxNa_ss+zca*JncxCa_ss);
   Knai = Knai0*exp(delta*Y[38]*F/(3.0*R*T));
   Knao = Knao0*exp((1.0-delta)*Y[38]*F/(3.0*R*T));
   P = eP/(1.0+H/Khp+Y[36]/Knap+Y[34]/Kxkur);


   a1 = k1p*pow(Y[36]/Knai, 3.0)/(pow(1.0+Y[36]/Knai, 3.0)+pow(1.0+Y[34]/Kki, 2.0)-1.0);
   b2 = k2m*pow(nao/Knao, 3.0)/(pow(1.0+nao/Knao, 3.0)+pow(1.0+ko/Kko, 2.0)-1.0);
   a3 = k3p*pow(ko/Kko, 2.0)/(pow(1.0+nao/Knao, 3.0)+pow(1.0+ko/Kko, 2.0)-1.0);
   b3 = k3m*P*H/(1.0+MgATP/Kmgatp);
   b4 = k4m*pow(Y[34]/Kki, 2.0)/(pow(1.0+Y[36]/Knai, 3.0)+pow(1.0+Y[34]/Kki, 2.0)-1.0);
   x1 = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
   x2 = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
   x3 = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
   x4 = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
   E1 = x1/(x1+x2+x3+x4);
   E2 = x2/(x1+x2+x3+x4);
   E3 = x3/(x1+x2+x3+x4);
   E4 = x4/(x1+x2+x3+x4);
   JnakNa = 3.0*(E1*a3-E2*b3);
   JnakK = 2.0*(E4*b1-E3*a1);
   INaK = Pnak*(zna*JnakNa+zk*JnakK);






   mLss = 1.0/(1.0+exp(-(Y[38]+42.85)/5.264));
   tmL = tm;
   dY[17] = (mLss-Y[17])/tmL;
   hLss = 1.0/(1.0+exp((Y[38]+87.61)/7.488));
   dY[15] = (hLss-Y[15])/thL;
   hLssp = 1.0/(1.0+exp((Y[38]+93.81)/7.488));
   dY[16] = (hLssp-Y[16])/thLp;
   fINaLp = 1.0/(1.0+KmCaMK/CaMKa);
   INaL = GNaL*(Y[38]-ENa)*Y[17]*((1.0-fINaLp)*Y[15]+fINaLp*Y[16]);
   INab = PNab*vffrt*(Y[36]*exp(vfrt)-nao)/(exp(vfrt)-1.0);


   IpCa = GpCa*Y[30]/(KmCap+Y[30]);



   ass = 1.0/(1.0+exp(-(Y[38]-14.34)/14.82));
   ta = 1.0515/(1.0/(1.2089*(1.0+exp(-(Y[38]-18.4099)/29.3814)))+3.5/(1.0+exp((Y[38]+100.0)/29.3814)));
   dY[24] = (ass-Y[24])/ta;
   iss = 1.0/(1.0+exp((Y[38]+43.94)/5.711));

   if (celltype == 1.0)
      delta_epi = 1.0-0.95/(1.0+exp((Y[38]+70.0)/5.0));
   else
      delta_epi = 1.0;

   tiF_b = 4.562+1.0/(0.3933*exp(-(Y[38]+100.0)/100.0)+0.08004*exp((Y[38]+50.0)/16.59));
   tiS_b = 23.62+1.0/(0.001416*exp(-(Y[38]+96.52)/59.05)+1.78e-8*exp((Y[38]+114.1)/8.079));
   tiF = tiF_b*delta_epi;
   tiS = tiS_b*delta_epi;
   AiF = 1.0/(1.0+exp((Y[38]-213.6)/151.2));
   AiS = 1.0-AiF;
   dY[26] = (iss-Y[26])/tiF;
   dY[28] = (iss-Y[28])/tiS;
   i = AiF*Y[26]+AiS*Y[28];
   assp = 1.0/(1.0+exp(-(Y[38]-24.34)/14.82));
   dY[25] = (assp-Y[25])/ta;
   dti_develop = 1.354+1.0e-4/(exp((Y[38]-167.4)/15.89)+exp(-(Y[38]-12.23)/0.2154));
   dti_recover = 1.0-0.5/(1.0+exp((Y[38]+70.0)/20.0));
   tiFp = dti_develop*dti_recover*tiF;
   tiSp = dti_develop*dti_recover*tiS;
   dY[27] = (iss-Y[27])/tiFp;
   dY[29] = (iss-Y[29])/tiSp;
   ip = AiF*Y[27]+AiS*Y[29];
   fItop = 1.0/(1.0+KmCaMK/CaMKa);
   Ito = Gto*(Y[38]-EK)*((1.0-fItop)*Y[24]*i+fItop*Y[25]*ip);


// the CaM parameters are set in the CaM parameter location.
/*
   data->p[98]  = 0.004375;
   data->p[99]  = 0.00092;
   data->p[100] = 2.75*0.004375;
   data->p[101] = 0.00075;
   data->p[102] = 0.0039375/15.0;
*/

/*
   Jupnp = upScale* 0.004375*Y[30]/(Y[30]+0.00092);
   Jupp  = upScale* 2.75*0.004375*Y[30]/(Y[30]+0.00075);
   fJupp = 1.0/(1.0+KmCaMK/CaMKa);
   Jleak = 0.0039375*Y[32]/15.0;
   Jup   = (1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
*/

   Jupnp = upScale* data->p[98] *Y[30]/(Y[30]+data->p[99]);
   Jupp  = upScale* data->p[100]*Y[30]/(Y[30]+data->p[101]);
   fJupp = 1.0/(1.0+KmCaMK/CaMKa);
   Jleak = data->p[102]*Y[32];
   Jup   = (1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;


   JdiffNa = (Y[37]-Y[36])/2.0;
   JdiffK = (Y[35]-Y[34])/2.0;
   Jdiff = (Y[33]-Y[30])/0.2;

   dY[36] = -(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap*cm/(F*vmyo)+JdiffNa*vss/vmyo;

// printf("Nai: %30.30f %30.30f\n",Y[36],3.0*INaK);

   dY[37] = -(ICaNa+3.0*INaCa_ss)*cm*Acap/(F*vss)-JdiffNa;

   dY[34] = -(Ito+IKr+IKs+IK1+IKb-Istim-2.0*INaK)*cm*Acap/(F*vmyo)+JdiffK*vss/vmyo;
//	printf("Ki: %30.30f\n",-2.0*INaK);
   dY[35] = -ICaK*cm*Acap/(F*vss)-JdiffK;

   Bcai   = 1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+Y[30], 2.0)+trpnmax*kmtrpn/pow(kmtrpn+Y[30], 2.0));
   Bcass  = 1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+Y[33], 2.0)+BSLmax*KmBSL/pow(KmBSL+Y[33], 2.0));
   Bcajsr = 1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+Y[31], 2.0));

   dY[30]        = Bcai*(-(IpCa+ICab-2.0*INaCa_i)*cm*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);
   fJrelp        = 1.0/(1.0+KmCaMK/CaMKa);
   Jrel          = (1.0-fJrelp)*Y[39]+fJrelp*Y[40];
   dY[33]        = Bcass*(-(ICaL-2.0*INaCa_ss)*cm*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);
   Jtr           = (Y[32]-Y[31])/100.0;
   dY[32]        = Jup-Jtr*vjsr/vnsr;
   dY[31]        = Bcajsr*(Jtr-Jrel);
   dY[38]        = -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim);
   Jrel_inf_temp = a_rel*-ICaL/(1.0+1.0*pow(1.5/Y[31], 8.0));

   if (celltype == 2.0)
      Jrel_inf = Jrel_inf_temp*1.7;
   else
      Jrel_inf = Jrel_inf_temp;

   tau_rel_temp = bt/(1.0+0.0123/Y[31]);

   if (tau_rel_temp < 0.001)
      tau_rel = 0.001;
   else
      tau_rel = tau_rel_temp;

   dY[39] = (Jrel_inf-Y[39])/tau_rel;
   Jrel_temp = a_relp*-ICaL/(1.0+pow(1.5/Y[31], 8.0));

   if (celltype == 2.0)
      Jrel_infp = Jrel_temp*1.7;
   else
      Jrel_infp = Jrel_temp;

   tau_relp_temp = btp/(1.0+0.0123/Y[31]);

   if (tau_relp_temp < 0.001)
      tau_relp = 0.001;
   else
      tau_relp = tau_relp_temp;

   dY[40] = (Jrel_infp-Y[40])/tau_relp;

/***********************************************************************/

 for(ii=0;ii<NEQ;ii++)
  Ith(ydot,ii+1) = dY[ii];

data->ICaL = ICaL;
data->IKr = IKr;
data->IK1 = IK1;
data->IKs = IKs;
data->INa = INa;
data->Ito = Ito;

data->dvdtf = dY[38];

 return(0);
}

/*
 * error weights. this does not help much with eps_3_0
 */

static int ewt(N_Vector u, N_Vector w, void *user_data)
{
  int i;
  for (i=1; i<=NEQ; i++) Ith(w,i) = 1.0/(RTOL * ABS(Ith(u,i)) + ATOL);  
  return(0);
}


static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
