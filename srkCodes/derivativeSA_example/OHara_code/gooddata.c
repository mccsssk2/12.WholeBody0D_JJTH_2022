/*
SK. Jan 14 2015.
This file is to be run on the MI simulation from ord.c
This file will sift through the simulation and output
the data that is within a reasonable range. This output
is then fed to the MI and Sobie indices.
*/


// my standard header and #defines, as of 22 Dec. 2013
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#define NUM_DATA 133

int main(int argc, char *argv[]){

double p[NUM_DATA];
double temp;
int i, j, k, ret;

FILE *output, *input;

input = fopen("together_temp.dat","r");
output = fopen("goodData.dat","w+");

j = 0;

while(1){
	for(i=0;i<NUM_DATA;i++)
		ret = fscanf(input,"%lf ",&p[i]);
		if(ret!=1) break;
		else{

/*
//						111      112      113                           114                     115                116
fprintf(output,"%d %f %f %f %f %f %f %f %f ",apd_counter,pcl,restingPotential[apd_counter],dvdtmax[apd_counter],maxPotential[apd_counter],apd90[apd_counter],apd30[apd_counter],maxcai[apd_counter],mincai[apd_counter]);
fprintf(output,"%f %f %f %f ",maxnai[apd_counter],minnai[apd_counter],maxki[apd_counter],minki[apd_counter]);
fprintf(output,"%f %f %f %f ",maxcam[apd_counter],mincam[apd_counter],maxcajsr[apd_counter],mincajsr[apd_counter]);
fprintf(output,"%f %f %f %f ",maxcansr[apd_counter],mincansr[apd_counter],maxcass[apd_counter],mincass[apd_counter]);
fprintf(output,"%f %f ",caimax_time[apd_counter],caimin_time[apd_counter]); // calculate Cai relaxation from here. last entry is 133
*/
/* sensibility constraints. I will analyse only good data. */
			if(((int)p[110]==147)&&                 /* unless the data reached 150 beats, I dont want it. */
			        (p[112]>-1000&&p[112]<10000) &&    /* resting */
				(p[113]>-200 &&p[113]<10000)&&    /* dvdtmax */
				(p[114]>-200 &&p[114]<10000)&&    /* vmax    */
				(p[115]>0    &&p[115]<10000)&&    /* APD90   */
				(p[116]>0    &&p[116]<10000)&&    /* APD30   */
				(p[117]>0    &&p[117]<10000)&&    /* maxCai  */
				(p[118]>0    &&p[118]<10000)&&    /* minCai   */
				(p[119]>0    &&p[119]<10000)&&     
				(p[120]>0    &&p[120]<10000)&&
				(p[121]>0    &&p[121]<10000)&&     
				(p[122]>0    &&p[122]<10000)&&
				(p[123]>0    &&p[123]<10000)&&     
				(p[124]>0    &&p[124]<10000)&&
				(p[125]>0    &&p[125]<10000)&&
				(p[126]>0    &&p[126]<10000)&&     
				(p[127]>0    &&p[127]<10000)&&
				(p[128]>0    &&p[128]<10000)&&     
				(p[129]>0    &&p[129]<10000)&&
				(p[130]>0    &&p[130]<10000)&&
				(p[131]>0    &&p[131]<200000)&&     
				(p[132]>0    &&p[132]<200000)
				){
					for(i=0;i<NUM_DATA;i++)
					ret = fprintf(output,"%40.40f ",p[i]);
					      fprintf(output,"\n");
				j++;
			}
		}
} // end of while
printf("Total number of data lines read: %d \n",j);


return 0;
} // end of main.
