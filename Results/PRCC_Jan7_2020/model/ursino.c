/* Sanjay R. Kharche.
The combined Heldt-dialysis model, nominally called Ursino model. This model combines
Heldt whole body circulation and Ursino dialysis unit into a 0D discription.
The names of these files are misnomers, if anything, they should be called Lim and Eu Bo Shim.
References:
The composite model is desciribed in this paper. In the pdfs/ directory, refs refer to references in this paper:
1) Ki Moo Lim, Sung Wook Choi, Byung Goo Min, and Eun Bo Shim.  Numerical Simulation of the Effect of Sodium Profile on Cardiovascular Response to Hemodialysis. Yonsei Med J 49(4):581 - 591, 2008.
2) Heldt 2002. Journal of Applied Physiology.
3) shear from Secomb.
4) detailed kidney as a generalization of single resisitance-capacitor.

Further, papers by Secomb have been used to implement myogenic response and shear.

Units in this model are:
time: seconds.
volume: Liters.
flow: ml/s.
pressure: mmHg.
concentration: mmol
resistance: mmHg*s/ml
capacitance: ml/mmHg
elastance:
*/
// my standard headers and functions.
#include "ursino.h"

int main (int argc, char *argv[])
{
// exit(-1);
void *cvode_mem;  // pointer to memory: the full state lives here.
realtype t, tout;
int iout, NOUT, retval, i;
char *str;
FILE *output, *cardiac_output_file;
UserData 	data; // instance pointer.
data 		= (UserData) malloc(sizeof *data); // now it is created. // allocated memory to pointer.

/* Create serial vector of length NEQ for I.C. and abstol */
N_Vector 	y_ursino = N_VNew_Serial(NEQ); // allocated memory to pointer.

Ith(y_ursino, 1)   = 25;    		// % Vic, 	units: L
Ith(y_ursino, 2)   = 11;    		// % Vis, 	units: L
Ith(y_ursino, 3)   = 3.25;  		// % Vpl, 	units: L
Ith(y_ursino, 4)   = 8.97;  		// % Pup, 	units: mmHg
Ith(y_ursino, 5)   = 11.32;  		// % Pk, 		units: mmHg. This is now 12 variables.
Ith(y_ursino, 6)   = 10.11;     // % Psp, 	units: mmHg
Ith(y_ursino, 7)   = 12.24;  		// % Pll, 	units: mmHg
Ith(y_ursino, 8)   = 3.81;  		// % Pab, 	units: mmHg
Ith(y_ursino, 9)   = -5;    		// % Pth, 	units: mmHg. // % is pth a variable, or a constant? 5 Sept. 2020.
Ith(y_ursino, 10)  = 10.0;  		// % Cl, 	units: ml/mmHg
Ith(y_ursino, 11)  = 20.0;  		// % Cr, 		units: ml/mmHg
Ith(y_ursino, 12)  = 12.78;  		// % Pl, left ventricular pressure, units: mmHg
Ith(y_ursino, 13)  = 91;  			// % Pa, aortic pressure, units: mmHg.
Ith(y_ursino, 14)  = 3.47;     	// % Psup, 		units: mmHg
Ith(y_ursino, 15)  = 3.31;     	// % Pinf, 		units: mmHg
Ith(y_ursino, 16)  = 2.60;   		// % Pr, 			units: mmHg. right ventricle pressure.
Ith(y_ursino, 17)  = 15.66;  		// % Ppa,			units: mmHg. Pulmonary artery pressure.
Ith(y_ursino, 18)  = 12.99;  		// % Ppv, 		units: mmHg. Pulmonary vein pressure.
Ith(y_ursino, 19)  = 100;   		// % Muic, 		units: mmol
Ith(y_ursino, 20)  = 250.0; 		// % Mnaic, 	units: mmol
Ith(y_ursino, 21)  = 3535.0; 		// % Mkic, 		units: mmol
Ith(y_ursino, 22)  = 84.0;  		// % Mclic, 	units: mmol
Ith(y_ursino, 23)  = 10.0;  		// % MHco3ic, units: mmol
Ith(y_ursino, 24)  = 100.0; 		// % Mhic,	 	units: mmol
Ith(y_ursino, 25)  = 0.0;   		// % Mpic, 		units: mmol
Ith(y_ursino, 26)  = 55.0;  		// % Muex, 		units: mmol
Ith(y_ursino, 27)  = 2130.0; 		// % Mnaex, 	units: mmol
Ith(y_ursino, 28)  = 75.0;  		// % Mkex, 		units: mmol
Ith(y_ursino, 29)  = 1470.0; 		// % Mclex, 	units: mmol
Ith(y_ursino, 30)  = 100.0; 		// % Mhco3ex, units: mmol
Ith(y_ursino, 31)  = 100.0; 		// % Mhex, 		units: mmol
Ith(y_ursino, 32)  = 0.0;   		// % Mpex, 		units: mmol

// new
Ith(y_ursino, 33)  = 1.0; 			// Cla variable, left atrial elastance.
Ith(y_ursino, 34)  = 2.0;   		//  Cra variable, right atrial elastance.

Ith(y_ursino, 35)  = 8.0;   		//  left atrial pressure initial condition, mmHg.
Ith(y_ursino, 36)  = 1.0;   		//  right atrial pressure initial condition, mmHg.

// all kidney pressures, including inlet right and left pressures.
for(i=37;i<=48;i++) Ith(y_ursino, i) = 15.0;

Ith(y_ursino, 49)  = 50.0;   //  Right Kidnet Inlet Pressure, mmHg.
Ith(y_ursino, 50)  = 50.0;   //  Left Kidnet Inlet Pressure, mmHg.

Ith(y_ursino, 51)  = 15.0;   //  Brachiocephalic Aortic Pressure
Ith(y_ursino, 52)  = 15.0;   //  Thoracic Aortic Pressure
Ith(y_ursino, 53)  = 15.0;   //  Abdominal Aortic Pressure

cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
CVodeInit(cvode_mem, f_ursino, 0.0, y_ursino);
CVodeSStolerances(cvode_mem, ATOL, RTOL);
CVDense(cvode_mem, NEQ);
CVodeSetMaxStep(cvode_mem,DELTAT);



for(i = 0; i<N_PARAMETER; i++) data->p_ursino[i] = 0.0;

#include "p_ursino.c"

//*******************************************************************************************************************************
//**************************** END OF BAROREFLEX PARAMETERS *********************************************************************
//*******************************************************************************************************************************


// opening various files for writing
// In future can combine multiple files to save computation time
// str = malloc(32*sizeof(char)); sprintf(str,"output%05d.dat", atoi(argv[1])); 				output = fopen(str,"w+"); 	free(str);

// JJ TH Nov 3
// calculate length of each buffer
int s_l = (int)(data->p_ursino[151]/DELTAT);
int p_l = (int)(data->p_ursino[143]/DELTAT);

// Declaring Buffers
double SNA_buffer[s_l], PNA_buffer[p_l];
for (int s = 0; s < s_l; s++){
	SNA_buffer[s] = 3.0;
}
for (int p = 0; p < p_l; p++){
	PNA_buffer[p] = 3.0;
}
double deltaHR, HR;


//***********************************************************************
// initialise iterators.
if(atoi(argv[54])==1)
NOUT = (int)(dialysisTime/DELTAT);
else
NOUT = (int)(shortTime/DELTAT);

iout = 0; tout = DELTAT;
double cardiac_output = 0.0; // this is time integral of Qlo (LV output flow) over one heart period. for now, it is p11 as the reflex is switched off 10 Oct 2020.

double flowpermin[6]		 		= {0.0};

double sysPressures[7] 			= {0.0};
double PresToCompare[7] 		= {0.0};
double diasPressures[7] 		= {0.0};
double maxFlows[7] 					= {0.0};
double FlowsToCompare[7] 		= {0.0};
double minFlows[7] 					= {0.0};


data->MAP_avg				= 0.0;
data->tau_avg[9]		= 0.0;
data->tau_avg[18]		= 0.0;
data->v1_avg				= 0.0;
data->v2_avg				= 0.0;
data->v3_avg				= 0.0;

// TIME LOOP STARTS HERE
for(iout=0; iout<NOUT;iout++){

  data->p_ursino[1] = tout; // pursino1 is time in seconds. time dependent paramter.

  // ***************************************************************************************************
  // Baroreflex Implementation JJ TH Nov 3
  #include "BR_Integration.c"

  // Update SNA_buffer
  for(int s = s_l-1; s>0; s--) SNA_buffer[s] = SNA_buffer[s-1];
  SNA_buffer[0] = data->p_ursino[107];

  // Update PNA_buffer
  for(int p = p_l-1; p>0; p--) PNA_buffer[p] = PNA_buffer[p-1];
  PNA_buffer[0] = data->p_ursino[108];

  #include "BR_Effector.c"

  deltaHR = (19.64*data->p_ursino[111]) - (17.95*data->p_ursino[112]) - (1.225*pow(data->p_ursino[111],2)) + (1.357*pow(data->p_ursino[112],2)) - (1.523*data->p_ursino[111]*data->p_ursino[112]);
  // HR = HR_0 + deltaHR; % eq. 18 (from tim's matlab. HR_0 is just p(148))
  HR = data->p_ursino[167] + deltaHR; // eq. 18 of paper ABC.
  data->p_ursino[11] = (60.0/HR);

  // Timing Parameters
  if(data->p_ursino[1] >= data->p_ursino[106] + data->p_ursino[11]){
  	data->p_ursino[106]  		= data->p_ursino[1];
  	data->Ts             		= 0.37	*sqrt(data->p_ursino[11]);
  	data->Tasys          		= 0.25	*sqrt(data->p_ursino[11]);
  	data->Tav            		= 0.19	*sqrt(data->p_ursino[11]);
  }
  // END OF BAROREFLEX OPERATIONS


  // ***************************************************************************************************

  CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data. You have time dependent data.
  retval = CVode(cvode_mem, tout, y_ursino, &t, CV_NORMAL); if(check_retval(&retval, "CVode", 1)) exit(-1);

// JJJ Shear Stress

// array of lengths for each flow. Might be better to choose reasonable lengths. But for now all set to Sanjay's constant parameter.
for (i = 0; i<45; i++){
  data->lengths[i] = data->p_ursino[101];
}

// mu is p[100].
for(i = 1; i<45; i++){
  data->radius[i]    	= pow(( (8.0 * data->lengths[i]  * data->p_ursino[100] )/ (M_PI * data->p_R[i]) ), 0.25 );
  data->tau[i]       	= 4.0 * data->p_ursino[100] * (*(data->flows[i]))  / (M_PI * pow(data->radius[i],3));
}

  // Configure output files:
  #include "output_ursino.c"


  tout = tout + DELTAT;

	data->MAP_avg 		-= data->MAP_avg/tout;
	data->MAP_avg 		+= data->MAP/tout;

	data->tau_avg[9] 	-= data->tau_avg[9]/tout;
	data->tau_avg[9] 	+= data->tau[9]/tout;

	data->tau_avg[18] -= data->tau_avg[18]/tout;
	data->tau_avg[18] += data->tau[18]/tout;

	data->v1_avg-= data->v1_avg/tout;
	data->v1_avg+= Ith(y_ursino, 1)/tout;

	data->v2_avg-= data->v2_avg/tout;
	data->v2_avg+= Ith(y_ursino, 2)/tout;

	data->v3_avg-= data->v3_avg/tout;
	data->v3_avg+= Ith(y_ursino, 3)/tout;

} // end of time loop.

// 	fclose(output);

  N_VDestroy_Serial(y_ursino);
  CVodeFree(&cvode_mem);
  free(data);

return 0; // the main must return a success exit code to the middleware.
}
