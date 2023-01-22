// our my parameters. *************************************************************************************************************
data->p_ursino[0] = 0.0; // you never use p0, but it must have a value, lest you confuse the program.
data->p_ursino[1]  = 0.0;            // time, seconds. data->p_Ursino[1] = 0.0;
data->p_ursino[2]  = 150.0;          // Meqic
data->p_ursino[3]  = 150.0;          // Meqex
// why p4 is a constant, and not coming out of the calculation?

// Qf is the rate fluid is removed to the dialyzer
// data->p_ursino[4]  = 0.2083*atof(argv[3]);          // 1/60; // Qf; TEMPORARY from ursino 2000 table 2, p. 209
data->p_ursino[4]  = 0.2083;          // 1/60; // Qf; TEMPORARY from ursino 2000 table 2, p. 209
// Qinf is the rate that fluid is injected from dialyzer
// data->p_ursino[5]  = 0.2*atof(argv[3]);           // Qinf: infused fluid
data->p_ursino[5]  = 0.2;           // Qinf: infused fluid

// data->p_ursino[56] = 3.0*atof(argv[3]); // Q_B: bulk blood flow through the dialyzer
data->p_ursino[56] = 3.0; // Q_B: bulk blood flow through the dialyzer

data->p_ursino[6]  = 0;              // cud: concentration of urea in dialysate; units: mmol/L
data->p_ursino[7]  = 142;            // cnad: concentration of Na in dialysate; units: mmol/L; ref. ursino 1997 p666
data->p_ursino[8]  = 2;              // ckd: concentration of K in dialysate; units: mmol
data->p_ursino[9]  = 0;              // ccld: concentration of Cl in dialysate; units: mmol
data->p_ursino[10] = 35;             // chco3d: concentration of HCO3 in dialysate; units: mmol.
// is pursino11 up for time dependent baroreflex modulation?
// revised pursino11 from 60/60 to 60/70 according to the PhysioNet program.
data->p_ursino[11] = 60.0/70.0;      // baseline duration of cardiac cycle (s) . Severi cycle length: 0.825 s per beat.
data->p_ursino[12] = 12.0/70.0;        // baseline respiration rate; units: s/breath; ref. Heldt 2002
data->p_ursino[13] = 0;              // cuinf: concentration of U in infused fluid; units: mmol/L
data->p_ursino[14] = 0;              // cnainf: concentration of Na in infused fluid; units: mmol/L
data->p_ursino[15] = 0;              // ckinf: concentration of K in infused fluid; units: mmol/L
data->p_ursino[16] = 0;              // cclinf: concentration of Cl in infused fluid; units: mmol/L
// p 17,18,19 need to be argv'd
data->p_ursino[17] = 2.67*atof(argv[3]);           // D_s: dialysance (or clearance) of solute
data->p_ursino[18] = 2.67*atof(argv[3]);           // D_U: dialysance (or clearance) of urea
data->p_ursino[19] = 0.13*atof(argv[3])/60.0;      // D_HCO3: dialysance (or clearance) of HCO3
data->p_ursino[20] = 0.94;           // F_p: plasma water fraction
data->p_ursino[21] = 0.72;           // F_R: RBC water fraction
data->p_ursino[22] = 1;              // gamma_u: fraction of red blood cell water that participates in the transfer through the dialyzer
data->p_ursino[23] = 1;              // R_DU: Donnan ratio for Urea.
data->p_ursino[24] = 0;              // gamma_Na: fraction of red blood cell water that participates in the transfer through the dialyzer
data->p_ursino[25] = 0;              // gamma_K: fraction of red blood cell water that participates in the transfer through the dialyzer
data->p_ursino[26] = 0;              // gamma_Cl: fraction of red blood cell water that participates in the transfer through the dialyzer
data->p_ursino[27] = 0;              // gamma_HCO3: fraction of red blood cell water that participates in the transfer through the dialyzer
data->p_ursino[28] = 25;             // k_Na: mass transfer coefficient for K; units (ml/s)
data->p_ursino[29] = 0.0704;         // beta_Na: mass transfer coefficient for K; units (n/a)
data->p_ursino[30] = 0.0667;         // k_K: mass transfer coefficient for K; units (ml/s)
data->p_ursino[31] = 28.2  ;         // beta_K: mass transfer coefficient for K; units (n/a)
data->p_ursino[32] = 13    ;         // k_U: mass transfer coefficient for U; units (ml/s)
data->p_ursino[33] = 1     ;         // beta_U: mass transfer coefficient for U; units (n/a)
data->p_ursino[34] = 0.004 ;         // k_f: water exchange coefficient; units (L^2 s^-1 mmol^-1)
data->p_ursino[35] = 2.45  ;         // E_is: elastance of the interstitial space; units (mmHg/L)
data->p_ursino[36] = 11    ;         // V_isn: basal volume of interstitial compartments; units (L)
data->p_ursino[37] = 3.25  ;         // V_pln: basal volume of blood plasma; units (L)
data->p_ursino[38] = 7.4   ;         // c_ppln: basal protein concentration in plasma; units (g/dl)
data->p_ursino[39] = 1.37  ;         // c_pisn: basal protein concentration in interstitial compartment; units (g/dl)
data->p_ursino[40] = 1.05;           // Gibbs Donnan ratio for anions. Ursino 2000, table 1.
data->p_ursino[41] = 0.95;           // Gibbs Donnan ratio for cations. Ursino 2000, table 1.
data->p_ursino[42] = 0.2 /61.0;     // eta_hco3: bicarbonate mass transfer coefficient. table 1, Ursino 2000. units: L/s.
data->p_ursino[43] = 0.03 /60.0;      // eta_h: hydrogen ion mass transfer coefficient. table 1, Ursino 2000. units: L/s.
data->p_ursino[44] = 0.4	;          // g_hco3: bicarbonate equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
data->p_ursino[45] = 3.5	;          // g_h: hydrogen ion equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
data->p_ursino[46] = 6.0 / 60.0;       // etaprime_r: reaction velocity, hco3 buffer. table 1, Ursino 2000. see reaction i. units: L^2/s/mmol
data->p_ursino[47] = pow(10.0,(-6.1));// kprime_a: dissociation constant; p209 of Ursino 2000, right column 3rd para.
data->p_ursino[48] = 6.0 / 60;       // etaprimeprime_r: reaction velocity, protein buffer. table 1, Ursino 2000. see reaction ii. units: L^2/s/mmol
data->p_ursino[49] = pow(10.0,(-7.4));// kprimeprime_a:dissociation constant; p209 of Ursino 2000, right column 3rd para.
data->p_ursino[50] = 1.2;            // cco2ic: concentration of CO2 in ic compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.
data->p_ursino[51] = 1.2;            // cco2ex: concentration of CO2 in ex compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.
data->p_ursino[52] = 4;              // cpic0: basal protein concentration in intracellular compartment; units mmol/L, Table 1 ursino 2000
data->p_ursino[53] = -5.97;          // pis0: basal pressure in is compartment; units: mmHg; eq. 10, appendix 1. Ursino 2000. table 1.
data->p_ursino[54] = 0.01;           // La: arterial capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008, table 1
data->p_ursino[55] = 0.062;          // Lv: venous capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008, table 1


/******************************************************************************************************************************************/
// cooling mechanism paramters.
data->p_ursino[57] = 37.5; // physiological temperature is 37.5C.
data->p_ursino[58] = atof(argv[2]); //37.5; // 36.0; // cool temperature. You are wanting a 10% increase in microvascular resistance. Our range is 33C to 36.5C.
data->p_ursino[59] = 2.961; // alpha in J appl physiol 102: 1329; 2007. This is the Q10. Also see scanned notes.
data->p_ursino[60] = 0.08401; // beta in J appl physiol 102: 1329; 2007. This is the exponent. Also see scanned notes.
// the alpha^{-beta*(Tcool - Tphys)} factor.
data->p_ursino[61]  = pow(data->p_ursino[59], - data->p_ursino[60] * (data->p_ursino[58] - data->p_ursino[57]));
// printf("P_61 = %f\n",data->p_ursino[61]);

/******************************************************************************/

data->p_ursino[100] = 3 * pow(10,-5); // blood dynamic viscosity. Letcher et al. Am J Med. 70: 1195-1202. 1981. mmHg-s.

// the 0.1 resistance is a placeholder - find justification when writing paper/doing simulations.
data->p_ursino[101] = 0.1; // length of all microvasculature in cm, microvasculature is with resistance more than 0.1.
data->p_ursino[102] = 10.0; // length of all large vessels in cm, microvasculature is with resistance less than 0.1.

// the ascending aorta does not have a resistance as yet. For now, assign some constant value of radius to the aorta.
// Literature value of ascending aorta diameter: 3 cm (Journal of the American College of Cardiology,  Volume 69, Issue 11 Supplement, March 2017, DOI: 10.1016/S0735-1097(17)35464-5; Paruchuri et al. Cardiology 2015;131:265-272).
data->p_ursino[103] = 0.90; // ascending aorta diameter in cm, this is a large value as compared to what using lengths to work our radius gives in this model. The observed value is closer to 1.5 cm.
data->p_ursino[104] = 0.007; // units PRU. Ascending aorta resistance according to table on p35 of Heldt thesis. Note that written in this way, it presumes certain things and will need to be generalised. This resistance is NOT feeding into the aorta flow yet. There is no aorta ODE yet. Written in this way, ascending aorta shear will depend on LV elastance.

data->p_ursino[105] = 0.0;  // this is loc_t passed back from RHS to the main function thru' userdata. loc_t is from 0 to p[11].




/*******************************************************************************
HEART CHAMBER COMPLIANCES
*******************************************************************************/
// units: ml/mmHg;
// Heldt 2002, page 1243, paragraph 1.
// JJJ TH NOV 3 Took these values from heldt appendix

data->Edias_lv  = 0.13*atof(argv[6]);
data->Esys_lv   = 2.5*atof(argv[7]);

data->Edias_rv  = 0.07*atof(argv[8]);
data->Esys_rv   = 1.3*atof(argv[9]);

data->Edias_la  = 0.5*atof(argv[10]);
data->Esys_la   = 0.61*atof(argv[11]);

data->Edias_ra  = 0.3*atof(argv[12]);
data->Esys_ra   = 0.74*atof(argv[13]);

// data->p_ursino[57] = 37.5; // physiological temperature is 37.5C.
// data->p_ursino[58] = atof(argv[2]);

if((   (int)floor(data->p_ursino[57]-data->p_ursino[58])   )==0){ // If p58 = physiological temp
  data->Esys_lv   = 2.5*atof(argv[7]);
  data->Esys_rv   = 1.3*atof(argv[9]);
}else if((   (int)floor(data->p_ursino[57]-data->p_ursino[58])   )==2){
  data->Esys_lv   = 2.5*atof(argv[7])*0.60;
  data->Esys_rv   = 1.3*atof(argv[9])*0.60;
}else{
  printf("Therapeutic Hypothermia value %f is invalid.\n", data->p_ursino[58]); exit(-1);
}

if(0==1){ // switch on for AF.
  data->Edias_la  = 0.5*atof(argv[10]);
  data->Esys_la   = data->Edias_la; // 0.61*atof(argv[11]);

  data->Edias_ra  = 0.3*atof(argv[12]);
  data->Esys_ra   = data->Edias_ra; // 0.74*atof(argv[13]);
}



/******************************************************************************
RESISTANCE PARAMETERS
******************************************************************************/
// units: (mmHg*s/ml);
// Table 2 of Heldt 2002.

// lowering of temperature increases resistances. p99 multiples all resistances. at 36C, the increase is approx. 14%.
// rescaling of low resistances caused integration to fail around t = 280.
// a bunch of resistances are time dependent.

// data->p_R[] = ; //
data->p_R[0]  = 0.0; // Empty
data->p_R[1]  = 0.06 * (data->p_ursino[61]);         // Rsup:  superior vena cava;
data->p_R[2]  = 0.01;                                // Rab:   abdominal vena cava;
data->p_R[3]  = 0.015;                               // Rinf:  inferior vena cava;
data->p_R[4]  =  0.005;                             // Rao: right heart's atrio-ventricular valve (tricuspid valve)
data->p_R[5]  =  0.003;                             // Rro: Resistance of right heart outlet;
data->p_R[6]  = 0.08 * (data->p_ursino[61]);        // Rp: Resistance of pulmonary arteries;
data->p_R[7]  = 0.01;                               // Rpv: Resistance of pulmonary veins;
data->p_R[8]  = 0.01;                               // Rmv: Left heart's mitral value
data->p_R[9]  = 0.006;                             // Rlo: Resistance of left heart outlet;
data->p_R[10] = 0.003* (data->p_ursino[61]);        // R_Brachiocephalic aorta
data->p_R[11] = 0.011* (data->p_ursino[61]);        // R_Thoracic Aorta
data->p_R[12] = 0.010* (data->p_ursino[61]);        // R_abdominal aorta

data->p_R[13] = 3.9 * (data->p_ursino[61])*atof(argv[29]);  //  Rup1: Resistance of upper body(1);
data->p_R[14] = 0.23 * (data->p_ursino[61]);                //  Rup2: Resistance of upper body (2);
data->p_R[15] = 3.0 * (data->p_ursino[61])*atof(argv[28]);  //  Rsp1: Resistance of splanchic circulation (1);
data->p_R[16] = 0.18 * (data->p_ursino[61]);                //  Rsp2: Resistance of splanchic circulation (2);
data->p_R[17] = 5.0* (data->p_ursino[61]) *atof(argv[14]);                         //  R_Inlet_right Kidneys
data->p_R[18] = 5.0* (data->p_ursino[61]) *atof(argv[15]);                         //  R_Inlet_left Kidneys

for(i = 0; i<12; i++){
  data->p_R[19+i] = 15.0 * (data->p_ursino[61]) * 12 * atof(argv[16+i]);   // R_kidi,1 to R_kidi,12
  if(fabs(atof(argv[16+i])) < pow(10, -3)){
    printf("R_kid1_perterb%d is meaningless.\n", 16+i); exit(-1);
  }}

for(i = 0; i<12; i++) data->p_R[31+i] = 0.3 * (data->p_ursino[61]) * 12;  // R_kido,1 to R_kido,12

data->p_R[43] = 3.6 * (data->p_ursino[61])*atof(argv[30]);                // Rll1: Resistance of legs (1);
data->p_R[44] = 0.3 * (data->p_ursino[61]);                               // Rll2: Resistance of legs (2);

/*******************************************************************************
COMPARTMENT COMPLIANCES
*******************************************************************************/
// units: mL/mmHg;
// Heldt 2002 AJP. p. 1242.

data->p_C[0] = 0.0;    // not used
data->p_C[1] = 15.0*atof(argv[53]);   // Csup: Capacitance of superior vena cava;
data->p_C[2] = 25.0*atof(argv[48]);   // Cab: Capacitance of abdominal veins;
data->p_C[3] = 2.0*atof(argv[52]);    // Cinf: Capacitance of inferior vena cava;
data->p_C[4] = 4.3*atof(argv[50]);    // Cpa: Capacitance of pulmonary arteries;
data->p_C[5] = 8.4*atof(argv[51]);    // Cpv: Capacitance of pulmonary veins;
data->p_C[6] = 0.9*2.0*atof(argv[49]);    // Ca: Capacitance of systemic artery, i.e. aorta;
data->p_C[7] = 0.13;                  // C_Brachiocephalic aorta
data->p_C[8] = 0.21;                  // C_Thoracic Aorta
data->p_C[9] = 0.10;                  // C_abdominal aorta
data->p_C[10] = 8.0*atof(argv[46]);   // Cup: Capacitance of upper body;
data->p_C[11] = 55.0*atof(argv[45]);  // Csp: Splanchnic capacitance;
data->p_C[12] = 5.0*atof(argv[31]);  // C_Inlet_right Kidneys
data->p_C[13] = 5.0*atof(argv[32]);  // C_Inlet_left Kidneys

for(i = 0; i<12; i++){                // Ckid1-12: Kidney capacitance;
  data->p_C[14+i] = (15.0/12.0)*atof(argv[33+i]);
  if(fabs(atof(argv[33+i])) < pow(10, -5)){
    printf("C_kid1_perterb%d is meaningless.\n", 33+i); exit(-1);
  }
}
data->p_C[26] = 19.0*atof(argv[47]);   // Cll: Legs venous capacitance;

/*******************************************************************************
BAROREFLEX PARAMETERS
*******************************************************************************/

data->p_ursino[106]   = 0.0;	//t_cardCycleInit: initiation time of the current cardiac cycle

data->p_ursino[107]   = 0; // SNA
data->p_ursino[108]   = 0; // PNA
data->p_ursino[109]   = 30; // Pbco2; units: mmHg
data->p_ursino[110]   = 87; // Pb02; units: mmHg
data->p_ursino[111]   = 0; // deltaHR_S
data->p_ursino[112]   = 0; // deltaHR_V
data->p_ursino[113]   = 0; // sigma_lv; units: mmHg/ml
data->p_ursino[114]   = 0; // sigma_rv; units: mmHg/ml
data->p_ursino[115]   = 0; // delta_sigma_V; units: ml
data->p_ursino[116]   = 0; // sigma_R; units: mmHg/(ml/s)

// ZPFV
data->p_ursino[117]   = 645;//650; // ZPFV_up: zero-pressure filling volume upper body veins; units: ml;
data->p_ursino[118]   = 30;//150; // ZPFV_kid;
data->p_ursino[119]   = 1146;//1300; // ZPFV_sp;
data->p_ursino[120]   = 716;//350; // ZPFV_ll
data->p_ursino[121]   = 79;//250; // ZPFV_ab
data->p_ursino[122]   = 33;//75; // ZPFV_inf
data->p_ursino[123]   = 16;//10; // ZPFV_sup

// Baroreflex INTEGRATOR Parameters
data->p_ursino[124]  = 0.001;   // tau_aff: units: seconds
data->p_ursino[125]  = 1;       // G_aff
data->p_ursino[126]  = 0.001;   // tau_c: Same as tau_aff,  units: s
data->p_ursino[127]  = 0.0205;  // S_p
data->p_ursino[128]  = 6;       // PNA_max, units: Hz
data->p_ursino[129]  = 0.6;     // PNA_min, units: Hz
data->p_ursino[130]  = -0.0138; // S_s
data->p_ursino[131]  = 4;       // SNA_max, units: Hz
data->p_ursino[132]  = 1.12;    // SNA_min, units: Hz
data->p_ursino[133]  = 13.8;    // k1
data->p_ursino[134]  = 0.182;   // k2
data->p_ursino[135]  = 828;     // k3
data->p_ursino[136]  = 1;       // k4
data->p_ursino[137]  = -18.118; // k5

// Baroreflex EFFECTOR Parameters
double G_factor = 0.8*atof(argv[5]);
data->p_ursino[138]  = 90; 			// G_k_s0, units: beats/min/Hz
data->p_ursino[139]  = 0.28; 		// k_k_s0
data->p_ursino[140]  = 3; 			// T_s, units: s
data->p_ursino[141]  = 60; 			// G_v0, units: beats/min/Hz 45
data->p_ursino[142]  = 0.4; 		// k_v0
data->p_ursino[143]  = 0.5; 		// T_v, units: s
data->p_ursino[144]  = 1.5; 		// tau_sigma_lv, units: s
data->p_ursino[145]  = 2; 			// T_e_lv, units: s
// All parameters with G_factor are baroreflex gains.
data->p_ursino[146]  = 0.45*G_factor; 		// G_eff_lv, units: mmHg/ml/Hz

data->p_ursino[147]  = 1.5; 		// tau_sigma_rv, units: s
data->p_ursino[148]  = 2; 			// T_e_rv, units: s
data->p_ursino[149]  = 0.282*G_factor; 	// G_eff_rv, units: mmHg/ml/Hz
data->p_ursino[150]  = 10; 			// tau_sigma_V, units: s
data->p_ursino[151]  = 5; 			// T_e_V, units: s
data->p_ursino[152]  = -275*G_factor; 		// G_eff_V, units: ml/Hz
data->p_ursino[153]  = 1.5; 		// tau_sigma_R, units: s
data->p_ursino[154]  = 3; 			// T_e_R, units: s
data->p_ursino[155]  = 0.28*G_factor;//0.30; 		// G_eff_R, units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21
data->p_ursino[156]  = 25; 			// tau_s, simplified from equation 23, units: s
data->p_ursino[157]  = 0.8; 		// tau_v, simplified from eq. 28; units: s
// state vectors for baroreflex
data->p_ursino[158]  	= 1; 			// x0_P_aff
data->p_ursino[159]  	= 1; 			// x0_temp1
data->p_ursino[160] 	= 1; 			// x1_deltaHR_s;
data->p_ursino[161] 	= 1; 			// x1_deltaHR_v;
data->p_ursino[162] 	= 1; 			// x1_sigma_lv;
data->p_ursino[163] 	= 1; 			// x1_sigma_rv;
data->p_ursino[164] 	= 1; 			// x1_sigma_V;
data->p_ursino[165] 	= 1; 			// x1_sigma_R;
data->p_ursino[166] 	= 0; 			// sigma_V


// JJJ Dec 30
// Dialysis increases intrinsic HR by 20%
if(atoi(argv[3])==0){     // No Dialysis
  data->p_ursino[167] 	= 70*atof(argv[4]); 		// HR0 units: bpm
}else if(atoi(argv[3])==1){  // Dialysis On
  data->p_ursino[167] 	= 70*atof(argv[4])*1.2; 		// HR0 units: bpm
}else{
  // data->p_ursino[167] 	= 70*atof(argv[4]); 		// HR0 units: bpm
  printf("Argv[3] must be 0 or 1\n"); exit(-1);
}

data->flows[0]  = NULL;
data->flows[1]  = &data->Qsup;
data->flows[2]  = &data->Qab;
data->flows[3]  = &data->Qinf;
data->flows[4]  = &data->Qrao;
data->flows[5]  = &data->Qro;
data->flows[6]  = &data->Qpa;
data->flows[7]  = &data->Qli;
data->flows[8]  = &data->Qlao;
data->flows[9]  = &data->Qlo;
data->flows[10] = &data->Q_ao_bra;
data->flows[11] = &data->Q_ao_tho;
data->flows[12] = &data->Q_ao_abd;
data->flows[13] = &data->Qupi;
data->flows[14] = &data->Qupo;
data->flows[15] = &data->Qsp1;
data->flows[16] = &data->Qsp2;
data->flows[17] = &data->QkRi;
data->flows[18] = &data->QkLi;
for(i = 0; i<24; i++){
  data->flows[19+i] = &data->Q_kidneys[i];
}
data->flows[43] = &data->Qll1;
data->flows[44] = &data->Qll2;

data->DialysisOnOff = atoi(argv[3]);
