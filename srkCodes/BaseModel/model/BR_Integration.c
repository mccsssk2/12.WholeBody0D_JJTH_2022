//  12 October, 2020
//  Author: Timothy Hunter
//  Lead Developer: Sanjay R. Kharche
//  Part of PM3
//  This code includes implementation of the baroreflex model from Lin et al.
//  2012.
//
//  Ref. DOI: 10.1177/0954411912451823
//  Constants are from Lin et al. 2012

// JJJ Dec 30
// Calculated now in output_ursino.c as (SysP + 2DiasP)/3

double MAP          =  data->MAP; // Ith(y_ursino, 13);

double Pbco2        =  data->p_ursino[109];
double Pbo2         =  data->p_ursino[110];

double x0_P_aff     =  data->p_ursino[158];
double x0_temp1     =  data->p_ursino[159];

//  Table 3 - Baroreflex control parameters
//  Afferent compartment:
double tau_aff      =  data->p_ursino[124]; //  units: seconds
double G_aff        =  data->p_ursino[125];

//  Central Compartment:
double tau_c        =  data->p_ursino[126]; //  Same as tau_aff,  units: s

//  Efferent compartment:
double S_p      =  data->p_ursino[127];
double PNA_max  =  data->p_ursino[128]; //  units: Hz
double PNA_min  =  data->p_ursino[129]; //  units: Hz
double S_s      =  data->p_ursino[130];
double SNA_max  =  data->p_ursino[131]; //  units: Hz
double SNA_min  =  data->p_ursino[132]; //  units: Hz

//  Table A3
double k1 =  data->p_ursino[133];
double k2 =  data->p_ursino[134];
double k3 =  data->p_ursino[135];
double k4 =  data->p_ursino[136];
double k5 =  data->p_ursino[137];

//  sys1 = c2d(tf(G_aff,[tau_aff 1]),DELTAT);
//  [A1, B1, C1, D1] = tf2ss(G_aff,[tau_aff 1]);
double A1 = exp((-1/tau_aff) * DELTAT);
double B1 = (1-A1)*tau_aff;
double C1 = G_aff/tau_aff;
double D1 = 0.0;
double x1_P_aff    = A1 * x0_P_aff + B1 * MAP;
double P_aff       = C1 * x0_P_aff + D1 * MAP;

//  Central Compartment
//  Pbco2 = 30; //  TEMPORARY
//  Pbo2 = 104;
double deltaMAP;
deltaMAP = 0.0; // in case none of the conditions below are satisfied.
if 		((Pbco2 > 40) && (Pbo2 < 104))			deltaMAP = k1 + k2 * Pbco2 + k3 / Pbo2	;
else if 	((Pbco2 <= 40) && (Pbo2 < 104))			deltaMAP = k1 + k2 * 40.0 + k3 / Pbo2		;
else if 	((Pbco2 > 40) && (Pbo2 >= 104))			deltaMAP = k1 + k2 * Pbco2 + k3 / 104.0	;


double P_demand = 90.0 + deltaMAP/100.0; //  Desrcibed in fig. 2 and on page 793

//  temp1 is part of eq. 13
//  sys2 = c2d(tf(1,[tau_c 1]),DELTAT); //  for reference, the tf is of the form 1/(tau_c*s + 1) with s as the laplace variable
//  [A2, B2, C2, D2] = tf2ss(sys2.Numerator{1},sys2.Denominator{1});
double A2 = exp((-1/tau_c) * DELTAT);
double B2 = (1-A2)*tau_c;
double C2 = 1/tau_c;
double D2 = 0;
double x1_temp1    = A2 * x0_temp1 + B2 * P_demand;
double temp1       = C2 * x0_temp1 + D2 * P_demand;

double P_error = temp1 - P_aff ;//  eq. 13


//  Efferent Compartment:
//  eq. 31 chemoreceptor operating point
double deltaG_SNA;
if (Pbco2 > 40)
    deltaG_SNA = k4 * Pbco2 + k5;
else
    deltaG_SNA = k4 * 40 + k5;


double k_s = (SNA_max - SNA_min) / (4 * S_s); //  eq. 15
double SNA = ((SNA_max + SNA_min * exp(P_error/k_s))/(1 + exp(P_error/k_s))) * (1 + deltaG_SNA/100); //  eq. 14

double S_v = S_p; //  paper uses inconsistent notation
double k_v = (PNA_max - PNA_min) / (4 * S_v); //  eq. 17
double PNA = (PNA_max - PNA_min * exp(P_error/k_v)) / (1 + exp(P_error/k_v)); //  eq. 16

data->p_ursino[158] = x1_P_aff;
data->p_ursino[159] = x1_temp1;
data->p_ursino[107] = SNA;
data->p_ursino[108] = PNA;



// fprintf(S_P_NA, "%f\t%f\t%f\n", tout,SNA, PNA);
