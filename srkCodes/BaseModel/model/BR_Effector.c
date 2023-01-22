// Effector fxn
// This function controlls the action on the effector sites of the
// baroreflex. Effector sites include HR (sympathetic and parasympathetic),
// venous tone (zero pressure filling volume), Arterial tone (resistance),
// and right and left ventricular contractility.
//
// references:
// Lin et al. DOI: 10.1177/0954411912451823
//
// "Autonomic control of cardiac pacemaker activity and atrioventricular
// transmission" Levy and Zeiske 1969
//
// all values and equations are from Lin unless otherwise stated
//

//double  DELTAT = data->p_ursino[86+19]; // Not needed


// State vectors
double x0_deltaHR_s    = data->p_ursino[141+19];
double x0_deltaHR_v    = data->p_ursino[142+19];
double x0_sigma_lv     = data->p_ursino[143+19];
double x0_sigma_rv     = data->p_ursino[144+19];
double x0_sigma_V      = data->p_ursino[145+19];
double x0_sigma_R      = data->p_ursino[146+19];

// Table 3 baroreflex control parameters
double G_k_s0     = data->p_ursino[119+19]; // units: beats/min/Hz
double k_k_s0     = data->p_ursino[120+19];
double T_s        = data->p_ursino[121+19]; // units: s

double G_v0       = data->p_ursino[122+19]; // units: beats/min/Hz 45
double k_v0       = data->p_ursino[123+19];
double T_v        = data->p_ursino[124+19]; // units: s

// Myocardial constriction
double tau_sigma_lv = data->p_ursino[125+19]; // units: s
double T_e_lv       = data->p_ursino[126+19]; // units: s
double G_eff_lv     = data->p_ursino[127+19]; // units: mmHg/ml/Hz

double tau_sigma_rv = data->p_ursino[128+19]; // units: s
double T_e_rv       = data->p_ursino[129+19]; // units: s
double G_eff_rv     = data->p_ursino[130+19]; // units: mmHg/ml/Hz

// Veins
double tau_sigma_V  = data->p_ursino[131+19]; // units: s
double T_e_V        = data->p_ursino[132+19]; // units: s
double G_eff_V      = data->p_ursino[133+19]; // units: ml/Hz

// Arterioles
double tau_sigma_R  = data->p_ursino[134+19]; // units: s
double T_e_R        = data->p_ursino[135+19]; // units: s
double G_eff_R      = data->p_ursino[136+19]; // units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21

double tau_s        = data->p_ursino[137+19]; // simplified from equation 23, units: s

double k_s2   = G_k_s0 / k_k_s0; // eq. 22
double G_s    = 1 * (1 - exp(-k_s2 * SNA_buffer[1-1])); // eq. 21 EDIT gain value changed to 1
//double  sys3 = c2d(tf(G_s,[tau_s 1]),DELTAT);
// [A3, B3, C3, D3] = tf2ss(sys3.Numerator{1},sys3.Denominator{1});
double A3 = exp((-1/tau_s) * DELTAT);
double B3 = (1-A3)*tau_s;
double C3 = G_s/tau_s;
double D3 = 0;
double x1_deltaHR_s    = A3 * x0_deltaHR_s + B3 * SNA_buffer[(int)(T_s/DELTAT)-1];
double deltaHR_s       = C3 * x0_deltaHR_s + D3 * SNA_buffer[(int)(T_s/DELTAT)-1];


double tau_v = data->p_ursino[138+19]; // simplified from eq. 28; units: s

double k_v2 = G_v0 / k_v0; // eq. 27
double G_v = 1 * (1 - exp(-k_v2 * PNA_buffer[1-1])); // eq. 26, EDIT gain value changed to 1 for more accurate description of HR used in
// sys4 = c2d(tf(G_v,[tau_v 1]),DELTAT);
// [A4, B4, C4, D4] = tf2ss(sys4.Numerator{1},sys4.Denominator{1});
double A4 = exp((-1/tau_v) * DELTAT);
double B4 = (1-A4)*tau_v;
double C4 = G_v/tau_v;
double D4 = 0;
double x1_deltaHR_v    = A4 * x0_deltaHR_v + B4 * PNA_buffer[(int)(T_v/DELTAT)-1];
double deltaHR_v       = C4 * x0_deltaHR_v + D4 * PNA_buffer[(int)(T_v/DELTAT)-1];


// // SNA Effector Sites
// again the paper uses inconsistent notation
double T_sigma_lv = T_e_lv;
double T_sigma_rv = T_e_rv;
double T_sigma_V = T_e_V;
double T_sigma_R = T_e_R;

// The following are all described by eq. 29
// Left ventricular contractility
// sys5 = c2d(tf(G_eff_lv,[tau_sigma_lv 1]),DELTAT);
// [A5, B5, C5, D5] = tf2ss(sys5.Numerator{1},sys5.Denominator{1});
double A5 = exp((-1/tau_sigma_lv) * DELTAT);
double B5 = (1-A5)*tau_sigma_lv;
double C5 = G_eff_lv/tau_sigma_lv;
double D5 = 0;
double x1_sigma_lv    = A5 * x0_sigma_lv + B5 * SNA_buffer[(int)(T_sigma_lv/DELTAT)-1];
double sigma_lv       = C5 * x0_sigma_lv + D5 * SNA_buffer[(int)(T_sigma_lv/DELTAT)-1];

// Right ventricular contractility
// sys6 = c2d(tf(G_eff_rv,[tau_sigma_rv 1]),DELTAT);
// [A6, B6, C6, D6] = tf2ss(sys6.Numerator{1},sys6.Denominator{1});
double A6 = exp((-1/tau_sigma_rv) * DELTAT);
double B6 = (1-A6)*tau_sigma_rv;
double C6 = G_eff_rv/tau_sigma_rv;
double D6 = 0;
double x1_sigma_rv    = A6 * x0_sigma_rv + B6 * SNA_buffer[(int)(T_sigma_rv/DELTAT)-1];
double sigma_rv       = C6 * x0_sigma_rv + D6 * SNA_buffer[(int)(T_sigma_rv/DELTAT)-1];

// Venous tone
//double  sys7 = c2d(tf(G_eff_V,[tau_sigma_V 1]),DELTAT);
// [A7, B7, C7,double  D7] = tf2ss(sys7.Numerator{1},sys7.Denominator{1});
double A7 = exp((-1/tau_sigma_V) * DELTAT);
double B7 = (1-A7)*tau_sigma_V;
double C7 = G_eff_V/tau_sigma_V;
double D7 = 0;
double x1_sigma_V    = A7 * x0_sigma_V + B7 * SNA_buffer[(int)(T_sigma_V/DELTAT)-1];
double sigma_V       = C7 * x0_sigma_V + D7 * SNA_buffer[(int)(T_sigma_V/DELTAT)-1];

// Arterial resistance
//double  sys8 = c2d(tf(G_eff_R,[tau_sigma_R 1]),DELTAT);
// [A8, B8, C8,double  D8] = tf2ss(sys8.Numerator{1},sys8.Denominator{1});
double A8 = exp((-1/tau_sigma_R) * DELTAT);
double B8 = (1-A8)*tau_sigma_R;
double C8 = G_eff_R/tau_sigma_R;
double D8 = 0;
double x1_sigma_R    = A8 * x0_sigma_R + B8 * SNA_buffer[(int)(T_sigma_R/DELTAT)-1];
double sigma_R       = C8 * x0_sigma_R + D8 * SNA_buffer[(int)(T_sigma_R/DELTAT)-1];


data->p_ursino[141 + 19]  = x1_deltaHR_s;
data->p_ursino[142 + 19]  = x1_deltaHR_v;
data->p_ursino[143 + 19]  = x1_sigma_lv;
data->p_ursino[144 + 19]  = x1_sigma_rv;
data->p_ursino[145 + 19]  = x1_sigma_V;
data->p_ursino[146 + 19]  = x1_sigma_R;

data->p_ursino[92 + 19]   = deltaHR_s;
data->p_ursino[93 + 19]   = deltaHR_v;
data->p_ursino[94 + 19]   = sigma_lv;
data->p_ursino[95 + 19]   = sigma_rv;
data->p_ursino[96 + 19]   = (data->p_ursino[147 + 19] - sigma_V)/DELTAT;
data->p_ursino[147 + 19]  = sigma_V;
data->p_ursino[97 + 19]   = sigma_R;
