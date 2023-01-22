// *****************************************************************************************************************
// my OUTPUTS
// START OF CARDIAC OUTPUT
// Initialize Variables to compare to in cardiac output
PresToCompare[0] = Ith(y_ursino, 13); FlowsToCompare[0] = data->Qlo; 	// Aorta
PresToCompare[1] = Ith(y_ursino, 17); FlowsToCompare[1] = data->Qpa; 	// Pulmonary
PresToCompare[2] = Ith(y_ursino, 4); 	FlowsToCompare[2] = data->Qupi; // upper
PresToCompare[3] = Ith(y_ursino, 6); 	FlowsToCompare[3] = data->Qsp1; // splanchnic
PresToCompare[4] = Ith(y_ursino, 49); FlowsToCompare[4] = data->QkRi; // R. Kid Inlet
PresToCompare[5] = Ith(y_ursino, 50); FlowsToCompare[5] = data->QkLi; // L. Kid Inlet
PresToCompare[6] = Ith(y_ursino, 7); 	FlowsToCompare[6] = data->Qll1; // Lower Body

if((data->p_ursino[105] >= DELTAT)&& (data->p_ursino[105]<data->p_ursino[11])) {
		cardiac_output = cardiac_output + (data->Qlo) * DELTAT;
		flowpermin[0] += data->Qpa*DELTAT;
		flowpermin[1] += data->Qupi*DELTAT;
		flowpermin[2] += data->QkRi*DELTAT;
		flowpermin[3] += data->QkLi*DELTAT;
		flowpermin[4] += data->Qsp1*DELTAT;
		flowpermin[5] += data->Qll1*DELTAT;

		for(i = 0; i<7; i++){
			if(PresToCompare[i] > sysPressures[i]) 				sysPressures[i] 	= PresToCompare[i];
			else if(PresToCompare[i] < diasPressures[i]) 	diasPressures[i] 	= PresToCompare[i];

			if(FlowsToCompare[i] > maxFlows[i]) 				maxFlows[i] = FlowsToCompare[i];
			else if(FlowsToCompare[i] < minFlows[i]) 		minFlows[i] = FlowsToCompare[i];
		}
}else if((data->p_ursino[105] < DELTAT) && (cardiac_output > 0.0)) {
	str = malloc(32*sizeof(char)); sprintf(str,"measurements%05d.dat", atoi(argv[1]));  cardiac_output_file = fopen(str,"w"); 	free(str);

// inputs: the paramters. NP is 167. p_ursino is not all perturbed. There are 51 that are.
fprintf(cardiac_output_file,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", data->p_ursino[167], G_factor, data->Edias_lv, data->Esys_lv, data->Edias_rv, data->Esys_rv, data->Edias_la, data->Esys_la, data->Edias_ra, data->Esys_ra);
	for(i = 0; i<27; i++) fprintf(cardiac_output_file, "%f\t", data->p_C[i]);
	for(i = 0; i<45; i++) fprintf(cardiac_output_file, "%f\t", data->p_R[i]);

// outputs.
	fprintf(cardiac_output_file, "%f\t%f\t%f\t", HR, cardiac_output*HR, data->p_ursino[167]);
	cardiac_output = 0.0;

	// 11-16) flow per minute to each systemic compartment
	for(int f = 0; f<6; f++){
		fprintf(cardiac_output_file, "%f\t", flowpermin[f]*HR);
		flowpermin[f] = 0.0;
	}

	data->MAP = (sysPressures[0] + 2.0*diasPressures[0])/3.0;
	// Each compartment starts at 17, 21, 25, 29, 33, 37, 41
	for(i = 0; i<7; i++){
		fprintf(cardiac_output_file, "%f\t%f\t", sysPressures[i], diasPressures[i]);
		fprintf(cardiac_output_file, "%f\t%f\t", maxFlows[i], minFlows[i]);

		// Reset variables
		sysPressures[i] 	= diasPressures[i] 	= PresToCompare[i];
		maxFlows[i] 			= minFlows[i] 			= FlowsToCompare[i];
	}

	fprintf(cardiac_output_file, "%f\t%f\t%f\t", data->MAP_avg, data->tau_avg[9], data->tau_avg[18]);
	fprintf(cardiac_output_file, "%f\t%f\t%f\t", data->v1_avg, data->v2_avg, data->v3_avg);

	fprintf(cardiac_output_file, "\n");
	fclose(cardiac_output_file);
}
// END OF CARDIAC OUTPUT

