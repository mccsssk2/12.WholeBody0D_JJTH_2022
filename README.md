# Repository.  

Whole body (human) lumped parameter modelling. 2019-2023.

# Background.  

As part of the hemodynamics project (April 2019 to March 2023), two junior developers undertook graduate studies under
my advice. While implementing lumped paramter boundary conditions, which are essentially blood vessel terminal resistances,
impedences, capacitances, a whole body human model described by ordinary differential equations (lumped parameters) 
was implemented based on existing literature. In addition, a suitable ODE description of a dialyser machine representing
a hemodynamic stress was coupled to the human model. The code is setup for sensitivity analysis of relevant parameters and 
hemodynamic variables.  

The outcomes of this part of the project were as listed below. I was senior or co-senior author on most of these
articles and peer reviewed conference proceedings. I was also responsible for finishing of Timothy's and Jermiah's thesis.
They have now progressed as scientific programmers in Canadian universities/companies. In addition to Jermiah and Tim,
I also nurtured experiental learning student (Clara Sun) who is in the process of progressing to medical school in Canada 
or USA.  
* Hunter et al. Hunter, T.J., Joseph, J.J., Anazodo, U., Kharche, S.R., McIntyre, C.W., Goldman, D. (2021). Computational Modelling of the Role of Atrial Fibrillation on Cerebral Blood Perfusion. In: Ennis, D.B., Perotti, L.E., Wang, V.Y. (eds) Functional Imaging and Modeling of the Heart. FIMH 2021. Lecture Notes in Computer Science(), vol 12738. Springer, Cham. https://doi.org/10.1007/978-3-030-78710-3_65 (peer reviewed international conference proceedings, SRK was senior/co-senior author).  
* Joseph, J.J., Lee, TY., Goldman, D., McIntyre, C.W., Kharche, S.R. (2021). The Role of Extra-Coronary Vascular Conditions that Affect Coronary Fractional Flow Reserve Estimation. In: Ennis, D.B., Perotti, L.E., Wang, V.Y. (eds) Functional Imaging and Modeling of the Heart. FIMH 2021. Lecture Notes in Computer Science(), vol 12738. Springer, Cham. https://doi.org/10.1007/978-3-030-78710-3_57 (peer reviewed international conference proceedings, SRK was senior/co-senior author).  
* Joseph, J.J.; Hunter, T.J.; Sun, C.; Goldman, D.; Kharche, S.R.; McIntyre, C.W. Using a Human Circulation Mathematical Model to Simulate the Effects of Hemodialysis and Therapeutic Hypothermia. Appl. Sci. 2022, 12, 307. https://doi.org/10.3390/app12010307 (SRK is co-senior author on this journal article).  
*  Hunter, T.J.; Joseph, J.J.; Anazodo, U.; Kharche, S.R.; McIntyre, C.W.; Goldman, D. Atrial Fibrillation and Anterior Cerebral Artery Absence Reduce Cerebral Perfusion: A De Novo Hemodynamic Model. Appl. Sci. 2022, 12, 1750. https://doi.org/10.3390/app12031750 (SRK is co-senior author on this journal article).  
* Joseph, J.J.; Sun, C.; Lee, T.-Y.; Goldman, D.; Kharche, S.R.; McIntyre, C.W. Structure (Epicardial Stenosis) and Function (Microvascular Dysfunction) That Influence Coronary Fractional Flow Reserve Estimation. Appl. Sci. 2022, 12, 4281. https://doi.org/10.3390/app12094281 (SRK is co-senior author on this journal article).  
* Jermiah Joseph. MSc graduate thesis. https://ir.lib.uwo.ca/cgi/viewcontent.cgi?article=11688&context=etd (SRK was advisor).  
* Timothy Hunter. MSc graduate thesis. https://ir.lib.uwo.ca/cgi/viewcontent.cgi?article=11510&context=etd (SRK was advisor).  

# Dependencies.  

This code in C language requires GNU C compiler and Sundials libraries.

# Install.  

Availability of compiler and Sundials.

# Source description.  

The program solves the human-dialyser circulation system using a higher order Backward Difference Fomulae ODE solver available in Sundials. Further, dynamic local sensivitiy is provided by use of CVODES in Sundials. Global sensitivity is obtained
by first performing model simulations for N instances. To do so, a Latin Hypercube Sampling of paramter space using the Meresenne Twister random number generator is done. The outputs are then provided to an uptaken lhs-prcc code. Figures showing inputs, outputs, and Results are provided in multiple sub-dirs of source/ directory.  

The driver function is ursino.c in the source/model directory. The human artetial network model and dialyser models are implemented in f_ursino.c . Autoregulation is provided in C files that start with BR. Other functions are implemented in helper C and header files. In case sensitivity analysis is to be performed, the source/randomPars/ directory will generate N instances of inputs, which is the input to the model and eventually to the partial rank correlation coefficient (PRCC, uptaken from Kuhl group in Stanford) code in sources/Results/PRCC_Jan/ .  

# Use.  

Editing of the provided makefile for library locations and compiler should generate an executable binary. If SA is desired, generate inputs using the program in randomPars/ and move the files to appropriate locations. A bash shell script may be provided to help.

# Maintainer.  

This code is provided as is to interested users.

# Provenance.  

This code can be further developed. It can be made modular. It can be extended to incorporate several known hemodynamic
mechnisms which may render it useful in experimental research. For example, the predictor of cerebral flow and pressure,
that is inaccessible to experiment, based on large vessel recordings may become feasible. While it has commercialisation potential in the mid-term (5 years), development of the package can eventually be deployed in following as examples:  
* Virtual patient simulation on the effects of disease (dialysis, loss of autoregulation, loss of tone) on hemodynamics.  
* A rapid and efficient way to provide clinical metrics, see Sheffield Insignio activities.  
* Interpretation of routine clinical measurements such as cerebral TCD and velocimetry.  
* complement in vivo experiments for data analysis.  
* Complement epidemiology data science.  

# Acknowledements.

This project was generously funded by CANARIE Inc. Canada (April 2020 to March 2023). We also thank the Kidney Unit in Lawson, London Ontario, Canada, for hosting this project.

# Licence.  

These codes are the collective IP of SR Kharche, Jermiah Joseph, Tim Hunter, grant investigators, Lawson, and Western University; whose approval may be required for use and deployment.  

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

