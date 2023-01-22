#!/bin/bash
#
# script to run job array for HDT.
# There are always 4 cases: control, HD, control+cool, control+HD+cool.
#
#
#
cd randomPars
make veryclean
make ra
mv ra ..
cd ..
cd model
make veryclean
make hdt
mv hdt ..
cd ..
#
#
# shortOrLong is 1000 seconds with value 0, and 6 hours with value 1.
shortOrLong=0
for run in 0 # 0 = no renal failure, 1 = homogenous, 2 = heterogenous one kid, 3 = heterogenous both
do
mkdir rundir${run}
for temperature in 35.5 37.5
do
	for dialysisonoff in 0 1
	do
		mkdir 	 	rundir${run}/dir_${dialysisonoff}_${temperature}
		cp hdt 		rundir${run}/dir_${dialysisonoff}_${temperature}
		cp ra	 	rundir${run}/dir_${dialysisonoff}_${temperature}
		cd 			rundir${run}/dir_${dialysisonoff}_${temperature}
		./ra 			${shortOrLong} ${temperature} ${dialysisonoff} ${run}
		parallel -j 4 < job_file > stderr 2>&1 &
		cd ../..
	done
done
done

