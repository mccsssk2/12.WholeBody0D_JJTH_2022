#!/bin/bash
#
# script to run job array for HDT.
# There are always 4 cases: control, HD, control+cool, control+HD+cool.
#
#
#
cd randomPars
make ra
cp ra ../
cd ..
cd model
make hdt
cp hdt ../
cd ..
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
		FILE="saFig.sh"
# the EOM has to be first character on line, else file is assumed not to be finished. this comment does go into the file, so moving it elsewhere.
/bin/cat <<EOM >$FILE
#!/bin/bash
#SBATCH --account=def-kharches
#SBATCH --mail-user=jjosep56@uwo.ca		# Email
#SBATCH --mail-type=ALL								# Email me notifications
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48 			# Graham has 32 cores, Cedar has 48 cores to a node.
#SBATCH --mem-per-cpu=1024M      		# memory; default unit is megabytes
#SBATCH --time=00-02:59          			# time (DD-HH:MM)
parallel -j 96 < job_file
EOM
sbatch $FILE		
#		parallel -j 4 < job_file > stderr 2>&1 &
		cd ../..
	done
done
done

