#!/bin/bash
#
#
mkdir mainFig1
mkdir mainFig1/dir_0_35.5/
mkdir mainFig1/dir_1_35.5/
mkdir mainFig1/dir_0_37.5/
mkdir mainFig1/dir_1_37.5/
# copy instance 0.
cp dir_0_35.5/*00000* mainFig1/dir_0_35.5/
cp dir_0_37.5/*00000* mainFig1/dir_0_37.5/
cp dir_1_35.5/*00000* mainFig1/dir_1_35.5/
cp dir_1_37.5/*00000* mainFig1/dir_1_37.5/
#
# copy 100 instances, so that TH can pick 20 he likes.
cp dir_0_35.5/output001* mainFig1/dir_0_35.5/
cp dir_0_37.5/output001* mainFig1/dir_0_37.5/
cp dir_1_35.5/output001* mainFig1/dir_1_35.5/
cp dir_1_37.5/output001* mainFig1/dir_1_37.5/
#
#
# organize cardiacoutput.
cd dir_0_35.5/
cat cardiac_output* >> allcardiacoutput.dat
cd ..
cd dir_1_35.5/
cat cardiac_output* >> allcardiacoutput.dat
cd ..
cd dir_0_37.5/
cat cardiac_output* >> allcardiacoutput.dat
cd ..
cd dir_1_37.5/
cat cardiac_output* >> allcardiacoutput.dat
cd ..
#
#
cp dir_0_35.5/allcardiacoutput.dat mainFig1/dir_0_35.5/
cp dir_1_35.5/allcardiacoutput.dat mainFig1/dir_1_35.5/
cp dir_0_37.5/allcardiacoutput.dat mainFig1/dir_0_37.5/
cp dir_1_37.5/allcardiacoutput.dat mainFig1/dir_1_37.5/
