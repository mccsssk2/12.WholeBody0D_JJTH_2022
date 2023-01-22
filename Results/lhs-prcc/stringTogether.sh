#!/bin/bash
#
#
rm   figSAFrontiers.png
mogrify -trim *.png
convert +append Frontiers_PRCC2.png Frontiers_PRCC1.png top.png
convert +append Frontiers_PRCC4.png Frontiers_PRCC5.png mid.png
convert +append Frontiers_PRCC7.png Frontiers_PRCC8.png bot.png
convert -append top.png mid.png bot.png figSAFrontiers.png
