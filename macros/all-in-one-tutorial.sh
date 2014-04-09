#! /bin/bash
#
#  Before attempting to run this, you need to edit gen_input_file.c and set the location
#  of your reduced trees.  You should also put in any selection modifications.
#  You can use unskimmed, skimmed, or skimmed and slimmed reduced trees as input.
#
#  This assumes you are running tcsh for your shell and have the path to root set up already.
#  Source this file from the parent directory (hbb_stats).  That is, do this
#
#    source macros/all-in-one-tutorial
#
#  The first steps (gen_input_file, build_hbb_workspace, and fitqual_plots) should run in
#  one or two minutes.  The limit and significance toys will take around 10 minutes.
#
#  Sample output and log files can be found here, if you want to compare yours to mine.
#
#  http://physics.ucr.edu/~owen/hbb-stats-demo-output/
#
#  Owen
#  July 16, 2013
#


mkdir -p logfiles-demo
mkdir -p outputfiles-demo


#-- generate a text input file with a fake dataset that has no signal with the sig250 model for the signal.
root -b -q 'gen_input_file3.c+("outputfiles-demo/input-metsig-nosig-250-4metbin.txt",250,0.,4,0.,"METsig")' 2>&1 | tee logfiles-demo/gen-input-metsig-nosig-250-4metbin.log
#-- generate a text input file with a fake dataset that includes some signal at the predicted cross sec with the sig250 model for the signal.
root -b -q 'gen_input_file3.c+("outputfiles-demo/input-metsig-withsig-250-4metbin.txt",250,1.,4,0.,"METsig")' 2>&1 | tee logfiles-demo/gen-input-metsig-withsig-250-4metbin.log



#-- build the workspace root files.
root -b -q 'build_hbb_workspace1.c+("outputfiles-demo/input-metsig-nosig-250-4metbin.txt","outputfiles-demo/ws-metsig-nosig-250-4metbin.root")' 2>&1 | tee logfiles-demo/build-ws-metsig-nosig-250-4metbin.log
root -b -q 'build_hbb_workspace1.c+("outputfiles-demo/input-metsig-withsig-250-4metbin.txt","outputfiles-demo/ws-metsig-withsig-250-4metbin.root")' 2>&1 | tee logfiles-demo/build-ws-metsig-withsig-250-4metbin.log





#-- Run the fit quality plots.  This is just to do a single fit with the observables set to the
#   MC values.  You can skip this step if you are only interested in the frequentist limit and/or significance
#   analysis below.
root -b -q 'fitqual_plots.c+("outputfiles-demo/ws-metsig-nosig-250-4metbin.root")' 2>&1 | tee logfiles-demo/fitqual-metsig-nosig-250-4metbin.log
root -b -q 'fitqual_plots.c+("outputfiles-demo/ws-metsig-withsig-250-4metbin.root")' 2>&1 | tee logfiles-demo/fitqual-metsig-withsig-250-4metbin.log




ntoy=2000
minmu=0.30
maxmu=2.00
npoints=10

#-- run frequentist toys for CLs limits.
root -b -q 'rundemo.c("outputfiles-demo/ws-metsig-nosig-250-4metbin.root",0,$maxmu,$minmu,$npoints,$ntoy)' 2>&1 | tee logfiles-demo/limit-toys-met-nosig-250-4metbin.log


#-- run the toys-based significance analysis.
root -b -q 'rundemosignif.c("outputfiles-demo/ws-metsig-withsig-250-4metbin.root",0,5000)  ' 2>&1 | tee logfiles-demo/signif-met-withsig-250-4metbin.log


