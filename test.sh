#!/bin/bash

WORKSPACEFILE=/tmp/exo/results/example_DS_04p000_combined_DS_04p000_model.root
WORKSPACENAME=combined
MODELCONFIGNAME=ModelConfig
DATASETNAME=obsData
PARAMNAME=mass
PARAMVALUE=123
WORKSPACETAG=test
OUTPUTFOLDER=validation
DOBLIND=0
CL=0.90

root -b -q runAsymptoticsCLs.C\(\"$WORKSPACEFILE\",\"$WORKSPACENAME\",\"$MODELCONFIGNAME\",\"$DATASETNAME\",\"$PARAMNAME\",$PARAMVALUE,\"$WORKSPACETAG\",\"$OUTPUTFOLDER\",$DOBLIND,$CL\)
root -b -q runSig.C\(\"$WORKSPACEFILE\",\"$WORKSPACENAME\",\"$MODELCONFIGNAME\",\"$DATASETNAME\",\"$PARAMNAME\",$PARAMVALUE,\"$WORKSPACETAG\",\"$OUTPUTFOLDER\",$DOBLIND\)
root -b -q getGlobalP0.C\(3.4,2\)
python compareHistos.py -o $WORKSPACEFILE -n $WORKSPACEFILE
python obtainHistosFromWS.py -i $WORKSPACEFILE -o validation_histos.root
#root -b -q getCorrMatrix.C\(\)
