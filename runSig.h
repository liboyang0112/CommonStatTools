#ifndef RUNSIG_H
#define RUNSIG_H

#include "TString.h"

void runSig(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
            TString paramName, Float_t paramValue, TString workspaceTag, TString outputFolder, Bool_t keepDataBlind,
            const char *asimovDataName = "asimovData_1", const char *conditionalSnapshot = "conditionalGlobs_1",
            const char *nominalSnapshot = "nominalGlobs", Bool_t doInjection = kFALSE, Float_t muInjection = 1, Int_t debugLevel = 2);

#endif
