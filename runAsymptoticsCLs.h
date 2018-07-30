#ifndef ASYMPTOTICSCLS_H
#define ASYMPTOTICSCLS_H

#include <TROOT.h>
#include "AsymptoticsCLsRunner.h"

void runAsymptoticsCLs(const char *inputFile, const char *workspaceName, const char *modelConfigName,
                       const char *dataName, TString paramName, Float_t paramValue, TString workspaceTag,
                       TString outputFolder, Bool_t keepDataBlind, Float_t CL = 0.95,
                       const char *asimovDataName = "asimovData_0", Bool_t doInjection = kFALSE,
                       Float_t muInjection = 1, Int_t debugLevel = 2);

#endif
