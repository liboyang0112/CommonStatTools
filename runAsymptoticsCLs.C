#include <TROOT.h>
#include "AsymptoticsCLsRunner.h"

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)
R__LOAD_LIBRARY(AsymptoticsCLsRunner.C+g)

void runAsymptoticsCLs(const char *inputFile, const char *workspaceName, const char *modelConfigName,
                       const char *dataName, TString paramName, Float_t paramValue, TString workspaceTag,
                       TString outputFolder, Bool_t keepDataBlind, Float_t CL = 0.95,
                       const char *asimovDataName = "asimovData_0", Bool_t doInjection = kFALSE,
                       Float_t muInjection = 1, Int_t debugLevel = 2)
{
   EXOSTATS::AsymptoticsCLsRunner limitRunner;
   limitRunner.setBlind(keepDataBlind);
   limitRunner.setInjection(doInjection);
   limitRunner.setInjectionStrength(muInjection);
   limitRunner.setDebugLevel(debugLevel);

   limitRunner.run(inputFile, workspaceName, modelConfigName, dataName, paramName, paramValue, workspaceTag,
                   outputFolder, CL, asimovDataName);
}
