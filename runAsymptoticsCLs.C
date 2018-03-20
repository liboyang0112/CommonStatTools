#include <TROOT.h>
#include "AsymptoticsCLsRunner.h"

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)
R__LOAD_LIBRARY(AsymptoticsCLsRunner.C+)

void runAsymptoticsCLs(const char *inputFile, const char *workspaceName, const char *modelConfigName,
                       const char *dataName, TString paramName, Float_t paramValue, TString workspaceTag,
                       TString outputFolder, Bool_t keepDataBlind, Float_t CL = 0.95, const char *asimovDataName = 0,
                       Bool_t doInjection = kFALSE, Float_t muInjection = 1)
{
   gROOT->ProcessLine(".L Minimization.C+");
   gROOT->ProcessLine(".L AsimovDataMaking.C+");
   gROOT->ProcessLine(".L AsymptoticsCLsRunner.C+");

   EXOSTATS::AsymptoticsCLsRunner limitRunner;
   limitRunner.setBlind(keepDataBlind);
   limitRunner.setInjection(doInjection);
   limitRunner.setInjectionStrength(muInjection);

   limitRunner.run(inputFile, workspaceName, modelConfigName, dataName, paramName, paramValue, workspaceTag,
                   outputFolder, CL, asimovDataName);
}
