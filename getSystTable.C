#include "HistFactoryInspector.h"

R__LOAD_LIBRARY(HistFactoryInspector.C+)

void getSystTable(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName, TString workspaceTag, TString outputFolder, TString samples, TString evalRegions, TString fitRegions, Bool_t doYields = kTRUE, Bool_t doImpacts = kTRUE)
{
  const TString rangeName = ""; // do not use

  EXOSTATS::HistFactoryInspector hf;
  hf.setInput(inputFile, workspaceName, modelConfigName, dataName, rangeName);
  hf.setEvalRegions(evalRegions);
  hf.setFitRegions(fitRegions);

  hf.getYields();
}
