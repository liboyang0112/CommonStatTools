/// \file
/// Macro to extract yield and impact tables from a HistFactory workspace
#include "HistFactoryInspector.h"

R__LOAD_LIBRARY(RooExpandedFitResult.C+)
R__LOAD_LIBRARY(HistFactoryInspector.C+)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)

/// Prints yield tables and nuisance parameter impact tables from a HistFactory workspace
///
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input ModelConfig
/// \param[in] dataName name of the dataset to fit
/// \param[in] workspaceTag prefix for the output files
/// \param[in] outputFolder path under which the output files will be stored; it will be created if it does not exist
/// \param[in] evalRegions comma-separated list of regions where yields/impacts should be evaluated
/// \param[in] fitRegions comma-separated list of regions where the fit must be performed
/// \param[in] doYields if \c kTRUE, calculate pre- and post-fit yields for all samples in the specified regions
/// \param[in] doImpacts if \c kTRUE, calculate nuisance parameter impacts on the sum of the yields of \c samplesForImpacts in the specified regions
/// \param[in] samplesForImpact comma-separated list of samples to evaluate nuisance parameter impacts for
/// \param[in] debugLevel (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
///
/// This function takes any HistFactory workspace and computes:
///   - yields of all samples before and after the fit to specific regions
///   - the impact of all nuisance parameters on the sum of the yields of specified samples after the fit to specific regions
/// 
/// Keep in mind that the likelihood will be fitted only in the regions specified by \c fitRegions, and that outputs
/// will be calculated only in the regions specified by \c evalRegions.
/// 
/// The function, as is, simply prints outputs on screen. If you want to parse them more conveniently, for example
/// inserting them in CSV tables, python dictionaries of LaTeX files, please see EXOSTATS::HistFactoryInspector::getYields()
/// and EXOSTATS::HistFactoryInspector::getYields() for details on how to read the output of these functions.
void getHFtables(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName, TString workspaceTag, TString outputFolder, TString evalRegions, TString fitRegions, Bool_t doYields = kTRUE, Bool_t doImpacts = kFALSE, TString samplesForImpact = "", Int_t debugLevel = 2)
{
  const TString rangeName = ""; // do not use

  EXOSTATS::HistFactoryInspector hf;
  hf.setDebugLevel(debugLevel);
  hf.setInput(inputFile, workspaceName, modelConfigName, dataName, rangeName);
  hf.setEvalRegions(evalRegions);
  hf.setFitRegions(fitRegions);

  if (doYields) hf.getYields();
  if (doImpacts) hf.getImpacts(samplesForImpact);
}
