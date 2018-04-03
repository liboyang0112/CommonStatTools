/// \file
/// Macro to run FitCrossCheckForLimits.C

#include <TROOT.h>

R__LOAD_LIBRARY(Minimization.C +)
R__LOAD_LIBRARY(AsimovDataMaking.C +)
R__LOAD_LIBRARY(FitCrossCheckForLimits.C +)

/// Runs the CLs limit using profile likelihood and asymptotic formulae
///
/// \param[in] algorithm to run
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input ModelConfig
/// \param[in] dataName name of the dataset to fit
/// \param[in] workspaceTag prefix for the output ROOT file
/// \param[in] outputFolder path under which the output ROOT file will be stored; it will be created if it does not
/// exist \param[in] doConditional if true, perform conditional fit to data \param[in] signalStrength value of the
/// parameter of interest for the conditional fit \param[in] numberOfSigmas number of sigmas for the definition of pull
/// \param[in] draw1DResponse
/// \param[in] createPostFitAsimov
/// \param[in] debugLevel (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
///
/// This function takes an input workspace and computes upper limits on the parameter of interest
/// (POI, which usually represents the signal strength) at a given confidence level, using the
/// profile-likelihood test statistics, the CLs technique and asymptotic formulae. Compared to ROOT's
/// default calculation (as implemented in $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C),
/// this calculation is faster, makes full usage of asymptotic formulae and does not rely on the user's
/// guesstimate of the range of the POI to scan.
///
/// The function creates, in the folder outputFolder/asymptotics/, a TFile named e.g. workspaceTag_[BLIND]_CL95.root,
/// which contains a TTree. This output TTree, named \c stats, contains a single entry, with standard branches
/// with limit results, a special branch called \c paramName which contains the value specified by \c paramValue,
/// and a branch per nuisance parameter with the values of each of them after the background-only and
/// signal-plus-background fit.

void runFitCrossCheckForLimits(Algs algorithm, const char *inputFile, const char *workspaceName,
                               const char *modelConfigName, const char *dataName, TString workspaceTag,
                               TString outputFolder, Bool_t doConditional, Float_t signalStrength,
                               Float_t numberOfSigmas, Bool_t draw1DResponse = kFALSE,
                               Bool_t createPostFitAsimov = kFALSE, Int_t debugLevel = 2)
{
   LimitCrossChecker fitcheck;
   fitcheck.setDebugLevel(2);
   fitcheck.FitCrossCheckForLimits(algorithm, signalStrength, numberOfSigmas, doConditional, inputFile, outputFolder,
                                   workspaceName, modelConfigName, dataName, draw1DResponse, createPostFitAsimov);
}
