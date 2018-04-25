/// \file
/// Macro to run the LimitCrossChecker

#include <TROOT.h>

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)
R__LOAD_LIBRARY(FitCrossCheckForLimits.C+)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)

/// Runs the LimitCrossChecker class over a workspace
///
/// \param[in] algorithm to run
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input ModelConfig
/// \param[in] dataName name of the dataset to fit
/// \param[in] outputFolder path under which the output will be stored
/// \param[in] doConditional if true, perform conditional fit to data
/// \param[in] signalStrength value of the parameter of interest for the conditional fit
/// \param[in] numberOfSigmas number of sigmas for the up/down variations in plots
/// \param[in] draw1DResponse plot the 1D response for each nuisance parameter
/// \param[in] createPostFitAsimov
/// \param[in] debugLevel (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
///
/// Implemented algorithms:
///
/// - PlotHistosBeforeFit,
/// - PlotMorphingControlPlots,
/// - PlotHistosAfterFitEachSubChannel,
/// - PlotHistosAfterFitGlobal,
/// - PlotsNuisanceParametersVSmu,
/// - PlotsStatisticalTest,
/// - FitToAsimov
///
/// For example, you might want to run:
///     
/// - \c PlotHistosAfterFitGlobal with mu=1 and \c doConditional set to \c false (only if not blind)
/// - \c PlotHistosAfterFitGlobal with mu=1 and \c doConditional set to \c true
/// - \c FitToAsimov and \c doConditional set to \c false
/// - \c FitToAsimov and \c doConditional set to \c true


void runFitCrossCheckForLimits(Algs algorithm, const char *inputFile, const char *workspaceName,
                               const char *modelConfigName, const char *dataName,
                               TString outputFolder, Bool_t doConditional, Float_t signalStrength,
                               Float_t numberOfSigmas, Bool_t draw1DResponse = kFALSE,
                               Bool_t createPostFitAsimov = kFALSE, Int_t debugLevel = 2)
{
   LimitCrossChecker fitcheck;
   fitcheck.setDebugLevel(2);

   outputFolder = outputFolder + "/crosschecks";
   fitcheck.run(algorithm, signalStrength, numberOfSigmas, doConditional, inputFile, outputFolder,
                                   workspaceName, modelConfigName, dataName, draw1DResponse, createPostFitAsimov);
}
