/// \file
/// Macro to run limit calculation using profile likelihood, CLs and asymptotic formulae

#include <TROOT.h>
#include "AsymptoticsCLsRunner.h"
#include "runAsymptoticsCLs.h"

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)
R__LOAD_LIBRARY(AsymptoticsCLsRunner.C+g)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)

/// Runs the CLs limit using profile likelihood and asymptotic formulae
///
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input ModelConfig
/// \param[in] dataName name of the dataset to fit
/// \param[in] paramName name of the parameter which characterizes this workspace (e.g. "m_Higgs"); it will be added as a branch of the output TTree
/// \param[in] paramValue value of the parameter which characterizes this workspace (e.g. 125.0); it will be stored in the corresponding branch of the output TTree
/// \param[in] workspaceTag prefix for the output ROOT file
/// \param[in] outputFolder path under which the output ROOT file will be stored; it will be created if it does not exist
/// \param[in] keepDataBlind if true, data are not used (blind analysis)
/// \param[in] CL confidence level for the calculation (e.g. 0.95 for 95% CL limits)
/// \param[in] asimovDataName name of the Asimov dataset to be used; if empty, will be created
/// \param[in] doInjection compute the limit after signal injection
/// \param[in] muInjection value of the parameter of interest to be used for injection
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
/// and a branch per nuisance parameter with the values of each of them after the background-only and signal-plus-background fit.
/// For example:
/// \code
///   root[0] stats->Print()
///   mass      : mass/F
///   CLb_med   : CLb_med/F
///   pb_med    : pb_med/F
///   CLs_med   : CLs_med/F
///   CLsplusb_med : CLsplusb_med/F
///   CLb_obs   : CLb_obs/F
///   pb_obs    : pb_obs/F
///   CLs_obs   : CLs_obs/F
///   CLsplusb_obs : CLsplusb_obs/F
///   obs_upperlimit : obs_upperlimit/F
///   inj_upperlimit : inj_upperlimit/F
///   exp_upperlimit : exp_upperlimit/F
///   exp_upperlimit_plus1 : exp_upperlimit_plus1/F
///   exp_upperlimit_plus2 : exp_upperlimit_plus2/F
///   exp_upperlimit_minus1 : exp_upperlimit_minus1/F
///   exp_upperlimit_minus2 : exp_upperlimit_minus2/F
///   fit_status : fit_status/F
///   mu_hat_obs : mu_hat_obs/F
///   mu_hat_exp : mu_hat_exp/F
///   param_alpha_G2_hat : param_alpha_G2_hat/F
///   param_alpha_bkgER_hat : param_alpha_bkgER_hat/F
///   param_alpha_bkg_NormSys_bkgCryo_hat :
///   param_alpha_bkg_NormSys_bkgInt_hat :
///   param_alpha_bkg_NormSys_bkgPMT_hat :
///   param_alpha_sigQy_hat : param_alpha_sigQy_hat/F
///   param_alpha_G2_med : param_alpha_G2_med/F
///   param_alpha_bkgER_med : param_alpha_bkgER_med/F
///   param_alpha_bkg_NormSys_bkgCryo_med :
///   param_alpha_bkg_NormSys_bkgInt_med :
///   param_alpha_bkg_NormSys_bkgPMT_med :
///   param_alpha_sigQy_med : param_alpha_sigQy_med/F
/// \endcode
///
///
/// Tip: the best way to profit from the TTree structure when one needs to run over different workspaces
/// representing similar physics processes (e.g. scanning the resonance mass of a given new physics model)
/// is to \c hadd all output files, so that the resulting TTree contains one entry per mass point. To do so,
/// make sure you give something meaningful as a workspaceTag parameter (e.g. "mDM_125gev")...

void runAsymptoticsCLs(const char *inputFile, const char *workspaceName, const char *modelConfigName,
                       const char *dataName, TString paramName, Float_t paramValue, TString workspaceTag,
                       TString outputFolder, Bool_t keepDataBlind, Float_t CL,
                       const char *asimovDataName, Bool_t doInjection,
                       Float_t muInjection, Int_t debugLevel)
{
   EXOSTATS::AsymptoticsCLsRunner limitRunner;
   limitRunner.setBlind(keepDataBlind);
   limitRunner.setInjection(doInjection);
   limitRunner.setInjectionStrength(muInjection);
   limitRunner.setDebugLevel(debugLevel);

   limitRunner.run(inputFile, workspaceName, modelConfigName, dataName, paramName, paramValue, workspaceTag,
                   outputFolder, CL, asimovDataName);
}
