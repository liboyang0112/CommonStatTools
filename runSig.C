/// \file
/// Macro to run p0 calculation using profile likelihood and asymptotic formulae

/*
Author: Aaron Armbruster
Date:   2012-06-01
Email:  armbrusa@umich.edu
Description:

Compute statistical significance with profile likelihood test stat.
Option for uncapped test stat is added (doUncap), as well as an option
to choose which mu value to profile observed data at before generating expected

*/

#include "runSig.h"

#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>
#include <RooDataSet.h>
#include <RooMinimizerFcn.h>
#include <RooNLLVar.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooRealSumPdf.h>
#include <Math/MinimizerOptions.h>
#include <RooMinimizer.h>

#include <TStopwatch.h>

#include "Minimization.h"
#include "AsimovDataMaking.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>

#include <iostream>
#include <iomanip>
#include <sstream>

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)

using namespace std;
using namespace RooFit;
using namespace RooStats;

/// Runs the p0 calculation using profile likelihood and asymptotic formulae
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
/// \param[in] asimovDataName name of the Asimov dataset to be used; if empty, will be created
/// \param[in] conditionalSnapshot name of the conditional snapshot
/// \param[in] nominalSnapshot name of the nominal nuisance parameter snapshot
/// \param[in] doInjection compute the limit after signal injection
/// \param[in] muInjection value of the parameter of interest to be used for injection
/// \param[in] debugLevel (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
///
/// This function takes an input workspace and computes the p-value of the background-only hypothesis.
///
/// The function creates, in the folder outputFolder/asymptotics/, a TFile named e.g. workspaceTag_[BLIND]_p0.root,
/// which contains a TTree. This output TTree, named \c p0, contains a single entry, with standard branches
/// with p-value calculation results, and a special branch called \c paramName which contains the value specified by \c paramValue.
/// For example:
/// \code
///   root[0] p0->Print()
///   mass      : mass/F
///   obs_sig   : obs_sig/F
///   obs_pval  : obs_pval/F
///   med_sig   : med_sig/F
///   med_pval  : med_pval/F
///   inj_sig   : inj_sig/F
///   inj_pval  : inj_pval/F
///   med_q0    : med_q0/F
///   inj_q0    : inj_q0/F
/// \endcode
///
///
/// Tip: the best way to profit from the TTree structure when one needs to run over different workspaces
/// representing similar physics processes (e.g. scanning the resonance mass of a given new physics model)
/// is to \c hadd all output files, so that the resulting TTree contains one entry per mass point.
void runSig(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
            TString paramName, Float_t paramValue, TString workspaceTag, TString outputFolder, Bool_t keepDataBlind,
            const char *asimovDataName, const char *conditionalSnapshot,
            const char *nominalSnapshot, Bool_t doInjection, Float_t muInjection, Int_t debugLevel)

{
   const UInt_t nCPU             = 1; // number of CPUs to be used in the fit
   const double       mu_profile_value = 1; // mu value to profile the obs data at before generating the expected
   const bool         doConditional    = 1 && !keepDataBlind; // do conditional expected data
   const bool         doUncap          = 1;                   // uncap p0
   const bool         doInj            = doInjection; // setup the poi for injection study (zero is faster if you're not)
   const double       muInj            = muInjection; // injection value for mu (if injection is requested)
   const bool         doObs            = 1 && !keepDataBlind; // compute median significance
   const bool         doMedian         = 1;                   // compute observed significance

   TStopwatch timer;
   timer.Start();

   TFile         f(inputFile);
   RooWorkspace *ws = (RooWorkspace *)f.Get(workspaceName);
   if (!ws) {
      cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
      return;
   }
   RooFIter   rfiter = ws->components().fwdIterator();
   RooAbsArg *arg;
   while ((arg = rfiter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
         arg->setAttribute("BinnedLikelihood");
      }
   }

   ModelConfig *mc = (ModelConfig *)ws->obj(modelConfigName);
   if (!mc) {
      cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
      return;
   }
   RooDataSet *data = (RooDataSet *)ws->data(dataName);
   if (!data) {
      cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
      return;
   }

   mc->GetNuisanceParameters()->Print("v");

   // RooNLLVar::SetIgnoreZeroEntries(1);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
   cout << "Setting max function calls" << endl;
   // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
   // RooMinimizer::SetMaxFunctionCalls(10000);

   ws->loadSnapshot("conditionalNuis_0");
   RooArgSet nuis(*mc->GetNuisanceParameters());

   RooRealVar *mu = (RooRealVar *)mc->GetParametersOfInterest()->first();

   RooAbsPdf *  pdf     = mc->GetPdf();

   string     condSnapshot(conditionalSnapshot);
   RooArgSet  nuis_tmp2 = *mc->GetNuisanceParameters();
   RooNLLVar *obs_nll   = doObs ? (RooNLLVar *)pdf->createNLL(*data, Constrain(nuis_tmp2)) : NULL;

   RooDataSet *asimovData1 = (RooDataSet *)ws->data(asimovDataName);
   if (!asimovData1) {
      cout << "Asimov data doesn't exist! Please, allow me to build one for you..." << endl;
      string mu_str, mu_prof_str;
      asimovData1  = EXOSTATS::makeAsimovData(ws, mc->GetName(), doConditional, obs_nll, 1, &mu_str, &mu_prof_str,
                                              mu_profile_value, true);
      condSnapshot = "conditionalGlobs" + mu_prof_str;
   }

   if (!doUncap)
      mu->setRange(0, 40);
   else
      mu->setRange(-40, 40);

   //   RooAbsPdf* pdf = mc->GetPdf();
   pdf = mc->GetPdf();

   RooArgSet  nuis_tmp1 = *mc->GetNuisanceParameters();
   RooNLLVar *asimov_nll =
      (RooNLLVar *)pdf->createNLL(*asimovData1, Constrain(nuis_tmp1), Offset(1), Optimize(2), NumCPU(nCPU, 1));

   // do asimov
   mu->setVal(1);
   mu->setConstant(0);
   if (!doInj) mu->setConstant(1);

   int    status, sign;
   double med_sig = 0, obs_sig = 0, inj_sig = 0, asimov_q0 = 0, obs_q0 = 0, inj_q0 = 0;

   if (doMedian) {
      ws->loadSnapshot(condSnapshot.c_str());
      if (doInj)
         ws->loadSnapshot("conditionalNuis_inj");
      else
         ws->loadSnapshot("conditionalNuis_1");
      mc->GetGlobalObservables()->Print("v");
      mu->setVal(0);
      mu->setConstant(1);
      status = EXOSTATS::minimize(asimov_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(asimov_nll);
         if (status >= 0) cout << "Success!" << endl;
      }
      double asimov_nll_cond = asimov_nll->getVal();

      mu->setVal(1);
      if (doInj)
         ws->loadSnapshot("conditionalNuis_inj");
      else
         ws->loadSnapshot("conditionalNuis_1");
      if (doInj) mu->setConstant(0);
      status = EXOSTATS::minimize(asimov_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(asimov_nll);
         if (status >= 0) cout << "Success!" << endl;
      }

      double asimov_nll_min = asimov_nll->getVal();
      asimov_q0             = 2 * (asimov_nll_cond - asimov_nll_min);
      if (doUncap && mu->getVal() < 0) asimov_q0 = -asimov_q0;

      sign    = int(asimov_q0 != 0 ? asimov_q0 / fabs(asimov_q0) : 0);
      med_sig = sign * sqrt(fabs(asimov_q0));

      ws->loadSnapshot(nominalSnapshot);
   }

   if (doObs) {

      ws->loadSnapshot("conditionalNuis_0");
      mu->setVal(0);
      mu->setConstant(1);
      status = EXOSTATS::minimize(obs_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(obs_nll);
         if (status >= 0) cout << "Success!" << endl;
      }
      double obs_nll_cond = obs_nll->getVal();

      // ws->loadSnapshot("ucmles");
      mu->setConstant(0);
      status = EXOSTATS::minimize(obs_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(obs_nll);
         if (status >= 0) cout << "Success!" << endl;
      }

      double obs_nll_min = obs_nll->getVal();

      obs_q0 = 2 * (obs_nll_cond - obs_nll_min);
      if (doUncap && mu->getVal() < 0) obs_q0 = -obs_q0;

      sign = int(obs_q0 == 0 ? 0 : obs_q0 / fabs(obs_q0));
      if (!doUncap && ((obs_q0 < 0 && obs_q0 > -0.1) || mu->getVal() < 0.001))
         obs_sig = 0;
      else
         obs_sig = sign * sqrt(fabs(obs_q0));
   }

   if (doInj) {

      string mu_str, mu_prof_str;
      double mu_inj = 1.;
      if (ws->var("ATLAS_norm_muInjection")) {
         mu_inj = ws->var("ATLAS_norm_muInjection")->getVal();
      } else {
         mu_inj = muInj; // for the mass point at the inj
      }
      RooDataSet *injData1 =
         EXOSTATS::makeAsimovData(ws, mc->GetName(), doConditional, obs_nll, 0, &mu_str, &mu_prof_str, 1, true, mu_inj, debugLevel);
      string globObsSnapName = "conditionalGlobs" + mu_prof_str;
      ws->loadSnapshot(globObsSnapName.c_str());
      RooNLLVar *inj_nll =
         (RooNLLVar *)pdf->createNLL(*injData1, Constrain(nuis_tmp2), Offset(1), Optimize(2), NumCPU(nCPU, 1));

      ws->loadSnapshot("conditionalNuis_0");
      mu->setVal(0);
      mu->setConstant(1);
      status = EXOSTATS::minimize(inj_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(inj_nll);
         if (status >= 0) cout << "Success!" << endl;
      }
      double inj_nll_cond = inj_nll->getVal();

      // ws->loadSnapshot("ucmles");
      mu->setConstant(0);
      status = EXOSTATS::minimize(inj_nll);
      if (status < 0) {
         cout << "Retrying with conditional snapshot at mu=1" << endl;
         ws->loadSnapshot("conditionalNuis_0");
         status = EXOSTATS::minimize(inj_nll);
         if (status >= 0) cout << "Success!" << endl;
      }

      double inj_nll_min = inj_nll->getVal();

      inj_q0 = 2 * (inj_nll_cond - inj_nll_min);
      if (doUncap && mu->getVal() < 0) inj_q0 = -inj_q0;

      sign = int(inj_q0 == 0 ? 0 : inj_q0 / fabs(inj_q0));
      if (!doUncap && ((inj_q0 < 0 && inj_q0 > -0.1) || mu->getVal() < 0.001))
         inj_sig = 0;
      else
         inj_sig = sign * sqrt(fabs(inj_q0));
   }

   f.Close();

   // prepare TTree
   TTree * tree          = new TTree("p0", "runSig");
   Float_t tree_point    = paramValue;
   Float_t tree_obs_sig  = obs_sig;
   Float_t tree_obs_pval = (1 - TMath::Erf(obs_sig / sqrt(2.))) / 2.;
   Float_t tree_med_sig  = med_sig;
   Float_t tree_med_pval = (1 - TMath::Erf(med_sig / sqrt(2.))) / 2.;
   Float_t tree_inj_sig  = inj_sig;
   Float_t tree_inj_pval = (1 - TMath::Erf(inj_sig / sqrt(2.))) / 2.;
   Float_t tree_med_q0   = asimov_q0;
   Float_t tree_inj_q0   = inj_q0;
   tree->Branch(paramName, &tree_point);
   tree->Branch("obs_sig", &tree_obs_sig);
   tree->Branch("obs_pval", &tree_obs_pval);
   tree->Branch("med_sig", &tree_med_sig);
   tree->Branch("med_pval", &tree_med_pval);
   tree->Branch("inj_sig", &tree_inj_sig);
   tree->Branch("inj_pval", &tree_inj_pval);
   tree->Branch("med_q0", &tree_med_q0);
   tree->Branch("inj_q0", &tree_inj_q0);
   tree->Fill();

   // print out
   if (!keepDataBlind) {
      cout << "Observed significance: " << tree_obs_sig << endl;
      cout << "Observed pValue: " << tree_obs_pval << endl;
   }
   if (med_sig) {
      cout << "Median test stat val: " << tree_med_q0 << endl;
      cout << "Median significance:   " << tree_med_sig << endl;
      cout << "Median pValue: " << tree_med_pval << endl;
   }
   if (inj_sig) {
      cout << "Injected test stat val: " << tree_inj_q0 << endl;
      cout << "Injected significance:   " << tree_inj_sig << endl;
      cout << "Injected pValue: " << tree_inj_pval << endl;
   }

   // save
   const TString blindness     = (keepDataBlind) ? "_BLIND" : "";
   const TString fullOutFolder = outputFolder + "/asymptotics/";
   system("mkdir -vp " + fullOutFolder);
   const TString outFileName = fullOutFolder + workspaceTag + blindness + "_p0.root";
   TFile         f_out(outFileName, "RECREATE");
   tree->SetDirectory(&f_out);
   f_out.Write();

   timer.Stop();
   timer.Print();
}
