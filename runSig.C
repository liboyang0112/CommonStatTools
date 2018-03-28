/*
Author: Aaron Armbruster
Date:   2012-06-01
Email:  armbrusa@umich.edu
Description:

Compute statistical significance with profile likelihood test stat.
Option for uncapped test stat is added (doUncap), as well as an option
to choose which mu value to profile observed data at before generating expected

*/

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooMinimizerFcn.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "Math/MinimizerOptions.h"
#include "RooMinimizer.h"

#include "TH1D.h"

#include "TStopwatch.h"

#include "Minimization.h"
#include "AsimovDataMaking.h"

//#include "macros/makeData.C"

#include "TFile.h"

#include <iostream>
#include <iomanip>
#include <sstream>

R__LOAD_LIBRARY(Minimization.C+)
R__LOAD_LIBRARY(AsimovDataMaking.C+)

using namespace std;
using namespace RooFit;
using namespace RooStats;

void runSig(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
            TString paramName, Float_t paramValue, TString workspaceTag, TString outputFolder, Bool_t keepDataBlind,
            const char *asimovDataName = "asimovData_1", const char *conditionalSnapshot = "conditionalGlobs_1",
            const char *nominalSnapshot = "nominalGlobs")

{
   const UInt_t nCPU             = 3; // number of CPUs to be used in the fit
   double       mu_profile_value = 1; // mu value to profile the obs data at before generating the expected
   bool         doConditional    = 1 && !keepDataBlind; // do conditional expected data
   bool         doUncap          = 1;                   // uncap p0
   bool         doInj            = 0; // setup the poi for injection study (zero is faster if you're not)
   bool         doObs            = 1 && !keepDataBlind; // compute median significance
   bool         doMedian         = 1;                   // compute observed significance

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

   const double mu_init = mu->getVal();
   RooAbsPdf *  pdf     = mc->GetPdf();

   string     condSnapshot(conditionalSnapshot);
   RooArgSet  nuis_tmp2 = *mc->GetNuisanceParameters();
   RooNLLVar *obs_nll   = doObs ? (RooNLLVar *)pdf->createNLL(*data, Constrain(nuis_tmp2)) : NULL;

   RooDataSet *asimovData1 = (RooDataSet *)ws->data(asimovDataName);
   RooRealVar *emb         = (RooRealVar *)mc->GetNuisanceParameters()->find("ATLAS_EMB");
   if (!asimovData1 || (string(inputFile).find("ic10") != string::npos && emb)) {
      if (emb) emb->setVal(0.7);
      cout << "Asimov data doesn't exist! Please, allow me to build one for you..." << endl;
      string mu_str, mu_prof_str;
      asimovData1  = EXOSTATS::makeAsimovData(ws, mc->GetName(), doConditional, obs_nll, 1, &mu_str, &mu_prof_str,
                                             mu_profile_value, true);
      condSnapshot = "conditionalGlobs" + mu_prof_str;

      // makeAsimovData(mc, true, ws, mc->GetPdf(), data, 0);
      // ws->Print();
      // asimovData1 = (RooDataSet*)ws->data("asimovData_1");
   }

   if (!doUncap)
      mu->setRange(0, 40);
   else
      mu->setRange(-40, 40);

   //   RooAbsPdf* pdf = mc->GetPdf();
   pdf = mc->GetPdf();

   RooArgSet  nuis_tmp1 = *mc->GetNuisanceParameters();
   RooNLLVar *asimov_nll =
      (RooNLLVar *)pdf->createNLL(*asimovData1, Constrain(nuis_tmp1), Offset(1), Optimize(2), NumCPU(nCPU, 3));

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
         mu_inj = mu_init; // for the mass point at the inj
      }
      RooDataSet *injData1 =
         EXOSTATS::makeAsimovData(ws, mc->GetName(), doConditional, obs_nll, 0, &mu_str, &mu_prof_str, 1, true, mu_inj);
      string globObsSnapName = "conditionalGlobs" + mu_prof_str;
      ws->loadSnapshot(globObsSnapName.c_str());
      RooNLLVar *inj_nll =
         (RooNLLVar *)pdf->createNLL(*injData1, Constrain(nuis_tmp2), Offset(1), Optimize(2), NumCPU(nCPU, 3));

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
   const TString outFileName = fullOutFolder + workspaceTag + blindness + "_sig.root";
   TFile         f_out(outFileName, "RECREATE");
   tree->SetDirectory(&f_out);
   f_out.Write();

   timer.Stop();
   timer.Print();
}
