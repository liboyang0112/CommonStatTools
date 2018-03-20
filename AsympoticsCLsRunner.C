#include <sstream>
#include <Math/MinimizerOptions.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>
#include <RooNLLVar.h>
#include <RooDataSet.h>
#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>
#include <RooRealVar.h>
#include <RooRealSumPdf.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "AsimovDataMaking.h"
#include "Minimization.h"
#include "AsympoticsCLsRunner.h"

using namespace RooFit;

void EXOSTATS::AsymptoticsCLsRunner::init()
{
   // band configuration
   betterBands = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
   betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
   profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero
   // other configuration
   defaultMinimizer    = "Minuit2"; // or "Minuit"
   defaultPrintLevel   = 1;         // Minuit print level
   defaultStrategy     = 2;         // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
   killBelowFatal      = 1;         // In case you want to suppress RooFit warnings further, set to 1
   doBlind             = 0;         // in case your analysis is blinded
   conditionalExpected = 1 && !doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
   doTilde             = 1;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
   doExp               = 1;             // compute expected limit
   doObs               = 1 && !doBlind; // compute observed limit
   doInj               = 0;             // compute expected limit after signal injection
   precision           = 0.005;         // % precision in mu that defines iterative cutoff
   verbose             = 0;             // 1 = very spammy
   usePredictiveFit    = 1;    // experimental, extrapolate best fit nuisance parameters based on previous fit results
   extrapolateSigma    = 1;    // experimantal, extrapolate sigma based on previous fits
   maxRetries          = 3;    // number of minimize(fcn) retries before giving up
   doPvals             = true; // perform pvalue calculation

   // don't touch!
   map_nll_muhat.clear();
   map_muhat.clear();
   map_data_nll.clear();
   map_snapshots.clear();
   map_nll_mu_sigma.clear();
   w             = nullptr;
   mc            = nullptr;
   data          = nullptr;
   firstPOI      = nullptr;
   asimov_0_nll  = nullptr;
   asimov_1_nll  = nullptr;
   obs_nll       = nullptr;
   nrMinimize    = 0;
   direction     = 1;
   global_status = 0;
   target_CLs    = 0.05;
   // range of firstPOI from ModelConfig mc
   firstPOIMax = 0;
   firstPOIMin = 0;
}

int EXOSTATS::AsymptoticsCLsRunner::minimize(RooNLLVar *nll)
{
   nrMinimize++;
   return EXOSTATS::minimize(nll, maxRetries);
}

void EXOSTATS::AsymptoticsCLsRunner::runAsymptoticsCLs(const char *infile, const char *workspaceName,
                                                       const char *modelConfigName, const char *dataName,
                                                       const char *asimovDataName, const char *conditionalSnapshot,
                                                       const char *nominalSnapshot, string folder, string mass,
                                                       double CL, bool doBetterBands)
{
   stringstream smass;
   smass << mass;

   conditionalSnapshot = ""; // warningless compile
   nominalSnapshot     = ""; // warningless compile

   runAsymptoticsCLs(infile, workspaceName, modelConfigName, dataName, asimovDataName, folder, smass.str(), CL,
                     doBetterBands, 0);
}

void EXOSTATS::AsymptoticsCLsRunner::runAsymptoticsCLs(const char *infile, const char *workspaceName,
                                                       const char *modelConfigName, const char *dataName,
                                                       const char *asimovDataName, string folder, string mass,
                                                       double CL, bool doBetterBands, double mu_inj)
{
   TStopwatch timer;
   timer.Start();
   betterBands = doBetterBands;
   if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
   // RooNLLVar::SetIgnoreZeroEntries(1);

   // check inputs
   TFile f(infile);
   w = (RooWorkspace *)f.Get(workspaceName);
   if (!w) {
      cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
      return;
   }
   RooFIter   rfiter = w->components().fwdIterator();
   RooAbsArg *arg;
   while ((arg = rfiter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
         arg->setAttribute("BinnedLikelihood");
      }
   }

   mc = (RooStats::ModelConfig *)w->obj(modelConfigName);
   if (!mc) {
      cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
      return;
   }
   firstPOI    = (RooRealVar *)mc->GetParametersOfInterest()->first();
   firstPOIMax = firstPOI->getMax();
   firstPOIMin = firstPOI->getMin();
   if (verbose) {
      cout << "runAsymptoticsCLs: get min and max of firstPOI " << endl;
      cout << "firstPOIMin " << firstPOIMin << endl;
      cout << "firstPOIMax " << firstPOIMax << endl;
   }

   data = (RooDataSet *)w->data(dataName);
   if (!data) {
      cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
      return;
   }

   // RooAbsPdf* pdf = mc->GetPdf();
   obs_nll                = createNLL(data); //(RooNLLVar*)pdf->createNLL(*data);
   map_snapshots[obs_nll] = "nominalGlobs";
   map_data_nll[data]     = obs_nll;
   w->saveSnapshot("nominalGlobs", *mc->GetGlobalObservables());
   w->saveSnapshot("nominalNuis", (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet()));

   global_status            = 0;
   RooDataSet *asimovData_0 = (RooDataSet *)w->data(asimovDataName);
   RooDataSet *asimovData_1 = EXOSTATS::makeAsimovData(w, modelConfigName, conditionalExpected, obs_nll, 1);
   if (!asimovData_0) {
      asimovData_0 = EXOSTATS::makeAsimovData(w, modelConfigName, conditionalExpected, obs_nll, 0);

      // asimovData_0 = makeAsimovData2((conditionalExpected ? obs_nll : (RooNLLVar*)nullptr), 0., 0.);
   }
   int asimov0_status = global_status;

   asimov_0_nll                = createNLL(asimovData_0); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
   asimov_1_nll                = createNLL(asimovData_1); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
   map_snapshots[asimov_0_nll] = "conditionalGlobs_0";
   map_data_nll[asimovData_0]  = asimov_0_nll;
   setMu(0);
   map_muhat[asimov_0_nll] = 0;
   saveSnapshot(asimov_0_nll, 0);
   w->loadSnapshot("conditionalNuis_0");
   w->loadSnapshot("conditionalGlobs_0");
   map_nll_muhat[asimov_0_nll] = asimov_0_nll->getVal();

   target_CLs = 1 - CL;

   double med_limit = -1;
   double med_muhat = -1;

   if (doExp) {
      cout << "Calculating Expected Limit" << endl;
      getLimit(asimov_0_nll, 1.0, med_limit, med_muhat);
   }

   int med_status = global_status;

   double inj_limit  = 0;
   int    inj_status = 0;

   if (doInj) {
      RooDataSet *asimovData_inj = EXOSTATS::makeAsimovData(w, modelConfigName, conditionalExpected, obs_nll, 0,
                                                            nullptr, nullptr, -999, true, mu_inj);

      int asimovinj_status = global_status;

      RooNLLVar *asimov_inj_nll     = createNLL(asimovData_inj); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
      map_snapshots[asimov_inj_nll] = "conditionalGlobs_0";
      map_data_nll[asimovData_inj]  = asimov_inj_nll;
      setMu(0);
      w->loadSnapshot("conditionalNuis_0");
      w->loadSnapshot("conditionalGlobs_0");

      target_CLs = 1 - CL;
      inj_limit  = getLimit(asimov_inj_nll, med_limit);
      inj_status = global_status;
   }

   double sigma = med_limit / sqrt(3.84); // pretty close
   double mu_up_p2_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - target_CLs * ROOT::Math::gaussian_cdf(2), 1) + 2);
   double mu_up_p1_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - target_CLs * ROOT::Math::gaussian_cdf(1), 1) + 1);
   double mu_up_n1_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - target_CLs * ROOT::Math::gaussian_cdf(-1), 1) - 1);
   double mu_up_n2_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - target_CLs * ROOT::Math::gaussian_cdf(-2), 1) - 2);

   double mu_up_p2 = mu_up_p2_approx;
   double mu_up_p1 = mu_up_p1_approx;
   double mu_up_n1 = mu_up_n1_approx;
   double mu_up_n2 = mu_up_n2_approx;

   firstPOI->setRange(-5 * sigma, 5 * sigma);
   map<int, int> N_status;
   if (betterBands && doExp) { // no better time than now to do this
      // find quantiles, starting with +2, since we should be at +1.96 right now

      double init_targetCLs = target_CLs;
      firstPOI->setRange(-5 * sigma, 5 * sigma);
      for (int N = 2; N >= -2; N--) {
         if (N < 0 && !betterNegativeBands) continue;
         if (N == 0) continue;
         target_CLs =
            2 * (1 - ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
         if (N < 0) direction = -1;

         // get the acual value
         double NtimesSigma =
            getLimit(asimov_0_nll, N * med_limit / sqrt(3.84)); // use N * sigma(0) as an initial guess
         N_status[N] += global_status;
         sigma = NtimesSigma / N;
         cout << endl;
         cout << "Found N * sigma = " << N << " * " << sigma << endl;

         string muStr, muStrPr;
         w->loadSnapshot("conditionalGlobs_0");
         double pr_val = NtimesSigma;
         if (N < 0 && profileNegativeAtZero) pr_val = 0;
         RooDataSet *asimovData_N =
            EXOSTATS::makeAsimovData(w, modelConfigName, 1, asimov_0_nll, NtimesSigma, &muStr, &muStrPr, pr_val, 0);
         // RooDataSet* asimovData_N = makeAsimovData2(asimov_0_nll, NtimesSigma, pr_val, &muStr, &muStrPr);

         RooNLLVar *asimov_N_nll     = createNLL(asimovData_N); //(RooNLLVar*)pdf->createNLL(*asimovData_N);
         map_data_nll[asimovData_N]  = asimov_N_nll;
         map_snapshots[asimov_N_nll] = "conditionalGlobs" + muStrPr;
         w->loadSnapshot(map_snapshots[asimov_N_nll].c_str());
         w->loadSnapshot(("conditionalNuis" + muStrPr).c_str());
         setMu(NtimesSigma);

         double nll_val = asimov_N_nll->getVal();
         saveSnapshot(asimov_N_nll, NtimesSigma);
         map_muhat[asimov_N_nll] = NtimesSigma;
         if (N < 0 && doTilde) {
            setMu(0);
            firstPOI->setConstant(1);
            nll_val = getNLL(asimov_N_nll);
         }
         map_nll_muhat[asimov_N_nll] = nll_val;

         target_CLs           = init_targetCLs;
         direction            = 1;
         double initial_guess = findCrossing(NtimesSigma / N, NtimesSigma / N, NtimesSigma);
         double limit         = getLimit(asimov_N_nll, initial_guess);
         N_status[N] += global_status;

         if (N == 2)
            mu_up_p2 = limit;
         else if (N == 1)
            mu_up_p1 = limit;
         else if (N == -1)
            mu_up_n1 = limit;
         else if (N == -2)
            mu_up_n2 = limit;
         // return;
      }
      direction  = 1;
      target_CLs = init_targetCLs;
   }

   w->loadSnapshot("conditionalNuis_0");
   firstPOI->setRange(firstPOIMin, firstPOIMax);
   double obs_limit = -1;
   double obs_muhat = -1;
   if (doObs) {
      cout << "Calculating Observed Limit" << endl;
      getLimit(obs_nll, med_limit, obs_limit, obs_muhat);
   }
   int obs_status = global_status;

   bool hasFailures = false;
   if (obs_status != 0 || med_status != 0 || asimov0_status != 0 || inj_status != 0) hasFailures = true;
   for (map<int, int>::iterator itr = N_status.begin(); itr != N_status.end(); itr++) {
      if (itr->second != 0) hasFailures = true;
   }
   if (hasFailures) {
      cout << "--------------------------------" << endl;
      cout << "Unresolved fit failures detected" << endl;
      cout << "Asimov0:  " << asimov0_status << endl;
      for (map<int, int>::iterator itr = N_status.begin(); itr != N_status.end(); itr++) {
         cout << "+" << itr->first << "sigma:  " << itr->first << endl;
      }
      cout << "Median:   " << med_status << endl;
      cout << "Injected: " << inj_status << endl;
      cout << "Observed: " << obs_status << endl;
      cout << "--------------------------------" << endl;
   }

   if (betterBands) cout << "Guess for bands" << endl;
   cout << "+2sigma:  " << mu_up_p2_approx << endl;
   cout << "+1sigma:  " << mu_up_p1_approx << endl;
   cout << "-1sigma:  " << mu_up_n1_approx << endl;
   cout << "-2sigma:  " << mu_up_n2_approx << endl;
   if (betterBands) {
      cout << endl;
      cout << "Correct bands" << endl;
      cout << "+2sigma:  " << mu_up_p2 << endl;
      cout << "+1sigma:  " << mu_up_p1 << endl;
      cout << "-1sigma:  " << mu_up_n1 << endl;
      cout << "-2sigma:  " << mu_up_n2 << endl;
   }

   cout << "Injected: " << inj_limit << endl;
   cout << "Median:   " << med_limit << endl;
   cout << "Observed: " << obs_limit << endl;
   cout << endl;

   // Pvalues
   double med_sig  = -1;
   double med_CLb  = -1;
   double med_CLsb = -1;
   double med_CLs  = -1;

   double obs_sig  = -1;
   double obs_CLb  = -1;
   double obs_CLsb = -1;
   double obs_CLs  = -1;

   if (doPvals) {
      getExpPvalue(med_CLb);
      med_CLb  = 1 - med_CLb;
      med_CLsb = 0.5;
      med_CLs  = med_CLsb / (med_CLb + 1e-9);

      getObsPvalue(0, obs_CLb);
      obs_CLb = 1 - obs_CLb;
      getObsPvalue(1, obs_CLsb);
      obs_CLs = obs_CLsb / (obs_CLb + 1e-9);

      cout << "************************" << endl;
      cout << "*-* Expected Pvalues *-*" << endl;
      cout << "************************" << endl;
      cout << "pb = " << 1 - med_CLb << endl;
      cout << "CLb = (1-pb) = " << med_CLb << endl;
      cout << "CLsb = " << med_CLsb << endl;
      cout << "CLs = " << med_CLs << endl;
      cout << "************************" << endl;
      cout << endl;
      cout << "************************" << endl;
      cout << "*-* Observed Pvalues *-*" << endl;
      cout << "************************" << endl;
      cout << "pb = " << 1 - obs_CLb << endl;
      cout << "CLb = (1-pb) = " << obs_CLb << endl;
      cout << "CLsb = " << obs_CLsb << endl;
      cout << "CLs = " << obs_CLs << endl;
      cout << "************************" << endl;
   }

   system(("mkdir -vp root-files/" + folder).c_str());

   stringstream fileName;
   fileName << "root-files/" << folder << "/" << mass << ".root";
   TFile  fout(fileName.str().c_str(), "recreate");
   TTree *tree = new TTree("stats", "runAsymptoticsCLs");

   Float_t tree_point;
   //  Float_t tree_null_pvalue;
   //  Float_t tree_null_pvalue_err;
   //  Float_t tree_alt_pvalue;
   Float_t tree_CLb_med;
   Float_t tree_CLs_med;
   Float_t tree_CLsplusb_med;
   Float_t tree_CLb_obs;
   Float_t tree_pb_obs;
   Float_t tree_pb_med;
   Float_t tree_CLs_obs;
   Float_t tree_CLsplusb_obs;
   //  Float_t tree_obs_lowerlimit;
   //  Float_t tree_obs_lowerlimit_err;
   Float_t tree_obs_upperlimit;
   // Float_t tree_obs_upperlimit_err;
   Float_t tree_exp_upperlimit;
   Float_t tree_exp_upperlimit_plus1;
   Float_t tree_exp_upperlimit_plus2;
   Float_t tree_exp_upperlimit_minus1;
   Float_t tree_exp_upperlimit_minus2;
   Float_t tree_muhat_obs;
   Float_t tree_muhat_exp;
   Float_t tree_fit_status;

   tree->Branch("point", &tree_point, "point/F");
   //  tree->Branch("null_pvalue", &tree_null_pvalue, "null_pvalue/F");
   //  tree->Branch("null_pvalue_err", &tree_null_pvalue_err, "null_pvalue_err/F");
   //  tree->Branch("alt_pvalue", &tree_alt_pvalue, "alt_pvalue/F");
   tree->Branch("CLb_med", &tree_CLb_med, "CLb_med/F");
   tree->Branch("pb_med", &tree_pb_med, "pb_med/F");
   tree->Branch("CLs_med", &tree_CLs_med, "CLs_med/F");
   tree->Branch("CLsplusb_med", &tree_CLsplusb_med, "CLsplusb_med/F");
   tree->Branch("CLb_obs", &tree_CLb_obs, "CLb_obs/F");
   tree->Branch("pb_obs", &tree_pb_obs, "pb_obs/F");
   tree->Branch("CLs_obs", &tree_CLs_obs, "CLs_obs/F");
   tree->Branch("CLsplusb_obs", &tree_CLsplusb_obs, "CLsplusb_obs/F");
   //  tree->Branch("obs_lowerlimit", &tree_obs_lowerlimit, "obs_lowerlimit/F");
   // tree->Branch("obs_lowerlimit_err", &tree_obs_lowerlimit_err, "obs_lowerlimit_err/F");
   tree->Branch("obs_upperlimit", &tree_obs_upperlimit, "obs_upperlimit/F");
   // tree->Branch("obs_upperlimit_err", &tree_obs_upperlimit_err, "obs_upperlimit_err/F");
   tree->Branch("exp_upperlimit", &tree_exp_upperlimit, "exp_upperlimit/F");
   tree->Branch("exp_upperlimit_plus1", &tree_exp_upperlimit_plus1, "exp_upperlimit_plus1/F");
   tree->Branch("exp_upperlimit_plus2", &tree_exp_upperlimit_plus2, "exp_upperlimit_plus2/F");
   tree->Branch("exp_upperlimit_minus1", &tree_exp_upperlimit_minus1, "exp_upperlimit_minus1/F");
   tree->Branch("exp_upperlimit_minus2", &tree_exp_upperlimit_minus2, "exp_upperlimit_minus2/F");
   tree->Branch("fit_status", &tree_fit_status, "fit_status/F");
   tree->Branch("mu_hat_obs", &tree_muhat_obs, "mu_hat_obs/F");
   tree->Branch("mu_hat_exp", &tree_muhat_exp, "mu_hat_exp/F");

   tree_point        = atof(mass.c_str());
   tree_CLb_med      = med_CLb;
   tree_pb_med       = 1 - med_CLb;
   tree_CLs_med      = med_CLs;
   tree_CLsplusb_med = med_CLsb;
   tree_CLb_obs      = obs_CLb;
   tree_pb_obs       = 1 - obs_CLb;
   tree_CLs_obs      = obs_CLs;
   tree_CLsplusb_obs = obs_CLsb;

   //  tree_obs_lowerlimit = xxx;
   //  tree_obs_lowerlimit_err = kv.second.obs_lowerlimit_err;
   tree_obs_upperlimit = obs_limit;
   //  tree_obs_upperlimit_err = -1;//kv.second.obs_upperlimit_err;
   tree_exp_upperlimit = med_limit;

   tree_exp_upperlimit_plus1  = mu_up_p1;
   tree_exp_upperlimit_plus2  = mu_up_p2;
   tree_exp_upperlimit_minus1 = mu_up_n1;
   tree_exp_upperlimit_minus2 = mu_up_n2;

   tree_muhat_obs  = obs_muhat;
   tree_muhat_exp  = med_muhat;
   tree_fit_status = global_status;
   tree->Fill();

   fout.Write();
   fout.Close();

   cout << "Finished with " << nrMinimize << " calls to minimize(nll)" << endl;
   timer.Print();
}

double EXOSTATS::AsymptoticsCLsRunner::getLimit(RooNLLVar *nll, double initial_guess)
{
   double upperLimit, muhat;
   getLimit(nll, initial_guess, upperLimit, muhat);
   return upperLimit;
}

void EXOSTATS::AsymptoticsCLsRunner::getLimit(RooNLLVar *nll, double initial_guess, double &upper_limit, double &mu_hat)
{

   upper_limit = -1;
   mu_hat      = -1;

   cout << "------------------------" << endl;
   cout << "Getting limit for nll: " << nll->GetName() << endl;
   // get initial guess based on muhat and sigma(muhat)
   firstPOI->setConstant(0);
   global_status = 0;

   if (nll == asimov_0_nll) {
      setMu(0);
      firstPOI->setConstant(1);
   }

   double muhat;
   if (map_nll_muhat.find(nll) == map_nll_muhat.end()) {
      double nll_val = getNLL(nll);
      muhat          = firstPOI->getVal();
      saveSnapshot(nll, muhat);
      map_muhat[nll] = muhat;
      if (muhat < 0 && doTilde) {
         setMu(0);
         firstPOI->setConstant(1);
         nll_val = getNLL(nll);
      }

      map_nll_muhat[nll] = nll_val;
   } else {
      muhat = map_muhat[nll];
   }

   if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
   double qmu, qmuA;
   double sigma_guess = getSigma(asimov_0_nll, firstPOI->getVal(), 0, qmu);
   double sigma_b     = sigma_guess;
   double mu_guess    = findCrossing(sigma_guess, sigma_b, muhat);
   double pmu         = calcPmu(qmu, sigma_b, mu_guess);
   double pb          = calcPb(qmu, sigma_b, mu_guess);
   double CLs         = calcCLs(qmu, sigma_b, mu_guess);
   double qmu95       = getQmu95(sigma_b, mu_guess);
   setMu(mu_guess);

   cout << "Initial guess:  " << mu_guess << endl;
   cout << "Sigma(obs):     " << sigma_guess << endl;
   cout << "Sigma(mu,0):    " << sigma_b << endl;
   cout << "muhat:          " << muhat << endl;
   cout << "qmu95:          " << qmu95 << endl;
   cout << "qmu:            " << qmu << endl;
   cout << "pmu:            " << pmu << endl;
   cout << "1-pb:           " << pb << endl;
   cout << "CLs:            " << CLs << endl;
   cout << endl;

   int                 nrDamping = 1;
   map<double, double> guess_to_corr;
   double              damping_factor = 1.0;
   // double damping_factor_pre = damping_factor;
   int    nrItr   = 0;
   double mu_pre  = muhat; // mu_guess-10*precision*mu_guess;
   double mu_pre2 = muhat;
   while (fabs(mu_pre - mu_guess) > precision * mu_guess * direction) {
      cout << "----------------------" << endl;
      cout << "Starting iteration " << nrItr << " of " << nll->GetName() << endl;
      // do this to avoid comparing multiple minima in the conditional and unconditional fits
      if (nrItr == 0)
         loadSnapshot(nll, muhat);
      else if (usePredictiveFit)
         doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
      else
         loadSnapshot(asimov_0_nll, mu_pre);

      sigma_guess = getSigma(nll, mu_guess, muhat, qmu);
      saveSnapshot(nll, mu_guess);

      if (nll != asimov_0_nll) {
         if (nrItr == 0)
            loadSnapshot(asimov_0_nll, map_nll_muhat[asimov_0_nll]);
         else if (usePredictiveFit) {
            if (nrItr == 1)
               doPredictiveFit(nll, map_nll_muhat[asimov_0_nll], mu_pre, mu_guess);
            else
               doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
         } else
            loadSnapshot(asimov_0_nll, mu_pre);

         sigma_b = getSigma(asimov_0_nll, mu_guess, 0, qmuA);
         saveSnapshot(asimov_0_nll, mu_guess);
      } else {
         sigma_b = sigma_guess;
         qmuA    = qmu;
      }

      double corr = damping_factor * (mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
      for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
         if (fabs(itr->first - (mu_guess - corr)) < direction * mu_guess * 0.02 &&
             fabs(corr) > direction * mu_guess * precision) {
            damping_factor *= 0.8;
            cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
            if (nrDamping++ > 10) {
               nrDamping      = 1;
               damping_factor = 1.0;
            }
            corr *= damping_factor;
            break;
         }
      }

      // subtract off the difference in the new and damped correction
      guess_to_corr[mu_guess] = corr;
      mu_pre2                 = mu_pre;
      mu_pre                  = mu_guess;
      mu_guess -= corr;

      pmu   = calcPmu(qmu, sigma_b, mu_pre);
      pb    = calcPb(qmu, sigma_b, mu_pre);
      CLs   = calcCLs(qmu, sigma_b, mu_pre);
      qmu95 = getQmu95(sigma_b, mu_pre);

      cout << "NLL:            " << nll->GetName() << endl;
      cout << "Previous guess: " << mu_pre << endl;
      cout << "Sigma(obs):     " << sigma_guess << endl;
      cout << "Sigma(mu,0):    " << sigma_b << endl;
      cout << "muhat:          " << muhat << endl;
      cout << "pmu:            " << pmu << endl;
      cout << "1-pb:           " << pb << endl;
      cout << "CLs:            " << CLs << endl;
      cout << "qmu95:          " << qmu95 << endl;
      cout << "qmu:            " << qmu << endl;
      cout << "qmuA0:          " << qmuA << endl;
      cout << "Precision:      " << direction * mu_guess * precision << endl;
      cout << "Correction:    " << (-corr < 0 ? " " : "") << -corr << endl;
      cout << "New guess:      " << mu_guess << endl;
      cout << endl;

      nrItr++;
      if (nrItr > 25) {
         cout << "Infinite loop detected in getLimit(). Please intervene." << endl;
         break;
      }
   }

   cout << "Found limit for nll " << nll->GetName() << ": " << mu_guess << endl;
   cout << "Finished in " << nrItr << " iterations." << endl;
   cout << "NLL is " << nll->getVal() << endl;
   cout << endl;
   upper_limit = mu_guess;
   mu_hat      = muhat;

   return;
}

double EXOSTATS::AsymptoticsCLsRunner::getSigma(RooNLLVar *nll, double mu, double muhat, double &qmu)
{
   qmu = getQmu(nll, mu);
   if (verbose) cout << "qmu = " << qmu << endl;
   if (mu * direction < muhat)
      return fabs(mu - muhat) / sqrt(qmu);
   else if (muhat < 0 && doTilde)
      return sqrt(mu * mu - 2 * mu * muhat * direction) / sqrt(qmu);
   else
      return (mu - muhat) * direction / sqrt(qmu);
}

void EXOSTATS::AsymptoticsCLsRunner::getExpPvalue(double &pb)
{

   bool isConst = firstPOI->isConstant();
   w->loadSnapshot("conditionalNuis_1");
   setMu(0);
   firstPOI->setConstant(1);

   const RooArgSet *np  = mc->GetNuisanceParameters();
   RooRealVar *     var = (RooRealVar *)(np->first());
   var->setVal(var->getVal() + 0.1);

   int status = minimize(asimov_1_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      w->loadSnapshot("conditionalNuis_0");
      status = minimize(asimov_1_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double asimov_nll_cond = asimov_1_nll->getVal();

   setMu(1);
   w->loadSnapshot("conditionalNuis_1");

   const RooArgSet *np2  = mc->GetNuisanceParameters();
   RooRealVar *     var2 = (RooRealVar *)(np2->first());
   var2->setVal(var2->getVal() + 0.1);

   status = minimize(asimov_1_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      w->loadSnapshot("conditionalNuis_0");
      status = minimize(asimov_1_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double asimov_nll_min = asimov_1_nll->getVal();
   double asimov_q0      = 2 * (asimov_nll_cond - asimov_nll_min);

   int    sign     = int(asimov_q0 != 0 ? asimov_q0 / fabs(asimov_q0) : 0);
   double med_sig  = sign * sqrt(fabs(asimov_q0));
   double med_pval = (1 - TMath::Erf(med_sig / sqrt(2))) / 2;

   cout << "****************************" << endl;
   cout << "Expected P-value Calculation " << endl;
   cout << "****************************" << endl;
   cout << "Asimov q0: " << asimov_q0 << endl;
   cout << "Med sig: " << med_sig << endl;
   cout << "Med pval: " << med_pval << endl;
   cout << "****************************" << endl;

   firstPOI->setConstant(isConst);

   pb = med_pval;

   return;
}

void EXOSTATS::AsymptoticsCLsRunner::getObsPvalue(double mu, double &pval)
{

   bool isConst = firstPOI->isConstant();
   w->loadSnapshot("conditionalNuis_0");
   setMu(mu);
   firstPOI->setConstant(1);

   //  const RooArgSet *np = mc->GetNuisanceParameters();
   //  RooRealVar* var = (RooRealVar*)(np->first());
   //  var->setVal(var->getVal() + 0.1);

   int status = minimize(obs_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=0" << endl;
      w->loadSnapshot("conditionalNuis_0");
      status = minimize(obs_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double obs_nll_cond = obs_nll->getVal();

   //  setMu(1);
   firstPOI->setConstant(0);
   //  w->loadSnapshot("conditionalNuis_1");

   //  const RooArgSet *np2 = mc->GetNuisanceParameters();
   //  RooRealVar* var2 = (RooRealVar*)(np2->first());
   //  var2->setVal(var2->getVal() + 0.1);

   status = minimize(obs_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      w->loadSnapshot("conditionalNuis_0");
      status = minimize(obs_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double muhat       = firstPOI->getVal();
   double obs_nll_min = obs_nll->getVal();
   double obs_q0      = 2 * (obs_nll_cond - obs_nll_min);
   if (muhat < 0) obs_q0 = 0;

   int    sign     = int(obs_q0 != 0 ? obs_q0 / fabs(obs_q0) : 0);
   double obs_sig  = sign * sqrt(fabs(obs_q0));
   double obs_pval = (1 - TMath::Erf(obs_sig / sqrt(2))) / 2;

   cout << "****************************" << endl;
   cout << "Observed P-value Calculation " << endl;
   cout << "****************************" << endl;
   cout << "Obs q0: " << obs_q0 << endl;
   cout << "Obs sig: " << obs_sig << endl;
   cout << "Obs pval: " << obs_pval << endl;
   cout << "****************************" << endl;

   firstPOI->setConstant(isConst);

   pval = obs_pval;
   return;
}

double EXOSTATS::AsymptoticsCLsRunner::getQmu(RooNLLVar *nll, double mu)
{
   double nll_muhat = map_nll_muhat[nll];
   bool   isConst   = firstPOI->isConstant();
   firstPOI->setConstant(1);
   setMu(mu);
   double nll_val = getNLL(nll);
   firstPOI->setConstant(isConst);
   // cout << "qmu = 2 * (" << nll_val << " - " << nll_muhat << ")" << endl;
   return 2 * (nll_val - nll_muhat);
}

void EXOSTATS::AsymptoticsCLsRunner::saveSnapshot(RooNLLVar *nll, double mu)
{
   stringstream snapshotName;
   snapshotName << nll->GetName() << "_" << mu;
   w->saveSnapshot(snapshotName.str().c_str(),
                   (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet()));
}

void EXOSTATS::AsymptoticsCLsRunner::loadSnapshot(RooNLLVar *nll, double mu)
{
   stringstream snapshotName;
   snapshotName << nll->GetName() << "_" << mu;
   w->loadSnapshot(snapshotName.str().c_str());
}

void EXOSTATS::AsymptoticsCLsRunner::doPredictiveFit(RooNLLVar *nll, double mu1, double mu2, double mu)
{
   if (fabs(mu2 - mu) < direction * mu * precision * 4) {
      loadSnapshot(nll, mu2);
      return;
   }

   // extrapolate to mu using mu1 and mu2 assuming nuis scale linear in mu
   const RooArgSet *nuis = mc->GetNuisanceParameters();
   if (nuis != 0) {
      int     nrNuis    = nuis->getSize();
      double *theta_mu1 = new double[nrNuis];
      double *theta_mu2 = new double[nrNuis];

      TIterator * itr = nuis->createIterator();
      RooRealVar *var;
      int         counter = 0;
      loadSnapshot(nll, mu1);
      while ((var = (RooRealVar *)itr->Next())) {
         theta_mu1[counter++] = var->getVal();
      }

      itr->Reset();
      counter = 0;
      loadSnapshot(nll, mu2);
      while ((var = (RooRealVar *)itr->Next())) {
         theta_mu2[counter++] = var->getVal();
      }

      itr->Reset();
      counter = 0;
      while ((var = (RooRealVar *)itr->Next())) {
         double m            = (theta_mu2[counter] - theta_mu1[counter]) / (mu2 - mu1);
         double b            = theta_mu2[counter] - m * mu2;
         double theta_extrap = m * mu + b;

         var->setVal(theta_extrap);
         counter++;
      }

      delete itr;
      delete[] theta_mu1;
      delete[] theta_mu2;
   }
}

RooNLLVar *EXOSTATS::AsymptoticsCLsRunner::createNLL(RooDataSet *_data)
{
   const RooArgSet *nuis = mc->GetNuisanceParameters();
   RooNLLVar *      nll;
   if (nuis != 0)
      nll = (RooNLLVar *)mc->GetPdf()->createNLL(*_data, Constrain(*nuis), NumCPU(4, 3), Optimize(2), Offset(true));
   else
      nll = (RooNLLVar *)mc->GetPdf()->createNLL(*_data, NumCPU(4, 3), Optimize(2), Offset(true));
   return nll;
}

double EXOSTATS::AsymptoticsCLsRunner::getNLL(RooNLLVar *nll)
{
   string snapshotName = map_snapshots[nll];
   if (snapshotName != "") w->loadSnapshot(snapshotName.c_str());
   minimize(nll);
   double val = nll->getVal();
   w->loadSnapshot("nominalGlobs");
   return val;
}

double EXOSTATS::AsymptoticsCLsRunner::findCrossing(double sigma_obs, double sigma, double muhat)
{
   double mu_guess  = muhat + ROOT::Math::gaussian_quantile(1 - target_CLs, 1) * sigma_obs * direction;
   int    nrItr     = 0;
   int    nrDamping = 1;

   map<double, double> guess_to_corr;
   double              damping_factor = 1.0;
   double              mu_pre         = mu_guess - 10 * mu_guess * precision;
   while (fabs(mu_guess - mu_pre) > direction * mu_guess * precision) {
      double sigma_obs_extrap = sigma_obs;
      double eps              = 0;
      if (extrapolateSigma) {
         // map<double, double> map_mu_sigma = map_nll_mu_sigma[nll];
      }

      mu_pre = mu_guess;

      double qmu95 = getQmu95(sigma, mu_guess);
      double qmu;
      qmu = 1. / sigma_obs_extrap / sigma_obs_extrap * (mu_guess - muhat) * (mu_guess - muhat);
      if (muhat < 0 && doTilde)
         qmu = 1. / sigma_obs_extrap / sigma_obs_extrap * (mu_guess * mu_guess - 2 * mu_guess * muhat);

      double dqmu_dmu = 2 * (mu_guess - muhat) / sigma_obs_extrap / sigma_obs_extrap - 2 * qmu * eps;

      double corr = damping_factor * (qmu - qmu95) / dqmu_dmu;
      for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
         if (fabs(itr->first - mu_guess) < direction * mu_guess * precision) {
            damping_factor *= 0.8;
            if (verbose)
               cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
            if (nrDamping++ > 10) {
               nrDamping      = 1;
               damping_factor = 1.0;
            }
            corr *= damping_factor;
            break;
         }
      }
      guess_to_corr[mu_guess] = corr;

      mu_guess = mu_guess - corr;
      nrItr++;
      if (nrItr > 100) {
         cout << "Infinite loop detected in findCrossing. Please intervene." << endl;
         exit(1);
      }
      if (verbose)
         cout << "mu_guess = " << mu_guess << ", mu_pre = " << mu_pre << ", qmu = " << qmu << ", qmu95 = " << qmu95
              << ", sigma_obs_extrap = " << sigma_obs_extrap << ", sigma = " << sigma
              << ", direction*mu*prec = " << direction * mu_guess * precision << endl;
   }

   return mu_guess;
}

void EXOSTATS::AsymptoticsCLsRunner::setMu(double mu)
{
   if (mu != mu) {
      cout << "ERROR::POI gave nan. Please intervene." << endl;
      exit(1);
   }
   if (mu > 0 && firstPOI->getMax() < mu) firstPOI->setMax(2 * mu);
   if (mu < 0 && firstPOI->getMin() > mu) firstPOI->setMin(2 * mu);
   firstPOI->setVal(mu);
}

double EXOSTATS::AsymptoticsCLsRunner::getQmu95_brute(double sigma, double mu)
{
   double step_size = 0.001;
   double start     = step_size;
   if (mu / sigma > 0.2) start = 0;
   for (double qmu = start; qmu < 20; qmu += step_size) {
      double CLs = calcCLs(qmu, sigma, mu);

      if (CLs < target_CLs) return qmu;
   }

   return 20;
}

double EXOSTATS::AsymptoticsCLsRunner::getQmu95(double sigma, double mu)
{
   double qmu95 = 0;
   // no sane man would venture this far down into |mu/sigma|
   double target_N = ROOT::Math::gaussian_cdf(1 - target_CLs, 1);
   if (fabs(mu / sigma) < 0.25 * target_N) {
      qmu95 = 5.83 / target_N;
   } else {
      map<double, double> guess_to_corr;
      double              qmu95_guess    = pow(ROOT::Math::gaussian_quantile(1 - target_CLs, 1), 2);
      int                 nrItr          = 0;
      int                 nrDamping      = 1;
      double              damping_factor = 1.0;
      double              qmu95_pre      = qmu95_guess - 10 * 2 * qmu95_guess * precision;
      while (fabs(qmu95_guess - qmu95_pre) > 2 * qmu95_guess * precision) {
         qmu95_pre = qmu95_guess;
         if (verbose) {
            cout << "qmu95_guess = " << qmu95_guess << endl;
            cout << "CLs = " << calcCLs(qmu95_guess, sigma, mu) << endl;
            cout << "Derivative = " << calcDerCLs(qmu95_guess, sigma, mu) << endl;
         }

         double corr =
            damping_factor * (calcCLs(qmu95_guess, sigma, mu) - target_CLs) / calcDerCLs(qmu95_guess, sigma, mu);
         for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
            if (fabs(itr->first - qmu95_guess) < 2 * qmu95_guess * precision) {
               damping_factor *= 0.8;
               if (verbose)
                  cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
               if (nrDamping++ > 10) {
                  nrDamping      = 1;
                  damping_factor = 1.0;
               }
               corr *= damping_factor;
            }
         }

         guess_to_corr[qmu95_guess] = corr;
         qmu95_guess                = qmu95_guess - corr;

         if (verbose) {
            cout << "next guess = " << qmu95_guess << endl;
            cout << "precision = " << 2 * qmu95_guess * precision << endl;
            cout << endl;
         }
         nrItr++;
         if (nrItr > 200) {
            cout << "Infinite loop detected in getQmu95. Please intervene." << endl;
            exit(1);
         }
      }
      qmu95 = qmu95_guess;
   }

   if (qmu95 != qmu95) {
      qmu95 = getQmu95_brute(sigma, mu);
   }
   if (verbose) cout << "Returning qmu95 = " << qmu95 << endl;

   return qmu95;
}

double EXOSTATS::AsymptoticsCLsRunner::calcCLs(double qmu_tilde, double sigma, double mu)
{
   double pmu = calcPmu(qmu_tilde, sigma, mu);
   double pb  = calcPb(qmu_tilde, sigma, mu);
   if (verbose) {
      cout << "pmu = " << pmu << endl;
      cout << "pb = " << pb << endl;
   }
   if (pb == 1) return 0.5;
   return pmu / (1 - pb);
}

double EXOSTATS::AsymptoticsCLsRunner::calcPmu(double qmu, double sigma, double mu)
{
   double pmu;
   if (qmu < mu * mu / (sigma * sigma) || !doTilde) {
      pmu = 1 - ROOT::Math::gaussian_cdf(sqrt(qmu));
   } else {
      pmu = 1 - ROOT::Math::gaussian_cdf((qmu + mu * mu / (sigma * sigma)) / (2 * fabs(mu / sigma)));
   }
   if (verbose)
      cout << "for pmu, qmu = " << qmu << ", sigma = " << sigma << ", mu = " << mu << ", pmu = " << pmu << endl;
   return pmu;
}

double EXOSTATS::AsymptoticsCLsRunner::calcPb(double qmu, double sigma, double mu)
{
   if (qmu < mu * mu / (sigma * sigma) || !doTilde) {
      return 1 - ROOT::Math::gaussian_cdf(fabs(mu / sigma) - sqrt(qmu));
   } else {
      return 1 - ROOT::Math::gaussian_cdf((mu * mu / (sigma * sigma) - qmu) / (2 * fabs(mu / sigma)));
   }
}

double EXOSTATS::AsymptoticsCLsRunner::calcDerCLs(double qmu, double sigma, double mu)
{
   double dpmu_dq  = 0;
   double d1mpb_dq = 0;

   if (qmu < mu * mu / (sigma * sigma)) {
      double zmu = sqrt(qmu);
      dpmu_dq    = -1. / (2 * sqrt(qmu * 2 * TMath::Pi())) * exp(-zmu * zmu / 2);
   } else {
      double zmu = (qmu + mu * mu / (sigma * sigma)) / (2 * fabs(mu / sigma));
      dpmu_dq    = -1. / (2 * fabs(mu / sigma)) * 1. / (sqrt(2 * TMath::Pi())) * exp(-zmu * zmu / 2);
   }

   if (qmu < mu * mu / (sigma * sigma)) {
      double zb = fabs(mu / sigma) - sqrt(qmu);
      d1mpb_dq  = -1. / sqrt(qmu * 2 * TMath::Pi()) * exp(-zb * zb / 2);
   } else {
      double zb = (mu * mu / (sigma * sigma) - qmu) / (2 * fabs(mu / sigma));
      d1mpb_dq  = -1. / (2 * fabs(mu / sigma)) * 1. / (sqrt(2 * TMath::Pi())) * exp(-zb * zb / 2);
   }

   double pb = calcPb(qmu, sigma, mu);
   return dpmu_dq / (1 - pb) - calcCLs(qmu, sigma, mu) / (1 - pb) * d1mpb_dq;
}

/*
RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, double mu_prof, string* mu_str, string*
mu_prof_str)
{
  RooNLLVar* conditioningNLL = nullptr;
  if (conditioningData)
  {
    conditioningNLL = (RooNLLVar*)mc->GetPdf()->createNLL(*conditioningData);
  }
  return makeAsimovData2(conditioningNLL, mu_true, mu_prof, mu_str, mu_prof_str);
}


RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, double mu_prof, string* mu_str, string*
mu_prof_str)
{
  if (mu_prof == -999) mu_prof = mu_true;
  bool doTest = 0;

  cout << "Creating asimov data at mu = " << mu_true << ", profiling at mu = " << mu_prof << endl;
  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();

  int test = 0;
  if (doTest) cout << "test = " << test++ << endl;

  int _printLevel = 0;

  stringstream muStr;
  muStr << setprecision(5);
  muStr << "_" << mu_true;
  if (mu_str) *mu_str = muStr.str();

  stringstream muStrProf;
  muStrProf << setprecision(5);
  muStrProf << "_" << mu_prof;
  if (mu_prof_str) *mu_prof_str = muStrProf.str();


  if (doTest) cout << "test = " << test++ << endl;
  const RooArgSet* globs = mc->GetGlobalObservables();
  const RooArgSet* nuis = mc->GetNuisanceParameters();
  const RooArgSet* obs = mc->GetObservables();
  const RooArgSet* pois = mc->GetParametersOfInterest();

  TIterator* gItr = globs->createIterator();
  TIterator* nItr = nuis->createIterator();
  RooRealVar* var;

  //cout << "test = " << test++ << endl;
  RooArgSet emptySet;
  RooArgSet params(*nuis);
  params.add(*globs);
  params.add(*pois);
  w->saveSnapshot("initial_params", params);

  if (doTest) cout << "test = " << test++ << endl;

//condition the MLEs
  if (conditioningNLL)
  {
    //get the conditional MLEs
    firstPOI->setVal(mu_prof);
    firstPOI->setConstant(1);
    minimize(conditioningNLL);
  }

  if (doTest) cout << "test = " << test++ << endl;
  w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*nuis);


//to find the conditional globs, do a fit to the constraint only pdf with the globs floating and the MLEs constant
  RooArgSet obsCopy = *obs;
  RooArgSet nuisCopy = *nuis;

  RooArgSet constraints(*mc->GetPdf()->getAllConstraints(obsCopy, nuisCopy));
  RooRealVar minusOne("minusOne","minusOne",-1);
  constraints.add(minusOne);
  RooProduct constrFunc("constrFunc","constrFunc",constraints);

  if (doTest) cout << "test = " << test++ << endl;
  while ((var = (RooRealVar*)gItr->Next()))
  {
    var->setConstant(false);
  }
  gItr->Reset();

  while ((var = (RooRealVar*)nItr->Next()))
  {
    var->setConstant(true);
  }
  nItr->Reset();

  minimize(&constrFunc);

  while ((var = (RooRealVar*)gItr->Next()))
  {
    var->setConstant(true);
  }
  gItr->Reset();

  while ((var = (RooRealVar*)nItr->Next()))
  {
    var->setConstant(false);
  }
  nItr->Reset();

  w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*globs);


  if (doTest) cout << "test = " << test++ << endl;



//make the asimov data
  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mc->GetObservables());

  RooRealVar* weightVar = nullptr;
  if (!(weightVar = w->var(weightName)))
  {
    w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
    weightVar = w->var(weightName);
  }
  obsAndWeight.add(*w->var(weightName));


  if (doTest) cout << "test = " << test++ << endl;

  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
  map<string, RooDataSet*> asimovDataMap;

  //try fix for sim pdf
  RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
  TIterator* iter = channelCat->typeIterator() ;
  RooCatType* tt = nullptr;
  int nrIndices = 0;
  int iFrame=0;
  while((tt=(RooCatType*) iter->Next())) {
    nrIndices++;
  }
  for (int i=0;i<nrIndices;i++){
    channelCat->setIndex(i);
    iFrame++;
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;

    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

    if (_printLevel >= 1)
    {
      obstmp->Print();
      cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
    }

    RooDataSet* obsDataUnbinned = new
RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
    RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
    double expectedEvents = pdftmp->expectedEvents(*obstmp);
    double thisNorm = 0;
    for(int jj=0; jj<thisObs->numBins(); ++jj){
      thisObs->setBin(jj);

      thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
      if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18))
obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
    }

    if (_printLevel >= 1)
    {
      obsDataUnbinned->Print();
      cout <<"sum entries "<<obsDataUnbinned->sumEntries()<<endl;
    }
    if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }


    asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;

    if (_printLevel >= 1)
    {
      cout << "channel: " << channelCat->getLabel() << ", data: ";
      obsDataUnbinned->Print();
      cout << endl;
    }
  }

  if (doTest) cout << "test = " << test++ << endl;
  RooDataSet* asimovData = new
RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  if (w->data(("asimovData"+muStr.str()).c_str()))
  {
    w->import(*asimovData, true);
  }
  else
  {
    w->import(*asimovData);
  }




  w->loadSnapshot("initial_params");
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

  if (doTest) cout << "test = " << test++ << endl;
  return asimovData;

}
*/

void EXOSTATS::AsymptoticsCLsRunner::doExpected(bool isExpected)
{
   if (isExpected) {
      doBlind             = 1;             // in case your analysis is blinded
      conditionalExpected = 1 && !doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
      doExp               = 1;             // compute expected limit
      doObs               = 1 && !doBlind; // compute observed limit
   } else {
      doBlind             = 0;             // in case your analysis is blinded
      conditionalExpected = 1 && !doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
      doExp               = 1;             // compute expected limit
      doObs               = 1 && !doBlind; // compute observed limit
   }
}

void EXOSTATS::AsymptoticsCLsRunner::doPvalues(bool calc)
{
   doPvals = calc;
}

void EXOSTATS::AsymptoticsCLsRunner::doBetterBands(bool isBetterBands)
{
   betterBands = isBetterBands;
}

void EXOSTATS::AsymptoticsCLsRunner::doInjection(bool injection)
{
   doInj = injection;
}
