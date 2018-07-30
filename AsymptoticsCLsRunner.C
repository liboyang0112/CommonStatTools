#include <sstream>
#include <stdexcept>
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
#include <TIterator.h>

#include "AsimovDataMaking.h"
#include "Minimization.h"
#include "AsymptoticsCLsRunner.h"

using namespace RooFit;
using namespace std;

//////////////////////////////////////////////////////////////////////////
/// Utility function to create a map of the values of all nuisance parameters
///
/// \param[in] mc pointer to the RooStats::ModelConfig from which the list of nuisance parameters must be taken
///
/// Stores names and values of all parameters in a map.
std::map<TString, Float_t> getParameterValuesMap(RooStats::ModelConfig *mc)
{
   std::map<TString, Float_t> result;
   auto                       np = mc->GetNuisanceParameters();
   if (np) {
      auto iter = np->createIterator();
      auto var  = iter->Next();
      while (var) {
         auto rrv = dynamic_cast<RooRealVar *>(var);
         if (rrv) {
            result[rrv->GetName()] = (Float_t)rrv->getVal();
         }
         var = iter->Next();
      }
   }

   return result;
}

EXOSTATS::AsymptoticsCLsRunner::AsymptoticsCLsRunner()
{
   reset();
}

EXOSTATS::AsymptoticsCLsRunner::~AsymptoticsCLsRunner() {}

//////////////////////////////////////////////////////////////////////////
/// Resets configuration
///
void EXOSTATS::AsymptoticsCLsRunner::reset()
{
   // band configuration
   m_betterBands = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
   m_betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
   m_profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero

   // other configuration
   m_defaultMinimizer    = "Minuit2"; // or "Minuit"
   m_defaultPrintLevel   = 1;         // Minuit print level
   m_defaultStrategy     = 2;         // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
   m_killBelowFatal      = 1;         // In case you want to suppress RooFit warnings further, set to 1
   m_doBlind             = 0;         // in case your analysis is blinded
   m_conditionalExpected = 1 && !m_doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
   m_doTilde             = 1;               // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
   m_doExp               = 1;               // compute expected limit
   m_doObs               = 1 && !m_doBlind; // compute observed limit
   m_doInj               = 0;               // compute expected limit after signal injection
   m_muInjection         = 1;               // mu value to be used for signal injection
   m_precision           = 0.005;           // % precision in mu that defines iterative cutoff
   m_debugLevel          = 0;               // 0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent
   m_usePredictiveFit    = 1;    // experimental, extrapolate best fit nuisance parameters based on previous fit results
   m_extrapolateSigma    = 1;    // experimantal, extrapolate sigma based on previous fits
   m_maxRetries          = 3;    // number of minimize(fcn) retries before giving up
   m_doPvals             = true; // perform pvalue calculation
   m_NumCPU              = 4;    // for the parallelisation of the likelihood calculation

   // don't touch!
   m_map_nll_muhat.clear();
   m_map_muhat.clear();
   m_map_data_nll.clear();
   m_map_snapshots.clear();
   m_map_nll_mu_sigma.clear();
   m_w             = nullptr;
   m_mc            = nullptr;
   m_data          = nullptr;
   m_firstPOI      = nullptr;
   m_asimov_0_nll  = nullptr;
   m_asimov_1_nll  = nullptr;
   m_obs_nll       = nullptr;
   m_nrMinimize    = 0;
   m_direction     = 1;
   m_global_status = 0;
   m_target_CLs    = 0.05;
   // range of m_firstPOI from ModelConfig m_mc
   m_firstPOIMax = 0;
   m_firstPOIMin = 0;
}

//////////////////////////////////////////////////////////////////////////
/// Wrapper around EXOSTATS::minimize(RooNLLVar *nll, Int_t m_maxRetries)
///
/// \param[in] nll pointer to the RooNLLVar to be minimized
///
/// This fnuction keeps track of the number of calls to the EXOSTATS::minimize() function
int EXOSTATS::AsymptoticsCLsRunner::minimize(RooNLLVar *nll)
{
   if (m_debugLevel == 0) cout << "Call #" << m_nrMinimize << " to minimize" << endl;
   m_nrMinimize++;
   return EXOSTATS::minimize(nll, m_maxRetries);
}

//////////////////////////////////////////////////////////////////////////
/// Runs profile likelihood limits with CLs and asymptotic approximation
///
/// \param[in] inputFile name of the input file containing the workspace
/// \param[in] workspaceName name of the workspace
/// \param[in] modelConfigName name of the ModelConfig object to be retrieved from the workspace
/// \param[in] dataName name of the dataset to be fitted
/// \param[in] paramName name of the parameter which characterizes this workspace (e.g. "m_Higgs"); it will be added as a branch of the output TTree
/// \param[in] paramValue value of the parameter which characterizes this workspace (e.g. 125.0); it will be stored in the corresponding branch of the output TTree
/// \param[in] workspaceTag prefix for the output ROOT file
/// \param[in] outputFolder path under which the output ROOT file will be stored; it will be created if it does not exist
/// \param[in] CL confidence level for the calculation (e.g. 0.95 for 95% CL limits)
/// \param[in] asimovDataName name of the Asimov dataset to be used; if empty, will be created
/// 
/// This function is a wrapper around the computeLimit() method. It allows standalone
/// usage of the class, since it deals both with input workspace retrieval and output
/// tree storage in a TFile. For more details on the output TTree, see computeLimit().
void EXOSTATS::AsymptoticsCLsRunner::run(const char *inputFile, const char *workspaceName, const char *modelConfigName,
                                         const char *dataName, TString paramName, float paramValue,
                                         TString workspaceTag, TString outputFolder, Double_t CL,
                                         const char *asimovDataName)
{
   // check inputs
   TFile         f(inputFile);
   RooWorkspace *workspace = (RooWorkspace *)f.Get(workspaceName);
   if (!workspace) {
      cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
      return;
   }

   TTree *tree = computeLimit(workspace, modelConfigName, dataName, paramName, paramValue, CL, asimovDataName);

   const TString blindness     = (m_doBlind) ? "_BLIND" : "";
   const TString clstring      = TString::Format("_CL%2.0f", CL * 100.0);
   const TString fullOutFolder = outputFolder + "/asymptotics/";
   system("mkdir -vp " + fullOutFolder);
   const TString outFileName = fullOutFolder + workspaceTag + blindness + clstring + ".root";
   TFile         f_out(outFileName, "RECREATE");
   tree->SetDirectory(&f_out);
   f_out.Write();
}

//////////////////////////////////////////////////////////////////////////
/// Runs profile likelihood limits with CLs and asymptotic approximation
///
/// \param[in] workspace pointer to the RooWorkspace object to be considered
/// \param[in] modelConfigName name of the ModelConfig object to be retrieved from the workspace
/// \param[in] dataName name of the dataset to be fitted
/// \param[in] paramName name of the parameter which characterizes this workspace (e.g. "m_Higgs"); it will be added as a branch of the output TTree
/// \param[in] paramValue value of the parameter which characterizes this workspace (e.g. 125.0); it will be stored in the corresponding branch of the output TTree
/// \param[in] CL confidence level for the calculation (e.g. 0.95 for 95% CL limits)
/// \param[in] asimovDataName name of the Asimov dataset to be used; if empty, will be created
/// \param[out] stats tree containing limit setting results
/// 
/// This method takes as an input a RooWorkspace, and specifications on which ModelConfig and dataset to be used, and
/// performs limit setting with a given CL.
/// The output TTree, named \c stats, contains a single entry, with standard branches with limit results, a special branch
/// called paramName which contains the value specified by paramValue, and a branch per nuisance parameter with the
/// values of each of them after the background-only and signal-plus-background fit.
/// For example:
/// \code
///   root[0] stats->Print()
///   mass      : mass/F                                                 *
///   CLb_med   : CLb_med/F                                              *
///   pb_med    : pb_med/F                                               *
///   CLs_med   : CLs_med/F                                              *
///   CLsplusb_med : CLsplusb_med/F                                      *
///   CLb_obs   : CLb_obs/F                                              *
///   pb_obs    : pb_obs/F                                               *
///   CLs_obs   : CLs_obs/F                                              *
///   CLsplusb_obs : CLsplusb_obs/F                                      *
///   obs_upperlimit : obs_upperlimit/F                                  *
///   inj_upperlimit : inj_upperlimit/F                                  *
///   exp_upperlimit : exp_upperlimit/F                                  *
///   exp_upperlimit_plus1 : exp_upperlimit_plus1/F                      *
///   exp_upperlimit_plus2 : exp_upperlimit_plus2/F                      *
///   exp_upperlimit_minus1 : exp_upperlimit_minus1/F                    *
///   exp_upperlimit_minus2 : exp_upperlimit_minus2/F                    *
///   fit_status : fit_status/F                                          *
///   mu_hat_obs : mu_hat_obs/F                                          *
///   mu_hat_exp : mu_hat_exp/F                                          *
///   param_alpha_G2_hat : param_alpha_G2_hat/F                          *
///   param_alpha_bkgER_hat : param_alpha_bkgER_hat/F                    *
///   param_alpha_bkg_NormSys_bkgCryo_hat :                              *
///   param_alpha_bkg_NormSys_bkgInt_hat :                               *
///   param_alpha_bkg_NormSys_bkgPMT_hat :                               *
///   param_alpha_sigQy_hat : param_alpha_sigQy_hat/F                    *
///   param_alpha_G2_med : param_alpha_G2_med/F                          *
///   param_alpha_bkgER_med : param_alpha_bkgER_med/F                    *
///   param_alpha_bkg_NormSys_bkgCryo_med :                              *
///   param_alpha_bkg_NormSys_bkgInt_med :                               *
///   param_alpha_bkg_NormSys_bkgPMT_med :                               *
///   param_alpha_sigQy_med : param_alpha_sigQy_med/F                    *
/// \endcode
TTree *EXOSTATS::AsymptoticsCLsRunner::computeLimit(RooWorkspace *workspace, const char *modelConfigName,
                                                    const char *dataName, TString paramName, float paramValue,
                                                    Double_t CL, const char *asimovDataName)
{
   TStopwatch timer;
   timer.Start();

   m_w          = workspace;
   m_target_CLs = 1 - CL;

   printOptionValues();

   if (m_killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(m_defaultMinimizer.c_str());
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(m_defaultStrategy);
   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(m_defaultPrintLevel);
   // RooNLLVar::SetIgnoreZeroEntries(1);

   RooFIter   rfiter = m_w->components().fwdIterator();
   RooAbsArg *arg;
   while ((arg = rfiter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
         arg->setAttribute("BinnedLikelihood");
      }
   }

   m_mc = (RooStats::ModelConfig *)m_w->obj(modelConfigName);
   if (!m_mc) {
      cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
      throw std::runtime_error("Invalid input");
      return nullptr;
   }
   m_firstPOI    = (RooRealVar *)m_mc->GetParametersOfInterest()->first();
   m_firstPOIMax = m_firstPOI->getMax();
   m_firstPOIMin = m_firstPOI->getMin();
   if (m_debugLevel == 0) {
      cout << "runAsymptoticsCLs: get min and max of m_firstPOI " << endl;
      cout << "firstPOIMin " << m_firstPOIMin << endl;
      cout << "firstPOIMax " << m_firstPOIMax << endl;
   }

   m_data = (RooDataSet *)m_w->data(dataName);
   if (!m_data) {
      cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
      throw std::runtime_error("Invalid input");
      return nullptr;
   }
   if (m_debugLevel == 0) {
      cout << "runAsymptoticsCLs: read data " << dataName << endl;
   }

   // RooAbsPdf* pdf = m_mc->GetPdf();
   m_obs_nll                  = createNLL(m_data); //(RooNLLVar*)pdf->createNLL(*m_data);
   m_map_snapshots[m_obs_nll] = "nominalGlobs";
   m_map_data_nll[m_data]     = m_obs_nll;
   m_w->saveSnapshot("nominalGlobs", *m_mc->GetGlobalObservables());
   m_w->saveSnapshot("nominalNuis", (m_mc->GetNuisanceParameters() ? *m_mc->GetNuisanceParameters() : RooArgSet()));
   if (m_debugLevel == 0) {
      cout << "runAsymptoticsCLs: saved nominal snapshots" << endl;
   }

   m_global_status = 0;
   if (m_debugLevel == 0) {
      cout << "runAsymptoticsCLs: preparing Asimov datasets" << endl;
   }
   RooDataSet *asimovData_0 = (RooDataSet *)m_w->data(asimovDataName);
   RooDataSet *asimovData_1 = EXOSTATS::makeAsimovData(m_w, modelConfigName, m_conditionalExpected, m_obs_nll, 1);
   if (!asimovData_0) {
      if (m_debugLevel == 0) {
         cout << "runAsymptoticsCLs: null hypothesis Asimov dataset does not exist, regenerating..." << endl;
      }
      asimovData_0 = EXOSTATS::makeAsimovData(m_w, modelConfigName, m_conditionalExpected, m_obs_nll, 0);

      // asimovData_0 = makeAsimovData2((m_conditionalExpected ? m_obs_nll : (RooNLLVar*)nullptr), 0., 0.);
   }
   int asimov0_status = m_global_status;

   if (m_debugLevel == 0) {
      cout << "runAsymptoticsCLs: creating NLL and mapping snapshots" << endl;
   }
   m_asimov_0_nll                  = createNLL(asimovData_0); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
   m_asimov_1_nll                  = createNLL(asimovData_1); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
   m_map_snapshots[m_asimov_0_nll] = "conditionalGlobs_0";
   m_map_data_nll[asimovData_0]    = m_asimov_0_nll;
   setMu(0);
   m_map_muhat[m_asimov_0_nll] = 0;
   saveSnapshot(m_asimov_0_nll, 0);
   m_w->loadSnapshot("conditionalNuis_0");
   m_w->loadSnapshot("conditionalGlobs_0");
   m_map_nll_muhat[m_asimov_0_nll] = m_asimov_0_nll->getVal();

   m_target_CLs = 1 - CL;

   double med_limit = -1;
   double med_muhat = -1;

   std::map<TString, Float_t> np_hatmed_map;
   if (m_doExp) {
      cout << "Calculating Expected Limit" << endl;
      getLimit(m_asimov_0_nll, 1.0, med_limit, med_muhat, np_hatmed_map);
   }

   int med_status = m_global_status;

   double inj_limit  = 0;
   int    inj_status = 0;

   if (m_doInj) {
      if (m_debugLevel == 0) {
         cout << "runAsymptoticsCLs: injecting signal with strength " << m_muInjection << endl;
      }
      RooDataSet *asimovData_inj = EXOSTATS::makeAsimovData(m_w, modelConfigName, m_conditionalExpected, m_obs_nll, 0,
                                                            nullptr, nullptr, -999, true, m_muInjection);

      RooNLLVar *asimov_inj_nll       = createNLL(asimovData_inj); //(RooNLLVar*)pdf->createNLL(*asimovData_0);
      m_map_snapshots[asimov_inj_nll] = "conditionalGlobs_0";
      m_map_data_nll[asimovData_inj]  = asimov_inj_nll;
      setMu(0);
      m_w->loadSnapshot("conditionalNuis_0");
      m_w->loadSnapshot("conditionalGlobs_0");

      m_target_CLs = 1 - CL;
      inj_limit    = getLimit(asimov_inj_nll, med_limit);
      inj_status   = m_global_status;
   }

   double sigma = med_limit / sqrt(3.84); // pretty close
   double mu_up_p2_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - m_target_CLs * ROOT::Math::gaussian_cdf(2), 1) + 2);
   double mu_up_p1_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - m_target_CLs * ROOT::Math::gaussian_cdf(1), 1) + 1);
   double mu_up_n1_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - m_target_CLs * ROOT::Math::gaussian_cdf(-1), 1) - 1);
   double mu_up_n2_approx =
      sigma * (ROOT::Math::gaussian_quantile(1 - m_target_CLs * ROOT::Math::gaussian_cdf(-2), 1) - 2);

   double mu_up_p2 = mu_up_p2_approx;
   double mu_up_p1 = mu_up_p1_approx;
   double mu_up_n1 = mu_up_n1_approx;
   double mu_up_n2 = mu_up_n2_approx;

   m_firstPOI->setRange(-5 * sigma, 5 * sigma);
   map<int, int> N_status;
   if (m_betterBands && m_doExp) { // no better time than now to do this
      // find quantiles, starting with +2, since we should be at +1.96 right now

      double init_targetCLs = m_target_CLs;
      m_firstPOI->setRange(-5 * sigma, 5 * sigma);
      for (int N = 2; N >= -2; N--) {
         if (N < 0 && !m_betterNegativeBands) continue;
         if (N == 0) continue;
         m_target_CLs =
            2 * (1 - ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
         if (N < 0) m_direction = -1;

         // get the acual value
         double NtimesSigma =
            getLimit(m_asimov_0_nll, N * med_limit / sqrt(3.84)); // use N * sigma(0) as an initial guess
         N_status[N] += m_global_status;
         sigma = NtimesSigma / N;
         cout << endl;
         cout << "Found N * sigma = " << N << " * " << sigma << endl;

         string muStr, muStrPr;
         m_w->loadSnapshot("conditionalGlobs_0");
         double pr_val = NtimesSigma;
         if (N < 0 && m_profileNegativeAtZero) pr_val = 0;
         RooDataSet *asimovData_N =
            EXOSTATS::makeAsimovData(m_w, modelConfigName, 1, m_asimov_0_nll, NtimesSigma, &muStr, &muStrPr, pr_val, 0);
         // RooDataSet* asimovData_N = makeAsimovData2(m_asimov_0_nll, NtimesSigma, pr_val, &muStr, &muStrPr);

         RooNLLVar *asimov_N_nll       = createNLL(asimovData_N); //(RooNLLVar*)pdf->createNLL(*asimovData_N);
         m_map_data_nll[asimovData_N]  = asimov_N_nll;
         m_map_snapshots[asimov_N_nll] = "conditionalGlobs" + muStrPr;
         m_w->loadSnapshot(m_map_snapshots[asimov_N_nll].c_str());
         m_w->loadSnapshot(("conditionalNuis" + muStrPr).c_str());
         setMu(NtimesSigma);

         double nll_val = asimov_N_nll->getVal();
         saveSnapshot(asimov_N_nll, NtimesSigma);
         m_map_muhat[asimov_N_nll] = NtimesSigma;
         if (N < 0 && m_doTilde) {
            setMu(0);
            m_firstPOI->setConstant(1);
            nll_val = getNLL(asimov_N_nll);
         }
         m_map_nll_muhat[asimov_N_nll] = nll_val;

         m_target_CLs         = init_targetCLs;
         m_direction          = 1;
         double initial_guess = findCrossing(NtimesSigma / N, NtimesSigma / N, NtimesSigma);
         double limit         = getLimit(asimov_N_nll, initial_guess);
         N_status[N] += m_global_status;

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
      m_direction  = 1;
      m_target_CLs = init_targetCLs;
   }

   m_w->loadSnapshot("conditionalNuis_0");
   m_firstPOI->setRange(m_firstPOIMin, m_firstPOIMax);
   double                     obs_limit = -1;
   double                     obs_muhat = -1;
   std::map<TString, Float_t> np_hat_map;
   if (m_doObs) {
      cout << "Calculating Observed Limit" << endl;
      getLimit(m_obs_nll, med_limit, obs_limit, obs_muhat, np_hat_map);
   }
   int obs_status = m_global_status;

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

   if (m_betterBands) cout << "Guess for bands" << endl;
   cout << "+2sigma:  " << mu_up_p2_approx << endl;
   cout << "+1sigma:  " << mu_up_p1_approx << endl;
   cout << "-1sigma:  " << mu_up_n1_approx << endl;
   cout << "-2sigma:  " << mu_up_n2_approx << endl;
   if (m_betterBands) {
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
   double med_CLb  = -1;
   double med_CLsb = -1;
   double med_CLs  = -1;

   double obs_CLb  = -1;
   double obs_CLsb = -1;
   double obs_CLs  = -1;

   if (m_doPvals) {
      getExpPvalue(med_CLb);
      med_CLb  = 1 - med_CLb;
      med_CLsb = 0.5;
      med_CLs  = med_CLsb / (med_CLb + 1e-9);

      if (m_doObs) {
         getObsPvalue(0, obs_CLb);
         obs_CLb = 1 - obs_CLb;
         getObsPvalue(1, obs_CLsb);
         obs_CLs = obs_CLsb / (obs_CLb + 1e-9);
      }

      cout << "************************" << endl;
      cout << "*-* Expected Pvalues *-*" << endl;
      cout << "************************" << endl;
      cout << "pb = " << 1 - med_CLb << endl;
      cout << "CLb = (1-pb) = " << med_CLb << endl;
      cout << "CLsb = " << med_CLsb << endl;
      cout << "CLs = " << med_CLs << endl;
      cout << "************************" << endl;
      cout << endl;

      if (m_doObs) {
         cout << "************************" << endl;
         cout << "*-* Observed Pvalues *-*" << endl;
         cout << "************************" << endl;
         cout << "pb = " << 1 - obs_CLb << endl;
         cout << "CLb = (1-pb) = " << obs_CLb << endl;
         cout << "CLsb = " << obs_CLsb << endl;
         cout << "CLs = " << obs_CLs << endl;
         cout << "************************" << endl;
      }
   }

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
   Float_t tree_inj_upperlimit;
   // Float_t tree_obs_upperlimit_err;
   Float_t tree_exp_upperlimit;
   Float_t tree_exp_upperlimit_plus1;
   Float_t tree_exp_upperlimit_plus2;
   Float_t tree_exp_upperlimit_minus1;
   Float_t tree_exp_upperlimit_minus2;
   Float_t tree_muhat_obs;
   Float_t tree_muhat_exp;
   Float_t tree_fit_status;

   tree->Branch(paramName, &tree_point);
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
   tree->Branch("inj_upperlimit", &tree_inj_upperlimit, "inj_upperlimit/F");
   // tree->Branch("obs_upperlimit_err", &tree_obs_upperlimit_err, "obs_upperlimit_err/F");
   tree->Branch("exp_upperlimit", &tree_exp_upperlimit, "exp_upperlimit/F");
   tree->Branch("exp_upperlimit_plus1", &tree_exp_upperlimit_plus1, "exp_upperlimit_plus1/F");
   tree->Branch("exp_upperlimit_plus2", &tree_exp_upperlimit_plus2, "exp_upperlimit_plus2/F");
   tree->Branch("exp_upperlimit_minus1", &tree_exp_upperlimit_minus1, "exp_upperlimit_minus1/F");
   tree->Branch("exp_upperlimit_minus2", &tree_exp_upperlimit_minus2, "exp_upperlimit_minus2/F");
   tree->Branch("fit_status", &tree_fit_status, "fit_status/F");
   tree->Branch("mu_hat_obs", &tree_muhat_obs, "mu_hat_obs/F");
   tree->Branch("mu_hat_exp", &tree_muhat_exp, "mu_hat_exp/F");

   for (auto &kv : np_hat_map) {
      const TString param = "param_" + kv.first + "_hat";
      tree->Branch(param, &kv.second, param + "/F");
   }
   for (auto &kv : np_hatmed_map) {
      const TString param = "param_" + kv.first + "_med";
      tree->Branch(param, &kv.second, param + "/F");
   }

   tree_point        = paramValue;
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
   tree_inj_upperlimit = inj_limit;
   //  tree_obs_upperlimit_err = -1;//kv.second.obs_upperlimit_err;
   tree_exp_upperlimit = med_limit;

   tree_exp_upperlimit_plus1  = mu_up_p1;
   tree_exp_upperlimit_plus2  = mu_up_p2;
   tree_exp_upperlimit_minus1 = mu_up_n1;
   tree_exp_upperlimit_minus2 = mu_up_n2;

   tree_muhat_obs  = obs_muhat;
   tree_muhat_exp  = med_muhat;
   tree_fit_status = m_global_status;
   tree->Fill();

   cout << "Finished with " << m_nrMinimize << " calls to minimize(nll)" << endl;
   timer.Print();

   return tree;
}

//////////////////////////////////////////////////////////////////////////
/// Wrapper around getLimit(RooNLLVar *nll, double initial_guess, double &upper_limit, double &mu_hat, std::map<TString, Float_t> &np_hat_map)
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] initial_guess initial guess for the parameter of interest best fit value
double EXOSTATS::AsymptoticsCLsRunner::getLimit(RooNLLVar *nll, double initial_guess)
{
   double                     upperLimit, muhat;
   std::map<TString, Float_t> dummy;
   getLimit(nll, initial_guess, upperLimit, muhat, dummy);
   return upperLimit;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates the upper limit on the signal strength
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] initial_guess initial guess for the parameter of interest best fit value (unconditional fit)
/// \param[in] upper_limit value of the upper limit on the parameter of interest
/// \param[in] mu_hat value of the best fit value of the parameter of interest
/// \param[in] np_hat_map reference to the map of nuisance parameter names vs values for the unconditional fit
void EXOSTATS::AsymptoticsCLsRunner::getLimit(RooNLLVar *nll, double initial_guess, double &upper_limit, double &mu_hat,
                                              std::map<TString, Float_t> &np_hat_map)
{

   upper_limit = -1;
   mu_hat      = -1;

   cout << "------------------------" << endl;
   cout << "Getting limit for nll: " << nll->GetName() << endl;
   // get initial guess based on muhat and sigma(muhat)
   m_firstPOI->setConstant(0);
   m_global_status = 0;

   if (nll == m_asimov_0_nll) {
      setMu(0);
      m_firstPOI->setConstant(1);
   }

   double                     muhat;
   std::map<TString, Float_t> nphatmap;
   if (m_map_nll_muhat.find(nll) == m_map_nll_muhat.end()) {
      double nll_val = getNLL(nll);
      muhat          = m_firstPOI->getVal();
      nphatmap       = getParameterValuesMap(m_mc);
      saveSnapshot(nll, muhat);
      m_map_muhat[nll] = muhat;
      if (muhat < 0 && m_doTilde) {
         setMu(0);
         m_firstPOI->setConstant(1);
         nll_val = getNLL(nll);
      }

      m_map_nll_muhat[nll] = nll_val;
   } else {
      muhat    = m_map_muhat[nll];
      nphatmap = getParameterValuesMap(m_mc);
   }

   if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
   double qmu, qmuA;
   double sigma_guess = getSigma(m_asimov_0_nll, m_firstPOI->getVal(), 0, qmu);
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
   while (fabs(mu_pre - mu_guess) > m_precision * mu_guess * m_direction) {
      cout << "----------------------" << endl;
      cout << "Starting iteration " << nrItr << " of " << nll->GetName() << endl;
      // do this to avoid comparing multiple minima in the conditional and unconditional fits
      if (nrItr == 0)
         loadSnapshot(nll, muhat);
      else if (m_usePredictiveFit)
         doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
      else
         loadSnapshot(m_asimov_0_nll, mu_pre);

      sigma_guess = getSigma(nll, mu_guess, muhat, qmu);
      saveSnapshot(nll, mu_guess);

      if (nll != m_asimov_0_nll) {
         if (nrItr == 0)
            loadSnapshot(m_asimov_0_nll, m_map_nll_muhat[m_asimov_0_nll]);
         else if (m_usePredictiveFit) {
            if (nrItr == 1)
               doPredictiveFit(nll, m_map_nll_muhat[m_asimov_0_nll], mu_pre, mu_guess);
            else
               doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
         } else
            loadSnapshot(m_asimov_0_nll, mu_pre);

         sigma_b = getSigma(m_asimov_0_nll, mu_guess, 0, qmuA);
         saveSnapshot(m_asimov_0_nll, mu_guess);
      } else {
         sigma_b = sigma_guess;
         qmuA    = qmu;
      }

      double corr = damping_factor * (mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
      for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
         if (fabs(itr->first - (mu_guess - corr)) < m_direction * mu_guess * 0.02 &&
             fabs(corr) > m_direction * mu_guess * m_precision) {
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
      cout << "Precision:      " << m_direction * mu_guess * m_precision << endl;
      cout << "Correction:    " << (-corr < 0 ? " " : "") << -corr << endl;
      cout << "New guess:      " << mu_guess << endl;
      cout << endl;

      nrItr++;
      if (nrItr > 25) {
         cout << "Infinite loop detected in getLimit(). Please intervene." << endl;
         throw std::runtime_error("Infinite loop");
         break;
      }
   }

   cout << "Found limit for nll " << nll->GetName() << ": " << mu_guess << endl;
   cout << "Finished in " << nrItr << " iterations." << endl;
   cout << "NLL is " << nll->getVal() << endl;
   cout << endl;
   upper_limit = mu_guess;
   mu_hat      = muhat;
   np_hat_map  = nphatmap;

   return;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates the sigma parameter
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] mu value of the parameter of interest
/// \param[in] best fit value of the parameter of interest (unconditional fit)
/// \param[in] reference to the value of the test statistics q_mu
/// \param[out] sigma
double EXOSTATS::AsymptoticsCLsRunner::getSigma(RooNLLVar *nll, double mu, double muhat, double &qmu)
{
   qmu = getQmu(nll, mu);
   if (m_debugLevel == 0) cout << "qmu = " << qmu << endl;
   if (mu * m_direction < muhat)
      return fabs(mu - muhat) / sqrt(qmu);
   else if (muhat < 0 && m_doTilde)
      return sqrt(mu * mu - 2 * mu * muhat * m_direction) / sqrt(qmu);
   else
      return (mu - muhat) * m_direction / sqrt(qmu);
}

//////////////////////////////////////////////////////////////////////////
/// Calculates the expected p-value of the background-only hypothesis
///
/// \param[in] pb reference to the p-value
void EXOSTATS::AsymptoticsCLsRunner::getExpPvalue(double &pb)
{

   bool isConst = m_firstPOI->isConstant();
   m_w->loadSnapshot("conditionalNuis_1");
   setMu(0);
   m_firstPOI->setConstant(1);

   const RooArgSet *np  = m_mc->GetNuisanceParameters();
   RooRealVar *     var = (RooRealVar *)(np->first());
   var->setVal(var->getVal() + 0.1);

   int status = minimize(m_asimov_1_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      m_w->loadSnapshot("conditionalNuis_0");
      status = minimize(m_asimov_1_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double asimov_nll_cond = m_asimov_1_nll->getVal();

   setMu(1);
   m_w->loadSnapshot("conditionalNuis_1");

   const RooArgSet *np2  = m_mc->GetNuisanceParameters();
   RooRealVar *     var2 = (RooRealVar *)(np2->first());
   var2->setVal(var2->getVal() + 0.1);

   status = minimize(m_asimov_1_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      m_w->loadSnapshot("conditionalNuis_0");
      status = minimize(m_asimov_1_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double asimov_nll_min = m_asimov_1_nll->getVal();
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

   m_firstPOI->setConstant(isConst);

   pb = med_pval;

   return;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates the observed p-value of the background-only hypothesis
///
/// \param[in] value of the parameter of interest to consider
/// \param[in] reference to the p-value
void EXOSTATS::AsymptoticsCLsRunner::getObsPvalue(double mu, double &pval)
{

   bool isConst = m_firstPOI->isConstant();
   m_w->loadSnapshot("conditionalNuis_0");
   setMu(mu);
   m_firstPOI->setConstant(1);

   //  const RooArgSet *np = m_mc->GetNuisanceParameters();
   //  RooRealVar* var = (RooRealVar*)(np->first());
   //  var->setVal(var->getVal() + 0.1);

   int status = minimize(m_obs_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=0" << endl;
      m_w->loadSnapshot("conditionalNuis_0");
      status = minimize(m_obs_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double obs_nll_cond = m_obs_nll->getVal();

   //  setMu(1);
   m_firstPOI->setConstant(0);
   //  m_w->loadSnapshot("conditionalNuis_1");

   //  const RooArgSet *np2 = m_mc->GetNuisanceParameters();
   //  RooRealVar* var2 = (RooRealVar*)(np2->first());
   //  var2->setVal(var2->getVal() + 0.1);

   status = minimize(m_obs_nll);
   if (status < 0) {
      cout << "Retrying with conditional snapshot at mu=1" << endl;
      m_w->loadSnapshot("conditionalNuis_0");
      status = minimize(m_obs_nll);
      if (status >= 0) cout << "Success!" << endl;
   }

   double muhat       = m_firstPOI->getVal();
   double obs_nll_min = m_obs_nll->getVal();
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

   m_firstPOI->setConstant(isConst);

   pval = obs_pval;
   return;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates the test statistics q_mu
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] mu value of the parameter of interest
double EXOSTATS::AsymptoticsCLsRunner::getQmu(RooNLLVar *nll, double mu)
{
   double nll_muhat = m_map_nll_muhat[nll];
   bool   isConst   = m_firstPOI->isConstant();
   m_firstPOI->setConstant(1);
   setMu(mu);
   double nll_val = getNLL(nll);
   m_firstPOI->setConstant(isConst);
   // cout << "qmu = 2 * (" << nll_val << " - " << nll_muhat << ")" << endl;
   return 2 * (nll_val - nll_muhat);
}

//////////////////////////////////////////////////////////////////////////
/// Saves a snapshot of all workspace parameters
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] mu value of the parameter of interest (used to name the snapshot)
void EXOSTATS::AsymptoticsCLsRunner::saveSnapshot(RooNLLVar *nll, double mu)
{
   stringstream snapshotName;
   snapshotName << nll->GetName() << "_" << mu;
   m_w->saveSnapshot(snapshotName.str().c_str(),
                     (m_mc->GetNuisanceParameters() ? *m_mc->GetNuisanceParameters() : RooArgSet()));
}

//////////////////////////////////////////////////////////////////////////
/// Loads the snapshot of all workspace parameters corresponding to the desired value of the parameter of interest
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] mu value of the parameter of interest
void EXOSTATS::AsymptoticsCLsRunner::loadSnapshot(RooNLLVar *nll, double mu)
{
   stringstream snapshotName;
   snapshotName << nll->GetName() << "_" << mu;
   m_w->loadSnapshot(snapshotName.str().c_str());
}

//////////////////////////////////////////////////////////////////////////
/// Performs the predictive fit
///
/// \param[in] nll pointer to the RooNLLVar to use
/// \param[in] mu1
/// \param[in] mu2
/// \param[in] mu
void EXOSTATS::AsymptoticsCLsRunner::doPredictiveFit(RooNLLVar *nll, double mu1, double mu2, double mu)
{
   if (fabs(mu2 - mu) < m_direction * mu * m_precision * 4) {
      loadSnapshot(nll, mu2);
      return;
   }

   // extrapolate to mu using mu1 and mu2 assuming nuis scale linear in mu
   const RooArgSet *nuis = m_mc->GetNuisanceParameters();
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

//////////////////////////////////////////////////////////////////////////
/// Creates the negative log-likelihood to be used in the fit
///
/// \param[in] _data pointer to the \c RooDataSet to compute the likelihood on
/// \param[out] nll pointer to the output \c RooNLLVar
///
/// Constraints are applied, as well as likelihood optimization. Multiple CPU running
/// is available if set (see setNumCPU()).
RooNLLVar *EXOSTATS::AsymptoticsCLsRunner::createNLL(RooDataSet *_data)
{
   const RooArgSet *nuis = m_mc->GetNuisanceParameters();
   RooNLLVar *      nll;
   if (nuis != 0)
      nll = (RooNLLVar *)m_mc->GetPdf()->createNLL(*_data, Constrain(*nuis), NumCPU(m_NumCPU, 3), Optimize(2),
                                                   Offset(true));
   else
      nll = (RooNLLVar *)m_mc->GetPdf()->createNLL(*_data, NumCPU(m_NumCPU, 3), Optimize(2), Offset(true));
   return nll;
}

//////////////////////////////////////////////////////////////////////////
/// Computes the best fit value of a negative log-likelihood
///
/// \param[in] nll pointer to the \c RooNLLVar to minimize
/// \param[out] val value of the negative log-likelihood after the fit
double EXOSTATS::AsymptoticsCLsRunner::getNLL(RooNLLVar *nll)
{
   string snapshotName = m_map_snapshots[nll];
   if (snapshotName != "") m_w->loadSnapshot(snapshotName.c_str());
   minimize(nll);
   double val = nll->getVal();
   m_w->loadSnapshot("nominalGlobs");
   return val;
}

//////////////////////////////////////////////////////////////////////////
/// Extrapolates the parameter of interest value
///
/// \param[in] sigma_obs
/// \param[in] sigma
/// \param[in] muhat
/// \param[out] mu_guess 
double EXOSTATS::AsymptoticsCLsRunner::findCrossing(double sigma_obs, double sigma, double muhat)
{
   double mu_guess  = muhat + ROOT::Math::gaussian_quantile(1 - m_target_CLs, 1) * sigma_obs * m_direction;
   int    nrItr     = 0;
   int    nrDamping = 1;

   map<double, double> guess_to_corr;
   double              damping_factor = 1.0;
   double              mu_pre         = mu_guess - 10 * mu_guess * m_precision;
   while (fabs(mu_guess - mu_pre) > m_direction * mu_guess * m_precision) {
      double sigma_obs_extrap = sigma_obs;
      double eps              = 0;
      if (m_extrapolateSigma) {
         // map<double, double> map_mu_sigma = m_map_nll_mu_sigma[nll];
      }

      mu_pre = mu_guess;

      double qmu95 = getQmu95(sigma, mu_guess);
      double qmu;
      qmu = 1. / sigma_obs_extrap / sigma_obs_extrap * (mu_guess - muhat) * (mu_guess - muhat);
      if (muhat < 0 && m_doTilde)
         qmu = 1. / sigma_obs_extrap / sigma_obs_extrap * (mu_guess * mu_guess - 2 * mu_guess * muhat);

      double dqmu_dmu = 2 * (mu_guess - muhat) / sigma_obs_extrap / sigma_obs_extrap - 2 * qmu * eps;

      double corr = damping_factor * (qmu - qmu95) / dqmu_dmu;
      for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
         if (fabs(itr->first - mu_guess) < m_direction * mu_guess * m_precision) {
            damping_factor *= 0.8;
            if (m_debugLevel == 0)
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
         throw std::runtime_error("Infinite loop");
         exit(1);
      }
      if (m_debugLevel == 0)
         cout << "mu_guess = " << mu_guess << ", mu_pre = " << mu_pre << ", qmu = " << qmu << ", qmu95 = " << qmu95
              << ", sigma_obs_extrap = " << sigma_obs_extrap << ", sigma = " << sigma
              << ", direction*mu*prec = " << m_direction * mu_guess * m_precision << endl;
   }

   return mu_guess;
}

//////////////////////////////////////////////////////////////////////////
/// Sets the parameter of interest to a given value, extending its range if needed.
///
/// \param[in] mu desided value of the parameter of interest
void EXOSTATS::AsymptoticsCLsRunner::setMu(double mu)
{
   if (mu != mu) {
      cout << "ERROR::POI gave nan. Please intervene." << endl;
      throw std::runtime_error("NaN encountered");
      exit(1);
   }
   if (mu > 0 && m_firstPOI->getMax() < mu) m_firstPOI->setMax(2 * mu);
   if (mu < 0 && m_firstPOI->getMin() > mu) m_firstPOI->setMin(2 * mu);
   m_firstPOI->setVal(mu);
}

//////////////////////////////////////////////////////////////////////////
/// Performs a brute-force calculation of q_mu95
///
/// \param[in] sigma
/// \param[in] mu
/// \param[out] qmu
double EXOSTATS::AsymptoticsCLsRunner::getQmu95_brute(double sigma, double mu)
{
   double step_size = 0.001;
   double start     = step_size;
   if (mu / sigma > 0.2) start = 0;
   for (double qmu = start; qmu < 20; qmu += step_size) {
      double CLs = calcCLs(qmu, sigma, mu);

      if (CLs < m_target_CLs) return qmu;
   }

   return 20;
}

//////////////////////////////////////////////////////////////////////////
/// Tries to calculate of q_mu95
///
/// \param[in] sigma
/// \param[in] mu
/// \param[out] qmu
///
/// Reverts to the brute-force calculation getQmu95_brute() if needed
double EXOSTATS::AsymptoticsCLsRunner::getQmu95(double sigma, double mu)
{
   double qmu95 = 0;
   // no sane man would venture this far down into |mu/sigma|
   double target_N = ROOT::Math::gaussian_cdf(1 - m_target_CLs, 1);
   if (fabs(mu / sigma) < 0.25 * target_N) {
      qmu95 = 5.83 / target_N;
   } else {
      map<double, double> guess_to_corr;
      double              qmu95_guess    = pow(ROOT::Math::gaussian_quantile(1 - m_target_CLs, 1), 2);
      int                 nrItr          = 0;
      int                 nrDamping      = 1;
      double              damping_factor = 1.0;
      double              qmu95_pre      = qmu95_guess - 10 * 2 * qmu95_guess * m_precision;
      while (fabs(qmu95_guess - qmu95_pre) > 2 * qmu95_guess * m_precision) {
         qmu95_pre = qmu95_guess;
         if (m_debugLevel == 0) {
            cout << "qmu95_guess = " << qmu95_guess << endl;
            cout << "CLs = " << calcCLs(qmu95_guess, sigma, mu) << endl;
            cout << "Derivative = " << calcDerCLs(qmu95_guess, sigma, mu) << endl;
         }

         double corr =
            damping_factor * (calcCLs(qmu95_guess, sigma, mu) - m_target_CLs) / calcDerCLs(qmu95_guess, sigma, mu);
         for (map<double, double>::iterator itr = guess_to_corr.begin(); itr != guess_to_corr.end(); itr++) {
            if (fabs(itr->first - qmu95_guess) < 2 * qmu95_guess * m_precision) {
               damping_factor *= 0.8;
               if (m_debugLevel == 0)
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

         if (m_debugLevel == 0) {
            cout << "next guess = " << qmu95_guess << endl;
            cout << "precision = " << 2 * qmu95_guess * m_precision << endl;
            cout << endl;
         }
         nrItr++;
         if (nrItr > 200) {
            cout << "Infinite loop detected in getQmu95. Please intervene." << endl;
            throw std::runtime_error("Infinite loop");
            exit(1);
         }
      }
      qmu95 = qmu95_guess;
   }

   if (qmu95 != qmu95) {
      qmu95 = getQmu95_brute(sigma, mu);
   }
   if (m_debugLevel == 0) cout << "Returning qmu95 = " << qmu95 << endl;

   return qmu95;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates CLs given a value of qmu_tilde, sigma and the parameter of interest
///
/// \param[in] qmu_tilde
/// \param[in] sigma
/// \param[in] mu
/// \param[out] CLs
///
/// Reverts to the brute-force calculation getQmu95_brute() if needed
double EXOSTATS::AsymptoticsCLsRunner::calcCLs(double qmu_tilde, double sigma, double mu)
{
   double pmu = calcPmu(qmu_tilde, sigma, mu);
   double pb  = calcPb(qmu_tilde, sigma, mu);
   if (m_debugLevel == 0) {
      cout << "pmu = " << pmu << endl;
      cout << "pb = " << pb << endl;
   }
   if (pb == 1) return 0.5;
   return pmu / (1 - pb);
}

//////////////////////////////////////////////////////////////////////////
/// Calculates p_mu
///
/// \param[in] qmu
/// \param[in] sigma
/// \param[in] mu
/// \param[out] p_mu
double EXOSTATS::AsymptoticsCLsRunner::calcPmu(double qmu, double sigma, double mu)
{
   double pmu;
   if (qmu < mu * mu / (sigma * sigma) || !m_doTilde) {
      pmu = 1 - ROOT::Math::gaussian_cdf(sqrt(qmu));
   } else {
      pmu = 1 - ROOT::Math::gaussian_cdf((qmu + mu * mu / (sigma * sigma)) / (2 * fabs(mu / sigma)));
   }
   if (m_debugLevel == 0)
      cout << "for pmu, qmu = " << qmu << ", sigma = " << sigma << ", mu = " << mu << ", pmu = " << pmu << endl;
   return pmu;
}

//////////////////////////////////////////////////////////////////////////
/// Calculates p_b
///
/// \param[in] qmu
/// \param[in] sigma
/// \param[in] mu
/// \param[out] p_b
double EXOSTATS::AsymptoticsCLsRunner::calcPb(double qmu, double sigma, double mu)
{
   if (qmu < mu * mu / (sigma * sigma) || !m_doTilde) {
      return 1 - ROOT::Math::gaussian_cdf(fabs(mu / sigma) - sqrt(qmu));
   } else {
      return 1 - ROOT::Math::gaussian_cdf((mu * mu / (sigma * sigma) - qmu) / (2 * fabs(mu / sigma)));
   }
}

//////////////////////////////////////////////////////////////////////////
/// Internal use
///
/// \param[in] qmu
/// \param[in] sigma
/// \param[in] mu
/// \param[out] dpmu_dq
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

//////////////////////////////////////////////////////////////////////////
/// Activates the better implementation of the calculation of limit bands
///
/// Improves bands by using a more appropriate Asimov dataset for those points. Recommended value is \c kTRUE
void EXOSTATS::AsymptoticsCLsRunner::setBetterBands(Bool_t value)
{
   m_betterBands = value;
}

//////////////////////////////////////////////////////////////////////////
/// Activates the better implementation of the calculation of limit bands also for the negative bands
///
/// Recommended value is \c kFALSE
void EXOSTATS::AsymptoticsCLsRunner::setBetterNegativeBands(Bool_t value)
{
   m_betterNegativeBands = value;
}

//////////////////////////////////////////////////////////////////////////
/// Profiles the Asimov for negative bands at zero
///
/// Recommended value is \c kFALSE
void EXOSTATS::AsymptoticsCLsRunner::setProfileNegativeAtZero(Bool_t value)
{
   m_profileNegativeAtZero = value;
}

//////////////////////////////////////////////////////////////////////////
/// Sets the name of the default minimizer
///
/// E.g. Minuit/Minuit2
void EXOSTATS::AsymptoticsCLsRunner::setDefaultMinimizer(TString value)
{
   m_defaultMinimizer = value;
}

//////////////////////////////////////////////////////////////////////////
/// Sets the Minuit print level
void EXOSTATS::AsymptoticsCLsRunner::setDefaultMinimizerPrintLevel(Int_t value)
{
   m_defaultPrintLevel = value;
}

//////////////////////////////////////////////////////////////////////////
/// Sets the Minuit minimization strategy.
///
/// 0 = fastest, least robust. 2 = slowest, most robust
void EXOSTATS::AsymptoticsCLsRunner::setDefaultMinimizerStrategy(Int_t value)
{
   m_defaultStrategy = value;
}

//////////////////////////////////////////////////////////////////////////
/// Suppresses all RooFit warnings below FATAL
void EXOSTATS::AsymptoticsCLsRunner::setRooFitWarningSuppression(Bool_t value)
{
   m_killBelowFatal = value;
}

//////////////////////////////////////////////////////////////////////////
/// Specifies the analysis is blinded
///
/// If set to \c kTRUE, observed limits are deactivated and the expected limits are not performed
/// by profiling nuisance parameters on the observed dataset
void EXOSTATS::AsymptoticsCLsRunner::setBlind(Bool_t value)
{
   m_doBlind             = value;
   m_conditionalExpected = m_conditionalExpected && !m_doBlind;
   m_doObs               = m_doObs && !m_doBlind;
}

//////////////////////////////////////////////////////////////////////////
/// Specifies the profiling mode for Asimov data
///
/// If set to \c kTRUE, conditional maximum-likelihood estimates of all parameters are used (i.e. they
/// are fitted to data in a background-only fit). If set to \c kFALSE, nominal parameter values are
/// used instead.
void EXOSTATS::AsymptoticsCLsRunner::setConditionalExpected(Bool_t value)
{
   m_conditionalExpected = value && !m_doBlind;
}

//////////////////////////////////////////////////////////////////////////
/// If set to \c kTRUE, binds the parameter of interest at zero and computes the q_mu_tilde calculation
void EXOSTATS::AsymptoticsCLsRunner::setTilde(Bool_t value)
{
   m_doTilde = value;
}

//////////////////////////////////////////////////////////////////////////
/// If set to kTRUE, compute only the expected limits.
void EXOSTATS::AsymptoticsCLsRunner::setExpected(Bool_t value)
{
   if (value) {
      m_doBlind             = 1;               // in case your analysis is blinded
      m_conditionalExpected = 1 && !m_doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
      m_doExp               = 1;               // compute expected limit
      m_doObs               = 1 && !m_doBlind; // compute observed limit
   } else {
      m_doBlind             = 0;               // in case your analysis is blinded
      m_conditionalExpected = 1 && !m_doBlind; // Profiling mode for Asimov data: 1 = conditional MLEs, 0 = nominal MLEs
      m_doExp               = 1;               // compute expected limit
      m_doObs               = 1 && !m_doBlind; // compute observed limit
   }
}

//////////////////////////////////////////////////////////////////////////
/// Compute the observed limit. 
///
/// Must be called \b after the setBlind() has been called.
void EXOSTATS::AsymptoticsCLsRunner::setObserved(Bool_t value)
{
   m_doObs = value && !m_doBlind;
}

//////////////////////////////////////////////////////////////////////////
/// Compute the limit after signal injection. 
///
/// See the method setInjectionStrength().
void EXOSTATS::AsymptoticsCLsRunner::setInjection(Bool_t value)
{
   m_doInj = value;
}

//////////////////////////////////////////////////////////////////////////
/// Configures the signal injection test
///
/// If and only if setInjection(kTRUE) was called, a signal injection test
/// will be performed, with the specified value of the parameter of interest
void EXOSTATS::AsymptoticsCLsRunner::setInjectionStrength(Double_t value)
{
   m_muInjection = value;
}

//////////////////////////////////////////////////////////////////////////
/// Precision to be used in the parameter of interest calculation iteration
void EXOSTATS::AsymptoticsCLsRunner::setPrecision(Double_t value)
{
   m_precision = value;
}

//////////////////////////////////////////////////////////////////////////
/// Debug level: 0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent
void EXOSTATS::AsymptoticsCLsRunner::setDebugLevel(Int_t value)
{
   m_debugLevel = value;
}

//////////////////////////////////////////////////////////////////////////
/// Experimental: extrapolate best fit nuisance parameters based on previous fit results
void EXOSTATS::AsymptoticsCLsRunner::setUsePredictiveFit(Bool_t value)
{
   m_usePredictiveFit = value;
}

//////////////////////////////////////////////////////////////////////////
/// Experimental: extrapolate sigma based on previous fits
void EXOSTATS::AsymptoticsCLsRunner::setExtrapolateSigma(Bool_t value)
{
   m_extrapolateSigma = value;
}

//////////////////////////////////////////////////////////////////////////
/// Set number of retries of the negative log-likelihood minimization before giving up
void EXOSTATS::AsymptoticsCLsRunner::setMaxRetries(Int_t value)
{
   m_maxRetries = value;
}

//////////////////////////////////////////////////////////////////////////
/// Perform the p-value calculation
void EXOSTATS::AsymptoticsCLsRunner::setCalculatePvalues(Bool_t value)
{
   m_doPvals = value;
}

//////////////////////////////////////////////////////////////////////////
/// Set the number of CPUs to be used for the parallelisation of the likelihood calculation
void EXOSTATS::AsymptoticsCLsRunner::setNumCPU(Int_t value)
{
   m_NumCPU = value;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getBetterBands()
{
   return m_betterBands;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getBetterNegativeBands()
{
   return m_betterNegativeBands;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getProfileNegativeAtZero()
{
   return m_profileNegativeAtZero;
}

TString EXOSTATS::AsymptoticsCLsRunner::getDefaultMinimizer()
{
   return m_defaultMinimizer;
}

Int_t EXOSTATS::AsymptoticsCLsRunner::getDefaultMinimizerPrintLevel()
{
   return m_defaultPrintLevel;
}

Int_t EXOSTATS::AsymptoticsCLsRunner::getDefaultMinimizerStrategy()
{
   return m_defaultStrategy;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getRooFitWarningSuppression()
{
   return m_killBelowFatal;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getBlind()
{
   return m_doBlind;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getConditionalExpected()
{
   return m_conditionalExpected;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getTilde()
{
   return m_doTilde;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getExpected()
{
   return m_doExp;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getObserved()
{
   return m_doObs;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getInjection()
{
   return m_doInj;
}

Double_t EXOSTATS::AsymptoticsCLsRunner::getInjectionStrength()
{
   return m_muInjection;
}

Double_t EXOSTATS::AsymptoticsCLsRunner::getPrecision()
{
   return m_precision;
}

Int_t EXOSTATS::AsymptoticsCLsRunner::getDebugLevel()
{
   return m_debugLevel;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getUsePredictiveFit()
{
   return m_usePredictiveFit;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getExtrapolateSigma()
{
   return m_extrapolateSigma;
}

Int_t EXOSTATS::AsymptoticsCLsRunner::getMaxRetries()
{
   return m_maxRetries;
}

Bool_t EXOSTATS::AsymptoticsCLsRunner::getCalculatePvalues()
{
   return m_doPvals;
}

Int_t EXOSTATS::AsymptoticsCLsRunner::getNumCPU()
{
   return m_NumCPU;
}

void EXOSTATS::AsymptoticsCLsRunner::printOptionValues()
{
   cout << "Settings:" << endl;
   cout << "  - betterBands is set to " << m_betterBands << endl;
   cout << "  - betterNegativeBands is set to " << m_betterNegativeBands << endl;
   cout << "  - profileNegativeAtZero is set to " << m_profileNegativeAtZero << endl;
   cout << "  - defaultMinimizer is set to " << m_defaultMinimizer << endl;
   cout << "  - defaultPrintLevel is set to " << m_defaultPrintLevel << endl;
   cout << "  - defaultStrategy is set to " << m_defaultStrategy << endl;
   cout << "  - killBelowFatal is set to " << m_killBelowFatal << endl;
   cout << "  - doBlind is set to " << m_doBlind << endl;
   cout << "  - conditionalExpected is set to " << m_conditionalExpected << endl;
   cout << "  - doTilde is set to " << m_doTilde << endl;
   cout << "  - doExp is set to " << m_doExp << endl;
   cout << "  - doObs is set to " << m_doObs << endl;
   cout << "  - doInj is set to " << m_doInj << endl;
   cout << "  - muInjection is set to " << m_muInjection << endl;
   cout << "  - precision is set to " << m_precision << endl;
   cout << "  - debugLevel is set to " << m_debugLevel << endl;
   cout << "  - usePredictiveFit is set to " << m_usePredictiveFit << endl;
   cout << "  - extrapolateSigma is set to " << m_extrapolateSigma << endl;
   cout << "  - maxRetries is set to " << m_maxRetries << endl;
   cout << "  - doPvals is set to " << m_doPvals << endl;
   cout << "  - NumCPU is set to " << m_NumCPU << endl;
   cout << "  - target_CLs is set to " << m_target_CLs << endl;
   cout << endl;
}
