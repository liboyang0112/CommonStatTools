#ifndef __EXOSTATS_ASYMPTOTICSCLSRUNNER_H__
#define __EXOSTATS_ASYMPTOTICSCLSRUNNER_H__

/*
Author: Aaron Armbruster
Date:   2012-05-25
Email:  armbrusa@umich.edu
Description: Script to run asymptotic CLs.

--------
00-01-00
-First version with updated bands

--------
00-01-01
-Fixed problem in asimov data creation that affected +1,2sigma bands

--------
00-01-02
-(Re)added support for non-sim pdfs (still need to be extended)
-Fixed default doFit arg of makeAsimovData
-Added better output for unresolved fit failures
-Improved retry loop for fit failures


/////////////////////
//////PREAMBLE///////
/////////////////////

The script uses an iterative process to find the crossing of qmu with the qmu95(mu/sigma) curve,
where qmu95(mu/sigma) is found assuming asymptotic formula for the distribution of the
test statistic f(qmu|mu') (arxiv 1007.1727) and of the test statistic qmu (or tilde)

The sequence is

mu_i+1 = mu_i - gamma_i*(mu_i - mu'_i)

where gamma_i is a dynamic damping factor used for convergence (nominal gamma_i = 1), and mu'_i is
determined by extrapolating the test statistic to the qmu95 curve assuming qmu is parabolic:

qmu'_i = (mu'_i - muhat)^2 / sigma_i^2 = qmu95(mu'_i / sigma_i)

where sigma_i is determined by computing qmu_i (not '):

sigma_i = (mu_i - muhat) / sqrt(qmu_i)

At the crossing qmu_N = qmu95 the assumption that qmu is a parabola goes away,
so we're not ultimately dependent on this assumption beyond its use in the asymptotic formula.

The sequence ends when the relative correction factor gamma*(mu_i - mu'_i) / mu_i is less than some
specified precision (0.005 by default)




///////////////////////////
//////AFTER RUNNING////////
///////////////////////////


The results will be printed as well as stored in a root file in the folder 'root-files/<folder>', where <folder>
is specified by you (default 'test').  The results are stored in a TTree.



//////////////////////////


This version is functionally fully consistent with the previous tag.

NOTE: The script runs significantly faster when compiled
*/

#include <string>
#include <map>
#include <Rtypes.h>

class RooNLLVar;
class RooDataSet;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}
class RooRealVar;

namespace EXOSTATS {
class AsymptoticsCLsRunner {
private:
   // band configuration
   Bool_t betterBands;
   Bool_t betterNegativeBands;
   Bool_t profileNegativeAtZero;
   // other configuration
   std::string defaultMinimizer;
   int         defaultPrintLevel;
   int         defaultStrategy;
   Bool_t      killBelowFatal;
   Bool_t      doBlind;
   Bool_t      conditionalExpected;
   Bool_t      doTilde;
   Bool_t      doExp;
   Bool_t      doObs;
   Bool_t      doInj;
   Double_t    precision;
   Bool_t      verbose;
   Bool_t      usePredictiveFit;
   Bool_t      extrapolateSigma;
   int         maxRetries;
   Bool_t      doPvals;

   // don't touch!
   std::map<RooNLLVar *, Double_t>                     map_nll_muhat;
   std::map<RooNLLVar *, Double_t>                     map_muhat;
   std::map<RooDataSet *, RooNLLVar *>                 map_data_nll;
   std::map<RooNLLVar *, std::string>                  map_snapshots;
   std::map<RooNLLVar *, std::map<Double_t, Double_t>> map_nll_mu_sigma;
   RooWorkspace *                                      w;
   RooStats::ModelConfig *                             mc;
   RooDataSet *                                        data;
   RooRealVar *                                        firstPOI;
   RooNLLVar *                                         asimov_0_nll;
   RooNLLVar *                                         asimov_1_nll;
   RooNLLVar *                                         obs_nll;
   int                                                 nrMinimize;
   int                                                 direction;
   int                                                 global_status;
   Double_t                                            target_CLs;
   // range of firstPOI from ModelConfig mc
   Double_t firstPOIMax;
   Double_t firstPOIMin;

public:
   void init();
   // main
   void runAsymptoticsCLs(const char *infile, const char *workspaceName, const char *modelConfigName,
                          const char *dataName, const char *asimovDataName, std::string folder, std::string mass,
                          Double_t CL, Bool_t betterBands, Double_t mu_inj = 1);

   // for backwards compatibility
   void runAsymptoticsCLs(const char *infile, const char *workspaceName = "combWS",
                          const char *modelConfigName = "ModelConfig", const char *dataName = "combData",
                          const char *asimovDataName      = "asimovData_0",
                          const char *conditionalSnapshot = "conditionalGlobs_0",
                          const char *nominalSnapshot = "nominalGlobs", std::string folder = "test",
                          std::string mass = "130", Double_t CL = 0.95, Bool_t betterBands = false);

protected:
   Int_t    minimize(RooNLLVar *nll);
   Double_t getLimit(RooNLLVar *nll, Double_t initial_guess = 0);
   void     getLimit(RooNLLVar *nll, Double_t initial_guess, Double_t &upper_limit, Double_t &muhat);
   Double_t getSigma(RooNLLVar *nll, Double_t mu, Double_t muhat, Double_t &qmu);
   Double_t getQmu(RooNLLVar *nll, Double_t mu);
   void     getExpPvalue(Double_t &pb);
   void     getObsPvalue(Double_t mu, Double_t &pv);

   void       saveSnapshot(RooNLLVar *nll, Double_t mu);
   void       loadSnapshot(RooNLLVar *nll, Double_t mu);
   void       doPredictiveFit(RooNLLVar *nll, Double_t mu1, Double_t m2, Double_t mu);
   RooNLLVar *createNLL(RooDataSet *_data);
   Double_t   getNLL(RooNLLVar *nll);
   Double_t   findCrossing(Double_t sigma_obs, Double_t sigma, Double_t muhat);
   void       setMu(Double_t mu);
   Double_t   getQmu95_brute(Double_t sigma, Double_t mu);
   Double_t   getQmu95(Double_t sigma, Double_t mu);
   Double_t   calcCLs(Double_t qmu_tilde, Double_t sigma, Double_t mu);
   Double_t   calcPmu(Double_t qmu_tilde, Double_t sigma, Double_t mu);
   Double_t   calcPb(Double_t qmu_tilde, Double_t sigma, Double_t mu);
   Double_t   calcDerCLs(Double_t qmu, Double_t sigma, Double_t mu);

   void doExpected(Bool_t isExpected);
   void doPvalues(Bool_t calc);
   void doBetterBands(Bool_t isBetterBands);
   void doInjection(Bool_t injection);
};

} // namespace EXOSTATS

#endif
