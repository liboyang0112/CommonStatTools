#ifndef __EXOSTATS_ASYMPTOTICSCLSRUNNER_H__
#define __EXOSTATS_ASYMPTOTICSCLSRUNNER_H__

/** \class EXOSTATS::AsymptoticsCLsRunner
A class to run profile likelihood limits using CLs and asymptotic formulae

This class runs profile likelihood limits using CLs and asymptotic formulae, based on Eur.Phys.J.C71:1554,2011
(arXiv:1007.1727).
Limit setting properties can be configured (and retrieved) using setters/getters methods.

A standalone function to perform limits based on a given input is provided, run().

\author Aaron Armbruster (original code), Valerio Ippolito (adaptation into a C++ class)
*/


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
#include <TString.h>

class RooNLLVar;
class RooDataSet;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}
class RooRealVar;
class TTree;

namespace EXOSTATS {
class AsymptoticsCLsRunner {
public:
   AsymptoticsCLsRunner();
   ~AsymptoticsCLsRunner();

   void reset();

   /// standalone-user-friendly function
   void run(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
            TString paramName, float paramValue, TString workspaceTag, TString outputFolder, Double_t CL = 0.95,
            const char *asimovDataName = 0);

   /// main function
   TTree *computeLimit(RooWorkspace *workspace, const char *modelConfigName, const char *dataName, TString paramName,
                       float paramValue, Double_t CL, const char *asimovDataName);

public:
   void     setBetterBands(Bool_t value);
   void     setBetterNegativeBands(Bool_t value);
   void     setProfileNegativeAtZero(Bool_t value);
   void     setDefaultMinimizer(TString value);
   void     setDefaultMinimizerPrintLevel(Int_t value);
   void     setDefaultMinimizerStrategy(Int_t value);
   void     setRooFitWarningSuppression(Bool_t value);
   void     setBlind(Bool_t value);
   void     setConditionalExpected(Bool_t value);
   void     setTilde(Bool_t value);
   void     setExpected(Bool_t value);
   void     setObserved(Bool_t value);
   void     setInjection(Bool_t value);
   void     setInjectionStrength(Double_t value);
   void     setPrecision(Double_t value);
   void     setDebugLevel(Int_t value);
   void     setUsePredictiveFit(Bool_t value);
   void     setExtrapolateSigma(Bool_t value);
   void     setMaxRetries(Int_t value);
   void     setCalculatePvalues(Bool_t value);
   void     setNumCPU(Int_t value);
   Bool_t   getBetterBands();
   Bool_t   getBetterNegativeBands();
   Bool_t   getProfileNegativeAtZero();
   TString  getDefaultMinimizer();
   Int_t    getDefaultMinimizerPrintLevel();
   Int_t    getDefaultMinimizerStrategy();
   Bool_t   getRooFitWarningSuppression();
   Bool_t   getBlind();
   Bool_t   getConditionalExpected();
   Bool_t   getTilde();
   Bool_t   getExpected();
   Bool_t   getObserved();
   Bool_t   getInjection();
   Double_t getInjectionStrength();
   Double_t getPrecision();
   Int_t    getDebugLevel();
   Bool_t   getUsePredictiveFit();
   Bool_t   getExtrapolateSigma();
   Int_t    getMaxRetries();
   Bool_t   getCalculatePvalues();
   Int_t    getNumCPU();
   void     printOptionValues();

protected:
   Int_t    minimize(RooNLLVar *nll);
   Double_t getLimit(RooNLLVar *nll, Double_t initial_guess = 0);
   void     getLimit(RooNLLVar *nll, Double_t initial_guess, Double_t &upper_limit, Double_t &muhat,
                     std::map<TString, Float_t> &np_hat_map);
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

private:
   // band configuration
   Bool_t m_betterBands;
   Bool_t m_betterNegativeBands;
   Bool_t m_profileNegativeAtZero;
   // other configuration
   std::string m_defaultMinimizer;
   int         m_defaultPrintLevel;
   int         m_defaultStrategy;
   Bool_t      m_killBelowFatal;
   Bool_t      m_doBlind;
   Bool_t      m_conditionalExpected;
   Bool_t      m_doTilde;
   Bool_t      m_doExp;
   Bool_t      m_doObs;
   Bool_t      m_doInj;
   Double_t    m_muInjection;
   Double_t    m_precision;
   Int_t       m_debugLevel;
   Bool_t      m_usePredictiveFit;
   Bool_t      m_extrapolateSigma;
   int         m_maxRetries;
   Bool_t      m_doPvals;
   Int_t       m_NumCPU;

   // don't touch!
   std::map<RooNLLVar *, Double_t>                     m_map_nll_muhat;
   std::map<RooNLLVar *, Double_t>                     m_map_muhat;
   std::map<RooDataSet *, RooNLLVar *>                 m_map_data_nll;
   std::map<RooNLLVar *, std::string>                  m_map_snapshots;
   std::map<RooNLLVar *, std::map<Double_t, Double_t>> m_map_nll_mu_sigma;
   RooWorkspace *                                      m_w;
   RooStats::ModelConfig *                             m_mc;
   RooDataSet *                                        m_data;
   RooRealVar *                                        m_firstPOI;
   RooNLLVar *                                         m_asimov_0_nll;
   RooNLLVar *                                         m_asimov_1_nll;
   RooNLLVar *                                         m_obs_nll;
   int                                                 m_nrMinimize;
   int                                                 m_direction;
   int                                                 m_global_status;
   Double_t                                            m_target_CLs;
   // range of m_firstPOI from ModelConfig m_mc
   Double_t m_firstPOIMax;
   Double_t m_firstPOIMin;
};

} // namespace EXOSTATS

#endif
