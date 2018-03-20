#include <Math/MinimizerOptions.h>
#include <RooMinimizerFcn.h>
#include <RooMinimizer.h>
#include <RooWorkspace.h>
#include <RooNLLVar.h>
#include <RooAbsReal.h>

#include "Minimization.h"

int EXOSTATS::minimize(RooNLLVar *nll, Int_t maxRetries, RooWorkspace *w, TString mu0Snapshot, TString nominalSnapshot,
                       Int_t debugLevel)
{
   RooAbsReal *fcn = (RooAbsReal *)nll;
   return EXOSTATS::minimize(fcn, maxRetries, w, mu0Snapshot, nominalSnapshot, debugLevel);
}

int EXOSTATS::minimize(RooAbsReal *fcn, Int_t maxRetries, RooWorkspace *w, TString mu0Snapshot, TString nominalSnapshot,
                       Int_t debugLevel)
{
   static int nrItr = 0;
   if (debugLevel == 0) {
      cout << "Starting minimization" << endl;
   }

   int              printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
   RooFit::MsgLevel msglevel   = RooMsgService::instance().globalKillBelow();
   if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   int          strat      = ROOT::Math::MinimizerOptions::DefaultStrategy();
   int          save_strat = strat;
   RooMinimizer minim(*fcn);
   minim.setStrategy(strat);
   minim.setPrintLevel(printLevel);
   minim.optimizeConst(2);

   int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                               ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

   // up the strategy
   if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
   }

   if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
   }

   // cout << "status is " << status << endl;

   // //switch minuit version and try again
   if (status != 0 && status != 1) {
      string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
      string newMinType;
      if (minType == "Minuit2")
         newMinType = "Minuit";
      else
         newMinType = "Minuit2";

      cout << "Switching minuit type from " << minType << " to " << newMinType << endl;

      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
      strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
      minim.setStrategy(strat);

      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

      if (status != 0 && status != 1 && strat < 2) {
         strat++;
         cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
         minim.setStrategy(strat);
         status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                 ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      }

      if (status != 0 && status != 1 && strat < 2) {
         strat++;
         cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
         minim.setStrategy(strat);
         status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                 ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      }

      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
   }

   if (status != 0 && status != 1) {
      nrItr++;
      if (nrItr > maxRetries) {
         nrItr = 0;
         cout << "WARNING::Fit failure unresolved with status " << status << endl;
         return status;
      } else {
         if (nrItr == 0) { // retry with mu=0 snapshot
            if (w)
               w->loadSnapshot(mu0Snapshot);
            else
               cout << "WARNING: workspace not provided, unable to set mu=0 snapshot; will simply retry as is" << endl;
            return minimize(fcn);
         } else if (nrItr == 1) { // retry with nominal snapshot
            if (w)
               w->loadSnapshot(nominalSnapshot);
            else
               cout << "WARNING: workspace not provided, unable to set nominal NP snapshot; will simply retry as is"
                    << endl;
            return minimize(fcn);
         }
      }
   }

   if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);

   if (nrItr != 0) cout << "Successful fit" << endl;
   nrItr = 0;
   return status;
}
