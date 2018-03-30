/// \file
/// Tools to create Asimov datasets
#include <RooDataSet.h>
#include <RooNLLVar.h>
#include <RooStats/ModelConfig.h>
#include <RooRealVar.h>
#include <sstream>
#include <iomanip>
#include <RooSimultaneous.h>
#include <RooCategory.h>

#include "AsimovDataMaking.h"
#include "Minimization.h"

using namespace RooFit;
using namespace std;

////////////////////////////////////////////////////
/// \param[in] initial reference to the starting RooArgSet, which contains pointers to RooAbsPdf 
/// \param[in] final reference to the result RooArgSet
/// \param[in] obs reference to the RooArgSet containing p.d.f. observables
/// \param[in] nuis reference to the RooArgSet containing p.d.f. nuisance parameters
/// \param[in] counter reference to a counter of calls to the function
void EXOSTATS::unfoldConstraints(RooArgSet &initial, RooArgSet &final, RooArgSet &obs, RooArgSet &nuis, int &counter)
{
   if (counter > 50) {
      cout << "ERROR::Couldn't unfold constraints!" << endl;
      cout << "Initial: " << endl;
      initial.Print("v");
      cout << endl;
      cout << "Final: " << endl;
      final.Print("v");
      throw std::runtime_error("Infinite loop");
      exit(1);
   }
   TIterator *itr = initial.createIterator();
   RooAbsPdf *pdf;
   while ((pdf = (RooAbsPdf *)itr->Next())) {
      RooArgSet nuis_tmp = nuis;
      RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
      // if (constraint_set.getSize() > 1)
      //{
      string className(pdf->ClassName());
      if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" &&
          className != "RooPoisson" && className != "RooBifurGauss") {
         counter++;
         EXOSTATS::unfoldConstraints(constraint_set, final, obs, nuis, counter);
      } else {
         final.add(*pdf);
      }
   }
   delete itr;
}

////////////////////////////////////////////////////
/// \param[in] w pointer to the RooWorkspace containing the likelihood
/// \param[in] w name of the ModelConfig containing the model specifications
/// \param[in] doConditional profile nuisance parameters using a conditioning negative log-likelihood
/// \param[in] conditioning_mll pointer to the negative log-likelihood to be used for profiling, if activated
/// \param[in] mu_val value of the parameter of interest (POI) at which the Asimov dataset must be created
/// \param[in] mu_str pointer to a string where \c mu_val will be stored
/// \param[in] mu_prof_str pointer to a string where the POI value used for profiling will be stored
/// \param[in] mu_val_profile value of the POI to be used for profiling; if not specified, \c mu_val is used
/// \param[in] doFit perform the fit
/// \param[in] mu_injection if > 0, specifies the POI value of the injected signal
/// \param[in] debugLevel debug level (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
/// \param[out] asimovData Asimov dataset
///
/// This function creates an Asimov dataset, given a model. The user can specify both which POI value
/// should be used for generating the dataset, and which POI value (and likelihood) should be used to 
/// evaluate the nuisance parameters before the generation is performed.
///
/// If \c mu_injection is > 0, a signal will be injected with the specified signal strength instead.
/// Note that, if the workspace contains a variable called \c ATLAS_norm_muInjection, its value is set to
/// \c mu_injection and the POI is not changed.
RooDataSet *EXOSTATS::makeAsimovData(RooWorkspace *w, TString modelConfigName, Bool_t doConditional,
                                     RooNLLVar *conditioning_nll, Double_t mu_val, std::string *mu_str,
                                     std::string *mu_prof_str, Double_t mu_val_profile, Bool_t doFit,
                                     Double_t mu_injection, Int_t debugLevel)
{
   if (mu_val_profile == -999) mu_val_profile = mu_val;

   if (debugLevel >= 0)
      cout << "Creating asimov data at mu = " << mu_val << ", profiling at mu = " << mu_val_profile << endl;

   // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // int strat = ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   // int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
   // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
   // RooMinuit::SetMaxIterations(10000);
   // RooMinimizer::SetMaxFunctionCalls(10000);

   ////////////////////
   // make asimov data//
   ////////////////////
   RooStats::ModelConfig *mc      = (RooStats::ModelConfig *)w->obj(modelConfigName);
   RooAbsPdf *            combPdf = mc->GetPdf();

   int _printLevel = debugLevel;

   stringstream muStr;
   muStr << setprecision(5);
   muStr << "_" << mu_val;
   if (mu_str) *mu_str = muStr.str();

   stringstream muStrProf;
   muStrProf << setprecision(5);
   muStrProf << "_" << mu_val_profile;
   if (mu_prof_str) *mu_prof_str = muStrProf.str();

   RooRealVar *mu = (RooRealVar *)mc->GetParametersOfInterest()->first(); // w->var("mu");
   mu->setVal(mu_val);

   RooArgSet mc_obs   = *mc->GetObservables();
   RooArgSet mc_globs = *mc->GetGlobalObservables();
   RooArgSet mc_nuis  = (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet());

   // pair the nuisance parameter to the global observable
   RooArgSet  mc_nuis_tmp = mc_nuis;
   RooArgList nui_list("ordered_nuis");
   RooArgList glob_list("ordered_globs");
   RooArgSet  constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
   RooArgSet  constraint_set;
   int        counter_tmp = 0;
   EXOSTATS::unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

   TIterator *cIter = constraint_set.createIterator();
   RooAbsArg *arg;
   while ((arg = (RooAbsArg *)cIter->Next())) {
      RooAbsPdf *pdf = (RooAbsPdf *)arg;
      if (!pdf) continue;
      //     cout << "Printing pdf" << endl;
      //     pdf->Print();
      //     cout << "Done" << endl;
      TIterator * nIter   = mc_nuis.createIterator();
      RooRealVar *thisNui = NULL;
      RooAbsArg * nui_arg;
      while ((nui_arg = (RooAbsArg *)nIter->Next())) {
         if (pdf->dependsOn(*nui_arg)) {
            thisNui = (RooRealVar *)nui_arg;
            break;
         }
      }
      delete nIter;

      // RooRealVar* thisNui = (RooRealVar*)pdf->getObservables();

      // need this incase the observable isn't fundamental.
      // in this case, see which variable is dependent on the nuisance parameter and use that.
      RooArgSet *components = pdf->getComponents();
      //     cout << "\nPrinting components" << endl;
      //     components->Print();
      //     cout << "Done" << endl;
      components->remove(*pdf);
      if (components->getSize()) {
         TIterator *itr1 = components->createIterator();
         RooAbsArg *arg1;
         while ((arg1 = (RooAbsArg *)itr1->Next())) {
            TIterator *itr2 = components->createIterator();
            RooAbsArg *arg2;
            while ((arg2 = (RooAbsArg *)itr2->Next())) {
               if (arg1 == arg2) continue;
               if (arg2->dependsOn(*arg1)) {
                  components->remove(*arg1);
               }
            }
            delete itr2;
         }
         delete itr1;
      }
      if (components->getSize() > 1) {
         cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
         return NULL;
      } else if (components->getSize() == 1) {
         thisNui = (RooRealVar *)components->first();
      }

      TIterator * gIter    = mc_globs.createIterator();
      RooRealVar *thisGlob = NULL;
      RooAbsArg * glob_arg;
      while ((glob_arg = (RooAbsArg *)gIter->Next())) {
         if (pdf->dependsOn(*glob_arg)) {
            thisGlob = (RooRealVar *)glob_arg;
            break;
         }
      }
      delete gIter;

      if (!thisNui || !thisGlob) {
         cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
         // return;
         continue;
      }

      if (_printLevel >= 1)
         cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " << thisGlob->GetName()
              << ", from constraint: " << pdf->GetName() << endl;

      nui_list.add(*thisNui);
      glob_list.add(*thisGlob);

      //     cout << "\nPrinting Nui/glob" << endl;
      //     thisNui->Print();
      //     cout << "Done nui" << endl;
      //     thisGlob->Print();
      //     cout << "Done glob" << endl;
   }
   delete cIter;

   // save the snapshots of nominal parameters, but only if they're not already saved
   w->saveSnapshot("tmpGlobs", *mc->GetGlobalObservables());
   w->saveSnapshot("tmpNuis", (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet()));
   if (!w->loadSnapshot("nominalGlobs")) {
      cout << "nominalGlobs doesn't exist. Saving snapshot." << endl;
      w->saveSnapshot("nominalGlobs", *mc->GetGlobalObservables());
   } else
      w->loadSnapshot("tmpGlobs");
   if (!w->loadSnapshot("nominalNuis")) {
      cout << "nominalNuis doesn't exist. Saving snapshot." << endl;
      w->saveSnapshot("nominalNuis", (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet()));
   } else
      w->loadSnapshot("tmpNuis");

   RooArgSet nuiSet_tmp(nui_list);

   mu->setVal(mu_val_profile);
   mu->setConstant(1);
   // int status = 0;
   if (doConditional && doFit) {
      EXOSTATS::minimize(conditioning_nll);
      // cout << "Using globs for minimization" << endl;
      // mc->GetGlobalObservables()->Print("v");
      // cout << "Starting minimization.." << endl;
      // RooAbsReal* nll;
      // if (!(nll = map_data_nll[combData])) nll = combPdf->createNLL(*combData, RooFit::Constrain(nuiSet_tmp));
      // RooMinimizer minim(*nll);
      // minim.setStrategy(0);
      // minim.setPrintLevel(1);
      // status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
      // ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()); if (status != 0)
      // {
      //   cout << "Fit failed for mu = " << mu->getVal() << " with status " << status << endl;
      // }
      // cout << "Done" << endl;

      // combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
   }
   mu->setConstant(0);
   mu->setVal(mu_val);

   // loop over the nui/glob list, grab the corresponding variable from the tmp ws, and set the glob to the value of the
   // nui
   int nrNuis = nui_list.getSize();
   if (nrNuis != glob_list.getSize()) {
      cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
      return NULL;
   }

   for (int i = 0; i < nrNuis; i++) {
      RooRealVar *nui  = (RooRealVar *)nui_list.at(i);
      RooRealVar *glob = (RooRealVar *)glob_list.at(i);

      // cout << "nui: " << nui << ", glob: " << glob << endl;
      // cout << "Setting glob: " << glob->GetName() << ", which had previous val: " << glob->getVal() << ", to
      // conditional val: " << nui->getVal() << endl;

      glob->setVal(nui->getVal());
   }

   // save the snapshots of conditional parameters
   // cout << "Saving conditional snapshots" << endl;
   // cout << "Glob snapshot name = " << "conditionalGlobs"+muStrProf.str() << endl;
   // cout << "Nuis snapshot name = " << "conditionalNuis"+muStrProf.str() << endl;
   w->saveSnapshot(("conditionalGlobs" + muStrProf.str()).c_str(), *mc->GetGlobalObservables());
   w->saveSnapshot(("conditionalNuis" + muStrProf.str()).c_str(),
                   (mc->GetNuisanceParameters() ? *mc->GetNuisanceParameters() : RooArgSet()));

   if (!doConditional) {
      w->loadSnapshot("nominalGlobs");
      w->loadSnapshot("nominalNuis");
   }

   if (_printLevel >= 1) cout << "Making asimov" << endl;
   // make the asimov data (snipped from Kyle)
   mu->setVal(mu_val);

   if (mu_injection > 0) {
      RooRealVar *norm_injection = w->var("ATLAS_norm_muInjection");
      if (norm_injection) {
         norm_injection->setVal(mu_injection);
      } else {
         mu->setVal(mu_injection);
      }
   }

   int iFrame = 0;

   const char *weightName = "weightVar";
   RooArgSet   obsAndWeight;
   // cout << "adding obs" << endl;
   obsAndWeight.add(*mc->GetObservables());
   // cout << "adding weight" << endl;

   RooRealVar *weightVar = NULL;
   if (!(weightVar = w->var(weightName))) {
      w->import(*(new RooRealVar(weightName, weightName, 1, 0, 10000000)));
      weightVar = w->var(weightName);
   }
   // cout << "weightVar: " << weightVar << endl;
   obsAndWeight.add(*w->var(weightName));

   // cout << "defining set" << endl;
   w->defineSet("obsAndWeight", obsAndWeight);

   //////////////////////////////////////////////////////
   //////////////////////////////////////////////////////
   //////////////////////////////////////////////////////
   //////////////////////////////////////////////////////
   //////////////////////////////////////////////////////
   // MAKE ASIMOV DATA FOR OBSERVABLES

   // dummy var can just have one bin since it's a dummy
   // if(w->var("ATLAS_dummyX"))  w->var("ATLAS_dummyX")->setBins(1);

   // cout <<" check expectedData by category"<<endl;
   // RooDataSet* simData=NULL;
   RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(mc->GetPdf());

   RooDataSet *asimovData;
   if (!simPdf) {
      // Get pdf associated with state from simpdf
      RooAbsPdf *pdftmp = mc->GetPdf(); // simPdf->getPdf(channelCat->getLabel()) ;

      // Generate observables defined by the pdf associated with this state
      RooArgSet *obstmp = pdftmp->getObservables(*mc->GetObservables());

      if (_printLevel >= 1) {
         obstmp->Print();
      }

      asimovData = new RooDataSet(("asimovData" + muStr.str()).c_str(), ("asimovData" + muStr.str()).c_str(),
                                  RooArgSet(obsAndWeight), WeightVar(*weightVar));

      RooRealVar *thisObs        = ((RooRealVar *)obstmp->first());
      double      expectedEvents = pdftmp->expectedEvents(*obstmp);
      double      thisNorm       = 0;
      for (int jj = 0; jj < thisObs->numBins(); ++jj) {
         thisObs->setBin(jj);

         thisNorm = pdftmp->getVal(obstmp) * thisObs->getBinWidth(jj);
         if (thisNorm * expectedEvents <= 0) {
            cout << "WARNING::Detected bin with zero expected events (" << thisNorm * expectedEvents
                 << ") ! Please check your inputs. Obs = " << thisObs->GetName() << ", bin = " << jj << endl;
         }
         if (thisNorm * expectedEvents > 0 && thisNorm * expectedEvents < pow(10.0, 18))
            asimovData->add(*mc->GetObservables(), thisNorm * expectedEvents);
      }

      if (_printLevel >= 1) {
         asimovData->Print();
         cout << "sum entries " << asimovData->sumEntries() << endl;
      }
      if (asimovData->sumEntries() != asimovData->sumEntries()) {
         cout << "sum entries is nan" << endl;
         throw std::runtime_error("NaN encountered");
         exit(1);
      }

      //((RooRealVar*)obstmp->first())->Print();
      // cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;

      w->import(*asimovData);

      if (_printLevel >= 1) {
         asimovData->Print();
         cout << endl;
      }
   } else {
      map<string, RooDataSet *> asimovDataMap;

      // try fix for sim pdf
      RooCategory *channelCat =
         (RooCategory *)&simPdf
            ->indexCat(); //(RooCategory*)w->cat("master_channel");//(RooCategory*) (&simPdf->indexCat());
      //    TIterator* iter = simPdf->indexCat().typeIterator() ;
      TIterator * iter      = channelCat->typeIterator();
      RooCatType *tt        = NULL;
      int         nrIndices = 0;
      while ((tt = (RooCatType *)iter->Next())) {
         nrIndices++;
      }
      for (int i = 0; i < nrIndices; i++) {
         channelCat->setIndex(i);
         iFrame++;
         // Get pdf associated with state from simpdf
         RooAbsPdf *pdftmp = simPdf->getPdf(channelCat->getLabel());

         // Generate observables defined by the pdf associated with this state
         RooArgSet *obstmp = pdftmp->getObservables(*mc->GetObservables());

         if (_printLevel >= 1) {
            obstmp->Print();
            cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
         }

         RooDataSet *obsDataUnbinned =
            new RooDataSet(Form("combAsimovData%d", iFrame), Form("combAsimovData%d", iFrame),
                           RooArgSet(obsAndWeight, *channelCat), WeightVar(*weightVar));
         RooRealVar *thisObs        = ((RooRealVar *)obstmp->first());
         double      expectedEvents = pdftmp->expectedEvents(*obstmp);
         double      thisNorm       = 0;
         for (int jj = 0; jj < thisObs->numBins(); ++jj) {
            thisObs->setBin(jj);

            thisNorm = pdftmp->getVal(obstmp) * thisObs->getBinWidth(jj);
            if (thisNorm * expectedEvents > 0 && thisNorm * expectedEvents < pow(10.0, 18))
               obsDataUnbinned->add(*mc->GetObservables(), thisNorm * expectedEvents);
         }

         if (_printLevel >= 1) {
            obsDataUnbinned->Print();
            cout << "sum entries " << obsDataUnbinned->sumEntries() << endl;
         }
         if (obsDataUnbinned->sumEntries() != obsDataUnbinned->sumEntries()) {
            cout << "sum entries is nan" << endl;
            throw std::runtime_error("NaN encountered");
            exit(1);
         }

         // ((RooRealVar*)obstmp->first())->Print();
         // cout << "pdf: " << pdftmp->GetName() << endl;
         // cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
         // cout << "-----" << endl;

         asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned; // tempData;

         if (_printLevel >= 1) {
            cout << "channel: " << channelCat->getLabel() << ", data: ";
            obsDataUnbinned->Print();
            cout << endl;
         }
      }

      asimovData = new RooDataSet(("asimovData" + muStr.str()).c_str(), ("asimovData" + muStr.str()).c_str(),
                                  RooArgSet(obsAndWeight, *channelCat), Index(*channelCat), Import(asimovDataMap),
                                  WeightVar(*weightVar));
      w->import(*asimovData);
   }

   if (mu_injection > 0) {
      RooRealVar *norm_injection = w->var("ATLAS_norm_muInjection");
      if (norm_injection) {
         norm_injection->setVal(0);
      }
   }

   // bring us back to nominal for exporting
   // w->loadSnapshot("nominalNuis");
   w->loadSnapshot("nominalGlobs");

   // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

   return asimovData;
}
