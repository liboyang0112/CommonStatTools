#include <TFile.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/RooStatsUtils.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooRealSumPdf.h>
#include <RooRealVar.h>
#include <RooProduct.h>
#include <TPRegexp.h>
#include <RooFitResult.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TObjString.h>
#include <RooGaussian.h>
#include <iostream>
#include <sstream>

#include "RooExpandedFitResult.h"
#include "HistFactoryInspector.h"

using namespace std;

////////////////////////////////////
/// A function to tokenize a string
///
/// \param[in] line the string
/// \param[in] delim its separator
/// \param[out] vtokens a vector of tokens
vector<TString> getTokens(TString line, TString delim)
{
   vector<TString> vtokens;
   TObjArray *          tokens = TString(line).Tokenize(delim); // delimiters
   if (tokens->GetEntriesFast()) {
      TIter       iString(tokens);
      TObjString *os = 0;
      while ((os = (TObjString *)iString())) {
         vtokens.push_back(os->GetString().Data());
      }
   }
   delete tokens;

   return vtokens;
}

////////////////////////////////////
/// A copy of the python \c join method
///
/// \param[in] separator the separator to add
/// \param[in] vec a vector of strings
/// \param[out] result a string containing all vector elements separated by a separator
TString join(TString separator, vector<TString> vec)
{
   TString result("");
   for (auto item : vec) {
      if (result == "")
         result += item;
      else
         result += separator + item;
   }

   return result;
}

////////////////////////////////////
/// Convert a double into a string with fixed number of digits after the decimal point
///
/// \param[in] x the number
/// \param[in] decDigits the number of digits
/// \param[out] ss a string
TString prd(const double x, const int decDigits)
{
   stringstream ss;
   ss << fixed;
   ss.precision(decDigits); // set # places after decimal
   ss << x;
   return ss.str();
}

EXOSTATS::HistFactoryInspector::HistFactoryInspector()
{
   m_debugLevel      = 2;
   m_inputFile       = "";
   m_workspaceName   = "";
   m_modelConfigName = "";
   m_dataName        = "";
   m_rangeName       = "";
   m_file            = nullptr;
   m_w               = nullptr;
   m_mc              = nullptr;
   m_simPdf          = nullptr;
   m_cat             = nullptr;
   m_prefitSnap      = "DUMMY";
   m_postfitSnap     = "DUMMY";
}

/// Debug level: 0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent
void EXOSTATS::HistFactoryInspector::setDebugLevel(Int_t level)
{
   m_debugLevel = level;
}

/// \param[in] inputFile name of the input file containing the workspace
/// \param[in] workspaceName name of the workspace
/// \param[in] modelConfigName name of the ModelConfig object to be retrieved from the workspace
/// \param[in] dataName name of the dataset to be fitted
/// \param[in] rangeName name of the observable range to restrict calculation to (must be defined in the workspace)
void EXOSTATS::HistFactoryInspector::setInput(const char *inputFile, const char *workspaceName,
                                              const char *modelConfigName, const char *dataName, TString rangeName)
{
   m_inputFile       = inputFile;
   m_workspaceName   = workspaceName;
   m_modelConfigName = modelConfigName;
   m_dataName        = dataName;
   m_rangeName       = rangeName;

   m_file        = new TFile(m_inputFile);
   m_w           = dynamic_cast<RooWorkspace *>(m_file->Get(m_workspaceName));
   m_mc          = dynamic_cast<RooStats::ModelConfig *>(m_w->obj(m_modelConfigName));
   m_simPdf      = dynamic_cast<RooSimultaneous *>(m_mc->GetPdf());
   m_cat         = m_w->cat(m_simPdf->indexCat().GetName()); // we need non-const access
   m_prefitSnap  = "myPrefitSnap";
   m_postfitSnap = "myPostfitSnap";

   // make snapshot of things before we do anything
   m_w->saveSnapshot(m_prefitSnap, *m_mc->GetPdf()->getParameters(*m_w->data(m_dataName)));
}

/// \param[in] regions vector of names of regions (channels) to compute yields/impacts in
void EXOSTATS::HistFactoryInspector::setEvalRegions(vector<TString> regions)
{
   m_evalRegions = regions;
   retrieveSampleNames();
}

/// \param[in] regions comma-separated list of names of regions (channels) to compute yields/impacts in
void EXOSTATS::HistFactoryInspector::setEvalRegions(TString regions)
{
   setEvalRegions(getTokens(regions, ","));
}

/// \param[in] regions vector of names of regions (channels) where the fit must be performed
void EXOSTATS::HistFactoryInspector::setFitRegions(vector<TString> regions)
{
   m_fitRegions = regions;
}

/// \param[in] regions comma-separated list of names of regions (channels) where the fit must be performed
void EXOSTATS::HistFactoryInspector::setFitRegions(TString regions)
{
   setFitRegions(getTokens(regions, ","));
}

//////////////////////////////////////////////////////////////////////////////
/// \param[in] asymErrors activate error calculation using asymmetric errors
/// \param[out] pair pre- and post-fit yield tables
///
/// The yield tables for all samples in activated evalRegions before and after fit are also returned. For example:
/// \code
/// auto result = hf.getYields();
/// auto prefitTable = result.first;
/// auto postfitTable = result.second;
/// for (auto kv: prefitTable) {
///   auto region = kv.first;
///   for (auto kv2: kv.second) {
///     auto sample = kv2.first;
///     auto yield = kv2.second.first;
///     auto error = kv2.second.secondo;
///     cout << "region " << region << " sample " << sample << ": yield " << yield << " +/- " << error << endl;
///     // in other workds, prefitTable["SR"]["Zmumu"].first is the Zmumu yield before fit, and .second is the error on
///     this number
///   }
/// }
/// \endcode
/// Tip: use these returned objects for fancy output formatting / usage of function output in ancillary code.
EXOSTATS::YieldTable EXOSTATS::HistFactoryInspector::getYields(Bool_t asymErrors)
{
   // prefit
   stringstream myCout; // let's print everything at the end of the job, to avoid being flooded by RooFit printouts
   myCout << "\n\n\nPRE-FIT YIELDS\n*****************\n\n";
   m_w->loadSnapshot(m_prefitSnap); // crucial!

   map<TString, std::map<TString, RooFormulaVar *>> RFV_map;

   YieldTable result;
   result.first = YieldTableElement();
   for (auto kv : m_samples) {
      auto reg = kv.first;
      myCout << "region: " << reg << endl;
      for (auto sample : kv.second) {
         myCout << "   - " << sample << ": ";
         vector<TString> vec = {sample}; // we want yields for each single sample :)
         RFV_map[reg][sample]     = retrieveYieldRFV(reg, vec);
         auto yieldRFV            = RFV_map[reg][sample];

         // retrieve floating parameters and set their errors to initial sensible values
         auto floatParList = getFloatParList(*m_simPdf, *m_mc->GetObservables());
         resetError(floatParList);

         RooExpandedFitResult emptyFitResult(getFloatParList(*m_simPdf, *m_mc->GetObservables()));
         const Double_t       rfv_val = yieldRFV->getVal();
         myCout << rfv_val << " +/- ";
         // const Double_t rfv_err = yieldRFV->getPropagatedError(emptyFitResult);
         const Double_t rfv_err = getPropagatedError(yieldRFV, emptyFitResult, asymErrors);
         myCout << rfv_err << endl;

         result.first[reg][sample].first  = rfv_val;
         result.first[reg][sample].second = rfv_err;
      }
   }

   // fit
   m_w->loadSnapshot(m_prefitSnap);
   RooFitResult *fitResult = fitPdfInRegions(m_fitRegions, kTRUE, kTRUE);

   // postfit
   myCout << "\n\n\nPOST-FIT YIELDS\n*****************\n\n";
   result.second = YieldTableElement();
   for (auto kv : m_samples) {
      auto reg = kv.first;
      myCout << "region: " << reg << endl;
      for (auto sample : kv.second) {
         myCout << "   - " << sample << ": ";
         auto yieldRFV = RFV_map[reg][sample]; // we re-use the one created for the pre-fit

         const Double_t rfv_val = yieldRFV->getVal();
         myCout << rfv_val << " +/- ";
         // const Double_t rfv_err = yieldRFV->getPropagatedError(*fitResult);
         const Double_t rfv_err = getPropagatedError(yieldRFV, *fitResult, asymErrors);
         myCout << rfv_err << endl;

         result.first[reg][sample].first  = rfv_val;
         result.first[reg][sample].second = rfv_err;
      }
   }

   // garbage collection
   for (auto kv : RFV_map) {
      for (auto kv2 : kv.second) {
         delete kv2.second;
      }
   }

   cout << myCout.str() << std::endl;

   return EXOSTATS::YieldTable();
}

EXOSTATS::ImpactTable EXOSTATS::HistFactoryInspector::getImpacts(vector<TString> samples)
{
   // prefit
   stringstream myCout; // let's print everything at the end of the job, to avoid being flooded by RooFit printouts
   myCout << "\n\n\nPRE-FIT IMPACTS\n*****************\n\n";
   m_w->loadSnapshot(m_prefitSnap); // crucial!

   map<TString, RooFormulaVar *> RFV_map; // key: region

   ImpactTable result;
   result.first = ImpactTableElement();
   map<TString, TString> summedSamples;
   for (auto kv : m_samples) {
      auto reg = kv.first;
      myCout << "region: " << reg << endl;

      // we will sum up all samples, among those requested, which are present in this regions
      vector<TString> vec;
      for (auto sample : kv.second) {
         if (find(samples.begin(), samples.end(), sample) != samples.end()) {
            vec.push_back(sample);

            if (summedSamples[reg] == "")
               summedSamples[reg] += sample;
            else
               summedSamples[reg] += " + " + sample;
         }
      }
      myCout << "   - " << summedSamples[reg] << ": yield = ";
      RFV_map[reg]  = retrieveYieldRFV(reg, vec);
      auto yieldRFV = RFV_map[reg];

      const Double_t rfv_val = yieldRFV->getVal();
      myCout << rfv_val << endl; // we don't retrieve the error

      for (auto np : getFreeParameters()) {
         auto var                     = getYieldUpDown(np, yieldRFV, kFALSE, kFALSE, kFALSE);
         result.first[reg][np].first  = var.first / rfv_val - 1;
         result.first[reg][np].second = var.second / rfv_val - 1;
         myCout << "        - " << np << ": (up, down) = (" << prd(result.first[reg][np].first * 100, 2) << "%, "
                << prd(result.first[reg][np].second * 100, 2) << "%)" << endl;
      }
   }

   // fit
   m_w->loadSnapshot(m_prefitSnap);
   RooFitResult *fitResult = fitPdfInRegions(m_fitRegions, kTRUE, kTRUE);
   if (m_debugLevel <= 1) fitResult->Print();
   m_w->saveSnapshot(m_postfitSnap,
                     *m_mc->GetPdf()->getParameters(*m_w->data(
                        m_dataName))); // we save it and reload it before evaluating all ingredients of post-fit impacts

   // postfit
   myCout << "\n\n\nPOST-FIT IMPACTS\n*****************\n\n";
   result.second = ImpactTableElement();
   for (auto kv : m_samples) {
      m_w->loadSnapshot(m_postfitSnap);
      auto reg = kv.first;
      myCout << "region: " << reg << endl;
      myCout << "   - " << summedSamples[reg] << ": yield = ";
      auto yieldRFV = RFV_map[reg];

      const Double_t rfv_val = yieldRFV->getVal();
      myCout << rfv_val << " +/- ";
      const Double_t rfv_err = getPropagatedError(yieldRFV, *fitResult, kTRUE);
      myCout << rfv_err << endl;

      for (auto np : getFreeParameters()) {
         m_w->loadSnapshot(m_postfitSnap);
         auto var                     = getYieldUpDown(np, yieldRFV, kFALSE, kTRUE, kFALSE);
         result.first[reg][np].first  = var.first / rfv_val - 1;
         result.first[reg][np].second = var.second / rfv_val - 1;
         myCout << "        - " << np << ": (up, down) = (" << prd(result.first[reg][np].first * 100, 2) << "%, "
                << prd(result.first[reg][np].second * 100, 2) << "%)" << endl;
      }
   }

   // garbage collection
   for (auto kv : RFV_map) {
      delete kv.second;
   }

   cout << myCout.str() << std::endl;

   return result;
}

EXOSTATS::ImpactTable EXOSTATS::HistFactoryInspector::getImpacts(TString samples)
{
   return getImpacts(getTokens(samples, ","));
}

void EXOSTATS::HistFactoryInspector::retrieveSampleNames()
{
   for (auto reg : m_evalRegions) retrieveSampleNames(reg);
}

void EXOSTATS::HistFactoryInspector::retrieveSampleNames(TString region)
{
   m_samples[region].clear();

   TPMERegexp re("L_x_([a-zA-Z0-9_]+)_" + region + "_.*"); // hardcoded in HistFactory

   retrieveRooProductNames(region);

   for (auto prodName : m_products[region]) {
      const Int_t matched = re.Match(prodName);
      if (matched != 2) throw runtime_error("Unable to parse RooProduct name");
      m_samples[region].push_back(re[1]);
   }
}

////////////////////////////////////
/// Retrieves the names of all p.d.f. components (samples) in a given region
///
/// \param[in] region name of the region
///
/// The rough structure of an HistFactory workspace is
/// \code
///   RooSimultaneous::simPdf
///     RooProdPdf::model_REGION1
///       RooGaussian::lumiConstraint
///       RooGaussian::alpha_SYST1Constraint (for overallsys)
///       RooPoisson::gamma_stat_REGION1_bin_BIN1_constraint (for stat uncertainty)
///       RooRealSumPdf::REGION1_model (sum of RooProduct objects, each of them weighted by [the same] bin width)
///         RooProduct::L_x_SAMPLE1_REGION1_overallSyst_x_Exp (term for sample 1, including uncertainties [RooHistFunc
///         and FlexibleInterpVar or PieceWiseInterpolation])
/// \endcode
/// Note that NormFactors are included in the L_x_blabla term (e.g. this term changes integral when SigXsecOverSM is
/// changed).
///
/// Here we simply want to retrieve the list of RooProducts (L_x_blabla) associated to a given region (i.e. to the
/// RooRealSumPdf representing that region).
void EXOSTATS::HistFactoryInspector::retrieveRooProductNames(TString region)
{
   m_products[region].clear();

   const TString  rrsPdfName = TString::Format("%s_model", region.Data()); // hardcoded in HistFactory
   RooRealSumPdf *rrsPdf = dynamic_cast<RooRealSumPdf *>(m_simPdf->getPdf(region)->getComponents()->find(rrsPdfName));

   vector<TString> result;

   TIterator *itr = rrsPdf->funcList().createIterator();

   TObject *comp = itr->Next();

   while (comp) {
      result.push_back(comp->GetName());
      comp = itr->Next();
   }

   m_products[region] = result;
}

/// Create a RooFormulaVar containing the number of events associated to a sum of RooProducts
///
/// \param[in] region channel (region) to consider
/// \param[in] components list of the components (samples) to be summed up
/// \param[out] form_frac pointer to the output RooFormulaVar
///
/// For a given list of samples in a region ('components'), we create a \c RooRealSumPdf representing
/// the sum of the corresponding components, and return its integral in the form of a RooFormulaVar
/// so that getVal() will tell us the number of events associated to this sum.
/// A range for the integral can also be specified
RooFormulaVar *EXOSTATS::HistFactoryInspector::retrieveYieldRFV(TString region, vector<TString> components)
{
   if (components.size() < 1) {
      throw runtime_error("component list is empty");
   }

   RooAbsPdf *regPdf = m_simPdf->getPdf(region);

   RooRealVar *obs =
      dynamic_cast<RooRealVar *>(m_mc->GetObservables()->find("obs_x_" + region)); // name hardcoded in HistFactory
   RooRealVar *binWidth = dynamic_cast<RooRealVar *>(
      regPdf->getVariables()->find("binWidth_obs_x_" + region + "_0")); // name hardcoded in HistFactory

   // retrieve components
   RooArgList compFuncList;
   RooArgList compCoefList;

   vector<TString> available = m_products[region];

   for (auto avail : available) {
      for (auto wanted : components) {
         const TString target = "_" + wanted + "_";

         if (avail.Contains(target)) {
            compFuncList.add(*(dynamic_cast<RooProduct *>(m_w->obj(avail))));
            compCoefList.add(*binWidth);
         }
      }
   }

   if (compFuncList.getSize() == 0 || compCoefList.getSize() == 0 || compFuncList.getSize() != compCoefList.getSize()) {
      throw runtime_error("something went wrong when fetching components");
      // return 0;
   }

   // create RRSPdf and integral
   const TString compList = join("_", components);

   const TString  compRRSname = "RRS_region_" + region + "_" + compList;
   RooRealSumPdf *compRRS     = new RooRealSumPdf(compRRSname, compRRSname, compFuncList, compCoefList);

   RooAbsReal *compFunc(nullptr);
   if (m_rangeName != "")
      compFunc = compRRS->createIntegral(RooArgSet(*obs), m_rangeName);
   else
      compFunc = compRRS->createIntegral(RooArgSet(*obs));

   // create RooFormulaVar

   TString compRFVname = "form_frac_region_" + region + "_" + compList;
   if (m_rangeName != "") compRFVname += "_" + m_rangeName;

   RooFormulaVar *form_frac = new RooFormulaVar("form_fracError", "@0", RooArgList(*compFunc));
   form_frac->SetNameTitle(compRFVname, compRFVname);

   return form_frac;
}

/// Sum up components in different regions
///
/// \param[in] regions channels (regions) to consider
/// \param[in] components list of the components (samples) to be summed up
/// \param[out] form_frac pointer to the output RooFormulaVar
///
/// See the single-region version of retrieveYieldRFV() for more details
RooFormulaVar *EXOSTATS::HistFactoryInspector::retrieveYieldRFV(vector<TString> regions,
                                                                vector<TString> components)
{
   RooArgList list;
   TString    form      = "";
   int        count     = 0;
   TString    finalName = "";
   for (TString region : regions) {
      auto RFV_for_this_region = retrieveYieldRFV(region, components);

      if (form.Length() == 0) {
         form += "@";
         finalName += RFV_for_this_region->GetName();
      } else {
         form += "+ @";
         finalName += "_plus_" + TString(RFV_for_this_region->GetName());
      }
      list.add(*RFV_for_this_region);

      form += count;
      count++;
   }

   auto result = new RooFormulaVar(finalName, form, list);

   return result;
}

/// Fits the pdf of the HistFactory workspace only in a subset of the available channels (regions)
///
/// \param[in] w pointer to the input RooWorkspace
/// \param[in] dataName name of the dataset to fit
/// \param[in] regions vector of names of channels (regions) to be considered in the fit
/// \param[in] saveResult if set to \c kTRUE, the RooFitResult is returned
/// \param[in] doMinos if set to \c kTRUE, Minos is run (much slower!)
/// \param[out] fitResult the fit result, if \c saveResult is \c kTRUE
///
/// To do the job, this function defines a new simultaneous PDF and a new dataset, including only fit regions.
RooFitResult *EXOSTATS::HistFactoryInspector::fitPdfInRegions(vector<TString> regions, Bool_t saveResult,
                                                              Bool_t doMinos)
{
   // TODO: use the better fitting technique used by HistFitter's Util::FitPdf
   // (in https://svnweb.cern.ch/trac/atlasphys/browser/Physics/SUSY/Analyses/HistFitter/trunk/src/Utils.cxx)

   // to do so, a temporary PDF and a temporary dataset have to be built
   RooAbsData *dataFull = m_w->data(m_dataName);

   RooSimultaneous *pdf  = m_simPdf;
   RooDataSet *     data = dynamic_cast<RooDataSet *>(dataFull);

   // determine useful terms
   vector<RooAbsPdf *>  pdfVec;
   vector<RooDataSet *> dataVec;

   for (auto region : regions) {
      if (m_cat->setLabel(region, kTRUE)) {
         throw runtime_error("Unknown region found");
      } else {
         RooDataSet *regionData = dynamic_cast<RooDataSet *>(
            data->reduce(TString::Format("%s==%s::%s", m_cat->GetName(), m_cat->GetName(), region.Data())));
         RooAbsPdf *regionPdf = pdf->getPdf(region);

         dataVec.push_back(regionData);
         pdfVec.push_back(regionPdf);
      }
   }

   if (dataVec.size() == 0 || pdfVec.size() == 0 || dataVec.size() != pdfVec.size() || dataVec.size() != regions.size())
      throw runtime_error("Error in specified regions");

   // merge terms
   const TString nickname = join("_", regions);
   data                   = dynamic_cast<RooDataSet *>(dataVec[0]->Clone("obsDataReduced_" + nickname));
   for (UInt_t i = 1; i < dataVec.size(); i++) data->append(*dataVec[i]);

   pdf = new RooSimultaneous("simPdfReduced_" + nickname, "simultaneous pdf reduced to regions " + nickname, *m_cat);
   for (UInt_t i = 0; i < pdfVec.size(); i++) pdf->addPdf(*pdfVec[i], regions[i]);

   // perform the fit (this is the part which should possibly be improved)

   RooArgSet *allParams = pdf->getParameters(data);
   RooStats::RemoveConstantParameters(allParams);
   const RooArgSet *globObs   = m_mc->GetGlobalObservables();
   RooFitResult *   fitResult = nullptr;
   if (saveResult)
      fitResult = pdf->fitTo(*data, RooFit::GlobalObservables(*globObs), RooFit::Minos(doMinos), RooFit::Save());
   else
      pdf->fitTo(*data, RooFit::GlobalObservables(*globObs), RooFit::Minos(doMinos));
   return fitResult;
}

/// Evaluate the effect on a \c RooFormulaVar of changing a nuisance parameter
///
/// \param[in] param name of the nuisance parameter (NP)
/// \param[in] w pointer to the input RooWorspace
/// \param[in] impact pointer to the RooFormulaVar used to calculate the impact
/// \param[in] useErrorVal if \c kTRUE, will use fitted errors as plus and minus NP variations; otherwise, will use +/-
/// 1 \param[out] result pair of impact of up and down variation on the provided RooFormulaVar
/// \param[in] doFit if \c kTRUE, will fit the PDF after fixing the NP to up/down one sigma
/// \param[in] doMinos if \c kTRUE, will use MINOS for the fit, if activated (very slow!)
///
/// Uses a given RooFormulaVar (typically the output of getComponent) to evaluate the impact of
/// varying a nuisance parameter up or down by one sigma. The variation is done either manually
/// by \c avg +/- 1 (\c useErrorVar set to \c false) or by \c avg +/- 1 sigma (\c useErrorVar set to \c true)
///
/// Note that, for OverallSys uncertainties and fit deactivated, the output of this function
/// must be ~identical to what specified in the HistFactory's XMLs.
pair<Double_t, Double_t> EXOSTATS::HistFactoryInspector::getYieldUpDown(TString param, RooFormulaVar *yield,
                                                                             Bool_t useErrorVar, Bool_t doFit,
                                                                             Bool_t doMinos)
{
   pair<Double_t, Double_t> result;

   auto         var         = m_w->var(param);
   const Bool_t wasConstant = var->isConstant();

   // NOTE: assumes that the nominal value of the NP is meaningful already, i.e.
   // that the standard fit with all NPs was just done, if appropriate
   // that's why nom/up/down values of the NP are calculated now...
   const Double_t par_nom  = var->getVal();
   const Double_t par_up   = par_nom + ((useErrorVar) ? var->getErrorHi() : 1);
   const Double_t par_down = par_nom + ((useErrorVar) ? var->getErrorLo() : -1);

   // fit only if requested
   var->setVal(par_up);
   if (doFit) {
      var->setConstant(kTRUE);
      fitPdfInRegions(m_fitRegions, kFALSE, doMinos);
   }
   result.first = yield->getVal();

   var->setVal(par_down);
   if (doFit) {
      var->setConstant(kTRUE);
      fitPdfInRegions(m_fitRegions, kFALSE, doMinos);
   }
   result.second = yield->getVal();

   // restore original value and constantness
   var->setVal(par_nom);
   var->setConstant(wasConstant);

   return result;
}

/// Get list of free parameters of a p.d.f.
///
/// \param[out] result vector of nuisance parameter names
///
/// Gets list of free parameters (i.e. non-constant parameters)
vector<TString> EXOSTATS::HistFactoryInspector::getFreeParameters()
{
   // TODO: use RooStats::RemoveConstantParameters instead of the manual thing...
   const RooArgSet *NPs = m_mc->GetNuisanceParameters();

   TIterator *itr = NPs->createIterator();

   TObject *np = itr->Next();

   vector<TString> result;

   while (np) {
      RooAbsArg *raa = dynamic_cast<RooAbsArg *>(np);
      if (m_debugLevel <= 1)
         cout << "NP named " << raa->GetName() << " has constant=" << raa->isConstant() << std::endl;
      if (raa->isConstant() == kFALSE) result.push_back(np->GetName());
      np = itr->Next();
   }

   return result;
}

/////////////////////////////
/// Adapted from HistFitter's Util::GetPropagatedError
Double_t EXOSTATS::HistFactoryInspector::getPropagatedError(RooAbsReal *var, const RooFitResult &fitResult,
                                                            const Bool_t doAsym)
{
   if (m_debugLevel <= 1) cout << " GPP for variable = " << var->GetName() << std::endl;
   // Clone self for internal use
   RooAbsReal *cloneFunc   = var; //(RooAbsReal*) var->cloneTree();
   RooArgSet * errorParams = cloneFunc->getObservables(fitResult.floatParsFinal());
   RooArgSet * nset        = cloneFunc->getParameters(*errorParams);

   // Make list of parameter instances of cloneFunc in order of error matrix
   RooArgList        paramList;
   const RooArgList &fpf = fitResult.floatParsFinal();
   vector<int>  fpf_idx;
   for (Int_t i = 0; i < fpf.getSize(); i++) {
      RooAbsArg *par = errorParams->find(fpf[i].GetName());
      if (par) {
         if (!par->isConstant()) {
            paramList.add(*par);
            fpf_idx.push_back(i);
         }
      }
   }

   vector<Double_t> plusVar, minusVar;

   TMatrixDSym V(fitResult.covarianceMatrix());

   for (Int_t ivar = 0; ivar < paramList.getSize(); ivar++) {

      RooRealVar &rrv = (RooRealVar &)fpf[fpf_idx[ivar]];

      int newI = fpf_idx[ivar];

      Double_t cenVal = rrv.getVal();
      Double_t errHes = sqrt(V(newI, newI));

      Double_t errHi  = rrv.getErrorHi();
      Double_t errLo  = rrv.getErrorLo();
      Double_t errAvg = (TMath::Abs(errLo) + TMath::Abs(errHi)) / 2.0;

      Double_t errVal = errHes;
      if (doAsym) {
         errVal = errAvg;
      }

      if (m_debugLevel <= 1)
         cout << " GPP:  par = " << rrv.GetName() << " cenVal = " << cenVal << " errSym = " << errHes
                   << " errAvgAsym = " << errAvg << endl;

      // Make Plus variation
      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal + errVal);
      plusVar.push_back(cloneFunc->getVal(nset));

      // Make Minus variation
      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal - errVal);
      minusVar.push_back(cloneFunc->getVal(nset));

      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal);
   }

   TMatrixDSym         C(paramList.getSize());
   vector<double> errVec(paramList.getSize());
   for (int i = 0; i < paramList.getSize(); i++) {
      int newII = fpf_idx[i];
      errVec[i] = sqrt(V(newII, newII));
      for (int j = i; j < paramList.getSize(); j++) {
         int newJ = fpf_idx[j];
         C(i, j)  = V(newII, newJ) / sqrt(V(newII, newII) * V(newJ, newJ));
         C(j, i)  = C(i, j);
      }
   }

   // Make vector of variations
   TVectorD F(plusVar.size());

   for (unsigned int j = 0; j < plusVar.size(); j++) {
      F[j] = (plusVar[j] - minusVar[j]) / 2;
   }

   if (m_debugLevel < 1) {
      F.Print();
      C.Print();
   }

   // Calculate error in linear approximation 1 variations and correlation coefficient
   Double_t sum = F * (C * F);

   if (m_debugLevel >= 1) cout << " GPP : sum = " << sqrt(sum) << std::endl;

   return sqrt(sum);
}

//////////////////////////////
// adapted from HistFitter
RooArgList EXOSTATS::HistFactoryInspector::getFloatParList(const RooAbsPdf &pdf, const RooArgSet &obsSet)
{
   RooArgList floatParList;

   const RooArgSet *pars = pdf.getParameters(obsSet);
   if (pars == 0) {
      return floatParList;
   }

   TIterator *iter = pars->createIterator();
   RooAbsArg *arg;
   while ((arg = (RooAbsArg *)iter->Next())) {
      if (arg->InheritsFrom("RooRealVar") && !arg->isConstant()) {
         floatParList.add(*arg);
      }
   }
   delete iter;

   return floatParList;
}

/////////////////////////////////////////////////////////////////
/// Find the input parameter (systematic) with the given name and shift that parameter (systematic) by 1-sigma if given;
/// otherwise set error to small number
///
/// Adapted from HistFitter
void EXOSTATS::HistFactoryInspector::resetError(const RooArgList &parList, const RooArgList &vetoList)

{
   /// For the given workspace,
   /// find the input systematic with
   /// the given name and shift that
   /// systematic by 1-sigma

   if (m_debugLevel <= 1)
      cout << " starting with workspace: " << m_w->GetName() << "   parList.getSize(): " << parList.getSize()
                << "  vetoList.size() = " << vetoList.getSize() << endl;

   TIterator *iter = parList.createIterator();
   RooAbsArg *arg;
   while ((arg = (RooAbsArg *)iter->Next())) {

      string UncertaintyName;
      if (arg->InheritsFrom("RooRealVar") && !arg->isConstant()) {
         UncertaintyName = arg->GetName();
      } else {
         continue;
      }

      if (vetoList.FindObject(UncertaintyName.c_str()) != 0) {
         continue;
      }

      RooRealVar *var = m_w->var(UncertaintyName.c_str());
      if (!var) {
         if (m_debugLevel <= 2)
            cout << "Could not find variable: " << UncertaintyName << " in workspace: " << m_w->GetName() << ": "
                      << m_w << endl;
      }

      // Initialize
      double val_hi  = FLT_MAX;
      double val_low = FLT_MIN;
      double sigma   = 0.;
      bool   resetRange(false);

      if (UncertaintyName == "") {
         if (m_debugLevel <= 2) cout << "No Uncertainty Name provided" << std::endl;
         throw - 1;
      }
      // If it is a standard (gaussian) uncertainty
      else if (string(UncertaintyName).find("alpha") != string::npos) {
         // Assume the values are +1, -1
         val_hi     = 1.0;
         val_low    = -1.0;
         sigma      = 1.0;
         resetRange = true;
      }
      // If it is Lumi:
      else if (UncertaintyName == "Lumi") {
         // Get the Lumi's constraint term:
         RooGaussian *lumiConstr = (RooGaussian *)m_w->pdf("lumiConstraint");
         if (!lumiConstr) {
            if (m_debugLevel <= 2)
               cout << "Could not find m_w->pdf('lumiConstraint') "
                         << " in workspace: " << m_w->GetName() << ": " << m_w
                         << " when trying to reset error for parameter: Lumi" << endl;
            continue;
         }
         // Get the uncertainty on the Lumi:
         RooRealVar *lumiSigma = (RooRealVar *)lumiConstr->findServer(0);
         sigma                 = lumiSigma->getVal();

         RooRealVar *nominalLumi = m_w->var("nominalLumi");
         double      val_nom     = nominalLumi->getVal();

         val_hi     = val_nom + sigma;
         val_low    = val_nom - sigma;
         resetRange = true;
      }
      // If it is a stat uncertainty (gamma)
      else if (string(UncertaintyName).find("gamma") != string::npos) {

         // Get the constraint and check its type:
         RooAbsReal *constraint     = (RooAbsReal *)m_w->obj((UncertaintyName + "_constraint").c_str());
         string ConstraintType = "";
         if (constraint != 0) {
            ConstraintType = constraint->IsA()->GetName();
         }

         if (ConstraintType == "") {
            if (m_debugLevel <= 1)
               cout << "Assuming parameter :" << UncertaintyName << ": is a ShapeFactor and so unconstrained"
                         << endl;
            continue;
         } else if (ConstraintType == "RooGaussian") {
            RooAbsReal *sigmaVar = (RooAbsReal *)m_w->obj((UncertaintyName + "_sigma").c_str());
            sigma                = sigmaVar->getVal();

            // Symmetrize shifts
            val_hi     = 1 + sigma;
            val_low    = 1 - sigma;
            resetRange = true;
         } else if (ConstraintType == "RooPoisson") {
            RooAbsReal *nom_gamma     = (RooAbsReal *)m_w->obj(("nom_" + UncertaintyName).c_str());
            double      nom_gamma_val = nom_gamma->getVal();

            sigma      = 1 / TMath::Sqrt(nom_gamma_val);
            val_hi     = 1 + sigma;
            val_low    = 1 - sigma;
            resetRange = true;
         } else {
            if (m_debugLevel <= 2)
               cout << "Strange constraint type for Stat Uncertainties: " << ConstraintType << std::endl;
            throw - 1;
         }

      } // End Stat Error
      else {
         // Some unknown uncertainty
         if (m_debugLevel <= 1) {
            cout << "Couldn't identify type of uncertainty for parameter: " << UncertaintyName
                      << ". Assuming a normalization factor." << endl;
            cout << "Setting uncertainty to 0.0001 before the fit for parameter: " << UncertaintyName << std::endl;
         }
         sigma      = 0.0001;
         val_low    = var->getVal() - sigma;
         val_hi     = var->getVal() + sigma;
         resetRange = false;
      }

      var->setError(abs(sigma));
      if (resetRange) {
         double minrange = var->getMin();
         double maxrange = var->getMax();
         double newmin   = var->getVal() - 6. * sigma;
         double newmax   = var->getVal() + 6. * sigma;
         if (minrange < newmin) var->setMin(newmin);
         if (newmax < maxrange) var->setMax(newmax);
      }

      if (m_debugLevel <= 1)
         cout << "Uncertainties on parameter: " << UncertaintyName << " low: " << val_low << " high: " << val_hi
                   << " sigma: " << sigma << " min range: " << var->getMin() << " max range: " << var->getMax()
                   << endl;

      // Done
   } // end loop

   delete iter;
}
