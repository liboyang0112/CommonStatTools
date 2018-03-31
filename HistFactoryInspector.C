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
#include <iostream>

#include "RooExpandedFitResult.h"
#include "HistFactoryInspector.h"

////////////////////////////////////
/// A function to tokenize a string
///
/// \param[in] line the string
/// \param[in] delim its separator
/// \param[out] vtokens a vector of tokens
std::vector<TString> getTokens(TString line, TString delim)
{
   std::vector<TString> vtokens;
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
TString join(TString separator, std::vector<TString> vec)
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

   m_file       = new TFile(m_inputFile);
   m_w          = dynamic_cast<RooWorkspace *>(m_file->Get(m_workspaceName));
   m_mc         = dynamic_cast<RooStats::ModelConfig *>(m_w->obj(m_modelConfigName));
   m_simPdf     = dynamic_cast<RooSimultaneous *>(m_mc->GetPdf());
   m_cat        = m_w->cat(m_simPdf->indexCat().GetName()); // we need non-const access
   m_prefitSnap = "myPrefitSnap";

   // make snapshot of things before we do anything
   m_w->saveSnapshot(m_prefitSnap, *m_mc->GetPdf()->getParameters(*m_w->data(m_dataName)));
}

/// \param[in] regions vector of names of regions (channels) to compute yields/impacts in
void EXOSTATS::HistFactoryInspector::setEvalRegions(std::vector<TString> regions)
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
void EXOSTATS::HistFactoryInspector::setFitRegions(std::vector<TString> regions)
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
   std::stringstream myCout; // let's print everything at the end of the job, to avoid being flooded by RooFit printouts
   myCout << "\n\n\nPRE-FIT\n*****************\n\n";
   m_w->loadSnapshot(m_prefitSnap); // crucial!

   std::map<TString, std::map<TString, RooFormulaVar *>> RFV_map;

   YieldTable result;
   result.first = YieldTableElement();
   for (auto kv : m_samples) {
      auto reg = kv.first;
      myCout << "region: " << reg << std::endl;
      for (auto sample : kv.second) {
         myCout << "   - " << sample << ": ";
         std::vector<TString> vec = {sample}; // we want yields for each single sample :)
         RFV_map[reg][sample]     = retrieveYieldRFV(reg, vec);
         auto yieldRFV            = RFV_map[reg][sample];

         RooExpandedFitResult emptyFitResult(getFloatParList(*m_simPdf, *m_mc->GetObservables()));
         const Double_t       rfv_val = yieldRFV->getVal();
         myCout << rfv_val << " +/- ";
         // const Double_t rfv_err = yieldRFV->getPropagatedError(emptyFitResult);
         const Double_t rfv_err = getPropagatedError(yieldRFV, emptyFitResult, asymErrors);
         myCout << rfv_err << std::endl;

         result.first[reg][sample].first  = rfv_val;
         result.first[reg][sample].second = rfv_err;
      }
   }

   // fit
   m_w->loadSnapshot(m_prefitSnap);
   RooFitResult *fitResult = fitPdfInRegions(m_fitRegions, kTRUE, kTRUE);

   // postfit
   myCout << "\n\n\nPOST-FIT\n*****************\n\n";
   result.second = YieldTableElement();
   for (auto kv : m_samples) {
      auto reg = kv.first;
      myCout << "region: " << reg << std::endl;
      for (auto sample : kv.second) {
         myCout << "   - " << sample << ": ";
         auto yieldRFV = RFV_map[reg][sample]; // we re-use the one created for the pre-fit

         const Double_t rfv_val = yieldRFV->getVal();
         myCout << rfv_val << " +/- ";
         // const Double_t rfv_err = yieldRFV->getPropagatedError(*fitResult);
         const Double_t rfv_err = getPropagatedError(yieldRFV, *fitResult, asymErrors);
         myCout << rfv_err << std::endl;

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

   std::cout << myCout.str() << std::endl;

   return EXOSTATS::YieldTable();
}

EXOSTATS::ImpactTable EXOSTATS::HistFactoryInspector::getImpacts(std::vector<TString> samples) {}

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
      if (matched != 2) throw std::runtime_error(TString::Format("Unable to parse RooProduct name"));
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

   std::vector<TString> result;

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
RooFormulaVar *EXOSTATS::HistFactoryInspector::retrieveYieldRFV(TString region, std::vector<TString> components)
{
   if (components.size() < 1) {
      throw std::runtime_error("component list is empty");
   }

   RooAbsPdf *regPdf = m_simPdf->getPdf(region);

   RooAbsData *simData = m_w->data(m_dataName);
   RooAbsData *regData =
      simData->reduce(TString::Format("%s==%s::%s", m_cat->GetName(), m_cat->GetName(), region.Data()));

   RooRealVar *obs =
      dynamic_cast<RooRealVar *>(m_mc->GetObservables()->find("obs_x_" + region)); // name hardcoded in HistFactory
   RooRealVar *binWidth = dynamic_cast<RooRealVar *>(
      regPdf->getVariables()->find("binWidth_obs_x_" + region + "_0")); // name hardcoded in HistFactory

   // retrieve components
   RooArgList compFuncList;
   RooArgList compCoefList;

   std::vector<TString> available = m_products[region];

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
      throw std::runtime_error("something went wrong when fetching components");
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
RooFormulaVar *EXOSTATS::HistFactoryInspector::retrieveYieldRFV(std::vector<TString> regions,
                                                                std::vector<TString> components)
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
RooFitResult *EXOSTATS::HistFactoryInspector::fitPdfInRegions(std::vector<TString> regions, Bool_t saveResult,
                                                              Bool_t doMinos)
{
   // TODO: use the better fitting technique used by HistFitter's Util::FitPdf
   // (in https://svnweb.cern.ch/trac/atlasphys/browser/Physics/SUSY/Analyses/HistFitter/trunk/src/Utils.cxx)

   // to do so, a temporary PDF and a temporary dataset have to be built
   RooAbsData *dataFull = m_w->data(m_dataName);

   RooSimultaneous *pdf  = m_simPdf;
   RooDataSet *     data = dynamic_cast<RooDataSet *>(dataFull);

   // determine useful terms
   std::vector<RooAbsPdf *>  pdfVec;
   std::vector<RooDataSet *> dataVec;

   for (auto region : regions) {
      if (m_cat->setLabel(region, kTRUE)) {
         throw std::runtime_error("Unknown region found");
      } else {
         RooDataSet *regionData = dynamic_cast<RooDataSet *>(
            data->reduce(TString::Format("%s==%s::%s", m_cat->GetName(), m_cat->GetName(), region.Data())));
         RooAbsPdf *regionPdf = pdf->getPdf(region);

         dataVec.push_back(regionData);
         pdfVec.push_back(regionPdf);
      }
   }

   if (dataVec.size() == 0 || pdfVec.size() == 0 || dataVec.size() != pdfVec.size() || dataVec.size() != regions.size())
      throw std::runtime_error("Error in specified regions");

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

/////////////////////////////
/// Adapted from HistFitter's Util::GetPropagatedError
Double_t EXOSTATS::HistFactoryInspector::getPropagatedError(RooAbsReal *var, const RooFitResult &fitResult,
                                                            const Bool_t doAsym)
{
   if (m_debugLevel >= 1) std::cout << " GPP for variable = " << var->GetName() << std::endl;
   // Clone self for internal use
   RooAbsReal *cloneFunc   = var; //(RooAbsReal*) var->cloneTree();
   RooArgSet * errorParams = cloneFunc->getObservables(fitResult.floatParsFinal());
   RooArgSet * nset        = cloneFunc->getParameters(*errorParams);

   // Make list of parameter instances of cloneFunc in order of error matrix
   RooArgList        paramList;
   const RooArgList &fpf = fitResult.floatParsFinal();
   std::vector<int>  fpf_idx;
   for (Int_t i = 0; i < fpf.getSize(); i++) {
      RooAbsArg *par = errorParams->find(fpf[i].GetName());
      if (par) {
         if (!par->isConstant()) {
            paramList.add(*par);
            fpf_idx.push_back(i);
         }
      }
   }

   std::vector<Double_t> plusVar, minusVar;

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

      if (m_debugLevel >= 1)
         std::cout << " GPP:  par = " << rrv.GetName() << " cenVal = " << cenVal << " errSym = " << errHes
                   << " errAvgAsym = " << errAvg << std::endl;

      // Make Plus variation
      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal + errVal);
      plusVar.push_back(cloneFunc->getVal(nset));

      // Make Minus variation
      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal - errVal);
      minusVar.push_back(cloneFunc->getVal(nset));

      ((RooRealVar *)paramList.at(ivar))->setVal(cenVal);
   }

   TMatrixDSym         C(paramList.getSize());
   std::vector<double> errVec(paramList.getSize());
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

   if (m_debugLevel >= 1) std::cout << " GPP : sum = " << sqrt(sum) << std::endl;

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
