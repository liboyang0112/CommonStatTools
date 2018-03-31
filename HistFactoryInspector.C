#include <TFile.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooRealSumPdf.h>
#include <RooRealVar.h>
#include <RooProduct.h>

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

void EXOSTATS::HistFactoryInspector::setInput(const char *inputFile, const char *workspaceName,
                                              const char *modelConfigName, const char *dataName, TString rangeName)
{
   m_inputFile       = inputFile;
   m_workspaceName   = workspaceName;
   m_modelConfigName = modelConfigName;
   m_dataName        = dataName;
   m_rangeName       = rangeName;

   m_file   = new TFile(m_inputFile);
   m_w      = dynamic_cast<RooWorkspace *>(m_file->Get(m_workspaceName));
   m_mc     = dynamic_cast<RooStats::ModelConfig *>(m_w->obj(m_modelConfigName));
   m_simPdf = dynamic_cast<RooSimultaneous *>(m_mc->GetPdf());
   m_cat    = dynamic_cast<const RooCategory *>(&(m_simPdf->indexCat()));
}

void EXOSTATS::HistFactoryInspector::setEvalRegions(std::vector<TString> regions)
{
   m_evalRegions = regions;
}

void EXOSTATS::HistFactoryInspector::setEvalRegions(TString regions)
{
   setEvalRegions(getTokens(regions, ","));
}

void EXOSTATS::HistFactoryInspector::setFitRegions(std::vector<TString> regions)
{
   m_fitRegions = regions;
}

void EXOSTATS::HistFactoryInspector::setFitRegions(TString regions)
{
   setFitRegions(getTokens(regions, ","));
}

EXOSTATS::YieldTable EXOSTATS::HistFactoryInspector::getYields() {}

EXOSTATS::ImpactTable EXOSTATS::HistFactoryInspector::getImpacts(std::vector<TString> samples) {}

EXOSTATS::ImpactTable EXOSTATS::HistFactoryInspector::getImpacts(TString samples)
{
   return getImpacts(getTokens(samples, ","));
}

void EXOSTATS::HistFactoryInspector::retrieveSampleNames()
{
   for (auto reg : m_evalRegions) retrieveSampleNames(reg);
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
void EXOSTATS::HistFactoryInspector::retrieveSampleNames(TString region)
{
   const TString  rrsPdfName = TString::Format("%s_model", region.Data()); // hardcoded in HistFactory
   RooRealSumPdf *rrsPdf = dynamic_cast<RooRealSumPdf *>(m_simPdf->getPdf(region)->getComponents()->find(rrsPdfName));

   std::vector<TString> result;

   TIterator *itr = rrsPdf->funcList().createIterator();

   TObject *comp = itr->Next();

   while (comp) {
      result.push_back(comp->GetName());
      comp = itr->Next();
   }

   m_samples[region] = result;
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

   std::vector<TString> available = m_samples[region];

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
