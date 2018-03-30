/// \file
/// A simple macro to evaluate pre- and post- fit impacts of nuisance parameters,
/// using HistFactory workspaces
/// \authors Valerio Ippolito, Emma Tolley - Harvard University

#include <RooWorkspace.h>
#include <RooRealSumPdf.h>
#include <TFile.h>
#include <RooFormulaVar.h>
#include <RooAbsData.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooFitResult.h>
#include <vector>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooProduct.h>
#include <RooStats/ModelConfig.h>
#include <RooDataSet.h>
#include <RooStats/RooStatsUtils.h>
#include <sstream>
#include <iomanip>

////////////////////////////////////
/// \class ParamFitResult
/// Helper class to hold parameter info
class ParamFitResult {
public:
   Double_t val;
   Double_t val_up;
   Double_t val_down;
   ParamFitResult() : val(0), val_up(0), val_down(0) {}
   ParamFitResult(const Double_t in_val, const Double_t in_val_up, const Double_t in_val_down)
      : val(in_val), val_up(in_val_up), val_down(in_val_down)
   {
   }
};

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

////////////////////////////////////
/// Given the nuisance parameter names, determines if it'c constrained or not
/// 
/// \param[in] input the nuisance parameter name
/// \param[out] result \c kTRUE if the parameter is constrained, \c kFALSE otherwise
///
/// This function could easily be implemented without relying on the name of
/// the parameter, but using RooFit methods only.
bool is_np_constrained(TString input)
{
   return input.Contains("alpha") || input.Contains("gamma");
}

////////////////////////////////////
/// Retrieves the names of all p.d.f. components (samples) in a given region
///
/// \param[in] region name of the region
/// \param[in] regionPdf pointer to the region p.d.f.
/// \param[out] result vector of names of likelihood components in that region
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
std::vector<TString> getAllComponentNamesInRegions(TString region, RooAbsPdf *regionPdf)
{
   const TString  rrsPdfName = TString::Format("%s_model", region.Data());
   RooRealSumPdf *rrsPdf     = dynamic_cast<RooRealSumPdf *>(regionPdf->getComponents()->find(rrsPdfName));

   std::vector<TString> result;

   TIterator *itr = rrsPdf->funcList().createIterator();

   TObject *comp = itr->Next();

   while (comp) {
      result.push_back(comp->GetName());
      comp = itr->Next();
   }

   return result;
}

/// Create a RooFormulaVar containing the number of events associated to a sum of RooProducts
///
/// \param[in] w pointer to the input RooWorspace
/// \param[in] components list of the components (samples) to be summed up
/// \param[in] region channel (region) to consider
/// \param[in] rangeName name of the variable range to be used for the integral (ignored if empty)
/// \param[in] dataName name of the dataset (used to determine observables)
/// \param[out] form_frac pointer to the output RooFormulaVar
///
/// For a given list of samples in a region ('components'), we create a \c RooRealSumPdf representing
/// the sum of the corresponding components, and return its integral in the form of a RooFormulaVar
/// so that getVal() will tell us the number of events associated to this sum.
/// A range for the integral can also be specified
RooFormulaVar *getComponent(RooWorkspace *w, std::vector<TString> components, TString region, TString rangeName = "", TString dataName = "obsData")
{
   if (components.size() < 1) {
      throw std::runtime_error("component list is empty");
   }

   RooCategory *cat = w->cat("channelCat"); // name hardcoded in HistFactory

   RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(w->pdf("simPdf")); // name hardcoded in HistFactory
   RooAbsPdf *      regPdf = simPdf->getPdf(region);

   RooAbsData *simData = w->data(dataName); // name hardcoded in HistFactory
   RooAbsData *regData = simData->reduce(TString("channelCat==channelCat::" + region)); // name hardcoded in HistFactory

   RooRealVar *obs      = dynamic_cast<RooRealVar *>(regPdf->getObservables(regData)->find("obs_x_" + region)); // name hardcoded in HistFactory
   RooRealVar *binWidth = dynamic_cast<RooRealVar *>(regPdf->getVariables()->find("binWidth_obs_x_" + region + "_0")); // name hardcoded in HistFactory

   // retrieve components

   RooArgList compFuncList;
   RooArgList compCoefList;

   std::vector<TString> available = getAllComponentNamesInRegions(region, regPdf);

   for (auto avail : available) {
      for (auto wanted : components) {
         const TString target = "_" + wanted + "_";

         if (avail.Contains(target)) {
            compFuncList.add(*(dynamic_cast<RooProduct *>(w->obj(avail))));
            compCoefList.add(*binWidth);
         }
      }
   }

   if (compFuncList.getSize() == 0 || compCoefList.getSize() == 0 || compFuncList.getSize() != compCoefList.getSize()) {
      throw std::runtime_error("something went wrong when fetching components");
      //return 0;
   }

   // create RRSPdf and integral
   const TString compList = join("_", components);

   const TString  compRRSname = "RRS_region_" + region + "_" + compList;
   RooRealSumPdf *compRRS     = new RooRealSumPdf(compRRSname, compRRSname, compFuncList, compCoefList);

   RooAbsReal *compFunc(nullptr);
   if (rangeName != "")
      compFunc = compRRS->createIntegral(RooArgSet(*obs), rangeName);
   else
      compFunc = compRRS->createIntegral(RooArgSet(*obs));

   // create RooFormulaVar

   TString compRFVname = "form_frac_region_" + region + "_" + compList;
   if (rangeName != "") compRFVname += "_" + rangeName;

   RooFormulaVar *form_frac = new RooFormulaVar("form_fracError", "@0", RooArgList(*compFunc));
   form_frac->SetNameTitle(compRFVname, compRFVname);

   return form_frac;
}

/// Sum up components in different regions
///
/// \param[in] w pointer to the input RooWorspace
/// \param[in] components list of the components (samples) to be summed up
/// \param[in] regions channels (regions) to consider
/// \param[in] rangeName name of the variable range to be used for the integral (ignored if empty)
/// \param[in] dataName name of the dataset (used to determine observables)
/// \param[out] form_frac pointer to the output RooFormulaVar
///
/// See the single-region version of getComponent() for more details
RooFormulaVar *getComponent(RooWorkspace *w, std::vector<TString> components, std::vector<TString> regions,
                            TString rangeName = "", TString dataName = "obsData")
{
   RooArgList list;
   TString    form  = "";
   int        count = 0;
   for (TString region : regions) {
      if (form.Length() == 0)
         form += "@";
      else
         form += "+ @";
      form += count;
      count++;
      list.add(*getComponent(w, components, region, rangeName, dataName));
   }
   return new RooFormulaVar("sum_var", form, list);
}

/// Get list of free parameters of a p.d.f.
///
/// \param[in] w pointer to the input RooWorspace
/// \param[in] pdfName name of the p.d.f.
/// \param[in] dataName name of the dataset (used to determine observables)
/// \param[out] result vector of nuisance parameter names
/// 
/// Gets list of free parameters (i.e. non-constant parameters) given a pdf and a dataset (which
/// is used to determine which ones of the pdf free parameters are actually observables).
std::vector<TString> getFreeParameters(RooWorkspace *w, TString pdfName = "simPdf", TString dataName = "obsData")
{
   // TODO: use RooStats::RemoveConstantParameters instead of the manual thing...
   RooArgSet *NPs = w->pdf(pdfName)->getParameters(*w->data(dataName));

   TIterator *itr = NPs->createIterator();

   TObject *np = itr->Next();

   std::vector<TString> result;

   while (np) {
      RooAbsArg *raa = dynamic_cast<RooAbsArg *>(np);
      std::cout << "NP named " << raa->GetName() << " has constant=" << raa->isConstant() << std::endl;
      if (raa->isConstant() == kFALSE) result.push_back(np->GetName());
      np = itr->Next();
   }

   return result;
}

/// Evaluate the impact of a nuisance parameter on a \c RooFormulaVar (without fit)
///
/// \param[in] param name of the nuisance parameter (NP)
/// \param[in] w pointer to the input RooWorspace
/// \param[in] impact pointer to the RooFormulaVar used to calculate the impact
/// \param[in] useErrorVal if \c kTRUE, will use fitted errors as plus and minus NP variations; otherwise, will use +/- 1
/// \param[out] result pair of impact of up and down variation on the provided RooFormulaVar
/// 
/// Uses a given RooFormulaVar (typically the output of getComponent) to evaluate the impact, without fit, of
/// varying a nuisance parameter up or down by one sigma. The variation is done either manually
/// by \c avg +/- 1 (\c useErrorVar set to \c false) or by \c avg +/- 1 sigma (\c useErrorVar set to \c true)
///
/// Note that, for OverallSys uncertainties, the output of this function must be ~identical to what specified in the HistFactory's
/// XMLs.
std::pair<Double_t, Double_t> getParameterImpactNoFit(TString param, RooWorkspace *w, RooFormulaVar *impact,
                                                      Bool_t useErrorVar = kFALSE)
{
   std::pair<Double_t, Double_t> result;

   // TODO: if useErrorVar = kTRUE, this makes sense only for alpha_ nuisance parameters
   const Double_t var_value      = (useErrorVar) ? (w->var(param)->getVal()) : 0;
   const Double_t var_value_up   = (useErrorVar) ? (var_value + w->var(param)->getErrorHi()) : 1;
   const Double_t var_value_down = (useErrorVar) ? (var_value + w->var(param)->getErrorLo()) : -1;

   // up
   w->var(param)->setVal(var_value_up);
   result.first = impact->getVal();

   // down
   w->var(param)->setVal(var_value_down);
   result.second = impact->getVal();

   return result;
}

/// Fits the pdf of the HistFactory workspace only in a subset of the available channels (regions)
///
/// \param[in] w pointer to the input RooWorkspace
/// \param[in] dataName name of the dataset to fit
/// \param[in] regions vector of names of channels (regions) to be considered in the fit
/// \param[in] saveResult if set to \c kTRUE, the RooFitResult is returned
/// \param[in] doMinos if set to \c kTRUE, Minos is run (much slower!)
///
/// To do the job, this function defines a new simultaneous PDF and a new dataset, including only fit regions.
RooFitResult *fitPdfInRegions(RooWorkspace *w, TString dataName, std::vector<TString> regions, bool saveResult = false,
                              Bool_t doMinos = kTRUE)
{
   // TODO: use the better fitting technique used by HistFitter's Util::FitPdf
   // (in https://svnweb.cern.ch/trac/atlasphys/browser/Physics/SUSY/Analyses/HistFitter/trunk/src/Utils.cxx)

   // NOTE: all pdf/object names are specified by HistFactory by default, that's why they are hardcoded

   // to do so, a temporary PDF and a temporary dataset have to be built
   RooSimultaneous *pdfFull  = dynamic_cast<RooSimultaneous *>(w->pdf("simPdf"));
   RooAbsData *     dataFull = w->data(dataName);
   RooCategory *    cat      = w->cat("channelCat");

   RooSimultaneous *pdf  = pdfFull;
   RooDataSet *     data = dynamic_cast<RooDataSet *>(dataFull);


   // determine useful terms
   std::vector<RooAbsPdf *>  pdfVec;
   std::vector<RooDataSet *> dataVec;

   for (auto region : regions) {
      if (cat->setLabel(region, kTRUE)) {
         throw std::runtime_error("Unknown region found");
      } else {
         RooDataSet *regionData =
            dynamic_cast<RooDataSet *>(data->reduce(TString("channelCat==channelCat::" + region)));
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

   pdf = new RooSimultaneous("simPdfReduced_" + nickname, "simultaneous pdf reduced to regions " + nickname, *cat);
   for (UInt_t i = 0; i < pdfVec.size(); i++) pdf->addPdf(*pdfVec[i], regions[i]);

   // perform the fit (this is the part which should possibly be improved)

   RooArgSet *allParams = pdf->getParameters(data);
   RooStats::RemoveConstantParameters(allParams);
   RooStats::ModelConfig *mc        = dynamic_cast<RooStats::ModelConfig *>(w->obj("ModelConfig"));
   const RooArgSet *      globObs   = mc->GetGlobalObservables();
   RooFitResult *         fitResult = 0;
   if (saveResult)
      fitResult = pdf->fitTo(*data, RooFit::GlobalObservables(*globObs), RooFit::Minos(doMinos), RooFit::Save());
   else
      pdf->fitTo(*data, RooFit::GlobalObservables(*globObs), RooFit::Minos(doMinos));
   return fitResult;
}

/// Obtain a table containing the impact, on the normalisation of a few samples in given regions, of nuisance parameters.
/// 
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input \c ModelConfig
/// \param[in] dataName name of the dataset
/// \param[in] workspaceTag prefix for the output ROOT file
/// \param[in] outputFolder path under which the output ROOT file will be stored; it will be created if it does not exist
/// \param[in] samples comma-separated list of samples ("components") to consider
/// \param[in] evaluationRegions comma-separated list of regions ("channels") for which impacts must be calculated
/// \param[in] fitRegions comma-separated list of regions to be used in the fit
///
/// Let us suppose you want to evaluate the impact of a nuisance parameter (NP) on the signal region yield of the sum of
/// all backgrounds, when the fit is performed using only control region data. You would then have
/// \code
/// samples = "bkg1,bkg2,bkg3"
/// evaluationRegions = "SR"
/// fitRegions = "CR1,CR2,CR3"
/// \endcode
/// and this function will print the impact of all nuisance parameters on the sum of the yields of bkg1, bkg2 and bkg3
/// in the region "SR".
///
/// This impact is defined for each NP as the change in the bkg1+bkg2+bkg3 yield evaluated in the SR, in the following way:
///   -# a fit to CR1+CR2+CR3 is performed, and the yield is saved in memory (let us denote it by \c std)
///   -# the NP is set to <bestfit>+<onesigma>, where <bestfit> and <onesigma> are determined by the fit above
///   -# the NP is set to constant and the fit is performed again
///   -# the new yield is saved in memory (let us denote it by \c up)
///   -# the "up" impact of this NP is then evaluated by the ratio between \c up and \c std
///   -# the procedure is repeated for the "down" variation
///
/// Output is also written in the folder outputFolder/np_impact/workspaceTag_systtable.tex
void getSystTable(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
                  TString workspaceTag, TString outputFolder, TString samples, TString evaluationRegions,
                  TString fitRegions)
{
   TFile *                f  = TFile::Open(inputFile);
   RooWorkspace *         w  = dynamic_cast<RooWorkspace *>(f->Get(workspaceName));
   RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig *>(w->obj(modelConfigName));

   outputFolder += "/np_impact/";
   system("mkdir -p " + outputFolder);
   FILE *outfile = fopen(outputFolder + workspaceTag + "_systtable.tex", "w"); // file for tex output

   // transform strings in vectors via tokenization
   std::vector<TString> samplesVec     = getTokens(samples, ",");
   std::vector<TString> fitRegionsVec  = getTokens(fitRegions, ",");
   std::vector<TString> evalRegionsVec = getTokens(evaluationRegions, ",");

   // this marvellous object is what allows us to know the yield for a given sample (or combination of samples) in a
   // region taking into account all scale factors, systematics...
   map<TString, RooFormulaVar *> impactRFV;
   for (auto reg : evalRegionsVec) impactRFV[reg] = getComponent(w, samplesVec, reg, dataName);

   // list of NPs (actually it can include more than the usual NPs: it contains the signal strength!)
   std::vector<TString> NPs = getFreeParameters(w);

   // save state-of-the-art pre-fit parameter settings
   const TString snapBeforeFit = "adummyname"; // sigh
   w->saveSnapshot(snapBeforeFit, *mc->GetPdf()->getParameters(*w->data(dataName)));

   std::stringstream myCout; // let's print everything at the end of the job

   ////////////////////
   // PRE-FIT IMPACT //
   ////////////////////

   w->loadSnapshot(snapBeforeFit); // crucial to do this every time

   // fetch raw yields in each region from the likelihood
   map<TString, Double_t> prefit_yield;
   for (auto reg : evalRegionsVec) prefit_yield[reg] = impactRFV[reg]->getVal();

   myCout << "\n\n\nCONSIDERING ONLY SAMPLES " << join(" ", samplesVec) << std::endl;
   myCout << "IMPACTS are given in %, in the form (DOWN, UP) where DOWN means NP goes DOWN by 1 SIGMA, etc"
          << std::endl;
   for (auto reg : evalRegionsVec) myCout << "NOMINAL " << reg << " EVENT YIELD: " << prefit_yield[reg] << std::endl;
   myCout << "\nPRE-FIT IMPACTS" << std::endl;
   myCout << left << setw(46) << setfill(' ') << "SYSTEMATIC";
   myCout << left << setw(20) << setfill(' ');
   for (TString reg : evalRegionsVec) myCout << left << setw(20) << setfill(' ') << reg;
   myCout << std::endl;

   for (auto np : NPs) {
      if (!is_np_constrained(np)) {
         std::cout << "WARNING: NP to be ignored in pre-fit printouts found, named " << np << std::endl;
         continue;
      }
      w->loadSnapshot(snapBeforeFit);

      // fetch impacted yields from the likelihood
      map<TString, std::pair<Double_t, Double_t>> nofit_impact;
      for (auto reg : evalRegionsVec) {
         nofit_impact[reg] = getParameterImpactNoFit(np, w, impactRFV[reg], kFALSE); // kFALSE = use +/- 1 instead of +/- err

         const Float_t down = 100. * (nofit_impact[reg].first / prefit_yield[reg] - 1);
         const Float_t up   = 100. * (nofit_impact[reg].second / prefit_yield[reg] - 1);
         myCout << left << setw(46) << setfill(' ') << np;
         myCout << left << setw(20) << setfill(' ') << "(" + prd(up, 2) + "%, " + prd(down, 2) + "%)";
         fprintf(outfile, " & ($%2.2f$, $%2.2f$)", up, down);
      }
      myCout << std::endl;
   }
   myCout << std::endl;

   /////////////////
   // FIT         //
   /////////////////

   // perform pdf fit in specified regions
   w->loadSnapshot(snapBeforeFit);
   RooFitResult *fitResult = fitPdfInRegions(w, dataName, fitRegionsVec, true);
   fitResult->Print();
   for (auto reg : evalRegionsVec) {
      const Double_t full_fit_impact = impactRFV[reg]->getVal();
      const Double_t full_fit_err    = impactRFV[reg]->getPropagatedError(*fitResult);

      myCout << "FIT VALUES" << std::endl;
      // save impact and param values
      map<TString, ParamFitResult> bestfit;
      for (auto np : NPs) {
         const Double_t param_bestfit_val = w->var(np)->getVal();
         const Double_t param_bestfit_up  = param_bestfit_val + w->var(np)->getErrorHi(); // NOTE: assumes one ran Minos
         const Double_t param_bestfit_down = param_bestfit_val + w->var(np)->getErrorLo();
         ParamFitResult param_bestfit(param_bestfit_val, param_bestfit_up, param_bestfit_down);
         bestfit[np] = param_bestfit;
         if (np.Contains("alpha"))
            myCout << left << setw(46) << setfill(' ') << np;
         else if (np.Contains("gamma"))
            myCout << left << setw(33) << setfill(' ') << np;
         else
            myCout << left << setw(7) << setfill(' ') << np;
         myCout << left << setw(10) << setfill(' ') << " = " + prd(param_bestfit_val, 3);
         myCout << left << setw(10) << setfill(' ') << " + " + prd(w->var(np)->getErrorHi(), 3);
         myCout << left << setw(10) << setfill(' ') << " - " + prd(w->var(np)->getErrorLo() * -1, 3);
         myCout << endl;
      }

      /////////////////////
      // POST-FIT IMPACT //
      /////////////////////

      myCout << "POST-FIT IMPACTS when fitting only in " << join(" ", fitRegionsVec) << std::endl;

      map<TString, float> np_error;
      for (auto np : NPs) {
         // make sure that we don't evaluate the post-fit impact of a parameter which
         // was not part of the original fit (e.g. SR-only NP, and CR-only fit)
         if (fitResult->floatParsFinal().index(np) < 0) {
            std::cout << "WARNING: not-fitted NP to be ignored in post-fit printouts found, named " << np << std::endl;
            continue;
         }

         // now, fix this parameter to +1 sigma (TODO: check also +1 absolute?)
         w->loadSnapshot(snapBeforeFit);
         w->var(np)->setVal(bestfit[np].val_up);
         w->var(np)->setConstant(kTRUE);
         fitPdfInRegions(w, dataName, fitRegionsVec, kTRUE, kFALSE); // VALERIO's speedup
         const Double_t plus1_fit_impact = impactRFV[reg]->getVal();

         // now, fix this parameter to -1 sigma (TODO: check also -1 absolute?)
         w->loadSnapshot(snapBeforeFit);
         w->var(np)->setVal(bestfit[np].val_down);
         w->var(np)->setConstant(kTRUE);
         fitPdfInRegions(w, dataName, fitRegionsVec, kTRUE, kFALSE); // VALERIO's speedup
         const Double_t minus1_fit_impact = impactRFV[reg]->getVal();

         // finally, evaluate impact!
         const float down = 100. * (minus1_fit_impact / full_fit_impact - 1);
         const float up   = 100. * (plus1_fit_impact / full_fit_impact - 1);

         // average up and down for total uncertainty calculation
         np_error[np] = 0.5 * ((plus1_fit_impact / full_fit_impact - 1) - (minus1_fit_impact / full_fit_impact - 1));

         // but print things only if this is not a free fit parameter!
         if (is_np_constrained(np)) {
            myCout << left << setw(46) << setfill(' ') << np;
            myCout << left << setw(20) << setfill(' ') << "(" + prd(up, 2) + "%, " + prd(down, 2) + "%)";
            myCout << std::endl;

            fprintf(outfile, "%s & ($%2.2f$, $%2.2f$) \\\\ \\noalign{\\smallskip} \n",
                    (np.ReplaceAll("_", "\\_")).Data(), up, down);

         } else {
            std::cout << "WARNING: free NP to be ignored in post-fit printouts found, named " << np << std::endl;
         }
      }

      myCout << std::endl << "{";
      for (auto np : NPs) myCout << " " << np_error[np];
      myCout << "}" << std::endl;

      // fit without nuisance parameters
      w->loadSnapshot(snapBeforeFit);

      for (auto np : NPs)
         if (is_np_constrained(np)) w->var(np)->setConstant(kTRUE);
      RooFitResult * frozen_fitResult  = fitPdfInRegions(w, dataName, fitRegionsVec, true);
      const Double_t frozen_fit_impact = impactRFV[reg]->getVal();
      const Double_t frozen_fit_err    = impactRFV[reg]->getPropagatedError(*frozen_fitResult);

      myCout << "FITTED EVENT YIELD: " << full_fit_impact << std::endl;
      myCout << "NO-NP   EVENT YIELD: " << frozen_fit_impact << std::endl;
      myCout << "Total [?] error: " << 100. * full_fit_err / full_fit_impact << "%" << std::endl;
      myCout << "Stat [?]  error: " << 100. * frozen_fit_err / full_fit_impact << "%" << std::endl;
      // myCout << "Syst  error: " << 100. * (full_fit_err - frozen_fit_err)/full_fit_impact << "%" << std::endl;

      TMatrixDSym CoVarMatrix(fitResult->covarianceMatrix());
      Double_t    total_error = 0;
      for (unsigned int i = 0; i < NPs.size(); i++) {
         for (unsigned int j = i; j < NPs.size(); j++) {
            const int idx_i = fitResult->floatParsFinal().index(NPs[i]);
            const int idx_j = fitResult->floatParsFinal().index(NPs[j]);
            if (idx_i < 0 || idx_j < 0) {
               std::cout << "WARNING: ignoring covariance matrix term " << NPs[i] << " " << idx_i << " " << NPs[j]
                         << " " << idx_j << std::endl;
               continue;
            }
            //	       myCout << NPs[i] << " impt = " <<  np_error[NPs[i]] << ", "
            //		      << NPs[j] << " impt = " <<  np_error[NPs[j]] << ", COV = " <<  CoVarMatrix(idx_i,idx_j) <<  ",
            // CORR
            //" << fitResult->correlation(NPs[i], NPs[j]) << std::endl;
            total_error += np_error[NPs[i]] * np_error[NPs[j]] * CoVarMatrix(idx_i, idx_j);
         }
      }
      myCout << "Total SYST error: " << sqrt(total_error) * 100 << "%" << endl;
      
      // we do this in preparation for the next region
      w->loadSnapshot(snapBeforeFit);
   }

   // printout (avoids MINUIT verbosity)
   std::cout << myCout.str() << std::endl;
   fclose(outfile);
}
