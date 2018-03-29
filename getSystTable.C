// A simple macro to evaluate pre- and post- fit impacts of nuisance parameters,
// using HistFactory workspaces
// Valerio Ippolito, Emma Tolley - Harvard University
// run with root -b -q getSystTable.C+

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

// helper class to hold param info
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

// string & output formatting

void printTableBegin(FILE *outfile)
{
   fputs("\\documentclass{article}\n\\usepackage[margin=0.5in]{geometry}\n", outfile);
   fputs("\\begin{document}\n\n\\begin{table}\n\\begin{center}\n\\setlength{\\tabcolsep}{0.0pc}\n", outfile);
   fputs("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lccccc}\n ", outfile);
   fputs("\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n", outfile);
   fputs("{\\bf Systematic} & SR & CR1mu & CR2mu & CR1e & CR2e \\\\ \n", outfile);
   fputs("\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n\n", outfile);
}

void printTableEnd(FILE *outfile)
{
   fputs("\\\\ \\hline\\noalign{\\smallskip}\n\\end{tabular*}\n\\end{center}\n\\end{table}\n\\end{document}\n\n\n",
         outfile);
}

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

TString prd(const double x, const int decDigits)
{
   stringstream ss;
   ss << fixed;
   ss.precision(decDigits); // set # places after decimal
   ss << x;
   return ss.str();
}

/////////////

bool is_np_constrained(TString input)
{
   return input.Contains("alpha") || input.Contains("gamma");
}

std::vector<TString> getAllComponentNamesInRegions(TString region, RooAbsPdf *regionPdf)
{
   // the rough structure of an HistFactory workspace is
   //   RooSimultaneous::simPdf
   //     RooProdPdf::model_REGION1
   //       RooGaussian::lumiConstraint
   //       RooGaussian::alpha_SYST1Constraint (for overallsys)
   //       RooPoisson::gamma_stat_REGION1_bin_BIN1_constraint (for stat uncertainty)
   //       RooRealSumPdf::REGION1_model (sum of RooProduct objects, each of them weighted by [the same] bin width)
   //         RooProduct::L_x_SAMPLE1_REGION1_overallSyst_x_Exp (term for sample 1, including uncertainties [RooHistFunc
   //         and FlexibleInterpVar or PieceWiseInterpolation])
   // note that NormFactors are included in the L_x_blabla term (e.g. this term changes integral when SigXsecOverSM is
   // changed)
   //
   // here we simply want to retrieve the list of RooProducts (L_x_blabla) associated to a given region (i.e. to the
   // RooRealSumPdf representing that region)
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

RooFormulaVar *getComponent(RooWorkspace *w, std::vector<TString> components, TString region, TString rangeName = "")
{
   // for a given list of samples in a region ('components'), we create a RooRealSumPdf representing
   // the sum of the corresponding components, and return its integral in the form of a RooFormulaVar
   // so that getVal() will tell us the number of events associated to this sum
   // a range for the integral can also be specified
   if (components.size() < 1) {
      throw std::runtime_error("component list is empty");
   }

   RooCategory *cat = w->cat("channelCat");

   RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(w->pdf("simPdf"));
   RooAbsPdf *      regPdf = simPdf->getPdf(region);

   RooAbsData *simData = w->data("obsData");
   RooAbsData *regData = simData->reduce(TString("channelCat==channelCat::" + region));

   RooRealVar *obs      = dynamic_cast<RooRealVar *>(regPdf->getObservables(regData)->find("obs_x_" + region));
   RooRealVar *binWidth = dynamic_cast<RooRealVar *>(regPdf->getVariables()->find("binWidth_obs_x_" + region + "_0"));

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
      return 0;
      throw std::runtime_error("something went wrong when fetching components");
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

RooFormulaVar *getComponent(RooWorkspace *w, std::vector<TString> components, std::vector<TString> regions,
                            TString rangeName = "")
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
      list.add(*getComponent(w, components, region, rangeName));
   }
   return new RooFormulaVar("sum_var", form, list);
   // return sum_var;

   /*RooFormulaVar *sum_var = 0;
   for (TString region: regions){
       RooFormulaVar* region_var = getComponent(w, components, region, rangeName);
       if (sum_var) sum_var = new RooFormulaVar("form_fracError", "@0 + @1", RooArgList(*region_var, *sum_var));
       else sum_var = region_var;
   }
   return sum_var;*/
}

std::vector<TString> getFreeParameters(RooWorkspace *w)
{
   // gets list of free parameters (i.e. non-constant parameters) given a pdf and a dataset (which
   // is used to determine which ones of the pdf free parameters are actually observables)
   // TODO: use RooStats::RemoveConstantParameters instead of the manual thing...
   RooArgSet *NPs = w->pdf("simPdf")->getParameters(*w->data("obsData"));

   TIterator *itr = NPs->createIterator();

   TObject *np = itr->Next();

   std::vector<TString> result;

   while (np) {
      RooAbsArg *raa = dynamic_cast<RooAbsArg *>(np);
      std::cout << "VALERIO: NP is " << raa->GetName() << " and constant=" << raa->isConstant() << std::endl;
      if (raa->isConstant() == kFALSE) result.push_back(np->GetName());
      np = itr->Next();
   }

   return result;
}

std::pair<Double_t, Double_t> getParameterImpactNoFit(TString param, RooWorkspace *w, RooFormulaVar *impact,
                                                      Bool_t useErrorVar = kFALSE)
{
   // uses a given RooFormulaVar (output of getComponent) to evaluate the impact, without fit, of
   // varying a nuisance parameter UP or DOWN by one sigma
   // variation is done either manually by +/- 1 (useErrorVar=false) or by avg +/- 1 sigma (useErrorVar=true)
   // NOTE: for OverallSys uncertainties, the output of this function must be ~identical to what's in the HistFactory's
   // XMLs
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

RooFitResult *fitPdfInRegions(RooWorkspace *w, std::vector<TString> regions, bool saveResult = false,
                              Bool_t doMinos = kTRUE)
{
   // fits simPdf only in a subset of regions
   // (strategy: define a new simultaneous PDF and a new dataset, including only fit regions)
   // TODO: use the better fitting technique used by HistFitter's Util::FitPdf
   // (in https://svnweb.cern.ch/trac/atlasphys/browser/Physics/SUSY/Analyses/HistFitter/trunk/src/Utils.cxx)

   // to do so, a temporary PDF and a temporary dataset have to be built
   RooSimultaneous *pdfFull  = dynamic_cast<RooSimultaneous *>(w->pdf("simPdf"));
   RooAbsData *     dataFull = w->data("obsData");
   RooCategory *    cat      = w->cat("channelCat");

   RooSimultaneous *pdf  = pdfFull;
   RooDataSet *     data = dynamic_cast<RooDataSet *>(dataFull);

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

void getSystTable(void)
{
   // a test function
   // to produce simple, easy-to-understand inputs:
   //  prepareHistFactor
   //  hist2workspace config/example.xml
   // TFile *f = TFile::Open("results/example_combined_GaussExample_model.root");

   ///////////////
   // SETUP     //
   ///////////////
   //  TString fullfitblind =
   //  "../HistFitter/results/SR700_Unblind_MJ135_independent/BkgOnly_L3p21ifb_combined_NormalMeasurement_model.root";
   TString fullfitblind = "/afs/cern.ch/user/e/etolley/work/public/MonoJet/workspaces/"
                          "testSimplifiedShapeFit_7Bins_Unblind_MJ135_DMA_750_10_gq0p25/"
                          "BkgOnly_L3p21ifb_combined_NormalMeasurement_model.root";
   TFile *       f = TFile::Open(fullfitblind);
   RooWorkspace *w = dynamic_cast<RooWorkspace *>(f->Get("combined"));

   // define binning, sample, and CR/SR configuration
   const std::vector<TString> low_bin_edges = {"250", "300", "350", "400", "500", "600", "700"};
   //    const std::vector<TString> low_bin_edges = {"700"};
   const std::vector<TString> samples = {"Znunu",   "Wmunu", "Wenu",    "Wtaunu", "Zmumu",   "Zee",
                                         "Ztautau", "top",   "diboson", "NCB",    "multijet"};
   const std::vector<TString> CRs     = {"CR1mu", "CR2mu", "CR1e"}; //, "CR2e"};
   const TString              SR      = "SR";
   // const std::vector<TString> CRs = { "CR1mu150", "CR2mu150", "CR1e150"};
   // const TString SR = "VR150";

   FILE *outfile = fopen("out_fit.tex", "w"); // file for tex output

   ///////////////////
   // END SETUP     //
   ///////////////////

   // automatically build fit regions and eval regions
   // fit region = all bins of all CRs
   // eval region = all bins of a given SR or CR
   std::vector<TString>                    fitRegions;
   std::map<TString, std::vector<TString>> evalRegions;
   for (TString bin_edge : low_bin_edges) {
      for (TString CR : CRs) {
         fitRegions.push_back(CR + "_" + bin_edge + "_cuts");
         evalRegions[CR].push_back(CR + "_" + bin_edge + "_cuts");
      }
      evalRegions[SR].push_back(SR + "_" + bin_edge + "_cuts");
   }

   // this marvellous object is what allows us to know the yield for a given sample (or combination of samples) in a
   // region taking into account all scale factors, systematics...
   map<TString, RooFormulaVar *> impactRFV;
   impactRFV[SR] = getComponent(w, samples, evalRegions[SR]);
   for (TString CR : CRs) impactRFV[CR] = getComponent(w, samples, evalRegions[CR]);

   // list of NPs (actually it can include more than the usual NPs)
   std::vector<TString> NPs = getFreeParameters(w); // VALERIO: contiene anche i mu

   // save state-of-the-art pre-fit parameter settings
   const TString snapBeforeFit = "pippo"; // sigh
   w->saveSnapshot(snapBeforeFit, *w->pdf("simPdf")->getParameters(*w->data("obsData")));

   std::stringstream myCout; // let's print everything at the end of the job

   ////////////////////
   // PRE-FIT IMPACT //
   ////////////////////

   w->loadSnapshot(snapBeforeFit); // crucial to do this every time

   // fetch raw SR/CR yields from the likelihood
   map<TString, Double_t> original_impact;
   original_impact[SR] = impactRFV[SR]->getVal();
   for (TString CR : CRs) original_impact[CR] = impactRFV[CR]->getVal();

   myCout << "\n\n\nCONSIDERING ONLY SAMPLES " << join(" ", samples) << "\n...IN BINS " << join(" ", low_bin_edges)
          << std::endl;
   myCout << "IMPACTS are given in %, in the form (DOWN, UP) where DOWN means NP goes DOWN by 1 SIGMA, etc"
          << std::endl;
   myCout << "\nNOMINAL " << SR << " EVENT YIELD: " << original_impact[SR] << std::endl;
   for (TString CR : CRs) myCout << "NOMINAL " << CR << " EVENT YIELD: " << original_impact[CR] << std::endl;
   myCout << "\nPRE-FIT IMPACTS" << std::endl;
   myCout << left << setw(46) << setfill(' ') << "SYSTEMATIC";
   myCout << left << setw(20) << setfill(' ') << SR;
   for (TString CR : CRs) myCout << left << setw(20) << setfill(' ') << CR;
   myCout << std::endl;

   printTableBegin(outfile);
   for (auto np : NPs) {
      if (!is_np_constrained(np)) {
         std::cout << "WARNING: NP to be ignored in pre-fit printouts found, named " << np << std::endl;
         continue;
      }
      w->loadSnapshot(snapBeforeFit);

      // fetch impacted SR/CR yields from the likelihood
      map<TString, std::pair<Double_t, Double_t>> nofit_impact;
      nofit_impact[SR] = getParameterImpactNoFit(np, w, impactRFV[SR], kFALSE);
      for (TString CR : CRs) nofit_impact[CR] = getParameterImpactNoFit(np, w, impactRFV[CR], kFALSE);

      float down = 100. * (nofit_impact[SR].first / original_impact[SR] - 1);
      float up   = 100. * (nofit_impact[SR].second / original_impact[SR] - 1);
      myCout << left << setw(46) << setfill(' ') << np;
      myCout << left << setw(20) << setfill(' ') << "(" + prd(up, 2) + "%, " + prd(down, 2) + "%)";
      fprintf(outfile, "%s & ($%2.2f$, $%2.2f$)", (np.ReplaceAll("_", "\\_")).Data(), up, down);

      for (TString CR : CRs) {
         down = 100. * (nofit_impact[CR].first / original_impact[CR] - 1);
         up   = 100. * (nofit_impact[CR].second / original_impact[CR] - 1);
         myCout << left << setw(20) << setfill(' ') << "(" + prd(up, 2) + "%, " + prd(down, 2) + "%)";
         fprintf(outfile, " & ($%2.2f$, $%2.2f$)", up, down);
      }
      fputs(" \\\\ \\noalign{\\smallskip} \n", outfile);
      myCout << std::endl;
   }
   myCout << std::endl;
   printTableEnd(outfile);

   /////////////////
   // FIT         //
   /////////////////

   // perform pdf fit in specified regions
   w->loadSnapshot(snapBeforeFit);
   RooFitResult *fitResult = fitPdfInRegions(w, fitRegions, true);
   fitResult->Print();
   const Double_t full_fit_impact = impactRFV[SR]->getVal();
   const Double_t full_fit_err    = impactRFV[SR]->getPropagatedError(*fitResult);

   myCout << "FIT VALUES" << std::endl;
   // save impact and param values
   map<TString, ParamFitResult> bestfit;
   for (auto np : NPs) {
      const Double_t param_bestfit_val  = w->var(np)->getVal();
      const Double_t param_bestfit_up   = param_bestfit_val + w->var(np)->getErrorHi(); // NOTE: assumes one ran Minos
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

   myCout << "POST-FIT IMPACTS when fitting only in " << join(" ", fitRegions) << std::endl;

   printTableBegin(outfile);

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
      fitPdfInRegions(w, fitRegions, kTRUE, kFALSE); // VALERIO's speedup
      const Double_t plus1_fit_impact = impactRFV[SR]->getVal();

      // now, fix this parameter to -1 sigma (TODO: check also -1 absolute?)
      w->loadSnapshot(snapBeforeFit);
      w->var(np)->setVal(bestfit[np].val_down);
      w->var(np)->setConstant(kTRUE);
      fitPdfInRegions(w, fitRegions, kTRUE, kFALSE); // VALERIO's speedup
      const Double_t minus1_fit_impact = impactRFV[SR]->getVal();

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

         fprintf(outfile, "%s & ($%2.2f$, $%2.2f$) \\\\ \\noalign{\\smallskip} \n", (np.ReplaceAll("_", "\\_")).Data(),
                 up, down);

      } else {
         std::cout << "WARNING: free NP to be ignored in post-fit printouts found, named " << np << std::endl;
      }
   }

   printTableEnd(outfile);
   myCout << std::endl << "{";
   for (auto np : NPs) myCout << " " << np_error[np];
   myCout << "}" << std::endl;

   // fit without nuisance parameters
   w->loadSnapshot(snapBeforeFit);

   for (auto np : NPs)
      if (is_np_constrained(np)) w->var(np)->setConstant(kTRUE);
   RooFitResult * frozen_fitResult  = fitPdfInRegions(w, fitRegions, true);
   const Double_t frozen_fit_impact = impactRFV[SR]->getVal();
   const Double_t frozen_fit_err    = impactRFV[SR]->getPropagatedError(*frozen_fitResult);

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
            std::cout << "WARNING: ignoring covariance matrix term " << NPs[i] << " " << idx_i << " " << NPs[j] << " "
                      << idx_j << std::endl;
            continue;
         }
         //	       myCout << NPs[i] << " impt = " <<  np_error[NPs[i]] << ", "
         //		      << NPs[j] << " impt = " <<  np_error[NPs[j]] << ", COV = " <<  CoVarMatrix(idx_i,idx_j) <<  ", CORR =
         //" << fitResult->correlation(NPs[i], NPs[j]) << std::endl;
         total_error += np_error[NPs[i]] * np_error[NPs[j]] * CoVarMatrix(idx_i, idx_j);
      }
   }
   myCout << "Total SYST error: " << sqrt(total_error) * 100 << "%" << endl;

   // printout (avoids MINUIT verbosity)
   std::cout << myCout.str() << std::endl;
   fclose(outfile);
}
