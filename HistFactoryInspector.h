#ifndef __HistFactoryInspector_h__
#define __HistFactoryInspector_h__

/// \class HistFactoryInspector
///
/// A class to inspect HistFactory workspaces, retrieving yields and nuisance parameter impacts
///
/// \author Valerio Ippolito

#include <vector>
#include <map>
#include <RooArgList.h>

class TFile;
class RooFormulaVar;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}
class RooSimultaneous;
class RooCategory;
class RooFitResult;

namespace EXOSTATS {
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> YieldTableElement;
typedef std::pair<YieldTableElement, YieldTableElement>                     YieldTable;
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> ImpactTableElement;
typedef std::pair<ImpactTableElement, ImpactTableElement>                   ImpactTable;

class HistFactoryInspector {
public:
   /// Constructor
   HistFactoryInspector();

   void setDebugLevel(Int_t level);

   /// Set the input
   void setInput(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
                 TString rangeName);

   /// Set the regions where yields/impacts are evaluated
   void setEvalRegions(std::vector<TString> regions);

   /// Set the regions where yields/impacts are evaluated
   void setEvalRegions(TString regions);

   /// Set the regions where fit is performed
   void setFitRegions(std::vector<TString> regions);

   /// Set the regions where fit is performed
   void setFitRegions(TString regions);

   /// Retrieves event yields before and after fit, for selected evalRegions, when the fit is performed in the selected
   /// fitRegions
   YieldTable getYields(Bool_t asymErrors = kTRUE);

   /// Calculates the impact of each nuisance parameter on the sum of the selected samples, when the fit is performed in
   /// the selected fitRegions
   ImpactTable getImpacts(std::vector<TString> samples);

   /// Calculates the impact of each nuisance parameter on the sum of the selected samples, when the fit is performed in
   /// the selected fitRegions
   ImpactTable getImpacts(TString samples);

protected:
   void           retrieveSampleNames();
   void           retrieveSampleNames(TString region);
   void           retrieveRooProductNames(TString region);
   RooFormulaVar *retrieveYieldRFV(std::vector<TString> regions, std::vector<TString> components);
   RooFormulaVar *retrieveYieldRFV(TString region, std::vector<TString> components);
   RooFitResult * fitPdfInRegions(std::vector<TString> regions, Bool_t saveResult = kFALSE, Bool_t doMinos = kTRUE);
   std::pair<Double_t, Double_t> getYieldUpDown(TString param, RooFormulaVar *yield, Bool_t useErrorVar, Bool_t doFit,
                                                Bool_t doMinos);
   std::vector<TString>          getFreeParameters();
   Double_t   getPropagatedError(RooAbsReal *var, const RooFitResult &fitResult, const Bool_t doAsym);
   RooArgList getFloatParList(const RooAbsPdf &pdf, const RooArgSet &obsSet);
   void       resetError(const RooArgList &parList, const RooArgList &vetoList = RooArgList());

private:
   Int_t                m_debugLevel;
   TString              m_inputFile;
   TString              m_workspaceName;
   TString              m_modelConfigName;
   TString              m_dataName;
   TString              m_rangeName;
   std::vector<TString> m_evalRegions;
   std::vector<TString> m_fitRegions;
   std::map<TString, std::vector<TString>>
                                           m_products; ///< name of the RooProduct representing each sample in each region
   std::map<TString, std::vector<TString>> m_samples; ///< name of the regions, inferred from \c m_products

   TFile *                m_file;
   RooWorkspace *         m_w;
   RooStats::ModelConfig *m_mc;
   RooSimultaneous *      m_simPdf;
   RooCategory *          m_cat;
   TString                m_prefitSnap;
   TString                m_postfitSnap;

   std::map<TString, RooFormulaVar *> m_yieldRFV;
};

} // namespace EXOSTATS

#endif
