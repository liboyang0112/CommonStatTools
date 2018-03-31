#ifndef __HistFactoryInspector_h__
#define __HistFactoryInspector_h__

/// \class HistFactoryInspector
/// A class to inspect HistFactory workspaces, retrieving yields and nuisance parameter impacts

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
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> YieldTable;
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> ImpactTable;

class HistFactoryInspector {
public:
   HistFactoryInspector();
   void setDebugLevel(Int_t level);
   void setInput(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
                 TString rangeName);
   void setEvalRegions(std::vector<TString> regions);
   void setEvalRegions(TString regions);
   void setFitRegions(std::vector<TString> regions);
   void setFitRegions(TString regions);
   YieldTable  getYields(Bool_t asymErrors = kTRUE);
   ImpactTable getImpacts(std::vector<TString> samples);
   ImpactTable getImpacts(TString samples);
   RooArgList  getFloatParList(const RooAbsPdf &pdf, const RooArgSet &obsSet);

protected:
   void           retrieveSampleNames();
   void           retrieveSampleNames(TString region);
   void           retrieveRooProductNames(TString region);
   RooFormulaVar *retrieveYieldRFV(std::vector<TString> regions, std::vector<TString> components);
   RooFormulaVar *retrieveYieldRFV(TString region, std::vector<TString> components);
   RooFitResult * fitPdfInRegions(std::vector<TString> regions, Bool_t saveResult = kFALSE, Bool_t doMinos = kTRUE);
   Double_t       getPropagatedError(RooAbsReal *var, const RooFitResult &fitResult, const Bool_t doAsym);

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

   std::map<TString, RooFormulaVar *> m_yieldRFV;
};

} // namespace EXOSTATS

#endif
