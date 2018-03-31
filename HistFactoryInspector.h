#ifndef __HistFactoryInspector_h__
#define __HistFactoryInspector_h__

/// \class HistFactoryInspector
/// A class to inspect HistFactory workspaces, retrieving yields and nuisance parameter impacts

#include <vector>
#include <map>
#include <TString.h>

class TFile;
class RooFormulaVar;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}
class RooSimultaneous;
class RooCategory;

namespace EXOSTATS {
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> YieldTable;
typedef std::map<TString, std::map<TString, std::pair<Double_t, Double_t>>> ImpactTable;

class HistFactoryInspector {
public:
   void setInput(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
                 TString rangeName);
   void setEvalRegions(std::vector<TString> regions);
   void setEvalRegions(TString regions);
   void setFitRegions(std::vector<TString> regions);
   void setFitRegions(TString regions);
   YieldTable  getYields();
   ImpactTable getImpacts(std::vector<TString> samples);
   ImpactTable getImpacts(TString samples);

protected:
   void           retrieveSampleNames();
   void           retrieveSampleNames(TString region);
   RooFormulaVar *retrieveYieldRFV(std::vector<TString> regions, std::vector<TString> components);
   RooFormulaVar *retrieveYieldRFV(TString region, std::vector<TString> components);

private:
   TString                                 m_inputFile;
   TString                                 m_workspaceName;
   TString                                 m_modelConfigName;
   TString                                 m_dataName;
   TString                                 m_rangeName;
   std::vector<TString>                    m_evalRegions;
   std::vector<TString>                    m_fitRegions;
   std::map<TString, std::vector<TString>> m_samples;

   TFile *                m_file;
   RooWorkspace *         m_w;
   RooStats::ModelConfig *m_mc;
   RooSimultaneous *      m_simPdf;
   const RooCategory *    m_cat;

   std::map<TString, RooFormulaVar *> m_yieldRFV;
};

} // namespace EXOSTATS

#endif
