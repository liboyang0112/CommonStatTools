#ifndef __FITCROSSCHECKFORLIMITS_H__
#define __FITCROSSCHECKFORLIMITS_H__

#include <TString.h>
#include <map>
#include <vector>
#include <list>

class RooFitResult;
namespace RooStats {
class ModelConfig;
}
class RooAbsPdf;
class RooAbsData;
class RooAbsReal;
class RooRealVar;
class TCanvas;
class TF1;
class TLatex;
class TDirectory;
class RooArgSet;
class TLegend;
class TH2D;
class TH1;
class TGraph;
class TH1F;
class TFile;
class TDirectory;

using namespace std;

enum Algs {
   PlotHistosBeforeFit = 0,
   PlotMorphingControlPlots,
   PlotHistosAfterFitEachSubChannel,
   PlotHistosAfterFitGlobal,
   PlotsNuisanceParametersVSmu,
   PlotsStatisticalTest,
   FitToAsimov
};

class LimitCrossChecker {
public:
   LimitCrossChecker();
   void resetParams();

   void setDebugLevel(Int_t value);
   void setDrawPlots(bool value);            ///< create eps & png files and creat a webpage
   void setPlotRelative(bool value);         ///< plot % shift of systematic
   void setDraw1DResponse(bool value);       ///< draw 1D response for each NP
   void setWritePostfitAsimData(bool value); ///< create and add to workspace Asimov data created from fit
   void setUseMinosError(bool value);        ///< compute minos error (if false : use minuit error)
   void setBlind(bool value);                ///< blind di-jet mass from 100-150 GeV
   void setMakePostFitPlots(bool value);     ///< takes a lot of time and the error band is wrong
   void setPullMaxAcceptable(double value);  ///< Threshold to consider a NP[central value] as suspicious
   void setErrorMinAcceptable(double value); ///< Threshold to consider a NP[error] as suspicious
   void setXAxisLabel(TString value);        ///< set what the x-axis of the distribution is
   void setNJobs(int value);                 ///< number of subjobs for parallel processing
   void setIJob(int value);                  ///< index of subjob: i = 0..n-1

   Int_t   getDebugLevel();
   bool    getDrawPlots();
   bool    getPlotRelative();
   bool    getDraw1DResponse();
   bool    getWritePostfitAsimData();
   bool    getUseMinosError();
   bool    getBlind();
   bool    getMakePostFitPlots();
   double  getPullMaxAcceptable();
   double  getErrorMinAcceptable();
   TString getXAxisLabel();
   int     getNJobs();
   int     getIJob();

public:
   // Global functions
   RooFitResult *FitPDF(RooStats::ModelConfig *model, RooAbsPdf *fitpdf, RooAbsData *fitdata,
                        TString minimType = "Minuit2");
   void          PlotHistosBeforeFit(double nSigmaToVary, double mu);
   void          PlotMorphingControlPlots();
   void          PlotHistosAfterFitEachSubChannel(bool IsConditionnal, double mu);
   void          PlotHistosAfterFitGlobal(bool IsConditionnal, double mu, bool isAsimov = false);
   void          PlotsNuisanceParametersVSmu();
   void          PlotsStatisticalTest(double mu_pe, double mu_hyp, int nToyMC = 100, int rndmSeed = 0);
   void Plot1DResponse(RooAbsReal *nll, RooRealVar *var, TString cname, TCanvas *can, TF1 *poly, bool IsFloating,
                       TLatex *latex, TDirectory *tdir, RooArgSet *SliceSet = 0);
   void Plot1DResponseNew(RooAbsReal *nll, RooRealVar *var, TString cname, TCanvas *can, bool IsFloating, TLatex *latex,
                          TDirectory *tdir, TString snapshotName);
   TTree *                         createObservableTree(RooStats::ModelConfig *model);
   void                            setObservableTreeValues(RooStats::ModelConfig *model, TTree *tree);
   double                          FindMuUpperLimit();
   void                            PrintModelObservables();
   void                            PrintNuisanceParameters();
   void                            PrintAllParametersAndValues(RooArgSet para);
   void                            PrintNumberOfEvents(RooAbsPdf *pdf);
   void                            FindConstants(RooAbsPdf *pdf);
   void                            PrintSubChannels();
   void                            PrintSuspiciousNPs();
   bool                            IsSimultaneousPdfOK();
   bool                            IsChannelNameOK();
   void                            SetAllNuisanceParaToSigma(double Nsigma);
   void                            SetAllStatErrorToSigma(double Nsigma);
   void                            SetNuisanceParaToSigma(RooRealVar *var, double Nsigma);
   void                            GetNominalValueNuisancePara();
   void                            SetNominalValueNuisancePara();
   void                            SetPOI(double mu);
   void                            SetStyle();
   void                            LegendStyle(TLegend *l);
   int                             GetPosition(RooRealVar *var, TH2D *corrMatrix);
   list<pair<RooRealVar *, float>> GetOrderedCorrelations(RooRealVar *var, RooFitResult *fitres);

   /// create the canvas and put stuff on it
   /// to be used when plotting the +/- 1 sigma shifts
   TCanvas *DrawShift(TString channel, TString var, TString comp, double mu, TH1 *d, TH1 *n, TH1 *p1s, TH1 *m1s);
   TH1F *   MakeHist(TString name, RooCurve *curve);
   TH1F *   ConvertGraphToHisto(TGraph *pGraph);
   TH2D *GetSubsetOfCorrMatrix(RooRealVar *var, list<pair<RooRealVar *, float>> &pairs, RooFitResult *fitres, int size);
   void  FillGraphIntoHisto(TGraph *pGraph, TH1F *pHisto);
   void  Initialize(const char *infile, const char *outputdir, const char *workspaceName, const char *modelConfigName,
                    const char *ObsDataName);
   void  Finalize(const char *infile);
   RooArgList getFloatParList(const RooAbsPdf &pdf, const RooArgSet &obsSet = RooArgSet());

   //======================================================
   // ================= Main function =====================
   //======================================================
   void PlotFitCrossChecks(const char *infile = "WorkspaceForTest1.root", const char *outputdir = "./results/",
                           const char *workspaceName = "combined", const char *modelConfigName = "ModelConfig",
                           const char *ObsDataName = "obsData");

   void run(const Algs algorithm = Algs::PlotHistosBeforeFit, float mu = 0, float sigma = 1,
                               bool IsConditional = true, const char *infile = "WorkspaceForTest1.root",
                               const char *outputdir = "./results/", const char *workspaceName = "combined",
                               const char *modelConfigName = "ModelConfig", const char *ObsDataName = "obsData",
                               bool draw1DResponse = false, bool createPostfitAsimov = false);

   struct NPContainer {
      TString NPname;
      double  NPvalue;
      double  NPerrorHi;
      double  NPerrorLo;
      TString WhichFit;
   };

private:
   // User configuration one
   Int_t   debugLevel;           ///< debug level
   bool    drawPlots;            ///< create eps & png files and creat a webpage
   bool    plotRelative;         ///< plot % shift of systematic
   bool    draw1DResponse;       ///< draw 1D response for each NP
   bool    writePostfitAsimData; ///< create and add to workspace Asimov data created from fit
   bool    UseMinosError;        ///< compute minos error (if false : use minuit error)
   bool    blind;                ///< blind di-jet mass from 100-150 GeV
   bool    makePostFitPlots;     ///< takes a lot of time and the error band is wrong
   double  PullMaxAcceptable;    ///< Threshold to consider a NP[central value] as suspicious
   double  ErrorMinAcceptable;   ///< Threshold to consider a NP[error] as suspicious
   TString xAxisLabel;           ///< set what the x-axis of the distribution is
   int     nJobs;                ///< number of subjobs for parallel processing
   int     iJob;                 ///< index of subjob: i = 0..n-1

   static bool comp_second_abs_decend(const pair<RooRealVar *, float> &i, const pair<RooRealVar *, float> &j)
   {
      return fabs(i.second) > fabs(j.second);
   }

   // not switches
   RooWorkspace *         w;
   RooStats::ModelConfig *mc;
   RooAbsData *           data;
   TFile *                outputfile;
   double                 LumiRelError;
   TDirectory *           MainDirSyst;
   TDirectory *           MainDirMorphing;
   TDirectory *           MainDirFitEachSubChannel;
   TDirectory *           MainDirFitGlobal;
   TDirectory *           MainDirModelInspector;
   TDirectory *           MainDirStatTest;
   TDirectory *           MainDirFitAsimov;
   map<string, double>    MapNuisanceParamNom;
   vector<NPContainer>    AllNPafterEachFit;
   TString                OutputDir;
};
#endif
