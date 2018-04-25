/// \file
/// Tools to visualize fit correlation matrices

#include <TFile.h>
#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <TIterator.h>
#include <RooCatType.h>
#include <RooStats/RooStatsUtils.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <RooCategory.h>
#include <RooFitResult.h>

#include "Minimization.h"

R__LOAD_LIBRARY(Minimization.C+)
//R__LOAD_LIBRARY(CommonStatTools/build/libExoStats.so)

/// Saves the fit correlation matrix to file
/// \param[in] inputFile name of the input file
/// \param[in] workspaceName name of the input workspace
/// \param[in] modelConfigName name of the input ModelConfig
/// \param[in] dataName name of the dataset to fit
/// \param[in] workspaceTag prefix for the output ROOT file
/// \param[in] outputFolder path under which the output ROOT file will be stored; it will be created if it does not exist
/// \param[in] outputFormat format of the output file (e.g. ".root")
/// \param[in] doReduced keep this to kTRUE
/// \param[in] debugLevel (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
void getCorrMatrix(const char *inputFile, const char *workspaceName, const char *modelConfigName, const char *dataName,
                   TString workspaceTag, TString outputFolder, TString outputFormat, Bool_t doReduced = kTRUE,
                   Int_t debugLevel = 2)
{

   using namespace std;
   using namespace RooFit;
   using namespace RooStats;

   gStyle->SetOptStat(0);

   TFile                  _file0(inputFile);
   RooWorkspace *         w        = (RooWorkspace *)_file0.Get(workspaceName);
   RooStats::ModelConfig *mc       = (RooStats::ModelConfig *)w->obj(modelConfigName);
   RooRealVar *           firstPOI = dynamic_cast<RooRealVar *>(mc->GetParametersOfInterest()->first());
   firstPOI->setVal(0.0);

   RooSimultaneous *simPdf = (RooSimultaneous *)(mc->GetPdf());
   RooAbsData *     data   = w->data(dataName);

   RooCategory *channelCat = (RooCategory *)(&simPdf->indexCat());
   TIterator *  iter       = channelCat->typeIterator();
   RooCatType * tt         = NULL;

   while ((tt = (RooCatType *)iter->Next())) {
      cout << "category : " << tt->GetName() << " " << endl;
      const TString      ttname = tt->GetName();
      RooAbsPdf * pdftmp = simPdf->getPdf(tt->GetName());
      RooArgSet * obstmp = pdftmp->getObservables(*mc->GetObservables());
      RooAbsData *datatmp =
         data->reduce(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()));

      RooArgSet *      constrainedParams = pdftmp->getParameters(*data);
      const RooArgSet *glbObs            = mc->GetGlobalObservables();

      RooRealVar *poi = (RooRealVar *)mc->GetParametersOfInterest()->first();
      cout << "Constant POI: " << poi->isConstant() << endl;
      cout << "Value of POI: " << poi->getVal() << endl;

      RooStats::RemoveConstantParameters(constrainedParams);
      RooFit::Constrain(*constrainedParams);

      RooAbsReal *nll =
         pdftmp->createNLL(*datatmp, Constrain(*constrainedParams), GlobalObservables(*glbObs), Offset(1));
      double nllval = nll->getVal();

      std::cout << "initial parameters" << std::endl;
      constrainedParams->Print("v");

      std::cout << "INITIAL NLL = " << nllval << std::endl;

      RooFitResult *r = nullptr;

      EXOSTATS::minimize(nll, 3, nullptr, "", "", debugLevel, kTRUE, &r);

      TCanvas *c = new TCanvas(ttname, ttname, 1300, 1200);
      c->SetRightMargin(0.11);
      c->SetLeftMargin(0.28);
      c->SetBottomMargin(0.25);
      c->SetTopMargin(0.03);
      // r->correlationHist()->Draw("colz");
      // c->Print( Form( "%s.eps", ttname.c_str() ) );

      if (!doReduced) continue;

      TH2D *    h2cm = (TH2D *)r->correlationHist();
      const int px   = h2cm->GetXaxis()->GetNbins();
      const int py   = h2cm->GetYaxis()->GetNbins();

      TH2D *h2ncm = new TH2D("CorrelationMatrixNew", "", px, 0., px, py, 0, py);

      double valz;
      string valx, valy;
      int    nx = 0, ny = 0;

      for (int i1 = 0; i1 < px; i1++) {
         valx = h2cm->GetXaxis()->GetBinLabel(i1 + 1);
         cout << "debug " << valx << endl;
         // Don't keep bin systematics
         if (valx.find("_bin_") < valx.length()) continue;
         // remove alpha
         unsigned long pos_s = valx.find("alpha_");
         if (pos_s < valx.length()) {
            valx.erase(pos_s, 6);
         }
         // remove auto
         unsigned long pos_a = valx.find("auto_");
         while (pos_a < valx.length()) {
            valx.erase(pos_a, 5);
            pos_a = valx.find("auto_", pos_a);
         }
         // remove file one
         unsigned long pos_1 = valx.find("_fileOne");
         while (pos_1 < valx.length()) {
            valx.erase(pos_1, 8);
            pos_1 = valx.find("_fileOne", pos_1);
         }
         unsigned long pos_2 = valx.find("_fileTwo");
         while (pos_2 < valx.length()) {
            valx.erase(pos_2, 8);
            pos_2 = valx.find("_fileTwo", pos_2);
         }

         nx++;
         ny = 0;
         //
         for (int i2 = 0; i2 < py; i2++) {
            valy = h2cm->GetYaxis()->GetBinLabel(i2 + 1);
            // Dont keep bin systematics
            if (valy.find("_bin_") < valy.length()) continue;
            // Remove alpha
            unsigned long pos_sy = valy.find("alpha_");
            if (pos_sy < valy.length()) {
               valy.erase(pos_sy, 6);
            }
            // remove auto
            unsigned long pos_ay = valy.find("auto_");
            while (pos_ay < valy.length()) {
               valy.erase(pos_ay, 5);
               pos_ay = valy.find("auto_", pos_ay);
            }
            // remove file one
            unsigned long pos_1y = valy.find("_fileOne");
            while (pos_1y < valy.length()) {
               valy.erase(pos_1y, 8);
               pos_1y = valy.find("_fileOne", pos_1y);
            }
            unsigned long pos_2y = valy.find("_fileTwo");
            while (pos_2y < valy.length()) {
               valy.erase(pos_2y, 8);
               pos_2y = valy.find("_fileTwo", pos_2y);
            }

            ny++;
            //
            valz = h2cm->GetBinContent(i1 + 1, i2 + 1);
            // cout<<"debug "<<i1<<" "<<i2<<" "<<nx<<" "<<ny<<" "<<valz<<endl;
            h2ncm->SetBinContent(nx, ny, valz);
            h2ncm->GetXaxis()->SetBinLabel(nx, valx.c_str());
            h2ncm->GetYaxis()->SetBinLabel(ny, valy.c_str());
         }
      }

      // c->SetBottomMargin( 0.14 );
      // new TCanvas();
      h2ncm->SetTitle("");
      h2ncm->GetXaxis()->SetRangeUser(0, nx);
      h2ncm->GetYaxis()->SetRangeUser(0, ny);
      h2ncm->GetXaxis()->SetLabelSize(0.017);
      h2ncm->GetYaxis()->SetLabelSize(0.02);
      h2ncm->GetZaxis()->SetLabelSize(0.028);
      h2ncm->GetZaxis()->SetRangeUser(-1, 1);
      h2ncm->Draw("colz");
      c->Modified();
      c->Update();

      const TString fullOutFolder = outputFolder + "/corrmatrix/";
      system("mkdir -vp " + fullOutFolder);
      const TString outFileName = fullOutFolder + workspaceTag + "_corrmatrix_" + ttname + outputFormat;

      c->SaveAs(outFileName);
      // c->SaveAs( Form( "%s_%s_red.png",inname.c_str(), ttname.c_str() ) );
   }

   return;
}
