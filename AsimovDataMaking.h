#ifndef __EXOSTATS_ASIMOVDATAMAKING_H__
#define __EXOSTATS_ASIMOVDATAMAKING_H__

/// \file
/// Tools to create Asimov datasets

#include <string>
#include <RooArgSet.h>

class RooDataSet;
class RooWorkspace;
class RooNLLVar;

namespace EXOSTATS {

/// Utility function to unfold constraints
void        unfoldConstraints(RooArgSet &initial, RooArgSet &final, RooArgSet &obs, RooArgSet &nuis, int &counter);

/// Creates an Asimov dataset
RooDataSet *makeAsimovData(RooWorkspace *w, TString modelConfigName, Bool_t doConditional, RooNLLVar *conditioning_nll,
                           Double_t mu_val, std::string *mu_str = nullptr, std::string *mu_prof_str = nullptr,
                           Double_t mu_val_profile = -999, Bool_t doFit = true, Double_t mu_injection = -1,
                           Int_t debugLevel = 2);
} // namespace EXOSTATS

#endif
