#ifndef __EXOSTATS_MINIMIZATION_H__
#define __EXOSTATS_MINIMIZATION_H__

#include <TString.h>

class RooAbsReal;
class RooNLLVar;
class RooWorkspace;

namespace EXOSTATS {
Int_t minimize(RooNLLVar *nll, Int_t maxRetries = 3, RooWorkspace *w = nullptr,
               TString mu0Snapshot = "conditionalNuis_0", TString nominalSnapshot = "nominalNuis",
               Int_t debugLevel = 2, Bool_t saveFitResult = kFALSE, RooFitResult **fitResult = nullptr);
Int_t minimize(RooAbsReal *fcn, Int_t maxRetries = 3, RooWorkspace *w = nullptr,
               TString mu0Snapshot = "conditionalNuis_0", TString nominalSnapshot = "nominalNuis",
               Int_t debugLevel = 2, Bool_t saveFitResult = kFALSE, RooFitResult **fitResult = nullptr);
} // namespace EXOSTATS

#endif
