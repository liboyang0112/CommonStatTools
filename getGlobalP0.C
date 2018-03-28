#include <iostream>
#include <TMath.h>
#include <RooStats/RooStatsUtils.h>

// macro for global p-value calculation
Double_t getGlobalP0(Double_t maximumSignificance, UInt_t nCrossings)
{

   // reference: https://cds.cern.ch/record/1375842/files/ATL-PHYS-PUB-2011-011.pdf
   // CMS: https://cds.cern.ch/record/1406347/files/HIG-11-032-pas.pdf
   // p (global) = p0 (local) + N*exp( - (q0 - q_ref)/2.)
   // q0: test-statistics for local max significance, q0=Z0*Z0 asymptotically
   // q_ref: reference point, eg. 1sigma
   // N: up-crossing points at reference point q_ref (down corssing if looking at p-value plot)

   const Double_t p0 = RooStats::SignificanceToPValue(maximumSignificance);
   const Double_t q0 = TMath::Power(maximumSignificance, 2);

   const Double_t pglobal = ((Double_t)nCrossings) * TMath::Exp(-(q0 - 1.) / 2.) + p0;

   const Double_t nsigglobal = RooStats::PValueToSignificance(pglobal);

   std::cout << " Local p-value :  " << p0 << " significance: " << maximumSignificance << std::endl;
   std::cout << " Global p-value:  " << pglobal << " significance: " << nsigglobal << std::endl;

   return pglobal;
}

Double_t getGlobalSignificance(Double_t maximumSignificance, UInt_t nCrossings)
{
   return RooStats::PValueToSignificance(getGlobalP0(maximumSignificance, nCrossings));
}
