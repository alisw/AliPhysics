/*
***********************************************************
  Implementation of AliReducedPairInfo class.
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDPAIRINFO_H
#include "AliReducedPairInfo.h"
#endif

ClassImp(AliReducedPairInfo)

const Char_t* AliReducedPairInfo::fgkDecayChannelNames[AliReducedPairInfo::kNMaxCandidateTypes][4] = 
{
   "GammaConv",               "#gamma #rightarrow e^{+} e^{-}",                     "e^{+}", "e^{-}",
   "K0sToPiPi",                    "K^{0}_{s} #rightarrow #pi^{+} #pi^{-}",           "#pi^{+}", "#pi^{-}",
   "Lambda0ToPPi",            "#Lambda^{0} #rightarrow p #pi^{-}",                  "p", "#pi^{-}",
   "ALambda0ToPPi",          "#bar{#Lambda} #rightarrow #bar{p} #pi^{+}", "#bar{p}", "#pi^{+}",
   "PhiToKK",                       "#phi #rightarrow K^{+} K^{-}",                            "K^{+}", "K^{-}",
   "JpsiToEE",                      "J/#psi #rightarrow e^{+} e^{-}",                           "e^{+}", "e^{-}",
   "Upsilon",                        "#Upsilon #rightarrow e^{+} e^{-}",                      "e^{+}", "e^{-}",
   "DplusToK0sPiplus",        "D^{+} #rightarrow K^{0}_{s} #pi^{+}",            "K^{0}_{s}", "#pi^{+}",
   "DplusToK0sKplus",         "D^{+} #rightarrow K^{0}_{s} K^{+}",               "K^{0}_{s}", "K^{+}",
   "DplusToPhiPiplus",         "D^{+} #rightarrow #phi #pi^{+}",                      "#phi", "#pi^{+}",
   "DminusToK0sPiminus",  "D^{-} #rightarrow K^{0}_{s} #pi^{-}",              "K^{0}_{s}", "#pi^{-}",
   "DminusToK0sKminus",   "D^{-} #rightarrow K^{0}_{s} K^{-}",                 "K^{0}_{s}", "K^{-}",
   "DminusToPhiPiminus",   "D^{-} #rightarrow #phi #pi^{-}",                        "#phi", "#pi^{-}",
   "ADzeroToKplusPiminus", "#bar{D^{0}} #rightarrow K^{+} #pi^{-}",       "K^{+}", "#pi^{-}",
   "DsplusToK0sKplus",         "D^{+}_{s} #rightarrow K^{0}_{s} K^{+}",     "K^{0}_{s}", "K^{+}",
   "DsminusToK0sKminus",   "D^{-}_{s} #rightarrow K^{0}_{s} K^{-}",        "K^{0}_{s}", "K^{-}", 
};

//_______________________________________________________________________________
AliReducedPairInfo::AliReducedPairInfo() :
  AliReducedBaseTrack(),
  fCandidateId(-1),
  fPairType(-1), 
  fPairTypeSPD(-1), 
  fLegIds(),
  fMass(),
  fLxy(0.0),
  fPsProper(0.0),
  fPointingAngle(0.0),
  fChisquare(0.0)
{
  //
  // Constructor
  //
  fLegIds[0] = 0; fLegIds[1] = 0;
  fMass[0]=-999.; fMass[1]=-999.; fMass[2]=-999.; fMass[3]=-999.;
}


//_______________________________________________________________________________
AliReducedPairInfo::AliReducedPairInfo(const AliReducedPairInfo &c) :
  AliReducedBaseTrack(c),
  fCandidateId(c.CandidateId()),
  fPairType(c.PairType()),
  fPairTypeSPD(c.PairTypeSPD()),
  fLegIds(),
  fMass(),
  fLxy(c.Lxy()),
  fPsProper(c.PsProper()),
  fPointingAngle(c.PointingAngle()),
  fChisquare(c.Chi2())
{
  //
  // copy constructor
  //
  fLegIds[0] = c.LegId(0);
  fLegIds[1] = c.LegId(1);
  fMass[0] = c.Mass(0); fMass[1] = c.Mass(1); fMass[2] = c.Mass(2); fMass[3] = c.Mass(3);
}


//_______________________________________________________________________________
AliReducedPairInfo::~AliReducedPairInfo() {
  //
  // destructor
  //
}
