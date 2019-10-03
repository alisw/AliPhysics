/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJet
// \helper class to handle jet objects
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "TMath.h"
#include "AliHFJet.h"

/// \cond CLASSIMP
ClassImp(AliHFJet);
/// \endcond

//________________________________________________________________
// Definitions of class AliHFJet

/// \cond CLASSIMP
ClassImp(AliHFJet);
/// \endcond


AliHFJet::AliHFJet():
  fID(-99.),
  fHFMeson(-99.),
  fPt(-99.),
  fEta(-99.),
  fPhi(-99.),
  fDeltaEta(-99.),
  fDeltaPhi(-99.),
  fDeltaR(-99.),
  fN(-99.),
  fZg(-99.),
  fRg(-99.)
{
}


//________________________________________________________________

AliHFJet::AliHFJet(const AliHFJet &source):
  fID(source.fID),
  fHFMeson(source.fHFMeson),
  fPt(source.fPt),
  fEta(source.fEta),
  fPhi(source.fPhi),
  fDeltaEta(source.fDeltaEta),
  fDeltaPhi(source.fDeltaPhi),
  fDeltaR(source.fDeltaR),
  fN(source.fN),
  fZg(source.fZg),
  fRg(source.fRg)
{
}

//________________________________________________________________
AliHFJet::~AliHFJet()
{
  //
  // Destructor
  //
}

//________________________________________________________________
// Reset the current object
void AliHFJet::Reset()
{
  fID =-99.;
  fHFMeson=-99.;
  fPt = -99.;
  fEta = -99.;
  fPhi = -99.;
  fDeltaEta = -99.;
  fDeltaPhi = -99.;
  fDeltaR = -99.;
  fN=-99.;
  fZg=-99.;
  fRg=-99.;

}
