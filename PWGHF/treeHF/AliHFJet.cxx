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
  fID(-9.),
  fHFMeson(-9.),
  fPt(-9.),
  fEta(-9.),
  fPhi(-9.),
  fDeltaEta(-9.),
  fDeltaPhi(-9.),
  fDeltaR(-9.),
  fN(-9.),
  fZg(-9.),
  fRg(-9.)
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
  fID =-9;
  fHFMeson=-9;
  fPt = -9;
  fEta = -9;
  fPhi = -9;
  fDeltaEta = -9;
  fDeltaPhi = -9;
  fDeltaR = -9;
  fN=-9;
  fZg=-9;
  fRg=-9;

}
