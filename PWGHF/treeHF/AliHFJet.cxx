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
  fID(-1.),
  fHFMeson(-1.),
  fPt(-1.),
  fEta(-1.),
  fPhi(-1.),
  fDeltaEta(-1.),
  fDeltaPhi(-1.),
  fDeltaR(-1.),
  fN(-1.),
  fZg(-1.),
  fRg(-1.),
  fNsd(-1.),
  fPt_splitting(-1.),
  fk0(-1.),
  fZk0(-1.),
  fRk0(-1.),
  fk1(-1.),
  fZk1(-1.),
  fRk1(-1.),
  fk2(-1.),
  fZk2(-1.),
  fRk2(-1.),
  fkT(-1.),
  fZkT(-1.),
  fRkT(-1.)
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
  fRg(source.fRg),
  fNsd(source.fNsd),
  fPt_splitting(source.fPt_splitting),
  fk0(source.fk0),
  fZk0(source.fZk0),
  fRk0(source.fRk0),
  fk1(source.fk1),
  fZk1(source.fZk1),
  fRk1(source.fRk1),
  fk2(source.fk2),
  fZk2(source.fZk2),
  fRk2(source.fRk2),
  fkT(source.fkT),
  fZkT(source.fZkT),
  fRkT(source.fRkT)
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
  fID =-1.;
  fHFMeson=-1.;
  fPt = -1.;
  fEta = -1.;
  fPhi = -1.;
  fDeltaEta = -1.;
  fDeltaPhi = -1.;
  fDeltaR = -1.;
  fN=-1.;
  fZg=-1.;
  fRg=-1.;
  fNsd=-1.;
  fPt_splitting=-1.;
  fk0=-1.;
  fZk0=-1.;
  fRk0=-1.;
  fk1=-1.;
  fZk1=-1.;
  fRk1=-1.;
  fk2=-1.;
  fZk2=-1.;
  fRk2=-1.;
  fkT=-1.;
  fZkT=-1.;
  fRkT=-1.;

}
