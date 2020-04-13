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
  fPt_mother(-1.),
  fk0(-1.),
  fk1(-1.),
  fk2(-1.),
  fkT(-1.)
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
  fPt_mother(source.fPt_mother),
  fk0(source.fk0),
  fk1(source.fk1),
  fk2(source.fk2),
  fkT(source.fkT)
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
  fPt_mother=-1.;
  fk0=-1.;
  fk1=-1.;
  fk2=-1.;
  fkT=-1.;

}
