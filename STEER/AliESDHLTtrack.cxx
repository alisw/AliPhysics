/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//           Implementation of the ESD HLT track class
//   ESD = Event Summary Data
//   HLT = High Level Trigger
//   This is the class to deal with during the phisical analysis of data
//-----------------------------------------------------------------

#include "TMath.h"
#include "AliESDHLTtrack.h"

ClassImp(AliESDHLTtrack)

AliESDHLTtrack::AliESDHLTtrack() : TObject()
{
  fNHits = 0;
  fMCid = 0;
  fWeight = 0;
  fFromMainVertex = kFALSE;
  fRowRange[0] = fRowRange[1] = 0;
  fSector = 0;
  fFirstPoint[0] = fFirstPoint[1] = fFirstPoint[2] = 0;
  fLastPoint[0] = fLastPoint[1] = fLastPoint[2] = 0;
  fQ = 0;
  fTanl = 0;
  fPsi = 0;
  fPt = 0;
  fPterr = 0;
  fPsierr = 0;
  fTanlerr = 0;
  fBinX = 0;
  fBinY = 0;
  fSizeX = 0;
  fSizeY = 0;
  fPID =0;
}

Double_t AliESDHLTtrack::GetP() const
{
  // Returns total momentum.  
  return TMath::Abs(GetPt())*sqrt(1. + GetTgl()*GetTgl());
}

Double_t AliESDHLTtrack::GetPseudoRapidity() const
{
  return 0.5 * TMath::Log((GetP() + GetPz()) / (GetP() - GetPz()));
}
