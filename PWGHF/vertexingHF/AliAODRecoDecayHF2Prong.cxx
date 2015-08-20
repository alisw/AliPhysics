/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour 2-prong decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"

ClassImp(AliAODRecoDecayHF2Prong)

//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong() :
  AliAODRecoDecayHF()
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayHF(vtx2,2,0,px,py,pz,d0,d0err)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayHF(vtx2,2,0,d0,d0err)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(const AliAODRecoDecayHF2Prong &source) :
  AliAODRecoDecayHF(source)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong &AliAODRecoDecayHF2Prong::operator=(const AliAODRecoDecayHF2Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF::operator=(source);

  return *this;
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF2Prong::SelectD0(const Double_t *cuts,Int_t &okD0,Int_t &okD0bar) 
  const {
//
// This function compares the D0 with a set of cuts:
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = dca [cm]
// cuts[2] = cosThetaStar 
// cuts[3] = pTK [GeV/c]
// cuts[4] = pTPi [GeV/c]
// cuts[5] = d0K [cm]   upper limit!
// cuts[6] = d0Pi [cm]  upper limit!
// cuts[7] = d0d0 [cm^2]
// cuts[8] = cosThetaPoint
//
// If the D0/D0bar doesn't pass the cuts it sets the weights to 0
// If neither D0 nor D0bar pass the cuts return kFALSE
//
  Double_t mD0,mD0bar,ctsD0,ctsD0bar;
  okD0=1; okD0bar=1;

  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

  if(PtProng(1) < cuts[3] || PtProng(0) < cuts[4]) okD0 = 0;
  if(PtProng(0) < cuts[3] || PtProng(1) < cuts[4]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(TMath::Abs(Getd0Prong(1)) > cuts[5] || 
     TMath::Abs(Getd0Prong(0)) > cuts[6]) okD0 = 0;
  if(TMath::Abs(Getd0Prong(0)) > cuts[6] ||
     TMath::Abs(Getd0Prong(1)) > cuts[5]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(GetDCA() > cuts[1]) { okD0 = okD0bar = 0; return kFALSE; }

  InvMassD0(mD0,mD0bar);
  if(TMath::Abs(mD0-mD0PDG)    > cuts[0]) okD0 = 0;
  if(TMath::Abs(mD0bar-mD0PDG) > cuts[0]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  CosThetaStarD0(ctsD0,ctsD0bar);
  if(TMath::Abs(ctsD0)    > cuts[2]) okD0 = 0;
  if(TMath::Abs(ctsD0bar) > cuts[2]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(Prodd0d0() > cuts[7]) { okD0 = okD0bar = 0; return kFALSE; }

  if(CosPointingAngle()   < cuts[8]) { okD0 = okD0bar = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF2Prong::SelectBtoJPSI(const Double_t *cuts,Int_t &okB)
  const {
//
// This function compares the Secondary JPSI candidates with a set of cuts:
//
// cuts[0] = inv. mass half width [GeV]
// cuts[1] = dca [cm]
// cuts[2] = cosThetaStar (negative electron)
// cuts[3] = pTP [GeV/c]
// cuts[4] = pTN [GeV/c]
// cuts[5] = d0P [cm]   upper limit!
// cuts[6] = d0N [cm]  upper limit!
// cuts[7] = d0d0 [cm^2]
// cuts[8] = cosThetaPoint
//
// If the candidate doesn't pass the cuts it sets the weight to 0
// and return kFALSE
//
  Double_t mJPsi,ctsJPsi;
  okB=1;

  Double_t mJPSIPDG = TDatabasePDG::Instance()->GetParticle(443)->Mass();

  if(PtProng(1) < cuts[3] || PtProng(0) < cuts[4]) okB = 0;
  if(!okB) return kFALSE;

  if(TMath::Abs(Getd0Prong(1)) > cuts[5] ||
     TMath::Abs(Getd0Prong(0)) > cuts[6]) okB = 0;
  if(!okB) return kFALSE;

  if(GetDCA() > cuts[1]) { okB = 0; return kFALSE; }

  mJPsi=InvMassJPSIee();
  if(TMath::Abs(mJPsi-mJPSIPDG)    > cuts[0]) okB = 0;
  if(!okB) return kFALSE;

  ctsJPsi=CosThetaStarJPSI();
  if(TMath::Abs(ctsJPsi)    > cuts[2]) okB = 0;
  if(!okB) return kFALSE;

  if(Prodd0d0() > cuts[7]) { okB = 0; return kFALSE; }

  if(CosPointingAngle()   < cuts[8]) { okB = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
