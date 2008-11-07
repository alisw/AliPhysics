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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour 4-prong decay
//
// Authors: G.E.Bruno Giuseppe.Bruno@to.infn.it, R.Romita Rossella.Romita@ba.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "TRandom.h" // for the time being

ClassImp(AliAODRecoDecayHF4Prong)

//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong() :
  AliAODRecoDecayHF(), 
  //fSigmaVert(0),
  fDist12toPrim(0),
  fDist23toPrim(0),
  fDist14toPrim(0),
  fDist34toPrim(0)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, //Double_t sigvert,
						 Double_t dist12,Double_t dist23,
						 Double_t dist14,Double_t dist34,
						 Short_t charge) :
  AliAODRecoDecayHF(vtx2,4,charge,px,py,pz,d0,d0err),
  // fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist23toPrim(dist23),
  fDist14toPrim(dist14),
  fDist34toPrim(dist34)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  SetDCAs(6,dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, //Double_t sigvert,
						 Double_t dist12,Double_t dist23, 
						 Double_t dist14,Double_t dist34, 
						 Short_t charge) :
  AliAODRecoDecayHF(vtx2,4,charge,d0,d0err),
  //fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist23toPrim(dist23),
  fDist14toPrim(dist14),
  fDist34toPrim(dist34)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  SetDCAs(6,dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(const AliAODRecoDecayHF4Prong &source) :
  AliAODRecoDecayHF(source),
  //fSigmaVert(source.fSigmaVert),
  fDist12toPrim(source.fDist12toPrim),
  fDist23toPrim(source.fDist23toPrim),
  fDist14toPrim(source.fDist14toPrim),
  fDist34toPrim(source.fDist34toPrim)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong &AliAODRecoDecayHF4Prong::operator=(const AliAODRecoDecayHF4Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF::operator=(source);

  fDist12toPrim= source.fDist12toPrim;
  fDist23toPrim= source.fDist23toPrim;
  fDist14toPrim= source.fDist14toPrim;
  fDist34toPrim= source.fDist34toPrim;
  //fSigmaVert= source.fSigmaVert;

  return *this;
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF4Prong::SelectD0(const Double_t *cuts,Int_t &okD0,Int_t &okD0bar)
  const {
//
// This function compares the D0 with a set of cuts:
// 
// to be implemented 
//
// cuts[1] = ... 
// cuts[2] = ...
// cuts[3] = ...
// cuts[4] = ...
//
// If candidate D0 does not pass the cuts return kFALSE
//

  okD0=0; okD0bar=0;
  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mD0=InvMassD0();
  if(TMath::Abs(mD0-mD0PDG)>cuts[0]) return kFALSE;

  //single track
  //if(TMath::Abs(PtProng(2)) < cuts[2] || TMath::Abs(Getd0Prong(2))<cuts[4])return kFALSE;//Pion2

  //DCA
  //for(Int_t i=0;i<3;i++) if(cuts[11]>0 && GetDCA(i)>cuts[11])return kFALSE;

  //2track cuts
  //if(fDist12toPrim<cuts[5] || fDist23toPrim<cuts[5])return kFALSE;
  //if(Getd0Prong(0)*Getd0Prong(1)<0 && Getd0Prong(2)*Getd0Prong(1)<0)return kFALSE;

  //sec vert
  //if(fSigmaVert>cuts[6])return kFALSE;

  //if(DecayLength()<cuts[7])return kFALSE;
  
  return kFALSE; // for the time being 
  // return kTRUE; // after implementation 
}
