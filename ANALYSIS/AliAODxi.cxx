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

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) Xi vertex class
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//     Purpose: Having observables for physics available for Xis
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliESD.h"
#include "AliAODv0.h"
#include "AliAODxi.h"

ClassImp(AliAODxi)

AliAODxi::AliAODxi() : AliAODv0() {
  //--------------------------------------------------------------------
  // Default constructor
  //--------------------------------------------------------------------
  fCharge                  = 0;
  fDecayVertexXiX          = 999;
  fDecayVertexXiY          = 999;
  fDecayVertexXiZ          = 999;
  fDcaXiDaughters          = 999;
  fDcaXiToPrimVertex       = 999;
  fDcaBachelorToPrimVertex = 999;
  fMomBachelorX            = 999;
  fMomBachelorY            = 999;
  fMomBachelorZ            = 999;

  fKeyBachelor             = 999;

  fChi2Xi                  = 999;
}

AliAODxi::AliAODxi(AliESDcascade* rXiVertex ,AliESD* rEvent){
  this->Fill(rXiVertex,rEvent);
}

AliAODxi::AliAODxi(const AliAODxi &rAliAODxi) : AliAODv0(rAliAODxi){
  //--------------------------------------------------------------------
  // Copy constructor
  //--------------------------------------------------------------------
  fDecayVertexXiX     = rAliAODxi.fDecayVertexXiX;
  fDecayVertexXiY     = rAliAODxi.fDecayVertexXiY;
  fDecayVertexXiZ     = rAliAODxi.fDecayVertexXiZ;
  fDcaXiDaughters     = rAliAODxi.fDcaXiDaughters;
  fDcaXiToPrimVertex  = rAliAODxi.fDcaXiToPrimVertex ;
  fDcaBachelorToPrimVertex = rAliAODxi.fDcaBachelorToPrimVertex;
  fMomBachelorX = rAliAODxi.fMomBachelorX;
  fMomBachelorY = rAliAODxi.fMomBachelorY;
  fMomBachelorZ = rAliAODxi.fMomBachelorZ;
  fKeyBachelor  = rAliAODxi.fKeyBachelor;
  fChi2Xi  = rAliAODxi.fChi2Xi;
}

AliAODxi& AliAODxi::operator=(const AliAODxi& rAliAODxi){
  //--------------------------------------------------------------------
  // Assignment overload
  //--------------------------------------------------------------------
  AliAODv0::operator=(rAliAODxi);
  this->fDecayVertexXiX     = rAliAODxi.fDecayVertexXiX;
  this->fDecayVertexXiY     = rAliAODxi.fDecayVertexXiY;
  this->fDecayVertexXiZ     = rAliAODxi.fDecayVertexXiZ;
  this->fDcaXiDaughters     = rAliAODxi.fDcaXiDaughters;
  this->fDcaXiToPrimVertex  = rAliAODxi.fDcaXiToPrimVertex;
  this->fDcaBachelorToPrimVertex = rAliAODxi.fDcaBachelorToPrimVertex;
  this->fMomBachelorX = rAliAODxi.fMomBachelorX;
  this->fMomBachelorY = rAliAODxi.fMomBachelorY;
  this->fMomBachelorZ = rAliAODxi.fMomBachelorZ;
  this->fKeyBachelor  = rAliAODxi.fKeyBachelor;
  this->fChi2Xi  = rAliAODxi.fChi2Xi;
  return *this;
}

AliAODxi::~AliAODxi(){
  //--------------------------------------------------------------------
  // Empty destructor
  //--------------------------------------------------------------------
}

void AliAODxi::Fill(AliESDcascade* rXiVertex ,AliESD* rEvent){
  // Fills the data memebers of the AOD
  Double_t tDecayVertexXi[3]; rXiVertex->GetXYZcascade(tDecayVertexXi[0],tDecayVertexXi[1],tDecayVertexXi[2]); 
  fDecayVertexXiX = tDecayVertexXi[0];
  fDecayVertexXiY = tDecayVertexXi[1];
  fDecayVertexXiZ = tDecayVertexXi[2];

  fDcaXiDaughters = rXiVertex->GetDcaXiDaughters();
  fDcaXiToPrimVertex = rXiVertex->GetD();

  Double_t tMomPos[3]; rXiVertex->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
  fMomPosX = tMomPos[0];
  fMomPosY = tMomPos[1];
  fMomPosZ = tMomPos[2];

  Double_t tMomNeg[3]; rXiVertex->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 
  fMomNegX = tMomNeg[0];
  fMomNegY = tMomNeg[1];
  fMomNegZ = tMomNeg[2];

  Double_t tMomBachelor[3]; rXiVertex->GetBPxPyPz(tMomBachelor[0],tMomBachelor[1],tMomBachelor[2]); 
  fMomBachelorX = tMomBachelor[0];
  fMomBachelorY = tMomBachelor[1];
  fMomBachelorZ = tMomBachelor[2];

  fKeyPos = (UInt_t)TMath::Abs(rXiVertex->GetPindex());// need to ask why Abs
  fKeyNeg = (UInt_t)TMath::Abs(rXiVertex->GetNindex());
  fKeyBachelor = (UInt_t)TMath::Abs(rXiVertex->GetBindex());

  AliESDtrack *pTrack=rEvent->GetTrack(fKeyPos);
  AliESDtrack *nTrack=rEvent->GetTrack(fKeyNeg);
  AliESDtrack *bTrack=rEvent->GetTrack(fKeyBachelor);

  Float_t tDcaPosToPrimVertex[2];
  if(pTrack) pTrack->GetImpactParameters(tDcaPosToPrimVertex[0],tDcaPosToPrimVertex[1]);
  else { tDcaPosToPrimVertex[0]=999.;  tDcaPosToPrimVertex[1]=999.;}

  fDcaPosToPrimVertex = TMath::Sqrt(tDcaPosToPrimVertex[0]*tDcaPosToPrimVertex[0]+tDcaPosToPrimVertex[1]*tDcaPosToPrimVertex[1]);

  Float_t tDcaNegToPrimVertex[2];
  if(nTrack) nTrack->GetImpactParameters(tDcaNegToPrimVertex[0],tDcaNegToPrimVertex[1]);
  else { tDcaNegToPrimVertex[0]=999.;  tDcaNegToPrimVertex[1]=999.;}


  fDcaNegToPrimVertex = TMath::Sqrt(tDcaNegToPrimVertex[0]*tDcaNegToPrimVertex[0]+tDcaNegToPrimVertex[1]*tDcaNegToPrimVertex[1]);

  Float_t tDcaBachelorToPrimVertex[2];
  if(bTrack) bTrack->GetImpactParameters(tDcaBachelorToPrimVertex[0],tDcaBachelorToPrimVertex[1]);
  else { tDcaBachelorToPrimVertex[0]=999.;  tDcaBachelorToPrimVertex[1]=999.;}

  fDcaBachelorToPrimVertex = TMath::Sqrt(tDcaBachelorToPrimVertex[0]*tDcaBachelorToPrimVertex[0]+tDcaBachelorToPrimVertex[1]*tDcaBachelorToPrimVertex[1]);

  fDcaV0Daughters    = rXiVertex->GetDcaV0Daughters();
  fDcaV0ToPrimVertex = rXiVertex->GetDcascade(rEvent->GetVertex()->GetXv(),
					      rEvent->GetVertex()->GetYv(),
					      rEvent->GetVertex()->GetZv());

  double tDecayVertexV0[3]; rXiVertex->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
  fDecayVertexV0X=tDecayVertexV0[0];
  fDecayVertexV0Y=tDecayVertexV0[1];
  fDecayVertexV0Z=tDecayVertexV0[2];
}

void AliAODxi::ResetXi(){
  // Sets the default values of the AOD data members
  fDecayVertexXiX          = 999;
  fDecayVertexXiY          = 999;
  fDecayVertexXiZ          = 999;
  fDcaXiDaughters          = 999;
  fDcaXiToPrimVertex       = 999;
  fDcaBachelorToPrimVertex = 999;
  fMomBachelorX = 999;
  fMomBachelorY = 999;
  fMomBachelorZ = 999;

  fKeyBachelor  = 999;
  fChi2Xi       = 999;
}
