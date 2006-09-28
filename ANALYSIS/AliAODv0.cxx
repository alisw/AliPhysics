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
//     Implementation of the Analysis Oriented Data (AOD) V0 vertex class
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//     Purpose: Having observables for physics available for V0s
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliESD.h"
#include "AliAODv0.h"

ClassImp(AliAODv0)

AliAODv0::AliAODv0() : TObject() {
  //--------------------------------------------------------------------
  // Default constructor
  //--------------------------------------------------------------------
  fDecayVertexV0X     = 999;
  fDecayVertexV0Y     = 999;
  fDecayVertexV0Z     = 999;
  fDcaV0Daughters     = 999;
  fDcaV0ToPrimVertex  = 999;
  fDcaPosToPrimVertex = 999;
  fDcaNegToPrimVertex = 999;
  fMomPosX = 999;
  fMomPosY = 999;
  fMomPosZ = 999;
  fMomNegX = 999;
  fMomNegY = 999;
  fMomNegZ = 999;

  fKeyPos  = 999;
  fKeyNeg  = 999;

  fChi2    = 999;
}

AliAODv0::AliAODv0(AliESDv0* rV0Vertex ,AliESD* rEvent){
  this->Fill(rV0Vertex,rEvent);
}


AliAODv0::AliAODv0(const AliAODv0& rAliAODv0) : TObject(rAliAODv0) {
  //--------------------------------------------------------------------
  // Copy constructor
  //--------------------------------------------------------------------
  fDecayVertexV0X     = rAliAODv0.fDecayVertexV0X;
  fDecayVertexV0Y     = rAliAODv0.fDecayVertexV0Y;
  fDecayVertexV0Z     = rAliAODv0.fDecayVertexV0Z;
  fDcaV0Daughters     = rAliAODv0.fDcaV0Daughters;
  fDcaV0ToPrimVertex  = rAliAODv0.fDcaV0ToPrimVertex ;
  fDcaPosToPrimVertex = rAliAODv0.fDcaPosToPrimVertex;
  fDcaNegToPrimVertex = rAliAODv0.fDcaNegToPrimVertex;
  fMomPosX = rAliAODv0.fMomPosX;
  fMomPosY = rAliAODv0.fMomPosY;
  fMomPosZ = rAliAODv0.fMomPosZ;
  fMomNegX = rAliAODv0.fMomNegX;
  fMomNegY = rAliAODv0.fMomNegY;
  fMomNegZ = rAliAODv0.fMomNegZ;

  fKeyNeg  = rAliAODv0.fKeyNeg;
  fKeyPos =  rAliAODv0.fKeyPos;

  fChi2    = rAliAODv0.fChi2;
}

AliAODv0& AliAODv0::operator=(const AliAODv0& rAliAODv0){
  //--------------------------------------------------------------------
  // Assignment overload
  //--------------------------------------------------------------------
  this->fDecayVertexV0X     = rAliAODv0.fDecayVertexV0X;
  this->fDecayVertexV0Y     = rAliAODv0.fDecayVertexV0Y;
  this->fDecayVertexV0Z     = rAliAODv0.fDecayVertexV0Z;
  this->fDcaV0Daughters     = rAliAODv0.fDcaV0Daughters;
  this->fDcaV0ToPrimVertex  = rAliAODv0.fDcaV0ToPrimVertex ;
  this->fDcaPosToPrimVertex = rAliAODv0.fDcaPosToPrimVertex;
  this->fDcaNegToPrimVertex = rAliAODv0.fDcaNegToPrimVertex;
  this->fMomPosX = rAliAODv0.fMomPosX;
  this->fMomPosY = rAliAODv0.fMomPosY;
  this->fMomPosZ = rAliAODv0.fMomPosZ;
  this->fMomNegX = rAliAODv0.fMomNegX;
  this->fMomNegY = rAliAODv0.fMomNegY;
  this->fMomNegZ = rAliAODv0.fMomNegZ;

  this->fKeyPos  = rAliAODv0.fKeyPos;
  this->fKeyNeg  = rAliAODv0.fKeyNeg;

  this->fChi2    = rAliAODv0.fChi2;
  return *this;
}

AliAODv0::~AliAODv0(){
  //--------------------------------------------------------------------
  // Empty destructor
  //--------------------------------------------------------------------
}


void AliAODv0::Fill(AliESDv0* rV0Vertex ,AliESD* rEvent){
  // Fills the data memebers of the AOD
  Double_t tDecayVertexV0[3]; rV0Vertex->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 
  fDecayVertexV0X = tDecayVertexV0[0];
  fDecayVertexV0Y = tDecayVertexV0[1];
  fDecayVertexV0Z = tDecayVertexV0[2];

  fDcaV0Daughters = rV0Vertex->GetDcaV0Daughters();
  fDcaV0ToPrimVertex = rV0Vertex->GetD();

  Double_t tMomPos[3]; rV0Vertex->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
  fMomPosX = tMomPos[0];
  fMomPosY = tMomPos[1];
  fMomPosZ = tMomPos[2];

  Double_t tMomNeg[3]; rV0Vertex->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 
  fMomNegX = tMomNeg[0];
  fMomNegY = tMomNeg[1];
  fMomNegZ = tMomNeg[2];

  fKeyPos = (UInt_t)TMath::Abs(rV0Vertex->GetPindex());// need to ask why Abs
  fKeyNeg = (UInt_t)TMath::Abs(rV0Vertex->GetNindex());

  AliESDtrack *pTrack=rEvent->GetTrack(fKeyPos);
  AliESDtrack *nTrack=rEvent->GetTrack(fKeyNeg);

  Float_t tDcaPosToPrimVertex[2];
  if(pTrack) pTrack->GetImpactParameters(tDcaPosToPrimVertex[0],tDcaPosToPrimVertex[1]);
  else { tDcaPosToPrimVertex[0]=999.;  tDcaPosToPrimVertex[1]=999.;}

  fDcaPosToPrimVertex = TMath::Sqrt(tDcaPosToPrimVertex[0]*tDcaPosToPrimVertex[0]+tDcaPosToPrimVertex[1]*tDcaPosToPrimVertex[1]);

  Float_t tDcaNegToPrimVertex[2];
  if(nTrack) nTrack->GetImpactParameters(tDcaNegToPrimVertex[0],tDcaNegToPrimVertex[1]);
  else { tDcaNegToPrimVertex[0]=999.;  tDcaNegToPrimVertex[1]=999.;}

  fDcaNegToPrimVertex = TMath::Sqrt(tDcaNegToPrimVertex[0]*tDcaNegToPrimVertex[0]+tDcaNegToPrimVertex[1]*tDcaNegToPrimVertex[1]);
}

void AliAODv0::ResetV0(){
  // Sets the default values of the AOD data members
  fDecayVertexV0X     = 999;
  fDecayVertexV0Y     = 999;
  fDecayVertexV0Z     = 999;
  fDcaV0Daughters     = 999;
  fDcaV0ToPrimVertex  = 999;
  fDcaPosToPrimVertex = 999;
  fDcaNegToPrimVertex = 999;
  fMomPosX = 999;
  fMomPosY = 999;
  fMomPosZ = 999;
  fMomNegX = 999;
  fMomNegY = 999;
  fMomNegZ = 999;

  fKeyPos  = 999;
  fKeyNeg  = 999;

  fChi2    = 999;
}
