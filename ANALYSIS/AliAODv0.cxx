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
//
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//
//     Purpose: Having observables for physics available for V0s
//-------------------------------------------------------------------------
#include <Riostream.h>
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
  fEvent   = 0;
}

AliAODv0::AliAODv0(AliESDv0* rV0Vertex ,AliESD* rEvent){
  this->Fill(rV0Vertex,rEvent);
}

// AliAODv0::~AliAODv0(){
// }

void AliAODv0::Fill(AliESDv0* rV0Vertex ,AliESD* rEvent){// Filling method
  fEvent=rEvent;
  Double_t tDecayVertexV0[3]; rV0Vertex->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 
  fDecayVertexV0X = tDecayVertexV0[0];
  fDecayVertexV0Y = tDecayVertexV0[1];
  fDecayVertexV0Z = tDecayVertexV0[2];

  fDcaV0Daughters = rV0Vertex->GetDcaDaughters();

  fDcaV0ToPrimVertex = rV0Vertex->GetD();


  Double_t tMomPos[3]; rV0Vertex->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
  fMomPosX = tMomPos[0];
  fMomPosY = tMomPos[1];
  fMomPosZ = tMomPos[2];

  Double_t tMomNeg[3]; rV0Vertex->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 
  fMomNegX = tMomNeg[0];
  fMomNegY = tMomNeg[1];
  fMomNegZ = tMomNeg[2];

  fKeyPos = TMath::Abs(rV0Vertex->GetPindex());// need to ask why Abs
  fKeyNeg = TMath::Abs(rV0Vertex->GetNindex());

  AliESDtrack *pTrack=fEvent->GetTrack(fKeyPos);
  AliESDtrack *nTrack=fEvent->GetTrack(fKeyNeg);

  Float_t tDcaPosToPrimVertex[2];
  if(pTrack) pTrack->GetImpactParameters(tDcaPosToPrimVertex[0],tDcaPosToPrimVertex[1]);
  else { tDcaPosToPrimVertex[0]=999.;  tDcaPosToPrimVertex[1]=999.;}

  fDcaPosToPrimVertex = TMath::Sqrt(tDcaPosToPrimVertex[0]*tDcaPosToPrimVertex[0]+tDcaPosToPrimVertex[1]*tDcaPosToPrimVertex[1]);

  Float_t tDcaNegToPrimVertex[2];
  if(nTrack) nTrack->GetImpactParameters(tDcaNegToPrimVertex[0],tDcaNegToPrimVertex[1]);
  else { tDcaNegToPrimVertex[0]=999.;  tDcaNegToPrimVertex[1]=999.;}

  fDcaNegToPrimVertex = TMath::Sqrt(tDcaNegToPrimVertex[0]*tDcaNegToPrimVertex[0]+tDcaNegToPrimVertex[1]*tDcaPosToPrimVertex[1]);
}

void AliAODv0::ResetV0(){// Reset method
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
