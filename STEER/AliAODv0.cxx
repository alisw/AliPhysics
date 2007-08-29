/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
//     Analysis Oriented Data (AOD) V0 vertex class
//     Authors: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//              G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//-------------------------------------------------------------------------

//#include "AliESDEvent.h"
//#include "AliESDv0.h"

#include "AliAODv0.h"

//#include "AliAODTrack.h"

ClassImp(AliAODv0)

  AliAODv0::AliAODv0() : 
    AliAODRecoDecay(),
    fDcaV0ToPrimVertex(999)
{
  //--------------------------------------------------------------------
  // Default constructor
  //--------------------------------------------------------------------
  fCharge  = 0;
  fNProngs = 2;
  fNDCA    = 1;
  fNPID    = 2;

  fDCA = new Float_t[fNDCA];
  fDCA[0] = 999;

  fPx = new Double_t[GetNProngs()];
  fPy = new Double_t[GetNProngs()];
  fPz = new Double_t[GetNProngs()];
  fPx[0] = 999;
  fPy[0] = 999;
  fPz[0] = 999;

  fPx[1] = 999;
  fPy[1] = 999;
  fPz[1] = 999;

  fd0 = new Double_t[GetNProngs()];
  fd0[0] = 999;
  fd0[1] = 999;
}

//AliAODv0::AliAODv0(AliESDv0* rV0Vertex ,AliESDEvent* rEvent) :
//  AliAODRecoDecay(),
//  fDcaV0ToPrimVertex(999)
//{
  //--------------------------------------------------------------------
  // Constructor via fill to be removed eventually
  //--------------------------------------------------------------------
//  fCharge  = 0;
//  fNProngs = 2;
//  fNDCA    = 1;
//  fNPID    = 2;

//  fDCA = new Float_t[fNDCA];

//   fPx = new Double_t[GetNProngs()];
//   fPy = new Double_t[GetNProngs()];
//   fPz = new Double_t[GetNProngs()];

//   fd0 = new Double_t[GetNProngs()];

//   this->Fill(rV0Vertex,rEvent);
// }

AliAODv0::AliAODv0(AliAODVertex* rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
	   Double_t *rMomPos, Double_t *rMomNeg, Double_t *rDcaDaughterToPrimVertex) :
  AliAODRecoDecay(rAODVertex,2,0,rDcaDaughterToPrimVertex),
  fDcaV0ToPrimVertex(rDcaV0ToPrimVertex)
{
  //--------------------------------------------------------------------
  // Constructor via setting each data member
  //--------------------------------------------------------------------
  fCharge  = 0;
  fNProngs = 2;
  fNDCA    = 1;
  fNPID    = 2;

  fDCA = new Float_t[fNDCA];

  fDCA[0] = rDcaV0Daughters;
  fDcaV0ToPrimVertex = rDcaV0ToPrimVertex;

  fPx = new Double_t[GetNProngs()];
  fPy = new Double_t[GetNProngs()];
  fPz = new Double_t[GetNProngs()];

  fPx[0] = rMomPos[0] ;
  fPy[0] = rMomPos[1];
  fPz[0] = rMomPos[2];

  fPx[1] = rMomNeg[0];
  fPy[1] = rMomNeg[1];
  fPz[1] = rMomNeg[2];
}

AliAODv0::AliAODv0(const AliAODv0& rAliAODv0) :
  AliAODRecoDecay(rAliAODv0),
  fDcaV0ToPrimVertex(rAliAODv0.fDcaV0ToPrimVertex)
 {
  //--------------------------------------------------------------------
  // Copy constructor
  //--------------------------------------------------------------------
}

AliAODv0& AliAODv0::operator=(const AliAODv0& rAliAODv0){
  //--------------------------------------------------------------------
  // Assignment overload
  //--------------------------------------------------------------------
  this->fDcaV0ToPrimVertex  = rAliAODv0.fDcaV0ToPrimVertex ;
  return *this;
}

AliAODv0::~AliAODv0(){
  //--------------------------------------------------------------------
  // Empty destructor
  //--------------------------------------------------------------------
}


// void AliAODv0::Fill(AliESDv0* rV0Vertex ,AliESDEvent* rEvent){

//   Double_t tDecayVertexV0[3]; rV0Vertex->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 
//   fSecondaryVtx->SetX(tDecayVertexV0[0]);
//   fSecondaryVtx->SetY(tDecayVertexV0[1]);
//   fSecondaryVtx->SetZ(tDecayVertexV0[2]);

//   Double_t lCovVtx[6];
//   rV0Vertex->GetPosCov(lCovVtx);
//   fSecondaryVtx->SetCovMatrix(lCovVtx);

//   fSecondaryVtx->SetChi2perNDF(rV0Vertex->GetChi2V0());
//   fSecondaryVtx->SetType(AliAODVertex::kV0);

//   UInt_t lKeyPos = (UInt_t)TMath::Abs(rV0Vertex->GetPindex());// need to ask why Abs
//   UInt_t lKeyNeg = (UInt_t)TMath::Abs(rV0Vertex->GetNindex());
//   fSecondaryVtx->AddDaughter(rEvent->GetTrack(lKeyPos));
//   fSecondaryVtx->AddDaughter(rEvent->GetTrack(lKeyNeg));

//   fDCA[0] = rV0Vertex->GetDcaV0Daughters();
//   fDcaV0ToPrimVertex = rV0Vertex->GetD();

//   Double_t tMomPos[3]; rV0Vertex->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
//   fPx[0] = tMomPos[0];
//   fPy[0] = tMomPos[1];
//   fPz[0] = tMomPos[2];

//   Double_t tMomNeg[3]; rV0Vertex->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 
//   fPx[1] = tMomNeg[0];
//   fPy[1] = tMomNeg[1];
//   fPz[1] = tMomNeg[2];

//   AliESDtrack *pTrack=rEvent->GetTrack(lKeyPos);
//   AliESDtrack *nTrack=rEvent->GetTrack(lKeyNeg);

//   Float_t tDcaPosToPrimVertex[2];
//   if(pTrack) pTrack->GetImpactParameters(tDcaPosToPrimVertex[0],tDcaPosToPrimVertex[1]);
//   else { tDcaPosToPrimVertex[0]=999.;  tDcaPosToPrimVertex[1]=999.;}
//   fd0[0] = TMath::Sqrt(tDcaPosToPrimVertex[0]*tDcaPosToPrimVertex[0]+tDcaPosToPrimVertex[1]*tDcaPosToPrimVertex[1]);

//   Float_t tDcaNegToPrimVertex[2];
//   if(nTrack) nTrack->GetImpactParameters(tDcaNegToPrimVertex[0],tDcaNegToPrimVertex[1]);
//   else { tDcaNegToPrimVertex[0]=999.;  tDcaNegToPrimVertex[1]=999.;}

//   fd0[1] = TMath::Sqrt(tDcaNegToPrimVertex[0]*tDcaNegToPrimVertex[0]+tDcaNegToPrimVertex[1]*tDcaNegToPrimVertex[1]);
// }

void AliAODv0::Fill(AliAODVertex *rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
		    Double_t *rMomPos, Double_t *rMomNeg, Double_t *rDcaDaughterToPrimVertex){

  this->SetSecondaryVtx(rAODVertex);

  fDCA[0] = rDcaV0Daughters;
  fDcaV0ToPrimVertex = rDcaV0ToPrimVertex;

  fPx[0] = rMomPos[0] ;
  fPy[0] = rMomPos[1];
  fPz[0] = rMomPos[2];

  fPx[1] = rMomNeg[0];
  fPy[1] = rMomNeg[1];
  fPz[1] = rMomNeg[2];

  fd0[0] = rDcaDaughterToPrimVertex[0];
  fd0[1] = rDcaDaughterToPrimVertex[1];
}

void AliAODv0::ResetV0(){

  fSecondaryVtx->SetX(999);
  fSecondaryVtx->SetY(999);
  fSecondaryVtx->SetZ(999);
  fSecondaryVtx->SetChi2perNDF(999);
  fSecondaryVtx->SetType(AliAODVertex::kUndef);

  Int_t lNumDaughters = fSecondaryVtx->GetNDaughters();
  for(Int_t iDaughterIndex = 0; iDaughterIndex<lNumDaughters;iDaughterIndex++)
    fSecondaryVtx->RemoveDaughter(fSecondaryVtx->GetDaughter(iDaughterIndex));

  fDCA[0] = 999;
  fDcaV0ToPrimVertex  = 999;

  fPx[0] = 999;
  fPy[0] = 999;
  fPz[0] = 999;

  fPx[1] = 999;
  fPy[1] = 999;
  fPz[1] = 999;

  fd0[0] = 999;
  fd0[1] = 999;
}

void AliAODv0::Print(Option_t* /*option*/) const {
  //
  // Print some information
  //
  AliAODRecoDecay::Print();
  printf("AliAODv0: invariant mass (k0s %.6f, lambda %.6f, anti-lambda %.6f) \n",MassK0Short(),MassLambda(),MassAntiLambda());
  printf("AliAODv0: dca (v0d %.6f, v0tpv %.6f, postpv %.6f, negtpv %.6f ) \n",DcaV0Daughters(),DcaV0ToPrimVertex(),DcaPosToPrimVertex(),DcaNegToPrimVertex());
  printf("AliAODv0: mom (ptot2 %.6f, pt2 %.6f, rapk0 %.6f, rapla %.6f ) \n",Ptot2V0(),Pt2V0(),RapK0Short(),RapLambda());
  printf("AliAODv0: cin (mpav0 %.6f, mnav0 %.6f, alpha %.6f, ptarm %.6f ) \n",MomPosAlongV0(),MomNegAlongV0(),AlphaV0(),PtArmV0());
  printf("AliAODv0: nrg (eppro %.6f, enpro %.6f, eppio %.6f, enpio %.6f ) \n",EPosProton(),ENegProton(),EPosPion(),ENegPion());

  return;
}
