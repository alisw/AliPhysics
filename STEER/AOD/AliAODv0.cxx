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

#include "AliAODv0.h"
#include "AliAODTrack.h"

ClassImp(AliAODv0)

  AliAODv0::AliAODv0() : 
    AliAODRecoDecay(),
    fDcaV0ToPrimVertex(999),
    fOnFlyStatus(kFALSE)
{
  //--------------------------------------------------------------------
  // Default constructor
  //--------------------------------------------------------------------
  fCharge  = 0;
  fNProngs = 2;
  fNDCA    = 1;
  fNPID    = 0; // used to be 2!

  fDCA = new Double_t[fNDCA];
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

AliAODv0::AliAODv0(AliAODVertex* rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
	   const Double_t *rMomPos, const Double_t *rMomNeg, Double_t *rDcaDaughterToPrimVertex) :
  AliAODRecoDecay(rAODVertex,2,0,rDcaDaughterToPrimVertex),
  fDcaV0ToPrimVertex(rDcaV0ToPrimVertex),
  fOnFlyStatus(kFALSE)
{
  //--------------------------------------------------------------------
  // Constructor via setting each data member
  //--------------------------------------------------------------------
  fCharge  = 0;
  fNProngs = 2;
  fNDCA    = 1;
  fNPID    = 0; // used to be 2!

  fDCA = new Double_t[fNDCA];

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
  fDcaV0ToPrimVertex(rAliAODv0.fDcaV0ToPrimVertex),
  fOnFlyStatus(rAliAODv0.fOnFlyStatus)
 {
  //--------------------------------------------------------------------
  // Copy constructor
  //--------------------------------------------------------------------
}

AliAODv0& AliAODv0::operator=(const AliAODv0& rAliAODv0){
  //--------------------------------------------------------------------
  // Assignment overload
  //--------------------------------------------------------------------
  if(this!=&rAliAODv0) {
    AliAODRecoDecay::operator=(rAliAODv0);
    this->fDcaV0ToPrimVertex  = rAliAODv0.fDcaV0ToPrimVertex ;
    this->fOnFlyStatus        = rAliAODv0.fOnFlyStatus;
  }
  return *this;
}

AliAODv0::~AliAODv0(){
  //--------------------------------------------------------------------
  // Empty destructor
  //--------------------------------------------------------------------
}

void AliAODv0::Fill(AliAODVertex *rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
		    const Double_t *rMomPos, const Double_t *rMomNeg, const Double_t *rDcaDaughterToPrimVertex){
  //--------------------------------------------------------------------
  // Filling with all needed info
  //--------------------------------------------------------------------
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
  //--------------------------------------------------------------------
  // Resetting all the info
  //--------------------------------------------------------------------
  GetSecondaryVtx()->SetChi2perNDF(999);
  GetSecondaryVtx()->RemoveCovMatrix();
  GetSecondaryVtx()->RemoveDaughters();
  GetSecondaryVtx()->SetParent((TObject*) 0x0);
  GetSecondaryVtx()->SetID(-1);
  GetSecondaryVtx()->SetPosition(999,999,999);
  GetSecondaryVtx()->SetType(AliAODVertex::kUndef);

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

Short_t AliAODv0::GetPosID() const {
  	AliAODTrack *posTrack = (AliAODTrack *) (this->GetSecondaryVtx()->GetDaughter(0));
	Short_t posID = posTrack->GetID();
	return posID;
}

Short_t AliAODv0::GetNegID() const {
  	AliAODTrack *negTrack = (AliAODTrack *) (this->GetSecondaryVtx()->GetDaughter(1));
	Short_t negID = negTrack->GetID();
	return negID;
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
