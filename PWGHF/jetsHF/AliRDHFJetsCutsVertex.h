#ifndef ALIRDHFJETSCUTSVERTEX_H
#define ALIRDHFJETSCUTSVERTEX_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFJetsCuts
// base class for cuts on AOD reconstructed heavy-flavour decays
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include <TString.h>

#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPidHF.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliRDHFJetsCuts.h"

class AliAODTrack;
class AliAODRecoDecayHF;
class AliESDVertex;

class AliRDHFJetsCutsVertex : public AliRDHFJetsCuts {
 public:


  AliRDHFJetsCutsVertex(const Char_t* name="RDHFJetsCutsVertex", const Char_t* title="");
  
  virtual ~AliRDHFJetsCutsVertex();
  
  AliRDHFJetsCutsVertex(const AliRDHFJetsCutsVertex& source);
  AliRDHFJetsCutsVertex& operator=(const AliRDHFJetsCutsVertex& source);

  Int_t IsVertexSelected(AliAODVertex* vert, AliAODEvent* aod, Double_t magzkG ,Double_t dispersion,Double_t massParticle=0.138);
  void IsElecInVert(AliAODEvent *aod, AliAODVertex* vtx, Bool_t &fFlagElec);
  static Double_t CosPointingAngle(AliAODVertex* vtx1, Double_t point[3]);

  //setters
  void SetNprongs(Int_t n){fNprongs=n;}
  void SetSecVtxWithKF(Bool_t flag=kFALSE){fSecVtxWithKF=flag;}
  void SetMinPtHardestTrack(Double_t minpt){fMinPtHardestTrack=minpt;}
  void SetImpParCut(Double_t d0){fImpPar=d0;}
  void SetDistPrimSec(Double_t d){fDistPrimSec=d;}
  void SetCospCut(Double_t c){fCosp=c;}
  void SetInvMassCut(Double_t mass){fInvMassCut=mass;}
  void SetSigmaVert(Double_t sig){fSigvert=sig;}
  void SetChi2(Double_t chi){fChi2=chi;}
  void SetIsElec(Bool_t eflag=kFALSE){fIsElec=eflag;}

  //getters
 
  Int_t GetNprongs(){return fNprongs;}
  Bool_t GetSecVtxWithKF(){return fSecVtxWithKF;}
  Double_t GetMinPtHardestTrack(){return fMinPtHardestTrack;}
  Double_t GetImpParCut(){return fImpPar;}
  Double_t GetDistPrimSec(){return fDistPrimSec;}
  Double_t GetCosp(){return fCosp;}
  Double_t GetInvMassCut(){return fInvMassCut;}
  Double_t GetSigmaVert(){return fSigvert;}
  Double_t GetVertChi2(){return fChi2;}
  Double_t GetIsElec(){return fIsElec;}
 
 protected:

  Int_t fNprongs;//number of prongs
  Bool_t fSecVtxWithKF;//secondary vertex finder w/ KF
  Double_t fMinPtHardestTrack;//minimum pT of the hardest track in the vertex
  Double_t fImpPar; //Impact parameter single tracks
  Double_t fDistPrimSec; //distance primary-secondary vertex
  Double_t fCosp;//cos theta point
  Double_t fInvMassCut;
  Double_t fSigvert;
  Double_t fChi2;
  Bool_t fIsElec;

  ClassDef(AliRDHFJetsCutsVertex,1);  // base class for cuts on AOD reconstructed heavy-flavour decays
};

#endif
