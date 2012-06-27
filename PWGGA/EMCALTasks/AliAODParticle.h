#ifndef ALIAODPARTICLE_H
#define ALIAODPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */
 
//_________________________________________________________________________
// Class to fill photon and track for correlation 

#include "TLorentzVector.h"

class AliAODParticle :public TLorentzVector{
  
 public:
  
  AliAODParticle() ;
  AliAODParticle(Double_t px,Double_t py,Double_t pz,Double_t E) ; 
  ~AliAODParticle(){} 

   Int_t    GetChargedSign() const { return fChargedSign ;}

   Float_t  GetLambda0(void)        const {return fL0      ;}
   Float_t  GetLambda1(void)        const {return fL1      ;}
   Int_t    GetModule(void)         const {return fModule  ;}
   Double_t GetDistBad(void)        const {return fBadDist ;}
   Int_t    GetNCells(void)         const {return fNCells  ;}
   Double_t GetTOF(void)            const {return fClusterTime;}
   Int_t    GetClusterID(void)      const {return fClusterID;}
   Int_t    GetAODClusterID(void)   const {return fAODClusterID;}

   Float_t  GetPhotonPairDTime(void)    const {return fPhotonPairDTime   ;} 
   Int_t    GetPhotonPairDModule(void)  const {return fPhotonPairDModule ;} 
   Float_t  GetPhotonPairAsy(void)      const {return fPhotonPairAsy     ;}
   Float_t  GetPhotonPairAngle(void)    const {return fPhotonPairAngle   ;}
   Int_t    GetPhotonPairID(Int_t i)    const {if(i ==0 ) return fPhotonPairID0;   else if(i ==1) return fPhotonPairID1; else return -999;}
   Int_t    GetAODPhotonPairID(Int_t i)    const {if(i ==0 ) return fAODPhotonPairID0;   else if(i ==1) return fAODPhotonPairID1; else return -999;}

   Bool_t   IsInSSA(void)          const {return kInSSA           ;}
   Bool_t   IsInTrackMatched(void) const {return kInTrackMatched  ;}
   Bool_t   IsInTOF(void)          const {return kInTOF           ;}
   Bool_t   IsLeading(void)        const {return kIsLeading       ;}
   Bool_t   IsIsolated(void)       const {return kIsIsolated      ;}

   Bool_t   IsPIDOK(const Int_t fPid) const ;

   void SetChargedSign(Int_t n){fChargedSign = n;}
   
   void SetLambdas(Float_t l0,Float_t l1){fL0=l0,fL1=l1;}
   void SetModule(Int_t mod){fModule = mod;}
   void SetDistBad(Double_t nbad){fBadDist = nbad;}
   void SetNCells(Int_t ncell){fNCells = ncell;}
   void SetTOF(Double_t time){fClusterTime = time;}
   void SetClusterID(Int_t id){fClusterID = id;}
   void SetAODClusterID(Int_t aodid){fAODClusterID = aodid;}
   void SetSSABit(Bool_t ssa){kInSSA = ssa;} 
   void SetTOFBit(Bool_t tof){kInTOF = tof;} 
   void SetTrackMatchedBit(Bool_t matched){kInTrackMatched = matched;}
   void SetDistToBad(Int_t dist){fBadDist=dist;} 
   void SetPhotonPairDTime(Float_t pairtime){fPhotonPairDTime = pairtime;}
   void SetPhotonPairDModule(Int_t pairmodule){fPhotonPairDModule = pairmodule;}
   void SetPhotonPairAsy(Float_t asy){fPhotonPairAsy = asy;}
   void SetPhotonPairAngle(Float_t angle){fPhotonPairAngle = angle;}
   void SetPhotonPairID(Int_t id0, Int_t id1){fPhotonPairID0 =id0, fPhotonPairID1=id1 ; }
   void SetAODPhotonPairID(Int_t aodid0, Int_t aodid1){fAODPhotonPairID0 =aodid0, fAODPhotonPairID1=aodid1 ; }
   void SetIsLeading(Bool_t leading){kIsLeading = leading;}
   void SetIsolated(Bool_t iso){kIsIsolated = iso;}
private:

  Int_t     fChargedSign ;   //sign of charged track

  Float_t   fL0               ;
  Float_t   fL1               ;
  Int_t     fModule           ; 
  Int_t     fBadDist          ;
  Int_t     fNCells           ;
  Double_t  fClusterTime      ;
  Int_t     fClusterID        ;
  Int_t     fAODClusterID        ;
  Bool_t    kInSSA            ;   
  Bool_t    kInTOF            ;    
  Bool_t    kInTrackMatched   ;   
  Double_t  fPhotonPairDTime  ;
  Int_t     fPhotonPairDModule;
  Float_t   fPhotonPairAsy    ;
  Float_t   fPhotonPairAngle  ;
  Int_t     fPhotonPairID0    ;
  Int_t     fAODPhotonPairID0    ;
  Int_t     fPhotonPairID1    ;
  Int_t     fAODPhotonPairID1    ;
  Bool_t    kIsLeading        ;
  Bool_t    kIsIsolated       ;

  ClassDef(AliAODParticle,2)

};

#endif // #ifdef ALIAODPARTICLE_H

  
