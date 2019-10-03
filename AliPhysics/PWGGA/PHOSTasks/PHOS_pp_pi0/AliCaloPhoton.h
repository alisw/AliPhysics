#ifndef ALICALOPHOTON_H
#define ALICALOPHOTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */
 
//_________________________________________________________________________
// Class to fill two-photon invariant mass hisograms
// to be used to extract pi0 raw yield.
//
//-- Author: Dmitri Peressounko (RRC "KI")
// This class contains all (minimal) necessary information about photon to 
// calculate invarint mass distr for pi0
// and for tagging and isolation analysis

class AliVCluster;

#include "TLorentzVector.h"

class AliCaloPhoton :public TLorentzVector{
  
 public:
  
  AliCaloPhoton() ;
  AliCaloPhoton(Double_t px,Double_t py,Double_t pz,Double_t E) ; 
  ~AliCaloPhoton(){} 

   const TLorentzVector * GetMomV2()const{return &fMomV2;}
   Int_t    DistToBad() const {return fBadDist ;}
   Double_t DistToBadfp() const {return fBadDistfp ;}
   Double_t EMCx(void)  const {return fX;}
   Double_t EMCy(void)  const {return fY;}
   Double_t EMCz(void)  const {return fZ;}
   Int_t    Module(void)const {return fModule;}
   Int_t    GetBC(void) const {return fBC;}
   Int_t    GetFiducialArea(void) const {return fFiducialArea ;}
   Int_t    GetIsolationTag(void) const {return fIsolationTag ;}
   Double_t GetLambda1(void) const {return fLambda0;}
   Double_t GetLambda2(void) const {return fLambda1;}
   Int_t    GetNCells() const { return fNCells ;} 
   Int_t    GetPrimary()const {return fPrimary;}
   Int_t    GetPrimaryAtVertex()  const {return fPrimaryAtVertex;}
   Double_t GetPartnerPt(void)    const {return fPartnerPt;}  
   Int_t    GetTagInfo(void) const {return fTagInfo;}
   Double_t GetTime(void)    const {return fTime ;}
   Double_t GetWeight(void)  const {return fWeight;}

   Int_t    IsConvertedPartner() const { if(fConvertedPartner == 1) return 1; else return 0; }
   Bool_t   IsCPVOK(void)   const {return fCpv;}
   Bool_t   IsCPV2OK(void)  const {return fCpv2;}
   Bool_t   IsDispOK(void)  const {return fDisp;}
   Bool_t   IsDisp2OK(void) const {return fDisp2;} //stricter cut
   Bool_t   IsIsolated(void)const {return fIsIsolated ;}
   Bool_t   IsPhoton() const {return fIsPhoton ;} //check if this particle is indeed photon (this bit is set with MC stack info
   Bool_t   IsPIDOK(const Int_t ipid) const ;
   Bool_t   IsTagged(void)  const {return fIsTagged ;} //check if this photon is tagged
   Bool_t   IsTagged(Int_t i,Int_t k) const {return fIsTagged_reg[i][k] ;} //check if this photon is tagged
   Bool_t   IsTOFOK(void)   const {return fTof;}
   Bool_t   IsTrig(void)    const{ return fTrig ; }
   Bool_t   IsntUnfolded(void)const{return fUnfolded;}

   Double_t GetNsigmaCPV()      {return fNsigmaCPV;}
   Double_t GetNsigmaFullDisp() {return fNsigmaFullDisp;}
   Double_t GetNsigmaCoreDisp() {return fNsigmaCoreDisp;}

   Double_t GetTOFCutEfficiency() {return fTOFCutEfficiency;}
   Int_t GetEmbeddedEventID() {return fEmbEventID;}

   //ConvertedPair bit is set for events when photon's FirstMother is not e+/e- but pi0, but after pi0 decayed
//there is conversion of one or both of the photons and results of their conversion are registered by PHOS.
//This process is marked as tagged photons but actually the energy of photons is changed and pi0 can't be
//correctly found.
   Int_t IsConverted(void) const { if(fConverted == 1) return 1; else return 0; }
//Converted bit is set if this photon originate from e+/e- conversion on medium
   Int_t IsPi0Decay(void) const { if(fPi0Decayflag == 1) return 1; else return 0; }
//Pi0Decayflag is set if this photon originate from pi0 decay
   void Pi0Decay(Int_t flag){ fPi0Decayflag=flag; }
   void Pi0Id(Int_t id){ fPi0Id=id; }
//Id of pi0 from which this photon is decayed (to check if 2 photons originate from the same pi0 or not)

   Int_t ComparePi0Ids( AliCaloPhoton *phot) { if(AliCaloPhoton::fPi0Id!=0 && (*phot).fPi0Id !=0 && AliCaloPhoton::fPi0Id == (*phot).fPi0Id) return 1; else return 0; }

   void SetBC(Int_t bc){fBC = bc;}
   void SetCluster(AliVCluster* cluster) { fCluster = cluster; }
   void SetConverted(Int_t flag){ fConverted=flag; }
   void SetConvertedPartner(Int_t flag){ fConvertedPartner=flag; }
   void SetCPVBit(Bool_t cpv){fCpv = cpv; }
   void SetCPV2Bit(Bool_t cpv){fCpv2 = cpv; }
   void SetDispBit(Bool_t chi2){fDisp = chi2 ;} 
   void SetDisp2Bit(Bool_t chi2){fDisp2 = chi2 ;} 
   void SetDistToBad(Int_t dist){fBadDist=dist;} 
   void SetDistToBadfp(Double_t dist){fBadDistfp=dist;} 
   void SetEMCx(Double_t x){fX = x ;} 
   void SetEMCy(Double_t y){fY = y ;} 
   void SetEMCz(Double_t z){fZ = z ;} 
   void SetFiducialArea(Int_t a){fFiducialArea=a ;}
   void SetIsolationTag(Int_t tag){fIsolationTag=tag ;}
   void SetIsolated(Bool_t bit){fIsIsolated=bit;}
   void SetLambdas(Double_t l1,Double_t l2){fLambda0=l1; fLambda1=l2;}
   void SetModule(Int_t mod){fModule = mod ;} 
   void SetMomV2(TLorentzVector * p){fMomV2=(*p);}
   void SetNCells(Int_t n){fNCells=n;}
   void SetPartnerPt(Double_t pt){fPartnerPt=pt;}
   void SetPCAPID(Bool_t pca){fPCA = pca;}
   void SetPhoton(Int_t flag){ fIsPhoton=flag; }
   void SetPrimary(Int_t label){fPrimary=label;}
   void SetPrimaryAtVertex(Int_t label){fPrimaryAtVertex=label;}
   void SetTagged(Bool_t bit){fIsTagged=bit;}
   void SetTagged(Bool_t bit,Int_t i,Int_t k){fIsTagged_reg[i][k]=bit;}
   void SetTagInfo(Int_t bits){fTagInfo=bits;}
   void SetTime(Double_t t) {fTime=t ;}
   void SetTOFBit(Bool_t tof){fTof = tof ;} 
   void SetTrig(Bool_t trig){fTrig=trig;}
   void SetUnfolded(Bool_t wasNotUnfolded){fUnfolded=wasNotUnfolded;} 
   void SetWeight(Double_t w){fWeight=w;}

   void SetNsigmaCPV(Double_t nsigma)      {fNsigmaCPV      = nsigma;}
   void SetNsigmaFullDisp(Double_t nsigma) {fNsigmaFullDisp = nsigma;}
   void SetNsigmaCoreDisp(Double_t nsigma) {fNsigmaCoreDisp = nsigma;}

   void SetTOFCutEfficiency(Double_t eff) {fTOFCutEfficiency = eff;}
   void SetEmbeddedEventID(Int_t id) {fEmbEventID = id;}

   AliVCluster* GetCluster() { return fCluster; }

private:
  AliCaloPhoton(const AliCaloPhoton&); // not implemented
  AliCaloPhoton& operator=(const AliCaloPhoton&);
  
  TLorentzVector fMomV2 ; //Alternative momentum
  Bool_t    fDisp ;   //Dispersion bit
  Bool_t    fDisp2 ;  //Strict Dispersion bit
  Bool_t    fTof ;    //TOF bit
  Bool_t    fCpv ;    //Charged bit
  Bool_t    fCpv2 ;   //Strict Charged bit
  Bool_t    fPCA ;    //Principal Component Analysis bit
  Bool_t    fTrig ;      //If this photon fired trigger
  Bool_t    fIsTagged;   //If it is tagged 
  Bool_t    fIsTagged_reg[10][20];   //If it is tagged 
  Bool_t    fIsIsolated ; //it is isolated
  Bool_t    fIsPhoton; //If it is really photon or not
  Bool_t    fUnfolded;  //True if was not unfolded
  Int_t     fModule ;   //Module number
  Int_t     fBC ;       //Bunch crossing number (BC=0 is main-main collision)
  Int_t     fBadDist ;  //Distance to bad module in module units
  Double_t   fBadDistfp ;  //Distance to the closest bad cell in cell units with floating point.
  Int_t     fNCells ;   //Number of cells in cluster
  Int_t     fFiducialArea ; //class of fiducial areas
  Int_t     fPi0Decayflag; //if this photon is from pi0 decay (from simulation)
  Int_t     fPi0Id;
  Int_t     fConverted; //If this photon originated from convertion on material (i.e. its primary is electron)
  Int_t	    fConvertedPartner;
  Int_t     fIsolationTag ;
  Int_t     fTagInfo ;
  Int_t     fPrimary;   //Primary entered PHOS
  Int_t     fPrimaryAtVertex;   //Primary at vertex
  Double_t  fX ;        //Cluster coordinates in ALICE ref system 
  Double_t  fY ;        //Cluster coordinates in ALICE ref system
  Double_t  fZ ;        //Cluster coordinates in ALICE ref system
  Double_t  fLambda0 ;  //Short and 
  Double_t  fLambda1 ;  //Long dispersion axis
  Double_t  fTime ;     //time of the cluster
  Double_t  fPartnerPt;
  Double_t  fWeight ;   //Weight of parent particle
  Double_t  fNsigmaCPV ;      //distance to a matched track in unit of sigma
  Double_t  fNsigmaFullDisp ; //shower dispersion of a full cluster in unit of sigma
  Double_t  fNsigmaCoreDisp ; //shower dispersion at a core of a cluster in unit of sigma
  Double_t  fTOFCutEfficiency; //TOF cut effciency at a cluster level
  Int_t fEmbEventID;//event ID used in only embedding analysis for re-shuffle
  AliVCluster* fCluster; //! Originating Cluster the Photon Candidate is based on

  ClassDef(AliCaloPhoton,11);

};

#endif // #ifdef ALICALOPHOTON_H

  
