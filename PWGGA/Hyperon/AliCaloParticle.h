#ifndef ALICALOPARTICLE_H
#define ALICALOPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */
 
//_________________________________________________________________________
// Class to fill two-photon invariant mass hisograms
// to be used to extract pi0 raw yield.
//
//-- Author: Dmitri Peressounko (RRC "KI"), Alexander Borissov (Wayne State University, USA)
// This class contains all (minimal) necessary information about photon to 
// calculate invarint mass distr for pi0
// and for tagging and isolation analysis


#include "TLorentzVector.h"

class AliCaloParticle :public TLorentzVector{
  
 public:
  
  AliCaloParticle() ;
  AliCaloParticle(TLorentzVector p) ; 
  AliCaloParticle(Double_t px,Double_t py,Double_t pz,Double_t E) ; 
  ~AliCaloParticle(){} 

   const TLorentzVector * GetMomV2()const{return &fMomV2;}
   
   Double_t EMCx(void)const {return fZ;}
   Double_t EMCy(void)const {return fZ;}
   Double_t EMCz(void)const {return fZ;}
   Int_t    Module(void)const{return fModule;}
   Int_t    DistToBad()const  {return fBadDist ;}
   Int_t    GetNCells()const { return fNCells ;}

   Bool_t   IsDispOK(void)const {return fDisp;}
   Bool_t   IsTOFOK(void)const {return fTof;}
   Bool_t   IsCPVOK(void)const {return fCpv;}
   Bool_t   IsCPV2OK(void)const {return fCpv2;}
   Int_t    IsIsolated(void)const{return fIsIsolated ;}
   Bool_t   IsTagged(void) const{return fIsTagged ;} //check if this photon is tagged
   Bool_t   IsTagged(Int_t i,Int_t k) const{return fIsTagged_reg[i][k] ;} //check if this photon is tagged
   Bool_t   IsPIDOK(const Int_t ipid) const ;
   Bool_t   IsPhoton()const {return fIsPhoton ;} //check if this particle is indeed photon (this bit is set with MC stack info
   Int_t    IsConvertedPartner(){ if(fConvertedPartner == 1) return 1; else return 0; }
//ConvertedPair bit is set for events when photon's FirstMother is not e+/e- but pi0, but after pi0 decayed
//there is conversion of one or both of the photons and results of their conversion are registered by PHOS.
//This process is marked as tagged photons but actually the energy of photons is changed and pi0 can't be
//correctly found.
   Int_t IsConverted(){ if(fConverted == 1) return 1; else return 0; }
//Converted bit is set if this photon originate from e+/e- conversion on medium
   Int_t IsPi0Decay(){ if(fPi0Decayflag == 1) return 1; else return 0; }
//Pi0Decayflag is set if this photon originate from pi0 decay
  
// void SetPtsum(Double_t Ptsum[10]){for(Int_t j=0; j<10; j++)  fPtsum[j]=Ptsum[j];}

void SetPtsum(Double_t *Ptsum){for(Int_t j=0; j<10; j++)  fPtsum[j]=Ptsum[j];}

   void Pi0Decay(Int_t flag){ fPi0Decayflag=flag; }
   void Pi0Id(Int_t id){ fPi0Id=id; }
//Id of pi0 from which this photon is decayed (to check if 2 photons originate from the same pi0 or not)


   void SetMomV2(TLorentzVector * p){fMomV2=(*p);}
   void SetNCells(Int_t n){fNCells=n;}
   void SetConverted(Int_t flag){ fConverted=flag; }
   Int_t ComparePi0Ids( AliCaloParticle *phot) { if(AliCaloParticle::fPi0Id!=0 && (*phot).fPi0Id !=0 && AliCaloParticle::fPi0Id == (*phot).fPi0Id) return 1; else return 0; }
   void SetConvertedPartner(Int_t flag){ fConvertedPartner=flag; }
   void SetPhoton(Int_t flag){ fIsPhoton=flag; }
   void SetDispBit(Bool_t chi2){fDisp = chi2 ;} 
   void SetTOFBit(Bool_t tof){fTof = tof ;} 
   void SetCPVBit(Bool_t cpv){fCpv = cpv; }
   void SetCPV2Bit(Bool_t cpv){fCpv2 = cpv; }
   void SetPCAPID(Bool_t pca){fPCA = pca;}
   void SetTrig(Bool_t trig){fTrig=trig;}
   void SetEMCx(Double_t x){fX = x ;} 
   void SetEMCy(Double_t y){fY = y ;} 
   void SetEMCz(Double_t z){fZ = z ;} 
   void SetModule(Int_t mod){fModule = mod ;} 
   void SetDistToBad(Int_t dist){fBadDist=dist;} 
   void SetTagged(Bool_t bit){fIsTagged=bit;}
   void SetTagged(Bool_t bit,Int_t i,Int_t k){fIsTagged_reg[i][k]=bit;}
   void SetIsolation(Int_t isolation){fIsol=isolation;}
   void SetPartnerPt(Double_t pt){fPartnerPt=pt;}
   void SetPrimary(Int_t label){fPrimary=label;}
   void SetXt(Double_t Xt){fXt=Xt;}

   void SetPpix(Double_t Ppix){fPpix=Ppix;}
   void SetPpiy(Double_t Ppiy){fPpiy=Ppiy;}
   void SetPpiz(Double_t Ppiz){fPpiz=Ppiz;}
   void SetPprx(Double_t Pprx){fPprx=Pprx;}
   void SetPpry(Double_t Ppry){fPpry=Ppry;}
   void SetPprz(Double_t Pprz){fPprz=Pprz;}

   void SetV0CosPointingAngle(Double_t value) { fV0CosPointingAngle = value; }
   void SetV0DCADaughters(Double_t value) { fV0DCADaughters = value; }
   void SetV0IPNeg(Double_t value) { fV0IPNeg = value ; }
   void SetV0IPPos(Double_t value) { fV0IPPos = value ; }
   void SetV0Radius( Double_t value) { fV0Radius = value; }
   void SetV0Chi2 ( Double_t value ) { fV0Chi2 = value; }
   void SetZConv ( Double_t value ) { fZConv = value; } 
   void SetArmPt (Double_t value ) { fArmPt = value; }
   void SetArmAlpha (Double_t value ) { fArmAlpha = value; }




   Double_t GetPtsum( Double_t* Ptsum ) ;
   Int_t GetPrimary(){return fPrimary;}
   Int_t GetIsolation(){return fIsol;}
   Double_t GetPartnerPt(){return fPartnerPt;} 
   Double_t GetXt(){return fXt;}

   Double_t GetPpix(){return fPpix;}
   Double_t GetPpiy(){return fPpiy;}
   Double_t GetPpiz(){return fPpiz;}
   Double_t GetPprx(){return fPprx;}
   Double_t GetPpry(){return fPpry;}
   Double_t GetPprz(){return fPprz;}

   Double_t GetV0CosPointingAngle() { return fV0CosPointingAngle; }
   Double_t GetV0DCADaughters() { return fV0DCADaughters; }
   Double_t GetV0IPNeg() { return fV0IPNeg; }
   Double_t GetV0IPPos() { return fV0IPPos; }
   Double_t GetV0Radius( ) { return fV0Radius; }
   Double_t GetV0Chi2 (  ) { return fV0Chi2; }
   Double_t GetZConv (  ) { return fZConv; }
   Double_t GetArmPt () { return fArmPt ; }
   Double_t GetArmAlpha () { return fArmAlpha; }


 
private:
  TLorentzVector fMomV2 ; //Alternative momentum
  Bool_t    fDisp ;   //Dispersion bit
  Bool_t    fTof ;    //TOF bit
  Bool_t    fCpv ;    //Charged bit
  Bool_t    fCpv2 ;   //Stronger charged bit
  Bool_t    fPCA ;    //Principal Component Analysis bit
  Bool_t    fTrig ;      //If this photon fired trigger
  Bool_t    fIsTagged;   //If it is tagged 
  Bool_t    fIsTagged_reg[10][20];   //If it is tagged 
  Bool_t    fIsPhoton; //If it is really photon or not
  Double_t  fX ;        //Cluster coordinates in ALICE ref system 
  Double_t  fY ;        //Cluster coordinates in ALICE ref system
  Double_t  fZ ;        //Cluster coordinates in ALICE ref system
  Int_t     fModule ;   //Module number
  Int_t     fBadDist ;  //Distance to bad module in module units
  Int_t     fNCells ;   //Number of cells in cluster
  Int_t     fPi0Decayflag; //if this photon is from pi0 decay (from simulation)
  Int_t     fPi0Id;
  Int_t     fConverted; //If this photon originated from convertion on material (i.e. its primary is electron)
  Int_t	    fConvertedPartner;
  Double_t  fPartnerPt;
  Int_t     fPrimary;   //Primary label
  Int_t     fIsIsolated ; //Isolation bits for different criteria
  Int_t     fIsol ;
  Double_t  fXt ; // scaling variable = 2Pt/sqrt(s)

  Double_t fPpix;
  Double_t fPpiy;
  Double_t fPpiz;
  Double_t fPprx;
  Double_t fPpry;
  Double_t fPprz;

  Double32_t  fPtsum[10] ; 

  
  //_______________________________________________
  Double_t fV0CosPointingAngle;
  Double_t fV0DCADaughters;
  Double_t fV0IPNeg;
  Double_t fV0IPPos;
  Double_t fV0Radius;
  Double_t fV0Chi2;
  Double_t fZConv;
  Double_t fArmPt;
  Double_t fArmAlpha;
  /*
	 Double_t fNSigmaPos;
	 Double_t fNSigmaNeg;
	 Double_t fLeastNbrCrossedRows;
	 */




  ClassDef(AliCaloParticle,1)

};

#endif // #ifdef ALICALOPHOTON_H

  
