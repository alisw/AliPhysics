#ifndef ALICALOPHOTON_H
#define ALICALOPHOTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */
 
//_________________________________________________________________________
// Class to fill two-photon invariant mass hisograms
// to be used to extract pi0 raw yield.
//
//-- Author: Dmitri Peressounko (RRC "KI")
// This class contains all (minimal) necessary information about photon to 
// calculate invarint mass distr for pi0
// and for tagging and isolation analysis


#include "TLorentzVector.h"

class AliCaloPhoton :public TLorentzVector{
  
 public:
  
  AliCaloPhoton() ;
  AliCaloPhoton(Double_t px,Double_t py,Double_t pz,Double_t E) ; 
  ~AliCaloPhoton(){} 

   Bool_t   IsDispOK(void)const {return fDisp;}
   Bool_t   IsTOFOK(void)const {return fTof;}
   Bool_t   IsCPVOK(void)const {return fCpv;}
   Bool_t   IsIsolated(void)const{return fIsIsolated ;}
   Bool_t   IsTagged(void) const{return fIsTagged ;}
   Double_t EMCx(void)const {return fZ;}
   Double_t EMCy(void)const {return fZ;}
   Double_t EMCz(void)const {return fZ;}
   Int_t    Module(void)const{return fModule;}
   Bool_t   IsPIDOK(const Int_t ipid) const ;
   Int_t    DistToBad()const  {return fBadDist ;} 
   

   void SetDispBit(Bool_t chi2){fDisp = chi2 ;} 
   void SetTOFBit(Bool_t tof){fTof = tof ;} 
   void SetCPVBit(Bool_t cpv){fCpv = cpv; }
   void SetPCAPID(Bool_t pca){fPCA = pca;}
   void SetTrig(Bool_t trig){fTrig=trig;}
   void SetEMCx(Double_t x){fX = x ;} 
   void SetEMCy(Double_t y){fY = y ;} 
   void SetEMCz(Double_t z){fZ = z ;} 
   void SetModule(Int_t mod){fModule = mod ;} 
   void SetDistToBad(Int_t dist){fBadDist=dist;} 
   void SetTagged(Bool_t bit){fIsTagged=bit;}
   void SetIsolated(Bool_t bit){fIsIsolated=bit;}
   
private:
  Bool_t    fDisp ;   //Dispersion bit
  Bool_t    fTof ;    //TOF bit
  Bool_t    fCpv ;    //Charged bit
  Bool_t    fPCA ;    //Principal Component Analysis bit
  Bool_t    fTrig ;      //If this photon fired trigger
  Bool_t    fIsTagged;   //If it is tagged 
  Bool_t    fIsIsolated ; //it is isolated
  Double_t  fX ;        //Cluster coordinates in ALICE ref system 
  Double_t  fY ;        //Cluster coordinates in ALICE ref system
  Double_t  fZ ;        //Cluster coordinates in ALICE ref system
  Int_t     fModule ;   //Module number
  Int_t     fBadDist ;  //Distance to bad module in module units

  ClassDef(AliCaloPhoton,1)

};

#endif // #ifdef ALICALOPHOTON_H

  
