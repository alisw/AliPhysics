#ifndef ALIPHOSBEAMTESTEVENT_H
#define ALIPHOSBEAMTESTEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class for PHOS Beam Test event header. Contanes all information 
//  about triggers, etc.    
//                  
//*-- Author: Maxim Volkov (RRC KI) & Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TObject.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSBeamTestEvent:public TObject {

public:
  AliPHOSBeamTestEvent() ;          // ctor

  virtual ~AliPHOSBeamTestEvent() ; // dtor

  Float_t   GetBeamEnergy(void) const {return fBeamEnergy ;}
  UInt_t *  GetUserVector(void)       {return fUserVector ;}
  UInt_t *  GetHeader(void)           {return fHeader ;}
  UShort_t  GetPattern(void) const    {return fPattern ;}
  UShort_t *GetScanning(void)         {return fScanning ;}
  UShort_t *GetCharge(void)           {return fCharge ;}
  UInt_t *  GetScaler(void)           {return fScaler ;}
  UShort_t *GetTDC(void)              {return fTDC2228 ;}

  void SetBeamEnergy(Float_t energy ){fBeamEnergy = energy ;}
  void SetUserVector(UInt_t * uv){
	  for(Int_t i=0;i<16;i++)fUserVector[i]=uv[i];}
  void SetHeader(UInt_t * h){
	  for(Int_t i=0;i<12;i++)fHeader[i]=h[i];}
  void SetPattern(UShort_t pat){fPattern=pat ;}
  void SetScanning(UShort_t * scan){
	 for(Int_t i=0;i<32;i++) fScanning[i]=scan[i] ;}
  void SetCharge(UShort_t *charg){
	  for(Int_t i=0;i<12;i++) fCharge[i]=charg[i] ;}
  void SetScaler(UInt_t * sc){
	 for(Int_t i=0;i<12;i++) fScaler[i]=sc[i] ;}
  void SetTDC(UShort_t * tdc) {
	  for(Int_t i=0;i<12;i++) fTDC2228[i]=tdc[i] ;}
private:
  Float_t  fBeamEnergy ;         //Beam energy 
  UInt_t   fUserVector[16] ;     //ZEBRA Event user vector
  UInt_t   fHeader[12] ;         //ZEBRA event header
  UInt_t   fScaler[12] ;         //Scalers, 1 module X 12 (4 byte) ch.
  UShort_t fPattern ;            //Trigger bit register
  UShort_t fScanning[32] ;       //Scanning ADCs,4 modulesX8=32 channels
  UShort_t fCharge[12] ;         //Charge ADCs, 1 module X 12 = 12 ch.
  UShort_t fTDC2228[32] ;        //LeCroy TDC 2228A, 4 module X 8 =32 ch

  ClassDef(AliPHOSBeamTestEvent,1)  // description 

};

#endif // ALIPHOSBEAMTESTEVENT_H
