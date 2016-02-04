#ifndef ALIADV1_H
#define ALIADV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : AD       //
//////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                  AD (ALICE Diffractive)  Detector                     //
//                                                                       //
//  This class contains the base procedures for the AD  detector         //
//  New geometry of 2014                                                 //
//  All comments should be sent to :                                     //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliAD.h"
#include "TGeoCompositeShape.h"
class AliADv1 : public AliAD {
public:
   
                        AliADv1();
                        AliADv1(const char *name, const char *title);
  virtual void   AddAlignableVolumes() const;
  virtual              ~AliADv1();

   
  virtual TString  Version() { return TString("v1"); }
  virtual   Int_t  IsVersion() const { return 1; }
  virtual    void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  // virtual    void  AddDigits(Int_t* track, Int_t module, Float_t time);
  virtual    void  MakeBranch(Option_t *option);
  virtual    void  CreateGeometry();
  virtual    void  Init();
  virtual    void  StepManager();
  virtual    void  DisableTunnelStruct() { fADCstruct = kFALSE; }
  virtual    void  KeepHistory() { fKeepHistory = kTRUE; }

  enum ADCPosition_t { kADCInTunnel, kADCInCavern, kADCInBoth};

protected:

  // functions for ADA and ADC
  void ReadADCFromEnv(void);
  TGeoCompositeShape * MakeShapeADCpadH(const Double_t W, const Double_t H, const Double_t dz);
  virtual    void  CreateAD();
private:
  // Position of ADC: In the Tunnel, In the Cavern, or in Both
  Bool_t      fADCstruct;
  ADCPosition_t fADCPosition;
  //! ADC Geometrical & Optical parameters :
  
  Double_t    fADCLightYield;       //! Lightyield in NE102
  Double_t    fADCPhotoCathodeEfficiency;

  //! ADA Geometrical & Optical parameters :
 
  Double_t    fADALightYield;       //! Lightyield in NE102
  Double_t    fADAPhotoCathodeEfficiency; 
  
  Bool_t      fKeepHistory; 


  AliADv1(const AliAD&); 
  AliADv1& operator = (const AliADv1&); 
  
  ClassDef(AliADv1, 1)  //!Class for the AD detector
   
};


#endif
