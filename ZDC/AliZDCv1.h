#ifndef ALIZDCV1_H
#define ALIZDCV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:ZDC      //
////////////////////////////////////////////////

#include "AliZDC.h"

//____________________________________________________________________________ 
class AliZDCv1 : public AliZDC {

public:
  AliZDCv1();
  AliZDCv1(const char *name, const char *title);
  virtual      ~AliZDCv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateBeamLine();
  virtual void  CreateZDC();
  virtual void  CreateMaterials();
  Int_t         Digitize(Int_t Det, Int_t Quad, Int_t Light);
  virtual void  FinishEvent();
  virtual void  MakeBranch(Option_t* opt);
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
  virtual void  Init();
  virtual void  InitTables();
  virtual void  StepManager();
  
  // Switching off the shower development in ZDCs
  void  NoShower(){fNoShower=1;}
  void  Shower()  {fNoShower=0;}
  
  // Digitization parameters setters and getters

  // ADC pedestal mean value
  void SetPedMean(Int_t Det, Int_t PMDet, Int_t PedMean)
       {fPedMean[Det][PMDet] = PedMean;}
  Float_t GetPedMean(Int_t Det, Int_t PMDet)
       {return fPedMean[Det][PMDet];} 
  // ADC pedestal width
  void SetPedSigma(Int_t Det, Int_t PMDet, Int_t PedSigma)
       {fPedSigma[Det][PMDet] = PedSigma;}
  Float_t GetPedSigma(Int_t Det, Int_t PMDet)
       {return fPedSigma[Det][PMDet];}
  // PM gain
  void SetPMGain(Int_t Det, Int_t PMDet, Int_t PMGain)
       {fPMGain[Det][PMDet] = PMGain;}
  Float_t GetPMGain(Int_t Det, Int_t PMDet)
       {return fPMGain[Det][PMDet];}
  // Conversion factor from charge to ADC channels
  //   F = 1.6E-19 / Resolution [Coulomb/ch]
  void SetADCRes(Int_t ADCRes) {fADCRes =  ADCRes;}
  Float_t GetADCRes() {return fADCRes;}
 
protected:
  // Sensitive media
  Int_t   fMedSensF1;         // Sensitive medium F1
  Int_t   fMedSensF2;         // Sensitive medium F2
  Int_t   fMedSensZP;         // Sensitive medium for ZP
  Int_t   fMedSensZN;         // Sensitive medium for ZN
  Int_t   fMedSensZEM;        // Sensitive medium for EM ZDC
  Int_t   fMedSensGR;         // Other sensitive medium
  Int_t   fMedSensPI;         // Beam pipe and magnet coils
  
  // Parameters for light tables
  Int_t   fNalfan;            // Number of Alfa neutrons
  Int_t   fNalfap;            // Number of Alfa protons
  Int_t   fNben;              // Number of beta neutrons
  Int_t   fNbep;              // Number of beta protons
  Float_t fTablen[4][90][18]; // Neutrons light table
  Float_t fTablep[4][90][28]; // Protons light table
  
  // Parameters for conversion of light yield in ADC channels
  Float_t fPedMean[3][5];     // ADC pedestal mean value
  Float_t fPedSigma[3][5];    // ADC pedestal width
  Float_t fPMGain[3][5];      // PM gain
  Float_t fADCRes;            // ADC conversion factor

public:
  //Flag for fast simulation (no shower)
  Int_t   fNoShower;
  
   ClassDef(AliZDCv1,1)  // Zero Degree Calorimeter version 1
}; 
 
#endif
