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
  virtual void  SDigits2Digits();
  virtual void  MakeBranch(Option_t* opt, char *file=0);
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
  virtual void  Init();
  virtual void  InitTables();
  virtual void  Hits2Digits(Int_t ntracks = 0);
  virtual void  StepManager();
  
//  // Digitization parameters setters and getters
//  // PM gain
//  void SetPMGain(Int_t Det, Int_t PMDet, Int_t PMGain)
//       {fPMGain[Det][PMDet] = PMGain;}
//  Float_t GetPMGain(Int_t Det, Int_t PMDet)
//       {return fPMGain[Det][PMDet];}
//  // Conversion factor from charge to ADC channels
//  //   F = 1.6E-19 / Resolution [Coulomb/ch]
//  void SetADCRes(Int_t ADCRes) {fADCRes =  ADCRes;}
//  Float_t GetADCRes() {return fADCRes;}
 
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
  Int_t   fNalfan;	      // Number of Alfa (neutrons)
  Int_t   fNalfap;	      // Number of Alfa (protons)
  Int_t   fNben;	      // Number of beta (neutrons)
  Int_t   fNbep;	      // Number of beta (protons)
  Float_t fTablen[4][90][18]; // Neutrons light table
  Float_t fTablep[4][90][28]; // Protons light table

  // Parameters for hadronic calorimeters geometry
  // NB -> parameters used in CreateZDC() and in StepManager()
  // (other parameters are defined in CreateZDC())
  Float_t fDimZP[3];	// Dimensions of proton detector
  Float_t fPosZN[3];   	// Position of neutron detector
  Float_t fPosZP[3]; 	// Position of proton detector
  Float_t fFibZN[3]; 	// Fibers for neutron detector
  Float_t fFibZP[3];  	// Fibers for proton detector

  // Parameters for EM calorimeter geometry
  // NB -> parameters used in CreateZDC() and in StepManager()
  // (other parameters are defined in CreateZDC())
  Float_t fPosZEM[3]; // Position of EM detector
  
//  // Parameters for conversion of light yield in ADC channels
//  Float_t fPMGain[3][5];      // PM gain
//  Float_t fADCRes;            // ADC conversion factor
  
   ClassDef(AliZDCv1,1)  // Zero Degree Calorimeter version 1
}; 
 
#endif
