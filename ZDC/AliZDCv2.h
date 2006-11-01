#ifndef ALIZDCV2_H
#define ALIZDCV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:ZDC      //
////////////////////////////////////////////////

#include "AliZDC.h"

//____________________________________________________________________________ 
class AliZDCv2 : public AliZDC {

public:
  AliZDCv2();
  AliZDCv2(const char *name, const char *title);
  virtual  ~AliZDCv2() {}
  virtual void  CreateGeometry();
  virtual void  CreateBeamLine();
  virtual void  CreateZDC();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule() const;
  virtual void  AddAlignableVolumes() const;
  virtual void  Init();
  virtual void  InitTables();
  virtual void  StepManager();
  
 
protected:

  // Sensitive media
  Int_t   fMedSensF1;         // Sensitive medium F1
  Int_t   fMedSensF2;         // Sensitive medium F2
  Int_t   fMedSensZP;         // Sensitive medium for ZP
  Int_t   fMedSensZN;         // Sensitive medium for ZN
  Int_t   fMedSensZEM;        // Sensitive medium for EM ZDC
  Int_t   fMedSensGR;         // Other sensitive medium
  Int_t   fMedSensPI;         // Beam pipe and magnet coils
  Int_t   fMedSensTDI;        // TDI Cu shielding 
  
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
  Float_t fDimZN[3];	// Dimensions of proton detector
  Float_t fDimZP[3];	// Dimensions of proton detector
  Float_t fPosZN[3];   	// Position of neutron detector
  Float_t fPosZP[3]; 	// Position of proton detector
  Float_t fFibZN[3]; 	// Fibers for neutron detector
  Float_t fFibZP[3];  	// Fibers for proton detector

  // Parameters for EM calorimeter geometry
  // NB -> parameters used in CreateZDC() and in StepManager()
  // (other parameters are defined in CreateZDC())
  Float_t fPosZEM[3]; // Position of EM detector
  Float_t fZEMLength; // ZEM length
  
  // Parameters for tracking studies
  Int_t fpLostIT, fpLostD1, fpLostTDI, fpDetected, fnDetected; // For spectators acceptance
  
   ClassDef(AliZDCv2,4)  // Zero Degree Calorimeter version 1
}; 
 
#endif
