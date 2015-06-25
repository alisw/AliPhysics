#ifndef ALIZDCV5_H
#define ALIZDCV5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set: ZDC     //
////////////////////////////////////////////////

#include "AliZDC.h"

//____________________________________________________________________________ 
class AliZDCv5 : public AliZDC {

public:
  AliZDCv5();
  AliZDCv5(const char *name, const char *title);
  virtual  ~AliZDCv5() {}
  virtual void  CreateGeometry();
  virtual void  CreateBeamLine();
  virtual void  CreateZDC();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  AddAlignableVolumes() const;
  virtual void  Init();
  virtual void  InitTables();
  virtual void  StepManager();
  
  void SetVCollSideCAperture(Float_t aperture)
  	{if(aperture<3.5) fVCollSideCAperture = aperture; 
	 else printf("\n\n AliZDCv5: WARNING! SideC TCTVB aperture set to max. value: 3.5 cm\n\n");}
  void SetVCollSideCApertureNeg(Float_t aperture)
  	{if(aperture<3.5) fVCollSideCApertureNeg = aperture; 
	 else printf("\n\n AliZDCv5: WARNING! SideC TCTVB aperture set to max. value: -3.5 cm\n\n");}
  void SetVCollSideCCentre(Float_t centre) {fVCollSideCCentreY = centre;}
  
  virtual void SetTCDDAperturePos(Float_t aperture) 
  	{if(aperture<=2.2) fTCDDAperturePos = aperture;
	 else printf("\n\n AliZDCv5: WARNING! TCDD pos. aperture set to max. value: 2.0 cm\n\n");}
  virtual void SetTCDDApertureNeg(Float_t aperture) 
    	{if(aperture<=2.4) fTCDDApertureNeg = aperture;
	 else printf("\n\n AliZDCv5: WARNING! TCDD neg. aperture set to max. value: -2.0 cm\n\n");}

  
  virtual void SetTDIAperturePos(Float_t aperture) 
  	{if(aperture<=6.) fTDIAperturePos = aperture;
	 else printf("\n\n AliZDCv5: WARNING! TDI pos. aperture set to max. value: 5.5 cm\n\n");}
  virtual void SetTDIApertureNeg(Float_t aperture) 
  	{if(aperture<=6.) fTDIApertureNeg = aperture;
	 else printf("\n\n AliZDCv5: WARNING! TDI neg. aperture set to max. value: -5.5 cm\n\n");}
  virtual void SetTDIConfiguration(Int_t configuration) 
  	{if(fTDIConfiguration>=0 && fTDIConfiguration<=2) fTDIConfiguration=configuration;
	 else printf("\n\n AliZDCv5: WARNING! TDI invalid configuration -> setting to 2\n\n");}

  void SetLumiLength(Float_t length) {fLumiLength = length;}
  
  void SetYZNC(Float_t yZNC) {fPosZNC[1] = yZNC;}
  void SetYZNA(Float_t yZNA) {fPosZNA[1] = yZNA;}
  
  void SetYZPC(Float_t yZPC) {fPosZPC[1] = yZPC;}
  void SetYZPA(Float_t yZPA) {fPosZPA[1] = yZPA;}
  
  void SetSwitchOnTrackreferences() {fSwitchOnTrackRef = kTRUE;}
 
protected:

  // Sensitive media
  Int_t   fMedSensF1;         // Sensitive medium F1
  Int_t   fMedSensF2;         // Sensitive medium F2
  Int_t   fMedSensZP;         // Sensitive medium for ZP
  Int_t   fMedSensZN;         // Sensitive medium for ZN
  Int_t   fMedSensZEM;        // Sensitive medium for EM ZDC
  Int_t   fMedSensGR;         // Other sensitive medium
  Int_t   fMedSensPI;         // Beam pipe and magnet coils
  Int_t   fMedSensTDI;        // Cu materials along beam pipe
  Int_t   fMedSensVColl;      // W jaws of vertical collimators
  Int_t   fMedSensLumi;       // luminometer medium
  
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
  Float_t fPosZNC[3];   // Position of neutron detector side C
  Float_t fPosZNA[3];   // Position of neutron detector side A  
  Float_t fPosZPC[3]; 	// Position of proton detector side C
  Float_t fPosZPA[3]; 	// Position of proton detector side A
  Float_t fFibZN[3]; 	// Fibers for neutron detector
  Float_t fFibZP[3];  	// Fibers for proton detector

  // Parameters for EM calorimeter geometry
  // NB -> parameters used in CreateZDC() and in StepManager()
  // (other parameters are defined in CreateZDC())
  Float_t fPosZEM[3]; // Position of EM detector
  Float_t fZEMLength; // ZEM length
  
  // Parameters for proton accepancy studies
  Int_t fpLostITC, fpLostD1C, fpcVCollC, fpDetectedC, fnDetectedC; // Side C
  Int_t fpLostITA, fpLostD1A, fpLostTDI, fpcVCollA, fpDetectedA, fnDetectedA; // Side A
  
  // Apertures to describe beam line elements variable apertures
  
  // Vertical collimator
  Float_t fVCollSideCAperture;    // Semi-aperture of TCTVB jaws pos. y dir.
  Float_t fVCollSideCApertureNeg; // Semi-aperture of TCTVB jaws neg. y dir (abs. value)
  Float_t fVCollSideCCentreY;     // Centre of TCTVB jaw apertures
  
  // TCDD
  Float_t fTCDDAperturePos;       // TCDD semi-aperture pos. y dir.
  Float_t fTCDDApertureNeg;       // TCDD semi-aperture neg. y dir. (abs. value)
  
  // TDI
  Float_t fTDIAperturePos;	  // TDI jaw semi-aperture pos. y dir.
  Float_t fTDIApertureNeg;	  // TDI jaw semi-aperture  neg. y dir. (abs. value)
  //
  Int_t fTDIConfiguration;	  // choose TDI design
  
  Float_t fLumiLength;  	  // Luminometer length
  Bool_t  fSwitchOnTrackRef;      // to switch on/off storing of track references
  
  ClassDef(AliZDCv5, 2)  // Zero Degree Calorimeter version 1
}; 
 
#endif
