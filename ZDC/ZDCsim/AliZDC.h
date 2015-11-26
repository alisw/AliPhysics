#ifndef ALIZDC_H
#define ALIZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and classes for set ZDC           //
////////////////////////////////////////////////

#include <TSystem.h>

#include "AliDetector.h"
#include "AliZDCTrigger.h"
#include "AliZDCChMap.h"

class AliZDCPedestals;
class AliZDCEnCalib;
class AliZDCTowCalib;
 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual       ~AliZDC();
  AliZDC(const AliZDC&);
  //
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  virtual Int_t IsVersion() const =0;
  virtual Float_t ZMin() const;	// Minimum overall dimension of the ZDC
  virtual Float_t ZMax() const;	// Maximum overall dimension of the ZDC
  virtual void  SetTreeAddress();
  virtual void  MakeBranch(Option_t* opt);
  virtual void  Hits2SDigits();
  virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
  virtual void  Digits2Raw();
  virtual Bool_t Raw2SDigits(AliRawReader* rawReader);
  Int_t   Pedestal(Int_t Detector, Int_t Quadrant, Int_t Res) const;
  Int_t   ADCch2Phe(Int_t Detector, Int_t Quadrant, Int_t ADCVal, Int_t Res) const;
  virtual void  StepManager() {}
    
  // Switching off the shower development in ZDCs
  void  NoShower(){fNoShower=1;}
  void  Shower()  {fNoShower=0;}
  
  virtual void SetVCollSideCAperture(Float_t /*aperture*/) {}
  virtual void SetVCollSideCApertureNeg(Float_t /*aperture*/) {}
  virtual void SetVCollSideCCentre(Float_t /*centre*/) {}
  
  virtual void SetVCollSideAAperture(Float_t /*aperture*/) {}
  virtual void SetVCollSideAApertureNeg(Float_t /*aperture*/) {}
  virtual void SetVCollSideACentre(Float_t /*centre*/) {}
  
  virtual void SetTCDDAperturePos(Float_t /*aperture*/) {}
  virtual void SetTCDDApertureNeg(Float_t /*aperture*/) {}
  
  virtual void SetTDIAperturePos(Float_t /*aperture*/) {}
  virtual void SetTDIApertureNeg(Float_t /*aperture*/) {}
  virtual void SetTDIConfiguration(Int_t /*configuration*/) {}
  
  virtual void SetLumiLength(Float_t /*length*/) {}
  
  virtual void SetYZNC(Float_t /*yZNC*/) {}
  virtual void SetYZNA(Float_t /*yZNA*/) {}
  virtual void SetYZPC(Float_t /*yZPC*/) {}
  virtual void SetYZPA(Float_t /*yZPA*/) {}
  
  virtual void SetSwitchOnTrackreferences() {}
  
  //Calibration methods 
  void    SetZDCCalibFName(const char *name);
  char*   GetZDCCalibFName() const {return (char*)fZDCCalibFName.Data();}
  AliZDCPedestals* GetPedCalib()   const  {return fPedCalib;}
  AliZDCEnCalib*   GetECalibData() const  {return fEnCalibData;}
  
  // Map from OCDB
  AliZDCChMap*     GetChMap() const;

  // Trigger
  virtual AliTriggerDetector* CreateTriggerDetector() const
  	{return new AliZDCTrigger();}

  void  SetSpectatorParam(int imod=1) {fSpectatorParam=imod;}
  void  SetSpectatorsTrack() {fSpectatorTracked=0;}
  Int_t AreSpectatorsTracked() const {return fSpectatorTracked;}
  void  SetBeamEnergy(Float_t beamEnergy) {fBeamEnergy = beamEnergy;}
  void  SetpAsystem() {fIspASystem = kTRUE;}
  void  SetRELDISGenerator() {fIsRELDISgen = kTRUE;}
  
  void  SetOnlyZEM() 		{fOnlyZEM=kTRUE;}
  void  SetMotherFinding()	{fFindMother=kTRUE;}

private:

  AliZDC& operator = (const AliZDC&);

protected:

  Int_t        fNoShower;		// Flag to switch off the shower	

  //Calibration data member 
  AliZDCPedestals* fPedCalib;		//! Pedestal data for ZDC
  AliZDCEnCalib*   fEnCalibData;	//! Energy data for ZDC
  AliZDCTowCalib*  fTowCalibData;	//! Equalization data for ZDC

  TString          fZDCCalibFName; 	// Name of the ZDC calibration data
 
  Int_t   fSpectatorTracked; // Are spectator tracked by generator? 0=NO
  Float_t fBeamEnergy;	     // beam energy from generator (AliGenZDC + RELDIS)
  Bool_t  fIspASystem;       // Configuring pA collisions (MC only)
  Bool_t  fIsRELDISgen;	     // Is RELDIS used as generator
  
  Bool_t  fOnlyZEM;	     // build only ZEM (no had. calorimeters!)
  Bool_t  fFindMother;	     // look for particle mothers in the stack in StepManager
  Int_t   fSpectatorParam;   // kinematic model fro spectators (1=from AliGenZDC(DEFAULT), 2=from HIJING)
  
  ClassDef(AliZDC,15)  	// Zero Degree Calorimeter base class
};
 
// Calibration
//_____________________________________________________________________________
inline void AliZDC::SetZDCCalibFName(const char *name)  
{ 
  fZDCCalibFName = name;        
  gSystem->ExpandPathName(fZDCCalibFName);
}


#endif
