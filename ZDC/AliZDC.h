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
  virtual void SetVCollSideCCentre(Float_t /*centre*/) {}
  
  virtual void SetVCollSideAAperture(Float_t /*aperture*/) {}
  virtual void SetVCollSideACentre(Float_t /*centre*/) {}
  
  virtual void SetLumiLength(Float_t /*length*/) {}
  
  virtual void SetYZNC(Float_t /*yZNC*/) {}
  virtual void SetYZNA(Float_t /*yZNA*/) {}
  virtual void SetYZPC(Float_t /*yZPC*/) {}
  virtual void SetYZPA(Float_t /*yZPA*/) {}

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

  
  void  SetSpectatorsTrack() {fSpectatorTracked=0;}
  Int_t SpectatorsTracked() const {return fSpectatorTracked;}

private:

  AliZDC& operator = (const AliZDC&);

protected:

  Int_t        fNoShower;		// Flag to switch off the shower	

  //Calibration data member 
  AliZDCPedestals* fPedCalib;		//! Pedestal data for ZDC
  AliZDCEnCalib*   fEnCalibData;	//! Energy data for ZDC
  AliZDCTowCalib*  fTowCalibData;	//! Equalization data for ZDC

  TString          fZDCCalibFName; 	// Name of the ZDC calibration data
 
  Int_t fSpectatorTracked; // Are spectator tracked by generator?
  
  ClassDef(AliZDC,10)  	// Zero Degree Calorimeter base class
};
 
// Calibration
//_____________________________________________________________________________
inline void AliZDC::SetZDCCalibFName(const char *name)  
{ 
  fZDCCalibFName = name;        
  gSystem->ExpandPathName(fZDCCalibFName);
}


#endif
