#ifndef ALIPHOSRAWDIGIPRODUCER_H
#define ALIPHOSRAWDIGIPRODUCER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

//This class produces PHOS digits of one event
//using AliPHOSRawDecoder. See cxx source for use case.

class AliPHOSRawDecoder;
class AliPHOSCalibData ;
class AliPHOSDigit ;
class AliPHOSRecoParam ;
class AliPHOSGeometry ;
class AliPHOSPulseGenerator;
#include "TObject.h"

class AliPHOSRawDigiProducer: public TObject {

public:

  AliPHOSRawDigiProducer() ;
  AliPHOSRawDigiProducer(const AliPHOSRecoParam* recoParam) ;
  AliPHOSRawDigiProducer(const AliPHOSRawDigiProducer &dp);
  AliPHOSRawDigiProducer& operator= (const AliPHOSRawDigiProducer &dp);
 
  virtual ~AliPHOSRawDigiProducer(); 

  void MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder);

protected:

  void GetCalibrationParameters() ; //Extract calibration parameters from DB
  void CleanDigits(TClonesArray* digits) ; //remove digits below threshold and bad ones
  
  Bool_t IsInEMC(AliPHOSDigit* digit) const ; //tests if digit belongs to EMC
  Bool_t IsInCPV(AliPHOSDigit* digit) const ;

  Double_t CalibrateE(Double_t amp, Int_t* relId, Bool_t isLowGain) ; //calibrate energy 
  Double_t CalibrateT(Double_t amp, Int_t* relId, Bool_t isLowGain) ; //calibrate time

private:
  Float_t fEmcMinE ;                 // minimum energy of digit to be included into cluster
  Float_t fCpvMinE ;                 // minimum energy of digit to be included into cluster
  Float_t fSampleQualityCut;         // Cut on sample shapes: 0: no samples; 1: default parameterization; 999: accept even obviously bad
  Int_t   fGlobalAltroOffset;        // Global ALTRO offset used in ZS runs
  Int_t   fGlobalAltroThreshold;     // Global ALTRO threshold used in ZS runs
  Int_t fEmcCrystals ;               //  number of EMC crystals
  AliPHOSGeometry * fGeom ;          //! PHOS geometry
  static AliPHOSCalibData * fgCalibData ;   //! Calibration database if avalable
  AliPHOSPulseGenerator * fPulseGenerator ; //! Class with pulse shape parameters

  ClassDef(AliPHOSRawDigiProducer,4)
};

#endif
