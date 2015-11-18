#ifndef ALIPHOSRAWDIGIPRODUCER_H
#define ALIPHOSRAWDIGIPRODUCER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

//This class produces PHOS digits of one event
//using AliPHOSRawFitter. See cxx source for use case.

class AliPHOSCalibData ;
class AliPHOSDigit ;
class AliPHOSGeometry ;
class AliPHOSPulseGenerator;
class AliRawReader;
class AliCaloRawStreamV3;
class AliPHOSRawFitterv0;

#include "AliAltroMapping.h"
#include "TObject.h"

class AliPHOSRawDigiProducer: public TObject {

public:

  AliPHOSRawDigiProducer() ;
  AliPHOSRawDigiProducer(AliRawReader *rawReader, AliAltroMapping **mapping = NULL);
  AliPHOSRawDigiProducer(const AliPHOSRawDigiProducer &dp);
  AliPHOSRawDigiProducer& operator= (const AliPHOSRawDigiProducer &dp);
 
  virtual ~AliPHOSRawDigiProducer(); 

  void MakeDigits(TClonesArray *digits, AliPHOSRawFitterv0* fitter);
  void MakeDigits(TClonesArray *digits, TClonesArray *tmpDigLG, AliPHOSRawFitterv0* fitter);

  void SetEmcMinAmp(Float_t emcMin) { fEmcMinE=emcMin; }
  void SetCpvMinAmp(Float_t cpvMin) { fCpvMinE=cpvMin; }
  void SetSampleQualityCut(Float_t qcut) { fSampleQualityCut=qcut; }
  void SetSubtractL1phase(Bool_t a=kTRUE){ fSubtractL1phase=a ; }

protected:

  void GetCalibrationParameters() ; //Extract calibration parameters from DB
  void CleanDigits(TClonesArray* digits) ; //remove digits below threshold and bad ones
  
  Bool_t IsInEMC(AliPHOSDigit* digit) const ; //tests if digit belongs to EMC
  Bool_t IsInCPV(AliPHOSDigit* digit) const ;

  Double_t CalibrateE(Double_t amp, Int_t* relId, Bool_t isLowGain) ; //calibrate energy 
  Double_t CalibrateT(Double_t amp, Int_t* relId, Bool_t isLowGain) ; //calibrate time

private:
  Bool_t  fSubtractL1phase ;         // To correct time for L1phase
  Float_t fEmcMinE ;                 // minimum energy of digit (ADC)
  Float_t fCpvMinE ;                 // minimum energy of digit (ADC)
  Float_t fSampleQualityCut;         // Cut on sample shapes: 0: no samples; 1: default parameterization; 999: accept even obviously bad
  Float_t fSampleToSec ;             // Conversion coeff from sample time step to seconds
  Int_t fEmcCrystals ;               //  number of EMC crystals
  AliPHOSGeometry * fGeom ;          //! PHOS geometry
  static AliPHOSCalibData * fgCalibData ;     //! Calibration database if avalable
  AliPHOSPulseGenerator   * fPulseGenerator ; //! Class with pulse shape parameters
  AliRawReader            * fRawReader;       //! Raw data reader
  AliCaloRawStreamV3      * fRawStream;       //! Calorimeter decoder of ALTRO format
  Int_t *fADCValuesLG;               //! Array of low-gain ALTRO samples
  Int_t *fADCValuesHG;               //! Array of high-gain ALTRO samples
  static const Int_t fgkSTUDDL = 20; //! DDL ID of the PHOS STU

  ClassDef(AliPHOSRawDigiProducer,9)
};

#endif
