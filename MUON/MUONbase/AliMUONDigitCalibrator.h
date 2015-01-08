/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONDigitCalibrator
/// \brief Class to calibrate the digits
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONDIGITCALIBRATOR_H
#define ALIMUONDIGITCALIBRATOR_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliMUONCalibrationData;
class AliMUONLogger;
class AliMUONVStore;
class AliMUONVDigitStore;
class AliMUONVDigit;
class AliMUONPadStatusMaker;
class AliMUONPadStatusMapMaker;
class AliMUONRecoParam;
class TExMap;

class AliMUONDigitCalibrator : public TObject
{
public:

  AliMUONDigitCalibrator(Int_t runNumber);
  
  AliMUONDigitCalibrator(const AliMUONCalibrationData& calib, 
                         const AliMUONRecoParam* recoParams);
  
  AliMUONDigitCalibrator(const AliMUONCalibrationData& calib, int b);
  
  virtual ~AliMUONDigitCalibrator();
  
  virtual void Calibrate(AliMUONVDigitStore& digitStore);

  Bool_t IsValidDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel, 
                      Int_t* statusMap=0x0) const;

  Float_t CalibrateDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                         Float_t adc, Float_t nsigmas=-1, 
                         Bool_t* isSaturated=0x0) const;
                     
  Int_t PadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;

  Int_t StatusMap(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;

private:    
  
  /// Not implemented
  AliMUONDigitCalibrator(const AliMUONDigitCalibrator& other);
  /// Not implemented
  AliMUONDigitCalibrator& operator=(const AliMUONDigitCalibrator& other);
  
  void Ctor(const AliMUONCalibrationData& calib,
            const AliMUONRecoParam* recoParams,
            Bool_t deferredInitialization=kTRUE);
  
private:
	AliMUONLogger* fLogger; //!< to log repeated messages
	AliMUONPadStatusMaker* fStatusMaker; //!< to build pad statuses
	AliMUONPadStatusMapMaker* fStatusMapMaker; //!< to build status map
	AliMUONVStore* fPedestals; //!< pedestal values
	AliMUONVStore* fGains; //!< gain values
	Int_t fApplyGains; //!< whether we should apply gains or not, capa or not...
	AliMUONVStore* fCapacitances; //!< capa values
	Double_t fNumberOfBadPads; //!< # of times we've rejected a bad pad
	Double_t fNumberOfPads; //!< # of pads we've seen
	Double_t fChargeSigmaCut; //!< number of sigmas to cut on
  UInt_t fMask; //!< mask used to cut bad channels
  
	static const Int_t fgkNoGain; //!< do not apply gain calib at all
	static const Int_t fgkGainConstantCapa; //!< apply gain (from OCDB) with constant capa
	static const Int_t fgkGain; //!< apply gain and capa (from OCDB)
	static const Int_t fgkInjectionGain; //!< apply injection gain (from OCDB)
	
  ClassDef(AliMUONDigitCalibrator,10) // Calibrate raw digit
};

#endif
