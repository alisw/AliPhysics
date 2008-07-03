#ifndef ALIMUONTRACKERCALIBRATEDDATAMAKER_H
#define ALIMUONTRACKERCALIBRATEDDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerCalibratedDataMaker
/// \brief Creator of calibrated AliMUONVTrackerData from AliRawReader
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVTRACKERDATAMAKER_H
#  include "AliMUONVTrackerDataMaker.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliRawReader;
class AliMUONCalibrationData;
class AliMUONDigitCalibrator;
class AliMUONDigitMaker;
class AliMUONVTrackerData;
class AliMUONVStore;
class AliMUONVDigitStore;

class AliMUONTrackerCalibratedDataMaker : public AliMUONVTrackerDataMaker
{
public:

  AliMUONTrackerCalibratedDataMaker(TRootIOCtor* ioCtor);
  
  AliMUONTrackerCalibratedDataMaker(Int_t runNumber,
                                    AliRawReader* reader, 
                                    const char* cdbpath=0x0,
                                    const char* calibMode=0x0,
                                    Bool_t histogram=kFALSE,
                                    Double_t xmin=0.0,
                                    Double_t xmax=4096.0,
                                    Bool_t useHPdecoder=kTRUE);
  
  AliMUONTrackerCalibratedDataMaker(AliRawReader* reader = 0x0, 
                                    const char* cdbpath=0x0,
                                    const char* calibMode=0x0,
                                    Bool_t histogram=kFALSE,
                                    Double_t xmin=0.0,
                                    Double_t xmax=4096.0,
                                    Bool_t useHPdecoder=kTRUE);

  virtual ~AliMUONTrackerCalibratedDataMaker();
  
  /// Whether we have a valid raw reader
  Bool_t IsValid() const { return fRawReader != 0x0; }
  
  /// Our data
  AliMUONVTrackerData* Data() const { return fAccumulatedData; }
  
  /// We can be run
  virtual Bool_t IsRunnable() const { return kTRUE; }
  
  /// Whether we are running or not
  virtual Bool_t IsRunning() const { return fIsRunning; }
  
  /// Set the running status
  virtual void SetRunning(Bool_t flag) { fIsRunning = flag; }
  
	Bool_t ProcessEvent();
	
  Bool_t NextEvent();
  
  void Print(Option_t* opt="") const;
  
  void Rewind();
  
  /// Get our source URI
  virtual TString Source() const { return fSource.Data(); }
  
  /// Set our source URI
  void SetSource(const char* source) { fSource = source; }
  
  /// Number of events seen
    Int_t NumberOfEvents() const { return fNumberOfEvents; }

  virtual Long64_t Merge(TCollection* list);
  
	void SetRawReader(AliRawReader* rawReader);
	
private:
  /// Not implemented
  AliMUONTrackerCalibratedDataMaker(const AliMUONTrackerCalibratedDataMaker& rhs);
  /// Not implemented
  AliMUONTrackerCalibratedDataMaker& operator=(const AliMUONTrackerCalibratedDataMaker& rhs);
  
  void Ctor(Int_t runNumber, const char* calibMode,
            Bool_t histogram, Double_t xmin, Double_t xmax);
  
private:
  AliRawReader* fRawReader; //!< reader of the data
	Bool_t fIsOwnerOfRawReader; //!< whether we must delete fRawReader or not
  AliMUONVTrackerData* fAccumulatedData; ///< data (owner)
  AliMUONVStore* fOneEventData; ///< data for one event (owner)
  TString fSource; ///< where the data comes from
  Bool_t fIsRunning; ///< whether we are running or are paused
  AliMUONDigitCalibrator* fDigitCalibrator; //!< digit calibrator (if calibrating data)
  AliMUONCalibrationData* fCalibrationData; ///< calibration data (if calibrating data)  
  TString fCDBPath; ///< OCDB path (if calibrating data)
  Int_t fNumberOfEvents; ///< number of events seen
  Bool_t fUseHPDecoder; ///< whether to use High Performance decoder or not
  static Int_t fgkCounter; ///< to count the number of instances
  
  ClassDef(AliMUONTrackerCalibratedDataMaker,3) // Producer of calibrated AliMUONVTrackerData from raw data
};

#endif
