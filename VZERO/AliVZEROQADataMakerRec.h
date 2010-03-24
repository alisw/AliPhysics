#ifndef ALIVZEROQADATAMAKERREC_H
#define ALIVZEROQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//  Produces the data needed to calculate the quality assurance  
//  All data must be mergeable objects
//  Handles ESDs and RAWs
//  Histos will be used for Raw Data control and monitoring

// --- ROOT system ---
class TH1F; 
class TH1I; 
class TObjArray; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliCDBManager;
class AliCDBStorage;
class AliVZEROCalibData;

class AliVZEROQADataMakerRec: public AliQADataMakerRec {

public:
  // Histograms for Raw data control
  enum HRawType_t {kPedestalInt0,kPedestalInt1,kPedestalCycleInt0,kPedestalCycleInt1,kPedestalTimeInt0,kPedestalTimeInt1
  		  ,kChargeEoIInt0,kChargeEoIInt1,kChargeEoITimeInt0,kChargeEoITimeInt1,kChargeEoICycleInt0,kChargeEoICycleInt1
		  ,kChargeEoIBBInt0,kChargeEoIBBInt1,kChargeEoIBGInt0,kChargeEoIBGInt1,kChargeVsClockInt0,kChargeVsClockInt1
		  ,kChargeMBBB0BG0Int0,kChargeMBBB0BG1Int0,kChargeMBBB1BG0Int0,kChargeMBBB1BG1Int0
		  ,kChargeMBBB0BG0Int1,kChargeMBBB0BG1Int1,kChargeMBBB1BG0Int1,kChargeMBBB1BG1Int1
		  ,kWidth,kWidthBB,kWidthBG,kHPTDCTime,kHPTDCTimeBB,kHPTDCTimeBG,kBBFlagVsClock,kBGFlagVsClock
		  ,kMultiV0A,kMultiV0C,kChargeV0A,kChargeV0C,kChargeV0 
		  ,kV0ATime,kV0CTime,kDiffTime
		  ,kRawMIPV0A,kRawMIPV0C,kRawMIPV0,kRawMIPChannel
		  ,kRawMeanChargePerRing,kRawMeanFlagPerRing,kRawDQMCharge,kRawDQMFlag} ;
	
 enum HESDType_t {kCellMultiV0A,kCellMultiV0C,kMIPMultiV0A,kMIPMultiV0C,kMIPMultiChannel
		  ,kBBFlag,kBGFlag,kChargeChannel,kTimeChannel
		  ,kESDV0ATime,kESDV0CTime,kESDDiffTime};

public:
  AliVZEROQADataMakerRec() ;            // constructor
  AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) ;   
  AliVZEROQADataMakerRec& operator = (const AliVZEROQADataMakerRec& qadm) ;
  virtual ~AliVZEROQADataMakerRec() {;} // destructor
  AliVZEROCalibData *GetCalibData() const;
  virtual void   InitRaws() ; 
  void SetTrendingUpdateTime(size_t time) {fTrendingUpdateTime = time;};
  
protected: 
  AliVZEROCalibData *fCalibData;        //! calibration data
   
private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;
  virtual void   InitESDs() ; 
  virtual void   InitDigits();  
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeDigits() ; 
  virtual void   MakeDigits(TTree* digitTree) ; 
  virtual void   StartOfDetectorCycle() ; 
  void AddTrendingEntry();

  Int_t   fEvent;                     // event index
  Int_t   fEven[64];                  // even charge integrators
  Int_t   fOdd[64];                   // odd charge intergators
  Float_t fADCmean[128];              // mean adc per integrator
  size_t fNTotEvents;                 // total number of events
  size_t fNSubEvents;                 // number of events used in trending histos
  size_t fTrendingUpdateEvent;        // event index of last update of the trending histos
  size_t fNTrendingUpdates;           // number of updates in trending histos
  size_t fTrendingUpdateTime;         // trending histos update time
  UInt_t fCycleStartTime;             // timestamp of QA start-of-cycle
  UInt_t fCycleStopTime;              // timestamp of QA end-of-cycle
  Double_t fMonitorRate;              // monitoring rate
  Double_t fChargePerRing[8];         // charge per ring
  Double_t fFlagPerRing[8];           // flag per ring
  Double_t fChargePerChannel[64];     // charge per channel
  Double_t fFlagPerChannel[64];       // flag per channel
  Double_t fMeanChargePerChannel[64]; // mean charge per channel
  Double_t fMeanFlagPerChannel[64];   // mean flag per channel

  ClassDef(AliVZEROQADataMakerRec,2)  // description 

};

#endif // AliVZEROQADATAMAKERREC_H
