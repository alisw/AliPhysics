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
class AliVZEROTriggerData;

class AliVZEROQADataMakerRec: public AliQADataMakerRec {

public:
  // Histograms for Raw data control
  enum HRawType_t {kPedestalInt0,kPedestalInt1
  		  ,kChargeEoI,kChargeEoIInt0,kChargeEoIInt1
		  ,kChargeEoIBBInt0,kChargeEoIBBInt1,kChargeEoIBGInt0,kChargeEoIBGInt1,kChargeVsClockInt0,kChargeVsClockInt1
		  ,kChargeMBBB0BG0Int0,kChargeMBBB0BG1Int0,kChargeMBBB1BG0Int0,kChargeMBBB1BG1Int0
		  ,kChargeMBBB0BG0Int1,kChargeMBBB0BG1Int1,kChargeMBBB1BG0Int1,kChargeMBBB1BG1Int1
		  ,kWidth,kWidthBB,kWidthBG,kHPTDCTime,kHPTDCTimeBB,kHPTDCTimeBG,kBBFlagVsClock,kBGFlagVsClock
		  ,kMultiV0A,kMultiV0C,kChargeV0A,kChargeV0C,kChargeV0 
		  ,kV0ATime,kV0CTime,kDiffTime
		  ,kRawMIPV0A,kRawMIPV0C,kRawMIPV0,kRawMIPChannel
		  ,kBBFlagsPerChannel, kTriggers,kTriggers2,kTimeV0AV0C
		  ,kCentrChargeV0AV0C};
	
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
  
protected: 
  AliVZEROCalibData *fCalibData;        //! calibration data
  AliVZEROTriggerData *fTriggerData;    //! trigger config data
   
private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;
  virtual void   InitESDs() ; 
  virtual void   InitDigits();  
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeDigits() ; 
  virtual void   MakeDigits(TTree* digitTree) ; 
  virtual void   StartOfDetectorCycle() ; 
  Float_t CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const;
  
 
  //  Int_t   fEvent;                     // event index
  Int_t   fEven[64];                  // even charge integrators
  Int_t   fOdd[64];                   // odd charge intergators
  Float_t fADCmean[128];              // mean adc per integrator
  //  size_t fNTotEvents;                 // total number of events
  //  size_t fNSubEvents;                 // number of events used in trending histos
  //  size_t fTrendingUpdateEvent;        // event index of last update of the trending histos
  //  size_t fNTrendingUpdates;           // number of updates in trending histos
  size_t fTrendingUpdateTime;         // trending histos update time
  UInt_t fCycleStartTime;             // timestamp of QA start-of-cycle
  UInt_t fCycleStopTime;              // timestamp of QA end-of-cycle
  Float_t            fTimeOffset[64]; //! HPTDC time offsets channel by channel
  TF1*               fTimeSlewing;    //! Function for time slewing correction

  ClassDef(AliVZEROQADataMakerRec,3)  // description 

};

#endif // AliVZEROQADATAMAKERREC_H
