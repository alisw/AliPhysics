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
		  ,kRawMIPV0A,kRawMIPV0C,kRawMIPV0,kRawMIPChannel} ;
	
 enum HESDType_t {kCellMultiV0A,kCellMultiV0C,kMIPMultiV0A,kMIPMultiV0C,kMIPMultiChannel
		  ,kBBFlag,kBGFlag,kChargeChannel,kTimeChannel
		  ,kESDV0ATime,kESDV0CTime,kESDDiffTime};

public:
  AliVZEROQADataMakerRec() ;            // constructor
  AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) ;   
  AliVZEROQADataMakerRec& operator = (const AliVZEROQADataMakerRec& qadm) ;
  virtual ~AliVZEROQADataMakerRec() {;} // destructor
  AliVZEROCalibData *GetCalibData() const;
  
protected: 
  AliVZEROCalibData *fCalibData;        //! calibration data
   
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;
  virtual void   InitESDs() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 

  Int_t   fEvent;
  Int_t   fEven[64];
  Int_t   fOdd[64];
  Float_t fADCmean[128];

  ClassDef(AliVZEROQADataMakerRec,1)  // description 

};

#endif // AliVZEROQADATAMAKERREC_H
