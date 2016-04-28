#ifndef ALIADQADATAMAKERREC_H
#define ALIADQADATAMAKERREC_H
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
class AliADCalibData;
class AliADRecoParam;
class AliADQAParam;

class AliADQADataMakerRec: public AliQADataMakerRec {

public:
  // Histograms for Raw data control
  enum HRawType_t {
  		   kChargeADA,kChargeADC,kChargeEoI,kChargeEoIBB,kChargeEoIBG,
		   kHPTDCTime,kHPTDCTimeBB,kHPTDCTimeBG,kWidth,
		   kHPTDCTimeRebin,kHPTDCTimeRebinBB,kHPTDCTimeRebinBG,
		   kBBFlagVsClock,kBBFlagVsClock_ADOR,kBGFlagVsClock,kBBFlagsPerChannel,kBGFlagsPerChannel,
		   kChargeVsClockInt0,kChargeVsClockInt1,kMaxChargeClock,
		   kNBBCoincADA,kNBBCoincADC,kNBGCoincADA,kNBGCoincADC,
		   kPedestalDiffInt0,kPedestalDiffInt1,
		   kChargeEoIInt0,kChargeEoIInt1,kChargeSaturation,
		   kNBBCoincCorr,kNBGCoincCorr,
		   kTriggers,kDecisions,
		   kMeanTimeADA,kMeanTimeADC,kMeanTimeDiff,kMeanTimeCorr,kMeanTimeSumDiff,
		   kPedestalInt0,kPedestalInt1,
		   kNEventsBBFlag,kNEventsBGFlag,
		   kFlagNoTime,kTimeNoFlag,
		   kWidthBB,kWidthBG,
		   kTimeSlewingADA,kTimeSlewingADC,kWidthSlewing,
		   kMultiADA,kMultiADC,kChargeAD, 
		   
		   kChargeADA_PC,kChargeADC_PC,
		   kTrend_TriggerChargeQuantileADA,kTrend_TriggerChargeQuantileADC,
		   
		   kPairTimeDiffMean,kPairTimeDiffRMS,
		   kNChargeCorrADA,
		   kNChargeCorrADC = kNChargeCorrADA + 28,
		   kNTimeCorrADA = kNChargeCorrADC + 28,
		   kNTimeCorrADC = kNTimeCorrADA + 28,
		   kNTimeDiffADA = kNTimeCorrADC + 28,
		   kNTimeDiffADC = kNTimeDiffADA + 28};
		   
		   
  enum HESDType_t {kCellMultiADA,kCellMultiADC,
		   kBBFlag,kBGFlag,kChargeChannel,kTimeChannel,
		   kESDADATime,kESDADCTime,kESDDiffTime,kESDADATimeVsCharge,kESDADCTimeVsCharge,
		   kESDADAPairTimeSumDiff,kESDADCPairTimeSumDiff};
	
public:
  AliADQADataMakerRec() ;            // constructor
  AliADQADataMakerRec(const AliADQADataMakerRec& qadm) ;   
  AliADQADataMakerRec& operator = (const AliADQADataMakerRec& qadm) ;
  virtual ~AliADQADataMakerRec() {;} // destructor
  AliADCalibData *GetCalibData() const;
  AliADQAParam *GetQAParam() const;
  virtual void   InitRaws() ; 
  
protected: 
  AliADCalibData *fCalibData;        //! calibration data
  AliADRecoParam *fRecoParam;
  AliADQAParam *fQAParam;
   
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
  
  Int_t   fEven[16];                  // even charge integrators
  Int_t   fOdd[16];                   // odd charge intergators
  Float_t fADCmean[32];               // mean adc per integrator
  size_t fTrendingUpdateTime;         // trending histos update time
  UInt_t fCycleStartTime;             // timestamp of QA start-of-cycle
  UInt_t fCycleStopTime;              // timestamp of QA end-of-cycle		      
  Float_t fADADist;     	      // Z position of ADA
  Float_t fADCDist;     	      // Z position of ADC
  UInt_t fOldRun;

  ClassDef(AliADQADataMakerRec,4)  // description 

};

#endif // AliADQADATAMAKERREC_H
