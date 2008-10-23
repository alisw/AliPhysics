#ifndef ALIVZEROQADataMakerRec_H
#define ALIVZEROQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
*/


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
  //Histograms for Raw data control
  enum HRawType_t {kPedestal_Int0,kPedestal_Int1,kPedestal_Cycle_Int0,kPedestal_Cycle_Int1,kPedestal_Time_Int0,kPedestal_Time_Int1
  			,kChargeEoI_Int0,kChargeEoI_Int1,kChargeEoI_Time_Int0,kChargeEoI_Time_Int1,kChargeEoI_Cycle_Int0,kChargeEoI_Cycle_Int1
			,kChargeEoI_BB_Int0,kChargeEoI_BB_Int1,kChargeEoI_BG_Int0,kChargeEoI_BG_Int1,kChargeVsClock_Int0,kChargeVsClock_Int1
			,kChargeMB_BB0_BG0_Int0,kChargeMB_BB0_BG1_Int0,kChargeMB_BB1_BG0_Int0,kChargeMB_BB1_BG1_Int0
			,kChargeMB_BB0_BG0_Int1,kChargeMB_BB0_BG1_Int1,kChargeMB_BB1_BG0_Int1,kChargeMB_BB1_BG1_Int1
			,kWidth,kWidth_BB,kWidth_BG,kHPTDCTime,kHPTDCTime_BB,kHPTDCTime_BG,kBBFlagVsClock,kBGFlagVsClock
			,kMultiV0A,kMultiV0C,kChargeV0A,kChargeV0C,kChargeV0 
			,kV0ATime,kV0CTime,kDiffTime
			,kRawMIPV0A,kRawMIPV0C,kRawMIPV0,kRawMIPChannel} ;
	
	enum HESDType_t {kCellMultiV0A,kCellMultiV0C,kMIPMultiV0A,kMIPMultiV0C,kMIPMultiChannel
			,kBBFlag,kBGFlag,kChargeChannel,kTimeChannel
			,kESDV0ATime,kESDV0CTime,kESDDiffTime};

public:
  AliVZEROQADataMakerRec() ;           // constructor
  AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) ;   
  AliVZEROQADataMakerRec& operator = (const AliVZEROQADataMakerRec& qadm) ;
  virtual ~AliVZEROQADataMakerRec() {;} // destructor
  AliVZEROCalibData *GetCalibData() const;
  
protected: 
  AliVZEROCalibData *fCalibData;  //! calibration data
   
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray * list) ;
  virtual void   InitESDs() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 

  Int_t   fEvent;
  Int_t   fEven[64];
  Int_t   fOdd[64];
  Float_t fADC_Mean[128];

  ClassDef(AliVZEROQADataMakerRec,1)  // description 

};

#endif // AliVZEROQADataMakerRec_H
