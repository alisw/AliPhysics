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

class AliADQADataMakerRec: public AliQADataMakerRec {

public:
  // Histograms for Raw data control
  enum HRawType_t {kPedestalInt0,kPedestalInt1,
  		   kChargeEoI,kChargeEoIInt0,kChargeEoIInt1,
		   kWidth,kHPTDCTime,
		   kMultiADA,kMultiADC,kChargeADA,kChargeADC,kChargeAD, 
		   kADATime,kADCTime,kDiffTime,kTimeADAADC,
		   kNCoincADA,kNCoincADC,kPairDiffTime};
	
public:
  AliADQADataMakerRec() ;            // constructor
  AliADQADataMakerRec(const AliADQADataMakerRec& qadm) ;   
  AliADQADataMakerRec& operator = (const AliADQADataMakerRec& qadm) ;
  virtual ~AliADQADataMakerRec() {;} // destructor
  AliADCalibData *GetCalibData() const;
  virtual void   InitRaws() ; 
  
protected: 
  AliADCalibData *fCalibData;        //! calibration data
   
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
  Float_t fADCmean[32];              // mean adc per integrator
  size_t fTrendingUpdateTime;         // trending histos update time
  UInt_t fCycleStartTime;             // timestamp of QA start-of-cycle
  UInt_t fCycleStopTime;              // timestamp of QA end-of-cycle
  Float_t            fTimeOffset[16]; //! HPTDC time offsets channel by channel
  TF1*               fTimeSlewing;    //! Function for time slewing correction

  ClassDef(AliADQADataMakerRec,4)  // description 

};

#endif // AliADQADATAMAKERREC_H
