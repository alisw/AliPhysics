#ifndef ALITRDPREPROCESSOROFFLINE_H
#define ALITRDPREPROCESSOROFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//
//    Class to create OCDB entries - processing the results of the OFFLINE calibration
//


#include "TNamed.h"
#include "AliTRDCalChamberStatus.h"
#include <THnSparse.h>
class TObjArray;
class AliTRDCalDet;
class TH2I;
class TProfile2D;
class AliTRDCalibraVdriftLinearFit;
class AliTRDCalibraExbAltFit;
class TH1I;
class TH2F;
class TString;


class AliTRDPreprocessorOffline:public TNamed { 
public:
  enum{ kGain = 0,
	kVdriftPHDet = 1,
	kVdriftPHPad = 2,
	kT0PHDet = 3,
	kT0PHPad = 4,
	kVdriftLinear = 5,
	kLorentzLinear = 6,
	kChamberStatus = 7,
	kPRF = 8,
	kExbAlt = 9,
	kNumCalibObjs = 10
  };
  enum { kGainNotEnoughStatsButFill = 2,
	 kVdriftNotEnoughStatsButFill = 4,
	 kGainNotEnoughStatsNotFill = 8,
	 kVdriftNotEnoughStatsNotFill = 16,
	 kTimeOffsetNotEnoughStatsNotFill = 32,
	 kExBErrorRange = 64};  
  enum { kGainErrorOld = 2,
	 kVdriftErrorOld = 4,
	 kExBErrorOld = 8,
	 kGainErrorRange = 16,
	 kVdriftErrorRange = 32,
	 kTimeOffsetErrorRange = 64,
	 kChamberStatusErrorRange = 128};  


  AliTRDPreprocessorOffline();
  virtual ~AliTRDPreprocessorOffline();

  Bool_t Init(const Char_t* fileName);
  void Process(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage);

  // settings
  void SetLinearFitForVdrift(Bool_t methodsecond) { fMethodSecond = methodsecond;};
  void SetNameList(TString nameList) { fNameList = nameList;};
  void     SetMinStatsVdriftT0PH(Int_t minStatsVdriftT0PH)           { fMinStatsVdriftT0PH = minStatsVdriftT0PH; }  
  void     SetMinStatsVdriftLinear(Int_t minStatsVdriftLinear)       { fMinStatsVdriftLinear = minStatsVdriftLinear; }  
  void     SetMinStatsGain(Int_t minStatsGain)                       { fMinStatsGain = minStatsGain; }  
  void     SetMinStatsPRF(Int_t minStatsPRF)                         { fMinStatsPRF = minStatsPRF; }  
  void     SetMinStatsChamberStatus(Int_t minStatsChamberStatus)     { fMinStatsChamberStatus = minStatsChamberStatus; }
  void     SetLimitValidateNoData(Int_t nodatavalidate)              { fNoDataValidate = nodatavalidate; };
  void     SetLimitValidateBadCalib(Int_t badcalibvalidate)          { fBadCalibValidate = badcalibvalidate; };
  void     SetBackCorrectGain(Bool_t backCorrectGain)                { fBackCorrectGain = backCorrectGain; }
  void     SetBackCorrectVdrift(Bool_t backCorrectVdrift)            { fBackCorrectVdrift = backCorrectVdrift; }
  void     SetNoExBUsedInReco(Bool_t noExBUsedInReco)                { fNoExBUsedInReco    = noExBUsedInReco;   };
  void     SetSwitchOnValidation(Bool_t switchOnValidation)          { fSwitchOnValidation = switchOnValidation;};
  void     SetRMSBadCalibratedGain(Double_t rms)                     { fRMSBadCalibratedGain = rms;};
  void     SetRMSBadCalibratedVdrift(Double_t rms)                   { fRMSBadCalibratedVdrift = rms;};
  void     SetRMSBadCalibratedExB(Double_t rms)                      { fRMSBadCalibratedExB = rms;};
  void     SetRobustFitDriftVelocity(Bool_t robustFitDriftVelocity)  { fRobustFitDriftVelocity = robustFitDriftVelocity;};
  void     SetRobustFitExbAlt(Bool_t robustFitExbAlt)                { fRobustFitExbAlt = robustFitExbAlt;};


  Bool_t GetLinearFitForVdrift() const { return fMethodSecond;};
  TString GetNameList() const { return fNameList;}; 

  // status
  Bool_t      CheckStatus(Int_t status, Int_t bitMask) const;
  void PrintStatus() const;

  Bool_t      IsGainNotEnoughStatsButFill() const 
    { return CheckStatus(fStatusNeg, kGainNotEnoughStatsButFill);  };
  Bool_t      IsGainNotEnoughStatsNotFill() const 
    { return CheckStatus(fStatusNeg, kGainNotEnoughStatsNotFill);  };
  Bool_t      IsVdriftNotEnoughStatsButFill() const 
    { return CheckStatus(fStatusNeg, kVdriftNotEnoughStatsButFill);  };
  Bool_t      IsVdriftNotEnoughStatsNotFill() const 
    { return CheckStatus(fStatusNeg, kVdriftNotEnoughStatsNotFill);  };
  Bool_t      IsTimeOffsetNotEnoughStatsNotFill() const 
    { return CheckStatus(fStatusNeg, kTimeOffsetNotEnoughStatsNotFill);  };
  Bool_t      IsExBErrorRange() const 
    { return CheckStatus(fStatusNeg, kExBErrorRange);  };
  
  Bool_t      IsGainErrorOld() const 
    { return CheckStatus(fStatusPos, kGainErrorOld);  };
  Bool_t      IsVdriftErrorOld() const 
    { return CheckStatus(fStatusPos, kVdriftErrorOld);  };
  Bool_t      IsExBErrorOld() const 
    { return CheckStatus(fStatusPos, kExBErrorOld);  };
  Bool_t      IsGainErrorRange() const 
    { return CheckStatus(fStatusPos, kGainErrorRange);  };
  Bool_t      IsVdriftErrorRange() const 
    { return CheckStatus(fStatusPos, kVdriftErrorRange);  };
  Bool_t      IsTimeOffsetErrorRange() const 
    { return CheckStatus(fStatusPos, kTimeOffsetErrorRange);  };
  Bool_t      IsChamberStatusErrorRange() const 
    { return CheckStatus(fStatusPos, kChamberStatusErrorRange);  };
  

  // Back corrections
  void SetCalDetGain(AliTRDCalDet *calDetGainUsed) {fCalDetGainUsed = calDetGainUsed;};
  void SetCalDetVdrift(AliTRDCalDet *calDetVdriftUsed);
  void SetCalDetVdriftExB(AliTRDCalDet *calDetVdriftUsed,AliTRDCalDet *calDetExBUsed) {fCalDetVdriftUsed = calDetVdriftUsed; fCalDetExBUsed = calDetExBUsed;};
  Bool_t SetCalDetGain(Int_t runNumber, Int_t version, Int_t subversion);
  Bool_t SetCalDetVdriftExB(Int_t runNumber, Int_t versionv, Int_t subversionv, Int_t versionexb, Int_t subversionexb);
  
  AliTRDCalDet *GetCalDetGain() const { return fCalDetGainUsed;};
  AliTRDCalDet *GetCalDetVdrift() const { return fCalDetVdriftUsed;};
  Int_t    GetFirstRunGainUsed() const                               { return fFirstRunGainUsed;       }
  Int_t    GetVersionGainUsed() const                                { return fVersionGainUsed;        }
  Int_t    GetSubVersionGainUsed() const                             { return fSubVersionGainUsed;     }
  Int_t    GetFirstRunVdriftUsed() const                             { return fFirstRunVdriftUsed;     }
  Int_t    GetVersionVdriftUsed() const                              { return fVersionVdriftUsed;      }
  Int_t    GetSubVersionVdriftUsed() const                           { return fSubVersionVdriftUsed;   }
  Int_t    GetFirstRunExBUsed() const                                { return fFirstRunExBUsed;        }
  Int_t    GetVersionExBUsed() const                                 { return fVersionExBUsed;         }
  Int_t    GetSubVersionExBUsed() const                              { return fSubVersionExBUsed;      }


  // Internal functions

  void CalibVdriftT0(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage="");
  void CalibExbAlt(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage="");
  void CalibGain(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage="");
  void CalibPRF(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage="");
  void CalibChamberStatus(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage="");

  Bool_t ReadStatusGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadGainGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadVdriftT0Global(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadVdriftLinearFitGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadExbAltFitGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadPRFGlobal(const Char_t* fileName="CalibObjects.root");

  Bool_t AnalyzeGain(); 
  Bool_t AnalyzeVdriftT0(); 
  Bool_t AnalyzeVdriftLinearFit(); 
  Bool_t AnalyzeExbAltFit();
  Bool_t AnalyzePRF();
  Bool_t AnalyzeChamberStatus(); 
  
  void CorrectFromDetGainUsed();
  void CorrectFromDetVdriftUsed();
  
  void UpdateOCDBT0(Int_t startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBVdrift(Int_t startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBExB(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath);
  void UpdateOCDBExBAlt(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath);
  void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBPRF(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBChamberStatus(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath);

  Bool_t ValidateGain();
  Bool_t ValidateVdrift();
  Bool_t ValidateExB();
  Bool_t ValidateT0();
  Bool_t ValidatePRF() const;
  Bool_t ValidateChamberStatus();

  Int_t    GetStatus() const;
  Int_t    GetStatusPos() const                                      { return fStatusPos;              }
  Int_t    GetStatusNeg() const                                      { return fStatusNeg;              }
 

  
 private:
  Bool_t fMethodSecond;                      // Second Method for drift velocity   
  TString fNameList;                         // Name of the list
  AliTRDCalDet *fCalDetGainUsed;             // CalDet used and to be corrected for
  AliTRDCalDet *fCalDetVdriftUsed;           // CalDet used and to be corrected for
  AliTRDCalDet *fCalDetExBUsed;              // CalDet used and to be corrected for
  TH2I *fCH2d;                               // Gain
  TProfile2D *fPH2d;                         // Drift velocity first method
  TProfile2D *fPRF2d;                        // PRF
  THnSparseI *fSparse;                       // chamberstatus
  AliTRDCalibraVdriftLinearFit *fAliTRDCalibraVdriftLinearFit; // Drift velocity second method
  AliTRDCalibraExbAltFit* fAliTRDCalibraExbAltFit; //ExB alternative method
  TH1I *fNEvents;                         // Number of events 
  TH2F *fAbsoluteGain;                    // Absolute Gain calibration
  TObjArray * fPlots;                     // array with some plots to check
  TObjArray * fCalibObjects;              // array with calibration objects 
  Int_t    fFirstRunGainUsed;             // first run gain used 
  Int_t    fVersionGainUsed;              // VersionGainUsed 
  Int_t    fSubVersionGainUsed;           // SubVersionGainUsed
  Int_t    fFirstRunVdriftUsed;           // FirstRunVdrift 
  Int_t    fVersionVdriftUsed;            // VersionVdriftUsed 
  Int_t    fSubVersionVdriftUsed;         // SubVersionVdriftUsed
  Int_t    fFirstRunExBUsed;              // FirstRunExB 
  Int_t    fVersionExBUsed;               // VersionExBUsed 
  Int_t    fSubVersionExBUsed;            // SubVersionExBUsed
  Bool_t   fNoExBUsedInReco;              // ExB not used yet in the reco
  Bool_t   fSwitchOnValidation;           // Validation
  Bool_t   fVdriftValidated;              // Vdrift validation
  Bool_t   fExBValidated;                 // ExB validation
  Bool_t   fT0Validated;                  // T0 validation
  Int_t    fMinStatsVdriftT0PH;           // MinStats VdriftT0
  Int_t    fMinStatsVdriftLinear;         // MinStats Vdrift Linear
  Int_t    fMinStatsGain;                 // MinStats Gain
  Int_t    fMinStatsPRF;                  // MinStats PRF
  Int_t    fMinStatsChamberStatus;        // MinStats ChamberStatus
  Bool_t   fBackCorrectGain;              // Back correction afterwards gain  
  Bool_t   fBackCorrectVdrift;            // Back correction afterwards vdrift
  Bool_t   fNotEnoughStatisticsForTheGain;// Take the chamber per chamber distribution from the default distribution
  Bool_t   fNotEnoughStatisticsForTheVdriftLinear;// Take the chamber per chamber distribution from the default distribution
  Int_t    fStatusNeg;                    // Info but ok
  Int_t    fStatusPos;                    // Problems
  Int_t    fBadCalib[18];                 // number of bad calibrated chambers per sm
  Int_t    fNoData[18];                   // number of  chambers w/o data per sm
  Int_t    fBadCalibValidate;             // validation limit for bad calibrated chambers
  Int_t    fNoDataValidate;               // validation limit for chamber w/o data (sm w/o data excluded)
  Double_t fRMSBadCalibratedGain;         // value to decide when it is bad calibrated 
  Double_t fRMSBadCalibratedVdrift;       // value to decide when it is bad calibrated 
  Double_t fRMSBadCalibratedExB;          // value to decide when it is bad calibrated 
  Bool_t   fRobustFitDriftVelocity;       // Robust fit for the drift velocity
  Bool_t   fRobustFitExbAlt;              // Robust fit for the exb alt 

  Int_t GetSubVersion(TString name) const;
  Int_t GetVersion(TString name) const;
  Int_t GetFirstRun(TString name) const;

  

private:
  AliTRDPreprocessorOffline& operator=(const AliTRDPreprocessorOffline&); // not implemented
  AliTRDPreprocessorOffline(const AliTRDPreprocessorOffline&); // not implemented
  ClassDef(AliTRDPreprocessorOffline,3)
};

#endif


