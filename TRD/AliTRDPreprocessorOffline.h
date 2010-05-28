#ifndef ALITRDPREPROCESSOROFFLINE_H
#define ALITRDPREPROCESSOROFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//
//    Class to create OCDB entries - processing the results of the OFFLINE calibration
//


#include "TNamed.h"
class TObjArray;
class TH2I;
class TProfile2D;
class AliTRDCalibraVdriftLinearFit;
class TH1I;
class TH2F;


class AliTRDPreprocessorOffline:public TNamed { 
public:
  enum{
    kGain = 0,
      kVdriftPHDet = 1,
      kVdriftPHPad = 2,
      kT0PHDet = 3,
      kT0PHPad = 4,
      kVdriftLinear = 5,
      kLorentzLinear = 6,
      kPRF = 7
      };   
  
  AliTRDPreprocessorOffline();
  virtual ~AliTRDPreprocessorOffline();

  void SetLinearFitForVdrift(Bool_t methodsecond) { fMethodSecond = methodsecond;};
  Bool_t GetLinearFitForVdrift() const { return fMethodSecond;};

  void CalibVdriftT0(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage="");
  void CalibGain(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage="");
  void CalibPRF(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage="");
  
  Bool_t ReadGainGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadVdriftT0Global(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadVdriftLinearFitGlobal(const Char_t* fileName="CalibObjects.root");
  Bool_t ReadPRFGlobal(const Char_t* fileName="CalibObjects.root");

  Bool_t AnalyzeGain(); 
  Bool_t AnalyzeVdriftT0(); 
  Bool_t AnalyzeVdriftLinearFit(); 
  Bool_t AnalyzePRF(); 
  
  void UpdateOCDBT0(Int_t startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBVdrift(Int_t startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateOCDBPRF(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);

  
private:
  Bool_t fMethodSecond;                   // Second Method for drift velocity   
  TH2I *fCH2d;                            // Gain
  TProfile2D *fPH2d;                      // Drift velocity first method
  TProfile2D *fPRF2d;                     // PRF
  AliTRDCalibraVdriftLinearFit *fAliTRDCalibraVdriftLinearFit; // Drift velocity second method
  TH1I *fNEvents;                         // Number of events 
  TH2F *fAbsoluteGain;                    // Absolute Gain calibration
  TObjArray * fPlots;                     // array with some plots to check
  TObjArray * fCalibObjects;              // array with calibration objects 


  

private:
  AliTRDPreprocessorOffline& operator=(const AliTRDPreprocessorOffline&); // not implemented
  AliTRDPreprocessorOffline(const AliTRDPreprocessorOffline&); // not implemented
  ClassDef(AliTRDPreprocessorOffline,1)
};

#endif
