#ifndef ALIT0QACHECKER_H
#define ALIT0QACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//


// --- ROOT system ---
class TFile ; 
class TH1F ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliT0QAChecker: public AliQACheckerBase {

public:
  AliT0QAChecker();    
  AliT0QAChecker(const AliT0QAChecker& qac);
  AliT0QAChecker& operator=(const AliT0QAChecker& qac);  
  // dtor
  virtual ~AliT0QAChecker();
  //Double_t CheckLaser(TObjArray *listrec ) const ;
  Double_t CheckRaw(TObjArray *listrec ) const ;
  Double_t CheckESD(TObjArray *listrec ) const ;
   
 
private:

  enum{
    kT0Fatal=-1,  ///< error is really serious
    kT0Error=0,   ///< normal error, i.e. something is wrong
    kT0Warning=1, ///< warning, i.e. might become an error later on
    kT0Info=2     ///< just so you know...
  };

  virtual void Check(Double_t * test, AliQAv1::ALITASK_t, TObjArray ** list, const AliDetectorRecoParam * recoParam) ;
  void EraseOldMessages(TH1* h) const;
  Double_t ConvertQualityFlagToDouble(int qualityFlag) const;
  Float_t GetMeanAboveThreshold(TH1F* hV, Float_t thr) const;
  void GetMeanAndRmsAroundMainMaximum(Float_t &meanHisto,Float_t &rmsHisto, TH1F *histo, int type) const;

  Float_t fMeanCFDFromGoodRunParam[24]; //mean CFD for each PMT from a good run
  Float_t fMeanLEDFromGoodRunParam[24]; //mean LED for each PMT from a good run
  Float_t fMeanQTCFromGoodRunParam[24]; //mean QTC for each PMT from a good run
  Float_t fCFDErrorThreshold; //CFD error threshold instead of the yellow band 
  Float_t fLEDErrorThreshold; //LED error threshold 
  Float_t fQTCErrorThreshold; //QTC error threshold 
  Float_t fRatioCFDEffLEDEffErrorThreshold; //ratio CFD to LED efficiency error threshold 
  Float_t fQTCEfficiencyErrorThreshold; //QTC efficiency error threshold 
  Int_t   fBCIDPeriodParam; // period 
  Int_t   fBCIDOffsetParam;//offset of TRM BCID 
  Int_t   fBCIDBandWidthParam; // tollerated deviation of BCID from diagonal 
  Float_t fTZeroAPlusCErrorThreshold; // constraint on the tzero vertex displacement in ps
  Float_t fTZeroAMinusCErrorThreshold; // constraint on the tzero time shift in ps 


 
  ClassDef(AliT0QAChecker,1)  // description 

};

#endif // AliT0QAChecker_H
