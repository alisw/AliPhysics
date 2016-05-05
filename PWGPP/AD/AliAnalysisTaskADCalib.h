// -*- C++ -*-
// $Id$

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            This task is for AD calibration
//-----------------------------------------------------------------

class TList;
class TTree;
class TH2;
class TH3;
class TF1;
class AliCDBStorage;
class AliCDBEntry;
class ADESDFriendUtils;

#include <TString.h>

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskADCalib : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskADCalib(const char *name="AliAnalysisTaskADCalib");
  virtual ~AliAnalysisTaskADCalib();
  
  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetBCRangeTail(Int_t bcMin, Int_t bcMax) {
    fBCRangeTail[0] = bcMin;
    fBCRangeTail[1] = bcMax;
  }
  void SetBCRangeExtrapolation(Int_t bcMin, Int_t bcMax) {
    fBCRangeExtrapolation[0] = bcMin;
    fBCRangeExtrapolation[1] = bcMax;
  }

  void  ProcessOutput(const Char_t *filename, AliCDBStorage* db, Int_t runNumber);
  Int_t GetStatus() const;

protected:
  Bool_t  FillHist(TString name, Double_t x, Double_t y);             // TH2
  Bool_t  FillHist(TString name, Double_t x, Double_t y, Double_t z); // TH3
  TString GetHistName (Int_t ch, Int_t bc, Bool_t integrator) const;  // for TH2
  TString GetHistName (Int_t ch, Bool_t integrator) const;            // for TH3
  TString GetFcnName  (Int_t ch, Int_t bc, Bool_t integrator) const;
  TString GetHistTitle(Int_t ch, Int_t bc, Bool_t integrator) const;  // for TH2
  TString GetHistTitle(Int_t ch, Bool_t integrator) const;            // for TH3
  Bool_t  MakeExtrapolationFit(TH2* h, TF1* f, Int_t ch, Int_t bc, Double_t &xMax);

private:
  TTree*       MakeSaturationCalibObject();
  AliCDBEntry* UpdateGainParameters(Int_t runNumber, TTree *t);

  // not implemented:
  AliAnalysisTaskADCalib(const AliAnalysisTaskADCalib&);
  AliAnalysisTaskADCalib& operator=(const AliAnalysisTaskADCalib&);

  Int_t fBCRangeTail[2];          // BC range for tail charge
  Int_t fBCRangeExtrapolation[2]; // BC range for extrapolated BCs

  Float_t fTimeResolution[2];     //! HPTDC time resolution per side

  ADESDFriendUtils *fADESDFriendUtils; //! AD ESD friend helper 
  TList            *fList;             //! output histograms
                                       // TH2s (charge in a BC vs. tail charge)
                                       // TH3s (time1,time2,tail charge)

  enum EStatusCode_t {
    kOk,
    kInputError,       // open file error, missing histos
    kDataError,        // problems with histo information
    kMeasurementError, // MakeCalibObject returns NULL
    kLowStatistics,    // too low statistics
    kStoreError        // problems storing OCDB
  } ;
  Int_t fStatus; //! calibration status (after ProcessOutput)

  ClassDef(AliAnalysisTaskADCalib, 3);
} ;
