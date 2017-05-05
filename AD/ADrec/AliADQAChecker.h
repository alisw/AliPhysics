// -*- C++ -*-
#ifndef ALIADQACHECKER_H
#define ALIADQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality of the data
  by comparing with reference data
  which should be loaded from QA ref DB
*/

// --- ROOT system ---
class TH2;
class TObjArray;
class TVirtualPad;
class TCanvas;
class TSpline3;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
#include "AliADQAParam.h"

class AliADQAChecker: public AliQACheckerBase {
public:
  AliADQAChecker();
  virtual ~AliADQAChecker(); // destructor

  virtual void Init(const AliQAv1::DETECTORINDEX_t det);

  const AliADQAParam *GetQAParam() const;

protected:
  virtual void Check(Double_t* test, AliQAv1::ALITASK_t index, TObjArray** list, const AliDetectorRecoParam* recoParam);
  Double_t CheckRaws(TObjArray* list) const;
  Double_t CheckPedestals(TObjArray* list) const;
  Double_t CheckEsds(TObjArray* list) const;

  virtual void MakeImage(TObjArray** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ;
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t* value) const ;

  TCanvas* CreatePads(TCanvas *) const;
  TVirtualPad* GetPadByName(TCanvas *, const char*) const;

  TSpline3* MakeTimeSlewingSpline(TH2* ) const;
private:
  AliADQAChecker(const AliADQAChecker& qac);            // not implemented
  AliADQAChecker &operator=(const AliADQAChecker& qac); // not implemented

  const AliADQAParam *fQAParam;  //!
  AliADQAParam  fQAParamDefault; //! used as fallback when OCDB object is not found

  Int_t   fLowEventCut; //! Minimum number of events required by the QA checker
  Float_t fORvsANDCut;  //! AD OR vs AD AND counters cut
  Float_t fBGvsBBCut;   //! AD beam-gas vs beam-beam counters cut

  ClassDef(AliADQAChecker, 9);
} ;

#endif // AliADQAChecker_H
