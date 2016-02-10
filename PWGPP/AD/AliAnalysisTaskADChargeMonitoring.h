// -*- C++ -*-
// $Id$

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            This task is for AD calibration
//-----------------------------------------------------------------

class TTree;

class ADESDFriendUtils;

#include <TString.h>

#include "AliADRawStream.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskADChargeMonitoring : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskADChargeMonitoring(const char *name="AliAnalysisTaskADChargeMonitoring");
  virtual ~AliAnalysisTaskADChargeMonitoring();

  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

protected:

private:
  // not implemented:
  AliAnalysisTaskADChargeMonitoring(const AliAnalysisTaskADChargeMonitoring&);
  AliAnalysisTaskADChargeMonitoring& operator=(const AliAnalysisTaskADChargeMonitoring&);

  TTree    *fTE;                 //!
  UInt_t    fTimeStamp;          //!
  UShort_t  fBC;                 //!
  ULong64_t fClassMask;          //!
  ULong64_t fClassMaskNext50;    //!
  UInt_t    fPSInt0;             //!
  UInt_t    fPSInt1;             //!
  Float_t   fTriggerCharges[AliADRawStream::kNChannels]; //!
  Bool_t    fBBFlags[AliADRawStream::kNChannels];         //!
  UInt_t    fScalers[AliADRawStream::kNScalers];         //!

  ADESDFriendUtils *fESDADfriendUtils; //!

  ClassDef(AliAnalysisTaskADChargeMonitoring, 2);
} ;
