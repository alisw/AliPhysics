#ifndef ALIANALYSISTASKPIDRESPONSE_H
#define ALIANALYSISTASKPIDRESPONSE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskPIDResponse.h 43642 2010-09-17 15:50:04Z wiechula $ */
// Author: Jens Wiechula, 24/02/2011

//==============================================================================
//
//
//
//
//==============================================================================

#include <TVectorDfwd.h>
#include <TString.h>

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliPIDResponse;
class AliVEvent;

class AliAnalysisTaskPIDResponse : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDResponse();
  AliAnalysisTaskPIDResponse(const char *name);
  virtual ~AliAnalysisTaskPIDResponse();

  void SetIsMC(Bool_t isMC=kTRUE)   { fIsMC=isMC; }
  void SetCachePID(Bool_t cachePID) { fCachePID=cachePID; }
  Bool_t GetCachePID() const { return fCachePID; }
  
  virtual void UserCreateOutputObjects();
  
  virtual void UserExec(Option_t */*option*/);

  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char* GetOADBPath() const { return fOADBPath.Data(); }
  void SetTuneOnData(Bool_t flag,Int_t recopass){fIsTunedOnData=flag;fRecoPassTuned=recopass;};
  void SetTuneOnDataMask(Int_t mask){fTunedOnDataMask=mask;};
  
  void SetUseTPCEtaCorrection(Bool_t useTPCEtaCorrection) { fUseTPCEtaCorrection = useTPCEtaCorrection; };
  Bool_t UseTPCEtaCorrection() const { return fUseTPCEtaCorrection; };

  void SetSpecialDetectorResponse(const char* det) { fSpecialDetResponse=det; }

private:
  Bool_t fIsMC;                        // If we run on MC data
  Bool_t fCachePID;                    // Cache PID values in transient object
  TString fOADBPath;                   // OADB path to use
  TString fSpecialDetResponse;         // Special detector response files for debugging
  
  AliPIDResponse *fPIDResponse;        //! PID response Handler
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  Int_t   fRecoPass;                   //! reconstruction pass

  Bool_t  fIsTunedOnData;              // flag to tune MC on data (dE/dx)
  Int_t   fTunedOnDataMask;            // mask to activate tuning on data on specific detectors
  Int_t   fRecoPassTuned;              // Reco pass tuned on data for MC
  
  Bool_t  fUseTPCEtaCorrection;        // Use TPC eta correction
  
  //
  void SetRecoInfo();
    
  AliAnalysisTaskPIDResponse(const AliAnalysisTaskPIDResponse &other);
  AliAnalysisTaskPIDResponse& operator=(const AliAnalysisTaskPIDResponse &other);
  
  ClassDef(AliAnalysisTaskPIDResponse,6)  // Task to properly set the PID response functions of all detectors
};
#endif
