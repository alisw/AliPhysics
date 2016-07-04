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
  Bool_t GetTunedOnData() const { return fIsTunedOnData; };
  void SetTuneOnDataMask(Int_t mask){fTunedOnDataMask=mask;};
  
  void SetUseTPCEtaCorrection(Bool_t useTPCEtaCorrection) { fUseTPCEtaCorrection = useTPCEtaCorrection; };
  Bool_t UseTPCEtaCorrection() const { return fUseTPCEtaCorrection; };
  
  void SetUseTPCMultiplicityCorrection(Bool_t useMultiplicityCorrection = kTRUE) { fUseTPCMultiplicityCorrection = useMultiplicityCorrection; };
  Bool_t UseTPCMultiplicityCorrection() const { return fUseTPCMultiplicityCorrection; };


  void SetUseTRDEtaCorrection(Bool_t useTRDEtaCorrection) { fUseTRDEtaCorrection = useTRDEtaCorrection; };
  Bool_t UseTRDEtaCorrection() const { return fUseTRDEtaCorrection; };
  void SetUseTRDClusterCorrection(Bool_t useTRDClusterCorrection) { fUseTRDClusterCorrection = useTRDClusterCorrection; };
  Bool_t UseTRDClusterCorrection() const { return fUseTRDClusterCorrection; };


  void SetSpecialDetectorResponse(const char* det) { fSpecialDetResponse=det; }
  void SetUserDataRecoPass(Int_t pass){fUserDataRecoPass=pass;};


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
  
  Bool_t fUseTPCEtaCorrection;          // Use TPC eta correction
  Bool_t fUseTPCMultiplicityCorrection; // Use TPC multiplicity correction
  Bool_t fUseTRDEtaCorrection;          // Use TRD eta correction
  Bool_t fUseTRDClusterCorrection;      // Use TRD cluster correction

  Int_t  fUserDataRecoPass;            // forced DATA reco pass
  
  //
  void SetRecoInfo();
    
  AliAnalysisTaskPIDResponse(const AliAnalysisTaskPIDResponse &other);
  AliAnalysisTaskPIDResponse& operator=(const AliAnalysisTaskPIDResponse &other);
  
  ClassDef(AliAnalysisTaskPIDResponse,9)  // Task to properly set the PID response functions of all detectors
};
#endif
