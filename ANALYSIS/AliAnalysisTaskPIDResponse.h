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

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliPIDResponse;
class TList;

class AliAnalysisTaskPIDResponse : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDResponse();
  AliAnalysisTaskPIDResponse(const char *name);
  virtual ~AliAnalysisTaskPIDResponse();

  void SetIsMC(Bool_t isMC=kTRUE) { fIsMC=isMC; }
  
  virtual void UserCreateOutputObjects();
  
  virtual void UserExec(Option_t */*option*/);

  void SetTOFTimeZeroType(Int_t type) { fTOFTimeZeroTypeUser=type; }
  void SetTOFres(Float_t res)         { fTOFres=res;               }
  
private:
  Bool_t fIsMC;                        //  If we run on MC data

  Int_t   fTOFTimeZeroTypeUser;        //  start time type for tof (ESD)
  Int_t   fTOFTimeZeroType;            //! default start time type for tof (ESD)
  Float_t fTOFres;                     //  TOF resolution
  
  AliPIDResponse *fPIDResponse;        //! PID response Handler
  TList                 *fListQA;      //! list with all QA objects
  TList                 *fListQAits;   //! List with ITS QA objects
  TList                 *fListQAtpc;   //! List with TPC QA objects
  TList                 *fListQAtrd;   //! List with TRD QA objects
  TList                 *fListQAtof;   //! List with TOF QA objects

  TString fBeamType;                   //! beam type (PP) or (PBPB)
  TString fLHCperiod;                  //! LHC period
  TString fMCperiodTPC;                //! corresponding MC period to use for the TPC splines
  Int_t   fRecoPass;                   //! reconstruction pass
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  
  TObjArray *fArrPidResponseMaster;    //  TPC pid splines

  void ExecNewRun();

  //qa object initialisation
  void SetupTTSqa();
  void SetupTPCqa();
  void SetupTRDqa();
  void SetupTOFqa();

  //
  void FillITSqa();
  void FillTPCqa();
  void FillTOFqa();
  
  //
  //setup parametrisations
  //
  void SetITSParametrisation();

  //TPC
  void SetTPCPidResponseMaster();
  void SetTPCParametrisation();

  //
  void SetRecoInfo();
  
  //helper functions
  TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeArbitraryBinning(const char* bins);
  
  
  AliAnalysisTaskPIDResponse(const AliAnalysisTaskPIDResponse &other);
  AliAnalysisTaskPIDResponse& operator=(const AliAnalysisTaskPIDResponse &other);
  
  ClassDef(AliAnalysisTaskPIDResponse,1)  // Task to properly set the PID response functions of all detectors
};
#endif
