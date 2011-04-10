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
class AliVEvent;

class AliAnalysisTaskPIDResponse : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDResponse();
  AliAnalysisTaskPIDResponse(const char *name);
  virtual ~AliAnalysisTaskPIDResponse();

  void SetIsMC(Bool_t isMC=kTRUE) { fIsMC=isMC; }
  
  virtual void UserCreateOutputObjects();
  
  virtual void UserExec(Option_t */*option*/);

  
private:
  Bool_t fIsMC;                        //  If we run on MC data
  
  AliPIDResponse *fPIDResponse;        //! PID response Handler
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  Int_t   fRecoPass;                   //! reconstruction pass
  
  //
  void SetRecoInfo();
    
  AliAnalysisTaskPIDResponse(const AliAnalysisTaskPIDResponse &other);
  AliAnalysisTaskPIDResponse& operator=(const AliAnalysisTaskPIDResponse &other);
  
  ClassDef(AliAnalysisTaskPIDResponse,1)  // Task to properly set the PID response functions of all detectors
};
#endif
