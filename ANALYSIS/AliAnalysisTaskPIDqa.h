#ifndef ALIANALYSISTASKPIDQA_H
#define ALIANALYSISTASKPIDQA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskPIDqa.h 43642 2010-09-17 15:50:04Z wiechula $ */
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
class AliVEvent;

class AliAnalysisTaskPIDqa : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDqa();
  AliAnalysisTaskPIDqa(const char *name);
  virtual ~AliAnalysisTaskPIDqa();

  virtual void UserCreateOutputObjects();
  
  virtual void UserExec(Option_t */*option*/);

  
private: 
  AliPIDResponse *fPIDResponse;             //! PID response Handler
  TList                 *fListQA;           //! list with all QA histograms
  TList                 *fListQAits;        //! List with ITS QA histograms
  TList                 *fListQAitsSA;      //! List with ITS SA QA histograms
  TList                 *fListQAitsPureSA;  //! List with ITS pure SA QA histograms
  TList                 *fListQAtpc;        //! List with TPC QA histograms
  TList                 *fListQAtrd;        //! List with TRD QA histograms
  TList                 *fListQAtof;        //! List with TOF QA histograms
  TList                 *fListQAemcal;      //! List with EMCAL QA histograms
  TList                 *fListQAhmpid;      //! List with EMCAL QA histograms
  TList                 *fListQAtofhmpid;   //! List with EMCAL QA histograms
  TList                 *fListQAtpctof;     //! List with combined PID from TPC + TOF

  
  void ExecNewRun();

  //qa object initialisation
  void SetupITSqa();
  void SetupTPCqa();
  void SetupTRDqa();
  void SetupTOFqa();
  void SetupEMCALqa();
  void SetupHMPIDqa();
  void SetupTOFHMPIDqa();
  void SetupTPCTOFqa();

  //
  void FillITSqa();
  void FillTPCqa();
  void FillTRDqa();
  void FillTOFqa();
  void FillEMCALqa();
  void FillHMPIDqa();
  void FillTOFHMPIDqa();
  void FillTPCTOFqa();
  
  //
  void SetRecoInfo();
  
  //helper functions
  TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeArbitraryBinning(const char* bins);
  
  
  AliAnalysisTaskPIDqa(const AliAnalysisTaskPIDqa &other);
  AliAnalysisTaskPIDqa& operator=(const AliAnalysisTaskPIDqa &other);
  
  ClassDef(AliAnalysisTaskPIDqa,1)  // Task to properly set the PID response functions of all detectors
};
#endif
