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
class AliESDv0KineCuts;

class AliAnalysisTaskPIDqa : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDqa();
  AliAnalysisTaskPIDqa(const char *name);
  virtual ~AliAnalysisTaskPIDqa();

  virtual void UserCreateOutputObjects();
  
  virtual void UserExec(Option_t */*option*/);

  
private: 
  AliPIDResponse *fPIDResponse;             //! PID response Handler
  AliESDv0KineCuts *fV0cuts;                //! ESD V0 cuts

  TObjArray *fV0electrons;                  //! array with pointer to identified particles from V0 decays (electrons)
  TObjArray *fV0pions;                      //! array with pointer to identified particles from V0 decays (pions)
  TObjArray *fV0kaons;                      //! array with pointer to identified particles from V0 decays (kaons)
  TObjArray *fV0protons;                    //! array with pointer to identified particles from V0 decays (ptotons)

  TList                 *fListQA;           //! list with all QA histograms
  TList                 *fListQAits;        //! List with ITS QA histograms
  TList                 *fListQAitsSA;      //! List with ITS SA QA histograms
  TList                 *fListQAitsPureSA;  //! List with ITS pure SA QA histograms
  TList                 *fListQAtpc;        //! List with TPC QA histograms
  TList                 *fListQAtpcBasic;   //! Sub-list with TPC QA histograms - basic
  TList                 *fListQAtpcMCtruth; //! Sub-list with TPC QA histograms - only MC truth identified particles
//  TList                 *fListQAtpcHybrid;  //! Sub-list with TPC QA histograms - the "hybrid" scenario -> not used and commented for now
//  TList                 *fListQAtpcOROChigh;//! Sub-list with TPC QA histograms - the "OROChigh" scenario -> not used and commented for now
  TList                 *fListQAtpcV0;      //! Sub-list with TPC QA histograms - V0s
  TList                 *fListQAtrd;        //! List with TRD QA histograms
  TList                 *fListQAtrdNsig;    //! List with TRD QA histograms for Nsigma approach
  TList                 *fListQAtrdNsigTPCTOF; //! List with TRD QA histograms for Nsigma approach after TPC and TOF selection
  TList                 *fListQAtof;        //! List with TOF QA histograms
  TList                 *fListQAt0;         //! List with T0 QA histograms
  TList                 *fListQAemcal;      //! List with EMCAL QA histograms
  TList                 *fListQAhmpid;      //! List with EMCAL QA histograms
  TList                 *fListQAtofhmpid;   //! List with EMCAL QA histograms
  TList                 *fListQAtpctof;     //! List with combined PID from TPC + TOF
  TList                 *fListQAV0;         //! List with V0 kine cuts QA histograms
  TList                 *fListQAinfo;       //! List with information about loaded splines etc.

  
  void ExecNewRun();

  //qa object initialisation
  void SetupITSqa();
  void SetupTPCqa(Bool_t fillMC, Bool_t fill11h, Bool_t fillV0);
  void SetupTRDqa();
  void SetupTOFqa();
  void SetupT0qa();
  void SetupEMCALqa();
  void SetupHMPIDqa();
  void SetupTOFHMPIDqa();
  void SetupTPCTOFqa();
  void SetupV0qa();
  void SetupQAinfo();

  //
  void FillV0PIDlist();
  void ClearV0PIDlist();
  //
  void FillITSqa();
  void FillTPCqa();
  void FillTRDqa();
  void FillTOFqa();
  void FillT0qa();
  void FillEMCALqa();
  void FillHMPIDqa();
  void FillTOFHMPIDqa();
  void FillTPCTOFqa();
  void FillQAinfo();

  // Adding TPC Histograms - called in SetupTPCqa
  void AddTPCHistogramsSignal(TList *sublist, const char *scenario, Int_t scnumber);
  void AddTPCHistogramsNsigma(TList *sublist, const char *scenario, Int_t scnumber);

  // Fill TPC Histograms - called in FillTPCqa
  void FillTPCHistogramsSignal(TList *sublist, Int_t scenario, AliVTrack *track, Int_t nTracks);
  void FillTPCHistogramsNsigma(TList *sublist, Int_t scenario, AliVTrack *track, Int_t nTracks);

  //
  void SetRecoInfo();
  
  //helper functions
  TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  TVectorD* MakeArbitraryBinning(const char* bins);
  
  
  AliAnalysisTaskPIDqa(const AliAnalysisTaskPIDqa &other);
  AliAnalysisTaskPIDqa& operator=(const AliAnalysisTaskPIDqa &other);
  
  ClassDef(AliAnalysisTaskPIDqa,2)  // Task to properly set the PID response functions of all detectors
};
#endif
