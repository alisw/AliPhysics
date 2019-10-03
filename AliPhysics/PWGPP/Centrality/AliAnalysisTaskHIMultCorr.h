//-*- Mode: C++ -*-

#ifndef ALIANALYSISTASKHIMULTCORR_H
#define ALIANALYSISTASKHIMULTCORR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Task for HI Multiplicity Correlation checks
// Authors: Jochen Thaeder <jochen@thaeder.de>

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "TTreeStream.h"

class AliMultiplicityCorrelations;
class AliESDEvent;
class AliMCEvent;
class AliKineTrackCuts;
class AliTriggerAnalysis;

class AliAnalysisTaskHIMultCorr : public AliAnalysisTaskSE {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisTaskHIMultCorr(const char *name = "AliAnalysisTaskHIMultCorr");
  virtual ~AliAnalysisTaskHIMultCorr();
  
  /*
   * ---------------------------------------------------------------------------------
   *                                  Task Methods
   * ---------------------------------------------------------------------------------
   */

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  /*
   * ---------------------------------------------------------------------------------
   *                                    Setter
   * ---------------------------------------------------------------------------------
   */
  void SetMaxVertexZ(Float_t vZ) {fMaxVertexZ = vZ;}
  void SetUseCentrality(Int_t cent) {fUseCentralitySel = cent;}
  void SetESDCuts(AliESDtrackCuts* cuts)  {fESDTrackCuts  = cuts;}
  void SetESDCuts2(AliESDtrackCuts* cuts) {fESDTrackCuts2 = cuts;}
  void SetIsMC() { fIsMC = kTRUE; }

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisTaskHIMultCorr(const AliAnalysisTaskHIMultCorr&); // not implemented
  AliAnalysisTaskHIMultCorr& operator=(const AliAnalysisTaskHIMultCorr&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                            Setup Methods - private
   * ---------------------------------------------------------------------------------
   */

  Bool_t SetupEvent();
  Int_t GetCentralityBin();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  TTreeSRedirector *fpcstream;
  
  Bool_t           fIsMC;                        //!
  
  TH1F            *fHStat;                       //! cut statistics

  TList           *fOutList;                     //! output data container

  AliESDEvent     *fESD;                         //! ESD object
  AliESDtrackCuts *fESDTrackCuts;                // ESD cuts  
  AliESDtrackCuts *fESDTrackCuts2;               // ESD cuts 2 
  
  Int_t            fUseCentralitySel;            //  if 0 use none, 1 use VZERO - 2 use SPD
  Int_t            fCentralityBin;               //  current centrality bin

  Int_t            fCentralitySPDBin;            //  SPD centrality bin
  Int_t            fCentralityVZEROBin;          //  VZERO centrality bin
  Float_t          fCentralitySPD;               //  SPD centrality
  Float_t          fCentralityVZERO;             //  VZERO centrality

  Float_t          fMaxVertexZ;                  //  maxVertexZ

  AliTriggerAnalysis          *fTriggerAnalysis; //! trigger analysis object;

  AliMultiplicityCorrelations *fCorrObj;         //! correlations object
  AliMultiplicityCorrelations *fCorrObjCent0;    //! correlations object - centrality  0 -   5
  AliMultiplicityCorrelations *fCorrObjCent1;    //! correlations object - centrality 70 -  80
  AliMultiplicityCorrelations *fCorrObjCent2;    //! correlations object - centrality 80 -  90

  ClassDef(AliAnalysisTaskHIMultCorr, 1);
};

#endif
