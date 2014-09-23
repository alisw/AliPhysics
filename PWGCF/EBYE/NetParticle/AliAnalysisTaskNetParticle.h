//-*- Mode: C++ -*-

#ifndef ALIANALYSISTASKNETPARTICLE_H
#define ALIANALYSISTASKNETPARTICLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/**
 * Class for NetParticle Distributions
 * -- AnalysisTask
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "THnSparse.h"

#include "AliAnalysisNetParticleHelper.h"
#include "AliAnalysisNetParticleEffCont.h"
#include "AliAnalysisNetParticleDCA.h"
#include "AliAnalysisNetParticleDistribution.h"
#include "AliAnalysisNetParticleQA.h"

class AliESDEvent;
class AliESDInputHandler;
class AliAODEvent;
class AliAODInputHandler;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDResponse;

class AliAnalysisTaskNetParticle : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskNetParticle(const char *name = "AliAnalysisTaskNetParticle");
  virtual ~AliAnalysisTaskNetParticle();
  
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
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize Task */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Setter
   * ---------------------------------------------------------------------------------
   */

  void SetIsMC()                             {fIsMC             = kTRUE;}
  void SetIsAOD(Bool_t b)                    {fIsAOD            = b;}
  
  void SetESDTrackCutMode(Int_t i)           {fESDTrackCutMode  = i;}
  void SetModeEffCreation(Int_t i)           {fModeEffCreation  = i;}
  void SetModeDCACreation(Int_t i)           {fModeDCACreation  = i;}
  void SetModeDistCreation(Int_t i)          {fModeDistCreation = i;}
  void SetModeQACreation(Int_t i)            {fModeQACreation   = i;}
 
  void SetEtaMax(Float_t f)                  {fEtaMax           = f;}
  void SetEtaMaxEff(Float_t f)               {fEtaMaxEff        = f;}
  void SetPtRange(Float_t f1, Float_t f2)    {fPtRange[0]       = f1; fPtRange[1]    = f2;}
  void SetPtRangeEff(Float_t f1, Float_t f2) {fPtRangeEff[0]    = f1; fPtRangeEff[1] = f2;}

  void SetTrackFilterBit(Int_t i)            {fAODtrackCutBit   = i;}
  
  void SetNetParticleHelper(AliAnalysisNetParticleHelper *helper) {
    if (fHelper) 
      delete fHelper;
    fHelper = helper;
  }

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisTaskNetParticle(const AliAnalysisTaskNetParticle&); // not implemented
  AliAnalysisTaskNetParticle& operator=(const AliAnalysisTaskNetParticle&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                            Setup/Reset Methods - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Setup Event */
  Int_t SetupEvent();

  /** Setup ESD Event */
  Int_t SetupESDEvent();

  /** Setup AOD Event */
  Int_t SetupAODEvent();

  /** Setup MC Event */
  Int_t SetupMCEvent();

  /** Reset Event */
  void ResetEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  AliAnalysisNetParticleHelper       *fHelper;  //  Helper class
  AliAnalysisNetParticleEffCont      *fEffCont; //! Efficiency and Contamination class
  AliAnalysisNetParticleDCA          *fDCA;     //! DCA class
  AliAnalysisNetParticleDistribution *fDist;    //! Distributions class
  AliAnalysisNetParticleQA           *fQA;      //! QA class

  // --- OutLists ----------------------------------------------------------
  TList              *fOutList;                 //! Ptr to output data container
  TList              *fOutListEff;              //! Ptr to output data container - efficiency
  TList              *fOutListCont;             //! Ptr to output data container - contamination
  TList              *fOutListDCA;              //! Ptr to output data container - DCA
  TList              *fOutListQA;               //! Ptr to output data container - QA

  // --- ESD only ----------------------------------------------------------
  AliESDEvent        *fESD;                     //! Ptr to ESD event
  AliESDInputHandler *fESDHandler;              //! Ptr to ESD Handler
  // -----------------------------------------------------------------------
  AliESDtrackCuts    *fESDTrackCutsBase;        //! ESD cuts - base settings
  AliESDtrackCuts    *fESDTrackCuts;            //! ESD cuts  
  AliESDtrackCuts    *fESDTrackCutsBkg;         //! ESD cuts for Bkg
  AliESDtrackCuts    *fESDTrackCutsEff;         //! ESD cuts for efficiency determination -> larger pt Range

  // --- AOD only ----------------------------------------------------------
  AliAODEvent        *fAOD;                     //! Ptr to AOD event
  AliAODInputHandler *fAODHandler;              //! Ptr to AOD Handler
  
  // --- Flags -------------------------------------------------------------
  Bool_t              fIsMC;                    //  Is MC event
  Bool_t              fIsAOD;                   //  analysis mode            : 0 = ESDs  | 1 = AODs
  Int_t               fESDTrackCutMode;         //  ESD track cut mode       : 0 = clean | 1 = dirty
  Int_t               fModeEffCreation ;        //  Correction creation mode : 1 = on    | 0 = off
  Int_t               fModeDCACreation;         //  DCA creation mode        : 1 = on    | 0 = off
  Int_t               fModeDistCreation;        //  Dist creation mode       : 1 = on    | 0 = off
  Int_t               fModeQACreation;          //  QA creation mode         : 1 = on    | 0 = off

  // --- MC only -----------------------------------------------------------
  AliMCEvent         *fMCEvent;                 //! Ptr to MC event
  AliStack           *fMCStack;                 //! Ptr to MC stack
  // -----------------------------------------------------------------------
  Float_t             fEtaMax;                  //  Max, absolut eta 
  Float_t             fEtaMaxEff;               //  Max, absolut eta for efficiency
  Float_t             fPtRange[2];              //  Array of pt [min,max]
  Float_t             fPtRangeEff[2];           //  Array of pt [min,max] for efficiency
  // -----------------------------------------------------------------------
  Int_t               fAODtrackCutBit;          //  Track filter bit for AOD tracks
  // =======================================================================

  
  ClassDef(AliAnalysisTaskNetParticle, 1);
};

#endif
