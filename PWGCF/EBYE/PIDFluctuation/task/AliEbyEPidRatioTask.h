#ifndef ALIEBYEPIDRATIOTASK_H
#define ALIEBYEPIDRATIOTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    // 
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "THnSparse.h"

#include "AliEbyEPidRatioHelper.h"
#include "AliEbyEPidRatioEffCont.h"
#include "AliEbyEPidRatioEffContExtra.h"
#include "AliEbyEPidRatioDCA.h"
#include "AliEbyEPidRatioPhy.h"
#include "AliEbyEPidRatioQA.h"

class AliESDEvent;
class AliESDInputHandler;
class AliAODEvent;
class AliAODInputHandler;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDResponse;

class AliEbyEPidRatioTask : public AliAnalysisTaskSE {

 public:

  AliEbyEPidRatioTask(const char *name = "AliEbyEPidRatioTask");
  virtual ~AliEbyEPidRatioTask();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  /** Initialize Task */
  Int_t Initialize();

  void SetIsMC()                             {fIsMC         = kTRUE;}
  void SetIsRatio(Int_t i);
  void SetIsPer()                            {fIsPer        = kTRUE;}
  void SetEffExtra()                         {fIsEffExtra   = kTRUE;}  

  void SetIsAOD(Bool_t b          )          {fIsAOD            = b;}
  void SetESDTrackCutMode(Int_t i )          {fESDTrackCutMode  = i;}
  void SetModeEffCreation(Int_t i )          {fModeEffCreation  = i;}
  void SetModeDCACreation(Int_t i )          {fModeDCACreation  = i;}
  void SetModeDistCreation(Int_t i)          {fModeDistCreation = i;}

  void SetModeQACreation(Int_t i  )          {fModeQACreation   = i;}
 
  void SetEtaMax(Float_t f        )          {fEtaMax           = f;}
  void SetEtaMaxEff(Float_t f     )          {fEtaMaxEff        = f;}
  void SetPtRange(Float_t f1, Float_t f2)    {fPtRange[0]       = f1; fPtRange[1]    = f2;}
  void SetPtRangeEff(Float_t f1, Float_t f2) {fPtRangeEff[0]    = f1; fPtRangeEff[1] = f2;}

  void SetTrackFilterBit(Int_t i)            {fAODtrackCutBit   = i;}
  
  void SetNetParticleHelper(AliEbyEPidRatioHelper *helper) {
    if (fHelper) 
      delete fHelper;
    fHelper = helper;
  }

 private:

  AliEbyEPidRatioTask(const AliEbyEPidRatioTask&); // not implemented
  AliEbyEPidRatioTask& operator=(const AliEbyEPidRatioTask&); // not implemented

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

  AliEbyEPidRatioHelper       *fHelper;  //  Helper class
  AliEbyEPidRatioEffCont      *fEffCont; //! Efficiency and Contamination class
  AliEbyEPidRatioEffContExtra *fEffContExtra; //! Efficiency and Contamination class
  AliEbyEPidRatioDCA          *fDCA;     //! DCA class
  AliEbyEPidRatioPhy          *fDist;    //! Distributions class
  AliEbyEPidRatioQA           *fQA;      //! QA class

  TList              *fOutList;                 //! Ptr to output data container
  TList              *fOutListEff;              //! Ptr to output data container - efficiency
  TList              *fOutListCont;             //! Ptr to output data container - contamination
  TList              *fOutListDCA;              //! Ptr to output data container - DCA
  TList              *fOutListQA;               //! Ptr to output data container - QA

  AliESDEvent        *fESD;                     //! Ptr to ESD event
  AliESDInputHandler *fESDHandler;              //! Ptr to ESD Handler
 
  AliESDtrackCuts    *fESDTrackCutsBase;        //! ESD cuts - base settings
  AliESDtrackCuts    *fESDTrackCuts;            //! ESD cuts  
  AliESDtrackCuts    *fESDTrackCutsBkg;         //! ESD cuts for Bkg
  AliESDtrackCuts    *fESDTrackCutsEff;         //! ESD cuts for efficiency determination -> larger pt Range

  AliAODEvent        *fAOD;                     //! Ptr to AOD event
  AliAODInputHandler *fAODHandler;              //! Ptr to AOD Handler
  
  Bool_t              fIsMC;                    //  Is MC event 
  Bool_t              fIsRatio;                 //  Is MC event b
  Bool_t              fIsPtBin;                 //  Is MC Event a 
  Bool_t              fIsDetectorWise;          // Is detectorwise
  Bool_t              fIsAOD;                   //  analysis mode            : 0 = ESDs  | 1 = AODs
  Bool_t              fIsSub;                   //  analysis mode      SS  
  Bool_t              fIsBS;                    //  analysis mode      BS
  Bool_t              fIsPer;                   //  analysis mode      PER
  Bool_t              fIsEffExtra; 
  Int_t               fESDTrackCutMode;         //  ESD track cut mode       : 0 = clean | 1 = dirty
  Int_t               fModeEffCreation ;        //  Correction creation mode : 1 = on    | 0 = off
  Int_t               fModeDCACreation;         //  DCA creation mode        : 1 = on    | 0 = off
  Int_t               fModeDistCreation;        //  Dist creation mode       : 1 = on    | 0 = off
  Int_t               fModeQACreation;          //  QA creation mode         : 1 = on    | 0 = off
  AliMCEvent         *fMCEvent;                 //! Ptr to MC event
  AliStack           *fMCStack;                 //! Ptr to MC stack
  
  Float_t             fEtaMax;                  //  Max, absolut eta 
  Float_t             fEtaMaxEff;               //  Max, absolut eta for efficiency
  Float_t             fPtRange[2];              //  Array of pt [min,max]
  Float_t             fPtRangeEff[2];           //  Array of pt [min,max] for efficiency
  
  Int_t               fAODtrackCutBit;          //  Track filter bit for AOD tracks
    
  ClassDef(AliEbyEPidRatioTask, 1);
};

#endif
