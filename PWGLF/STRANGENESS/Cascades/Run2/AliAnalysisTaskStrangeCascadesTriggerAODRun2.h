/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//
// Used bits of code from AliAnalysisTaskStrangenessVsMultiplicityRun2
//					 and  AliAnalysisTaskTPCTOFCascade
//
// --- Romain Schotter
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskStrangeCascadesTriggerAODRun2_H
#define AliAnalysisTaskStrangeCascadesTriggerAODRun2_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCascadeResult.h"
#include "AliV0Result.h"
#include "AliESDtrackCuts.h"

#include "AliEventContainer.h"
#include "AliCascadeContainer.h"
#include "AliV0Container.h"
#include "AliTrackInfoContainer.h"

class AliAnalysisTaskStrangeCascadesTriggerAODRun2 : public AliAnalysisTaskSE  {
public:
                            AliAnalysisTaskStrangeCascadesTriggerAODRun2();
                            AliAnalysisTaskStrangeCascadesTriggerAODRun2(const char *name, Bool_t lSaveCascades, Bool_t lSaveV0s, Bool_t lSaveTracks);
    virtual                 ~AliAnalysisTaskStrangeCascadesTriggerAODRun2();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
	
    //---------------------------------------------------------------------------------------
    //For cascade
    void        AddConfiguration( AliCascadeResult *lCascadeResult );
    //---------------------------------------------------------------------------------------
    //For V0s
    void        AddConfiguration( AliV0Result      *lV0Result      );
    //---------------------------------------------------------------------------------------
    //For Tracks
    void        AddConfiguration( AliESDtrackCuts  *ltrackCuts     );
    //---------------------------------------------------------------------------------------    
    void        SetSaveCascades             (Bool_t use = kTRUE)        { fSaveCascades = use ;}
    void        SetSaveV0s                  (Bool_t use = kTRUE)        { fSaveV0s = use ;}
    void        SetSaveTracks            	(Bool_t use = kTRUE)        { fSaveTracks = use ;}
    void        SetTrigOnCascade            (Bool_t use = kTRUE)        { fTrigOnCascade = use ;}
    
    void        SetMaxPrimVtxR2D            (Double_t val)              { fkMaxPVR2D = val ;}
    void        SetMaxPrimVtxZ              (Double_t val)              { fkMaxPVZ = val ;}
    
    void        SetSaveAddCascConfig        (Bool_t use = kTRUE)        { fCascSaveAddConfig = use ;}
    void        SetSaveAddV0Config          (Bool_t use = kTRUE)        { fV0SaveAddConfig = use ;}
    void        SetSaveAddTrackConfig   	(Bool_t use = kTRUE)        { fTrackSaveAddConfig = use ;}
    
    void        SetMinNbrCrossedRows        (Int_t val)                 { fMinNbrCrossedRows = val;}
    void        SetMinPtToSave              (Double_t val)              { fMinPtToSave = val;}
    void        SetMaxPtToSave              (Double_t val)              { fMaxPtToSave = val;}
    void        SetMaxAbsEta                (Double_t val)              { fMaxAbsEta = val;}
    void        SetMaxAbsRap                (Double_t val)              { fMaxAbsRap = val;}
    void		SetAODFilterBit				(UInt_t val)				{ fAODFilterBit	= 1 << val ;}
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TTree                       *fOutputTree;           //! Output Tree, containing everything in the event
    TClonesArray                *fArrayCascade;         //! Array of Cascade containers
    TClonesArray                *fArrayV0;              //! Array of V0 containers
    TClonesArray                *fArrayTrack;           //! Array of Track containers
    
    AliEventContainer           *fAnalysisEvent;        //!
    AliCascadeContainer         *fAnalysisCascade;      //! Dummy container of Cascade ; only used for cuts
    AliV0Container              *fAnalysisV0;           //! Dummy container of V0      ; only used for cuts
    AliTrackInfoContainer       *fAnalysisTrack;        //! Dummy container of Track   ; only used for cuts
    
    AliPIDResponse      *fPIDResponse;          //! PID response object
    AliAnalysisUtils    *fUtils;                //! analysis utils (for MV pileup selection)
    
    //---> Flags controlling TTree outputs
    Bool_t      fSaveCascades;
    Bool_t      fSaveV0s;
    Bool_t      fSaveTracks;
    Bool_t      fTrigOnCascade;// default = kTRUE 
    //---> Variables controlling PV selections
    Double_t    fkMaxPVR2D;
    Double_t    fkMaxPVZ;
    //---> Variables controlling cascade and V0 default selections
    Int_t       fMinNbrCrossedRows;
    Double_t    fMinPtToSave;
    Double_t    fMaxPtToSave;
    Double_t    fMaxAbsEta;
    Double_t    fMaxAbsRap;
	UInt_t		fAODFilterBit;
    //---> Flags controlling cascade and V0 custom selections
    Bool_t          fCascSaveAddConfig;
    Bool_t          fV0SaveAddConfig;
    Bool_t          fTrackSaveAddConfig;
    TObjArray       fCascadeResult;                // list of cascade cuts
    TObjArray       fV0Result;                     // list of V0 cuts
    AliESDtrackCuts *fTrackCuts;               	   // Primary track selections
    
    AliESDtrackCuts* fTrackCuts2010;
    AliESDtrackCuts* fTrackCuts2011;
    AliESDtrackCuts* fTrackCutsTPCRefit;
    AliESDtrackCuts* fTrackCutsV0;
    AliESDtrackCuts* fTrackCuts2011Sys;
    AliESDtrackCuts* fTrackCutsHybrid_kNone;
    AliESDtrackCuts* fTrackCutsHybrid_kOff;


//===========================================================================================
//   Histograms
//===========================================================================================
    TList*                  fOutputList;    //! output list
    TH2D*                   fHistEventCounter;        //! dummy histogram

    AliAnalysisTaskStrangeCascadesTriggerAODRun2(const AliAnalysisTaskStrangeCascadesTriggerAODRun2&); // not implemented
    AliAnalysisTaskStrangeCascadesTriggerAODRun2& operator=(const AliAnalysisTaskStrangeCascadesTriggerAODRun2&); // not implemented

    ClassDef(AliAnalysisTaskStrangeCascadesTriggerAODRun2, 3);
};

//_____________________________________________________________________________
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2::AddConfiguration( AliCascadeResult *lCascadeResult )
{
    //
    // Add a new cut set for a new criterion for cascade selection.
    // A user can add as many as he wants, and each one corresponds either
    // to a Xi or a Omega.
    // An extra security is applied, to be sure we are working with a Xi or a Omega
    // A little paranoid...

    if 
    (    lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiMinus     || 
         lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiPlus      ||
         lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaMinus  ||
         lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaPlus   
    ){
        fCascadeResult.AddLast(lCascadeResult);
    }
    else
    {
        AliInfo(Form("AliCascadeResult '%s' does not correspond to a Xi or Omega", lCascadeResult->GetName()));
        return;
    }

}

//_____________________________________________________________________________
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2::AddConfiguration( AliV0Result *lV0Result )
{
    //
    // Add a new cut set for a new criterion for V0 selection.
    // A user can add as many as he wants, and each one corresponds either
    // to a K0s or a Lambda or AntiLambda.
    // An extra security is applied, to be sure we are working with aa K0s or a Lambda or AntiLambda
    // A little paranoid...

    if 
    (    lV0Result->GetMassHypothesis() == AliV0Result::kK0Short     || 
         lV0Result->GetMassHypothesis() == AliV0Result::kLambda      ||
         lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda  
    ){
        fV0Result.AddLast(lV0Result);
    }
    else
    {
        AliInfo(Form("AliV0Result '%s' does not correspond to a K0s or Lambda or AntiLambda", lV0Result->GetName()));
        return;
    }

}

//_____________________________________________________________________________
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2::AddConfiguration( AliESDtrackCuts *ltrackCuts )
{
    //
    // Define a cut set for primary tracks selection.
    //
    fTrackCuts = ltrackCuts;
}

#endif
