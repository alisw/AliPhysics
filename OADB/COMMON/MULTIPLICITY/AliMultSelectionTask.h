/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
// --- AD(A,C) and fTriggerChargeA(C) added by Tatiana Drozhzhova
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliMultSelectionTask_H
#define AliMultSelectionTask_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
class TObject;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliESDAD;
class AliMultVariable;
class AliMultEstimator;
class AliMultSelection;
class AliMultSelectionCuts;
class AliOADBMultSelection; 


//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliMultSelectionTask : public AliAnalysisTaskSE {
public:
    
    AliMultSelectionTask();
    AliMultSelectionTask(const char *name);
    virtual ~AliMultSelectionTask();
    
    //Static Event Selection Functions 
    static Bool_t IsSelectedTrigger                        (AliVEvent* event, AliVEvent::EOfflineTriggerTypes trigType = AliVEvent::kMB);
    static Bool_t IsINELgtZERO                         (AliVEvent *event);
    static Bool_t IsAcceptedVertexPosition             (AliVEvent *event);
    static Bool_t IsNotPileupSPD                       (AliVEvent *event);
    static Bool_t IsNotPileupSPDInMultBins             (AliVEvent *event);
    static Bool_t HasNoInconsistentSPDandTrackVertices (AliVEvent *event);
    
    //Cannot be static: requires AliAnalysisUtils Object (why not static?) 
    Bool_t IsNotPileupMV (AliVEvent *event);
    Bool_t PassesTrackletVsCluster (AliVEvent *event);
    
    //Setup Run if needed (depends on run number!)     
    Int_t SetupRun( const AliVEvent* const esd );

    void SetSaveCalibInfo( Bool_t lVar ) { fkCalibration = lVar; } ;
    void SetAddInfo      ( Bool_t lVar ) { fkAddInfo     = lVar; } ;
    void SetFilterMB     ( Bool_t lVar ) { fkFilterMB    = lVar; } ;
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //---------------------------------------------------------------------------------------

private:
    
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms
    TTree  *fTreeEvent;              //! Output Tree, Events

    //Objects Controlling Task Behaviour
    Bool_t fkCalibration; //if true, save Calibration object
    Bool_t fkAddInfo;     //if true, save info
    Bool_t fkFilterMB;    //if true, save only kMB events
    Bool_t fkAttached; //if true, has already attached to ESD (AOD)
    
    AliESDtrackCuts *fESDtrackCuts;
    AliAnalysisUtils *fUtils;         // analysis utils
    
    //To perform event selection
    //AliMultSelectionCuts *fMultSelectionCuts; //To perform event selections
    
    //===========================================================================
    //   Variables for Multiplicity Determination
    //===========================================================================
    AliMultVariable *fAmplitude_V0A;
    AliMultVariable *fAmplitude_V0C;
    AliMultVariable *fAmplitude_V0Apartial;
    AliMultVariable *fAmplitude_V0Cpartial;
    AliMultVariable *fAmplitude_V0AEq;
    AliMultVariable *fAmplitude_V0CEq;
    AliMultVariable *fAmplitude_OnlineV0A;
    AliMultVariable *fAmplitude_OnlineV0C;
    //Integer Variables 
    AliMultVariable *fnSPDClusters;
    AliMultVariable *fnTracklets; //tracklet estimator
    AliMultVariable *fRefMultEta5; //tracklet estimator
    AliMultVariable *fRefMultEta8; //tracklet estimator
    //AD Related
    AliMultVariable *fMultiplicity_ADA;
    AliMultVariable *fMultiplicity_ADC;
    
    // A.T.
    Float_t fAmplitude_V0A1;   //!
    Float_t fAmplitude_V0A2;   //!
    Float_t fAmplitude_V0A3;   //!
    Float_t fAmplitude_V0A4;   //!
    Float_t fAmplitude_V0C1;   //!
    Float_t fAmplitude_V0C2;   //!
    Float_t fAmplitude_V0C3;   //!
    Float_t fAmplitude_V0C4;   //!

    Int_t fRunNumber;                    

    //Event Characterization Variables - optional
    Bool_t fEvSel_VtxZCut;                  //!
    Bool_t fEvSel_IsNotPileup;              //!
    Bool_t fEvSel_IsNotPileupMV;            //!
    Bool_t fEvSel_IsNotPileupInMultBins;    //!
    Bool_t fEvSel_Triggered;                //!
    Bool_t fEvSel_INELgtZERO;               //! //done with SPD tracklets
    Bool_t fEvSel_HasNoInconsistentVertices;//!
    Bool_t fEvSel_PassesTrackletVsCluster;  //! 

    //Other Selections: more dedicated filtering to be studied!

    // A.T.
    Int_t   fnSPDClusters0;            //!
    Int_t   fnSPDClusters1;            //!
    Float_t fEvSel_VtxZ; //! the actual value

    // A.T.
    AliESDtrackCuts* fTrackCuts;  //! optional track cuts

    Float_t  fZncEnergy;          //!  ZNC Energy
    Float_t  fZpcEnergy;          //!  ZPC Energy
    Float_t  fZnaEnergy;          //!  ZNA Energy
    Float_t  fZpaEnergy;          //!  ZPA Energy
    Float_t  fZem1Energy;         //!  ZEM1 Energy
    Float_t  fZem2Energy;         //!  ZEM2 Energy
    Double_t fZnaTower;           //! common PMT of ZNA
    Double_t fZncTower;           //! common PMT of ZNC
    Double_t fZpaTower;           //! common PMT of ZPA
    Double_t fZpcTower;           //! common PMT of ZPC
    Bool_t   fZnaFired;
    Bool_t   fZncFired;
    Bool_t   fZpaFired;
    Bool_t   fZpcFired;
    
    Int_t    fNTracks;             //!  no. tracks
    Int_t fCurrentRun; 

    //Histograms / Anything else as needed
    TH1D *fHistEventCounter; //!
    
    //AliMultSelection Framework
    AliOADBMultSelection *oadbMultSelection;
    AliMultSelectionCuts *fMultCuts;
    AliMultSelection     *fSelection;
    AliMultInput         *fInput;

    AliMultSelectionTask(const AliMultSelectionTask&);            // not implemented
    AliMultSelectionTask& operator=(const AliMultSelectionTask&); // not implemented

    ClassDef(AliMultSelectionTask, 1);
};

#endif
