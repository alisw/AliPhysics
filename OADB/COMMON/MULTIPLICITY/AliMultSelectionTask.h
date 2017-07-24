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
class TProfile;
class TVector3;
class THnSparse;
class TObject;
class TRandom3;

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
class AliMultInput;
class AliMultSelection;
class AliMultSelectionCuts;
class AliOADBMultSelection;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliMultSelectionTask : public AliAnalysisTaskSE {
public:
    
    AliMultSelectionTask();
    AliMultSelectionTask(const char *name, TString lExtraOptions = "", Bool_t lCalib = kFALSE, Int_t lNDebugEstimators = 1);
    virtual ~AliMultSelectionTask();
    
    //Static Event Selection Functions 
    static Bool_t IsSelectedTrigger                    (AliVEvent* event, AliVEvent::EOfflineTriggerTypes lCheckedTrig);
    static Bool_t IsINELgtZERO                         (AliVEvent *event);
    static Bool_t IsAcceptedVertexPosition             (AliVEvent *event);
    static Bool_t IsNotPileupSPD                       (AliVEvent *event);
    static Bool_t IsNotPileupSPDInMultBins             (AliVEvent *event);
    static Bool_t HasNoInconsistentSPDandTrackVertices (AliVEvent *event);
    static Bool_t IsNotAsymmetricInVZERO               (AliVEvent *event);
    static Bool_t IsNotIncompleteDAQ                   (AliVEvent *event);
    static Bool_t HasGoodVertex2016                    (AliVEvent *event);
    
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fkTrigger = trigType;}
    
    //Get Period name (can be static)
    TString GetPeriodNameByLPM(); //try userInfo first
    TString GetPeriodNameByPath( const TString lPath ) const; //no input required, will have all info in globals...
    TString GetPeriodNameByRunNumber()  const; //no input required, use fCurrentRun
    Bool_t CheckOADB( TString lProdName ) const;
    
    //Check MC type
    Bool_t IsHijing() const;
    Bool_t IsDPMJet() const; 
 
    void CreateEmptyOADB(); //In case we really didn't get anything ...
    
    //Cannot be static: requires AliAnalysisUtils Object (why not static?) 
    Bool_t IsNotPileupMV           (AliVEvent *event);
    Bool_t PassesTrackletVsCluster (AliVEvent *event);
    
    //Setup Run if needed (depends on run number!)     
    Int_t SetupRun( const AliVEvent* const esd );
    
    //removed to avoid accidental usage!
    //void SetSaveCalibInfo( Bool_t lVar ) { fkCalibration = lVar; } ;
    void SetAddInfo      ( Bool_t lVar ) { fkAddInfo     = lVar; } ;
    void SetFilterMB     ( Bool_t lVar ) { fkFilterMB    = lVar; } ;
    void SetDebug        ( Bool_t lVar ) { fkDebug       = lVar; } ;
    void SetNDebug       ( Int_t  lVar ) { fNDebug       = lVar; } ;
    void SetHighMultQABinning( Bool_t lVar ) { fkHighMultQABinning = lVar; }
    
    //override for getting estimator definitions from different OADB file
    //FIXME: should preferably be protected, extra functionality required
    void SetAlternateOADBforEstimators      ( TString lFile ){ fAlternateOADBForEstimators      = lFile.Data(); }
    void SetAlternateOADBFullManualBypass   ( TString lFile ){ fAlternateOADBFullManualBypass   = lFile.Data(); }
    void SetAlternateOADBFullManualBypassMC ( TString lFile ){ fAlternateOADBFullManualBypassMC = lFile.Data(); }
    
    //Customize AliMultSelection object name
    void SetStoredObjectName ( TString lObjName ){ fStoredObjectName = lObjName.Data(); }
    
    //Default Setters
    void SetUseDefaultCalib   ( Bool_t lVar ){ fkUseDefaultCalib = lVar; }
    Bool_t GetUseDefaultCalib () const { return fkUseDefaultCalib; }
    
    void SetUseDefaultMCCalib ( Bool_t lVar ){ fkUseDefaultMCCalib = lVar; }
    Bool_t GetUseDefaultMCCalib () const { return fkUseDefaultMCCalib; }
    
    //Calibration mode downscaling for manageable output
    void SetDownscaleFactor ( Double_t lDownscale ) { fDownscaleFactor = lDownscale; }
    
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
    Bool_t fkAttached;    //if true, has already attached to ESD (AOD)
    Bool_t fkHighMultQABinning; //if true, use narrow binning for percentile histograms
    
    //Debug Options
    Bool_t fkDebug;       //if true, saves percentiles in TTree for debugging
    Bool_t fkDebugAliCentrality; //if true, adds V0M percentiles from AliCentrality in TTree
    Bool_t fkDebugAliPPVsMultUtils; //if true, adds V0M percentiles from AliCentrality in TTree
    Bool_t fkDebugIsMC; //if true, adds some MC info for cross-checks (needs MC)
    
    //Default options
    Bool_t fkUseDefaultCalib; //if true, allow for default data calibration
    Bool_t fkUseDefaultMCCalib; //if true, allow for default scaling factor in MC
    
    //Downscale factor:
    //-> if smaller than unity, reduce change of accepting a given event for calib tree
    Double_t fDownscaleFactor;
    TRandom3 *fRand; //PRNG (MT) for random downscaling
    
    //Trigger selection
    AliVEvent::EOfflineTriggerTypes fkTrigger; //kMB, kINT7, etc as needed
    
    TString fAlternateOADBForEstimators;
    TString fAlternateOADBFullManualBypass;
    TString fAlternateOADBFullManualBypassMC;
    
    //Object name for attaching to ESD/AOD
    TString fStoredObjectName;
    
    AliESDtrackCuts *fESDtrackCuts;
    AliAnalysisUtils *fUtils;         // analysis utils
    
    //To perform event selection
    //AliMultSelectionCuts *fMultSelectionCuts; //To perform event selections
    
    //===========================================================================
    //   Variables for Multiplicity Determination
    //===========================================================================
    AliMultVariable *fAmplitude_V0A;
    AliMultVariable *fAmplitude_V0A1;
    AliMultVariable *fAmplitude_V0A2;
    AliMultVariable *fAmplitude_V0A3;
    AliMultVariable *fAmplitude_V0A4;
    AliMultVariable *fAmplitude_V0C;
    AliMultVariable *fAmplitude_V0C1;
    AliMultVariable *fAmplitude_V0C2;
    AliMultVariable *fAmplitude_V0C3;
    AliMultVariable *fAmplitude_V0C4;
    AliMultVariable *fAmplitude_V0Apartial;
    AliMultVariable *fAmplitude_V0Cpartial;
    AliMultVariable *fAmplitude_V0AEq;
    AliMultVariable *fAmplitude_V0CEq;
    AliMultVariable *fAmplitude_OnlineV0A;
    AliMultVariable *fAmplitude_OnlineV0C;
    AliMultVariable *fAmplitude_V0AADC;
    AliMultVariable *fAmplitude_V0CADC;

    //Integer Variables 
    AliMultVariable *fnSPDClusters;
    AliMultVariable *fnSPDClusters0;
    AliMultVariable *fnSPDClusters1;
    AliMultVariable *fnTracklets; //tracklet estimator
    AliMultVariable *fnTracklets08; //tracklet estimator
    AliMultVariable *fnTracklets15; //tracklet estimator
    AliMultVariable *fRefMultEta5; //tracklet estimator
    AliMultVariable *fRefMultEta8; //tracklet estimator
    //AD Related
    AliMultVariable *fMultiplicity_ADA;
    AliMultVariable *fMultiplicity_ADC;
    //ZDC Related
    AliMultVariable *fZncEnergy;
    AliMultVariable *fZpcEnergy;
    AliMultVariable *fZnaEnergy;
    AliMultVariable *fZpaEnergy;
    AliMultVariable *fZem1Energy;
    AliMultVariable *fZem2Energy;
    //ZDC Tower info
    AliMultVariable *fZnaTower;
    AliMultVariable *fZncTower;
    AliMultVariable *fZpaTower;
    AliMultVariable *fZpcTower;
    
    //Event selection snippet for VtxZ as AliMultVariable
    AliMultVariable *fEvSel_VtxZ;

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
    Bool_t fEvSel_IsNotAsymmetricInVZERO;   //!
    Bool_t fEvSel_IsNotIncompleteDAQ;       //!
    Bool_t fEvSel_HasGoodVertex2016;        //!
    //Other Selections: more dedicated filtering to be studied!

    // A.T.
    Int_t   fnContributors; //! 'classical' number of contributors for vertex

    // A.T.
    AliESDtrackCuts* fTrackCuts;        // optional track cuts
    AliESDtrackCuts* fTrackCutsGlobal2015;  // optional track cuts
    AliESDtrackCuts* fTrackCutsITSsa2010; // optional track cuts
    
    AliMultVariable *fZnaFired;
    AliMultVariable *fZncFired;
    AliMultVariable *fZpaFired;
    AliMultVariable *fZpcFired;
    
    AliMultVariable *fNTracks;             //!  no. tracks
    AliMultVariable *fNTracksGlobal2015;             //!  no. tracks (2015 Global track cuts)
    AliMultVariable *fNTracksGlobal2015Trigger;             //!  no. tracks (2015 glob. + TOF-based selection for trigger event)
    AliMultVariable *fNTracksITSsa2010;                     //!  no. tracks ITSsa (2010 ITSsa track cuts)
    
    Int_t fCurrentRun;
    
    Float_t fQuantiles[100]; //! percentiles
    Int_t fEvSelCode; //Final code in event selection 
    Int_t fNDebug; // number of percentiles
    
    Float_t fAliCentralityV0M; //! percentiles from AliCentrality (for debugging) 
    Float_t fPPVsMultUtilsV0M; //! percentiles from AliPPVsMultUtils (for debugging)
    
    //Data needed for Monte Carlo
    Int_t fMC_NColl;
    Int_t fMC_NPart;
    Int_t fMC_NchV0A;
    Int_t fMC_NchV0C;
    Int_t fMC_NchEta05;
    Int_t fMC_NchEta08;
    Int_t fMC_NchEta10;
    
    //Histograms / Anything else as needed
    TH1D *fHistEventCounter; //!
    TH2D *fHistEventSelections; //! For keeping track of 
    
    //Simple QA histograms
    TH1D *fHistQA_V0M;
    TH1D *fHistQA_V0A;
    TH1D *fHistQA_V0C;
    TH1D *fHistQA_CL0; 
    TH1D *fHistQA_CL1;
    TH1D *fHistQA_SPDClusters;
    TH1D *fHistQA_SPDTracklets;
    TH1D *fHistQA_ZNA;
    TH1D *fHistQA_ZNC;
    TH1D *fHistQA_ZNApp;
    TH1D *fHistQA_ZNCpp;
    TProfile *fHistQA_TrackletsVsV0M; 
    TProfile *fHistQA_TrackletsVsCL0; 
    TProfile *fHistQA_TrackletsVsCL1; 
    
    TH1D *fHistQASelected_V0M;
    TH1D *fHistQASelected_V0A;
    TH1D *fHistQASelected_V0C;
    TH1D *fHistQASelected_CL0; 
    TH1D *fHistQASelected_CL1;
    TH1D *fHistQASelected_SPDClusters;
    TH1D *fHistQASelected_SPDTracklets;
    TH1D *fHistQASelected_ZNA;
    TH1D *fHistQASelected_ZNC;
    TH1D *fHistQASelected_ZNApp;
    TH1D *fHistQASelected_ZNCpp;
    TProfile *fHistQASelected_TrackletsVsV0M;
    TProfile *fHistQASelected_TrackletsVsCL0; 
    TProfile *fHistQASelected_TrackletsVsCL1; 

    TProfile *fHistQASelected_NTracksGlobalVsV0M; 
    TProfile *fHistQASelected_NTracksGlobalVsCL0; 
    TProfile *fHistQASelected_NTracksGlobalVsCL1; 
    TProfile *fHistQASelected_PtGlobalVsV0M; 
    TProfile *fHistQASelected_PtGlobalVsCL0; 
    TProfile *fHistQASelected_PtGlobalVsCL1; 

    TProfile *fHistQASelected_NTracksITSsaVsV0M;
    TProfile *fHistQASelected_NTracksITSsaVsCL0;
    TProfile *fHistQASelected_NTracksITSsaVsCL1;
    TProfile *fHistQASelected_PtITSsaVsV0M;
    TProfile *fHistQASelected_PtITSsaVsCL0;
    TProfile *fHistQASelected_PtITSsaVsCL1;

    //AliMultSelection Framework
    AliOADBMultSelection *fOadbMultSelection;
    AliMultInput         *fInput;

    AliMultSelectionTask(const AliMultSelectionTask&);            // not implemented
    AliMultSelectionTask& operator=(const AliMultSelectionTask&); // not implemented

    ClassDef(AliMultSelectionTask, 3);
    //3 - extra QA histograms
};

#endif
