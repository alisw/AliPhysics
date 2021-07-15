/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskStrangeCascadesTriggerAODRun2MC_H
#define AliAnalysisTaskStrangeCascadesTriggerAODRun2MC_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliRsnEvent.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"
#include "AliRsnMiniResonanceFinder.h"
#include "AliRsnCutSet.h"
#include "AliRsnMiniPair.h"
#include "AliCascadeResult.h"
#include "AliV0Result.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskStrangeCascadesTriggerAODRun2MC : public AliAnalysisTaskSE  {
public:
                            AliAnalysisTaskStrangeCascadesTriggerAODRun2MC();
                            AliAnalysisTaskStrangeCascadesTriggerAODRun2MC(const char *name, Bool_t lSaveV0s, Bool_t lSaveRsn, Bool_t lSavePrimaries, Bool_t lUseEventMixing, Bool_t lSaveMCGenPart);
    virtual                 ~AliAnalysisTaskStrangeCascadesTriggerAODRun2MC();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void            FinishTaskOutput();
    
    //---------------------------------------------------------------------------------------
    Float_t     GetDCAz(AliAODTrack *lTrack);
    Float_t     GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );
    //---------------------------------------------------------------------------------------
    Double_t    MyPhi(Double_t lPx, Double_t lPy) const;
    Double_t    MyTheta(Double_t lPt, Double_t lPz) const;
    //---------------------------------------------------------------------------------------
    Bool_t      CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2);
    //---------------------------------------------------------------------------------------
    //For cascade
    void        AddConfiguration( AliCascadeResult *lCascadeResult );
    //---------------------------------------------------------------------------------------
    //For V0s
    void        AddConfiguration( AliV0Result      *lV0Result      );
    //---------------------------------------------------------------------------------------
    //For Primaries
    void        AddConfiguration( AliESDtrackCuts  *ltrackCuts     );
    //---------------------------------------------------------------------------------------
    //For resonances
    Int_t       AddTrackCuts(AliRsnCutSet *cuts);
    Int_t       AddResonanceFinder(AliRsnMiniResonanceFinder *f);
    Int_t       GetNumberOfTrackCuts() { return fRsnTrackCuts.GetEntries(); }
    
    Bool_t      EventsMatch(AliRsnMiniEvent* event1, AliRsnMiniEvent *event2);
    Int_t       FillPair(AliRsnMiniEvent *miniEventRef, AliRsnMiniEvent *fRsnMiniEvent, TObjArray &lPairList, Bool_t refFirst);
    Double_t    ComputeSpherocity();
    Double_t    ComputeAngle();
    //---------------------------------------------------------------------------------------    
    void        SetSaveV0s                  (Bool_t use = kTRUE)        { fSaveV0s = use ;}
    void        SetSaveResonances           (Bool_t use = kTRUE)        { fSaveRsn = use ;}
    void        SetSavePrimaries            (Bool_t use = kTRUE)        { fSavePrimaries = use ;}
    void        SetUseEventMixing           (Bool_t use = kFALSE)       { fUseEventMixing = use ;}
    void        SetSaveMCGenPart            (Bool_t use = kFALSE)       { fSaveMCGenPart = use ;}
    void        SetTrigOnCascade            (Bool_t use = kTRUE)        { fTrigOnCascade = use ;}
    
    void        SetNMix                     (Int_t val)                 { fNMix = val;}
    void        SetMaxDiffVz                (Double_t val)              { fMaxDiffVz = val;}
    void        SetMaxDiffMult              (Double_t val)              { fMaxDiffMult = val;}
    void        SetMaxDiffAngle             (Double_t val)              { fMaxDiffAngle = val;}
    void        SetPairCuts                 (AliRsnCutSet *set)         { fRsnPairCuts = set;}
    
    void        SetMaxPrimVtxR2D            (Double_t val)              { fkMaxPVR2D = val ;}
    void        SetMaxPrimVtxZ              (Double_t val)              { fkMaxPVZ = val ;}
    
    void        SetSaveAddCascConfig        (Bool_t use = kTRUE)        { fCascSaveAddConfig = use ;}
    void        SetSaveAddV0Config          (Bool_t use = kTRUE)        { fV0SaveAddConfig = use ;}
    void        SetSaveAddPrimTrackConfig   (Bool_t use = kTRUE)        { fPrimariesSaveAddConfig = use ;}
    
    void        SetMinNbrCrossedRows        (Int_t val)                 { fMinNbrCrossedRows = val;}
    void        SetMinPtToSave              (Double_t val)              { fMinPtToSave = val;}
    void        SetMaxPtToSave              (Double_t val)              { fMaxPtToSave = val;}
    void        SetMaxAbsEta                (Double_t val)              { fMaxAbsEta = val;}
    void        SetMaxAbsRap                (Double_t val)              { fMaxAbsRap = val;}
    void        SetRejectCascKinkDaughters  (Bool_t use = kTRUE)        { fRejectCascKink = use ;}
    void        SetRejectV0KinkDaughters    (Bool_t use = kTRUE)        { fRejectV0Kink = use ;}
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TTree  *fTreeCascade;              //! Output Tree, Cascades
    TTree  *fTreeV0;                   //! Output Tree, V0s
    TTree  *fTreeRsn;                  //! Output Tree, Resonances
    TTree  *fTreeRsnMixed;             //! Output Tree, Mixed-event resonances
    TTree  *fTreePrimTrack;            //! Output Tree, Primary tracks
    TTree  *fTreeMCParticle;           //! Output Tree, MC particle
    
    TTree  *fDummyTree;                //! Dummy output Tree, Needed for event mixing on resonances
    
    ULong64_t fDummyVarEventNumber;    //!
    
    AliPIDResponse      *fPIDResponse;          //! PID response object
    AliAnalysisUtils    *fUtils;                //! analysis utils (for MV pileup selection)
    
    AliRsnEvent         *fRsnEvent;             //! Rsn Event object
    AliRsnMiniEvent     *fRsnMiniEvent;         //! Rsn Mini-Event object
    TObjArray           fResonanceFinders;      // list of AliRsnMiniResonanceFinder objects
    TObjArray           fRsnTrackCuts;          // list of single track cuts
    AliRsnCutSet        *fRsnPairCuts;          // cuts on the pair
    
    TString     fPairName;          // Name of the pair (or the resonance finder)
    Double_t    fMotherMass;        // nominal resonance mass
    Char_t      fCharge[2];         // required track charge
    RSNPID      fDaughter[2];       // species of resonance daughters
    Int_t       fCutID[2];          // ID of cut set used to select tracks
    
    //---> Flags controlling TTree outputs
    Bool_t      fSaveV0s;
    Bool_t      fSaveRsn;
    Bool_t      fSavePrimaries;
    Bool_t      fUseEventMixing; // kTRUE = Perform event mixing
    Bool_t      fSaveMCGenPart;
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
    Bool_t      fRejectCascKink;
    Bool_t      fRejectV0Kink;
    //---> Flags controlling cascade and V0 custom selections
    Bool_t          fCascSaveAddConfig;
    Bool_t          fV0SaveAddConfig;
    Bool_t          fPrimariesSaveAddConfig;
    TObjArray       fCascadeResult;                // list of cascade cuts
    TObjArray       fV0Result;                     // list of V0 cuts
    AliESDtrackCuts *fPrimTrackCuts;               // Primary track selections
    //---> Variables controlling event mixing
    Int_t       fNMix;          // Nbr of events to use for event mixing
    Double_t    fMaxDiffVz;     // Max diff of the PVz between two events (event mixing)
    Double_t    fMaxDiffMult;   // Max diff in multiplity (percentiles) between two events (event mixing)
    Double_t    fMaxDiffAngle;  // Max diff in plane angle between two events (event mixing)

//===========================================================================================
//   Variables for Cascade Tree
//===========================================================================================
    
    //-----------BASIC-INFO---------------------------
    Int_t fTreeCascVarCharge;//!
    Double_t fTreeCascVarMassAsXi;//!
    Double_t fTreeCascVarMassAsOmega;//!
    Double_t fTreeCascVarPtot;//!
    Double_t fTreeCascVarV0Ptot;//!
    Double_t fTreeCascVarPt;//!
    Double_t fTreeCascVarV0Pt;//!
    Double_t fTreeCascVarRapXi;//!
    Double_t fTreeCascVarRapOmega;//!
    Double_t fTreeCascVarBachEta;//!
    Double_t fTreeCascVarPosEta;//!
    Double_t fTreeCascVarNegEta;//!
    Double_t fTreeCascVarPhi;//!
    Double_t fTreeCascVarTheta;//!
    //-----------INFO-FOR-CUTS------------------------
    Double_t fTreeCascVarAlpha;//!
    Double_t fTreeCascVarPtArm;//!
    Double_t fTreeCascVarAlphaV0;//!
    Double_t fTreeCascVarPtArmV0;//!
    Double_t fTreeCascVarDCACascDau;//!
    Double_t fTreeCascVarDCABachToPV;//!
    Double_t fTreeCascVarDCAV0Dau;//!
    Double_t fTreeCascVarDCAV0ToPV;//!
    Double_t fTreeCascVarDCAPosToPV;//!
    Double_t fTreeCascVarDCANegToPV;//!
    Double_t fTreeCascVarCascCosPA;//!
    Double_t fTreeCascVarCascCosPASpecial;//!
    
    Double_t fTreeCascVarCascRadius;//!
    Double_t fTreeCascVarV0Mass;//!
    Double_t fTreeCascVarV0MassAsLambda;//!
    Double_t fTreeCascVarV0MassAsAntiLambda;//!
    Double_t fTreeCascVarV0Radius;//!
    Double_t fTreeCascVarV0CosPA;//!
    Double_t fTreeCascVarWrongCosPA;//!
    Double_t fTreeCascVarDCABachToBaryon;//!
    Int_t fTreeCascVarLeastNbrCrossedRows;//!
    Int_t fTreeCascVarLeastNbrClusters;//!
    Double_t fTreeCascVarLeastRatioCrossedRowsOverFindable;//!
    Double_t fTreeCascVarNbrCrossedRowsOverLength;//!
    Double_t fTreeCascVarMaxChi2PerCluster;//!
    Double_t fTreeCascVarMinTrackLength;  //!  
    //-----------MULTIPLICITY-INFO--------------------
    Double_t fTreeCascVarCentrality;          //! VZERO Percentile             - AliMultSelection
    Double_t fTreeCascVarVZEROMultSel;        //! VZERO Multiplicity           - AliMultSelection
    Int_t fTreeCascVarVZEROMultSig;        //! VZERO Multiplicity           - VZER0 signal
    Int_t fTreeCascVarVZEROMultSigCorr;    //! Corrected VZERO Multiplicity - VZER0 signal
    Int_t fTreeCascVarSPDMult;             //! SPD Multiplicity             - SPD Tracklets
    UInt_t fTreeCascVar_TriggerMask; //! save full info for checking later
    Bool_t fTreeCascVarMVPileupFlag; //!
    Bool_t fTreeCascVarOOBPileupFlag; //!
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t fTreeCascVarDistOverTotMom;//!
    Double_t fTreeCascVarV0DistOverTotMom;//!
    //------------------------------------------------
    //---------------PID-TPC-INFO---------------------
    Float_t fTreeCascVarBachNSigmaPion;//!
    Float_t fTreeCascVarBachNSigmaKaon;//!
    Float_t fTreeCascVarPosNSigmaProton;//!
    Float_t fTreeCascVarPosNSigmaPion;//!
    Float_t fTreeCascVarNegNSigmaProton;//!
    Float_t fTreeCascVarNegNSigmaPion;//!
    //---------------PID-TOF-INFO---------------------
    Float_t fTreeCascVarBachTOFNSigmaPion;//!
    Float_t fTreeCascVarBachTOFNSigmaKaon;//!
    Float_t fTreeCascVarPosTOFNSigmaProton;//!
    Float_t fTreeCascVarPosTOFNSigmaPion;//!
    Float_t fTreeCascVarNegTOFNSigmaProton;//!
    Float_t fTreeCascVarNegTOFNSigmaPion;//!
    //---------------PID-ITS-INFO---------------------
    Float_t fTreeCascVarBachITSNSigmaPion;//!
    Float_t fTreeCascVarBachITSNSigmaKaon;//!
    Float_t fTreeCascVarPosITSNSigmaProton;//!
    Float_t fTreeCascVarPosITSNSigmaPion;//!
    Float_t fTreeCascVarNegITSNSigmaProton;//!
    Float_t fTreeCascVarNegITSNSigmaPion;//!
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreeCascVarPosdEdx;//!
    Double_t fTreeCascVarNegdEdx;//!
    Double_t fTreeCascVarBachdEdx;//!
    Double_t fTreeCascVarPosPIDForTracking;//!
    Double_t fTreeCascVarNegPIDForTracking;//!
    Double_t fTreeCascVarBachPIDForTracking;//!
    //------------------------------------------------
    Double_t fTreeCascVarChi2Cascade;//!
    Double_t fTreeCascVarChi2CascadePerNDF;//!
    Double_t fTreeCascVarChi2V0;//!
    //------------------------------------------------
    ULong64_t fTreeCascVarBachTrackStatus;//!
    ULong64_t fTreeCascVarPosTrackStatus; //!
    ULong64_t fTreeCascVarNegTrackStatus; //!
    //------------------------------------------------
    Double_t fTreeCascVarBachDCAz; //!
    Double_t fTreeCascVarPosDCAz; //!
    Double_t fTreeCascVarNegDCAz; //!
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreeCascVarBachPx; //!
    Double_t fTreeCascVarBachPy; //!
    Double_t fTreeCascVarBachPz; //!
    Double_t fTreeCascVarPosPx; //!
    Double_t fTreeCascVarPosPy; //!
    Double_t fTreeCascVarPosPz; //!
    Double_t fTreeCascVarNegPx; //!
    Double_t fTreeCascVarNegPy; //!
    Double_t fTreeCascVarNegPz; //!
    //------------------------------------------------
    Double_t fTreeCascVarPrimVertexX; //!
    Double_t fTreeCascVarPrimVertexY; //!
    Double_t fTreeCascVarPrimVertexZ; //!
    Double_t fTreeCascVarCascadeDecayX; //!
    Double_t fTreeCascVarCascadeDecayY; //!
    Double_t fTreeCascVarCascadeDecayZ; //!
    Double_t fTreeCascVarV0DecayX; //!
    Double_t fTreeCascVarV0DecayY; //!
    Double_t fTreeCascVarV0DecayZ; //!
    //------------------------------------------------
    Int_t fTreeCascVarBachIndex; //!
    Int_t fTreeCascVarPosIndex; //!
    Int_t fTreeCascVarNegIndex; //!
    //------------------------------------------------
    Int_t fTreeCascVarBachLabel;//!
    Int_t fTreeCascVarPosLabel;//!
    Int_t fTreeCascVarNegLabel;//!
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreeCascVarBachITSClusters0;//!
    Bool_t fTreeCascVarBachITSClusters1;//!
    Bool_t fTreeCascVarBachITSClusters2;//!
    Bool_t fTreeCascVarBachITSClusters3;//!
    Bool_t fTreeCascVarBachITSClusters4;//!
    Bool_t fTreeCascVarBachITSClusters5;//!
    
    Bool_t fTreeCascVarPosITSClusters0;//!
    Bool_t fTreeCascVarPosITSClusters1;//!
    Bool_t fTreeCascVarPosITSClusters2;//!
    Bool_t fTreeCascVarPosITSClusters3;//!
    Bool_t fTreeCascVarPosITSClusters4;//!
    Bool_t fTreeCascVarPosITSClusters5;//!
    
    Bool_t fTreeCascVarNegITSClusters0;//!
    Bool_t fTreeCascVarNegITSClusters1;//!
    Bool_t fTreeCascVarNegITSClusters2;//!
    Bool_t fTreeCascVarNegITSClusters3;//!
    Bool_t fTreeCascVarNegITSClusters4;//!
    Bool_t fTreeCascVarNegITSClusters5;//!
    
    //------------------------------------------------
    Bool_t fTreeCascVarPosITSSharedClusters0;//!
    Bool_t fTreeCascVarPosITSSharedClusters1;//!
    Bool_t fTreeCascVarPosITSSharedClusters2;//!
    Bool_t fTreeCascVarPosITSSharedClusters3;//!
    Bool_t fTreeCascVarPosITSSharedClusters4;//!
    Bool_t fTreeCascVarPosITSSharedClusters5;//!
    
    Bool_t fTreeCascVarNegITSSharedClusters0;//!
    Bool_t fTreeCascVarNegITSSharedClusters1;//!
    Bool_t fTreeCascVarNegITSSharedClusters2;//!
    Bool_t fTreeCascVarNegITSSharedClusters3;//!
    Bool_t fTreeCascVarNegITSSharedClusters4;//!
    Bool_t fTreeCascVarNegITSSharedClusters5;//!
    
    Bool_t fTreeCascVarBachITSSharedClusters0;//!
    Bool_t fTreeCascVarBachITSSharedClusters1;//!
    Bool_t fTreeCascVarBachITSSharedClusters2;//!
    Bool_t fTreeCascVarBachITSSharedClusters3;//!
    Bool_t fTreeCascVarBachITSSharedClusters4;//!
    Bool_t fTreeCascVarBachITSSharedClusters5;//!

    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreeCascVarBachTOFExpTDiff; //!
    Double_t fTreeCascVarPosTOFExpTDiff; //!
    Double_t fTreeCascVarNegTOFExpTDiff; //!
    
    Double_t fTreeCascVarBachTOFSignal; //!
    Double_t fTreeCascVarPosTOFSignal; //!
    Double_t fTreeCascVarNegTOFSignal; //!
    
    Int_t   fTreeCascVarBachTOFBCid; //!
    Int_t   fTreeCascVarPosTOFBCid; //!
    Int_t   fTreeCascVarNegTOFBCid; //!
    
    //Kink tagging
    Bool_t fTreeCascVarBachIsKink;//!
    Bool_t fTreeCascVarPosIsKink;//!
    Bool_t fTreeCascVarNegIsKink;//!
    
    //Cowboy/sailor studies
    Bool_t  fTreeCascVarIsCowboy;   //! store if V0 is cowboy-like or sailor-like in XY plane
    Double_t fTreeCascVarCowboyness; //! negative -> cowboy, positive -> sailor
    Bool_t  fTreeCascVarIsCascadeCowboy;   //! store if V0 is cowboy-like or sailor-like in XY plane
    Double_t fTreeCascVarCascadeCowboyness; //! negative -> cowboy, positive -> sailor
    
    //Event info
    Double_t fTreeCascVarAmplitudeV0A; //!
    Double_t fTreeCascVarAmplitudeV0C; //!
    Int_t   fTreeCascVarClosestNonEmptyBC; //!
    
    Double_t fTreeCascVarMagField;//!
    Int_t fTreeCascVarRunNumber;//!
    ULong64_t fTreeCascVarEventNumber; //!
    
    //------------MONTE-CARLO-INFO------------------
    Double_t fTreeCascVarBachEtaMC;//!
    Double_t fTreeCascVarPosEtaMC;//!
    Double_t fTreeCascVarNegEtaMC;//!
    
    Double_t fTreeCascVarBachPxMC; //!
    Double_t fTreeCascVarBachPyMC; //!
    Double_t fTreeCascVarBachPzMC; //!
    Double_t fTreeCascVarPosPxMC; //!
    Double_t fTreeCascVarPosPyMC; //!
    Double_t fTreeCascVarPosPzMC; //!
    Double_t fTreeCascVarNegPxMC; //!
    Double_t fTreeCascVarNegPyMC; //!
    Double_t fTreeCascVarNegPzMC; //!
    
    Int_t fTreeCascVarBachMotherLabel;//!
    Int_t fTreeCascVarPosMotherLabel;//!
    Int_t fTreeCascVarNegMotherLabel;//!
    
    Int_t fTreeCascVarBachGrandMotherLabel;//!
    Int_t fTreeCascVarPosGrandMotherLabel;//!
    Int_t fTreeCascVarNegGrandMotherLabel;//!
    
    Int_t fTreeCascVarBachPID;//!
    Int_t fTreeCascVarPosPID;//!
    Int_t fTreeCascVarNegPID;//!
    
    Int_t fTreeCascVarBachMotherPID;//!
    Int_t fTreeCascVarPosMotherPID;//!
    Int_t fTreeCascVarNegMotherPID;//!
    
    Int_t fTreeCascVarBachGrandMotherPID;//!
    Int_t fTreeCascVarPosGrandMotherPID;//!
    Int_t fTreeCascVarNegGrandMotherPID;//!
    
    Int_t fTreeCascVarPID;//!
    Double_t fTreeCascVarRapMC;//!
    Double_t fTreeCascVarInvMassMC;//!
    Double_t fTreeCascVarTrueMassMC;//!
    Double_t fTreeCascVarPtotMC;//!
    Double_t fTreeCascVarPtMC;//!
    Double_t fTreeCascVarPhiMC;//!
    Double_t fTreeCascVarThetaMC;//!
    
    Int_t fTreeCascVarV0PID;//!
    Double_t fTreeCascVarV0MassMC;//!
    Double_t fTreeCascVarV0TrueMassMC;//!
    Double_t fTreeCascVarV0PtotMC;//!
    Double_t fTreeCascVarV0PtMC;//!
    
    Double_t fTreeCascVarCascadeDecayXMC;//!
    Double_t fTreeCascVarCascadeDecayYMC;//!
    Double_t fTreeCascVarCascadeDecayZMC;//!
    
    Double_t fTreeCascVarV0DecayXMC;//!
    Double_t fTreeCascVarV0DecayYMC;//!
    Double_t fTreeCascVarV0DecayZMC;//!
    
    Double_t fTreeCascVarPrimVertexXMC;//!
    Double_t fTreeCascVarPrimVertexYMC;//!
    Double_t fTreeCascVarPrimVertexZMC;//!

    Int_t fTreeCascVarIsPhysicalPrimary;//!
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBach;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryPos;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryNeg;//!
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBachMother;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryPosMother;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryNegMother;//!
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBachGrandMother;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryPosGrandMother;//!
    Bool_t fTreeCascVarIsPhysicalPrimaryNegGrandMother;//!
    
//===========================================================================================
//   Variables for V0 Tree
//===========================================================================================
    //-----------BASIC-INFO---------------------------
    Double_t fTreeV0VarMassAsK0s; //!
    Double_t fTreeV0VarMassAsLambda; //!
    Double_t fTreeV0VarMassAsAntiLambda; //!
    Double_t fTreeV0VarRapK0Short; //!
    Double_t fTreeV0VarRapLambda; //!
    Double_t fTreeV0VarPosEta; //!
    Double_t fTreeV0VarNegEta; //!
    Double_t fTreeV0VarPhi; //!
    Double_t fTreeV0VarTheta; //!
    Double_t fTreeV0VarPtot; //!
    Double_t fTreeV0VarPt; //!
    
    //-----------INFO-FOR-CUTS------------------------
    Double_t fTreeV0VarAlpha; //!
    Double_t fTreeV0VarPtArm;//!
    Double_t fTreeV0VarDCAV0Dau; //!
    Double_t fTreeV0VarDCAV0ToPV; //!
    Double_t fTreeV0VarDCAPosToPV; //!
    Double_t fTreeV0VarDCANegToPV; //!
    Double_t fTreeV0VarCosPA; //!
    Double_t fTreeV0VarRadius; //!
    
    Int_t   fTreeV0VarLeastNbrCrossedRows;//!
    Int_t   fTreeV0VarLeastNbrClusters;//!
    Double_t fTreeV0VarLeastRatioCrossedRowsOverFindable;//!
    Double_t fTreeV0VarMaxChi2PerCluster; //!
    Double_t fTreeV0VarMinTrackLength; //!
    //-----------MULTIPLICITY-INFO--------------------
    Double_t fTreeV0VarCentrality;          //! VZERO Percentile             - AliMultSelection
    Double_t fTreeV0VarVZEROMultSel;        //! VZERO Multiplicity           - AliMultSelection
    Int_t fTreeV0VarVZEROMultSig;        //! VZERO Multiplicity           - VZER0 signal
    Int_t fTreeV0VarVZEROMultSigCorr;    //! Corrected VZERO Multiplicity - VZER0 signal
    Int_t fTreeV0VarSPDMult;             //! SPD Multiplicity             - SPD Tracklets
    UInt_t fTreeV0Var_TriggerMask; //! save full info for checking later
    Bool_t fTreeV0VarMVPileupFlag; //!
    Bool_t fTreeV0VarOOBPileupFlag; //!
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t fTreeV0VarDistOverTotMom;//!
    //---------------PID-TPC-INFO---------------------
    Float_t fTreeV0VarPosNSigmaProton;//!
    Float_t fTreeV0VarPosNSigmaPion;//!
    Float_t fTreeV0VarNegNSigmaProton;//!
    Float_t fTreeV0VarNegNSigmaPion;//!
    //---------------PID-TOF-INFO---------------------
    Float_t fTreeV0VarPosTOFNSigmaProton;//!
    Float_t fTreeV0VarPosTOFNSigmaPion;//!
    Float_t fTreeV0VarNegTOFNSigmaProton;//!
    Float_t fTreeV0VarNegTOFNSigmaPion;//!
    //---------------PID-ITS-INFO---------------------
    Float_t fTreeV0VarPosITSNSigmaProton;//!
    Float_t fTreeV0VarPosITSNSigmaPion;//!
    Float_t fTreeV0VarNegITSNSigmaProton;//!
    Float_t fTreeV0VarNegITSNSigmaPion;//!
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreeV0VarPosdEdx;//!
    Double_t fTreeV0VarNegdEdx;//!
    Double_t fTreeV0VarPosPIDForTracking; //! uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Double_t fTreeV0VarNegPIDForTracking; //!
    //------------------------------------------------
    Double_t fTreeV0VarChi2V0;         //!
    //------------------------------------------------
    ULong64_t fTreeV0VarPosTrackStatus; //!
    ULong64_t fTreeV0VarNegTrackStatus; //!
    //------------------------------------------------
    Int_t fTreeV0VarPosLabel;//!
    Int_t fTreeV0VarNegLabel;//!
    //------------------------------------------------
    Double_t fTreeV0VarPosDCAz; //!
    Double_t fTreeV0VarNegDCAz; //!
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreeV0VarPosPx; //!
    Double_t fTreeV0VarPosPy; //!
    Double_t fTreeV0VarPosPz; //!
    Double_t fTreeV0VarNegPx; //!
    Double_t fTreeV0VarNegPy; //!
    Double_t fTreeV0VarNegPz; //!
    //------------------------------------------------
    Double_t fTreeV0VarPrimVertexX; //!
    Double_t fTreeV0VarPrimVertexY; //!
    Double_t fTreeV0VarPrimVertexZ; //!
    Double_t fTreeV0VarV0DecayX; //!
    Double_t fTreeV0VarV0DecayY; //!
    Double_t fTreeV0VarV0DecayZ; //!
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreeV0VarPosITSClusters0;//!
    Bool_t fTreeV0VarPosITSClusters1;//!
    Bool_t fTreeV0VarPosITSClusters2;//!
    Bool_t fTreeV0VarPosITSClusters3;//!
    Bool_t fTreeV0VarPosITSClusters4;//!
    Bool_t fTreeV0VarPosITSClusters5;//!
    
    Bool_t fTreeV0VarNegITSClusters0;//!
    Bool_t fTreeV0VarNegITSClusters1;//!
    Bool_t fTreeV0VarNegITSClusters2;//!
    Bool_t fTreeV0VarNegITSClusters3;//!
    Bool_t fTreeV0VarNegITSClusters4;//!
    Bool_t fTreeV0VarNegITSClusters5;//!
    //------------------------------------------------
    Bool_t fTreeV0VarPosITSSharedClusters0;//!
    Bool_t fTreeV0VarPosITSSharedClusters1;//!
    Bool_t fTreeV0VarPosITSSharedClusters2;//!
    Bool_t fTreeV0VarPosITSSharedClusters3;//!
    Bool_t fTreeV0VarPosITSSharedClusters4;//!
    Bool_t fTreeV0VarPosITSSharedClusters5;//!
    
    Bool_t fTreeV0VarNegITSSharedClusters0;//!
    Bool_t fTreeV0VarNegITSSharedClusters1;//!
    Bool_t fTreeV0VarNegITSSharedClusters2;//!
    Bool_t fTreeV0VarNegITSSharedClusters3;//!
    Bool_t fTreeV0VarNegITSSharedClusters4;//!
    Bool_t fTreeV0VarNegITSSharedClusters5;//!
    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreeV0VarNegTOFExpTDiff; //!
    Double_t fTreeV0VarPosTOFExpTDiff; //!
    Double_t fTreeV0VarNegTOFSignal; //!
    Double_t fTreeV0VarPosTOFSignal; //!
    Int_t   fTreeV0VarNegTOFBCid; //!
    Int_t   fTreeV0VarPosTOFBCid; //! 
    
    //Kink tag
    Bool_t fTreeV0VarPosIsKink;//!
    Bool_t fTreeV0VarNegIsKink;//!
    
    //Cowboy/sailor studies
    Bool_t fTreeV0VarIsCowboy; //! store if V0 is cowboy-like or sailor-like in XY plane
    
    //Event info
    Double_t fTreeV0VarAmplitudeV0A; //!
    Double_t fTreeV0VarAmplitudeV0C; //!
    Int_t   fTreeV0VarClosestNonEmptyBC; //!
    
    Double_t fTreeV0VarMagField;//!
    Int_t fTreeV0VarRunNumber;//!
    ULong64_t fTreeV0VarEventNumber; //!
    
    //------------MONTE-CARLO-INFO------------------
    Int_t fTreeV0VarPID;//! PDG Code
    Double_t fTreeV0VarRapMC;//! Perfect Y
    Double_t fTreeV0VarPtotMC;//!        
    Double_t fTreeV0VarPtMC;//!        
    Double_t fTreeV0VarInvMassMC;//!
    Double_t fTreeV0VarTrueMassMC;//!    
    Double_t fTreeV0VarPhiMC;//! 
    Double_t fTreeV0VarThetaMC;//!   
    Double_t fTreeV0VarPosEtaMC;//!
    Double_t fTreeV0VarNegEtaMC;//!
    
    Double_t fTreeV0VarPosPxMC;//!   
    Double_t fTreeV0VarPosPyMC;//!    
    Double_t fTreeV0VarPosPzMC;//!
    Double_t fTreeV0VarNegPxMC;//!   
    Double_t fTreeV0VarNegPyMC;//!    
    Double_t fTreeV0VarNegPzMC;//!
    
    Int_t fTreeV0VarPosMotherLabel;//!
    Int_t fTreeV0VarNegMotherLabel;//!
    Int_t fTreeV0VarPosGrandMotherLabel;//!
    Int_t fTreeV0VarNegGrandMotherLabel;//!
    
    Int_t fTreeV0VarPosPID;//!
    Int_t fTreeV0VarNegPID;//!
    Int_t fTreeV0VarPosMotherPID;//!
    Int_t fTreeV0VarNegMotherPID;//!
    Int_t fTreeV0VarPosGrandMotherPID;//!
    Int_t fTreeV0VarNegGrandMotherPID;//!
    
    Int_t fTreeV0VarMotherPID;//!
    Double_t fTreeV0VarMotherPt;//!
    Double_t fTreeV0VarMotherRapMC;//!
    
    Double_t fTreeV0VarDecayXMC;//!
    Double_t fTreeV0VarDecayYMC;//!
    Double_t fTreeV0VarDecayZMC;//!
    
    Double_t fTreeV0VarPrimVertexXMC;//!
    Double_t fTreeV0VarPrimVertexYMC;//!
    Double_t fTreeV0VarPrimVertexZMC;//!
    
    Int_t fTreeV0VarIsPhysicalPrimary;//!
    Int_t fTreeV0VarMotherIsPhysicalPrimary;//!
    
    Bool_t fTreeV0VarIsPhysicalPrimaryPos;//!
    Bool_t fTreeV0VarIsPhysicalPrimaryNeg;//!
    Bool_t fTreeV0VarIsPhysicalPrimaryPosMother;//!
    Bool_t fTreeV0VarIsPhysicalPrimaryNegMother;//!
    Bool_t fTreeV0VarIsPhysicalPrimaryPosGrandMother;//!
    Bool_t fTreeV0VarIsPhysicalPrimaryNegGrandMother;//!

//===========================================================================================
//   Variables for Resonance Tree
//===========================================================================================
    Int_t       fTreeRsnVarCutIDrsn; //!
    Double_t    fTreeRsnVarPx; //!
    Double_t    fTreeRsnVarPy; //!
    Double_t    fTreeRsnVarPz; //!
    Double_t    fTreeRsnVarInvMass; //!
    Bool_t      fTreeRsnVarPassesOOBPileupCut; //! ITSRefit AND/OR TOF hit (that is, at least ITS or TOF hit)
    ULong64_t   fTreeRsnVarEventNumber; //!
    //-----------MULTIPLICITY-INFO--------------------
    Double_t    fTreeRsnVarCentrality;          //! VZERO Percentile             - AliMultSelection
    Double_t    fTreeRsnVarVZEROMultSel;        //! VZERO Multiplicity           - AliMultSelection
    Int_t       fTreeRsnVarVZEROMultSig;        //! VZERO Multiplicity           - VZER0 signal
    Int_t       fTreeRsnVarVZEROMultSigCorr;    //! Corrected VZERO Multiplicity - VZER0 signal
    Int_t       fTreeRsnVarSPDMult;             //! SPD Multiplicity             - SPD Tracklets
    UInt_t      fTreeRsnVar_TriggerMask; //! save full info for checking later
    Bool_t      fTreeRsnVarMVPileupFlag; //!
    Bool_t      fTreeRsnVarOOBPileupFlag; //!
    //Event info
    Double_t    fTreeRsnVarAmplitudeV0A; //!
    Double_t    fTreeRsnVarAmplitudeV0C; //!
    Int_t       fTreeRsnVarClosestNonEmptyBC; //!
    Double_t    fTreeRsnVarMagField;//!
    Int_t       fTreeRsnVarRunNumber;//!
    //------------------------------------------------
    Double_t    fTreeRsnVarPrimVertexX; //!
    Double_t    fTreeRsnVarPrimVertexY; //!
    Double_t    fTreeRsnVarPrimVertexZ; //!
    //------------MONTE-CARLO-INFO------------------
    Bool_t      fTreeRsnVarDauFromSameMother;//!
    Int_t       fTreeRsnVarPID;//!
    Int_t       fTreeRsnVarDau1PID;//!
    Int_t       fTreeRsnVarDau2PID;//!
    Double_t    fTreeRsnVarPxMC;//!
    Double_t    fTreeRsnVarPyMC;//!
    Double_t    fTreeRsnVarPzMC;//!
    Double_t    fTreeRsnVarInvMassMC;//!
    Double_t    fTreeRsnVarPrimVertexXMC;//!
    Double_t    fTreeRsnVarPrimVertexYMC;//!
    Double_t    fTreeRsnVarPrimVertexZMC;//!
//===========================================================================================
//   Variables for Mixed Resonance Tree
//===========================================================================================
    Int_t       fTreeRsnFoundMixEvts;//!
    
    Int_t       fTreeRsnMixedVarCutIDrsn; //!
    Double_t    fTreeRsnMixedVarPx; //!
    Double_t    fTreeRsnMixedVarPy; //!
    Double_t    fTreeRsnMixedVarPz; //!
    Double_t    fTreeRsnMixedVarInvMass; //!
    Bool_t      fTreeRsnMixedVarPassesOOBPileupCut; //! ITSRefit AND/OR TOF hit (that is, at least ITS or TOF hit)
    ULong64_t   fTreeRsnMixedVarEventNumber; //!
    //------------MONTE-CARLO-INFO------------------
    Int_t       fTreeRsnMixedVarDau1PID;//!
    Int_t       fTreeRsnMixedVarDau2PID;//!
    Double_t    fTreeRsnMixedVarPxMC;//!
    Double_t    fTreeRsnMixedVarPyMC;//!
    Double_t    fTreeRsnMixedVarPzMC;//!
    Double_t    fTreeRsnMixedVarInvMassMC;//!
    
//===========================================================================================
//   Variables for Primary tracks Tree
//===========================================================================================
    //-----------BASIC-INFO---------------------------
    Int_t fTreePrimVarCharge; //!
    Double_t fTreePrimVarRapPion; //!
    Double_t fTreePrimVarRapProton; //!
    Double_t fTreePrimVarRapKaon; //!
    Double_t fTreePrimVarEta; //!
    Double_t fTreePrimVarTheta; //!
    Double_t fTreePrimVarPhi; //!
    Double_t fTreePrimVarPtot; //!
    Double_t fTreePrimVarPt; //!
    
    Double_t fTreePrimVarDCAxyToPV;//!
    Double_t fTreePrimVarDCAzToPV;//!
    
    Int_t fTreePrimVarNbrCrossedRows;//!
    Double_t fTreePrimVarRatioCrossedRowsOverFindable;//!
    Int_t fTreePrimVarNbrClusters;//!
    Double_t fTreePrimVarNbrCrossedRowsOverLength;//!
    Double_t fTreePrimVarFractionSharedTPCClusters;//!
    Double_t fTreePrimVarITSChi2PerCluster;//!
    Double_t fTreePrimVarTPCChi2PerCluster;//!
    Double_t fTreePrimVarTrackLength;  //!  
    //-----------MULTIPLICITY-INFO--------------------
    Double_t fTreePrimVarCentrality;          //! VZERO Percentile             - AliMultSelection
    Double_t fTreePrimVarVZEROMultSel;        //! VZERO Multiplicity           - AliMultSelection
    Int_t fTreePrimVarVZEROMultSig;        //! VZERO Multiplicity           - VZER0 signal
    Int_t fTreePrimVarVZEROMultSigCorr;    //! Corrected VZERO Multiplicity - VZER0 signal
    Int_t fTreePrimVarSPDMult;             //! SPD Multiplicity             - SPD Tracklets
    UInt_t fTreePrimVar_TriggerMask; //! save full info for checking later
    Bool_t fTreePrimVarMVPileupFlag; //!
    Bool_t fTreePrimVarOOBPileupFlag; //!
    //---------------PID-TPC-INFO---------------------
    Float_t fTreePrimVarNSigmaPion;//!
    Float_t fTreePrimVarNSigmaKaon;//!
    Float_t fTreePrimVarNSigmaProton;//!
    //---------------PID-TOF-INFO---------------------
    Float_t fTreePrimVarTOFNSigmaPion;//!
    Float_t fTreePrimVarTOFNSigmaKaon;//!
    Float_t fTreePrimVarTOFNSigmaProton;//!
    //---------------PID-ITS-INFO---------------------
    Float_t fTreePrimVarITSNSigmaPion;//!
    Float_t fTreePrimVarITSNSigmaKaon;//!
    Float_t fTreePrimVarITSNSigmaProton;//!
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreePrimVardEdx;//!
    Double_t fTreePrimVarPIDForTracking;//!
    //------------------------------------------------
    Double_t fTreePrimVarPrimVertexX; //!
    Double_t fTreePrimVarPrimVertexY; //!
    Double_t fTreePrimVarPrimVertexZ; //!
    //------------------------------------------------
    ULong64_t fTreePrimVarTrackStatus;//!
    //------------------------------------------------
    Int_t fTreePrimVarLabel;//!
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreePrimVarPx; //!
    Double_t fTreePrimVarPy; //!
    Double_t fTreePrimVarPz; //!
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreePrimVarITSClusters0;//!
    Bool_t fTreePrimVarITSClusters1;//!
    Bool_t fTreePrimVarITSClusters2;//!
    Bool_t fTreePrimVarITSClusters3;//!
    Bool_t fTreePrimVarITSClusters4;//!
    Bool_t fTreePrimVarITSClusters5;//!
    //------------------------------------------------
    Bool_t fTreePrimVarITSSharedClusters0;//!
    Bool_t fTreePrimVarITSSharedClusters1;//!
    Bool_t fTreePrimVarITSSharedClusters2;//!
    Bool_t fTreePrimVarITSSharedClusters3;//!
    Bool_t fTreePrimVarITSSharedClusters4;//!
    Bool_t fTreePrimVarITSSharedClusters5;//!

    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreePrimVarTOFExpTDiff; //!
    Double_t fTreePrimVarTOFSignal; //!
    Int_t   fTreePrimVarTOFBCid; //!
    //Event info
    Double_t fTreePrimVarAmplitudeV0A; //!
    Double_t fTreePrimVarAmplitudeV0C; //!
    Int_t   fTreePrimVarClosestNonEmptyBC; //!
    
    //Kink tagging
    Bool_t fTreePrimVarIsKink;//!
    
    Double_t fTreePrimVarMagField;//!
    Int_t fTreePrimVarRunNumber;//!
    ULong64_t fTreePrimVarEventNumber; //!
    
    //------------MONTE-CARLO-INFO------------------
    Int_t fTreePrimVarPID;//!
    Double_t fTreePrimVarRapMC;//!
    Double_t fTreePrimVarEtaMC;//!
    Double_t fTreePrimVarPtotMC;//!
    Double_t fTreePrimVarPtMC;//!
    Double_t fTreePrimVarPxMC;//!
    Double_t fTreePrimVarPyMC;//!
    Double_t fTreePrimVarPzMC;//!
    Double_t fTreePrimVarPhiMC;//!
    Double_t fTreePrimVarThetaMC;//!

    Double_t fTreePrimVarPrimVertexXMC;//!
    Double_t fTreePrimVarPrimVertexYMC;//!
    Double_t fTreePrimVarPrimVertexZMC;//!

    Bool_t fTreePrimVarIsPhysicalPrimary;//!
    
//===========================================================================================
//   Variables for MC tracks Tree
//===========================================================================================
    Int_t       fTreeMCVarPID;//!
    Int_t       fTreeMCVarCharge;//!
    Double_t    fTreeMCVarRap;//!
    Double_t    fTreeMCVarEta;//!
    Double_t    fTreeMCVarPtot;//!
    Double_t    fTreeMCVarPt;//!
    Double_t    fTreeMCVarPx;//!
    Double_t    fTreeMCVarPy;//!
    Double_t    fTreeMCVarPz;//!
    Double_t    fTreeMCVarTheta;//!
    Double_t    fTreeMCVarPhi;//!
    ULong64_t   fTreeMCVarEventNumber;//!
    
//===========================================================================================
//   Histograms
//===========================================================================================
    TList*                  fOutputList;    //! output list
    TH1D*                   fHistEventCounter;        //! dummy histogram

    AliAnalysisTaskStrangeCascadesTriggerAODRun2MC(const AliAnalysisTaskStrangeCascadesTriggerAODRun2MC&); // not implemented
    AliAnalysisTaskStrangeCascadesTriggerAODRun2MC& operator=(const AliAnalysisTaskStrangeCascadesTriggerAODRun2MC&); // not implemented

    ClassDef(AliAnalysisTaskStrangeCascadesTriggerAODRun2MC, 1);
};

//_____________________________________________________________________________
/// Add a new AliRsnMiniResonanceFinder object.
///
/// A user can add as many as he wants, and each one corresponds
/// to one of the available bits in the AliRsnMiniParticle mask.
/// The only check is the following: if a ResonanceFinder set with the same name
/// as the argument is there, this is not added.
///
/// \param f Pointer to the AliRsnMiniResonanceFinder object
/// \return Cut ID for the ResonanceFinder f.
///
inline Int_t AliAnalysisTaskStrangeCascadesTriggerAODRun2MC::AddResonanceFinder(AliRsnMiniResonanceFinder* f)
{
    TObject *obj = fResonanceFinders.FindObject(f->GetName());
    Int_t v = 0;

    if (obj) {
        AliInfo(Form("A ResonanceFinder named '%s' already exists", f->GetName()));
        v = fResonanceFinders.IndexOf(obj) + GetNumberOfTrackCuts();
    } else {
        fResonanceFinders.AddLast(f);
        v = fResonanceFinders.IndexOf(f) + GetNumberOfTrackCuts();
        f->SetResonanceCutID(v);
    }

    return v;
}

//_____________________________________________________________________________
/// Add a new cut set for a new criterion for track selection.
/// A user can add as many as (s)he wants, and each one corresponds
/// to one of the available bits in the AliRsnMiniParticle mask.
/// The only check is the following: if a cut set with the same name
/// as the argument is there, this is not added.
///
/// \param cuts Pointer to an AliRsnCutSet object
/// \return A value that is the array position of this set.
///
inline Int_t AliAnalysisTaskStrangeCascadesTriggerAODRun2MC::AddTrackCuts(AliRsnCutSet *cuts)
{
//
// Add a new cut set for a new criterion for track selection.
// A user can add as many as he wants, and each one corresponds
// to one of the available bits in the AliRsnMiniParticle mask.
// The only check is the following: if a cut set with the same name
// as the argument is there, this is not added.
// Return value is the array position of this set.
//

    TObject *obj = fRsnTrackCuts.FindObject(cuts->GetName());
    Int_t v = 0;

    if (obj) {
        AliInfo(Form("A cut set named '%s' already exists", cuts->GetName()));
        v = fRsnTrackCuts.IndexOf(obj);
    } else {
        fRsnTrackCuts.AddLast(cuts);
        v = fRsnTrackCuts.IndexOf(cuts);
    }

    for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
        AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
        if(f) f->IncrementResonanceCutID();
    }
        
    return v;
}

//_____________________________________________________________________________
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2MC::AddConfiguration( AliCascadeResult *lCascadeResult )
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
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2MC::AddConfiguration( AliV0Result *lV0Result )
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
inline void AliAnalysisTaskStrangeCascadesTriggerAODRun2MC::AddConfiguration( AliESDtrackCuts *ltrackCuts )
{
    //
    // Define a cut set for primary tracks selection.
    //
    fPrimTrackCuts = ltrackCuts;
}

#endif
