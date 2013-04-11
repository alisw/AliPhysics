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
//
// Task for Heavy Flavour Electron Analysis
// Fills a single-inclusive electron pt-spectrum
// For further information see implementation file
//
#ifndef ALIANALYSISTASKHFE_H
#define ALIANALYSISTASKHFE_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

#ifndef ROOT_TBits
#include <TBits.h>
#endif

class AliAnalysisUtils;
class AliESDtrackCuts;
class AliHFEcontainer;
class AliHFEcollection;
class AliHFEcuts;
class AliHFEextraCuts;
class AliHFEelecbackground;
class AliHFENonPhotonicElectron;
class AliHFEmcQA;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliHFEsecVtx;
class AliHFEsignalCuts;
class AliHFEvarManager;
class AliHFEtaggedTrackAnalysis;
class AliCFManager;
class AliMCEvent;
class AliOADBContainer;
class AliAODMCHeader;
class AliVEvent;
class AliVParticle;
class AliTriggerAnalysis;
class TH1I; 
class TList;
class TClonesArray;

class AliAnalysisTaskHFE : public AliAnalysisTaskSE{
  public:
    enum{
      kPIDqa = 0,
      kMCqa =1 
    };
    enum{
      kPriVtx = 0,
      kSecVtx = 1,
      kIsElecBackGround = 2,
      kPostProcess = 3,
      kDEstep = 4,
      kTaggedTrackAnalysis = 5,
      kNonPhotonicElectron = 6
    };
    enum CreationProcess_t{
      kSignalCharm = 0,
      kSignalBeauty = 1,
      kGammaConv = 2,
      kOther = 3
    };
    enum{
      kBgPtBins = 44,
      kElecBgSpecies = 6,
      kCentBins = 11,
      kBgLevels = 3
    };
    typedef enum{
      kpp = 0,
      kpPb = 1,
      kPbPb = 2
    } ECollisionSystem_t;

    AliAnalysisTaskHFE();
    AliAnalysisTaskHFE(const char * name);
    AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref);
    AliAnalysisTaskHFE& operator=(const AliAnalysisTaskHFE &ref);
    virtual void Copy(TObject &o) const;
    virtual ~AliAnalysisTaskHFE();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    virtual Bool_t IsEventInBinZero();

    Bool_t IsQAOn(Int_t qaLevel) const { return TESTBIT(fQAlevel, qaLevel); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };
    Bool_t HasMCData() const { return TestBit(kHasMCdata); }
    Bool_t Ispp() const { return fCollisionSystem.TestBitNumber(kpp); }
    Bool_t IsPbPb() const { return fCollisionSystem.TestBitNumber(kPbPb); }
    Bool_t IspPb() const { return fCollisionSystem.TestBitNumber(kpPb); }
    Bool_t IsHeavyIon() const { return IsPbPb() || IspPb(); }
    Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };

    // Get Components for configuration
    AliHFEvarManager *GetVarManager() const { return fVarManager; }
    AliHFEpidQAmanager *GetPIDQAManager() const { return fPIDqa; }
    AliHFEpid *GetPID() const { return fPID; }
    AliHFENonPhotonicElectron *GetHFEBackgroundSubtraction() const { return fBackgroundSubtraction; }

    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetTaggedTrackCuts(AliHFEcuts * const cuts) { fTaggedTrackCuts = cuts; }
    void SetCleanTaggedTrack(Bool_t clean) { fCleanTaggedTrack = clean; };
    void SetVariablesTRDTaggedTrack(Bool_t variablesTRD) { fVariablesTRDTaggedTrack = variablesTRD; };
    void SetHFECutsPreselect(AliESDtrackCuts * const cuts) { fCutspreselect = cuts; };
    void SetHFEElecBackGround(AliHFEelecbackground * const elecBackGround) { fElecBackGround = elecBackGround; };
    void SetHFEBackgroundSubtraction(AliHFENonPhotonicElectron * const backgroundSubtraction) { fBackgroundSubtraction = backgroundSubtraction; };
    void SetQAOn(Int_t qaLevel) { SETBIT(fQAlevel, qaLevel); };
    void SwitchOnPlugin(Int_t plug);
    void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
    void SetFillSignalOnly(Bool_t signalOnly) { fFillSignalOnly = signalOnly; }
   
    void SetFillNoCuts(Bool_t fillNoCuts) { fFillNoCuts = fillNoCuts; }
    void SetApplyCutAOD(Bool_t applyCutAOD)   { fApplyCutAOD = applyCutAOD; }
    void SetRemovePileUp(Bool_t removePileUp) { fRemovePileUp = removePileUp; }
    void SetRemoveFirstEventInChunk() {fRemoveFirstEvent = kTRUE;}
    void SetPIDPreselect(AliHFEpid * const cuts) { fPIDpreselect = cuts; };
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetTRDTrigger(Bool_t activateTRDTrigger, Int_t trdtrigger) {fTRDTrigger=activateTRDTrigger; fWhichTRDTrigger=trdtrigger;};
    void SetCollisionSystem(ECollisionSystem_t system){
      fCollisionSystem.Clear();
      fCollisionSystem.SetBitNumber(system, kTRUE);
    }
    void SetppAnalysis(){
      fCollisionSystem.SetBitNumber(kpPb, kFALSE); 
      fCollisionSystem.SetBitNumber(kPbPb, kFALSE); 
      fCollisionSystem.SetBitNumber(kpp, kTRUE); 
    }
    void SetpPbAnalysis() {
      fCollisionSystem.SetBitNumber(kpp, kFALSE); 
      fCollisionSystem.SetBitNumber(kPbPb, kFALSE); 
      fCollisionSystem.SetBitNumber(kpPb, kTRUE); 
    }
    void SetPbPbAnalysis() { 
      fCollisionSystem.SetBitNumber(kpp, kFALSE); 
      fCollisionSystem.SetBitNumber(kpPb, kFALSE); 
      fCollisionSystem.SetBitNumber(kPbPb, kTRUE); 
    };
    void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator; }
    void SetPbPbUserCentralityLimit(Bool_t isPbPbUserBinning = kFALSE){fPbPbUserCentralityBinning = isPbPbUserBinning; };
    void SetPbPbUserCentralityArray(Int_t icentr, Float_t valuecentr) {fCentralityLimits[icentr] = valuecentr;};
    void SetPPMultiBinAnalysis(Bool_t isppMultiBin) { fisppMultiBin = isppMultiBin; };
    void SetNonHFEsystematics(Bool_t isSystematics) {fisNonHFEsystematics = isSystematics; };
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
    void SetBackGroundFactorsFunction(const TF1 * const backGroundFactorsFunction, Int_t centralitybin=0){  
      fkBackGroundFactorArray[centralitybin]=backGroundFactorsFunction;
      fBackGroundFactorApply=kTRUE;
      SetBit(kBackgroundInitialized);
    };
    void SetElecBackGroundFactors(Int_t iPt, Int_t iType, Int_t iCent, Int_t iError, Double_t elecBackGroundFactor) {fElecBackgroundFactor[iError][iCent][iType][iPt] = elecBackGroundFactor; };
    void SetBinLimits(Int_t iPt, Double_t momentum){fBinLimit[iPt] = momentum;};
    void PrintStatus() const;
    Bool_t ReadCentrality();
    void RejectionPileUpVertexRangeEventCut();  
    void SelectSpecialTrigger(const Char_t *trgclust, Int_t runMin = 0, Int_t runMax = 999999999); 
    void SetDebugStreaming() {SetBit(kTreeStream);};
    Bool_t CheckTRDTrigger(AliESDEvent *ev);
    void DrawTRDTrigger(AliESDEvent *ev);

  private:
    enum{
      kHasMCdata = BIT(19),
      kAODanalysis = BIT(20),
      kBackgroundInitialized = BIT(21),
      kTreeStream = BIT(22)
    };

    Bool_t FillProductionVertex(const AliVParticle * const track) const;
    void MakeParticleContainer();
    void MakeEventContainer();
    void InitHistoITScluster();
    void InitContaminationQA();
    const Char_t *GetSpecialTrigger(Int_t run);
    void ProcessMC();
    void ProcessESD();
    void ProcessAOD();
    Int_t GetITSMultiplicity(AliVEvent *ev);
    Bool_t PreSelectTrack(AliESDtrack *track) const;
    Bool_t ProcessMCtrack(AliVParticle *track);
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    AliAODMCHeader *fAODMCHeader;         // ! MC info AOD
    TClonesArray *fAODArrayMCInfo;        // ! MC info particle AOD
    ULong_t fQAlevel;                     // QA level
    UShort_t fPlugins;                    // Enabled Plugins
    TBits fCollisionSystem;              // Collision System;
    Bool_t fFillSignalOnly;               // Fill container only with MC Signal Tracks
    Bool_t fFillNoCuts;                   // Fill container before any cut
    Bool_t fApplyCutAOD;                  // Apply the analysis cut for AOD tracks
    Bool_t fBackGroundFactorApply;        // Apply Background Function Subtraction,   MF: To be removed when transition to OADB container is finished
    Bool_t fRemovePileUp;                 // Remove Pile Up
    Bool_t fIdentifiedAsPileUp;           // Identified as pile-up
    Bool_t fIdentifiedAsOutInz;           // Out Of Range in z
    Bool_t fPassTheEventCut;              // Pass The Event Cut
    Bool_t fRejectKinkMother;             // Reject Kink Mother
    Bool_t fisppMultiBin;                 // pp Multiplicity Bin analysis
    Bool_t fPbPbUserCentralityBinning;    // PbPb user centrality binning
    Bool_t fRemoveFirstEvent;             // Remove first event from chunk
    Bool_t fisNonHFEsystematics;          // Non-HFE background systematics analysis
    AliOADBContainer *fSpecialTrigger;    // Special trigger selection
    Int_t   fCentralityF;                 // Centrality bin
    Float_t fCentralityPercent;           // Centrality percentile
    TString fCentralityEstimator;         // Centrality Estimator
    Float_t fContributors;                // Contributors
    Double_t fWeightBackGround;            // weight background function
    Double_t fVz;                         // z position of the primary vertex
    const TF1  *fkBackGroundFactorArray[12];   // Array of BackGround factors for each centrality bin, bin0 = min bias
    Double_t fElecBackgroundFactor[kBgLevels][kCentBins][kElecBgSpecies][kBgPtBins];     // Electron background factors
    Double_t fBinLimit[kBgPtBins+1];      // Electron pt bin edges
    Float_t fCentralityLimits[12];        // Limits for centrality bins
    AliHFEcontainer *fContainer;          //! The HFE container
    AliHFEvarManager *fVarManager;        // The var manager as the backbone of the analysis
    AliHFEsignalCuts *fSignalCuts;        //! MC true signal (electron coming from certain source) 
    AliCFManager *fCFM;                   //! Correction Framework Manager
    AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation
    AliHFEpid *fPID;                      // PID
    AliHFEpidQAmanager *fPIDqa;           // PID QA
    AliHFEpid *fPIDpreselect;             // PID oject for pre-selected tracks (without QA)
    AliHFEcuts *fCuts;                    // Cut Collection
    AliHFEcuts *fTaggedTrackCuts;         // Cut Collection for V0 tagged tracks
    Bool_t fCleanTaggedTrack;             // Loose cleaning of the V0 tagged tracks electron
    Bool_t fVariablesTRDTaggedTrack;      // Take the variables at the TRD for the V0 tagged tracks electron
    AliAnalysisUtils *fAnalysisUtils;     // Utility object to remove the first event of a chunk from the analysis
    AliESDtrackCuts *fCutspreselect;      // Cut Collection for pre-selected tracks
    AliHFEsecVtx *fSecVtx;                //! Secondary Vertex Analysis
    AliHFEelecbackground *fElecBackGround;//! Background analysis
    AliHFEmcQA *fMCQA;                    //! MC QA
    AliHFEtaggedTrackAnalysis *fTaggedTrackAnalysis;     //!Analyse V0-tagged tracks
    AliHFEextraCuts *fExtraCuts;          //! temporary implementation for IP QA
    AliHFENonPhotonicElectron *fBackgroundSubtraction; // Background subtraction
    Bool_t fTRDTrigger;                   // Check if event is TRD triggered event
    Int_t  fWhichTRDTrigger;               // Select type of TRD trigger
    //-----------QA and output---------------
    TList *fQA;                           //! QA histos for the cuts
    TList *fOutput;                       //! Container for Task Output
    TList *fHistMCQA;                     //! Output container for MC QA histograms 
    TList *fHistSECVTX;                   //! Output container for sec. vertexing results
    TList *fHistELECBACKGROUND;           //! Output container for electron background analysis
    AliHFEcollection *fQACollection;      //! Tasks own QA collection
    //---------------------------------------

    ClassDef(AliAnalysisTaskHFE, 2)       // The electron Analysis Task
};
#endif

