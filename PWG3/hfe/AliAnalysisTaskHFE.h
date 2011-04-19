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

class AliHFEcontainer;
class AliHFEcollection;
class AliHFEcuts;
class AliHFEelecbackground;
class AliHFEmcQA;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliHFEsecVtx;
class AliHFEsignalCuts;
class AliHFEvarManager;
class AliHFEtaggedTrackAnalysis;
class AliCFManager;
class AliVEvent;
class AliMCEvent;
class AliVParticle;
class AliTriggerAnalysis;
class TH1I; 
class TList;

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
      kTaggedTrackAnalysis = 5
    };
    enum CreationProcess_t{
      kSignalCharm = 0,
      kSignalBeauty = 1,
      kGammaConv = 2,
      kOther = 3
    };
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
    Bool_t IsPbPb() const { return TestBit(kBeamType); }
    Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };

    // Get Components for configuration
    AliHFEvarManager *GetVarManager() const { return fVarManager; }
    AliHFEpidQAmanager *GetPIDQAManager() const { return fPIDqa; }
    AliHFEpid *GetPID() const { return fPID; }

    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetTaggedTrackCuts(AliHFEcuts * const cuts) { fTaggedTrackCuts = cuts; }
    void SetCleanTaggedTrack(Bool_t clean) { fCleanTaggedTrack = clean; };
    void SetVariablesTRDTaggedTrack(Bool_t variablesTRD) { fVariablesTRDTaggedTrack = variablesTRD; };
    void SetHFECutsPreselect(AliHFEcuts * const cuts) { fCutspreselect = cuts; };
    void SetHFEElecBackGround(AliHFEelecbackground * const elecBackGround) { fElecBackGround = elecBackGround; };
    void SetQAOn(Int_t qaLevel) { SETBIT(fQAlevel, qaLevel); };
    void SwitchOnPlugin(Int_t plug);
    void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
    void SetFillSignalOnly(Bool_t signalOnly) { fFillSignalOnly = signalOnly; }
    void SetRemovePileUp(Bool_t removePileUp) { fRemovePileUp = removePileUp; }
    void SetPIDPreselect(AliHFEpid * const cuts) { fPIDpreselect = cuts; };
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetPbPbAnalysis(Bool_t isPbPb = kFALSE) { SetBit(kBeamType, isPbPb); };
    void SetBackGroundFactorsFunction(TF1 * const backGroundFactorsFunction, Int_t centralitybin=0)
    {  fBackGroundFactorArray[centralitybin]=backGroundFactorsFunction;
       fBackGroundFactorApply=kTRUE;};
    void PrintStatus() const;
    Bool_t ReadCentrality();
    void RejectionPileUpVertexRangeEventCut();  
    void SelectSpecialTrigger(const Char_t *trgclust){ fHasSpecialTriggerSelection = kTRUE; fSpecialTrigger = trgclust; }
 
  private:
    enum{
      kHasMCdata = BIT(19),
      kAODanalysis = BIT(20),
      kBeamType = BIT(21)
    };

    Bool_t FillProductionVertex(const AliVParticle * const track) const;
    void MakeParticleContainer();
    void MakeEventContainer();
    void InitPIDperformanceQA();
    void InitContaminationQA();
    void ProcessMC();
    void ProcessESD();
    void ProcessAOD();
    Bool_t PreSelectTrack(AliESDtrack *track) const;
    Bool_t ProcessMCtrack(AliVParticle *track);
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    ULong_t fQAlevel;                     // QA level
    UShort_t fPlugins;                    // Enabled Plugins
    Bool_t fFillSignalOnly;               // Fill container only with MC Signal Tracks
    Bool_t fBackGroundFactorApply;        // Apply Background Function Subtraction
    Bool_t fRemovePileUp;                 // Remove Pile Up
    Bool_t fIdentifiedAsPileUp;           // Identified as pile-up
    Bool_t fIdentifiedAsOutInz;           // Out Of Range in z
    Bool_t fPassTheEventCut;              // Pass The Event Cut
    Bool_t fHasSpecialTriggerSelection;   // Select special triggered events
    TString fSpecialTrigger;              // Special trigger selection
    Float_t fCentralityF;                 // Centrality
    Float_t fContributors;                // Contributors
    Double_t fWeightBackGround;            // weight background function
    Double_t fVz;                         // z position of the primary vertex
    TF1  *fBackGroundFactorArray[12];     // Array of BackGround factors for each centrality bin, bin0 = min bias
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
    AliHFEcuts *fCutspreselect;           // Cut Collection for pre-selected tracks
    AliHFEsecVtx *fSecVtx;                //! Secondary Vertex Analysis
    AliHFEelecbackground *fElecBackGround;//! Background analysis
    AliHFEmcQA *fMCQA;                    //! MC QA
    AliHFEtaggedTrackAnalysis *fTaggedTrackAnalysis;     //!Analyse V0-tagged tracks

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

