#ifndef ALIANALYSISTASKJETSERVICES_H
#define ALIANALYSISTASKJETSERVICES_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************

#include  "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF

////////////////
class AliJetHeader;
class AliESDEvent;
class AliAODEvent;
class AliAODJet;
class AliGenPythiaEventHeader;
class AliCFManager;

class TList;
class TChain;
class TH1F;
class TH2F;
class TH3F;
class TProfile;



class AliAnalysisTaskJetServices : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetServices();
    AliAnalysisTaskJetServices(const char* name);
    virtual ~AliAnalysisTaskJetServices() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetZVertexCut(Float_t f){fVtxZCut = f;}
    virtual Bool_t Notify();

    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}
    virtual void SetRunRange(Float_t fLo,Float_t fUp){fRunRange[0] = fLo;fRunRange[1] = fUp;}
    virtual void SetRealData(Bool_t b){fRealData = b;}
    virtual void SetUsePhysicsSelection(Bool_t b){fUsePhysicsSelection = b;}
    Bool_t IsEventSelectedESD(AliESDEvent* esd);
    Bool_t IsEventPileUpESD(AliESDEvent* esd);
    Bool_t IsEventCosmicESD(AliESDEvent* esd);
    Bool_t IsEventSelectedAOD(AliAODEvent* aod);
    enum { kAllTriggered = 0,kTriggeredVertex,kTriggeredVertexIn,kSelectedALICE,kSelectedALICEVertexIn,kSelected,kConstraints};



 private:

    AliAnalysisTaskJetServices(const AliAnalysisTaskJetServices&);
    AliAnalysisTaskJetServices& operator=(const AliAnalysisTaskJetServices&);

    Bool_t        fUseAODInput;        // take jet from input AOD not from ouptu AOD
    Bool_t        fUsePhysicsSelection;// decide wether we take into account physicsselction task
    Bool_t        fRealData;           // true for real data to allow correct trigger slection
    UInt_t        fSelectionInfoESD;      // slection info bit mask
    Float_t       fAvgTrials;          // Average number of trials
    Float_t       fVtxXMean;          // mean x for cuts
    Float_t       fVtxYMean;          // mean y for cuts
    Float_t       fVtxZMean;          // mean z for cuts
    Float_t       fVtxRCut;           // vtx cut in R
    Float_t       fVtxZCut;           // vtzx cut in Z
    Float_t       fRunRange[2];        // only important for real data for 
    TProfile*     fh1Xsec;             // pythia cross section and trials
    TH1F*         fh1Trials;           // trials are added
    TH1F*         fh1PtHard;           // Pt har of the event...       
    TH1F*         fh1PtHardTrials;     // Number of trials 
    TH1F*         fh1SelectionInfoESD; // Masks that satisfy fSelectionInfo
    TH2F*         fh2TriggerCount;     // number of fire triggers in each case
    TH2F*         fh2ESDTriggerCount;  // number of fire triggers in each case
    TH2F*         fh2TriggerVtx;       // vtx. position vs. trigger decision
    TH2F*         fh2ESDTriggerVtx;  // vtx. position vs. trigger decision 
    TH2F*         fh2ESDTriggerRun;  // fired triggers vs. run number
    TH2F*         fh2VtxXY;          // XY position of VTX were available
    TList *fHistList; // Output list
   
    ClassDef(AliAnalysisTaskJetServices,5)
};
 
#endif
