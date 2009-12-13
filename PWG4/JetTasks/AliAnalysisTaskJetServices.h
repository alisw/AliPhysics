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
    virtual void SetZVertexCut(Float_t f){fZVtxCut = f;}
    virtual Bool_t Notify();

    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}

    enum { kAllTriggered = 0,kTriggeredSPDVertex,kTriggeredVertexIn,kSelected,kConstraints};

 private:

    AliAnalysisTaskJetServices(const AliAnalysisTaskJetServices&);
    AliAnalysisTaskJetServices& operator=(const AliAnalysisTaskJetServices&);

    Bool_t        fUseAODInput;        // take jet from input AOD not from ouptu AOD
    Float_t       fAvgTrials;          // Average number of trials
    Float_t       fZVtxCut;          // Average number of trials
    TProfile*     fh1Xsec;             // pythia cross section and trials
    TH1F*         fh1Trials;           // trials are added
    TH1F*         fh1PtHard;           // Pt har of the event...       
    TH1F*         fh1PtHardTrials;     // Number of trials 
    TH2F*         fh2TriggerCount;     // number of fire triggers in each case
    TH2F*         fh2ESDTriggerCount;  // number of fire triggers in each case
    TH2F*         fh2TriggerVtx;       // vtx. position vs. trigger decision
    TH2F*         fh2ESDTriggerVtx;  // vtx. position vs. trigger decision 

    TList *fHistList; // Output list
   
    ClassDef(AliAnalysisTaskJetServices,1)
};
 
#endif
