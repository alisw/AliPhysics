#ifndef ALIANALYSISTASKJETSPECTRUM2_H
#define ALIANALYSISTASKJETSPECTRUM2_H
 
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



class AliAnalysisTaskJetSpectrum2 : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetSpectrum2();
    AliAnalysisTaskJetSpectrum2(const char* name);
    virtual ~AliAnalysisTaskJetSpectrum2() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();


    virtual void SetExternalWeight(Float_t f){fExternalWeight = f;}
    virtual void SetUseExternalWeightOnly(Bool_t b){fUseExternalWeightOnly = b;}
    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}
    virtual void SetLimitGenJetEta(Bool_t b){fLimitGenJetEta = b;}
    virtual void SetRecEtaWindow(Float_t f){fRecEtaWindow = f;}
    virtual void SetAnalysisType(Int_t i){fAnalysisType = i;}
    virtual void SetBranchGen(const char* c){fBranchGen = c;}
    virtual void SetBranchRec(const char* c){fBranchRec = c;}
    virtual void SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
    virtual void SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
    virtual void SetFilterMask(UInt_t i){fFilterMask = i;}
    // use for the CF


    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackUndef = 0, kTrackAODIn, kTrackAODOut, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged };
    enum {kAnaMC =  0x1, kAnaMCESD = 0x2};
    enum {kMaxJets = 4};
    enum {kMaxCorrelation =  3};
    
    // 
    // 0 all jets
    // 1 all jet in eta window
    // 2 all jets with partner
    // 3 all jets in eta window with partner
    // 4 all jets with partner in eta window
    enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};


 private:

    AliAnalysisTaskJetSpectrum2(const AliAnalysisTaskJetSpectrum2&);
    AliAnalysisTaskJetSpectrum2& operator=(const AliAnalysisTaskJetSpectrum2&);

    void MakeJetContainer();
    Int_t GetListOfTracks(TList *list,Int_t type);

    AliJetHeader *fJetHeaderRec;
    AliJetHeader *fJetHeaderGen;
    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD
    THnSparseF   *fhnJetContainer[kMaxStep*2];   // like particle container in corrfw with different steps need AliCFContainer with Scale(), and clone() to do the same
    THnSparseF   *fhnCorrelation;           // response matrix for unfolding 


    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fBranchGen;  // AOD brnach for genereated

    Bool_t        fUseAODInput;           // take jet from input AOD not from ouptu AOD
    Bool_t        fUseExternalWeightOnly; // use only external weight
    Bool_t        fLimitGenJetEta;        // Limit the eta of the generated jets
    UInt_t        fFilterMask;             // filter bit for slecected tracks
    Int_t         fAnalysisType;          // Analysis type 
    Int_t         fTrackTypeRec;          // type of tracks used for FF 
    Int_t         fTrackTypeGen;          // type of tracks used for FF 
    Float_t       fAvgTrials;             // Average nimber of trials
    Float_t       fExternalWeight;        // external weight
    Float_t       fRecEtaWindow;          // eta window used for corraltion plots between rec and gen 

    TProfile*     fh1Xsec;   // pythia cross section and trials
    TH1F*         fh1Trials; // trials are added
    TH1F*         fh1PtHard;  // Pt har of the event...       
    TH1F*         fh1PtHardNoW;  // Pt har of the event without weigt       
    TH1F*         fh1PtHardTrials;  // Number of trials 
    TH1F*         fh1NGenJets;      // nr of gen jets
    TH1F*         fh1NRecJets;      // nr of rec jets
    
    TH1F*         fh1PtRecIn[kMaxJets];  // Jet pt for all this info is also in the THNsparse      
    TH1F*         fh1PtGenIn[kMaxJets];  // Jet pt with corellated generated jet    

    TH2F*         fh2FragRec[kMaxJets];     // fragmentation function
    TH2F*         fh2FragLnRec[kMaxJets];   // fragmetation in xi

    TH2F*         fh2FragGen[kMaxJets];     // fragmentation function
    TH2F*         fh2FragLnGen[kMaxJets];   // fragmetation in xi


    TList *fHistList; // Output list
   

    ClassDef(AliAnalysisTaskJetSpectrum2, 1) // Analysis task for standard jet analysis
};
 
#endif
