#ifndef ALIANALYSISTASKJETHBOM_H
#define ALIANALYSISTASKJETHBOM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// Task used for the correction of detector effects for background fluctuations in jet spectra by the HBOM method
// *******************************************

#include  "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "TH1F.h"
#include "TH2D.h"


////////////////
class AliJetHeader;
class AliESDEvent;
class AliAODEvent;
class AliAODExtension;
class AliAODJet;
class AliGenPythiaEventHeader;
class AliCFManager;
class AliAODJetEventBackground;
class AliJetFinder;
class TList;
class TChain;
class TH3F;
class TProfile;
class TRandom3;
class TRefArray;
class TClonesArray;
class TF1;

class AliAnalysisTaskJetHBOM : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetHBOM();
    AliAnalysisTaskJetHBOM(const char* name);
    virtual ~AliAnalysisTaskJetHBOM();
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();

    

    virtual void SetAODTrackInput(Bool_t b){fUseAODTrackInput = b;}
    virtual void SetAODMCInput(Bool_t b){fUseAODMCInput = b;}
    virtual void SetEventSelection(Bool_t b){fEventSelection = b;}
    virtual void SetTrackEtaWindow(Float_t f){fTrackEtaWindow = f;}
    virtual void SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
    virtual void SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
    virtual void SetTrackPtCut(Float_t x){fTrackPtCut = x;}
    virtual void SetCentralityCut(Float_t xLo,Float_t xUp){fCentCutLo = xLo; fCentCutUp = xUp;}
    virtual void SetFilterMask(UInt_t i,Int_t iType = 0){fFilterMask = i;
      fFilterType = iType;}
    virtual void SetJetTypes(UInt_t i){fJetTypes = i;}
    virtual void SetVtxCuts(Float_t z,Float_t r = 1){fVtxZCut = z; fVtxR2Cut = r *r;}    
    virtual void SetBackgroundBranch(const char* c){fBackgroundBranch = c;}
    virtual const char* GetBackgroundBranch(){return fBackgroundBranch.Data();}    
    virtual void SetNSkipLeadingRan(Int_t x){fNSkipLeadingRan = x;}
    virtual void SetNSkipLeadingCone(Int_t x){fNSkipLeadingCone = x;}
    virtual void SetNRandomCones(Int_t x){fNRandomCones = x;}
    virtual void SetRandConePos(Double_t eta, Double_t phi){randCone_pos=1;randCone_Eta=eta;randCone_Phi=phi;}

    virtual void SetfNHBOM(Int_t x){fNHBOM = x;};
    virtual void SetEfficiencyPt(TH1F *h){fh1efficiencyPt = (TH1F*)h->Clone("h1efficiencyPt");}
    virtual void SetEfficiencyPhi(TH2D *h){fh2efficiencyPhi = (TH2D*)h->Clone("h2efficiencyPhi");}

    virtual void SetJetOutputBranch(const char *c){fNonStdBranch = c;}
    virtual const char* GetJetOutputBranch(){return fNonStdBranch.Data();}
    virtual void SetJetOutputFile(const char *c){fNonStdFile = c;}
    virtual const char* GetJetOutputFile(){return fNonStdFile.Data();}
    virtual void SetMaxTrackPtInJet(Float_t x){fMaxTrackPtInJet = x;}
    virtual void SetJetOutputMinPt(Float_t x){fJetOutputMinPt = x;}

    //Setters for detector level effects
    virtual void SetSmearResolution(Bool_t b){fUseTrMomentumSmearing = b;} 
    virtual void SetDiceEfficiency(Bool_t b){fUseDiceEfficiency = b;} 
    virtual void SetMomentumResolutionHybrid(TProfile *p1, TProfile *p2, TProfile *p3);
    virtual void SetEfficiencyHybrid(TH1 *h1, TH1 *h2, TH1 *h3);

    Double_t GetMomentumSmearing(Int_t cat, Double_t pt);
    void FitMomentumResolution();


    // for Fast Jet
    fastjet::JetAlgorithm        GetAlgorithm()         const {return fAlgorithm;}
    fastjet::Strategy            GetStrategy()          const {return fStrategy;}
    fastjet::RecombinationScheme GetRecombScheme()      const {return fRecombScheme;}
    fastjet::AreaType            GetAreaType()          const {return fAreaType;}
    // Setters
    void SetRparam(Double_t f)                           {fRparam = f;}
    void SetAlgorithm(fastjet::JetAlgorithm f)           {fAlgorithm = f;}
    void SetStrategy(fastjet::Strategy f)                {fStrategy = f;}
    void SetRecombScheme(fastjet::RecombinationScheme f) {fRecombScheme = f;}
    void SetAreaType(fastjet::AreaType f)                {fAreaType = f;}
    void SetGhostArea(Double_t f) {fGhostArea = f;}
    void SetActiveAreaRepeats(Int_t f) {fActiveAreaRepeats = f;}
    void SetGhostEtamax(Double_t f) {fGhostEtamax = f;}



    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackUndef = 0, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance, kTrackAODextra, kTrackAODextraonly};
    enum {kMaxJets = 4};
    enum {kMaxCorrelation =  3};
    enum {kMaxRadius =       5};
    enum {kMaxCent =         4};
    enum {kJet = 1<<0,
	  kJetRan = 1<<1,	  
	  kRC = 1<<2,
	  kRCRan = 1<<3
    };
    

 private:

    AliAnalysisTaskJetHBOM(const AliAnalysisTaskJetHBOM&);
    AliAnalysisTaskJetHBOM& operator=(const AliAnalysisTaskJetHBOM&);

    Int_t GetListOfTracks(TList *list,Int_t type);

    AliAODEvent     *fAOD;                // ! where we take the jets from can be input or output AOD
    AliAODExtension *fAODExtension;       // ! AOD extension in case we write a non-sdt branch to a separate file and the aod is standard
    TRefArray       *fRef;               // ! trefarray for track references within the jet
    Bool_t        fUseAODTrackInput;      // take track from input AOD not from ouptu AOD
    Bool_t        fUseAODMCInput;         // take MC from input AOD not from ouptu AOD
    Bool_t        fEventSelection;        // use the event selection of this task, otherwise analyse all
    UInt_t        fFilterMask;            // filter bit for slecected tracks
    UInt_t        fFilterMaskBestPt;      // filter bit to mark jets with high quality leading tracks

    UInt_t        fFilterType;            // filter type 0 = all, 1 = ITSTPC, 2 = TPC
    UInt_t        fJetTypes;              // 1<<0 regular jets, 1<<1 << randomized event 1<<2 random cones 1<<3 random cones randomiuzed event
    Int_t         fTrackTypeRec;          // type of tracks used for FF 
    Int_t         fTrackTypeGen;          // type of tracks used for FF 
    Int_t         fNSkipLeadingRan;       // number of leading tracks to be skipped in the randomized event
    Int_t         fNSkipLeadingCone;      // number of leading jets to be for the random cones
    Int_t         fNRandomCones;          // number of generated random cones
    Bool_t        randCone_pos;           // use fixed position for random cones
    Double_t      randCone_Eta;           // eta for random Cone at fixed position
    Double_t      randCone_Phi;           // phi for random Cone at fixed position
    Int_t         fNHBOM;                 // number of detector runs
    Float_t       fTrackEtaWindow;        // eta window used for corraltion plots between rec and gen 
    Float_t       fTrackPtCut;            // minimum track pt to be accepted
    Float_t       fJetOutputMinPt;        // minimum p_t for jets to be written out
    Float_t       fMaxTrackPtInJet;       // maximum track pt within a jet for flagging...
    //    Float_t       fJetTriggerPtCut;       // minimum jwt pT for AOD to be written
    Float_t       fVtxZCut;               // zvtx cut
    Float_t       fVtxR2Cut;              // R vtx cut (squared) 
    Float_t       fCentCutUp;             // upper limit on centrality
    Float_t       fCentCutLo;             // lower limit on centrality
    // output configurartion
    TString       fNonStdBranch;      // the name of the non-std branch name, if empty no branch is filled
    TString       fBackgroundBranch;  // name of the branch used for background subtraction
    TString       fNonStdFile;        // The optional name of the output file the non-std branch is written to

    //Detector level effects
    TProfile *fMomResH1; // Momentum resolution from TrackQA Hybrid Category 1
    TProfile *fMomResH2; // Momentum resolution from TrackQA Hybrid Category 2
    TProfile *fMomResH3; // Momentum resolution from TrackQA Hybrid Category 3
    TF1 *fMomResH1Fit; //fit
    TF1 *fMomResH2Fit; //fit
    TF1 *fMomResH3Fit; //fit

    TH1      *fhEffH1;        // Efficiency for Spectra Hybrid Category 1
    TH1      *fhEffH2;        // Efficiency for Spectra Hybrid Category 2
    TH1      *fhEffH3;        // Efficiency for Spectra Hybrid Category 3
    Bool_t    fUseTrMomentumSmearing;     // Apply momentum smearing on track level
    Bool_t    fUseDiceEfficiency;         // Apply efficiency on track level by dicing

    // Fast jet
    Double_t fRparam;                  // fastjet distance parameter
    fastjet::JetAlgorithm fAlgorithm; //fastjet::kt_algorithm
    fastjet::Strategy fStrategy;  //= fastjet::Best;
    fastjet::RecombinationScheme fRecombScheme; // = fastjet::BIpt_scheme;
    fastjet::AreaType fAreaType;  // fastjet area type
    Double_t fGhostArea;          // fasjet ghost area
    Int_t fActiveAreaRepeats;     // fast jet active area repeats
    Double_t fGhostEtamax;        // fast jet ghost area

    Double_t background; //background rho in the event

    TClonesArray  *fTCARandomConesOut;    //! TCA of output jets in randomized event

    TRandom3*     fRandom;   //! random number generator
    TProfile*     fh1Xsec;   //! pythia cross section and trials
    TH1F*         fh1Trials; //! trials are added
    TH1F*         fh1PtHard;  //! Pt har of the event...       
    TH1F*         fh1PtHardNoW;  //! Pt har of the event without weigt       
    TH1F*         fh1PtHardTrials;  //! Number of trials 

    TH1F*         fh1Nch;            //! charged particle mult
    TH1F*         fh1CentralityPhySel;          // ! centrality of anaylsed events 
    TH1F*         fh1Centrality;                // ! centrality of selected events 
    TH1F*         fh1DeltapT;        // pT of random Cone - background energy
    TH1F*         fh1Rho;            //background rho
    TH1F*         fh1RhoSigma;       //fluctuation of the background
    TH1F*         fh1PtRandCone;     //pT of random Cone

    TH1F*         fh1efficiencyPt;          //here efficiency is stored
    TH2D*         fh2efficiencyPhi;         //here efficiency is stored

    TH1F*         fh1ZPhySel;          // ! zvtx of anaylsed events 
    TH1F*         fh1Z;                // ! zvtx of selected events 


    TList *fHistList; //!leading tracks to be skipped in the randomized event Output list
   

    ClassDef(AliAnalysisTaskJetHBOM, 1) 
};
 
#endif
