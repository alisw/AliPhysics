#ifndef ALIANALYSISTASKJETCLUSTER_H
#define ALIANALYSISTASKJETCLUSTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// task used for comapring different jets D parmaters from fastjet 
// *******************************************

#include  "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

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
class TH1F;
class TH3F;
class TProfile;



class AliAnalysisTaskJetCluster : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetCluster();
    AliAnalysisTaskJetCluster(const char* name);
    virtual ~AliAnalysisTaskJetCluster() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();

    virtual void SetUseGlobalSelection(Bool_t b){fUseGlobalSelection = b;}
    virtual void SetAODTrackInput(Bool_t b){fUseAODTrackInput = b;}
    virtual void SetAODMCInput(Bool_t b){fUseAODMCInput = b;}
    virtual void SetRecEtaWindow(Float_t f){fRecEtaWindow = f;}
    virtual void SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
    virtual void SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
    virtual void SetFilterMask(UInt_t i){fFilterMask = i;}

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

    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackUndef = 0, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
    enum {kMaxJets = 4};
    enum {kMaxCorrelation =  3};
    enum {kMaxRadius =       5};
    

 private:

    AliAnalysisTaskJetCluster(const AliAnalysisTaskJetCluster&);
    AliAnalysisTaskJetCluster& operator=(const AliAnalysisTaskJetCluster&);

    Int_t GetListOfTracks(TList *list,Int_t type);

    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD

    Bool_t        fUseAODTrackInput;      // take track from input AOD not from ouptu AOD
    Bool_t        fUseAODMCInput;         // take MC from input AOD not from ouptu AOD
    Bool_t        fUseGlobalSelection;    // Limit the eta of the generated jets
    UInt_t        fFilterMask;             // filter bit for slecected tracks
    Int_t         fTrackTypeRec;          // type of tracks used for FF 
    Int_t         fTrackTypeGen;          // type of tracks used for FF 
    Float_t       fAvgTrials;             // Average nimber of trials
    Float_t       fExternalWeight;        // external weight
    Float_t       fRecEtaWindow;          // eta window used for corraltion plots between rec and gen 
    // Fast jet
    Double_t fRparam;
    fastjet::JetAlgorithm fAlgorithm; //fastjet::kt_algorithm
    fastjet::Strategy fStrategy;  //= fastjet::Best;
    fastjet::RecombinationScheme fRecombScheme; // = fastjet::BIpt_scheme;
    fastjet::AreaType fAreaType; 

    TProfile*     fh1Xsec;   // pythia cross section and trials
    TH1F*         fh1Trials; // trials are added
    TH1F*         fh1PtHard;  // Pt har of the event...       
    TH1F*         fh1PtHardNoW;  // Pt har of the event without weigt       
    TH1F*         fh1PtHardTrials;  // Number of trials 

    TH1F*         fh1NJetsRec; // number of reconstructed jets
    TH1F*         fh1NConstRec;// number of constiutens in leading jet
    TH1F*         fh1NConstLeadingRec;// number of constiutens in leading jet
    TH1F*         fh1PtJetsRecIn;  // Jet pt for all jets
    TH1F*         fh1PtJetsLeadingRecIn;  // Jet pt for all jets
    TH1F*         fh1PtJetConstRec;// pt of constituents
    TH1F*         fh1PtJetConstLeadingRec;// pt of constituents
    TH1F*         fh1PtTracksRecIn;  // track pt for all tracks
    TH1F*         fh1PtTracksLeadingRecIn;  // track pt for all tracks

    // Randomized track histos
    TH1F*         fh1NJetsRecRan; // number of reconstructed jets from randomized
    TH1F*         fh1NConstRecRan;// number of constiutens in leading jet
    TH1F*         fh1PtJetsLeadingRecInRan;  // Jet pt for all jets
    TH1F*         fh1NConstLeadingRecRan;// number of constiutens in leading jet
    TH1F*         fh1PtJetsRecInRan;  // Jet pt for all jets

    TH1F*         fh1PtTracksGenIn;  // track pt for all tracks
    TH1F*         fh1Nch;            // charged particle mult

    TH2F*         fh2NRecJetsPt;            // Number of found jets above threshold
    TH2F*         fh2NRecTracksPt;          // Number of found tracks above threshold
    TH2F*         fh2NConstPt;           // number of constituents vs. pt
    TH2F*         fh2NConstLeadingPt;           // number of constituents vs. pt
    TH2F*         fh2JetPhiEta;             // jet phi eta
    TH2F*         fh2LeadingJetPhiEta;      // leading jet phi eta
    TH2F*         fh2JetEtaPt;              // leading jet eta
    TH2F*         fh2LeadingJetEtaPt;              // leading jet eta
    TH2F*         fh2TrackEtaPt;              // track eta
    TH2F*         fh2LeadingTrackEtaPt;       // leading track eta
    TH2F*         fh2JetsLeadingPhiEta;     // jet phi eta
    TH2F*         fh2JetsLeadingPhiPt;      // jet correlation with leading jet
    TH2F*         fh2TracksLeadingPhiEta;   // track correlation with leading track
    TH2F*         fh2TracksLeadingPhiPt;    // track correlation with leading track
    TH2F*         fh2TracksLeadingJetPhiPt; // track correlation with leading Jet
    TH2F*         fh2JetsLeadingPhiPtW;      // jet correlation with leading jet
    TH2F*         fh2TracksLeadingPhiPtW;   // track correlation with leading track
    TH2F*         fh2TracksLeadingJetPhiPtW; // track correlation with leading Jet
    TH2F*         fh2NRecJetsPtRan;            // Number of found jets above threshold
    TH2F*         fh2NConstPtRan;           // number of constituents vs. pt
    TH2F*         fh2NConstLeadingPtRan;           // number of constituents vs. pt
    TH2F*         fh2PtNch;               // p_T of cluster vs. multiplicity,
    TH2F*         fh2PtNchRan;            // p_T of cluster vs. multiplicity,random
    TH2F*         fh2PtNchN;               // p_T of cluster vs. multiplicity, weigthed with constituents
    TH2F*         fh2PtNchNRan;            // p_T of cluster vs. multiplicity, weigthed with constituents random
    TH2F*         fh2TracksLeadingJetPhiPtRan; // track correlation with leading Jet
    TH2F*         fh2TracksLeadingJetPhiPtWRan; // track correlation with leading Jet
    TList *fHistList; // Output list
   

    ClassDef(AliAnalysisTaskJetCluster, 3) // Analysis task for standard jet analysis
};
 
#endif
