#ifndef ALIANALYSISTASKFASTEMBEDDING_H
#define ALIANALYSISTASKFASTEMBEDDING_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTaskSE.h"

class AliAODv0;
class AliAODVertex;
class AliPIDResponse;
class AliESDEvent;
class AliAODEvent;
class TTree;
class TFile;
class TChain;
class TObjArray;
class TObjString;
class TRandom3;
class TH1F;
class TH2F;
class TString;
class TList;
class TProfile;
class AliAODMCHeader;
class AliAODJet;

class AliAnalysisTaskFastEmbedding : public AliAnalysisTaskSE {

public:
   
   AliAnalysisTaskFastEmbedding();
   AliAnalysisTaskFastEmbedding(const char *name);
   AliAnalysisTaskFastEmbedding(const AliAnalysisTaskFastEmbedding &copy);
   AliAnalysisTaskFastEmbedding& operator=(const AliAnalysisTaskFastEmbedding &o);
   virtual ~AliAnalysisTaskFastEmbedding();

   virtual void UserCreateOutputObjects();
   virtual void LocalInit() { Init(); }
   virtual void Init();
   virtual Bool_t UserNotify();
   virtual void UserExec(Option_t*);
   virtual void Terminate(Option_t */*option*/);

   enum { kTrackUndef =0, kOnFly, kOnFlyPID, kOnFlydEdx, kOnFlyPrim, kOffl, kOfflPID, kOffldEdx, kOfflPrim };  
   enum { kK0, kLambda, kAntiLambda }; 


   void SetAODPath(TString path) {fAODPath = path;}
   void SetJetFriends(TString str){fStrJetFriends = str;}
   void SetArrayOfAODPaths(TObjArray* arr) {fAODPathArray = arr;}
   void SetArrayOfAODEntries(TArrayI* arr) {fAODEntriesArray = arr;}
   void SetAODEntriesSum(Int_t i){ fAODEntriesSum = i;}
   void SetAODEntriesMax(Int_t i){ fAODEntriesMax = i;}
   
   virtual void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; }
   virtual void     SetMinContribVtx(Int_t n) { fMinContribVtx = n; }
   virtual void     SetVtxZMin(Float_t z) { fVtxZMin = z; }
   virtual void     SetVtxZMax(Float_t z) { fVtxZMax = z; }
   virtual void     SetEvtClassMin(Int_t evtClass) { fEvtClassMin = evtClass; }
   virtual void     SetEvtClassMax(Int_t evtClass) { fEvtClassMax = evtClass; }
   virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
   virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
   virtual void     SetNInputTracksMin(Int_t nTr) { fNInputTracksMin = nTr; }
   virtual void     SetNInputTracksMax(Int_t nTr) { fNInputTracksMax = nTr; }
   
   void SetTrackBranch(TString name) {fTrackBranch = name;}
   void SetMCparticlesBranch(TString name) {fMCparticlesBranch = name;}
   void SetJetBranch(TString name) {fJetBranch = name;}

   void SetEmbedMode(Int_t m) {fEmbedMode = m;}
   void SetEPMode(Int_t m)  { fEPMode = m;}
   Int_t GetEmbedMode() const {return fEmbedMode;} 
   void SetEvtSelecMode(Int_t s) {fEvtSelecMode = s;}
   Int_t GetEvtSelecMode() const {return fEvtSelecMode;}

   void SetEvtSelJetPtRange(Float_t minPt, Float_t maxPt) {fEvtSelMinJetPt = minPt; fEvtSelMaxJetPt = maxPt;}
   void SetEvtSelJetEtaRange(Float_t minEta, Float_t maxEta) {fEvtSelMinJetEta = minEta; fEvtSelMaxJetEta = maxEta;}
   void SetEvtSelJetPhiRange(Float_t minPhi, Float_t maxPhi) {fEvtSelMinJetPhi = minPhi; fEvtSelMaxJetPhi = maxPhi;}
   void SetEvtSelJetMinLConstPt(Float_t minLPt)              {fEvtSelJetMinLConstPt = minLPt;}
   void SetEffExtra(Float_t effextra) {fExtraEffPb = effextra;}   
   void SetDiceMapEff(Bool_t dice) {fDiceMapEff = dice;}
   void SetLoadDiceTrackMapRootFile(TString path="$ALICE_ROOT/OADB/PWGJE/Efficiency/Efficiency_LHC11a2aj_Cent0_v1.root") {fPathDiceTrackMap=path;}
   void SetToyNumberOfTrackRange(Int_t minN = 1, Int_t maxN = 1){ fToyMinNbOfTracks = minN, fToyMaxNbOfTracks = maxN; }
   void SetToyTrackRanges(Double_t minPt = 50., Double_t maxPt = 50., Double_t ptDistr=0,
   Double_t minEta = -.5, Double_t maxEta = .5,
   Double_t minPhi = 0., Double_t maxPhi = 2*TMath::Pi())
   {
      fToyMinTrackPt = minPt; fToyMaxTrackPt = maxPt; fToyDistributionTrackPt = ptDistr;
      fToyMinTrackEta = minEta; fToyMaxTrackEta = maxEta;
      fToyMinTrackPhi = minPhi; fToyMaxTrackPhi = maxPhi;}
   void SetToyFilterMap(UInt_t f) {fToyFilterMap = f;}
   void SetTrackFilterMap(UInt_t f) {fTrackFilterMap = f;}

   //AZ
  
   virtual void SetK0Type(Int_t i){ fK0Type = i; }
   virtual void SetLaType(Int_t i){ fLaType = i; }
   virtual void SetALaType(Int_t i){ fALaType = i; }
   virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }
   virtual void   SetQAMode(Bool_t qa) { fQAMode = qa; }

   Float_t  GetFFRadius() const { return fFFRadius; }
   Bool_t DaughterTrackCheck(AliAODv0* v0, Int_t& nnum, Int_t& pnum);
   Bool_t ApplyV0Cuts(AliAODv0* v0, const Int_t type, const Int_t particletype, AliAODVertex* primVertex, AliAODEvent* aod);
   //void GetTracksInCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, const Double_t radius, Double_t& sumPt);


   //
 
 
  static Float_t GetPtHard(Bool_t bSet=kFALSE, Float_t newValue = 0.);
   
   virtual Int_t      GetPtHardBin(Double_t ptHard);
   virtual Bool_t     PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);
   Bool_t             HasMinLConstPt(AliAODJet* jet);
   virtual void  LoadDiceTrackMapRootFile();
   virtual void SetEfficiencyMap(TH1 *h1);

   // embedding modes
   enum {kAODFull=0, kAODJetTracks, kAODJet4Mom, kToyTracks};
   // event selection from AOD
   enum {kEventsAll=0, kEventsJetPt};


private:

   AliESDEvent    *fESD;        //! ESD object
   Bool_t fQAMode;              
   AliAODEvent    *fAODout;     //! AOD out
   AliAODEvent    *fAODevent;   //! AOD in
   AliAODEvent    *fAODeventJets; //! AOD in jets
   TTree          *fAODtree;    //! AODin tree
   TTree          *fAODtreeJets;//! AODin tree
   TFile          *fAODfile;    //! AODin file
   TFile          *fAODfileJets; //! AODin file jets (friends)
   AliAODMCHeader *mcHeader;    //! mc header
   Double_t fFFRadius;          // jet cone size


  Int_t fK0Type;      
  Int_t fLaType;      
  Int_t fALaType;   

  TList* fListK0s;
  TList* fListLa;
  TList* fListALa;
  /* TList* fListK0sCone;
  TList* fListLaCone;
  TList* fListALaCone;
  */

  AliPIDResponse *fPIDResponse;
	                           // PID AZ
   TRandom3       *rndm;        //! random nummer generator
   Int_t          fInputEntries; // total nb. of events (for this subjob)

   TObjArray *fAODPathArray;    // array of paths of AOD in file
   TArrayI   *fAODEntriesArray; // array of entries of AODs 
   TString    fAODPath;         // path of AOD in file
   TString    fStrJetFriends;   // AliAODFriends name 
   Int_t      fAODEntries;      // entries of AOD
   Int_t      fAODEntriesSum;   // sum of all entries of AODs
   Int_t      fAODEntriesMax;   // maximum entries of AODs

   AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline triggers to accept
   Int_t   fMinContribVtx; // minimum number of track contributors for primary vertex
   Float_t fVtxZMin;	      // lower bound on vertex z
   Float_t fVtxZMax;	      // upper bound on vertex z
   Int_t   fEvtClassMin;	// lower bound on event class
   Int_t   fEvtClassMax;	// upper bound on event class
   Float_t fCentMin;	      // lower bound on centrality
   Float_t fCentMax;	      // upper bound on centrality
   Int_t   fNInputTracksMin;  // lower bound of nb. of input tracks
   Int_t   fNInputTracksMax;  // upper bound of nb. of input tracks
   
   TString fTrackBranch;       // name of branch for extra tracks in AOD out
   TString fMCparticlesBranch; // name of branch for extra mcparticles in AOD out
   TString fJetBranch;         // name of branch for extra jets AOD in

   Int_t fFileId;   // nb. of file from the list
   Int_t fAODEntry; // entry of extra AOD
   Int_t fCountEvents; // count processed events in this file

   Int_t fEmbedMode;     // embedding mode: kAODFull=0, kAODJetTracks=1, kAODJet4Mom=2, kToyTracks=3
   Int_t fEPMode;        // event plane dependent embedding: fEPMode = 0 (not selected), fEPMode = 1 (in-plane), fEPMode = 2 (out-of-plane)
   Int_t fEvtSelecMode;  // event selection criterion: kEventsAll=0, kEventsJetPt=1

   // event selection from AOD
   Float_t fEvtSelMinJetPt;       // minimum pt of the leading jet
   Float_t fEvtSelMaxJetPt;       // maximum pt of the leading jet
   Float_t fEvtSelMinJetEta;      // minimum eta of the leading jet
   Float_t fEvtSelMaxJetEta;      // maximum eta of the leading jet
   Float_t fEvtSelMinJetPhi;      // minimum phi of the leading jet
   Float_t fEvtSelMaxJetPhi;      // maximum phi of the leading jet
   Float_t fEvtSelJetMinLConstPt; // minimum leading constituent pt of the leading jet
   Double_t fExtraEffPb;          //extra efficiency PbPb      
   Bool_t   fDiceMapEff;           //if kTRUE,dice tracks according to a map in a user histogram
 
   TString fPathDiceTrackMap;             // path to root file containing eff. maps
   TH1    *fhEffH1;               //histogram containing the efficiency map 

  // settings for toy "track generation"
   Int_t    fToyMinNbOfTracks;             // minimum nb. of tracks per event
   Int_t    fToyMaxNbOfTracks;             // maximum nb. of tracks per event
   Float_t  fToyMinTrackPt;                // minimum track pT
   Float_t  fToyMaxTrackPt;                // maximum track pT
   Float_t  fToyDistributionTrackPt;       // distribution of track pt
   Float_t  fToyMinTrackEta;               // minimum eta of tracks
   Float_t  fToyMaxTrackEta;               // maximum eta of tracks
   Float_t  fToyMinTrackPhi;               // minimum phi of tracks
   Float_t  fToyMaxTrackPhi;               // maximum phi of tracks
   UInt_t   fToyFilterMap;                 // filter map of tracks
   UInt_t   fTrackFilterMap;               // filter map of tracks for QA plots
   
   Int_t         fNPtHard;      // nb. of pT hard bins
   Double_t      fPtHard;       // pT hard
   Int_t         fPtHardBin;    // pT hard bin
   TClonesArray* fAODJets;      //! array of jets from aod
   Int_t         fNevents;      // number of events in aod
   Float_t       fXsection;     // average xsection of the event
   Float_t       fAvgTrials;    // average number of trials per event


   // histos
   TList *fHistList;          //  list of histograms
   TH1I  *fHistEvtSelection;  //! stastic of event selection
   TProfile *fh1Xsec;         //! cross-section
   TH1F  *fh1Trials;          //! nb. of trials (simulation)
   TH1F  *fh1TrialsEvtSel;    //! nb. of trials (event selection, e.g. jet pT)
   TH2F  *fh2PtHard;          //! pT hard bin
   TH2F  *fh2PtHardEvtSel;    //! pT hard bin (event selection)
   TH2F  *fh2PtHardTrials;     //! pT hard bin, weighted by nb. of trials
   
   // qa histos
   TH1F  *fh1TrackPt;         //! track pt
   TH2F  *fh2TrackEtaPhi;     //! track eta-phi
   TH1F  *fh1TrackN;          //! nb. of tracks
 
 
   //AZ:
   TH1F  *fh1V0Pt;         //! track pt
   TH2F  *fh2V0EtaPhi;     //! track eta-phi
   //  TH1F  *fh1V0N;          //! nb. of tracks
   TH1F  *fh1K0Pt;         //! track pt
   TH2F  *fh2K0EtaPhi;     //! track eta-phi
   //TH2F  *fh2K0sPtJetPtCone;     //! K0s candidate Pt (inside cone) vs Jet Pt
   TH1F  *fh1LaPt;         //! track pt
   TH2F  *fh2LaEtaPhi;     //! track eta-phi
   //TH2F  *fh2LaPtJetPtCone;     //! La candidate Pt (inside cone) vs Jet Pt
   TH1F  *fh1ALaPt;         //! track pt
   TH2F  *fh2ALaEtaPhi;     //! track eta-phi
   //TH2F  *fh2ALaPtJetPtCone;     //! ALa candidate Pt (inside cone) vs Jet Pt
   //

   TH1F  *fh1JetPt;           //! jet pt
   TH2F  *fh2JetEtaPhi;       //! jet eta-phi
   TH1F  *fh1JetRotPhi;       //! delta rotated jet to event plane phi
   TH1F  *fh1JetN;            //! nb. of jets
   TH1F  *fh1MCTrackPt;       //! MC track pt
   TH2F  *fh2MCTrackEtaPhi;   //! MC track eta-phi
   TH1F  *fh1MCTrackN;        //! nb. of MC tracks
   TH1I  *fh1AODfile;         //! used AOD files from AODPathArray
   TH2I  *fh2AODevent;        //! selected events in AODs
   TH1F  *fh1EP2;             //! 2nd order harmonic event plane ESD event 
   TH1F  *fh1EP3;             //! 3rd order harmonic event plane ESD event

   Int_t GetJobID();    // get job id (sub-job id on the GRID)
   Int_t SelectAODfile();
   Int_t OpenAODfile(Int_t trial = 0);


   ClassDef(AliAnalysisTaskFastEmbedding, 6);
};

#endif

