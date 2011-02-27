#ifndef ALIANALYSISTASKFASTEMBEDDING_H
#define ALIANALYSISTASKFASTEMBEDDING_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTaskSE.h"

class AliAODEvent;
class TTree;
class TFile;
class TChain;
class TObjArray;
class TObjString;
class TRandom3;
class TH1F;
class TH2F;

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
        virtual void UserExec(Option_t*);
	virtual void Terminate(Option_t */*option*/);

	void SetAODPath(TString path) {fAODPath = path;}
	void SetArrayOfAODPaths(TObjArray* arr) {fAODPathArray = arr;}
        void SetTrackBranch(TString name) {fTrackBranch = name;}
        void SetMCparticlesBranch(TString name) {fMCparticlesBranch = name;}
        void SetJetBranch(TString name) {fJetBranch = name;}

	void SetEmbedMode(Int_t m) {fEmbedMode = m;}
	Int_t GetEmbedMode() {return fEmbedMode;} 
	void SetEvtSelecMode(Int_t s) {fEvtSelecMode = s;}
	Int_t GetEvtSelecMode() {return fEvtSelecMode;}

	void SetEvtSelJetPtRange(Float_t minPt, Float_t maxPt) {fEvtSelMinJetPt = minPt; fEvtSelMaxJetPt = maxPt;}
	void SetEvtSelJetEtaRange(Float_t minEta, Float_t maxEta) {fEvtSelMinJetEta = minEta; fEvtSelMaxJetEta = maxEta;}
	void SetEvtSelJetPhiRange(Float_t minPhi, Float_t maxPhi) {fEvtSelMinJetPhi = minPhi; fEvtSelMaxJetPhi = maxPhi;}
	
        void SetToyNumberOfTrackRange(Int_t minN = 1, Int_t maxN = 1){ fToyMinNbOfTracks = minN, fToyMaxNbOfTracks = maxN; }
	void SetToyTrackRanges(Double_t minPt = 50., Double_t maxPt = 50., Double_t ptDistr=0,
		Double_t minEta = -.5, Double_t maxEta = .5,
		Double_t minPhi = 0., Double_t maxPhi = 2*TMath::Pi())
                {
	        fToyMinTrackPt = minPt; fToyMaxTrackPt = maxPt; fToyDistributionTrackPt = ptDistr;
	        fToyMinTrackEta = minEta; fToyMaxTrackEta = maxEta;
		fToyMinTrackPhi = minPhi; fToyMaxTrackPhi = maxPhi;}
        void SetToyFilterMap(UInt_t f) {fToyFilterMap = f;}


        // embedding modes
	enum {kAODFull=0, kAODJetTracks, kAODJet4Mom, kToyTracks};
	// event selection from AOD
	enum {kEventsAll=0, kEventsJetPt};


    private:

        AliAODEvent* fAODout;    //! AOD out
	AliAODEvent* fAODevent;  //! AOD in
	TTree* fAODtree;         //! AODin tree
        TFile* fAODfile;         //! AODin file
	TRandom3* rndm;           //! random nummer generator

	TObjArray* fAODPathArray;  // array of paths of AOD in file
	TString fAODPath;  // path of AOD in file

        TString fTrackBranch; // name of branch for extra tracks in AOD out
        TString fMCparticlesBranch; // name of branch for extra mcparticles in AOD out
        TString fJetBranch; // name of branch for extra jets AOD in

        Int_t fEntry; // entry of extra AOD

	Int_t fEmbedMode;
	Int_t fEvtSelecMode;

	// event selection from AOD
	Float_t fEvtSelMinJetPt;       // minimum pt of the leading jet
	Float_t fEvtSelMaxJetPt;       // maximum pt of the leading jet
        Float_t fEvtSelMinJetEta;      // minimum eta of the leading jet
        Float_t fEvtSelMaxJetEta;      // maximum eta of the leading jet
        Float_t fEvtSelMinJetPhi;      // minimum phi of the leading jet
        Float_t fEvtSelMaxJetPhi;      // maximum phi of the leading jet
        
         
        // settings for toy "track generation"
        Int_t    fToyMinNbOfTracks;             // minimum nb. of tracks per event
        Int_t    fToyMaxNbOfTracks;             // maximum nb. of tracks per event
	Double_t fToyMinTrackPt;                // minimum track pT
	Double_t fToyMaxTrackPt;                // maximum track pT
        Double_t fToyDistributionTrackPt;       // distribution of track pt
	Double_t fToyMinTrackEta;               // minimum eta of tracks
	Double_t fToyMaxTrackEta;               // maximum eta of tracks
	Double_t fToyMinTrackPhi;               // minimum phi of tracks
	Double_t fToyMaxTrackPhi;               // maximum phi of tracks
	UInt_t fToyFilterMap;                   // filter map of tracks


        // qa histos
        TList *fHistList;          //  list of histograms
        TH1F  *fh1TrackPt;         //! track pt
        TH2F  *fh2TrackEtaPhi;     //! track eta-phi
        TH1F  *fh1TrackN;          //! nb. of tracks
        TH1F  *fh1JetPt;           //! jet pt
        TH2F  *fh2JetEtaPhi;       //! jet eta-phi
        TH1F  *fh1JetN;            //! nb. of jets
        TH1F  *fh1MCTrackPt;       //! MC track pt
        TH2F  *fh2MCTrackEtaPhi;   //! MC track eta-phi
        TH1F  *fh1MCTrackN;        //! nb. of MC tracks
        TH1I  *fh1AODfile;         //! used AOD files from AODPathArray
		

	Int_t GetJobID();    // get job id (sub-job id on the GRID)
        Int_t SelectAODfile();
        Int_t OpenAODfile(Int_t trial = 0);


	ClassDef(AliAnalysisTaskFastEmbedding, 4);
};

#endif

