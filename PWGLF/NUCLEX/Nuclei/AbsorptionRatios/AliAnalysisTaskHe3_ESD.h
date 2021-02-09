#ifndef AliAnalysisTaskHe3_ESD_cxx
#define AliAnalysisTaskHe3_ESD_cxx

// example of analysis task

class TH1F;
class TH1I;
class TH1;
class AliVEvent;
class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class TProfile;
class TList;
class AliPIDResponse; // main class for PID analysis

#include "AliAnalysisTaskSE.h"
#include <iostream>

using namespace std;

class AliAnalysisTaskHe3_ESD : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskHe3_ESD();
	AliAnalysisTaskHe3_ESD(const char *name);
	virtual ~AliAnalysisTaskHe3_ESD();
    
	// Double_t Psi_pair(AliESDtrack *neg, AliESDtrack *pos);
    void   UserCreateOutputObjects();
    void   UserExec(Option_t *option);
    void   Terminate(Option_t *);
    
    void     SetTrackCuts(AliESDtrackCuts* const cuts) { fTrackCuts = cuts; }
    void     SetFilterBit(Int_t fBit)            {fFilterBit = fBit;}
    void     SetLowPCut(Float_t lowP)            {fLowPCut = lowP;}
    void     SetHighPCut(Float_t highP)          {fHighPCut = highP;}
    void     SetEtaCut(Float_t EtaCut)           {fEtaCut = EtaCut;}
    void     SetMinNITSCl(Int_t MinNclIts)       {fMinClIts = MinNclIts;}
    void     SetMaxDCAxyPreCut(Float_t maxDCAxy) {fMaxDCAxyCut = maxDCAxy;}
    void     SetMaxDCAxyFinal(Float_t maxDCAxy)  {fMaxDCAxyAna = maxDCAxy;}
    void     SetMaxDCAz(Float_t maxDCAz)         {fMaxDCAz = maxDCAz;}
    // PID
    void     SetMaxTPCnSigma(Float_t maxTPCnSigma)       {fMaxTPCnSigma = maxTPCnSigma;}
    void     SetMaxTOFnSigma(Float_t maxTOFnSigma)       {fMaxTOFnSigma = maxTOFnSigma;}
    void     SetUseTOFPidCut(Bool_t useTOFPidCut)        {fUseTOFPidCut = useTOFPidCut;}
    void     SetMomForTOFanaProt(Float_t momTOFprot)     {fMomTOFProt = momTOFprot;}
    void     SetMomForTOFanaDeut(Float_t momTOFdeut)     {fMomTOFDeut = momTOFdeut;}
	void	 SetAnalyseAllParticles(Bool_t doAnalyse)	 {kAnalyseAllParticles = doAnalyse;} // Call this function with argument kTRUE in the AddTask macro in order to run the analysis for protons and deuterons as well. 
    
    void     CreateHistosTrack(vector<TH1*> &histos);
    void     FillHistosTrack(vector<TH1*> &histos, AliESDtrack *track);
    Float_t  GetTOFBeta(AliESDtrack *);
    Float_t  GetMass2TOF(Float_t, AliESDtrack *);
    
private:
    // transient members are not streamed
    AliVEvent   *fVEvent;           //! vEvent
    AliAODEvent	*fAODEvent;         //! AOD event
	AliESDEvent *fESDEvent;			//! ESD event
    TList       *fOutputList;       //! Output list
    TList       *fOutputEvent;      //! Event list
    TList       *fOutputProtons;    //! protons list
    TList       *fOutputAProtons;   //! anti-protons list
    TList       *fOutputDeuterons;  //! deuterons list
    TList       *fOutputADeuterons; //! anti-deuterons list
	TList		*fOutputHe3;		//! helium 3 list
	TList		*fOutputAHe3;		//! anti-helium 3 list
    
    TProfile    *fHistTrackCuts;    //! track cuts config
    
    // Event histograms
    TH1I        *fEventStat;       //! Event statistics
    TH1F        *fVertexZ;         //! Vertex Z coordinate
    TH1F        *fVtxContrib;      //! Vertex contributors
    TH1F        *fSPDVtxResol;     //! SPD vertex resolution
    TH2F        *fVtxDisplacement; //! displacement btw track and spd vertices
    TH1F        *fMultV0;          //! V0 multiplicity
    TH1F        *fCentV0M;         //! V0M centrality
	TH1F		*fCentV0MZoomed;	//! V0M centrality Zoomed to 0-1%
	TH1F		*fNch;				//! Number of charged tracks at mid rapidity per event histogram
	TH1F		*fNchHeader;		//! -||- but from AliAODheader
	TH2F		*fHe3TOFCutVariation; //!
	TH2F		*fAHe3TOFCutVariation; //!

	// TNtuple
	
	TNtuple		*fNtupleHe3;		//! Ntuple for He3 Candidates
	TNtuple		*fNtupleAHe3;		//! Ntuple for AHe3 Candidates
    
    // track histograms
	//
	// to do: change the comment to not stream (//!) and re run to investigate why this wa set up in this way. 
    std::vector<TH1*> fHistsProton; // vector of proton histograms
    std::vector<TH1*> fHistsAProton; // vector of anti-proton histograms
    
    std::vector<TH1*> fHistsDeuteron; // vector of proton histograms
    std::vector<TH1*> fHistsADeuteron; // vector of anti-proton histograms

	std::vector<TH1*> fHistsHe3;		// vector of helium 3 histograms
	std::vector<TH1*> fHistsAHe3;		// vector of anti helium 3 histograms
    
    AliPIDResponse *fPIDResponse;   //! PID response object
    
    // persistent members are streamed (copied/stored)
    AliESDtrackCuts *fTrackCuts; // Track cuts
    Int_t         fFilterBit;
    Float_t       fLowPCut;
    Float_t       fHighPCut;
    Float_t       fEtaCut;
    Int_t         fMinClIts;
    Float_t       fMaxDCAxyCut;
    Float_t       fMaxDCAxyAna;
    Float_t       fMaxDCAz;
    // PID
    Float_t       fMaxTPCnSigma;
    Float_t       fMaxTOFnSigma;
    Bool_t        fUseTOFPidCut;
    Float_t       fMomTOFProt;
    Float_t       fMomTOFDeut;
	Bool_t		  kAnalyseAllParticles;		
    
    enum {kSelectedEvents=0, kINT7selected, kDAQincomplete, kV0timing, kClusterTrackletCut, kVertexZ, kVtxNcontrib, kSPDPileUp, kVtxDisplace, kVtxRes, kNbinsEvent};
    
	AliAnalysisTaskHe3_ESD(const AliAnalysisTaskHe3_ESD&); // not implemented
	AliAnalysisTaskHe3_ESD& operator=(const AliAnalysisTaskHe3_ESD&); // not implemented
    
	ClassDef(AliAnalysisTaskHe3_ESD, 1);
};

#endif
