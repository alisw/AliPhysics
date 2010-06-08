#ifndef AliAnalysisTaskFragFuncBB_h
#define AliAnalysisTaskFragFuncBB_h

class AliESDEvent;
class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFragFuncBB : public AliAnalysisTaskSE {

	public:
	
	AliAnalysisTaskFragFuncBB();
	AliAnalysisTaskFragFuncBB(const char *name);
	virtual ~AliAnalysisTaskFragFuncBB();

	virtual void   UserCreateOutputObjects();
	virtual void   Init();
	virtual void   LocalInit() {Init();}
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t* );
	
	virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
	virtual void   SetUseGlobalSelection(Bool_t b){fUseGlobalSelection = b;}
	virtual void   SetAODJetInput(Bool_t b){fUseAODJetInput = b;}
	virtual void   SetAODTrackInput(Bool_t b){fUseAODTrackInput = b;}
	virtual void   SetAODMCInput(Bool_t b){fUseAODMCInput = b;}
	virtual void   SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
	virtual void   SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
	virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
	virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}
	

	private:

	// Consts
	enum{ fgkMaxJets=4 }; // max. nb. of stored jets
	
	enum {kTrackUndef=0, kTrackAOD, kTrackKineAll, kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};

	// 
	Int_t   GetListOfTracks(TList *list, Int_t type);

	//
	AliESDEvent* fESD;
	AliAODEvent* fAOD;
	AliMCEvent*  fMCEvent;
	
	TString fBranchRecJets; // branch name for reconstructed jets
	TString fBranchGenJets; // branch name for generated jets

        Int_t   fLeadingRecJet; //
        Int_t   fLeadingGenJet; //
	
	Int_t  	fTrackTypeRec;       // type of reconstructed tracks
	Int_t  	fTrackTypeGen;       // type of generated tracks
	UInt_t 	fFilterMask;	     // filter bit for selected tracks
	Bool_t 	fUseGlobalSelection; // use selection of JetHelperTask
	Bool_t  fUseAODJetInput;	 // take jet from input AOD not from ouptut AOD
	Bool_t  fUseAODTrackInput;	 // take track from input AOD not from ouptut AOD
	Bool_t  fUseAODMCInput;		 // take MC from input AOD not from ouptut AOD
	
	Float_t fJetRadius;  // jet radius, used for fragmentation function 
	Float_t	fEtaMaxPart; // eta acceptance range for particles
	Float_t	fEtaMaxJets; // eta acceptance range for jets

	
	// Histograms
	TList		*fHistList;		  // List of histos

        TH1F            *fh1_eventSelection;
        TH1F		*fh1_vertexNContributors;
        TH1F		*fh1_vertexZ;

        // all jets 
	TH1F		*fh1_recJets_Et;
	TH1F		*fh1_genJets_Et;
	TH2F		*fh2_recJets_EtaPhi;
	TH2F		*fh2_genJets_EtaPhi;

        // jets, without acceptance cuts
	TH1F		*fh1_recJetsWoC_Et;
	TH1F		*fh1_genJetsWoC_Et;
	TH2F		*fh2_recJetsWoC_EtaPhi;
	TH2F		*fh2_genJetsWoC_EtaPhi;

        // leading jets
	TH1F		*fh1_recLJets_Et;
	TH1F		*fh1_genLJets_Et;
	TH2F		*fh2_recLJets_EtaPhi;
	TH2F		*fh2_genLJets_EtaPhi;
	
        // FF, all jets
	TH2F		*fh2_recFF_JetEt;
	TH2F		*fh2_genFF_JetEt;
	TH2F		*fh2_recHumpBacked_JetEt;
	TH2F		*fh2_genHumpBacked_JetEt;


        // FF, all jets without acceptance cuts
	TH2F		*fh2_recFF_JetWoCEt;
	TH2F		*fh2_genFF_JetWoCEt;
	TH2F		*fh2_recHumpBacked_JetWoCEt;
	TH2F		*fh2_genHumpBacked_JetWoCEt;

        // FF, leading jets
	TH2F		*fh2_recFF_LJetEt;
	TH2F		*fh2_genFF_LJetEt;
	TH2F		*fh2_recHumpBacked_LJetEt;
	TH2F		*fh2_genHumpBacked_LJetEt;

 
        // particles
	TH1F		*fh1_recPart_Pt;
	TH1F		*fh1_genPart_Pt;
	TH2F		*fh2_recPart_EtaPhi;
	TH2F		*fh2_genPart_EtaPhi;

	TH1F		*fh1_recJetPart_Pt;
	TH1F		*fh1_genJetPart_Pt;
	TH2F		*fh2_recJetPart_EtaPhi;
	TH2F		*fh2_genJetPart_EtaPhi;

        TH2F            *fh2_recJetPart_RJetPt;
        TH2F            *fh2_genJetPart_RJetPt;

	TH1F		*fh1_recJetWoCPart_Pt;
	TH1F		*fh1_genJetWoCPart_Pt;
	TH2F		*fh2_recJetWoCPart_EtaPhi;
	TH2F		*fh2_genJetWoCPart_EtaPhi;

        TH2F            *fh2_recJetWoCPart_RJetPt;
        TH2F            *fh2_genJetWoCPart_RJetPt;

	TH1F		*fh1_recLJetPart_Pt;
	TH1F		*fh1_genLJetPart_Pt;
	TH2F		*fh2_recLJetPart_EtaPhi;
	TH2F		*fh2_genLJetPart_EtaPhi;

        TH2F            *fh2_recLJetPart_RJetPt;
        TH2F            *fh2_genLJetPart_RJetPt;

        TH3F		*fh3_recPart_EtaPhiPt;
	
	ClassDef(AliAnalysisTaskFragFuncBB, 1);
};

#endif
