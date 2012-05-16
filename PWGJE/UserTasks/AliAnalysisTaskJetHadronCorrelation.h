#ifndef AliAnalysisTaskJetHadronCorrelation_cxx
#define AliAnalysisTaskJetHadronCorrelation_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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


#include "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF  
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>

class AliAnalysisTaskJetHadronCorrelation : public AliAnalysisTaskSE 
{
	public:
		AliAnalysisTaskJetHadronCorrelation();
		AliAnalysisTaskJetHadronCorrelation(const char *name);
		virtual ~AliAnalysisTaskJetHadronCorrelation() {;}

		// Implementation of interface methods
		virtual void   UserCreateOutputObjects();
		virtual void   Init();
		virtual Bool_t Notify();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);
		virtual void   SetDebug(Int_t debug = 0) {fDebug = debug;}
		virtual void   SetAlgorithm(TString jf="ANTIKT"){JFAlg=jf;}
		virtual void   SetRadius(Float_t radius=0.4){Radius=radius;}
		virtual void   SetFilterMask(UInt_t filter=256){Filtermask=filter;}
		virtual void   SetBackSubMode(Int_t backM=0){BackM=backM;}
		virtual void   SetTrackPtCut(Float_t tPtcut=0){TrackPtcut=tPtcut;}
		virtual void   SetSkipCone(Int_t skipCone=0){SkipCone=skipCone;}
		virtual void   SetMC(Bool_t ismc=true){IsMC=ismc;}
		virtual void   FinishTaskOutput();
		//		virtual Float_t GetTotalEvents(const char* currFile);
		virtual Double_t DeltaPhi(Double_t phi1,Double_t phi2);


		enum {kNPTBINS=10};

		// 0 all jets
		// 1 all jet in eta window
		// 2 all jets with partner
		// 3 all jets in eta window with partner
		// 4 all jets with partner in eta window
		enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};



	private:
		AliAnalysisTaskJetHadronCorrelation(const AliAnalysisTaskJetHadronCorrelation &det); // not implemented
		AliAnalysisTaskJetHadronCorrelation& operator=(const AliAnalysisTaskJetHadronCorrelation &det); // not implemented

		Int_t GetListOfJets(TList *list,TClonesArray *jarray,Int_t type);
		Bool_t JetSelected(AliAODJet *jet);

		Bool_t        fUseAODInput; // read jets from input AOD
		Bool_t        fFillAOD;     // option to fill AOD branch
		TString       fJetBranch;   // jet branch to read
		TString       fNonStdFile;

		AliAODEvent   *fAODIn;      // AOD event
		AliAODEvent   *fAODOut;      // AOD event
		AliAODExtension *fAODExtension;
		TString         JFAlg;
		Float_t         Radius;
		UInt_t          Filtermask;
		Int_t           BackM;
		Float_t         TrackPtcut;
		Int_t           SkipCone;
		Bool_t          IsMC;


		Float_t       fxsec;
		Float_t       ftrial;
		Float_t       fJetRecEtaWindow;       // eta window for rec jets
		Float_t       fMinJetPt;              // limits the jet p_T in addition to what already is done in the jet finder, this is important for jet matching for JF with lo threshold

		TList        *fHistList; // Output list
    Int_t        fIfiles;//!count no. of files


		TH1F         *fH1Events;
		TProfile     *fH1Xsec;
		TH1F         *fH1Trials;
		//for Reconstructed Jet (Data&MC)
		TH1F         *fH1JetAKT04_pt                ;
		TH1F         *fH1leadJetAKT04_pt            ;
		TH1F         *fH1leadJetAKT04_pt_dijet      ;
		TH1F         *fH1subJetAKT04_pt_dijet       ;
		TH2F         *fH2JetsJetAKT04_dphi          ;
		TH2F         *fH2JetsJetAKT04_deta          ;
		TH2F         *fH2JetsJetAKT04_Aj            ;
		TH2F         *fH2JetsJetAKT04_pt            ;

		TH1F         *fH1AKT04_ndiJ_ediv            [5];
		TH1F         *fH1JetHadronAKT04_dphi_ediv           [5][5];
		TH1F         *fH1JetHadronAKT04_dphi_tptweight_ediv [5][5];
		TH1F         *fH1JetHadronAKT04_dphi_tJptweight_ediv[5][5];

		//for Generated Jet (MC)

		TH1F         *fH1JetMCAKT04_pt                ;
		TH1F         *fH1leadJetMCAKT04_pt            ;
		TH1F         *fH1leadJetMCAKT04_pt_dijet      ;
		TH1F         *fH1subJetMCAKT04_pt_dijet       ;
		TH2F         *fH2JetsJetMCAKT04_dphi          ;
		TH2F         *fH2JetsJetMCAKT04_deta          ;
		TH2F         *fH2JetsJetMCAKT04_Aj            ;
		TH2F         *fH2JetsJetMCAKT04_pt            ;
		TH1F         *fH1leadJetMCAKT04_dphiResolution[5];
		TH1F         *fH1subJetMCAKT04_dphiResolution [5];

		ClassDef(AliAnalysisTaskJetHadronCorrelation, 13); // Analysis task for standard dijet analysis
};

#endif
