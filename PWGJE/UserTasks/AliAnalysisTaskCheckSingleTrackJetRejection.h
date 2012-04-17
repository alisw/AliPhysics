#ifndef AliAnalysisTaskCheckSingleTrackJetRejection_cxx
#define AliAnalysisTaskCheckSingleTrackJetRejection_cxx

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

class AliAnalysisTaskCheckSingleTrackJetRejection : public AliAnalysisTaskSE 
{
	public:
		AliAnalysisTaskCheckSingleTrackJetRejection();
		AliAnalysisTaskCheckSingleTrackJetRejection(const char *name);
		virtual ~AliAnalysisTaskCheckSingleTrackJetRejection() {;}

		// Implementation of interface methods
		virtual void    UserCreateOutputObjects();
		virtual void    Init();
		virtual Bool_t  Notify();
		virtual void    UserExec(Option_t *option);
		virtual void    Terminate(Option_t *);
		virtual void    SetDebug(Int_t debug = 0) {fDebug = debug;}
		virtual void    SetAlgorithm(const char *jf="ANTIKT"){JFAlg=jf;}
		virtual void    SetRadius(Float_t radius=0.4){Radius=radius;}
		virtual void    SetFilterMask(UInt_t filter=256){Filtermask=filter;}
		virtual void    SetBackSubMode(Int_t backM=0){BackM=backM;}
		virtual void    SetTrackPtCut(Float_t tPtcut=0){TrackPtcut=tPtcut;}
		virtual void    SetSkipCone(Int_t skipCone=0){SkipCone=skipCone;}
		virtual void    SetMC(Bool_t ismc=kFALSE){IsMC=ismc;}
    virtual void    FinishTaskOutput();
    virtual Bool_t  PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);
    virtual Float_t GetTotalEvents(const char* currFile);
    virtual Double_t DeltaPhi(Double_t phi1,Double_t phi2);
    virtual Double_t DeltaR(Double_t phi1,Double_t phi2,Double_t eta1,Double_t eta2);


		enum {kNPTBINS=10};

		// 0 all jets
		// 1 all jet in eta window
		// 2 all jets with partner
		// 3 all jets in eta window with partner
		// 4 all jets with partner in eta window
		enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};



	private:
		AliAnalysisTaskCheckSingleTrackJetRejection(const AliAnalysisTaskCheckSingleTrackJetRejection &det); // not implemented
		AliAnalysisTaskCheckSingleTrackJetRejection& operator=(const AliAnalysisTaskCheckSingleTrackJetRejection &det); // not implemented

    Bool_t  JetSelected(AliAODJet *jet);

		Bool_t        fUseAODInput; // read jets from input AOD
		Bool_t        fFillAOD;     // option to fill AOD branch
		TString       fNonStdFile;
		TString       fJetBranch;   // jet branch to read

		AliAODEvent   *fAODIn;      // AOD event
		AliAODEvent   *fAODOut;      // AOD event
		AliAODExtension *fAODExtension;
		TString JFAlg;
		Float_t Radius;
		UInt_t Filtermask;
		Int_t BackM;
		Float_t TrackPtcut;
		Int_t SkipCone;
		Bool_t        IsMC;

		TList        *fHistList; // Output list
		Float_t       fxsec;
		Float_t       ftrial;
		Float_t       fJetRecEtaWindow;       // eta window for rec jets
		Float_t       fMinJetPt;              // limits the jet p_T in addition to what already is done in the jet finder, this is important for jet matching for JF with lo threshold

		TProfile     *fH1Xsec;
		TH1F         *fH1Trials;
		TH1F         *fH1Events; 	

		TH1F         *fH1jetMCAKT04_pt         [6];
		TH2F         *fH2jetMCAKT04_Jetpt_maxpt;
		TH1F         *fH1jetAKT04_pt           [6];     
		TH2F         *fH2jetAKT04_Jetpt_maxpt  ;

		TH2F         *fH2jetMCAKT04_Eratio       [6];
		TH1F         *fH1jetMCAKT04_match        [6];



		ClassDef(AliAnalysisTaskCheckSingleTrackJetRejection, 13); // Analysis task for standard dijet analysis
};

#endif
