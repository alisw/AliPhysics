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
		virtual void   SetDebug(Int_t debug = 0)        {fDebug = debug;}
		virtual void   SetAlgorithm(TString jf="ANTIKT"){JFAlg=jf;}
		virtual void   SetRadius(Float_t radius=0.4)    {Radius=radius;}
		virtual void   SetFilterMask(UInt_t filter=256) {Filtermask=filter;}
		virtual void   SetBackSubMode(Int_t backM=0)    {BackM=backM;}
		virtual void   SetTrackPtCut(Float_t tPtcut=0)  {TrackPtcut=tPtcut;}
		virtual void   SetSkipCone(Int_t skipCone=0)    {SkipCone=skipCone;}
		virtual void   SetMC(Bool_t ismc=true)          {IsMC=ismc;}
		virtual void   SetJetEScale  (Float_t JEScale=1.){JetEScale=JEScale;}
		virtual void   SetTrackEScale(Float_t TEScale=1.){TrackEScale=TEScale;}
		virtual void   FinishTaskOutput();

		//enum {kNPTBINS=10};

		// 0 all jets
		// 1 all jet in eta window
		// 2 all jets with partner
		// 3 all jets in eta window with partner
		// 4 all jets with partner in eta window
		//enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};



	private:
		AliAnalysisTaskJetHadronCorrelation(const AliAnalysisTaskJetHadronCorrelation &det); // not implemented
		AliAnalysisTaskJetHadronCorrelation& operator=(const AliAnalysisTaskJetHadronCorrelation &det); // not implemented
		Double_t DeltaPhi(Double_t phi1,Double_t phi2);

		Bool_t          fUseAODInput; // read jets from input AOD
		TString         fJetBranch;   // jet branch to read
		TString         fNonStdFile;

		AliAODEvent     *fAODIn;       // AOD event
		AliAODEvent     *fAODOut;      // AOD event
		AliAODExtension *fAODExtension;
		TString         JFAlg;
		Float_t         Radius;
		UInt_t          Filtermask;
		Int_t           BackM;
		Float_t         TrackPtcut;
		Int_t           SkipCone;
		Bool_t          IsMC;
		Float_t         JetEScale;
		Float_t         TrackEScale;


		Float_t         fxsec;
		Float_t         ftrial;

		TList           *fHistList; // Output list
    Int_t           fIfiles;//!count no. of files


		TH1F            *fH1Events;
		TProfile        *fH1Xsec;
		TH1F            *fH1Trials;

		TH1F            *fH1Track_pt          ;
		TH1F            *fH1Track_phi         ;
		TH1F            *fH1Track_eta         ;
		TH1F            *fH1MCTrack_pt        ;
		TH1F            *fH1MCTrack_phi       ;
		TH1F            *fH1MCTrack_eta       ;
		TH1F            *fH1MCPrimTrack_pt    ;
		TH1F            *fH1MCPrimTrack_phi   ;
		TH1F            *fH1MCPrimTrack_eta   ;
		TH1F            *fH1Jet_pt            ;
		TH1F            *fH1Jet_phi           ;
		TH1F            *fH1Jet_eta           ;
		TH1F            *fH1leadJet_pt        ;
		TH1F            *fH1leadJet_pt_dijet  ;
		TH1F            *fH1subJet_pt_dijet   ;
		TH1F            *fH1JetMC_pt          ;
		TH1F            *fH1leadJetMC_pt      ;
		TH1F            *fH1leadJetMC_pt_dijet;
		TH1F            *fH1subJetMC_pt_dijet ;
		TH2F            *fH2JetsJet_dphi      ;
		TH2F            *fH2JetsJet_deta      ;
		TH2F            *fH2JetsJet_Aj        ;
		TH2F            *fH2JetsJet_pt        ;
		TH2F            *fH2JetsJetMC_dphi    ;
		TH2F            *fH2JetsJetMC_deta    ;
		TH2F            *fH2JetsJetMC_Aj      ;
		TH2F            *fH2JetsJetMC_pt      ;

		TH2F         *fH2Mult_Mtrack          ;
		TH2F         *fH2Mult_Mlead           ;
		TH2F         *fH2Mult_Mjet            ;
		TH2F         *fH2Mult_Njet            ;
		TH2F         *fH2Mult_Aj              ;
		TH2F         *fH2Mlead_Aj             ;
		TH2F         *fH2Jet_pt_Mlead         ;
		TH2F         *fH2Jet_pt_Munder        ;

		TH2F         *fH2leadJetMCptResolution ;
		TH2F         *fH2TrackMCptResolution   ;
		TH2F         *fH2TrackMCptEfficiency   ;
		TH2F         *fH2AjCorrelation_MCRec   ;
		TH2F         *fH2MleadCorrelation_MCRec;

		TH1F         *fH1ndiJ_ediv                     [5];
		TH1F         *fH1Aj                            [5];
		TH1F         *fH1Mlead                         [5];

		TH1F         *fH1leadJetMC_dphiResolution      [5];
		TH1F         *fH1subJetMC_dphiResolution       [5];
		TH1F         *fH1leadJetMC_Efficiency          [5];
		TH1F         *fH1subJetMC_Efficiency           [5];

		TH1F         *fH1JetHadron_dphi_ediv             [5][5];
		TH1F         *fH1JetHadron_dphi_tptweight_ediv   [5][5];
		TH1F         *fH1JetHadron_dphi_tJptweight_ediv  [5][5];
		TH1F         *fH1JetHadronMC_dphi_ediv           [5][5];
		TH1F         *fH1JetHadronMC_dphi_tptweight_ediv [5][5];
		TH1F         *fH1JetHadronMC_dphi_tJptweight_ediv[5][5];
		TH1F         *fH1JetHadronMCPrim_dphi_ediv           [5][5];
		TH1F         *fH1JetHadronMCPrim_dphi_tptweight_ediv [5][5];
		TH1F         *fH1JetHadronMCPrim_dphi_tJptweight_ediv[5][5];
		//TH1F         *fH1JetHadronMCIdeal_dphi_ediv             [5][5];
		//TH1F         *fH1JetHadronMCIdeal_dphi_tptweight_ediv   [5][5];
		//TH1F         *fH1JetHadronMCIdeal_dphi_tJptweight_ediv  [5][5];

		TH1F         *fH1ndiJ_2040Mlead                         [3];
		TH1F         *fH1ndiJ_2040Aj                            [3];
		TH1F         *fH1JetHadron_dphi_tptweight2040_Mleaddep  [3][5];
		TH1F         *fH1JetHadron_dphi_tptweight2040_Ajdep     [3][5];
		TH1F         *fH1JetHadronMC_dphi_tptweight2040_Mleaddep[3][5];
		TH1F         *fH1JetHadronMC_dphi_tptweight2040_Ajdep   [3][5];
		TH1F         *fH1JetHadronMCPrim_dphi_tptweight2040_Mleaddep[3][5];
		TH1F         *fH1JetHadronMCPrim_dphi_tptweight2040_Ajdep   [3][5];
		//TH1F         *fH1JetHadronMCIdeal_dphi_tptweight2040_Mleaddep[3][5];
		//TH1F         *fH1JetHadronMCIdeal_dphi_tptweight2040_Ajdep   [3][5];

		ClassDef(AliAnalysisTaskJetHadronCorrelation, 17); // Analysis task for JetHadronCorrelation
};

#endif
