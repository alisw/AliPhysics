#ifndef ALIANALYSISTASKEMCALJETBJETTAGGINGIP_H
#define ALIANALYSISTASKEMCALJETBJETTAGGINGIP_H

// $Id$
class TF1;
class TH1;
class TH2;


class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliHFJetTaggingIP;
class AliAODRecoDecayHF3Prong;
class AliVertexerTracks;
class AliAODRecoDecay;
class AliAODVertex;
class AliESDv0;
class AliAODv0;
class AliHFJetsTagging;
class AliRDHFJetsCuts;
class AliPIDResponse;
class AliPIDCombined;
class AliAODMCHeader;
class AliAnalysisUtils;
class TRandom3;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliHFJetsTagging.h"
#include <TTree.h>

class AliAnalysisTaskEmcalJetBJetTaggingIP : public AliAnalysisTaskEmcalJet
{
public:
	// ENUM definitions
	enum MatchingType{kNoMatching = 0,kGeometrical = 1};
	enum ExendedEFlavourTag { kBeautyJet = 1 << 7, kCharmJet = 1 << 8, kLFgJet = 1 << 9 };
	enum EQualityClass { kQtyStandard, kQtyVeryGood, kQtyGood, kQtyMedium, kQtyBad };

	Double_t                    fMatchingPar1;                           // matching parameter for jet1-jet2 matching
	Double_t                    fMatchingPar2;                           // matching parameter for jet2-jet1 matching
	Double_t                    fMinJetMCPt;                             // minimum jet MC pt

	//Constructors and standard function definitions
	AliAnalysisTaskEmcalJetBJetTaggingIP();
	AliAnalysisTaskEmcalJetBJetTaggingIP(const char* name);
	virtual ~AliAnalysisTaskEmcalJetBJetTaggingIP(){;}
	void UserCreateOutputObjects();
	virtual Bool_t Run();
	//Standard setters
	void SetMC(Bool_t val = kTRUE){fIsMC = val;};
	void SetUseAtlasSignTC(Bool_t val = kTRUE){fUseAtlasSignTCCalculation = kTRUE; };
	void SetEventSelectionMethod(int i = 0){fUseEventSelection = i;};
	void SetJetSelectionMethod(int i = 0){fUseJetSelection = i;};
	void SetJetTaggerMCMethod(int i = 0){fUseMCTagger = i;};
	void SetUseMVRejectionUtils(Bool_t val = kTRUE){fUseMVRejectionTools = val; };
	void SetUseFunctionalBunchCut(bool val, TF1* f){fUseFunctionalBunchCut = val;fFunctionalBunchCut = f;};
	//Standard getters
	virtual AliRDHFJetsCuts* GetJetCutsHF();
	virtual AliHFJetTaggingIP* GetTaggerIP(){return fTrackCountingTagger;}
	virtual Bool_t RunQATracksEvent();
	virtual Bool_t RunQATracksJet(const AliEmcalJet* jet);
	void SetDoTrackQA(Bool_t val = kTRUE){	fIsTrackQA = val;};
	void SetDoTrackQAConstituent(Bool_t val = kTRUE){	fIsTrackQAConstituent = val;    };
	void SetUseCorrectedPt(Bool_t val = kTRUE){fUseCorrectedJetPt = val;};
	void SetDoBackgroundFluctuations(Bool_t val = kTRUE){fDoRandomCones = val; };
	void EnableVertexingMassTemplates(Bool_t val = 0){fVetexingMassFitTest = val;};
	void EnablePtRelTemplates(Bool_t val = 0){fVetexingPtRelTest = val;};
	void SetUseSignificance(Bool_t val = kTRUE){fUseSignificance = kTRUE;}
	void SetUse3DsIP(Bool_t val = kTRUE){fUse3DsIP = kTRUE;};
	void SetUseConversions(Bool_t val = kTRUE){fUseConversions = kTRUE;};
	virtual void SelectPtHardBin(const char* bin = "", const char* flavour = "", Bool_t val = kTRUE)
	{fSelectPtHard = val;
	fPtHardBin = bin;
	fPtHardFlavour = flavour;};
	virtual void SetPtRelMomentumCuts(Double_t JetPt = 5., Double_t ElectronPt = 1.){
		fpTRelSoftElectronPt = JetPt;
		fpTRelSoftElectronPt = ElectronPt;};
	virtual void SetMCMatchingRadius(double radius=2.0){fMCMatchingRadius = radius;};

protected:
	//Event selection methods
	virtual Bool_t IsEventSelected();
	virtual Bool_t IsEventSelectedLegacy(AliAODEvent* aev);
	virtual Bool_t IsEventSelectedpA(AliAODEvent* aev);
	virtual Bool_t IsEventSelectedpp(AliAODEvent* aev);
	virtual Bool_t IsEventSelectedHF(AliAODEvent* aev);
	//Jet selection methods
	virtual Bool_t IsJetSelected(const AliEmcalJet* jet);
	virtual Bool_t IsJetSelectedLegacy(const AliEmcalJet* jet);
	virtual Bool_t IsJetSelectedHF(const AliEmcalJet* jet);
	//Helpers
	virtual Bool_t AddTagJet(AliEmcalJet* jet);
	virtual Bool_t IsQuality(const AliAODTrack* track, EQualityClass qclass);
	virtual Bool_t IsV0DaughterRadius(const AliAODTrack* track, Double_t& Radius);
	virtual Bool_t IsElectronTPC(double nSigmaTPCelectron, double nSigmaTOFElectron);
	virtual Double_t GetDeltaPtRandomCone();
	virtual Double_t GetPtCorrected(const AliEmcalJet* jet);
	virtual Bool_t FindVertexNProngSimple(const AliEmcalJet* jet, AliAODVertex*& vtx, Int_t& nProng);
	virtual Double_t GetVertexInvariantMass(AliAODVertex* vtx, Double_t massParticle);
	virtual Double_t GetElectronPIDnSigmaTPC(const AliVTrack* track);
	virtual Double_t GetPionPIDnSigmaTPC(const AliVTrack* track);
	virtual Double_t GetPtRel(const AliVTrack* track, const AliEmcalJet* jet);
	virtual Bool_t TrackIsFromConversion(const AliVTrack* track, Bool_t useMC);
	//Jet matching
	virtual Bool_t DoJetMatching();
	virtual void DoJetLoop();
	virtual void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching);
	virtual void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
	virtual void FillTCTree(Double_t  sIP[3], Bool_t  isTagged[3],AliEmcalJet * jet,AliEmcalJet * jetMC,int flavour, int qtyclass );
	MatchingType                fMatching;                               // matching type
	// Calls for special purpose analysis
	virtual TList* AddHistsPtRelTemplates();
	virtual TList* AddHistsVtxMassTemplates();
	virtual void ProcessPtRelTemplateAnalysis(const AliEmcalJet* jet,
			double corrected_jet_pt,
			int mcflavour = 0,
			double n3value = -999.);
	virtual void ProcessVtxMassTemplateAnalysis(const AliEmcalJet* jet,
			double corrected_jet_pt,
			int mcflavour = 0,
			double n3value = -999.);
	virtual void JetTrackLoop(const AliEmcalJet* jet, Double_t n3tag, Int_t MCtag);
	virtual void FillVertexingHists(const AliEmcalJet* jet, AliAODVertex* vtx, Int_t nProng, Int_t tag = 0, Int_t hist = 0);
	virtual void AddHistTH1(TH1** hist,	const char* histname,const char* title,	const char* titlex,const char* titley,Int_t nBinsX,
			Double_t minX,Double_t maxX,Bool_t setSumw2,TList* container);
	virtual void AddHistTH2(TH2** hist,const char* histname,const char* title,const char* titlex,const char* titley,Int_t nBinsX,
			Double_t minX,Double_t maxX,Int_t nBinsY,Double_t minY,Double_t maxY,Bool_t setSumw2,
			TList* container);
	// Jet containers
	AliJetContainer* fJetsCont;              //!Jets
	AliJetContainer* fJetsContMC;            //!Jets MC
	AliParticleContainer* fTracksCont;       //!Tracks
	AliHFJetsTagging* fTaggingHFClass;       //!Helper
	AliHFJetTaggingIP* fTrackCountingTagger; //!Track counting class
	TRandom3* fRandom;                       //! Random cone input
	AliAnalysisUtils* fUtils;                //! AliAnlalysisUtils helper object;
	AliRDHFJetsCuts* fJetCutsHF;             // HF jet cuts
	AliPIDResponse* fPIDResponse;            //!Combined bayesian pid;
	AliPIDCombined* fPIDCombined;            //!Combined bayesian pid;
	AliAODMCHeader* fMCHeader;               //!AOD MC Header;
	Bool_t fSelectPtHard;                    //!Select MC AOD pT hard bin;
	const char* fPtHardBin;                  //!pT hard  bin selection string
	const char* fPtHardFlavour;              //!pT hard flavour selection string
private:
	TClonesArray* fMCparticles;         //! Monte Carlo particle stack
	AliAODVertex* fCurrentNProngVertex; //!NProngVertex for mass templates
	Int_t fCurrentNProngs;              // current number of prongs
	Bool_t fIsMC;                       // Is Monte Carlo Event
	Bool_t fIsTrackQA;                  // Run track QA for analysis note plots
	Bool_t fIsTrackQAConstituent;       // Run constituent track QA for analysis note plots
	Bool_t fUseCorrectedJetPt;          // Subtract average rho
	Bool_t fDoRandomCones;              // Get fluctuations
	Bool_t fVetexingMassFitTest;        //
	Bool_t fVetexingPtRelTest;          //
	Bool_t fUseAtlasSignTCCalculation;  //
	Bool_t fUseSignificance;
	Bool_t fUse3DsIP;
	Bool_t fUseFunctionalBunchCut;
	Bool_t fUseMVRejectionTools;
	Bool_t fUseConversions;
	TF1* fFunctionalBunchCut;                      // needs to be streamed
	Bool_t fCalculateSoftElectronVtxInvariantMass; // Calculate inv. vtx mass for vertices with electrons
	Int_t fUseEventSelection;                      // Event selection method
	Int_t fUseJetSelection;                        // Jet selection method
	Int_t fUseMCTagger;                            // MC tagger method selection
	Int_t fNumberOfsIPBins;                        //
	Double_t fSigmaTPCElectronLow;                 //
	Double_t fSigmaTPCElectronHigh;                //
	Double_t fSigmaITSElectronLow;                 //
	Double_t fSigmaITSElectronHigh;                //
	Double_t fSigmaTOFElectronLow;                 //
	Double_t fSigmaTOFElectronHigh;                //
	Double_t fpTRelSoftElectronPt;                 //
	Double_t fpTRelSoftElectronJetPt;              //
	Double_t fCurrentWeight;//
	Double_t fXsec;//
	Double_t fTrials;//
	Double_t fMCMatchingRadius;//

	TList* fCurrentTrackContainerSecVtx; //!
	TList * flist_module_eventselection;//!
	TList * flist_module_jets;//!
	TList * flist_module_trackcounting;//!
	TList * flist_module_unfolding;//!
	TList * flist_module_trackqajet;//!
	TList * flist_module_trackqa;//!
	TList * flist_module_mc;//!



	TTree * fTCTree; //!
	Float_t fCorrVariable_JetPT;
	Float_t fCorrVariable_sIP1;
	Float_t fCorrVariable_sIP2;
	Float_t fCorrVariable_sIP3;
	Float_t fCorrVariable_tPt1;
	Float_t fCorrVariable_tPt2;
	Float_t fCorrVariable_tPt3;


	Float_t fCorrVariable_JetPTMC;
	Float_t fCorrVariableMC_tPt1MC;
	Float_t fCorrVariableMC_tPt2MC;
	Float_t fCorrVariableMC_tPt3MC;

	Float_t fCorrVariableMC_tPt1MotherMC;
	Float_t fCorrVariableMC_tPt2MotherMC;
	Float_t fCorrVariableMC_tPt3MotherMC;
	Float_t fCorrVariableMC_tPt1MothersMotherMC;
	Float_t fCorrVariableMC_tPt2MothersMotherMC;
	Float_t fCorrVariableMC_tPt3MothersMotherMC;


	Int_t fCorrVariableMC_Pdg1;
	Int_t fCorrVariableMC_Pdg2;
	Int_t fCorrVariableMC_Pdg3;

	Int_t fCorrVariableMC_PdgMothersMother1;
	Int_t fCorrVariableMC_PdgMothersMother2;
	Int_t fCorrVariableMC_PdgMothersMother3;

	Int_t fCorrVariableMC_PdgMother1;
	Int_t fCorrVariableMC_PdgMother2;
	Int_t fCorrVariableMC_PdgMother3;
	Int_t 	fCorrVariable_flv;//[2]
	Int_t 	fCorrVariable_qty;//[2]

	// Histograms

	TH1* fhist_Events;                                             //! Event selection statistics
	TH2* fhist_Jets;                                               //!   Jet selection statistics
	TH2* fhist_MonteCarloFlavour;                                  //!
	TH2* fhists_SPD_cluster_vs_tracklet_correlation;               //! monitor_plot_for pileup rejection
	TH2* fhists_SPD_cluster_vs_tracklet_correlation_PostSelection; //! monitor_plot_for pileup rejection
	TH2* fhist_Tracks_Eta_Phi;                                     //! Regular hybrid tracks
	TH2* fhist_Tracks_Eta_Phi_Bit4;                                //! Regular ITS tracks
	TH2* fhist_Tracks_Eta_Phi_Bit9;                                //! Regular complementary tracks
	TH1* fhist_QualityClasses;                                     //! Tracks in quality class X
	TH2* fhist_QualityClasses_sIP[4][2];                           //![4][2]] prim. vs sec. track sip
	TH2* fhist_QualityClasses_Eta_Phi[4][2];                       //![4][2]] prim. vs sec. track eta phi
	TH2* fhist_TC_sIP_Pt[3][4][5];                                 //![3][4][5] Track counting output
	TH2* fhist_TC_sIP_zIP[3][4][5];                                 //![3][4][5] Track counting output
	TH2* fhist_TC_zIP_Pt[3][4][5];                                 //![3][4][5] Track counting output
	TH2* fhist_TC_Eta_Phi[3][4][5];                                //![3][4][5] Track counting eta
	TH2* fhist_TC_sIP_Pt_Conversions[3][4][5];  				  //![3][4][5] Track counting output e+/- from conversions
	TH2* fhist_TC_Eta_Phi_Conversions[3][4][5]; 				 //![3][4][5] Track counting eta phi e+/- from conversions
	TH2* fhist_Jet_Eta_Phi;                     					//! Jet eta-phi  distribution
	TH2* fhist_Jet_Nconst_Pt;                   //! Jet pT vs N constituents
	TH1* fhist_Jet_Pt;                          //! Jet pT  distribution
	TH1* fhist_Jet_Background_Fluctuation;      //! Background Fluctuations
	TH1* fhist_Rho;                             //! Average background momentum density

	TH1* fhist_parton_genjet_dR;  //! Delta R mc jet matched parton
	TH2* fhist_parton_genjet_pT;  //! p_T mc jet matched parton
	TH1* fhist_parton_genjet_Eta; //! \eta mc jet matched parton
	TH1* fhist_parton_genjet_Phi; //! \phi mc jet matched parton

	TH1* fhist_jet_pt_mc; //! momentum of inclusive monte carlo jets to calculate kinematic efficiency

	TH2* fhist_momentum_response[4][2]; //![4][2] momentum response matrix

	// Vertexing
	TH2* fHist_MassVsJetPtHE[4];      //! Inv Vtx Mass vs. Jet pT all  Prong;
	TH2* fHist_MassVsJetPtHP[4];      //! Inv Vtx Mass vs. Jet pT all  Prong;
	TH2* fHist_MassVsJetPtHEMC[4][4]; //! Inv Vtx Mass vs. Jet pT all  Prong;
	TH2* fHist_MassVsJetPtHPMC[4][4]; //! Inv Vtx Mass vs. Jet pT all  Prong;
	// PTrel electron method
	TH2* fHistTPCnSigmaElectron;     //! electron n Sigma TPC
	TH2* fHistTPCnSigmaPion;         //! pion n Sigma TPC
	TH2* fHistpTrelElectron[4];      //! electron pTrel for inclusive jets , n=3 tagged jets,
	TH2* fHistpTrelElectronMC[4][4]; //! electron pTrel for inclusive jets , n=3 tagged jets,


	AliAnalysisTaskEmcalJetBJetTaggingIP(const AliAnalysisTaskEmcalJetBJetTaggingIP&);            // not implemented
	AliAnalysisTaskEmcalJetBJetTaggingIP& operator=(const AliAnalysisTaskEmcalJetBJetTaggingIP&); // not implemented

	ClassDef(AliAnalysisTaskEmcalJetBJetTaggingIP, 111) // jet sample analysis task
};
#endif
