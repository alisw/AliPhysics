#ifndef ALIANALYSISTASKEMCALJETBJETTAGGINGIP_H
#define ALIANALYSISTASKEMCALJETBJETTAGGINGIP_H

// $Id$
class TF1;
class TH1;
class TH2;
class THnSparse;

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

class AliAnalysisUtils;
class TRandom3;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliHFJetsTagging.h"
class AliAnalysisTaskEmcalJetBJetTaggingIP : public AliAnalysisTaskEmcalJet {
public:
  enum ExendedEFlavourTag{kBeautyJet = 1<<7,kCharmJet  = 1<<8,kLFgJet    = 1<<9};
  enum EQualityClass{kQtyStandard ,kQtyVeryGood,kQtyGood,kQtyMedium,kQtyBad};  
  
  AliAnalysisTaskEmcalJetBJetTaggingIP();
  AliAnalysisTaskEmcalJetBJetTaggingIP(const char *name);
  virtual ~AliAnalysisTaskEmcalJetBJetTaggingIP();

  void    UserCreateOutputObjects();
  virtual Bool_t Run();
  virtual Bool_t RunQATracksEvent();
  virtual Bool_t RunQATracksJet(const AliEmcalJet * jet);
  void   SetMC(Bool_t val = kTRUE){fIsMC=val;};
  void   SetDoTrackQA(Bool_t val = kTRUE){fIsTrackQA=val;};
  void   SetDoTrackQAConstituent(Bool_t val = kTRUE){fIsTrackQAConstituent=val;};
  void   SetUseCorrectedPt(Bool_t val = kTRUE){fUseCorrectedJetPt=val;};
  void   SetDoBackgroundFluctuations (Bool_t val = kTRUE){fDoRandomCones=val;};
  void   EnableVertexingMassTemplates(Bool_t val = 0){fVetexingMassFitTest=val;};

  void   SetEventSelectionMethod(int i = 0){fUseEventSelection=i;};
  void   SetJetSelectionMethod(int i = 0){fUseJetSelection=i;};
  void   SetJetTaggerMCMethod(int i = 0){fUseMCTagger=i;};
  virtual AliRDHFJetsCuts * GetJetCutsHF();
protected:
  virtual Bool_t IsEventSelected();
  virtual Bool_t IsEventSelectedLegacy( AliAODEvent * aev);
  virtual Bool_t IsEventSelectedpA( AliAODEvent * aev);
  virtual Bool_t IsEventSelectedpp( AliAODEvent * aev);
  virtual Bool_t IsEventSelectedHF( AliAODEvent * aev);

  virtual Bool_t IsJetSelected(const AliEmcalJet * jet);
  virtual Bool_t IsJetSelectedLegacy(const AliEmcalJet * jet);
  virtual Bool_t IsJetSelectedHF(const AliEmcalJet * jet);
  virtual Bool_t AddTagJet(AliEmcalJet * jet);
  virtual Bool_t IsQuality(const AliAODTrack *track ,EQualityClass qclass );
  virtual Bool_t IsV0DaughterRadius(const AliAODTrack *track ,Double_t &Radius);
  virtual Double_t GetDeltaPtRandomCone();
  virtual Double_t GetPtCorrected(const AliEmcalJet * jet);
  virtual Bool_t FindVertexNProngSimple(const AliEmcalJet * jet, AliAODVertex * &vtx, Int_t &nProng);
  virtual Double_t GetVertexInvariantMass(AliAODVertex *vtx,Double_t massParticle);

  // Jet containers
  AliJetContainer            *fJetsCont;                   //!Jets
  AliJetContainer            *fJetsContMC;                 //!Jets MC
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliHFJetsTagging           *fTaggingHFClass;             //!Helper
  AliHFJetTaggingIP          *fTrackCountingTagger;        //!Track counting class
  TRandom3                   *fRandom;                     //! Random cone input
  AliAnalysisUtils           *fUtils;                      //! AliAnlalysisUtils helper object;
  AliRDHFJetsCuts            *fJetCutsHF;                  // HF jet cuts
private:
  TClonesArray * fMCparticles; //! Monte Carlo particle stack 
  Bool_t fIsMC;// Is Monte Carlo Event
  Bool_t fIsTrackQA;// Run track QA for analysis note plots
  Bool_t fIsTrackQAConstituent;// Run constituent track QA for analysis note plots
  Bool_t fUseCorrectedJetPt;// Subtract average rho
  Bool_t fDoRandomCones;//Get fluctuations
  Bool_t fVetexingMassFitTest; //
  Int_t  fUseEventSelection;// Event selection method
  Int_t  fUseJetSelection;// Jet selection method
  Int_t  fUseMCTagger;// MC tagger method selection

  //Histograms


  TH1 * fhist_Events;//! Event selection statistics
  TH2 * fhist_Jets;//!   Jet selection statistics
  TH2 * fhist_MonteCarloFlavour;//!
  TH2 * fhists_SPD_cluster_vs_tracklet_correlation;//! monitor_plot_for pileup rejection
  TH2 * fhists_SPD_cluster_vs_tracklet_correlation_PostSelection;//! monitor_plot_for pileup rejection
  TH2 * fhist_Tracks_Eta_Phi;//! Regular hybrid tracks
  TH2 * fhist_Tracks_Eta_Phi_Bit4;//! Regular ITS tracks
  TH2 * fhist_Tracks_Eta_Phi_Bit9;//! Regular complementary tracks
  TH1 * fhist_QualityClasses;//! Tracks in quality class X
  TH2 * fhist_QualityClasses_sIP[4][2];//![4][2]] prim. vs sec. track sip
  TH2 * fhist_QualityClasses_Eta_Phi[4][2];//![4][2]] prim. vs sec. track eta phi
  TH2 * fhist_TC_sIP_Pt[3][4][5];//![3][4][5] Track counting output
  TH2 * fhist_TC_Eta_Phi[3][4][5];//![3][4][5] Track counting eta phi
  TH2 * fhist_Jet_Eta_Phi;//! Jet eta-phi  distribution
  TH2 * fhist_Jet_Nconst_Pt;//! Jet pT vs N constituents
  TH1 * fhist_Jet_Pt;//! Jet pT  distribution
  TH1 * fhist_Jet_Background_Fluctuation;//! Background Fluctuations
  TH1 * fhist_Rho;//! Average background momentum density

  TH1 * fhist_parton_genjet_dR;//! Delta R mc jet matched parton
  TH2 * fhist_parton_genjet_pT;//! p_T mc jet matched parton
  TH1 * fhist_parton_genjet_Eta;//! \eta mc jet matched parton
  TH1 * fhist_parton_genjet_Phi;//! \phi mc jet matched parton  
  
  TH2 * fhist_momentum_response[4][2];//![4][2] momentum response matrix
  
  //Vertexing

  //TH2 * fHist_nProngsVsJetPt; //! number of Prongs vs Jet pT
  TH2 * fHist_2Prong_MassVsJetPt; //! Inv Vtx Mass vs. Jet pT 2 Prong;
  TH2 * fHist_3Prong_MassVsJetPt; //! Inv Vtx Mass vs. Jet pT 3 Prong;
  TH2 * fHist_4Prong_MassVsJetPt; //! Inv Vtx Mass vs. Jet pT 4 Prong;
  TH2 * fHist_5Prong_MassVsJetPt; //! Inv Vtx Mass vs. Jet pT 5 Prong;
  TH2 * fHist_MassVsJetPtHE[4]; //! Inv Vtx Mass vs. Jet pT all  Prong;
  TH2 * fHist_MassVsJetPtHP[4]; //! Inv Vtx Mass vs. Jet pT all  Prong;
  virtual void FillVertexingHists(const AliEmcalJet *jet, AliAODVertex * vtx, Int_t nProng,Int_t tag=0);
  virtual void AddHistTH1 (TH1 **hist,const char* histname, const char * title, const char *titlex, const char *titley, Int_t nBinsX, Double_t minX,Double_t maxX, Bool_t setSumw2,TList * container);
  virtual void AddHistTH2 (TH2 **hist,const char* histname,const char * title,const char *titlex, const char *titley,Int_t nBinsX, Double_t minX,Double_t maxX,Int_t nBinsY, Double_t minY,Double_t maxY,Bool_t setSumw2,TList * container);
  AliAnalysisTaskEmcalJetBJetTaggingIP(const AliAnalysisTaskEmcalJetBJetTaggingIP&);            // not implemented
  AliAnalysisTaskEmcalJetBJetTaggingIP &operator=(const AliAnalysisTaskEmcalJetBJetTaggingIP&); // not implemented


  ClassDef(AliAnalysisTaskEmcalJetBJetTaggingIP,102) // jet sample analysis task
};
#endif








