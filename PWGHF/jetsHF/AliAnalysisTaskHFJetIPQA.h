#ifndef ALIANALYSISTASKJETIPQA_H
#define ALIANALYSISTASKJETIPQA_H
#include "AliHFJetsTagging.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "TGraph.h"
class AliEmcalJet;
class AliRDHFJetsCuts;
class AliAODVertex;
class AliAODTrack;
class TList;
class TH1D;
class TH2D;
class AliHFJetsTagging;
class TParticle;
class TClonesArray;
class AliAODMCParticle;
class AliMCEvent;
class AliESDEvent;
class AliESDtrack;
class TGraph;
class AliAnalysisUtils;
class TRandom3;
class AliTriggerAnalysis;
class THnSparse;
class AliOADBContainer;
class AliEmcalList;
class AliVertexerTracks;
class TGraph;
#include "TMatrixD.h"
#include "TF1.h"
#include "AliESDtrackCuts.h"
#include <vector>
#include <utility>
#include <map>
class AliAnalysisTaskHFJetIPQA: public AliAnalysisTaskEmcalJet
{
public:
  //STATIC ENUM DEFINITIONS
  enum EPileup {kNoPileupSelection,kRejectPileupEvent,kRejectTracksFromPileupVertex};
  enum ERejBits {kNotSelTrigger,kNoVertex,kTooFewVtxContrib,kVertexChi2NDF,kZVtxOutFid,kPileupSPD,kOutsideCentrality,kVertexZContrib,kPhysicsSelection,kNoContributors,kDeltaVertexZ,kNoVertexTracks,kVertexZResolution,kMVPileup,kSPDClusterCut,kZVtxSPDOutFid};
  enum TTypeImpPar {kXY,kXYSig,kXYZ,kXYZSig,kXYZSigmaOnly,kZSig};
  enum EParticleType  {bPi0=111,bEta=221,bEtaPrime=331,bPhi=333,bRho=113,bOmega=223,bSigma0=3212,bK0s=310,bLambda=3122,bPi=211,bProton=2212,bKaon=321,bOmegaBaryon=3334,
    bAntiOmegaBaryon=-3334,bXiBaryon=3312,bAntiXiBaryon=-3312,bD0=411,bDPlus=421,bDStarPlus=413,bDSPlus=431,bK0l=130,bSigmaPlus = 3222,bRhoPlus=213,
    bBPlus = 521,bB0 = 511,bLambdaB =5122,bLambdaC=4122,bBStarPlus=523,bK0S892 = 313,bK0S892plus = 323};
  enum EParticleArrayIdx
  {bIdxPi0=0,bIdxEta=1,bIdxEtaPrime=2,bIdxPhi=3,bIdxRho=4,bIdxOmega=5,bIdxK0s=6,bIdxLambda=7,bIdxPi=8,bIdxProton=9,bIdxKaon=10,bIdxD0=11,bIdxDPlus=12,
    bIdxDStarPlus=13,bIdxDSPlus=14,bIdxLambdaC=15,bIdxBPlus = 16,bIdxB0 = 17,bIdxLambdaB = 18,bIdxBStarPlus=19};
  //UTILITY STRUCT DEFINITIONS
  struct SJetIpPati {
    SJetIpPati(Double_t v1, Double_t v2, Bool_t b, Bool_t c): first(v1),second(v2),is_electron(b),is_fromB(c){}
    Double_t first; // to be compatible with std::pair
    Double_t second;// to be compatible with std::pair
    Bool_t   is_electron; // added for electron contribution check
    Bool_t   is_fromB; // added for electron contribution check
  };
  //FUNCTION DEFINITIONS
  AliAnalysisTaskHFJetIPQA();
  AliAnalysisTaskHFJetIPQA(const char *name);
 // AliAnalysisTaskHFJetIPQA(const AliHFJetsTagging&){}
  virtual ~AliAnalysisTaskHFJetIPQA(){;}
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Run();
  virtual Bool_t IsEventSelected();
  virtual Bool_t Notify();
  virtual Bool_t IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits);
  void DoJetProbabilityAnalysis(Bool_t val=true){fDoJetProbabilityAnalysis=val;}
  void SetESDCuts (AliESDtrackCuts  *cuts =NULL){fESDTrackCut =  new AliESDtrackCuts(*cuts);}
  void SetRunESD (Bool_t val = kTRUE){fESD = val;}
  virtual AliRDHFJetsCuts* GetJetCutsHF(){return fJetCutsHF;}
  void SetUseMonteCarloWeighingLinus(TH1F *Pi0 ,TH1F *Eta,TH1F *EtaP,TH1F *Rho,TH1F *Phi,TH1F *Omega,TH1F *K0s,TH1F *Lambda,TH1F *ChargedPi,
                                     TH1F *ChargedKaon,TH1F *Proton,TH1F *D0,TH1F *DPlus,TH1F *DStarPlus,
                                     TH1F *DSPlus,TH1F *LambdaC,TH1F *BPlus,TH1F *B0,TH1F *LambdaB,TH1F *BStarPlus);
  Bool_t SetResFunction( TGraph * f = 0x0, Int_t j=0);
  void SetAODBContainer(AliOADBContainer* cont);
  void EnableCorrectionSamplingMode(Bool_t val=true){fCorrrectionSamplingMode=val;}//
  void setN_ITSClusters_Input_global( Int_t value);
  void RotateTracksAroundYAxis(double angle, double shiftx);
  void RotateMatrixY(TMatrixD &m, double ang_rad);
  void PartiallyRotateMatrixY(TMatrixD &m, double ang_rad);
  void RotateVectorY(TVector3 &v, double ang_rad);
  void DisableCompositionCorrection(Bool_t val = kTRUE){fDisableWeightingMC = val;}
private:
  void DoJetLoop(); //jet matching function 2/4
  void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Int_t matching=0);
  void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
  void SmearTrackHybrid(AliVTrack * track);
  void FillHist(const char * name,Double_t x ,Double_t w);
  void FillHist(const char * name,Double_t x, Double_t y,Double_t w);
  void IncHist(const char * name,Int_t bin);
  void SubtractMean (Double_t val[2],AliVTrack *track);
  Bool_t CalculateTrackImpactParameter(AliAODTrack * track,Double_t *impar, Double_t * cov); // Removes track from Vertex calculation first
  Bool_t CalculateTrackImpactParameter(AliESDtrack * track,Double_t *impar, Double_t * cov,Bool_t useTRUEvtx=false); // Removes track from Vertex calculation first
  Bool_t CalculateTrackImpactParameter(AliVTrack * track,Double_t *impar, Double_t * cov);
  Bool_t CalculateTrackImpactParameterTruth(AliAODTrack * track,Double_t *impar, Double_t * cov); // calculates DCA on MC particle/event information
  Bool_t CalculateTrackImpactParameterTruth(AliESDtrack * track,Double_t *impar, Double_t * cov); // calculates DCA on MC particle/event information
  Bool_t CalculateJetSignedTrackImpactParameter(AliAODTrack * track,AliEmcalJet * jet ,Double_t *impar, Double_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength);
  Bool_t CalculateJetSignedTrackImpactParameter(AliESDtrack * track,AliEmcalJet * jet ,Double_t *impar, Double_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength);
  Bool_t CalculateJetSignedTrackImpactParameter(AliVTrack * track,AliEmcalJet * jet ,Double_t *impar, Double_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength);
  Bool_t IsV0PhotonFromBeamPipeDaughter(const AliAODTrack* track);
  Bool_t IsV0PhotonFromBeamPipeDaughter(const AliESDtrack* track);
  Bool_t IsTrackAccepted(AliVTrack* track,Int_t n=6);
  Bool_t MatchJetsGeometricDefault(); //jet matching function 1/4
  Bool_t ParticleIsPossibleSource(Int_t pdg);
  Bool_t IsSelectionParticle( AliAODMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t IsSelectionParticle( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t IsSelectionParticleALICE( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t IsSelectionParticleStrange( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t IsSelectionParticleMeson( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t IsSelectionParticleOmegaXiSigmaP( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx );
  Bool_t IsSecondaryFromWeakDecay( AliAODMCParticle * particle ) ;
  Bool_t IsSecondaryFromWeakDecay( AliMCParticle * particle ) ;
  Bool_t IsTruePrimary	(AliMCParticle * mcpart);
  Bool_t GetBMesonWeight( AliAODMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t GetBMesonWeight( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  );
  Bool_t GetESDITSMODULEINFO(AliESDtrack * track);
  Bool_t IsPromptDMeson(AliAODMCParticle * part );
  Bool_t IsPromptDMeson(AliMCParticle * part );
  Bool_t IsPromptBMeson(AliAODMCParticle * part );
  Bool_t IsPromptBMeson(AliMCParticle * part );
  static Bool_t mysort(const SJetIpPati& i, const SJetIpPati& j);
  Int_t IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg);
  Int_t GetRunNr(AliVEvent * event){return event->GetRunNumber();}
  Double_t CalculatePSTrack(Double_t sign, Double_t significance ,Double_t trackPt,Int_t trclass);
  Double_t CalculateJetProb(AliEmcalJet * jet);//!
  Double_t GetValImpactParameter(TTypeImpPar type,Double_t *impar, Double_t * cov);
  Double_t GetMonteCarloCorrectionFactor(AliVTrack* track,Int_t &pCorr_indx);
  Double_t GetWeightFactor( AliAODMCParticle * mcpart,Int_t &pCorr_indx);
  Double_t GetWeightFactor( AliMCParticle * mcpart,Int_t &pCorr_indx);
  Double_t GetArmenteros(AliESDv0 * v0 , Int_t pidneg,Int_t pidpos ,Double_t &alpha);
  Double_t GetPsiPair(AliESDv0 * v0);
  Double_t GetPtCorrected(const AliEmcalJet* jet);
  Double_t GetPtCorrectedMC(const AliEmcalJet *jet);
  //Functions to allow jet probability/TC System 8 efficiency estimation
  Bool_t IsJetTaggedTC(int n =0 ,double thres = 0.1);
  Bool_t IsJetTaggedJetProb(double thresProb = 0.90);


  void GetUDGResolutionFunctionHists(AliVTrack * track,AliEmcalJet * jet);
  AliAODMCParticle* GetMCTrack( const AliAODTrack* _track);
  TH1 *  AddHistogramm(const char * name,const char * title,Int_t x,Double_t xlow,Double_t xhigh, Int_t y=0,Double_t ylow=0,Double_t yhigh=0);
  TH1D * GetHist1D(const char * name){return (TH1D*)fOutput2->FindObject(name);}
  TH2D * GetHist2D(const char * name){return (TH2D*)fOutput2->FindObject(name);}

  TGraph * fGraphMean;//!
  TGraph * fGraphSigmaData;//!
  TGraph * fGraphSigmaMC;//!
  TGraph * fGraphXi; //!
  TGraph * fGraphOmega;
  TGraph * fK0Star; //!
  TGraph * fPhi; //!
  TGraph * fGeant3FlukaProton;//!
  TGraph * fGeant3FlukaAntiProton;//!
  TGraph * fGeant3FlukaLambda;//!
  TGraph * fGeant3FlukaAntiLambda;//!
  TGraph * fGeant3FlukaKMinus;//!
  TList   * fOutput2;//!
  TClonesArray     *fMCArray;//!
  AliRDHFJetsCuts  *fJetCutsHF;//
  AliOADBContainer *fAODBcont;//
  AliMCEvent       *fMCEvent;//!;
  AliESDtrackCuts  *fESDTrackCut;//
  AliAnalysisUtils *fUtils;//!
  AliVertexerTracks *fVertexer;//!
  Bool_t fDoJetProbabilityAnalysis;
  Bool_t fDisableWeightingMC;
  Bool_t fESD ;
  Bool_t fMcEvtSampled;//
  Bool_t fCorrrectionSamplingMode;//
  Int_t  fItsClustersInputGlobal;
  Double_t fBackgroundFactorLinus[21][498]; //[21][498]FineBinned correction factors up 0.1-25 GeV/c first value below last above 0.05 binwidth
  std::vector <Double_t > fEtaSEvt;//!
  std::vector <Double_t > fPhiSEvt;//!
  std::vector <Double_t > fEtaBEvt;//!
  std::vector <Double_t > fPhiBEvt;//!
  std::vector <Double_t > fEtaCEvt;//!
  std::vector <Double_t > fPhiCEvt;//!
  std::vector <Double_t > fEtaUdsgEvt;//!
  std::vector <Double_t > fPhiUdsgEvt;//!

  TGraph fResolutionFunction [5];//[5] ie 5 * n Pt bins
  TH2D * fh2dAcceptedTracksEtaPhiPerLayer[6];//![6]
  ClassDef(AliAnalysisTaskHFJetIPQA, 9	)
};
#endif

