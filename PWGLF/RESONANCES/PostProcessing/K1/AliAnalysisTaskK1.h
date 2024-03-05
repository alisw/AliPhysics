#ifndef AliAnalysisTaskK1_H
#define AliAnalysisTaskK1_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>
#include <tuple>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class AliAODMCHeader;
class AliPIDResponse;
class AliESDtrackCuts;

struct SK1Daughter
{
  enum
  {
    kPrimary = BIT(0),
    kSecondary = BIT(1)
  };
  int eventID;
  char centrality;
  float zVertex;
  int id;
  int motherID;
  int daughter1;
  int daughter2;
  float pt;
  float eta;
  int pdg;
  int flag;
};

struct RK1Daughter
{
  enum
  {
    kPositive = BIT(0),
    kHasTOF = BIT(1)
  };
  int eventID;
  char centrality;
  float zVertex;
  int id;
  int pdg;
  int motherID;
  float px;
  float py;
  float pz;
  float eta;
  Double32_t dcaxy;
  Double32_t dcaz;
  float tofNsigmaPi;
  float tpcNsigmaPi;
  float tofNsigmaKa;
  float tpcNsigmaKa;
  unsigned char flag;
};

class AliAnalysisTaskK1 : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskK1();
  AliAnalysisTaskK1(const char *name, Bool_t MCcase);
  virtual ~AliAnalysisTaskK1();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetFillQAPlot(Bool_t input) { fFillQAPlot = input; }

  void SetMixing(Bool_t setmixing) { fSetMixing = setmixing; }
  void SetnMix(Int_t nMix) { fnMix = nMix; }
  void SetIsPrimaryMC(Bool_t isprimarymc) { fIsPrimaryMC = isprimarymc; }
  void SetINEL(Bool_t input) { fIsINEL = input; }
  void SetHighMult(Bool_t input) { fIsHM = input; }
  void SetFillTree(Bool_t fillTree) { fFillTree = fillTree; }
  void SetSkipFillHistos(Bool_t skipHisto = kTRUE) { fSkipFillingHistogram = skipHisto; }

  // Setter for cut variables
  void SetFilterbitTracks(Double_t lParameter) { fFilterBit = lParameter; }

  void SetMaxTPCnSigPrimaryPion(Double_t lParameter) { fTPCNsigPrimaryPionCut = lParameter; }
  void SetMaxTOFnSigPrimaryPion(Double_t lParameter) { fTOFNsigPrimaryPionCut = lParameter; }
  void SetMaxEtaPrimaryPion(Double_t lParameter) { fPrimaryPionEtaCut = lParameter; }
  void SetMaxVertexZPrimaryPion(Double_t lParameter) { fPrimaryPionZVertexCut = lParameter; }
  void SetMaxVertexXYsigPrimaryPion(Double_t lParameter) { fPrimaryPionXYVertexSigmaCut = lParameter; }

  void SetMaxTPCnSigSecondaryPion(Double_t lParameter) { fTPCNsigSecondaryPionCut = lParameter; }
  void SetMaxTOFnSigTOFSecondaryPion(Double_t lParameter) { fTOFNsigSecondaryPionCut = lParameter; }
  void SetMaxEtaSecondaryPion(Double_t lParameter) { fSecondaryPionEtaCut = lParameter; }
  void SetMaxVertexZSecondaryPion(Double_t lParameter) { fSecondaryPionZVertexCut = lParameter; }
  void SetMaxVertexXYsigSecondaryPion(Double_t lParameter) { fSecondaryPionXYVertexSigmaCut = lParameter; }

  void SetMaxTPCnSigKaon(Double_t lParameter) { fTPCNsigKaonCut = lParameter; }
  void SetMaxTOFnSigKaon(Double_t lParameter) { fTOFNsigKaonCut = lParameter; }
  void SetMaxEtaKaon(Double_t lParameter) { fKaonEtaCut = lParameter; }
  void SetMaxVertexZKaon(Double_t lParameter) { fKaonZVertexCut = lParameter; }
  void SetMaxVertexXYsigKaon(Double_t lParameter) { fKaonXYVertexSigmaCut = lParameter; }

  void SetModeK892orRho(Bool_t lParameter) { fSetModeK892orRho = lParameter; }
  void SetSecondaryMassWindow(Double_t lParameter) { fSecondaryMassWindowCut = lParameter; }
  void SetSecondaryRapidityCut(Double_t lParameter) { fSecondaryRapCut = lParameter; }

  void SetK1RapidityCutHigh(Double_t lParameter) { fK1YCutHigh = lParameter; }
  void SetK1RapidityCutLow(Double_t lParameter) { fK1YCutLow = lParameter; }
  void SetK1OAMin(Double_t lParameter) { fMinK1OA = lParameter; }
  void SetK1OAMax(Double_t lParameter) { fMaxK1OA = lParameter; }
  void SetK1PairAssymMin(Double_t lParameter) { fMinPairAsym = lParameter; }
  void SetK1PairAssymMax(Double_t lParameter) { fMaxPairAsym = lParameter; }
  void SetK1AnotherSecnarioMassCutMin(Double_t lParameter) { fMinK1MassCutAnotherScenario = lParameter; }
  void SetK1AnotherSecnarioMassCutMax(Double_t lParameter) { fMaxK1MassCutAnotherScenario = lParameter; }
  void SetK1PiKaMassCutMin(Double_t lParameter) { fMinK1PiKa = lParameter; }
  void SetK1PiKaMassCutMax(Double_t lParameter) { fMaxK1PiKa = lParameter; }
  void SetSecondarypTCutMin(Double_t lParameter) { fMinSecondarypTCut = lParameter; }
  void SetSecondarypTCutMax(Double_t lParameter) { fMaxSecondarypTCut = lParameter; }

  Bool_t FillTrackPools();
  void FillHistograms();
  void FillMCinput(AliMCEvent *fMCEvent, int Fillbin = 0);
  void FillTrackToEventPool();
  Bool_t IsTrueK1(UInt_t v0, UInt_t pion, UInt_t BkgCheck = 0);
  Double_t GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type);
  Double_t GetTOFnSigma(AliVTrack *track, AliPID::EParticleType type);
  void PostAllData();
  void GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3]);
  std::tuple<Float_t, Float_t, std::array<Float_t, 3>> GetImpactParametersWrapper(AliVTrack *track);
  bool IsTrackAccepted(AliVTrack *track);
  Int_t TrackSelection(AliVTrack *track, Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t nTPCNSigKaon, Float_t nTOFNSigKaon, Float_t lpT, Float_t lDCAz, Float_t lDCAr, Float_t lEta);
  bool IsPassPrimaryPionSelection(Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr);
  bool IsPassSecondaryPionSelection(Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr);
  bool IsPassKaonSelection(Float_t nTPCNSigKaon, Float_t nTOFNSigKaon, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr);
  void SetCutOpen();

  // helper
  THnSparseD *CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt = "");
  Long64_t FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w = 1.);
  TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
  TAxis AxisVar(TString name, std::vector<Double_t> bin);
  TAxis AxisStr(TString name, std::vector<TString> bin);

  AliEventCuts fEventCut; // Event cuts

private:
  typedef std::vector<AliVTrack *> tracklist;
  typedef std::deque<tracklist> eventpool;
  typedef std::vector<std::vector<eventpool>> mixingpool;

  AliESDtrackCuts *fTrackCuts;  //!
  AliPIDResponse *fPIDResponse; //!

  AliVEvent *fEvt;              //!
  AliMCEvent *fMCEvent;         //!
  AliAODMCHeader *fAODMCHeader; //!<  MC info AOD
  TList *fList;                 //!
  TClonesArray *fMCArray;       //!
  AliAODVertex *fVertex;        //!
  TTree *fNanoTree;             ///<  Output nanoAOD ttree
  TTree *fNanoMCTree;           ///<  Output nanoAOD MC ttree
  RK1Daughter fRecK1daughter;   ///<  Reconstructed nucleus
  SK1Daughter fSimK1part;       ///<  Simulated nucleus

  //// Histograms
  THnSparseD *fHn5DK1Data;  //!
  THnSparseD *fHn5DK1MC;    //!
  THnSparseD *fHn2DEvtNorm; //!
  //// QA plots
  TH1F *hMultiplicity; //!
  // All tracks
  TH1D *hEtaTrack_before;    //!
  TH1D *hDCAPVTrack_before;  //!
  TH1D *hDCArPVTrack_before; //!
  TH1D *hPtTrack_before;     //!
  // Primary pions
  TH1D *hEtaTrack_ppion;            //!
  TH1D *hDCAPVTrack_ppion;          //!
  TH1D *hDCArPVTrack_ppion;         //!
  TH1D *hPtTrack_ppion;             //!
  TH2D *hTPCPIDTrack_ppion;         //!
  TH2D *hTPCPIDTrackNsigVspT_ppion; //!
  // Secondary pions
  TH1D *hEtaTrack_spion;            //!
  TH1D *hDCAPVTrack_spion;          //!
  TH1D *hDCArPVTrack_spion;         //!
  TH1D *hPtTrack_spion;             //!
  TH2D *hTPCPIDTrack_spion;         //!
  TH2D *hTPCPIDTrackNsigVspT_spion; //!
  // Kaons
  TH1D *hEtaTrack_kaon;            //!
  TH1D *hDCAPVTrack_kaon;          //!
  TH1D *hDCArPVTrack_kaon;         //!
  TH1D *hPtTrack_kaon;             //!
  TH2D *hTPCPIDTrack_kaon;         //!
  TH2D *hTPCPIDTrackNsigVspT_kaon; //!
  // K1
  TH1D *hK1OA;                          //!
  TH1D *hK1PairAsymm;                   //!
  TH2D *hInvMass_k892_rho;              //!
  TH2D *hInvMass_Secondary_PiKa;        //!
  TH1D *hK1SecondarypT;                 //!
  TH1D *hK1OA_cut;                      //!
  TH1D *hK1PairAsymm_cut;               //!
  TH2D *hInvMass_k892_rho_cut;          //!
  TH2D *hInvMass_Secondary_PiKa_cut;    //!
  TH1D *hK1SecondarypT_cut;             //!
  TH1D *hK1OA_MCTrue;                   //!
  TH1D *hK1PairAsymm_MCTrue;            //!
  TH2D *hInvMass_k892_rho_MCTrue;       //!
  TH2D *hInvMass_Secondary_PiKa_MCTrue; //!
  TH1D *hK1SecondarypT_MCTrue;          //!

  Bool_t fIsAOD;                //!
  Bool_t fIsNano;               //!
  Bool_t fSetMixing;            //
  Bool_t fFillQAPlot;           //
  Bool_t fIsMC;                 //
  Bool_t fIsPrimaryMC;          //
  Bool_t fFillTree;             //
  Bool_t fIsINEL;               //
  Bool_t fIsHM;                 //
  Bool_t fSkipFillingHistogram; //
  Bool_t fSetModeK892orRho;     // True: K892, False: Rho mode.

  mixingpool fEMpool; //!
  TAxis fBinCent;     //!
  TAxis fBinZ;        //!
  Double_t fPosPV[3]; //!
  Double_t fMagField; //!

  Double_t fCent;               //!
  unsigned long fCustomEventID; //!
  Int_t fnMix;                  //!
  Int_t fCentBin;               //!
  Int_t fZbin;                  //!

  // Pion cuts
  UInt_t fFilterBit;                     //
  Double_t fTPCNsigPrimaryPionCut;       //
  Double_t fTOFNsigPrimaryPionCut;       //
  Double_t fPrimaryPionEtaCut;           //
  Double_t fPrimaryPionZVertexCut;       //
  Double_t fPrimaryPionXYVertexSigmaCut; //

  Double_t fTPCNsigSecondaryPionCut;       //
  Double_t fTOFNsigSecondaryPionCut;       //
  Double_t fSecondaryPionEtaCut;           //
  Double_t fSecondaryPionZVertexCut;       //
  Double_t fSecondaryPionXYVertexSigmaCut; //

  Double_t fTPCNsigKaonCut;       //
  Double_t fTOFNsigKaonCut;       //
  Double_t fKaonEtaCut;           //
  Double_t fKaonZVertexCut;       //
  Double_t fKaonXYVertexSigmaCut; //

  Double_t fSecondaryMassWindowCut; //
  Double_t fSecondaryRapCut;        //

  Double_t fK1YCutHigh;                  //
  Double_t fK1YCutLow;                   //
  Double_t fMinK1OA;                     //
  Double_t fMaxK1OA;                     //
  Double_t fMinPairAsym;                 //
  Double_t fMaxPairAsym;                 //
  Double_t fMinK1MassCutAnotherScenario; //
  Double_t fMaxK1MassCutAnotherScenario; //
  Double_t fMinK1PiKa;                   //
  Double_t fMaxK1PiKa;                   //
  Double_t fMinSecondarypTCut;           //
  Double_t fMaxSecondarypTCut;           //

  // Good track arrays
  std::vector<UInt_t> fGoodPrimaryPionArray;
  std::vector<UInt_t> fGoodSecondaryPionArray;
  std::vector<UInt_t> fGoodKaonArray;
  std::vector<UInt_t> fGoodTracksArray;

  unsigned long long int fTracks; //!

  ClassDef(AliAnalysisTaskK1, 3);
  // 1: first implementation
  // 2: update Tree method and simplify the code
  // 3: Apply the robust daugther selection for both scenarios and add the secondary pT cut
};

#endif
