#ifndef AliAnalysisTaskK1_H
#define AliAnalysisTaskK1_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

struct SK1Resonance
{
  int pdg_pion1;
  int pdg_pion2;
  int pdg_pion3;
  int flag;
};
struct RK1Resonance
{
  // first pion
  float pt_pion1;
  float eta_pion1;
  bool hasTOF_pion1;
  Double32_t dcaxy_pion1;
  Double32_t dcaz_pion1;
  float tofNsigma_pion1;
  float tpcNsigma_pion1;
  // second pion
  float pt_pion2;
  float eta_pion2;
  bool hasTOF_pion2;
  Double32_t dcaxy_pion2;
  Double32_t dcaz_pion2;
  float tofNsigma_pion2;
  float tpcNsigma_pion2;
  // kaon
  float pt_kaon;
  float eta_kaon;
  bool hasTOF_kaon;
  Double32_t dcaxy_kaon;
  Double32_t dcaz_kaon;
  float tofNsigma_kaon;
  float tpcNsigma_kaon;
  // event info
  float centrality;
  float vertexZ;
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
  void SetSkipFillHistos(Bool_t skipHisto) { fSkipFillingHistogram = skipHisto; }

  // Setter for cut variables
  void SetFilterbitTracks(Double_t lParameter) { fFilterBit = lParameter; }

  void SetMaxNsigPrimaryPion(Double_t lParameter) { fTPCNsigPrimaryPionCut = lParameter; }
  void SetMaxEtaPrimaryPion(Double_t lParameter) { fPrimaryPionEtaCut = lParameter; }
  void SetMaxVertexZPrimaryPion(Double_t lParameter) { fPrimaryPionZVertexCut = lParameter; }
  void SetMaxVertexXYsigPrimaryPion(Double_t lParameter) { fPrimaryPionXYVertexSigmaCut = lParameter; }

  void SetMaxNsigSecondaryPion(Double_t lParameter) { fTPCNsigSecondaryPionCut = lParameter; }
  void SetMaxEtaSecondaryPion(Double_t lParameter) { fSecondaryPionEtaCut = lParameter; }
  void SetMaxVertexZSecondaryPion(Double_t lParameter) { fSecondaryPionZVertexCut = lParameter; }
  void SetMaxVertexXYsigSecondaryPion(Double_t lParameter) { fSecondaryPionXYVertexSigmaCut = lParameter; }

  void SetMaxNsigKaon(Double_t lParameter) { fTPCNsigKaonCut = lParameter; }
  void SetMaxEtaKaon(Double_t lParameter) { fKaonEtaCut = lParameter; }
  void SetMaxVertexZKaon(Double_t lParameter) { fKaonZVertexCut = lParameter; }
  void SetMaxVertexXYsigKaon(Double_t lParameter) { fKaonXYVertexSigmaCut = lParameter; }

  void SetMaxMassWindowK892(Double_t lParameter) { fK892MassWindowCut = lParameter; }
  void SetMaxRapidityCutK892(Double_t lParameter) { fK892RapCut = lParameter; }

  void SetK1RapidityCutHigh(Double_t lParameter) { fK1YCutHigh = lParameter; }
  void SetK1RapidityCutLow(Double_t lParameter) { fK1YCutLow = lParameter; }
  void SetK1OAMin(Double_t lParameter) { fMinK1OA = lParameter; }
  void SetK1OAMax(Double_t lParameter) { fMaxK1OA = lParameter; }
  void SetK1PairAssymMin(Double_t lParameter) { fMinPairAsym = lParameter; }
  void SetK1PairAssymMax(Double_t lParameter) { fMaxPairAsym = lParameter; }
  void SetK1PiPiMassCutMin(Double_t lParameter) { fMinK1PiPi = lParameter; }
  void SetK1PiPiMassCutMax(Double_t lParameter) { fMaxK1PiPi = lParameter; }
  void SetK1PiKaMassCutMin(Double_t lParameter) { fMinK1PiKa = lParameter; }
  void SetK1PiKaMassCutMax(Double_t lParameter) { fMaxK1PiKa = lParameter; }

  Bool_t FillTrackPools();
  void FillHistograms();
  void FillTree();
  void FillMCinput(AliMCEvent *fMCEvent, int Fillbin = 0);
  void FillTrackToEventPool();
  Bool_t IsTrueK1(UInt_t v0, UInt_t pion, UInt_t BkgCheck = 0);
  double GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type);
  void GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3]);
  Int_t trackSelection(AliVTrack* track, Float_t& nTPCNSigPion, Float_t& nTPCNSigKaon, Float_t lpT, Float_t lDCAz, Float_t lDCAr, Float_t lEta);
  void SetCutOpen();

  // helper
  THnSparse *CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt = "");
  THnSparse *CreateTHnSparse(TString name, TString title, TString templ, Option_t *opt = "");
  Long64_t FillTHnSparse(TString name, std::vector<Double_t> x, Double_t w = 1.);
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

  AliVEvent *fEvt;           //!
  AliMCEvent *fMCEvent;      //!
  THistManager *fHistos;     //!
  AliAODVertex *fVertex;     //!
  TTree *fRTree;             ///<  Output reconstructed ttree
  TTree *fSTree;             ///<  Output simulated ttree
  RK1Resonance fRecK1;         ///<  Reconstructed resonance
  SK1Resonance fSimK1;         ///<  Simulated resonance
  TClonesArray *fMCArray;    //!

  Bool_t fIsAOD;       //!
  Bool_t fIsNano;      //!
  Bool_t fSetMixing;   //
  Bool_t fFillQAPlot;  //
  Bool_t fIsMC;        //
  Bool_t fIsPrimaryMC; //
  Bool_t fFillTree;    //
  Bool_t fIsINEL;      //
  Bool_t fIsHM;        //
  Bool_t fSkipFillingHistogram; //

  mixingpool fEMpool; //!
  TAxis fBinCent;     //!
  TAxis fBinZ;        //!
  Double_t fPosPV[3]; //!
  Double_t fMagField; //!

  Double_t fCent; //!
  Int_t fnMix;    //!
  Int_t fCentBin; //!
  Int_t fZbin;    //!

  // Pion cuts
  UInt_t fFilterBit;                     //
  Double_t fTPCNsigPrimaryPionCut;       //
  Double_t fPrimaryPionEtaCut;           //
  Double_t fPrimaryPionZVertexCut;       //
  Double_t fPrimaryPionXYVertexSigmaCut; //

  Double_t fTPCNsigSecondaryPionCut;       //
  Double_t fSecondaryPionEtaCut;           //
  Double_t fSecondaryPionZVertexCut;       //
  Double_t fSecondaryPionXYVertexSigmaCut; //

  Double_t fTPCNsigKaonCut;       //
  Double_t fKaonEtaCut;           //
  Double_t fKaonZVertexCut;       //
  Double_t fKaonXYVertexSigmaCut; //

  Double_t fK892MassWindowCut; //
  Double_t fK892RapCut;        //

  Double_t fK1YCutHigh;        //
  Double_t fK1YCutLow;         //
  Double_t fMinK1OA;           //
  Double_t fMaxK1OA;           //
  Double_t fMinPairAsym;       //
  Double_t fMaxPairAsym;       //
  Double_t fMinK1PiPi;       //
  Double_t fMaxK1PiPi;       //
  Double_t fMinK1PiKa;       //
  Double_t fMaxK1PiKa;       //

  // Good track arrays
  std::vector<UInt_t> fGoodPrimaryPionArray;
  std::vector<UInt_t> fGoodSecondaryPionArray;
  std::vector<UInt_t> fGoodKaonArray;

  ClassDef(AliAnalysisTaskK1, 1);
  // 1: first implementation
};

#endif
