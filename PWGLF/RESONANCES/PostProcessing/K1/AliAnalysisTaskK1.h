#ifndef AliAnalysisTaskK1_H
#define AliAnalysisTaskK1_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliResoNanoEvent.h"
#include "AliResoNanoTrack.h"
class AliAODMCHeader;
class AliPIDResponse;
class AliESDtrackCuts;

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
  void FillMCinput(AliMCEvent *fMCEvent, int Fillbin = 0);
  void FillTrackToEventPool();
  Bool_t IsTrueK1(UInt_t v0, UInt_t pion, UInt_t BkgCheck = 0);
  Double_t GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type);
  Double_t GetTOFnSigma(AliVTrack *track, AliPID::EParticleType type);
  void GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3]);
  Int_t trackSelection(AliVTrack *track, Float_t &nTPCNSigPion, Float_t &nTPCNSigKaon, Float_t lpT, Float_t lDCAz, Float_t lDCAr, Float_t lEta);
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
  TClonesArray *fNanoEvents;    ///<  events for nanoAOD
  TClonesArray *fNanoTracks;    ///<  tracks for nanoAOD

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
  TH1D *hK1OA;                    //!
  TH1D *hK1PairAsymm;             //!
  TH2D *hInvMass_piK_pipi;        //!
  TH2D *hInvMass_piK_pika;        //!
  TH1D *hK1OA_cut;                //!
  TH1D *hK1PairAsymm_cut;         //!
  TH2D *hInvMass_piK_pipi_cut;    //!
  TH2D *hInvMass_piK_pika_cut;    //!
  TH1D *hK1OA_MCTrue;             //!
  TH1D *hK1PairAsymm_MCTrue;      //!
  TH2D *hInvMass_piK_pipi_MCTrue; //!
  TH2D *hInvMass_piK_pika_MCTrue; //!

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

  Double_t fK892MassWindowCut; //
  Double_t fK892RapCut;        //

  Double_t fK1YCutHigh;  //
  Double_t fK1YCutLow;   //
  Double_t fMinK1OA;     //
  Double_t fMaxK1OA;     //
  Double_t fMinPairAsym; //
  Double_t fMaxPairAsym; //
  Double_t fMinK1PiPi;   //
  Double_t fMaxK1PiPi;   //
  Double_t fMinK1PiKa;   //
  Double_t fMaxK1PiKa;   //

  // Good track arrays
  std::vector<UInt_t> fGoodPrimaryPionArray;
  std::vector<UInt_t> fGoodSecondaryPionArray;
  std::vector<UInt_t> fGoodKaonArray;
  std::vector<UInt_t> fGoodTracksArray;

  ClassDef(AliAnalysisTaskK1, 1);
  // 1: first implementation
};

#endif
