#ifndef AliAnalysisTaskNFactorialMoments_H
#define AliAnalysisTaskNFactorialMoments_H

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "TNtuple.h"
#include "TH3F.h"

const Int_t mMBins = 40, mDim = 2, mQs = 6, mPtmax = 5;

class TH1F;
class TH1D;
class TH2D;
class TH2F;
class TList;
class TNtuple;
class TString;
class TList;
class AliAODEvent;
class AliMCEvent;
class AliPIDResponse;
class AliPIDCombined;

class AliAnalysisTaskNFactorialMoments : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskNFactorialMoments();
  AliAnalysisTaskNFactorialMoments(const char* name);
  virtual ~AliAnalysisTaskNFactorialMoments();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  // Analysis Setters
  void SetTwoTrackCuts(Double_t deta, Double_t dphi, Bool_t twotrack, Bool_t QA)
  {
    this->fdeta = deta;
    this->fdphi = dphi;
    this->flag2Track = twotrack;
    this->flag2TrackQA = QA;
  }
  void SetSharingFraction(Double_t fshfr, Bool_t sharity)
  {
    this->fSharedFraction = fshfr;
    this->flagSharity = sharity;
  }
  void SetVtxCut(Double_t fVx, Double_t fVy, Double_t fVz)
  {
    this->fVxMax = fVx;
    this->fVyMax = fVy;
    this->fVzMax = fVz;
  }
  void SetIsMC(Bool_t isMC, TString gen)
  {
    this->flagMC = isMC;
    this->fGenName = gen;
  }
  void SetCentLim(Int_t f_minCent, Int_t f_maxCent)
  {
    this->minCent = f_minCent;
    this->maxCent = f_maxCent;
  }
  void SetEta(Double_t f_minEta, Double_t f_maxEta)
  {
    this->minEta = f_minEta;
    this->maxEta = f_maxEta;
  }
  void Setfbit(Int_t fb) { this->filterBit = fb; }
  void SetYear(TString year) { this->fYear = year; }
  void SetPSbinning(TArrayI MB, TArrayI NB, Int_t Mm)
  {
    this->mMBin2 = MB;
    this->mNBin2 = NB;
    this->Mmax = Mm;
  }
  void SetPtArray(TArrayD f_ptArray, Int_t nbins)
  {
    this->ptarray = f_ptArray;
    this->mPtBins = nbins;
  }
  void SetRejectElectrons(Bool_t fRejectElectron)
  {
    this->flagRejEls = fRejectElectron;
  }
  void SetRejectPileup(Bool_t RejectPileup) { this->flagPileup = RejectPileup; }
  void SetBfield(Int_t bf) { this->fBfield = bf; }
  void SetDCAXYRangeMax(Double_t dcaxy) { this->fDCAxyMax = dcaxy; }
  void SetDCAZRangeMax(Double_t dcaz) { this->fDCAzMax = dcaz; }
  void SetITSClusterCut(Double_t itscut) { this->fITSCls = itscut; }
  void SetTPCClusterCut(Double_t tpccut) { this->fTPCCls = tpccut; }
  void SetnTPCrossedrows(Double_t tpcrowcut) { this->fTPCRows = tpcrowcut; }
  void SetSelfAffine(Bool_t selfaffine) { this->flagSelfAff = selfaffine; }
  void SetSharedCls(Double_t fSharedCls, Double_t fSharedrows, Double_t fFindable, Bool_t fSharedMean)
  {
    this->fSharedClsMax = fSharedCls;
    this->fSharedRowsMax = fSharedrows;
    this->fFindableClsMin = fFindable;
    this->flagShClPro = fSharedMean;
  }

 protected:
  AliEventCuts fEventCuts;
  void FillTrackInfo();
  void FillMCTrackInfo();
  float GetDPhiStar(float phi1, float pt1, float charge1, float phi2, float pt2,
                    float charge2, float radius, float bSign);
  Bool_t GetParticleID(AliAODTrack* track, Bool_t fQA);
  float CalculateDPhiStar(float phi1, float eta1, float pt1, Int_t charge1,
                          float phi2, float eta2, float pt2, Int_t charge2,
                          float bSign);
  float SharedClusterFraction(TBits&, TBits&, TBits&, TBits&);
  void GetPtBin(Double_t);
  void CalculateNFMs(TH2D* h1[mPtmax][mMBins], Bool_t mcGen);
  void DataPosting();
  void ResetHistograms();

 private:
  // Event related (Data and MC)
  AliAODEvent* fAOD;    // AOD object
  AliMCEvent* fMCEvent; // MC event
  Bool_t flagMC;        // Flag for MC
  Bool_t flagPileup;    // Flag for pile-up rejection
  TString fGenName;     // Generator name
  TString fYear;        // Year of data taking
  Int_t counter;
  Double_t fdeta, fdphi, fSharedFraction;
  Bool_t flag2Track;
  Bool_t flagSharity;
  Bool_t flagRejEls;
  Bool_t flagSelfAff;
  Double_t fITSCls, fTPCCls, fTPCRows;
  Bool_t flagShClPro;
  Int_t fBfield;
  Bool_t flag2TrackQA;
  Int_t mPtBins;

  // Output lists
  TList* fHistList;
  TList* fNtupleList[mPtmax];
  TList* fNtupleListGen[mPtmax];

  // Objects retrieved from the input handler
  AliPIDResponse* fPIDResponse;
  AliAODMCHeader* mcHeader;

  // Tuples for Output
  TNtuple* fntpMBin[mPtmax][mMBins];    //! Tuples for output
  TNtuple* fntpMBinGen[mPtmax][mMBins]; //! Tuples for output

  // Filter Bit
  Int_t filterBit;

  // Variables for analysis
  Double_t fVxMax, fVyMax, fVzMax, minCent, maxCent, minEta, maxEta;
  Double_t fSharedClsMax, fSharedRowsMax, fFindableClsMin;
  float fDCAxyMax, fDCAzMax;
  Int_t Mmax;
  Int_t ptbin1, ptbin2, ptbin3, ptbin4;
  TArrayD ptarray;
  TArrayI mMBin2, mNBin2;
  std::vector<Bool_t> ptbin;    // vector of Bool_tean for pt bin
  std::vector<TBits> clusmap;   // vector of TBits for cluster map
  std::vector<TBits> sharedmap; // vector of TBits for shared map

  Double_t mfield;
  Double_t nsigmaTPC[5];
  Double_t nsigmaTOF[5];
  AliPIDCombined* fPIDCombined;

  // Histograms
  TH2D* fdEtadPhiBef[mPtmax];
  TH2D* fdEtadPhiAf[mPtmax];
  TH2D* fdEtadPhiChSame[mPtmax][2][3];
  TH2D* fdEtadPhiChDiff[mPtmax][2];
  TH2D* fdEtadPhiPtOrd[mPtmax][4][3][2];
  TH1D* fHistdEta;
  TH1D* fHistdPhi;
  TH2D* fHistQAPID[13];
  TH1D* fHistPDG[2];
  TH1F* fHistnTPCcls;
  TH1F* fHistnTPCcrossedrows;
  TH1F* fHistnITScls;
  TH1F* fHistnchi2ITScls;
  TH1F* fHistnchi2TPCcls;
  TH1F* fHistDCAxy;
  TH1F* fHistDCAz;
  TH2F* fHistDCAxypT;
  TH2F* fHistDCAzpT;
  TH1F* fHistNShCls[2];
  TH2F* fHistNShClsFra[2];
  TH2F* fHistNShClsFravspt[mPtmax + 1];
  TH2F* fHistNFoundClsFra[2];
  TH1F* fHistNFcls[2];
  // For Data and also for Generated MC (if flagMC = true)
  TH1F* fPtBin[mPtmax];             //! Pt spectrum
  TH1F* fEtaBin[mPtmax];            //! Eta spectrum
  TH1F* fPhiBin[mPtmax];            //! Phi spectrum
  TH1F* fMultBin[mPtmax];           //! Histogram to register event multiplicity
  TH2D* fEtaPhiBin[mPtmax][mMBins]; //! Eta-Phi for 4 different bins  distribution

  TH1F* fPtBinGen[mPtmax];             //! Pt spectrum
  TH1F* fEtaBinGen[mPtmax];            //! Eta spectrum
  TH1F* fPhiBinGen[mPtmax];            //! Phi spectrum
  TH1F* fMultBinGen[mPtmax];           //! Histogram to register event multiplicity
  TH2D* fEtaPhiBinGen[mPtmax][mMBins]; //! Eta-Phi for 4 different bins  distribution

  // QA hists
  TH1F* fHistQAVx;     //!
  TH1F* fHistQAVy;     //!
  TH1F* fHistQAVz;     //!
  TH1D* fEventCounter; //! Histogram to track events etc
  TH1D* fTrackCounter; //! Histogram to track tracks etc
  TH1F* fHistQACent;   //! Histogram for centrality Distribution
  TH1F* fHistQAEta[2]; //!
  TH1F* fHistQAPhi[2]; //!

  AliAnalysisTaskNFactorialMoments(
    const AliAnalysisTaskNFactorialMoments&); // not implemented
  AliAnalysisTaskNFactorialMoments&
    operator=(const AliAnalysisTaskNFactorialMoments&); // not implemented

  ClassDef(AliAnalysisTaskNFactorialMoments, 1);
};

#endif