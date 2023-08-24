#ifndef AliAnalysisTaskNFactorialMoments_H
#define AliAnalysisTaskNFactorialMoments_H

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

const int M = 40, Q = 6, fPt = 4, D = 2;

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
  void SetTwoTrackCuts(double deta, double dphi, bool twotrack)
  {
    this->fdeta = deta;
    this->fdphi = dphi;
    this->ftwotrack = twotrack;
  }
  void SetSharingFraction(double fshfr, bool sharity)
  {
    this->fSharedFraction = fshfr;
    this->fsharity = sharity;
  }
  void SetVtxCut(double fVx, double fVy, double fVz)
  {
    this->fVxMax = fVx;
    this->fVyMax = fVy;
    this->fVzMax = fVz;
  }
  void SetIsMC(bool isMC, TString gen)
  {
    this->fismc = isMC;
    this->fGenName = gen;
  }
  void SetCentLim(int f_minCent, int f_maxCent)
  {
    this->minCent = f_minCent;
    this->maxCent = f_maxCent;
  }
  void SetEta(double f_minEta, double f_maxEta)
  {
    this->minEta = f_minEta;
    this->maxEta = f_maxEta;
  }
  void Setfbit(int fb) { this->filterBit = fb; }
  void SetYear(TString year) { this->fYear = year; }
  void SetPSbinning(TArrayI MB, TArrayI NB, int Mm)
  {
    this->Mbins = MB;
    this->Nbins = NB;
    this->Mmax = Mm;
  }
  void SetPtArray(TArrayD f_ptArray) { this->ptarray = f_ptArray; }
  void SetRejectElectrons(bool fRejectElectron)
  {
    this->fRejectEls = fRejectElectron;
  }
  void SetRejectPileup(bool RejectPileup) { this->fpileup = RejectPileup; }
  void SetBfield(int bf) { this->fBfield = bf; }
  void SetDCAXYRangeMax(double dcaxy) { this->fDCAxyMax = dcaxy; }
  void SetDCAZRangeMax(double dcaz) { this->fDCAzMax = dcaz; }
  void SetITSClusterCut(double itscut) { this->fITSCls = itscut; }
  void SetTPCClusterCut(double tpccut) { this->fTPCCls = tpccut; }
  void SetnTPCrossedrows(double tpcrowcut) { this->fTPCRows = tpcrowcut; }
  void SetSelfAffine(bool selfaffine) { this->fSelfAffine = selfaffine; }
  void SetPIDAnalysis(bool prAnalysis, double nSigCut)
  {
    this->fProtonAnalysis = prAnalysis;
    this->nSigmaPrCut = nSigCut;
  }
  void SetSharedCls(double fSharedCls, double fSharedrows, double fFindable, bool fSharedMean)
  {
    this->fSharedClsMax = fSharedCls;
    this->fSharedRowsMax = fSharedrows;
    this->fFindableClsMin = fFindable;
    this->fSharedFrMean = fSharedMean;
  }

 protected:
  AliEventCuts fEventCuts;
  void FillTrackInfo();
  void FillMCTrackInfo();
  float GetDPhiStar(float phi1, float pt1, float charge1, float phi2, float pt2,
                    float charge2, float radius, float bSign);
  bool GetParticleID(AliAODTrack* track, bool fQA);
  float dphistarcalculation(float phi1, float eta1, float pt1, int charge1,
                            float phi2, float eta2, float pt2, int charge2,
                            float bSign);
  float SharedClusterFraction(TBits&, TBits&, TBits&, TBits&);
  void GetPtBin(double);
  void CalculateNFMs(TH2D* h1[M], TH2D* h2[M], TH2D* h3[M], TH2D* h4[M],
                     bool mcGen);
  bool ConfigureTrack();
  void DataPosting();
  void ResetHistograms();

 private:
  // Event related (Data and MC)
  AliAODEvent* fAOD;    // AOD object
  AliMCEvent* fMCEvent; // MC event
  bool fismc;           // Flag for MC
  bool fpileup;         // Flag for pile-up rejection
  TString fGenName;     // Generator name
  TString fYear;        // Year of data taking
  int counter;
  double fdeta, fdphi, fSharedFraction;
  bool ftwotrack;
  bool fsharity;
  bool fRejectEls;
  bool fSelfAffine;
  double fITSCls, fTPCCls, fTPCRows;
  bool fSharedFrMean;
  int fBfield;
  // int fPt;              // Number of pt bins

  // Output lists
  TList* fOutHList;
  TList* fQAList;
  TList* fQAList2;
  TList* fNtupleListBin1;
  TList* fNtupleListBin2;
  TList* fNtupleListBin3;
  TList* fNtupleListBin4;
  TList* fNtupleListBin[fPt];
  TList* fNtupleListBin1Gen;
  TList* fNtupleListBin2Gen;
  TList* fNtupleListBin3Gen;
  TList* fNtupleListBin4Gen;

  // Objects retrieved from the input handler
  AliPIDResponse* fPIDResponse;
  AliAODMCHeader* mcHeader;

  // Tuples for Output
  TNtuple* fntpMBin1[M];     //!
  TNtuple* fntpMBin2[M];     //!
  TNtuple* fntpMBin3[M];     //!
  TNtuple* fntpMBin4[M];     //!
  TNtuple* fntpMBin[M][fPt]; //!
  TNtuple* fntpMBin1Gen[M];  //!
  TNtuple* fntpMBin2Gen[M];  //!
  TNtuple* fntpMBin3Gen[M];  //!
  TNtuple* fntpMBin4Gen[M];  //!

  // Filter Bit
  int filterBit;

  // Variables for analysis
  double fVxMax, fVyMax, fVzMax, minCent, maxCent, minEta, maxEta;
  double fSharedClsMax, fSharedRowsMax, fFindableClsMin;
  float fDCAxyMax, fDCAzMax;
  int Mmax;
  int ptbin1, ptbin2, ptbin3, ptbin4;
  TArrayD ptarray;
  TArrayI Mbins, Nbins;
  std::vector<bool> ptbin;      // vector of boolean for pt bin
  std::vector<TBits> clusmap;   // vector of TBits for cluster map
  std::vector<TBits> sharedmap; // vector of TBits for shared map

  Double_t mfield;
  Double_t nsigmaTPC[5];
  Double_t nsigmaTOF[5];
  Bool_t fProtonAnalysis;
  Double_t nSigmaPrCut;
  AliPIDCombined* fPIDCombined;

  // Histograms
  TH2D* fHistbeforeHBT[4];
  TH2D* fHistafterHBT[4];
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
  TH1F* fHistnsharedcls[2];
  TH2F* fHistnshclsfra[2];
  TH2F* fHistnshclsfravspt[5];
  TH2F* fHistnfoundcls[2];
  TH1F* fHistnfcls[2];
  // For Data and also for Generated MC (if fismc = true)
  TH1F* fHistPtBin[4];   //! Pt spectrum
  TH1F* fEtaBin[4];      //! Eta spectrum
  TH1F* fPhiBin[4];      //! Phi spectrum
  TH1F* fHistMulBin[4];  //! Histogram to register event multiplicity
  TH2D* fHEtaPhiBin1[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin2[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin3[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin4[M]; //! Eta-Phi for 10 different bins  distribution
  // Extra Histograms for  Reconstucted MC
  TH1F* fHistPtBinGen[4];   //! Pt spectrum
  TH1F* fEtaBinGen[4];      //! Eta spectrum
  TH1F* fPhiBinGen[4];      //! Phi spectrum
  TH1F* fHistMulBinGen[4];  //! Histogram to register event multiplicity
  TH2D* fHEtaPhiBin1Gen[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin2Gen[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin3Gen[M]; //! Eta-Phi for 10 different bins  distribution
  TH2D* fHEtaPhiBin4Gen[M]; //! Eta-Phi for 10 different bins  distribution

  // QA hists
  TH1F* fHistQAVx;      //!
  TH1F* fHistQAVy;      //!
  TH1F* fHistQAVz;      //!
  TH1D* fEventCounter;  //! Histogram to track events etc
  TH1D* fTrackCounter;  //! Histogram to track tracks etc
  TH1F* fHistQACent;    //! Histogram for centrality Distribution
  TH1F* fHistQAEta[2];  //!
  TH1F* fHistQAPhi[2];  //!
  TH1F* fHistQAMult;    //!
  TH1F* fHistQAMultwpc; //!

  AliAnalysisTaskNFactorialMoments(
    const AliAnalysisTaskNFactorialMoments&); // not implemented
  AliAnalysisTaskNFactorialMoments&
    operator=(const AliAnalysisTaskNFactorialMoments&); // not implemented

  ClassDef(AliAnalysisTaskNFactorialMoments, 1);
};

#endif
