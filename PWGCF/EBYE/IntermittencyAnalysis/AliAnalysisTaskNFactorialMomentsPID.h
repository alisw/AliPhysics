#ifndef AliAnalysisTaskNFactorialMomentsPID_H
#define AliAnalysisTaskNFactorialMomentsPID_H

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "TNtuple.h"
#include "TH3F.h"

const Int_t maxBins = 40, maxDim = 2, maxqs = 6, maxPtBins = 5;

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

class AliAnalysisTaskNFactorialMomentsPID : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskNFactorialMomentsPID();
  AliAnalysisTaskNFactorialMomentsPID(const char* name);
  virtual ~AliAnalysisTaskNFactorialMomentsPID();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  // Analysis Setters
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
  void SetRejectPileup(Bool_t RejectPileup) { this->flagPileup = RejectPileup; }
  void SetBfield(Int_t bf) { this->fBfield = bf; }
  void SetDCAXYRangeMax(Double_t dcaxy) { this->fDCAxyMax = dcaxy; }
  void SetDCAZRangeMax(Double_t dcaz) { this->fDCAzMax = dcaz; }
  void SetITSClusterCut(Double_t itscut) { this->fITSCls = itscut; }
  void SetTPCClusterCut(Double_t tpccut) { this->fTPCCls = tpccut; }
  void SetnTPCrossedrows(Double_t tpcrowcut) { this->fTPCRows = tpcrowcut; }
  void SetSelfAffine(Bool_t selfaffine) { this->flagSelfAff = selfaffine; }
  void SetSharedCls(Double_t fSharedCls, Double_t fSharedrows, Double_t fFindable)
  {
    this->fSharedClsMax = fSharedCls;
    this->fSharedRowsMax = fSharedrows;
    this->fFindableClsMin = fFindable;
  }

 protected:
  AliEventCuts fEventCuts;
  void FillTrackInfo();
  void FillMCTrackInfo();
  void GetPtBin(Double_t);
  void CalculateNFMs(TH2D* h1[maxPtBins][maxBins], Bool_t mcGen);
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
  Bool_t flagSelfAff;
  Double_t fITSCls, fTPCCls, fTPCRows;
  Int_t fBfield;
  Int_t mPtBins;

  // Output lists
  TList* fHistList;
  TList* fNtupleList[maxPtBins];
  TList* fNtupleListCorr[maxPtBins];
  TList* fNtupleListGen[maxPtBins];
  TList* fNtupleListBin1; // Ntuple to store Fqes
  TList* fNtupleListBin2; // Ntuple to store Fqes
  TList* fNtupleListBin3; // Ntuple to store Fqes
  TList* fNtupleListBin4;

  TList* fWCNtupleListBin1; // Ntuple to store Fqes
  TList* fWCNtupleListBin2; // Ntuple to store Fqes
  TList* fWCNtupleListBin3; // Ntuple to store Fqes
  TList* fWCNtupleListBin4; // Ntuple to store Fqes

  // Objects retrieved from the input handler
  AliPIDResponse* fPIDResponse;
  AliAODMCHeader* mcHeader;

  // Tuples for Output
  TNtuple* fntpMBin[maxPtBins][maxBins];     //! Tuples for output
  TNtuple* fntpMBinCorr[maxPtBins][maxBins]; //! Tuples for output

  TNtuple* fntpMBinGen[maxPtBins][maxBins]; //! Tuples for output

  TNtuple* fWCntpMBin1[maxBins]; //!
  TNtuple* fWCntpMBin2[maxBins]; //!
  TNtuple* fWCntpMBin3[maxBins]; //!
  TNtuple* fWCntpMBin4[maxBins]; //!

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
  std::vector<Bool_t> ptbin; // vector of Bool_tean for pt bin

  Double_t mfield;
  Double_t nsigmaTPC[5];
  Double_t nsigmaTOF[5];
  AliPIDCombined* fPIDCombined;

  // Histograms
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
  TH2F* fHistNShClsFravspt;
  TH2F* fHistNFoundClsFra[2];
  TH1F* fHistNFcls[2];
  // For Data and also for Generated MC (if flagMC = true)
  TH1F* fPtBin[maxPtBins];              //! Pt spectrum
  TH1F* fEtaBin[maxPtBins];             //! Eta spectrum
  TH1F* fPhiBin[maxPtBins];             //! Phi spectrum
  TH1F* fMultBin[maxPtBins];            //! Histogram to register event multiplicity
  TH2D* fEtaPhiBin[maxPtBins][maxBins]; //! Eta-Phi for 4 different bins  distribution

  TH1F* fPtBinGen[maxPtBins];              //! Pt spectrum
  TH1F* fEtaBinGen[maxPtBins];             //! Eta spectrum
  TH1F* fPhiBinGen[maxPtBins];             //! Phi spectrum
  TH1F* fMultBinGen[maxPtBins];            //! Histogram to register event multiplicity
  TH2D* fEtaPhiBinGen[maxPtBins][maxBins]; //! Eta-Phi for 4 different bins  distribution

  TH2D* fhMapEtaPhiBin1M[maxBins]; //!
  TH2D* fhMapEtaPhiBin2M[maxBins]; //!
  TH2D* fhMapEtaPhiBin3M[maxBins]; //!
  TH2D* fhMapEtaPhiBin4M[maxBins]; //!

  // QA hists
  TH1F* fHistQAVx;     //!
  TH1F* fHistQAVy;     //!
  TH1F* fHistQAVz;     //!
  TH1D* fEventCounter; //! Histogram to track events etc
  TH1D* fTrackCounter; //! Histogram to track tracks etc
  TH1F* fHistQACent;   //! Histogram for centrality Distribution
  TH1F* fHistQAEta[2]; //!
  TH1F* fHistQAPhi[2]; //!
  Bool_t zeromult;

  AliAnalysisTaskNFactorialMomentsPID(
    const AliAnalysisTaskNFactorialMomentsPID&); // not implemented
  AliAnalysisTaskNFactorialMomentsPID&
    operator=(const AliAnalysisTaskNFactorialMomentsPID&); // not implemented

  ClassDef(AliAnalysisTaskNFactorialMomentsPID, 1);
};

#endif
