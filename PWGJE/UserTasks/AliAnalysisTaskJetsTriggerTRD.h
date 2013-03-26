#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"

#define ID(x) x, #x

class TList;

class AliAnalysisTaskJetsTriggerTRD :
  public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskJetsTriggerTRD(const char *name = "jets_trg_trd");
  ~AliAnalysisTaskJetsTriggerTRD();

  // analysis operations
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *option);

  // task configuration
  void  SetNumberOfJetPtBins(Int_t n) { fNoJetPtBins = n; }
  Int_t GetNumberOfJetPtBins() const { return fNoJetPtBins; }

  void    SetJetPtBinMax(Float_t ptmax) { fJetPtBinMax = ptmax; }
  Float_t GetJetPtBinMax() const { return fJetPtBinMax; }

  void SetJetBranchName(const char* const branchName) { strncpy(fJetBranchName, branchName, fgkStringLength-1); }
  const char* GetJetBranchName() const { return fJetBranchName; }

  // histograms
  enum Hist_t {
      kHistStat = 0,
      kHistNoJets,
      kHistTrackGTU,
      kHistJetPt, kHistJetPtITS, kHistJetPt3x3,
      kHistJetPtEMC, kHistJetPtHJT,
      kHistJetPtNoTracks3,
      kHistLast
  };

  // statistics
  enum Stat_t {
      kStatSeen = 1,
      kStatUsed,
      kStatMB,
      kStatLast
  };

  // trigger conditions
  enum Trigger_t { 
      kTrgMB = 0,
      kTrgInt,
      kTrgInt78,
      kTrgHJT,
      kTrgEMC,
      kTrgLast
  };	   

protected:
  UInt_t fTriggerMask;		// internal representation of trigger conditions

  Bool_t DetectTriggers();
  void   MarkTrigger(Trigger_t trg) { fTriggerMask |= (1 << trg); }
  Bool_t IsTrigger(Trigger_t trg) const { return (fTriggerMask & (1 << trg)); }

  // output objects
  TList *fOutputList;		// list of output objects

  // histogram management
  TH1  *fHist[kHistLast];	//! pointers to histogram
  const char *fShortTaskId;	//! short identifier for the task

  TH1*&  GetHistogram(Hist_t hist, const Int_t idx = 0) { return fHist[hist + idx]; }

  TH1*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax, Int_t binType = 1);
  TH2*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax, Int_t binType = 1);
  TH3*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax,
                      Int_t zbins, Float_t zmin, Float_t zmax, Int_t binType = 1);

  void    FillH1(Hist_t hist, Float_t x, Float_t weight = 1., Int_t idx = 0)
  { GetHistogram(hist, idx)->Fill(x, weight); }
  void    FillH2(Hist_t hist, Float_t x, Float_t y, Float_t weight = 1., Int_t idx = 0)
  { ((TH2*) GetHistogram(hist, idx))->Fill(x, y, weight); }
  void    FillH3(Hist_t hist, Float_t x, Float_t y, Float_t z, Float_t weight = 1., Int_t idx = 0)
  { ((TH3*) GetHistogram(hist, idx))->Fill(x, y, z, weight); }

  // task configuration
  Int_t    fNoJetPtBins;                // number of bins for jet pt
  Float_t  fJetPtBinMax;                // max jet pt (GeV) in histograms

  static const Int_t fgkStringLength = 100; // max length for the jet branch name
  char fJetBranchName[fgkStringLength];     // jet branch name

  // not implemented
  AliAnalysisTaskJetsTriggerTRD(const AliAnalysisTaskJetsTriggerTRD &rhs);
  AliAnalysisTaskJetsTriggerTRD& operator=(const AliAnalysisTaskJetsTriggerTRD &rhs);

  ClassDef(AliAnalysisTaskJetsTriggerTRD, 1);
};
