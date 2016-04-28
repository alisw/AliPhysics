#ifndef ALIANALYSISTASKJETSEVSHAPE
#define ALIANALYSISTASKJETSEVSHAPE

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

#include "AliLog.h"

#include "AliVParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

#define ID(x) x, #x
#define LAB(x) x + 1, #x

class TList;
class TClonesArray;
class AliVTrack;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

class AliAnalysisTaskJetsEvshape :
  public AliAnalysisTaskEmcalJet
{
public:
  AliAnalysisTaskJetsEvshape(const char *name = "jets_trg_trd");
  virtual ~AliAnalysisTaskJetsEvshape();

  // analysis operations
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *option);
  virtual void   PrintTask(Option_t *option, Int_t indent) const;

  // output lists
  enum OutputList_t {
    kOutputEmcal = 1,
    kOutputTask = 2
  };

  // histograms
  enum Hist_t {
      kHistStat = 0,
      kHistJetPt,
      kHistMult,
      kHistJetPtVsMult,
      kHistLast
  };

  // statistics
  enum Stat_t {
      kStatSeen = 1,
      kStatTrg,
      kStatEvCuts,
      kStatUsed,
      kStatLast
  };

protected:
  // from EMCAL framework
  virtual void   ExecOnce();
  virtual Bool_t FillHistograms();
  virtual Bool_t Run();

  // task internal
  AliMCEvent  *fMCEvent; //!
  AliESDEvent *fESDEvent; //!
  AliAODEvent *fAODEvent; //!

  Int_t fRunNumber; //! current run number

  Bool_t PrepareEvent();
  Bool_t CleanUpEvent();

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters

  // output objects
  TList *fOutputList;		//! list of output objects

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

private:
  // not implemented
  AliAnalysisTaskJetsEvshape(const AliAnalysisTaskJetsEvshape &rhs);
  AliAnalysisTaskJetsEvshape& operator=(const AliAnalysisTaskJetsEvshape &rhs);

  ClassDef(AliAnalysisTaskJetsEvshape, 1);
};
#endif
