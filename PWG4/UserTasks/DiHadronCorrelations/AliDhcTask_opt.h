// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011


#ifndef AliDhcTask_cxx
#define AliDhcTask_cxx

class TFormula;
class TH1;
class TH2;
class TObjArray;
class TObject;
class TProfile2D;
class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class KiddiePoolManager;

#include "AliAnalysisTaskSE.h"
#include "KiddiePoolClasses.h"

class AliDhcTask : public AliAnalysisTaskSE {
 public:
  AliDhcTask() : 
    AliAnalysisTaskSE(), fVerbosity(0), fEtaMax(1), fZVtxMax(10), fPtMin(0.25), fPtMax(15),
    fESD(0), fAOD(0), fOutputList(0), fHistPt(0), fHEvt(0), fHTrk(0), fHPtAss(0), 
    fHPtTrg(0), fHCent(0), fHZvtx(0), fNbins(0), fHSs(0), fHMs(0), fIndex(0), 
    fCentrality(99), fZVertex(99), fEsdTrackCutsTPCOnly(0), fPoolMgr(0) {}

  AliDhcTask(const char *name);
  virtual ~AliDhcTask() {}
  
  void         SetVerbosity(Int_t v)                  { fVerbosity = v;         }
  void         SetPtRange(Double_t min, Double_t max) { fPtMin=min; fPtMax=max; }
  void         SetEtaMax(Double_t eta)                { fEtaMax = eta;          }
  void         SetZvtx(Double_t zvtx)                 { fZVtxMax = zvtx;        }

 protected:
  enum ePairingScheme {kSameEvt, kDiffEvt};
  enum eDataType      {kESD, kAOD};

  void         BookHistos();
  void         InitEventMixer();
  MiniEvent*   GetESDTrax() const;
  MiniEvent*   GetAODTrax() const;
  Bool_t       VertexOk(TObject* obj) const;
  Double_t     DeltaPhi(Double_t phia, Double_t phib, 
			Double_t rangeMin = -TMath::Pi()/2, 
			Double_t rangeMax = 3*TMath::Pi()/2) const;
  Int_t        Correlate(const MiniEvent &arr1, const MiniEvent &arr2, 
			 Int_t pairing = kSameEvt, Double_t weight = 1.);
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

 private:
  Int_t        fVerbosity;       //  0 = silence
  Double_t     fEtaMax;          //  Max |eta| cut (cm)
  Double_t     fZVtxMax;         //  Max |z| cut (cm)
  Double_t     fPtMin;           //  Min pt cut
  Double_t     fPtMax;           //  Max pt cut
  AliESDEvent *fESD;             //! ESD object
  AliAODEvent *fAOD;             //! AOD object
  TList       *fOutputList;      //! Output list
  TH1F        *fHistPt;          //! Pt spectrum
  TH2         *fHEvt;            //! Cent, vtx, etc.
  TH2         *fHTrk;            //! Phi, Eta, etc.
  TH1         *fHPtAss;          //! Pt ass 
  TH1         *fHPtTrg;          //! Pt trg
  TH1         *fHCent;           //! Centrality
  TH1         *fHZvtx;           //! Zvertex
  Int_t        fNbins;           //! Number of histogram bins
  TH2        **fHSs;             //! Same-evt correlations
  TH2        **fHMs;             //! Diff-evt correlations
  TFormula    *fIndex;           //! Index for histograms
  TProfile2D **fMeanPtTrg;       //! Mean pt trig 
  TProfile2D **fMeanPtAss;       //! Mean pt ass
  TProfile2D **fMean2PtTrg;      //! RMS pt trig 
  TProfile2D **fMean2PtAss;      //! RMS pt ass
  Double_t     fCentrality;      //! V0M for now
  Double_t     fZVertex;         //! Of current event
  AliESDtrackCuts   *fEsdTrackCutsTPCOnly; //! Track cuts
  KiddiePoolManager *fPoolMgr;             //! Event mixer

  AliDhcTask(const AliDhcTask&);            // not implemented
  AliDhcTask &operator=(const AliDhcTask&); // not implemented

  ClassDef(AliDhcTask, 2);
};

#endif
