// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011, updated Mar 2012

#ifndef AliDhcTask_cxx
#define AliDhcTask_cxx

#include "AliAnalysisTaskSE.h"
#include "AliPool.h"
#include "THn.h"

class TFormula;
class TH1;
class TH2;
class TH3;
class TAxis;
class TObjArray;
class TObject;
class TProfile2D;
class AliAODEvent;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDMuonTrack;
class AliAODTrack;
class AliEvtPoolManager;


class AliDhcTask : public AliAnalysisTaskSE {
 public:
  AliDhcTask();
  AliDhcTask(const char *name);
  virtual ~AliDhcTask();
  
  void         SetCentBins(TAxis *bins)               { fBCent=bins;              }
  void         SetCentMethod(const char *name)        { fCentMethod = name;       }
  void         SetCentMixBins(TAxis *bins)            { fMixBCent=bins;           }
  void         SetDEtaDPhiBins(Int_t nbe, Int_t nbp)  { fNBdeta=nbe; fNBdphi=nbp; }
  void         SetDoWeights(Bool_t b)                 { fDoWeights = b;           }
  void         SetFillMuons(Bool_t b)                 { fFillMuons = b;           }
  void         SetPtTACrit(Bool_t b)                  { fPtTACrit = b;            }
  void         SetEtaMax(Double_t eta)                { fEtaMax = eta;            }
  void         SetPoolSize(Int_t p)                   { fPoolSize = p;            }
  void         SetPtABins(TAxis *bins)                { fBPtA=bins;               }
  void         SetPtRange(Double_t min, Double_t max) { fPtMin=min; fPtMax=max;   }
  void         SetPtTBins(TAxis *bins)                { fBPtT=bins;               }
  void         SetTrackDepth(Int_t t)                 { fTrackDepth = t;          }
  void         SetTracksName(const char *n)           { fTracksName = n;          }
  void         SetVerbosity(Int_t v)                  { fVerbosity = v;           }
  void         SetZVtxBins(TAxis *bins)               { fBZvtx=bins;              }
  void         SetZVtxMixBins(TAxis *bins)            { fMixBZvtx=bins;           }
  void         SetZvtx(Double_t zvtx)                 { fZVtxMax = zvtx;          }
  void         SetHEffT(THnF *h)                      { fHEffT=h;                 }
  void         SetHEffA(THnF *h)                      { fHEffA=h;                 }
  void         SetAnaMode(Int_t iAna);
  enum eAnaMode       {kHH, kMuH, kHMu, kMuMu, kPSide, kASide};

 protected:
  enum ePairingScheme {kSameEvt, kDiffEvt};
  enum eDataType      {kESD, kAOD};

  void         BookHistos();
  void         InitEventMixer();
  void         GetESDTracks(MiniEvent*);
  void         GetAODTracks(MiniEvent*);
  Bool_t       VertexOk(TObject* obj) const;
  Bool_t       IsGoodMUONtrack(AliESDMuonTrack &track);
  Bool_t       IsGoodMUONtrack(AliAODTrack &track);
  Double_t     DeltaPhi(Double_t phia, Double_t phib,
			Double_t rangeMin = -TMath::Pi()/2, 
			Double_t rangeMax = 3*TMath::Pi()/2) const;
  Int_t        Correlate(const MiniEvent &arr1, const MiniEvent &arr2, Int_t pairing = kSameEvt);
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

 private:
  Int_t              fVerbosity;       //  0 = silence
  Double_t           fEtaMax;          //  Max |eta| cut for standard ESD or AOD analysis
  Double_t           fZVtxMax;         //  Max |z| cut (cm)
  Double_t           fPtMin;           //  Min pt cut
  Double_t           fPtMax;           //  Max pt cut
  Int_t              fTrackDepth;      //  #tracks to fill pool
  Int_t              fPoolSize;        //  Maximum number of events
  TString            fTracksName;      //  name of track collection
  Bool_t             fDoWeights;       //  if true weight with 1/N per event
  Bool_t             fFillMuons;       //  fill the muon tracks into the mini event
  Bool_t             fPtTACrit;        //  use the pTT > pTA criterion?
  Double_t           fEtaTLo;          //  Min eta for triggers
  Double_t           fEtaTHi;          //  Max eta for triggers
  Double_t           fEtaALo;          //  Min eta for associated
  Double_t           fEtaAHi;          //  Max eta for associated
  AliESDEvent       *fESD;             //! ESD object
  AliAODEvent       *fAOD;             //! AOD object
  TList             *fOutputList;      //! Output list
  TH2               *fHEvt;            //! Cent vs vtx.
  TH2               *fHTrk;            //! Phi vs Eta
  TH1               *fHPtAss;          //! Pt ass 
  TH1               *fHPtTrg;          //! Pt trg
  TH1               *fHPtTrgEvt;       //! Pt trg per event for weighting
  TH3               *fHPtTrgNorm1S;    //! Pt trg same events in cent. and zvtx bins, method 1
  TH3               *fHPtTrgNorm1M;    //! Pt trg mixed events in cent. and zvtx bins, method 1
  TH3               *fHPtTrgNorm2S;    //! Pt trg same events in cent. and zvtx bins, method 2
  TH3               *fHPtTrgNorm2M;    //! Pt trg mixed events in cent. and zvtx bins, method 2
  TH1               *fHCent;           //! Centrality
  TH1               *fHZvtx;           //! Zvertex
  Int_t              fNbins;           //! Number of histogram bins
  TH2              **fHSs;             //! Same-evt correlations
  TH2              **fHMs;             //! Diff-evt correlations
  TH1              **fHPts;            //! Pt distributions
  TH3               *fHQAT;            //! trigger particle distribution for QA
  TH3               *fHQAA;            //! associated particle distribution for QA
  TFormula          *fIndex;           //! Index for histograms
  Double_t           fCentrality;      //! V0M for now
  Double_t           fZVertex;         //! Of current event
  AliESDtrackCuts   *fEsdTPCOnly;      //! Track cuts
  AliEvtPoolManager *fPoolMgr;         //! Event mixer
  TString            fCentMethod;      //  centrality selection method
  Int_t              fNBdeta;          //  no. deta bins
  Int_t              fNBdphi;          //  no. dphi bins
  TAxis             *fBPtT;            //  ptt binning
  TAxis             *fBPtA;            //  pta binning
  TAxis             *fBCent;           //  centrality binning
  TAxis             *fBZvtx;           //  zvtx binning
  TAxis             *fMixBCent;        //  centrality binning for mixing
  TAxis             *fMixBZvtx;        //  zvtx binning for mixing
  THnF              *fHEffT;           //  efficiency for trigger particles
  THnF              *fHEffA;           //  efficiency for associate particles

  AliDhcTask(const AliDhcTask&);            // not implemented
  AliDhcTask &operator=(const AliDhcTask&); // not implemented

  ClassDef(AliDhcTask, 3);
};

#endif
