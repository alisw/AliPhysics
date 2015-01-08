// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011, updated Mar 2012

#ifndef AliDhcTask_cxx
#define AliDhcTask_cxx

#include "AliAnalysisTaskSE.h"
#include "AliPool.h"
#include "AliMuonTrackCuts.h"

class TFormula;
class TH1;
class TH2;
class TH3;
class THn;
class TAxis;
class TObjArray;
class TObject;
class TProfile2D;
class AliAODEvent;
class AliAODTrack;
class AliAnalysisUtils;
class AliESDEvent;
class AliESDMuonTrack;
class AliESDtrackCuts;
class AliEvtPoolManager;

class AliDhcTask : public AliAnalysisTaskSE {
 public:
  AliDhcTask();
  AliDhcTask(const char *name, Bool_t def=kFALSE);
  virtual ~AliDhcTask();
  
  void         SetAllTAHists(Bool_t b)                { fAllTAHists = b;          }
  void         SetAnaMode(Int_t iAna);
  void         SetCentBins(TAxis *bins)               { fBCent=bins;              }
  void         SetCentMethod(const char *name)        { fCentMethod = name;       }
  void         SetCentMixBins(TAxis *bins)            { fMixBCent=bins;           }
  void         SetCheckVertex(Bool_t b)               { fCheckVertex = b;         }
  void         SetClassName(const char *n)            { fClassName = n;           }
  void         SetDEtaDPhiBins(Int_t nbe, Int_t nbp)  { fNBdeta=nbe; fNBdphi=nbp; }
  void         SetDoFillSame(Bool_t b)                { fDoFillSame = b;          }
  void         SetDoMassCut(Bool_t b)                 { fDoMassCut = b;           }
  void         SetDoWeights(Bool_t b)                 { fDoWeights = b;           }
  void         SetEtaMax(Double_t eta)                { fEtaMax = eta;            }
  void         SetEtaTRange(Double_t eL, Double_t eH) { fEtaTLo=eL; fEtaTHi=eH;   }
  void         SetEtaARange(Double_t eL, Double_t eH) { fEtaALo=eL; fEtaAHi=eH;   }
  void         SetFillMuons(Bool_t b)                 { fFillMuons = b;           }
  void         SetHEffA(THnF *h)                      { fHEffA=h;                 }
  void         SetHEffT(THnF *h)                      { fHEffT=h;                 }
  void         SetMixInEtaT(Bool_t b)                 { fMixInEtaT = b;           }
  void         SetMuonUtils(Bool_t b = kTRUE,
                            UInt_t mask = AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchApt ) {
    fUseMuonUtils = b;
    fMuonCutMask = mask;
  }
  void         SetOmitFirstEv(Bool_t b)               { fOmitFirstEv = b;         }
  void         SetPoolSize(Int_t p)                   { fPoolSize = p;            }
  void         SetPtABins(TAxis *bins)                { fBPtA=bins;               }
  void         SetPtRange(Double_t min, Double_t max) { fPtMin=min; fPtMax=max;   }
  void         SetPtTACrit(Bool_t b)                  { fPtTACrit = b;            }
  void         SetPtTBins(TAxis *bins)                { fBPtT=bins;               }
  void         SetRequireMuon(Bool_t b, Double_t l=0.5, Double_t h=4.0) {
    fRequireMuon = b;
    fReqPtLo = l;
    fReqPtHi = h;
  }
  void         SetTrackDepth(Int_t t)                 { fTrackDepth = t;          }
  void         SetTracksName(const char *n)           { fTracksName = n;          }
  void         SetTriggerMatch(Bool_t b)              { fTriggerMatch = b;        }
  void         SetVerbosity(Int_t v)                  { fVerbosity = v;           }
  void         SetZVtxBins(TAxis *bins)               { fBZvtx=bins;              }
  void         SetZVtxMixBins(TAxis *bins)            { fMixBZvtx=bins;           }
  void         SetZvtx(Double_t zvtx)                 { fZVtxMax = zvtx;          }
  void         PrintDhcSettings();

  enum eAnaMode       {kHH=1, kMuH, kHMu, kMuMu, kPSide, kASide};

 protected:
  enum ePairingScheme {kSameEvt, kDiffEvt};
  enum eDataType      {kESD, kAOD};

  void         BookHistos();
  void         InitEventMixer();
  void         GetESDTracks(MiniEvent*);
  void         GetAODTracks(MiniEvent*);
  Bool_t       VertexOk() const;
  Bool_t       HasMuonESD();
  Bool_t       HasMuonAOD();
  Bool_t       IsGoodMUONtrack(AliESDMuonTrack &track);
  Bool_t       IsGoodMUONtrack(AliAODTrack &track);
  Bool_t       IsTrigger(Double_t eta, Double_t pt);
  Bool_t       IsAssociated(Double_t eta, Double_t pt);
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
  Bool_t             fRequireMuon;     //  only run on events with a muon track
  Double_t           fReqPtLo;         //  require a muon above this pt
  Double_t           fReqPtHi;         //  and below this pt
  Bool_t             fPtTACrit;        //  use the pTT > pTA criterion?
  Bool_t             fAllTAHists;      //  create all pTT,pTA combination hists, even t<a?
  Bool_t             fMixInEtaT;       //  mix in bins of eta_T instead of z_vertex
  Bool_t             fUseMuonUtils;    //  use muon cuts from PWG/muon/AliMuonTrackCuts ?
  UInt_t             fMuonCutMask;     //  muon cut mask for above
  AliMuonTrackCuts  *fMuonTrackCuts;   //  muon track cut object
  Double_t           fEtaTLo;          //  Min eta for triggers
  Double_t           fEtaTHi;          //  Max eta for triggers
  Double_t           fEtaALo;          //  Min eta for associated
  Double_t           fEtaAHi;          //  Max eta for associated
  Bool_t             fOmitFirstEv;     //  if true skip first event in chunk
  Bool_t             fDoFillSame;      //  If true fill same event immediately (not waiting for pool)
  Bool_t             fDoMassCut;       //  If true cut on invariant mass
  Bool_t             fCheckVertex;     //  switch to flase for MC generator-level analysis
  TString            fClassName;       //  If not null only process events with given class
  TString            fCentMethod;      //  centrality selection method
  Int_t              fNBdeta;          //  no. deta bins
  Int_t              fNBdphi;          //  no. dphi bins
  Bool_t             fTriggerMatch;    //  muon trigger match
  TAxis             *fBPtT;            //  ptt binning
  TAxis             *fBPtA;            //  pta binning
  TAxis             *fBCent;           //  centrality binning
  TAxis             *fBZvtx;           //  zvtx binning
  TAxis             *fMixBCent;        //  centrality binning for mixing
  TAxis             *fMixBZvtx;        //  zvtx binning for mixing
  THnF              *fHEffT;           //  efficiency for trigger particles
  THnF              *fHEffA;           //  efficiency for associate particles
  AliESDEvent       *fESD;             //! ESD object
  AliAODEvent       *fAOD;             //! AOD object
  TList             *fOutputList;      //! Output list
  TH2               *fHEvt;            //! Cent vs vtx.
  TH2               *fHEvtWTr;         //! Cent vs vtx. for events with at least one track
  TH2               *fHTrk;            //! Phi vs Eta
  TH2               *fHPoolReady;      //! Check how many Jobs start mixing
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
  TH1              **fHSMass;          //! Mass distributions same
  TH1              **fHMMass;          //! Mass distributions mixed
  TH3               *fHQATp;           //! positive trigger particle distribution for QA
  TH3               *fHQAAp;           //! associated particle distribution for QA
  TH3               *fHQATpCorr;       //! correction applied trigger distribution
  TH3               *fHQAApCorr;       //! correction applied associated distribution
  TH3               *fHQATm;           //! the same for negative particles
  TH3               *fHQAAm;           //!
  TH3               *fHQATmCorr;       //!
  TH3               *fHQAAmCorr;       //!
  TH2               *fHPtCentT;        //! trigger pT vs centrality
  TH2               *fHPtCentA;        //! associated pT vs centrality
  TFormula          *fIndex;           //! index for histograms
  Double_t           fCentrality;      //! centrality of current event
  Double_t           fZVertex;         //! z vertex of current event
  AliESDtrackCuts   *fEsdTPCOnly;      //! track cuts
  AliEvtPoolManager *fPoolMgr;         //! event mixer
  AliAnalysisUtils  *fUtils;           //! analysis utils

  AliDhcTask(const AliDhcTask&);            // not implemented
  AliDhcTask &operator=(const AliDhcTask&); // not implemented

  ClassDef(AliDhcTask, 9)
};

#endif
