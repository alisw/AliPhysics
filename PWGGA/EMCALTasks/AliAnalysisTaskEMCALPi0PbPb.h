#ifndef AliAnalysisTaskEMCALPi0PbPb_h
#define AliAnalysisTaskEMCALPi0PbPb_h

// $Id$

class TAxis;
class TClonesArray;
class TH1;
class TH2;
class TNtuple;
class TObjArray;
class AliAODCaloCells;
class AliAODCaloCluster;
class AliAODEvent;
class AliAODTrack;
class AliAODVertex;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDCaloCells;
class AliESDCaloCluster;
class AliESDEvent;
class AliESDTrack;
class AliESDVertex;
class AliESDtrackCuts;
class AliMCEvent;
class AliMCParticle;
class AliStaHeader;
class AliStaVertex;

#include "AliAnalysisTaskSE.h"
#include "AliStaObjects.h"

class AliAnalysisTaskEMCALPi0PbPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPi0PbPb();
  AliAnalysisTaskEMCALPi0PbPb(const char *name);
  virtual ~AliAnalysisTaskEMCALPi0PbPb(); 
  
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  void         SetAsymMax(Double_t asymMax)                   { fAsymMax       = asymMax;   }
  void         SetCentrality(const char *n)                   { fCentVar       = n;         }
  void         SetCentralityRange(Double_t from, Double_t to) { fCentFrom=from; fCentTo=to; }
  void         SetClusName(const char *n)                     { fClusName      = n;         }
  void         SetDoAfterburner(Bool_t b)                     { fDoAfterburner = b;         }
  void         SetDoPhysicsSelection(Bool_t b)                { fDoPSel        = b;         }
  void         SetDoTrackMatWithGeom(Bool_t b)                { fDoTrMatGeom   = b;         }
  void         SetEmbedMode(Bool_t b)                         { fEmbedMode     = b;         }
  void         SetFillNtuple(Bool_t b)                        { fDoNtuple      = b;         }
  void         SetGeoName(const char *n)                      { fGeoName       = n;         }
  void         SetGeoUtils(AliEMCALGeometry *geo)             { fGeom          = geo;       }
  void         SetIsoDist(Double_t d)                         { fIsoDist       = d;         }
  void         SetL0TimeRange(Int_t l, Int_t h)               { fMinL0Time=l; fMaxL0Time=h; }
  void         SetMarkCells(const char *n)                    { fMarkCells     = n;         }
  void         SetMcMode(Bool_t b)                            { fMcMode        = b;         }
  void         SetMinClusEnergy(Double_t e)                   { fMinE          = e;         }
  void         SetMinEcc(Double_t ecc)                        { fMinEcc        = ecc;       }
  void         SetMinErat(Double_t erat)                      { fMinErat       = erat;      }
  void         SetMinNClustersPerTrack(Double_t m)            { fMinNClusPerTr = m;         }
  void         SetNminCells(Int_t n)                          { fNminCells     = n;         }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)           { fPrimTrCuts    = c;         }
  void         SetPrimTracksName(const char *n)               { fPrimTracksName = n;        }
  void         SetRecoUtils(AliEMCALRecoUtils *reco)          { fReco          = reco;      }
  void         SetTrClassNames(const char *n)                 { fTrClassNames  = n;         }
  void         SetTrackCuts(AliESDtrackCuts *c)               { fTrCuts        = c;         }
  void         SetTrainMode(Bool_t b)                         { fTrainMode     = b;         }
  void         SetTrigName(const char *n)                     { fTrigName      = n;         }
  void         SetUseQualFlag(Bool_t b)                       { fUseQualFlag   = b;         }
  void         SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }

 protected:
  virtual void CalcCaloTriggers();
  virtual void CalcClusterProps();
  virtual void CalcPrimTracks();
  virtual void CalcMcInfo();
  virtual void CalcTracks();
  virtual void ClusterAfterburner();
  virtual void FillCellHists();
  virtual void FillClusHists();
  virtual void FillNtuple();
  virtual void FillOtherHists();
  virtual void FillPionHists();
  virtual void FillMcHists();
  virtual void FillTrackHists();
  void         FillVertex(AliStaVertex *v, const AliESDVertex *esdv);
  void         FillVertex(AliStaVertex *v, const AliAODVertex *aodv);
  Double_t     GetCellIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2)                        const;
  Double_t     GetCellIsoNxM(Double_t cEta, Double_t cPhi, Int_t N, Int_t M)                              const;
  Double_t     GetCellEnergy(const AliVCluster *c)    const;
  Double_t     GetMaxCellEnergy(const AliVCluster *c) const { Short_t id=-1; return GetMaxCellEnergy(c,id); }
  Double_t     GetMaxCellEnergy(const AliVCluster *c, Short_t &id)                                        const;
  Double_t     GetSecondMaxCellEnergy(AliVCluster *clus, Short_t &id)                                     const;
  Int_t        GetNCells(const AliVCluster *c, Double_t emin=0.)                                          const;
  Int_t        GetNCells(Int_t sm, Double_t emin=0.)                                                      const;
  void         GetSigma(const AliVCluster *c, Double_t &sigmaMax, Double_t &sigmaMin)                     const;
  void         GetSigmaEtaEta(const AliVCluster *c, Double_t &sigmaEtaEta, Double_t &sigmaPhiPhi)         const;
  Double_t     GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2, Double_t pt=0.)       const;
  Double_t     GetTrackIsoStrip(Double_t cEta, Double_t cPhi, Double_t dEta=0.015, Double_t dPhi=0.3, Double_t pt=0.)       const;
  Bool_t       IsShared(const AliVCluster *c)                                                             const;
  Bool_t       IsIdPartOfCluster(const AliVCluster *c, Short_t id)                                        const;
  void         PrintDaughters(const AliVParticle *p, const TObjArray *arr, Int_t level=0)                 const;
  void         PrintDaughters(const AliMCParticle *p, const AliMCEvent *arr, Int_t level=0)               const;
  void         PrintTrackRefs(AliMCParticle *p)                                                           const;
  void         ProcessDaughters(AliVParticle *p, Int_t index, const TObjArray *arr);
  void         ProcessDaughters(AliMCParticle *p, Int_t index, const AliMCEvent *arr);

    // input members
  TString                fCentVar;                // variable for centrality determination
  Double_t               fCentFrom;               // min centrality (def=0)
  Double_t               fCentTo;                 // max centrality (def=100)
  Double_t               fVtxZMin;                // min primary vertex z (def=-10cm)
  Double_t               fVtxZMax;                // max primary vertex z (def=+10cm)
  Bool_t                 fUseQualFlag;            // if true use quality flag for centrality
  TString                fClusName;               // cluster branch name (def="")
  Bool_t                 fDoNtuple;               // if true write out ntuple
  Bool_t                 fDoAfterburner;          // if true run after burner
  Double_t               fAsymMax;                // maximum energy asymmetry (def=1)
  Int_t                  fNminCells;              // minimum number of cells attached to cluster (def=1)
  Double_t               fMinE;                   // minimum cluster energy (def=0.1 GeV/c)
  Double_t               fMinErat;                // minimum emax/ec ratio (def=0)
  Double_t               fMinEcc;                 // minimum eccentricity (def=0)
  TString                fGeoName;                // geometry name (def = EMCAL_FIRSTYEARV1)
  Double_t               fMinNClusPerTr;          // minimum number of cluster per track (def=50)
  Double_t               fIsoDist;                // isolation distance (def=0.2)
  TString                fTrClassNames;           // trigger class names
  AliESDtrackCuts       *fTrCuts;                 // track cuts
  AliESDtrackCuts       *fPrimTrCuts;             // track cuts
  TString                fPrimTracksName;         // name of track collection (if "" use branch)
  Bool_t                 fDoTrMatGeom;            // track matching including geometry
  Bool_t                 fTrainMode;              // train mode with minimal number of resources
  TString                fMarkCells;              // list of mark cells to monitor
  Int_t                  fMinL0Time;              // minimum accepted time for trigger
  Int_t                  fMaxL0Time;              // maximum accepted time for trigger
  Bool_t                 fMcMode;                 // monte carlo mode
  Bool_t                 fEmbedMode;              // embedding mode
  AliEMCALGeometry      *fGeom;                   // geometry utils
  AliEMCALRecoUtils     *fReco;                   // reco utils
  TString                fTrigName;               // trigger name
  Bool_t                 fDoPSel;                 // if false then accept all events
    // derived members (ie with ! after //)
  Bool_t                 fIsGeoMatsSet;           //!indicate that geo matrices are set 
  ULong64_t              fNEvs;                   //!accepted events 
  TList                 *fOutput;                 //!container of output histograms
  TObjArray             *fTrClassNamesArr;        //!array of trig class names  
  AliESDEvent           *fEsdEv;                  //!pointer to input esd event
  AliAODEvent           *fAodEv;                  //!pointer to input aod event
  const TObjArray       *fRecPoints;              //!pointer to rec points (AliAnalysisTaskEMCALClusterizeFast)
  const TClonesArray    *fDigits;                 //!pointer to digits     (AliAnalysisTaskEMCALClusterizeFast)
  TObjArray             *fEsdClusters;            //!pointer to esd clusters
  AliESDCaloCells       *fEsdCells;               //!pointer to esd cells
  TObjArray             *fAodClusters;            //!pointer to aod clusters
  AliAODCaloCells       *fAodCells;               //!pointer to aod cells
  TAxis                 *fPtRanges;               //!pointer to pt ranges
  TObjArray             *fSelTracks;              //!pointer to selected tracks
  TObjArray             *fSelPrimTracks;          //!pointer to selected primary tracks
    // ntuple
  TTree                 *fNtuple;                 //!pointer to ntuple
  AliStaHeader          *fHeader;                 //!pointer to header
  AliStaVertex          *fPrimVert;               //!pointer to primary vertex
  AliStaVertex          *fSpdVert;                //!pointer to SPD vertex
  AliStaVertex          *fTpcVert;                //!pointer to TPC vertex
  TClonesArray          *fClusters;               //!pointer to clusters
  TClonesArray          *fTriggers;               //!pointer to triggers
  TClonesArray          *fMcParts;                //!pointer to mc particles
    // histograms
  TH1                   *fHCuts;                  //!histo for cuts
  TH1                   *fHVertexZ;               //!histo for vtxz
  TH1                   *fHVertexZ2;              //!histo for vtxz after vtx cuts
  TH1                   *fHCent;                  //!histo for cent
  TH1                   *fHCentQual;              //!histo for cent after quality flag cut
  TH1                   *fHTclsBeforeCuts;        //!histo for trigger classes before cuts
  TH1                   *fHTclsAfterCuts;         //!histo for trigger classes after cuts

    // histograms for cells
  TH2                  **fHColuRow;               //!histo for cell column and row
  TH2                  **fHColuRowE;              //!histo for cell column and row weight energy
  TH1                  **fHCellMult;              //!histo for cell multiplicity in module
  TH1                   *fHCellE;                 //!histo for cell energy
  TH1                   *fHCellH;                 //!histo for highest cell energy
  TH1                   *fHCellM;                 //!histo for mean cell energy (normalized to hit cells)
  TH1                   *fHCellM2;                //!histo for mean cell energy (normalized to all cells)
  TH1                  **fHCellFreqNoCut;         //!histo for cell frequency without cut
  TH1                  **fHCellFreqCut100M;       //!histo for cell frequency with cut 100MeV
  TH1                  **fHCellFreqCut300M;       //!histo for cell frequency with cut 300MeV
  TH1                  **fHCellFreqE;             //!histo for cell frequency weighted with energy
  TH1                  **fHCellCheckE;            //!histo for cell E distribution for given channels
    // histograms for clusters
  TH1                   *fHClustEccentricity;     //!histo for cluster eccentricity
  TH2                   *fHClustEtaPhi;           //!histo for cluster eta vs. phi
  TH2                   *fHClustEnergyPt;         //!histo for cluster energy vs. pT
  TH2                   *fHClustEnergySigma;      //!histo for cluster energy vs. variance over long axis 
  TH2                   *fHClustSigmaSigma;       //!histo for sigma vs. lambda_0 comparison
  TH2                   *fHClustNCellEnergyRatio; //!histo for cluster n cells vs. energy ratio
  TH2			*fHClustEnergyNCell;      //!histo for cluster energy vs. cluster n cells
    // histograms for primary tracks
  TH1			*fHPrimTrackPt;           //!histo for primary track pt
  TH1			*fHPrimTrackEta;          //!histo for primary track eta
  TH1			*fHPrimTrackPhi;           //!histo for primary track phi
    // histograms for track matching
  TH1                   *fHMatchDr;               //!histo for dR track cluster matching
  TH1                   *fHMatchDz;               //!histo for dZ track cluster matching
  TH1                   *fHMatchEp;               //!histo for E/p track cluster matching
    // histograms for pion candidates
  TH2                   *fHPionEtaPhi;            //!histo for pion eta vs. phi
  TH2                   *fHPionMggPt;             //!histo for pion mass vs. pT
  TH2                   *fHPionMggAsym;           //!histo for pion mass vs. asym
  TH2                   *fHPionMggDgg;            //!histo for pion mass vs. opening angle
  TH1                   *fHPionInvMasses[21];     //!histos for invariant mass plots 
    // histograms for MC

 private:
  AliAnalysisTaskEMCALPi0PbPb(const AliAnalysisTaskEMCALPi0PbPb&);            // not implemented
  AliAnalysisTaskEMCALPi0PbPb &operator=(const AliAnalysisTaskEMCALPi0PbPb&); // not implemented

  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 13) // Analysis task for neutral pions in Pb+Pb
};
#endif
