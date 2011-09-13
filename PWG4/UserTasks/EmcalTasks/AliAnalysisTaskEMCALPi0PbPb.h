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
  void         SetDoTrackMatWithGeom(Bool_t b)                { fDoTrMatGeom   = b;         }
  void         SetFillNtuple(Bool_t b)                        { fDoNtuple      = b;         }
  void         SetGeoName(const char *n)                      { fGeoName       = n;         }
  void         SetGeoUtils(AliEMCALGeometry *geo)             { fGeom          = geo;       }
  void         SetIsoDist(Double_t d)                         { fIsoDist       = d;         }
  void         SetL0TimeRange(Int_t l, Int_t h)               { fMinL0Time=l; fMaxL0Time=h; }
  void         SetMarkCells(const char *n)                    { fMarkCells     = n;         }
  void         SetMcMode(Bool_t b)                            { fMcMode        = b;         }
  void         SetEmbedMode(Bool_t b)                         { fEmbedMode     = b;         }
  void         SetMinClusEnergy(Double_t e)                   { fMinE          = e;         }
  void         SetMinEcc(Double_t ecc)                        { fMinEcc        = ecc;       }
  void         SetMinErat(Double_t erat)                      { fMinErat       = erat;      }
  void         SetMinNClustersPerTrack(Double_t m)            { fMinNClusPerTr = m;         }
  void         SetNminCells(Int_t n)                          { fNminCells     = n;         }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)           { fPrimTrCuts    = c;         }
  void         SetRecoUtils(AliEMCALRecoUtils *reco)          { fReco          = reco;      }
  void         SetTrClassNames(const char *n)                 { fTrClassNames  = n;         }
  void         SetTrackCuts(AliESDtrackCuts *c)               { fTrCuts        = c;         }
  void         SetTrainMode(Bool_t b)                         { fTrainMode     = b;         }
  void         SetUseQualFlag(Bool_t b)                       { fUseQualFlag   = b;         }
  void         SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }
  void         SetDoPhysicsSelection(Bool_t b)                { fDoPSel        = b;         }

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
  Int_t        GetNCells(const AliVCluster *c, Double_t emin=0.)                                          const;
  void         GetSigma(const AliVCluster *c, Double_t &sigmaMax, Double_t &sigmaMin)                     const;
  void         GetSigmaEtaEta(const AliVCluster *c, Double_t &sigmaEtaEta, Double_t &sigmaPhiPhi)         const;
  Double_t     GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2, Double_t pt=0.)       const;
  Double_t     GetTrackIsoStrip(Double_t cEta, Double_t cPhi, Double_t dEta=0.015, Double_t dPhi=0.3, Double_t pt=0.)       const;
  Double_t     GetTrigEnergy(const AliVCluster *c)                                                        const;
  Bool_t       IsShared(const AliVCluster *c)                                                             const;
  void         PrintDaughters(const AliVParticle *p, const TObjArray *arr, Int_t level=0)                 const;
  void         PrintDaughters(const AliMCParticle *p, const AliMCEvent *arr, Int_t level=0)               const;
  void         PrintTrackRefs(AliMCParticle *p)                                                           const;
  void         ProcessDaughters(AliVParticle *p, Int_t index, const TObjArray *arr);
  void         ProcessDaughters(AliMCParticle *p, Int_t index, const AliMCEvent *arr);
  Double_t     GetSecondMaxCell(AliVCluster *clus);

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
  Bool_t                 fDoTrMatGeom;            // track matching including geometry
  Bool_t                 fTrainMode;              // train mode with minimal number of resources
  TString                fMarkCells;              // list of mark cells to monitor
  Int_t                  fMinL0Time;              // minimum accepted time for trigger
  Int_t                  fMaxL0Time;              // maximum accepted time for trigger
  Bool_t                 fMcMode;                 // monte carlo mode
  Bool_t                 fEmbedMode;              // embedding mode
  AliEMCALGeometry      *fGeom;                   // geometry utils
  AliEMCALRecoUtils     *fReco;                   // reco utils
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
  Int_t                  fNAmpInTrigger;          //!number of cells to keep trigger statistic
  Float_t               *fAmpInTrigger;           //!amplitude for calo cells which are part of trigger
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

  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 11) // Analysis task for neutral pions in Pb+Pb
};
#endif

#ifndef AliStaObjs_h
#define AliStaObjs_h
class AliStaHeader
{
 public:
  AliStaHeader() : fRun(0), fOrbit(0), fPeriod(0), fBx(0), fL0(0), fL1(0), fL2(0),
                   fTrClassMask(0), fTrCluster(0), fOffTriggers(0), fFiredTriggers(),
                   fTcls(0), fV0And(0), fV0Cent(0), fV0(0), fCl1Cent(0), fCl1(0), fTrCent(0), fTr(0),
                   fCqual(-1), fPsi(0), fPsiRes(0), fNSelTr(0), fNSelPrimTr(0), fNSelPrimTr1(0),
                   fNSelPrimTr2(0), fNCells(0), fNCells0(0), fNCells01(0), fNCells03(0), 
                   fNCells1(0), fNCells2(0), fNCells5(0), fNClus(0), fNClus1(0), fNClus2(0), fNClus5(0), 
                   fMaxCellE(0), fMaxClusE(0) {;}

  ULong64_t     GetEventId() const {
                  return (((ULong64_t)fPeriod << 36) |
                          ((ULong64_t)fOrbit  << 12) |
                          (ULong64_t)fBx); 
                }
  virtual ~AliStaHeader() {;}

 public:
  Int_t         fRun;            //         run number
  UInt_t        fOrbit;          //         orbit number
  UInt_t        fPeriod;         //         period number
  UShort_t      fBx;             //         bunch crossing id
  UInt_t        fL0;             //         l0 trigger bits
  UInt_t        fL1;             //         l1 trigger bits
  UShort_t      fL2;             //         l2 trigger bits
  ULong64_t     fTrClassMask;    //         trigger class mask
  UChar_t       fTrCluster;      //         trigger cluster mask
  UInt_t        fOffTriggers;    //         fired offline triggers for this event
  TString       fFiredTriggers;  //         string with fired triggers
  UInt_t        fTcls;           //         custom trigger definition
  Bool_t        fV0And;          //         V0AND (from AliTriggerAnalysis)
  Double32_t    fV0Cent;         //[0,0,16] v0 cent
  Double32_t    fV0;             //[0,0,16] v0 result used for cent 
  Double32_t    fCl1Cent;        //[0,0,16] cl1 cent
  Double32_t    fCl1;            //[0,0,16] cl1 result used for cent 
  Double32_t    fTrCent;         //[0,0,16] tr cent
  Double32_t    fTr;             //[0,0,16] tr result used for cent 
  Int_t         fCqual;          //         centrality quality
  Double32_t    fPsi;            //[0,0,16] event-plane angle
  Double32_t    fPsiRes;         //[0,0,16] event-plane ange resolution
  UShort_t      fNSelTr;         //         # selected tracks         
  UShort_t      fNSelPrimTr;     //         # selected tracks (primary)
  UShort_t      fNSelPrimTr1;    //         # selected tracks (primary) pt > 1 GeV/c
  UShort_t      fNSelPrimTr2;    //         # selected tracks (primary) pt > 2 GeV/c
  UShort_t      fNCells;         //         # cells
  UShort_t      fNCells0;        //         # cells > 0.45 GeV
  UShort_t      fNCells01;       //         # cells > 0.1  GeV
  UShort_t      fNCells03;       //         # cells > 0.3  GeV
  UShort_t      fNCells1;        //         # cells > 1    GeV
  UShort_t      fNCells2;        //         # cells > 2    GeV
  UShort_t      fNCells5;        //         # cells > 5    GeV
  UShort_t      fNClus;          //         # clus
  UShort_t      fNClus1;         //         # clus > 1 GeV
  UShort_t      fNClus2;         //         # clus > 2 GeV
  UShort_t      fNClus5;         //         # clus > 5 GeV
  Double32_t    fMaxCellE;       //[0,0,16] maximum cell energy
  Double32_t    fMaxClusE;       //[0,0,16] maximum clus energy

  ClassDef(AliStaHeader,4) // Header class
};

class AliStaVertex
{
 public:
  AliStaVertex(Double_t x=0, Double_t y=0, Double_t z=0) : fVx(x), fVy(y), fVz(z), fVc(-1), fDisp(0), fZres(0),
                                                           fChi2(0), fSt(0), fIs3D(0), fIsZ(0) {;}
  virtual ~AliStaVertex() {;}

 public:
  Double_t      fVx;          //[0,0,16] vertex x
  Double_t      fVy;          //[0,0,16] vertex y
  Double_t      fVz;          //[0,0,16] vertex z
  Double_t      fVc;          //[0,0,16] number of contributors to vertex
  Double_t      fDisp;        //[0,0,16] dispersion
  Double_t      fZres;        //[0,0,16] z-resolution
  Double_t      fChi2;        //[0,0,16] chi2 of fit
  Bool_t        fSt;          //         status bit
  Bool_t        fIs3D;        //         is vertex from 3D
  Bool_t        fIsZ;         //         is vertex from Z only

  ClassDef(AliStaVertex,1) // Vertex class
};

class AliStaCluster : public TObject
{
 public:
    AliStaCluster() : TObject(), fE(0), fR(0), fEta(0), fPhi(0), fN(0), fN1(0), fN3(0), fIdMax(0), fEmax(0), fTmax(0), 
                      fE2max(0),fDbc(-1), fDisp(-1), fM20(0), fM02(0), fEcc(0), fSig(0), fSigEtaEta(0), fSigPhiPhi(0),
                      fIsTrackM(0), fTrDz(0), fTrDr(-1), fTrEp(0), fTrDedx(0), fTrIso(0), fTrIso1(0), fTrIso2(0),  
                      fTrIsoD1(0), fTrIso1D1(0), fTrIso2D1(0), fTrIsoD3(0), fTrIso1D3(0), fTrIso2D3(0),fTrIsoStrip(0),
                      fCeIso(0), fCeIso1(0), fCeIso3(0), fCeIso4x4(0), fCeIso5x5(0), fCeCore(0), fCeIso3x22(0), 
                      fIsTrigM(0), fTrigE(-1), fTrigMaskE(-1), fIsShared(0), fMcLabel(-1), fEmbE(0) {;}

 public:
  Double32_t    fE;                //[0,0,16] energy
  Double32_t    fR;                //[0,0,16] radius (cylinder)
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  UChar_t       fN;                //         number of cells
  UChar_t       fN1;               //         number of cells > 100 MeV
  UChar_t       fN3;               //         number of cells > 300 MeV
  UShort_t      fIdMax;            //         id maximum cell
  Double32_t    fEmax;             //[0,0,16] energy of maximum cell
  Double32_t    fTmax;             //[0,0,16] time of maximum cell
  Double32_t    fE2max;            //[0,0,16] energy of second maximum cell
  Double32_t    fDbc;              //[0,0,16] distance to nearest bad channel
  Double32_t    fDisp;             //[0,0,16] cluster dispersion, for shape analysis
  Double32_t    fM20;              //[0,0,16] 2-nd moment along the main eigen axis
  Double32_t    fM02;              //[0,0,16] 2-nd moment along the second eigen axis
  Double32_t    fEcc;              //[0,0,16] eccentricity
  Double32_t    fSig;              //[0,0,16] sigma
  Double32_t    fSigEtaEta;        //[0,0,16] sigma eta-eta
  Double32_t    fSigPhiPhi;        //[0,0,16] sigma phi-phi
  Bool_t        fIsTrackM;         //         if true then track values are set
  Double32_t    fTrDz;             //[0,0,16] dZ to nearest track
  Double32_t    fTrDr;             //[0,0,16] dR to nearest track (in x,y)
  Double32_t    fTrEp;             //[0,0,16] E/P to nearest track 
  Double32_t    fTrDedx;           //[0,0,16] dE/dx (TPC signal) to nearest track 
  Double32_t    fTrIso;            //[0,0,16] track isolation
  Double32_t    fTrIso1;           //[0,0,16] track isolation (pt>1GeV/c)
  Double32_t    fTrIso2;           //[0,0,16] track isolation (pt>2GeV/c)
  Double32_t    fTrIsoD1;          //[0,0,16] track isolation, iso dist 0.25
  Double32_t    fTrIso1D1;         //[0,0,16] track isolation (pt>1GeV/c), iso dist 0.1
  Double32_t    fTrIso2D1;         //[0,0,16] track isolation (pt>2GeV/c), iso dist 0.1
  Double32_t    fTrIsoD3;          //[0,0,16] track isolation, iso dist 0.3
  Double32_t    fTrIso1D3;         //[0,0,16] track isolation (pt>1GeV/c), iso dist 0.3
  Double32_t    fTrIso2D3;         //[0,0,16] track isolation (pt>2GeV/c), iso dist 0.3
  Double32_t    fTrIsoStrip;       //[0,0,16] track isolation strip, dEtaXdPhi=0.015x0.3
  Double32_t    fCeIso;            //[0,0,16] cell isolation in R=0.20
  Double32_t    fCeIso1;           //[0,0,16] cell isolation in R=0.10
  Double32_t    fCeIso3;          //[0,0,16] cell isolation in  R=0.30
  Double32_t    fCeIso4x4;         //[0,0,16] cell isolation in  4x4 cells
  Double32_t    fCeIso5x5;         //[0,0,16] cell isolation in  5x5 cells
  Double32_t    fCeCore;           //[0,0,16] cell content in R=0.05 
  Double32_t    fCeIso3x22;        //[0,0,16] cell isolation in rectangular strip of dEtaXdPhi=0.042x0.308
  Bool_t        fIsTrigM;          //         if true then trigger values are set
  Double32_t    fTrigE;            //[0,0,16] trigger tower energy
  Double32_t    fTrigMaskE;        //[0,0,16] masked trigger tower energy
  Bool_t        fIsShared;         //         =true then extends across more than one super module
  Short_t       fMcLabel;          //         index of closest MC particle
  Double32_t    fEmbE;             //[0,0,16] sum of energy of embedded (MC) cells in cluster

  ClassDef(AliStaCluster,6) // Cluster class
};

class AliStaTrigger : public TObject
{
 public:
  AliStaTrigger() : TObject(), fE(0), fEta(0), fPhi(0), fAmp(0), fMinTime(0), fMaxTime(0) {}

 public:
  Double32_t    fE;                //[0,0,16] energy
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fAmp;              //[0,0,16] amplitude
  Short_t       fMinTime;          //        minimum L0 "time"
  Short_t       fMaxTime;          //        maximum L0 "time"

  ClassDef(AliStaTrigger,1) // Trigger class
};

class AliStaPart : public TObject
{
 public:
  AliStaPart() : TObject(), fPt(0), fEta(0), fPhi(0), fVR(0), fVEta(0), fVPhi(0), fPid(0), fMo(-1), fDet(-2), 
                 fLab(-1), fNs(0) { memset(fDs,-1,sizeof(Short_t)*99); }

  Int_t         OnEmcal() const { return (fDet==8);  }
  Int_t         IsSim()   const { return (fDet!=-2); }
    
 public:
  Double32_t    fPt;               //[0,0,16] pt
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fVR;               //[0,0,16] prod r (cylinder)
  Double32_t    fVEta;             //[0,0,16] prod eta
  Double32_t    fVPhi;             //[0,0,16] prod phi
  Short_t       fPid;              //         pid
  Short_t       fMo;               //         index of mother
  Short_t       fDet;              //         detector in which particle left trace (8 for EMCAL, see AliTrackReference.h)
    // the following must be filled before first usage
  Short_t       fLab;              //!        label (index in array)
  Short_t       fNs;               //!        number of daughters
  Short_t       fDs[99];           //!        daughters

  ClassDef(AliStaPart,1) // Particle class
};
#endif
