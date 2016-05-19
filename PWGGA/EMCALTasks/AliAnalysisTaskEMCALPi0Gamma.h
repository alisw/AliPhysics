#ifndef AliAnalysisTaskEMCALPi0Gamma_h
#define AliAnalysisTaskEMCALPi0Gamma_h

// $Id:
//
// Analysis task for neutral pions (into two gammas), and for direct photons by subtraction method
//
// Author: B. Sahlmueller, based on code by C. Loizides

class TAxis;
class TClonesArray;
class TH1;
class TH2;
class TF1;
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
class AliGenHijingEventHeader;
class AliGenPythiaEventHeader;
class AliGenEventHeader;

#include "AliAnalysisTaskSE.h"
#include "AliStaObjects.h"


// class to store emc hits (for now, very simple)
class EmcHit {
  
  TLorentzVector thishit;
  Short_t hittype; // 100 for bg event/real event, 1 for added pi0, 2 for added eta, -1 for not needed; from primary pi0: 101, from secondary pi0: 102, from K0: 103, from material: 104

  Int_t imo; // index of original mother in monte carlo stack
  Int_t pid; // particle ID
  Double_t weight; // weight from mother particle
  Bool_t bclean; // clean if only one contributor
  
public:
  //virtual ~EmcHit();
  EmcHit();
  friend class EmcEvent;
  friend class AliAnalysisTaskEMCALPi0Gamma;
  
  void Print(){Printf("E=%.2f, type=%d, MoID=%d, PID=%d, w=%.3f",thishit.E(),hittype,imo,pid,weight);   }
};

// class to store old events
class EmcEvent {
  
  //    Int_t fCenPercent;
  //    Int_t fVtx;
	Float_t TrigPhi; // phi of highest pT hit on EMCal
  Float_t TrigTheta; // eta of highest pT hit ...
  
  const static int nMaxHit = 800;
  
  int nHits;
  EmcHit hit[nMaxHit];
  
public:
  EmcEvent();
  EmcEvent(const EmcEvent &obj);
  //virtual ~EmcEvent();
  //    void SetGlobalInfo(const Int_t&, const Int_t&, const Int_t&, const Int_t&, const Double_t&, const Double_t&);
  void SetGlobalInfo(const Int_t&, const Float_t&, const Float_t&);
  int evsize() {return nHits;}
  void Reset();
  void Print();
  
  friend class AliAnalysisTaskEMCALPi0Gamma;
};


class AliAnalysisTaskEMCALPi0Gamma : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALPi0Gamma();
  AliAnalysisTaskEMCALPi0Gamma(const char *name);
  virtual ~AliAnalysisTaskEMCALPi0Gamma();
  
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);
  
  void         SetAsymMax1(Double_t asymMax)                  { fAsymMax1      = asymMax;   }
  void         SetAsymMax2(Double_t asymMax)                  { fAsymMax2      = asymMax;   }
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
  void         SetM02Cut(Double_t m02)                        { fM02           = m02;         }
  void         SetMinEcc(Double_t ecc)                        { fMinEcc        = ecc;       }
  void         SetMinErat(Double_t erat)                      { fMinErat       = erat;      }
  void         SetMinNClustersPerTrack(Double_t m)            { fMinNClusPerTr = m;         }
  void         SetNminCells(Int_t n)                          { fNminCells     = n;         }
  void         SetTrackMatchSimple(Bool_t b)                  { fDoTrMtSmpl    = b;         }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)           { fPrimTrCuts    = c;         }
  void         SetPrimTracksName(const char *n)               { fPrimTracksName = n;        }
  void         SetRecoUtils(AliEMCALRecoUtils *reco)          { fReco          = reco;      }
  void         SetTrClassNames(const char *n)                 { fTrClassNames  = n;         }
  void         SetTrackCuts(AliESDtrackCuts *c)               { fTrCuts        = c;         }
  void         SetTrainMode(Bool_t b)                         { fTrainMode     = b;         }
  void         SetTrigName(const char *n)                     { fTrigName      = n;         }
  void         SetUseQualFlag(Bool_t b)                       { fUseQualFlag   = b;         }
  void         SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }
  void         SetRotateMixed(Bool_t b)                       { fRotateMixed   = b;         }
  void         SetAddedSignal(Bool_t b)                       { fAddedSignal   = b;         }
  void         SetDataPeriod(Int_t b)                         { fDataPeriod   = b;         }
  void         SetDoManualRecal(Bool_t b)                     { fDoManualRecal = b;         }
  void         SetDoCalibRun(Bool_t b)                     { fCalibRun = b;         }
  void         SetDnDpT(Int_t i, Double_t par0, Double_t par1, Double_t par2, Double_t par3, Double_t par4);
  
protected:
  
  virtual void CalcMcInfo();
  virtual void ClusterAfterburner();
  void AddMixEvent(const Int_t, const Int_t, const Int_t, Int_t&, const Float_t&, const Float_t&);
  Double_t FillClusHists(Float_t&, Float_t&);
  void FillMixHists(const Int_t, const Int_t, const Int_t, const Double_t, const Double_t);
  virtual void FillNtuple();
  virtual void FillOtherHists();
  virtual void FillPionHists();
  virtual void FillMcHists();
  virtual void FillTrackHists();
  void         FillVertex(AliStaVertex *v, const AliESDVertex *esdv);
  void         FillVertex(AliStaVertex *v, const AliAODVertex *aodv);
  void GetMulClass(Int_t&);
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
  Int_t        GetModuleNumber(AliVCluster * cluster)                                                     const;
  void         FillCellQAHists(AliVCluster *);
  Double_t     PrivateEnergyRecal(Double_t energy, Int_t iCalib);
  // spectral shape
  Double_t CalcWeight(Double_t pt,Double_t eta, Int_t i);

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
  Double_t               fAsymMax1;                // maximum energy asymmetry (def=1)
  Double_t               fAsymMax2;                // maximum energy asymmetry (def=1)
  Double_t               fAsymMax3;                // maximum energy asymmetry (def=1)
  Int_t                  fNminCells;              // minimum number of cells attached to cluster (def=1)
  Double_t               fMinE;                   // minimum cluster energy (def=0.1 GeV/c)
  Double_t               fM02;                    // maximum M02
  Double_t               fMinErat;                // minimum emax/ec ratio (def=0)
  Double_t               fMinEcc;                 // minimum eccentricity (def=0)
  Bool_t                 fCalibRun;               // fill energy calibration histograms (def = 0)
  Bool_t                 fDoManualRecal;          // do manual recalibration here? (def = 0)
  Bool_t                 fDoTrMtSmpl;             // use built in track matching? (def=0)
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
//  Bool_t                 fAddedSignal;            // added signals in MC?
  Bool_t                 fIsGeoMatsSet;           // indicate that geo matrices are set
  Bool_t                 fRotateMixed;            // rotates the events before mixing such that the highest pT cluster has same phi as in real event
  Bool_t                 fAddedSignal;            // added signals in MC?
  Int_t                  fDataPeriod;             // which period(s)

  // derived members (ie with ! after //)
  ULong64_t              fNEvs;                   //!accepted events
  TList                 *fOutput;                 //!container of output histograms
  TObjArray             *fTrClassNamesArr;        //!array of trig class names
  AliESDEvent           *fEsdEv;                  //!pointer to input esd event
  AliAODEvent           *fAodEv;                  //!pointer to input aod event
  //AliGenPythiaEventHeader* PythiaGenHeader;     //!pointer to Pythia event header Philipp
  //AliGenHijingEventHeader* hijingGenHeader;     //!pointer to Hijing event header Philipp
  AliGenEventHeader * eventHeader;                //!pointer to Generated event header Evi
  AliGenEventHeader * pythiaHeader;                //!pointer to Generated event header Evi
  AliGenEventHeader * addedPi0Header;               //!pointer to Added event header Evi
  AliGenEventHeader * addedEtaHeader;               //!pointer to Added event header Evi
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
  
  // QA profile histograms
  TProfile              *fHMeanClusterEnergy;    //!profile for cluster energy vs. runnumber
  TProfile              *fHMeanClusterNumber;    //!profile for number of clusters in each event
  
  // histogram for cells
  TH2                   *fHCellIndexEnergy;      //!histo for cell energy vs cell number
  
  // histograms for clusters
  TH1                   *fHClusters;                  //!histo for cuts
  TH1                   *fHClustAllEtaPhi;        //!histo for all clusters eta and phi
  TH1                   *fHClustNoEvt;            //!histo for number of clusters in event
  TH1                   *fHClustAccEvt;            //!histo for number of clusters after cuts in event
  TH1                   *fHClustEccentricity;     //!histo for cluster eccentricity
  TH2                   *fHClustEtaPhi;           //!histo for cluster eta vs. phi
  TH2                   *fHClustEnergyPt;         //!histo for cluster energy vs. pT
  TH2                   *fHClustEnergyPtDCal;         //!histo for cluster energy vs. pT
  TH2                   *fHClustEnergySM;         //!histo for cluster energy vs. Supermodule
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
  TH2                   *fHAddPionEtaPt;            //!histo for pion eta vs. pt
  TH2                   *fHAddPionEtaPtWgt;            //!histo for pion eta vs. pt
  TH2                   *fHPyPionEtaPt;            //!histo for pion eta vs. pt
  TH1                   *fHdr;                    //!histo for dR of pairs
  TH2                   *fHPionMggPt;             //!histo for pion mass vs. pT
  TH2                   *fHPionMggAsym;           //!histo for pion mass vs. asym
  TH2                   *fHPionPtAsym;           //!histo for pion pT vs. asym
  TH2                   *fHPionMggDgg;            //!histo for pion mass vs. opening angle
  TH2                   *fHPionInvMasses;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesSym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAsym;     //!histos for invariant mass plots
  TH2                   *fHPriPionInvMasses;     //!histos for primary pion invariant mass
  TH2                   *fHSecPionInvMasses;     //!histos for secondary pion invariant mass
  TH2                   *fHK0PionInvMasses;     //!histos for K0 pion invariant mass
  TH2                   *fHMatPionInvMasses;     //!histos for material pion invariant mass
  TH2                   *fHPionInvMassesAdd1NoWgt; //!histos for added signal pions without weight
  TH2                   *fHPionInvMassesAdd1;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd1Sym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd1Asym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd1Mult;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd1MultSym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd1MultAsym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1Sym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1Asym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1Mult;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1MultSym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesGamAdd1MultAsym;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesAdd2;     //!histos for invariant mass plots
  // mixing
  TH2                   *fHPionInvMassesMix;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMix1;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMix2;     //!histos for invariant mass plots
  // DCal
  TH2                   *fHPionInvMassesDCal;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMixDCal;     //!histos for invariant mass plots

  // DCal+EMCal
  TH2                   *fHPionInvMassesEMCalDCal;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMixEMCalDCal;     //!histos for invariant mass plots

  // calibration
  TH2                   *fHPionInvMassesEMCalCalib;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMixEMCalCalib;     //!histos for invariant mass plots

  // calibration
  TH2                   *fHPionInvMassesDCalCalib;     //!histos for invariant mass plots
  TH2                   *fHPionInvMassesMixDCalCalib;     //!histos for invariant mass plots

  
  // quick histo for J/Psi
  //TH2                   *fHJPInvMasses;            //!histo for inv mass in JPsi region
  // primary pions
  TH2                   *fHPrimPionInvMasses;       //!histos for invariant mass plots
  TH2                   *fHPrimPionInvMassesAsym;       //!histos for invariant mass plots
  
  // conversion info
  TH1                   * fHConversionPoint;   //!histo for conversion position in XY
  
  // histograms for MC
  TH1                   *fHWgt;         //!histo for weight of particles
  TH1                   *fHPionTruthPt;       //!histo for pT from MC pion
  TH1                   *fHPionTruthPtIn;    //!histo for pT for MC pion in eta range
  TH1                   *fHPionTruthPtAcc;    //!histo for pT for MC pion in acceptance
  TH1                   *fHEtaTruthPt;       //!histo for pT from MC eta
  TH1                   *fHEtaTruthPtIn;    //!histo for pT for MC eta in eta range
  TH1                   *fHEtaTruthPtAcc;    //!histo for pT for MC eta in acceptance
  TH1                   *fHGamTruthPt;       //!histo for pT from MC gamma
  TH1                   *fHGamTruthPtIn;    //!histo for pT for MC gamma in eta range
  TH1                   *fHGamTruthPtAcc;    //!histo for pT for MC gamma in acceptance
  
  TH1                   *fHPionTruthPtAdd;       //!histo for pT from MC pion
  TH1                   *fHPionTruthPtInAdd;    //!histo for pT for MC pion in eta range
  TH1                   *fHPionTruthPtAccAdd;    //!histo for pT for MC pion in acceptance
  TH1                   *fHEtaTruthPtAdd;       //!histo for pT from MC eta
  TH1                   *fHEtaTruthPtInAdd;    //!histo for pT for MC eta in eta range
  TH1                   *fHEtaTruthPtAccAdd;    //!histo for pT for MC eta in acceptance
  //TH1                   *fHGamTruthPtAdd;       //!histo for pT from MC gamma
  //TH1                   *fHGamTruthPtInAdd;    //!histo for pT for MC gamma in eta range
  //TH1                   *fHGamTruthPtAccAdd;    //!histo for pT for MC gamma in acceptance

  // for MC particle stuff
  TH2                   *fHMCpartfrac; //!histo for fraction of MC of cluster and for fraction of highest MC to all MC contributors
  TH2                   *fHECluEMC;    //!histo for MC energy vs. cluster energy
  TH2                   *fHECluEMCAddPi0;    //!histo for MC energy vs. cluster energy
  TH2                   *fHECluEMCAddEta;    //!histo for MC energy vs. cluster energy
//  TH2                   *fHRecTrue;   //!histo for MC/Clu energy vs. cluster energy
//  TH2                   *fHRecTrueAddPi0;   //!histo for MC/Clu energy vs. cluster energy
//  TH2                   *fHRecTrueAddEta;   //!histo for MC/Clu energy vs. cluster energy
  TH2                   *fHECluEMCnofull;    //!histo for MC energy vs. cluster energy
  TH2                   *fHECluEMCnofullAdd;    //!histo for MC energy vs. cluster energy
  TH2                   *fHECluEMCelectron;  //!histo for MC energy vs. cluster energy for electrons
  TH2                   *fHECluEMCpion;  //!histo for MC energy vs. cluster energy for pions
  TH2                   *fHECluEMCkaon;  //!histo for MC energy vs. cluster energy for kaons
  TH2                   *fHECluEMCother;  //!histo for MC energy vs. cluster energy for others
  TH2                   *fHECluEMCpi0single; //!histo for MC pt vs. cluster pt for merged pi0
  TH2                   *fHNMothers; //!histo to store number of mother iterations to primary
  
  
  // more histograms
//  TH1		    		  *fHMixRotation;        //! histo to show how much the mixed events were rotated in phi
//  TH1             *fHCorrection;         //! histo to show the correction factor
//  TH2             *fHPionSm;             //! histo to see the change of the pair mass due to the correction
  
  // store start and end of pythia particles, added signals, ...
  Int_t ipymin; // first pythia particle index
  Int_t ipymax; // last pythia particle index
  Int_t ipi0min; // first added pi0 particle index
  Int_t ipi0max; // last added pi0 particle index
  Int_t ietamin; // first added eta particle index
  Int_t ietamax; // last added eta particle index
  
  const static int nMulClass =   5;
  const static int nZClass   =   3;
  const static int nPtClass = 1;
  int iEvt[nMulClass][nZClass][nPtClass];
  const static int nEvt      =   30; // mixing "depth"
  
  EmcEvent evt;
  EmcEvent EmcEventList[nMulClass][nZClass][nPtClass][nEvt];
  
  EmcEvent thisEvent;



private:
  AliAnalysisTaskEMCALPi0Gamma(const AliAnalysisTaskEMCALPi0Gamma&);            // not implemented
  AliAnalysisTaskEMCALPi0Gamma &operator=(const AliAnalysisTaskEMCALPi0Gamma&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALPi0Gamma, 15) // Analysis task for neutral pions in Pb+Pb
};


#endif


