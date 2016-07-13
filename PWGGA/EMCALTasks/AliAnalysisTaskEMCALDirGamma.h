#ifndef AliAnalysisTaskEMCALDirGamma_h
#define AliAnalysisTaskEMCALDirGamma_h

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
//class AliV0ReaderV1;
class AliGenHijingEventHeader;
class AliGenPythiaEventHeader;
class AliGenEventHeader;
class AliPIDResponse;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#include "AliStaObjects.h"

class AliESDv0KineCuts;

class AliAnalysisTaskEMCALDirGamma : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALDirGamma();
  AliAnalysisTaskEMCALDirGamma(const char *name);
  virtual ~AliAnalysisTaskEMCALDirGamma();
  
  
  
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
  void         SetfDoConvAna(Bool_t b)                        { fDoConvAna     = b;         }
  void         SetGeoName(const char *n)                      { fGeoName       = n;         }
  void         SetGeoUtils(AliEMCALGeometry *geo)             { fGeom          = geo;       }
  void         SetIsoDist(Double_t d)                         { fIsoDist       = d;         }
  void         SetL0TimeRange(Int_t l, Int_t h)               { fMinL0Time=l; fMaxL0Time=h; }
  void         SetMarkCells(const char *n)                    { fMarkCells     = n;         }
  void         SetMcMode(Bool_t b)                            { fMcMode        = b;         }
  void         SetMinClusEnergy(Double_t e)                   { fMinE          = e;         }
 // void         SetM02Cut(Double_t m02)                        { fM02           = m02;         }
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
  void         SetDoManualRecal(Bool_t b)                     { fDoManualRecal = b;         }
  void         SetDnDpT(Int_t i, Double_t par0, Double_t par1, Double_t par2, Double_t par3, Double_t par4);
  
  
  
  
protected:
  
  virtual void CalcMcInfo();
  virtual void ClusterAfterburner();
  Double_t FillClusHists(Float_t&, Float_t&);
  virtual void FillNtuple();
  virtual void FillOtherHists();
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
  Bool_t                 fDoConvAna;              // if true do conversion photon analysis
  Bool_t                 fDoAfterburner;          // if true run after burner
  Double_t               fAsymMax1;                // maximum energy asymmetry (def=1)
  Double_t               fAsymMax2;                // maximum energy asymmetry (def=1)
  Double_t               fAsymMax3;                // maximum energy asymmetry (def=1)
  Int_t                  fNminCells;              // minimum number of cells attached to cluster (def=1)
  Double_t               fMinE;                   // minimum cluster energy (def=0.1 GeV/c)
  //Double_t               fM02;                    // maximum M02
  Double_t               fM02min;                    // minimum M02
  Double_t               fM02max;                    // maximum M02
  Double_t               fMinErat;                // minimum emax/ec ratio (def=0)
  Double_t               fMinEcc;                 // minimum eccentricity (def=0)
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
  

  // derived members (ie with ! after //)
  //AliV0ReaderV1 	*fV0Reader;               //!v0 reader
  TClonesArray 		*fReaderGammas;           //!for conversion gammas
  ULong64_t              fNEvs;                   //!accepted events
  TList                 *fOutput;                 //!container of output histograms
  TObjArray             *fTrClassNamesArr;        //!array of trig class names
  AliESDEvent           *fEsdEv;                  //!pointer to input esd event
  AliAODEvent           *fAodEv;                  //!pointer to input aod event
  //AliESDEvent           *fesdTOF;                  //!pointer to input esd event  
  AliPIDResponse *fPIDResponse;     //! PID response object
  AliESDtrackCuts *fESDtrackCuts; //! Track cuts TPC
  
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
  TTree					*ftrcuts;					//! tree with cut data
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
  
  // histograms for clusters
  TH1                   *fHClustNoEvt;            //!histo for number of clusters in event
  TH1                   *fHClustAccEvt;            //!histo for number of clusters after cuts in event
  TH1                   *fHClustEccentricity;     //!histo for cluster eccentricity
  TH2                   *fHClustEtaPhi;           //!histo for cluster eta vs. phi  
  TH2                   *fHEMCalModule0;           //!histo for cluster in module 0
  
  TH2                   *fHClustEtaPhiRaw;           //!histo for cluster eta vs. phi    
  TH2                   *fHv0TrackEtaPhi;           //!histo for cluster eta vs. phi  
  TH2                   *fHv0ClustEtaPhi;           //!histo for cluster eta vs. phi  
  TH2                   *fHv0TrackEtaPhi2;           //!histo for cluster eta vs. phi  
  TH2                   *fHv0ClustEtaPhi2;           //!histo for cluster eta vs. phi  
     
  TH2                   *fHClustEnergyPt;         //!histo for cluster energy vs. pT
  TH2                   *fHClustEnergySM;         //!histo for cluster energy vs. Supermodule
  TH2                   *fHClustEnergySigma;      //!histo for cluster energy vs. variance over long axis
  TH2                   *fHClustSigmaSigma;       //!histo for sigma vs. lambda_0 comparison
  TH2                   *fHClustEtaM02;				//!histo for eta vs m02
  TH2                   *fHClustPhiM02;				//!histo for phi vs m02  
  TH1                   *fHv0TrackPtEMCal;		//!histo for track pT of electrons
  TH1                   *fHv0TrackPt;		//!histo for track pT of electrons  
  TH1                   *fHClustETrackP; //!histo for ratio of cluster energy and track p
  
  
  TH2                   *fHClustNCellEnergyRatio; //!histo for cluster n cells vs. energy ratio
  TH2					*fHClustEnergyNCell;      //!histo for cluster energy vs. cluster n cells
  TH2					*fHClustEnergyNCellPion;      //!histo for cluster energy vs. cluster n cells  
  TH2					*fHClustEnergyNCellPhoton;      //!histo for cluster energy vs. cluster n cells  
  TH2					*fHClustNCellM02Photon; 	//! histo for cluster ncell vs. m02
  TH2					*fHClustNCellM02Pion; 	//! histo for cluster ncell vs. m02
                  		
  TH2					*fHClustEnergyNCellRaw;      //!histo for cluster energy vs. cluster n cells
  TH2					*fHClustEnergyNCellPionRaw;      //!histo for cluster energy vs. cluster n cells  
  TH2					*fHClustEnergyNCellPhotonRaw;      //!histo for cluster energy vs. cluster n cells   
  TH2					*fHClustNCellM02PhotonRaw; 	//! histo for cluster ncell vs. m02
  TH2					*fHClustNCellM02PionRaw; 	//! histo for cluster ncell vs. m02  
  
  
  // histograms for conversion photons
  TH2                 *fHConvEnergyPt;
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
  
  // conversion
  TH2                   *fHPionEtaPhiConv;            //!histo for pion eta vs. phi
  TH2                   *fHPionMggPtConv;             //!histo for pion mass vs. pT
  TH2                   *fHPionMggAsymConv;           //!histo for pion mass vs. asym
  TH2                   *fHPionMggDggConv;            //!histo for pion mass vs. opening angle
  TH1                   *fHPionInvMassesConv;     //!histos for invariant mass plots
  // mixing
  TH2                   *fHPionInvMassesConvMix;     //!histos for invariant mass plots
  // conversion conversion
  TH2                   *fHPionEtaPhiConvConv;            //!histo for pion eta vs. phi
  TH2                   *fHPionMggPtConvConv;             //!histo for pion mass vs. pT
  TH2                   *fHPionMggAsymConvConv;           //!histo for pion mass vs. asym
  TH2                   *fHPionMggDggConvConv;            //!histo for pion mass vs. opening angle
  TH2                   *fHPionInvMassesConvConv;     //!histos for invariant mass plots
  // mixing conv conv
  TH2                   *fHPionInvMassesConvConvMix;     //!histos for invariant mass plots
  
  // conversion info
  TH1                   * fHConversionPoint;   //!histo for conversion position in XY
  
  // histograms for MC
  TH1                   *fHWgt;         //!histo for weight of particles
  TH1                   *fHPionTruthPt;       //!histo for pT from MC pion
  TH1                   *fHPionTruthPtIn;    //!histo for pT for MC pion in eta range
  TH1                   *fHPionTruthPtAcc;    //!histo for pT for MC pion in acceptance
  TH1                   *fHPionTruthPtConvAcc;    //!histo for pT for MC pion in acceptance
  TH1                   *fHEtaTruthPt;       //!histo for pT from MC eta
  TH1                   *fHEtaTruthPtIn;    //!histo for pT for MC eta in eta range
  TH1                   *fHEtaTruthPtAcc;    //!histo for pT for MC eta in acceptance
  TH1                   *fHEtaTruthPtConvAcc;    //!histo for pT for MC eta in acceptance
  TH1                   *fHGamTruthPt;       //!histo for pT from MC gamma
  TH1                   *fHGamTruthPtIn;    //!histo for pT for MC gamma in eta range
  TH1                   *fHGamTruthPtAcc;    //!histo for pT for MC gamma in acceptance
  
  TH1                   *fHPionTruthPtAdd;       //!histo for pT from MC pion
  TH1                   *fHPionTruthPtInAdd;    //!histo for pT for MC pion in eta range
  TH1                   *fHPionTruthPtAccAdd;    //!histo for pT for MC pion in acceptance
  TH1                   *fHPionTruthPtConvAccAdd;    //!histo for pT for MC pion in acceptance
  TH1                   *fHEtaTruthPtAdd;       //!histo for pT from MC eta
  TH1                   *fHEtaTruthPtInAdd;    //!histo for pT for MC eta in eta range
  TH1                   *fHEtaTruthPtAccAdd;    //!histo for pT for MC eta in acceptance
  TH1                   *fHEtaTruthPtConvAccAdd;    //!histo for pT for MC eta in acceptance
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
  
  // cluster studies
  TH2                   *fHClustEnergyM02Gamma;        //!histo for cluster energy vs. M02 for direct photons
  TH2                   *fHClustEnergyM02Pi0;        //!histo for cluster energy vs. M02 for pi0 photons
  TH2                   *fHClustEnergyM02Pion;        //!histo for cluster energy vs. M02 for charged pions 
  TH2                   *fHClustEnergyM02PionTM;	//!histo for cluster energy vs. M02 for charged pions 
  TH2                   *fHClustEnergyM02AllTM;		//!histo for cluster energy vs. M02 for charged pions 
  TH2                   *fHClustEnergyM02GammaAll;        //!histo for cluster energy vs. M02 for all gammas
  TH2                   *fHClustEnergyM02All;		//!histo for cluster energy vs. M02 for all cluster
  TH2                   *fHClustEnergyM02v0;		//!histo for cluster energy vs. M02 for v0 tagged electrons
  TH2                   *fHClustEnergyEPv0;  	//!histo for cluster energy vs. M02 for v0 tagged electrons
  TH1                   *fHClustM02Gamma0TM;	//!histo for m02 1d histogram
  TH1                   *fHClustM02Gamma1TM;	//!histo for m02 1d histogram
  TH1                   *fHClustM02Gamma2TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Gamma3TM;	//!histo for m02 1d histogram
  TH1                   *fHClustM02Gamma4TM;	//!histo for m02 1d histogram 
  TH1                   *fHClustM02Gamma5TM;	//!histo for m02 1d histogram   
  TH1                   *fHClustM02Gamma6TM;	//!histo for m02 1d histogram    
  TH1                   *fHClustM02Gamma7TM;	//!histo for m02 1d histogram   
  TH1                   *fHClustM02Pion0TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Pion1TM;	//!histo for m02 1d histogram
  TH1                   *fHClustM02Pion2TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Pion3TM;	//!histo for m02 1d histogram 
  TH1                   *fHClustM02Pion4TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Pion5TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Pion6TM;	//!histo for m02 1d histogram  
  TH1                   *fHClustM02Pion7TM;	//!histo for m02 1d histogram   
      
 TH2                   *fHClustEnergyM02GammaRaw;        //!histo for cluster energy vs. M02 for direct photons
 TH2                   *fHClustEnergyM02Pi0Raw;        //!histo for cluster energy vs. M02 for pi0 photons
 TH2                   *fHClustEnergyM02PionRaw;        //!histo for cluster energy vs. M02 for charged pions 
 TH2                   *fHClustEnergyM02GammaAllRaw;        //!histo for cluster energy vs. M02 for all gammas
 TH2                   *fHClustM02M20NoGammaRaw;        //!histo for cluster energy vs. M02 for charged pions 
 TH2                   *fHClustM02M20GammaAllRaw;        //!histo for cluster energy vs. M02 for all gammas 
 TH2                   *fHClustEnergyM02AllRaw;		//!histo for cluster energy vs. M02 for all cluster
 TH2                   *fHClustEnergyM02GammaSmallCut; //!histo for cluster energy vs. M02 for Gammas after some cuts
 TH2                   *fHClustEnergyM02PionSmallCut; //!histo for cluster energy vs. M02 for Pions after some cuts
  
 
 TH1		      *fHv0electrons;		//! histo to show v0 tagged electrons
 TH1		      *fHPtSpecAll;        //! histo to show Pt Spectrum for all cluster
 TH1 			  *fHPtSpecGamma;  //! histo to show Pt Spectrum for all photons
 TH1		      *fHPtSpecPion;        //! histo to show Pt Spectrum for charged pions  
 TH1		      *fHPtSpecElectron;		//! histo to show Pt Spectrum for electrons
 TH1		      *fHPtSpecMyon;		//! histo to show Pt Spectrum for myons  
 TH1		      *fHPtSpecProton;		//! histo to show Pt Spectrum for protons  
 TH1		      *fHPtSpecNeutron;		//! histo to show Pt Spectrum for neutrons  
 TH1		      *fHPtSpecKaon;		//! histo to show Pt Spectrum for charged kaons
 TH1		      *fHPtSpecKaon0;		//! histo to show Pt Spectrum for neutral kaons
 TH1		      *fHPtSpecNoGamma;		//! histo to show Pt Spectrum for particles that are not gammas
 TH1		      *fHPtSpecCharged; 		//! histo to show Pt Spectrum for charged particles
 TH1		      *fHPtSpecElectronTM;		//! histo to show Pt Spectrum after Trackamtching
 TH1		      *fHPtSpecPionTM;		//! histo to show Pt Spectrum after Trackamtching
 TH1		      *fHPtSpecElectronNoTM;		//! histo to show Pt Spectrum after Trackamtching
 TH1		      *fHPtSpecPionNoTM;		//! histo to show Pt Spectrum after Trackamtching 
 TH1		      *fHPtSpecGammaNoM02; 	//! histo to show Pt Spectrum after Trackamtching
 TH1		      *fHPtSpecPionNoM02; 	//! histo to show Pt Spectrum after Trackamtching
 TH1		      *fHPtSpecGammaM02; 	//! histo to show Pt Spectrum after M02Cut
 TH1		      *fHPtSpecPionM02; 	//! histo to show Pt Spectrum after M02Cut

 TH1          *fHPtSpecM02Cut0Pion; //! histo to study cut
 TH1          *fHPtSpecM02Cut1Pion; //! histo to study cut
 TH1          *fHPtSpecM02Cut0Photon; //! histo to study cut
 TH1          *fHPtSpecM02Cut1Photon; //! histo to study cut
  
  TH1          *fHPtSpecTrackCut0Pion; //! histo to study cut
  TH1          *fHPtSpecTrackCut1Pion; //! histo to study cut
  TH1          *fHPtSpecTrackCut0Photon; //! histo to study cut
  TH1          *fHPtSpecTrackCut1Photon; //! histo to study cut

  TH1		      *fHPtSpecAllRaw;        //! histo to show Pt Spectrum for all cluster
 TH1 			  *fHPtSpecGammaRaw;  //! histo to show Pt Spectrum for all photons
 TH1		      *fHPtSpecPionRaw;        //! histo to show Pt Spectrum for charged pions  
 TH1		      *fHPtSpecElectronRaw;		//! histo to show Pt Spectrum for electrons
 TH1		      *fHPtSpecMyonRaw;		//! histo to show Pt Spectrum for myons  
 TH1		      *fHPtSpecProtonRaw;		//! histo to show Pt Spectrum for protons  
 TH1		      *fHPtSpecNeutronRaw;		//! histo to show Pt Spectrum for neutrons  
 TH1		      *fHPtSpecKaonRaw;		//! histo to show Pt Spectrum for kaons
 TH1		      *fHPtSpecKaon0Raw;		//! histo to show Pt Spectrum for kaons
 TH1		      *fHPtSpecNoGammaRaw; //! histo to show Pt Spectrum for kaons
 TH1		      *fHPtSpecChargedRaw; 		//! histo to show Pt Spectrum for charged particles
 
  
  TH1 				*fHPtSpecSysEnergy1;  //! effieciency 
  TH1 				*fHPtSpecSysEnergy2;  //! effieciency   
  TH1 				*fHPtSpecSysEnergy3;  //! effieciency 
  TH1 				*fHPtSpecSysEnergy4;  //! effieciency 
  TH1 				*fHPtSpecSysEnergy5;  //! effieciency     
  TH1 				*fHPtSpecSysNcell1;  //! effieciency 
  TH1 				*fHPtSpecSysNcell2;  //! effieciency   
  TH1 				*fHPtSpecSysNcell3;  //! effieciency 
  TH1 				*fHPtSpecSysBorder1;  //! effieciency 
  TH1 				*fHPtSpecSysBorder2;  //! effieciency   
  TH1 				*fHPtSpecSysBorder3;  //! effieciency  
 
  
    
  TH1 			  *fHPtSpecGammaCompare;  //! histo to show Pt Spectrum for all photons
  TH1		      *fHPtSpecPionCompare;        //! histo to show Pt Spectrum for charged pions   
  TH1		      *fHPtSpecCompare;        //! histo to show Pt Spectrum all particles
  
  
  TH1		      *fHPtSpecEffParticle;  //! histo to show Pt Spectrum for all particles in EmCal Eta/Phi
  TH1		      *fHPtSpecEffPhoton;  //! histo to show Pt Spectrum for photons in EmCal Eta/Phi  
  TH1		      *fHPtSpecAccPhoton;		//! histo to show Pt Spectrum for photons overall  
  TH1		      *fHPtSpecEtaPhoton;		//! histo to show Pt Spectrum for photons overall 
  TH1		      *fHPtSpecPhiPhoton;		//! histo to show Pt Spectrum for photons overall 
  TH2 			  *fHGenEtaPhi;				//! histo to show Pt Spectrum for photons overall 
  TH1		      *fHPtSpecEffCluster;  	//! histo to show Pt Spectrum for cluster photons in EmCal Eta/Phi 
  TH1		      *fHPtSpecEffNeutron; 	//!histo for neutron/antineutron efficiency correction
  TH1		      *fHPtSpecConversion;  //!histo for conversion correction
  TH1		      *fHPtSpecConversionNot; //!histo for conversion correction
  TH1		      *fHPtSpecEffPhotonEta5;  //! histo to show Pt Spectrum for photons in EmCal Eta (+- 0.5)/Phi  
  TH1		      *fHPtSpecEffPhotonEta4;  //! histo to show Pt Spectrum for photons in EmCal Eta (+- 0.4)/Phi  
  TH1		      *fHPtSpecEffPhotonEta3;  //! histo to show Pt Spectrum for photons in EmCal Eta (+- 0.3)/Phi  
  TH1		      *fHPtSpecEffPhotonEta2;  //! histo to show Pt Spectrum for photons in EmCal Eta (+- 0.2)/Phi  
  TH1		      *fHPtSpecEffPhotonEta1;  //! histo to show Pt Spectrum for photons in EmCal Eta (+- 0.1)/Phi  
		  
  TH1		      *fHPtSpecDecayPi0;	//! histo to show Pt Spectrum for photons from pi0 decay
  TH1		      *fHPtSpecDecayEta;	//! histo to show Pt Spectrum for photons from eta decay
  TH1		      *fHPtSpecDecayOmega;	//! histo to show Pt Spectrum for photons from omega decay 
  TH1		      *fHPtSpecDecayEtap;	//! histo to show Pt Spectrum for photons from eta prime decay     
  
  
  
  TH1		      *fHCutVariationM02Photon; //!cutvariations
  TH1		      *fHCutVariationM02Pion; //!cutvariations
  TH1		      *fHCutVariationM02PhotonTest; //!cutvariations
  TH1		      *fHCutVariationM02PionTest; //!cutvariations  
  
  TH1		      *fHM02Photon;//!cutvariations
  TH1		      *fHM02Pion; //!cutvariations
 
  
  TH1		      *fHCutVariationPion; 			//!cutvariations
  TH1		      *fHCutVariationPionM021; 			//!cutvariations
  TH1		      *fHCutVariationPionM022; 			//!cutvariations
  TH1		      *fHCutVariationPionM023; 			//!cutvariations
  TH1		      *fHCutVariationPionM024; 			//!cutvariations
  TH1		      *fHCutVariationPionM025; 			//!cutvariations
  TH1		      *fHCutVariationPionM026; 			//!cutvariations  
  TH1		      *fHCutVariationPionM027; 			//!cutvariations    
  

  TH1		      *fHCutVariationPionEnergy1; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy2; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy3; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy4; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy5; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy6; 			//!cutvariations  
  TH1		      *fHCutVariationPionEnergy7; 			//!cutvariations  
   

  TH1		      *fHCutVariationPionNcell1; 			//!cutvariations   
  TH1		      *fHCutVariationPionNcell2; 			//!cutvariations   
  TH1		      *fHCutVariationPionNcell3; 			//!cutvariations   
  TH1		      *fHCutVariationPionNcell4; 			//!cutvariations   
  TH1		      *fHCutVariationPionNcell5; 			//!cutvariations        
  TH1		      *fHCutVariationPionNcell6; 			//!cutvariations            
  TH1		      *fHCutVariationPionNcell7; 			//!cutvariations   
  
  TH1		      *fHCutVariationPhoton; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM021; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM022; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM023; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM024; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM025; 			//!cutvariations
  TH1		      *fHCutVariationPhotonM026; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonM027; 			//!cutvariations    
  

  TH1		      *fHCutVariationPhotonEnergy1; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy2; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy3; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy4; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy5; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy6; 			//!cutvariations  
  TH1		      *fHCutVariationPhotonEnergy7; 			//!cutvariations  
   
  TH1		      *fHCutVariationPhotonNcell1; 			//!cutvariations   
  TH1		      *fHCutVariationPhotonNcell2; 			//!cutvariations   
  TH1		      *fHCutVariationPhotonNcell3; 			//!cutvariations   
  TH1		      *fHCutVariationPhotonNcell4; 			//!cutvariations   
  TH1		      *fHCutVariationPhotonNcell5; 			//!cutvariations        
  TH1		      *fHCutVariationPhotonNcell6; 			//!cutvariations            
  TH1		      *fHCutVariationPhotonNcell7; 			//!cutvariations   
  
      
  // more histograms
  TH1		      *fHMixRotation;        //! histo to show how much the mixed events were rotated in phi
  TH1             *fHCorrection;         //! histo to show the correction factor
  TH2             *fHPionSm;             //! histo to see the change of the pair mass due to the correction
  TH1		      *fHParticles;        //! histo to show all particles that formed a cluster   
  TH1		      *fHParticlesTM;	  //! histo to show all particles that formed a cluster 
  TH1		      *fHParticlesNoTM;	  //! histo to show all particles that formed a cluster  
  
  TH1		      *fHParticlesRaw;        //! histo to show all particles that formed a cluster  
  TH1		      *fHParticlesEff;        //! histo to show all particles that formed a cluster    
  TH1		      *fHParticlesCompare;        //! histo to show all particles that formed a cluster   
  TH1		      *fHParticlesCompare2;        //! histo to show all particles that formed a cluster    
  TH1		      *fHParticlesCompare3;        //! histo to show all particles that formed a cluster        
  TH1		      *fHParticleR;  //! histo to show range from electron vertices
  TH1		      *fHParticleRcut;  //! histo to show range from electron vertices afet cut  
  TH1		      *fHParticleRcutcut;  //! histo to show range from electron vertices afet cut
  TH1		      *fHParticleRcutcutcut;  //! histo to show range from electron vertices afet cut     
  TH1		      *fHMotherR;  //! histo to show range from electron vertices  
  TH1		      *fHMotherR2;  //! histo to show range from electron vertices 
  TH1		      *fHClustNcell;  //! histo with the number of cells per cluster
  TH1		      *fHClustNcellPhoton;  //! histo with the number of cells per cluster  
  TH1		      *fHClustNcellNoPhoton;  //! histo with the number of cells per cluster 
  TH1		      *fHClustNcellPhotonCut;  //! histo with the number of cells per cluster after cut
  TH1		      *fHClustNcellNoPhotonCut;  //! histo with the number of cells per cluster after cut  
  TH1		      *fHGammaMIP; //! histo with the number gammas with low Ec
  TH1		      *fHHadronMIP; //! histo with the number hadrons with low Ec
    
  TH1		      *fHPtSpecElectronMerge; //! asdfasdaf
  TH2			  *fHClustElectronZR;		//! asdfasdaf
  TH1		      *fHclusterTOFdistance; //! histo to show distance from emcal cluster to tof cluster
  TH1		      *fHistTOF; //! histo to show distance from emcal cluster to tof cluster
  
  
  // store start and end of pythia particles, added signals, ...
  Int_t ipymin; // first pythia particle index
  Int_t ipymax; // last pythia particle index
  Int_t ipi0min; // first added pi0 particle index
  Int_t ipi0max; // last added pi0 particle index
  Int_t ietamin; // first added eta particle index
  Int_t ietamax; // last added eta particle index
  
  const static int nMulClass =   4;
  const static int nZClass   =   4;
  const static int nPtClass = 1;
  int iEvt[nMulClass][nZClass][nPtClass];
  const static int nEvt      =   30; // mixing "depth"
  
private:
  AliAnalysisTaskEMCALDirGamma(const AliAnalysisTaskEMCALDirGamma&);            // not implemented
  AliAnalysisTaskEMCALDirGamma &operator=(const AliAnalysisTaskEMCALDirGamma&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALDirGamma, 14) // Analysis task for neutral pions in Pb+Pb
	    AliESDv0KineCuts *fV0cuts;                //! ESD V0 cuts
};



#endif


