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
class AliEMCALGeoUtils;
class AliEMCALRecoUtils;
class AliESDCaloCells;
class AliESDCaloCluster;
class AliESDEvent;
class AliESDTrack;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPi0PbPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPi0PbPb(const char *name=0);
  virtual ~AliAnalysisTaskEMCALPi0PbPb(); 
  
  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  void         SetAsymMax(Double_t asymMax)                   { fAsymMax = asymMax;         }
  void         SetCentrality(const char *name)                { fCentVar = name;            }
  void         SetCentralityRange(Double_t from, Double_t to) { fCentFrom=from; fCentTo=to; }
  void         SetClusName(const char *name)                  { fClusName = name;           }
  void         SetDoAfterburner(Bool_t b)                     { fDoAfterburner = b;         }
  void         SetDoTrackMatWithGeom(Bool_t b)                { fDoTrackMatWithGeom = b;    }
  void         SetDoTrackVtxConstrain(Bool_t b)               { fDoConstrain = b;           }
  void         SetFillNtuple(Bool_t b)                        { fDoNtuple = b;              }
  void         SetGeoName(const char *n)                      { fGeoName = n;               }
  void         SetIsoDist(Double_t d)                         { fIsoDist = d;               }
  void         SetMinClusEnergy(Double_t e)                   { fMinE = e;                  }
  void         SetMinEcc(Double_t ecc)                        { fMinEcc = ecc;              }
  void         SetMinErat(Double_t erat)                      { fMinErat = erat;            }
  void         SetMinNClustersPerTrack(Double_t mct)          { fMinNClustPerTrack = mct;   }
  void         SetMinPtPerMatchedTrack(Double_t mpt)          { fMinPtPerTrack = mpt;       }
  void         SetNminCells(Int_t n)                          { fNminCells = n;             }
  void         SetTrClassNames(const char *n)                 { fTrClassNames = n;          }
  void         SetTrackCuts(AliESDtrackCuts *c)               { fTrCuts = c;                }
  void         SetUseQualFlag(Bool_t b)                       { fUseQualFlag = b;           }
  void         SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }

 protected:
  virtual void CalcClusterProps();
  virtual void CalcTracks();
  virtual void ClusterAfterburner();
  virtual void FillCellHists();
  virtual void FillClusHists();
  virtual void FillPionHists();
  virtual void FillOtherHists();
  Double_t     GetCellIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2)                const;
  Double_t     GetMaxCellEnergy(AliVCluster *c)                                                   const;
  Int_t        GetNCells(AliVCluster *c, Double_t emin=0.)                                        const;
  void         GetSigma(AliVCluster *c, Double_t &sigmaMax, Double_t &sigmaMin)                   const;
  Double_t     GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2)               const;

  class ClusProps {
    public:
      ClusProps() : fTrIndex(-1), fTrDz(-1), fTrDr(-1), fTrDist(-1), fTrEp(0), 
                    fTrIso(0), fTrLowPtIso(0), fCellIso(0) {}
      void Reset() { fTrIndex=-1; fTrDz=-1; fTrDr=-1; fTrDist=-1; fTrEp=0; fTrIso=0; fTrLowPtIso=0; fCellIso=0; }
      Int_t    fTrIndex;
      Double_t fTrDz;
      Double_t fTrDr;
      Double_t fTrDist;
      Double_t fTrEp;
      Double_t fTrIso;
      Double_t fTrLowPtIso;
      Double_t fCellIso;
  };
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
  Double_t               fMinNClustPerTrack;      // minimum number of cluster per track (def=50)
  Double_t               fMinPtPerTrack;          // minimum pT per track (def=0.25 GeV/c)
  Double_t               fIsoDist;                // isolation distance (def=0.2)
  TString                fTrClassNames;           // trigger class names
  AliESDtrackCuts       *fTrCuts;                 // track cuts
  Bool_t                 fDoTrackMatWithGeom;     // track matching including geometry
  Bool_t                 fDoConstrain;            // if true constrain tracks to vertex 

    // derived members (ie with ! after //)
  ULong64_t              fNEvs;                   //!accepted events 
  AliEMCALGeoUtils      *fGeom;                   //!geometry utils
  AliEMCALRecoUtils     *fReco;                   //!geometry utils
  TList                 *fOutput;                 //!container of output histograms
  TObjArray             *fTrClassNamesArr;        //!array of trig class names  
  AliESDEvent           *fEsdEv;                  //!pointer to input esd event
  AliAODEvent           *fAodEv;                  //!pointer to input aod event
  TObjArray             *fRecPoints;              //!pointer to rec points (AliAnalysisTaskEMCALClusterizeFast)
  TObjArray             *fEsdClusters;            //!pointer to esd clusters
  AliESDCaloCells       *fEsdCells;               //!pointer to esd cells
  TObjArray             *fAodClusters;            //!pointer to aod clusters
  AliAODCaloCells       *fAodCells;               //!pointer to aod cells
  TAxis                 *fPtRanges;               //!pointer to pt ranges
  TNtuple               *fNtuple;                 //!pointer to ntuple
  TObjArray             *fSelTracks;              //!pointer to selected tracks
  ClusProps              fClusProps[1000];        //!array of cluster properties
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
  TH2                   *fHClustNCellEnergyRatio; //!histo for cluster n tow vs. energy ratio
    // histograms for pion candidates
  TH2                   *fHPionEtaPhi;            //!histo for pion eta vs. phi
  TH2                   *fHPionMggPt;             //!histo for pion mass vs. pT
  TH2                   *fHPionMggAsym;           //!histo for pion mass vs. asym
  TH2                   *fHPionMggDgg;            //!histo for pion mass vs. opening angle
  TH1                   *fHPionInvMasses[21];     //!histos for invariant mass plots 

 private:
  AliAnalysisTaskEMCALPi0PbPb(const AliAnalysisTaskEMCALPi0PbPb&);            // not implemented
  AliAnalysisTaskEMCALPi0PbPb &operator=(const AliAnalysisTaskEMCALPi0PbPb&); // not implemented

  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 5); // Analysis task for neutral pions in Pb+Pb
};
#endif
