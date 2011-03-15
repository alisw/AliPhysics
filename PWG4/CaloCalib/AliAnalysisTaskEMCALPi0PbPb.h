#ifndef AliAnalysisTaskEMCALPi0PbPb_cxx
#define AliAnalysisTaskEMCALPi0PbPb_cxx

// $Id$

class TAxis;
class TH1F;
class TH2F;
class TClonesArray;
class TObjArray;
class AliAODCaloCells;
class AliAODCaloCluster;
class AliAODEvent;
class AliEMCALGeoUtils;
class AliESDCaloCells;
class AliESDCaloCluster;
class AliESDEvent;

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
  void         SetUseQualFlag(Bool_t b)                       { fUseQualFlag = b;           }
  void         SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }

 protected:
  virtual void FillCellHists();
  virtual void FillClusHists();
  virtual void FillPionHists();
  Double_t     GetMaxCellEnergy(AliVCluster *c);
  Double_t     GetSigmaMax(AliVCluster *c);

    // input members
  Double_t               fAsymMax;              // energy asymmetry max (def=1)
  TString                fCentVar;              // variable for centrality determination
  Double_t               fCentFrom;             // min centrality (def=0)
  Double_t               fCentTo;               // max centrality (def=100)
  Double_t               fVtxZMin;              // min primary vertex z (def=-10cm)
  Double_t               fVtxZMax;              // max primary vertex z (def=+10cm)
  Bool_t                 fUseQualFlag;          // if true use quality flag for centrality
  TString                fClusName;             // cluster branch name (def="")
     // derived members (ie with ! after //)
  AliEMCALGeoUtils      *fGeom;                 //! geometry utils
  TList                 *fOutput;               //!container of output histograms
  AliESDEvent           *fEsdEv;                //!pointer to input esd event
  AliAODEvent           *fAodEv;                //!pointer to input aod event
  TObjArray             *fRecPoints;            //!pointer to rec points (AliAnalysisTaskEMCALClusterizeFast)
  TObjArray             *fEsdClusters;          //!pointer to esd clusters
  AliESDCaloCells       *fEsdCells;             //!pointer to esd cells
  TObjArray             *fAodClusters;          //!pointer to aod clusters
  AliAODCaloCells       *fAodCells;             //!pointer to aod cells
  TAxis                 *fPtRanges;             //!pointer to pt ranges
    // histograms
  TH1F                  *fHCuts;                //!histo for cuts
  TH1F                  *fHVertexZ;             //!histo for vtxz
  TH1F                  *fHVertexZ2;            //!histo for vtxz after vtx cuts
  TH1F                  *fHCent;                //!histo for cent
  TH1F                  *fHCentQual;            //!histo for cent after quality flag cut
    // histograms for cells
  TH2F                 **fHColuRow;             //!histo for cell column and row
  TH2F                 **fHColuRowE;            //!histo for cell column and row weight energy
  TH1F                 **fHCellMult;            //!histo for cell multiplicity in module
  TH1F                  *fHCellE;               //!histo for cell energy
  TH1F                  *fHCellH;               //!histo for highest cell energy
    // histograms for clusters
  TH2F                  *fHClustEtaPhi;         //!histo for cluster eta vs. phi
  TH2F                  *fHClustEnergyPt;       //!histo for cluster energy vs. pT
  TH2F                  *fHClustEnergySigma;    //!histo for cluster energy vs. variance over long axis 
  TH2F                  *fHClustSigmaSigma;     //!histo for sigma vs. lambda_0 comparison
  TH2F                  *fHClustNTowEnergyRatio;//!histo for cluster n tow vs. energy ratio
    // histograms for pion candidates
  TH2F                  *fHPionEtaPhi;          //!histo for pion eta vs. phi
  TH2F                  *fHPionMggPt;           //!histo for pion mass vs. pT
  TH2F                  *fHPionMggAsym;         //!histo for pion mass vs. asym
  TH1F                  *fHPionInvMasses[20];   //!histos for invariant mass plots 

 private:
  AliAnalysisTaskEMCALPi0PbPb(const AliAnalysisTaskEMCALPi0PbPb&);            // not implemented
  AliAnalysisTaskEMCALPi0PbPb &operator=(const AliAnalysisTaskEMCALPi0PbPb&); // not implemented

  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 1); // Analysis task for neutral pions in Pb+Pb
};
#endif
