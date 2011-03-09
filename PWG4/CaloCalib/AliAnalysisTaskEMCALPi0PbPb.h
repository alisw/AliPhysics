#ifndef AliAnalysisTaskEMCALPi0PbPb_cxx
#define AliAnalysisTaskEMCALPi0PbPb_cxx

// $Id$

class TH1F;
class TH2F;
class TClonesArray;
class TObjArray;
class AliAODCaloCells;
class AliAODCaloCluster;
class AliAODEvent;
class AliESDCaloCells;
class AliESDCaloCluster;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPi0PbPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPi0PbPb();
  AliAnalysisTaskEMCALPi0PbPb(const char *name);
  virtual ~AliAnalysisTaskEMCALPi0PbPb(); 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetCentrality(const char *name)                { fCentVar = name;            }
  void SetCentralityRange(Double_t from, Double_t to) { fCentFrom=from; fCentTo=to; }
  void SetClusName(const char *name)                  { fClusName = name;           }
  void SetVertexRange(Double_t z1, Double_t z2)       { fVtxZMin=z1; fVtxZMax=z2;   }

 protected:
    // input members
  TString                fCentVar;              // variable for centrality determination
  Double_t               fCentFrom;             // min centrality (def=0)
  Double_t               fCentTo;               // max centrality (def=100)
  Double_t               fVtxZMin;              // min primary vertex z (def=-7cm)
  Double_t               fVtxZMax;              // max primary vertex z (def=+7cm)
  TString                fClusName;             // cluster branch name (def="")
    // derived members (ie with ! after //)
  TList                 *fOutput;               //!container of output histograms
  AliESDEvent           *fEsdEv;                //!pointer to input esd event
  AliAODEvent           *fAodEv;                //!pointer to input aod event
  const TObjArray       *fEsdClusters;          //!pointer to esd clusters
  const AliESDCaloCells *fEsdCells;             //!pointer to esd cells
  const TObjArray       *fAodClusters;          //!pointer to aod clusters
  const AliAODCaloCells *fAodCells;             //!pointer to aod cells
    // histograms

 private:
  AliAnalysisTaskEMCALPi0PbPb(const AliAnalysisTaskEMCALPi0PbPb&);            // not implemented
  AliAnalysisTaskEMCALPi0PbPb &operator=(const AliAnalysisTaskEMCALPi0PbPb&); // not implemented

  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 1); // Analysis task for neutral pions in Pb+Pb
};

#endif
