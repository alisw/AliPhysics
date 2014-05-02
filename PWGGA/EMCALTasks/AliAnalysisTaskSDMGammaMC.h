#ifndef ALIANALYSISTASKSDMGAMMAMC_H
#define ALIANALYSISTASKSDMGAMMAMC_H

class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TH3D;
class TNtuple;
class TList;
class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDCaloCluster;
class AliAODCaloCluster;
class AliMCEvent;
class AliMCParticle;
class AliEMCALGeometry;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskSDMGammaMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSDMGammaMC();
  AliAnalysisTaskSDMGammaMC(const char *name);
  virtual ~AliAnalysisTaskSDMGammaMC();
    
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
    
  void         SetMcMode(Bool_t b) { fMcMode = b; }
  void         SetRecalScheme(Int_t kRecalibrator) { fRecalibrator = kRecalibrator; }
  void         SetFidPhiMinMax(Double_t kPhimin, Double_t kPhimax){ fPhimin = kPhimin; fPhimax = kPhimax; }
  void         SetFidEtaMinMax(Double_t kEtamin, Double_t kEtamax){ fEtamin = kEtamin; fEtamax = kEtamax; }
  
 private:
  static const int zvtx_bins = 16;
  static const int mult_bins = 25;
  static const unsigned int poolDepth = 25;
  
  Int_t GetMultBin(Int_t mult);
  Int_t GetZvtxBin(Double_t vertZ);
  Int_t isGoodEsdCluster(AliESDCaloCluster* esdclust);
  Int_t isGoodAodCluster(AliAODCaloCluster* aodclust);
  Double_t getDeltaPhi(TLorentzVector p1, TLorentzVector p2);
  Double_t getDeltaEta(TLorentzVector p1, TLorentzVector p2);
  Double_t PrivateEnergyRecal(Double_t energy, Int_t iCalib);

  TList           *fOutput;        //! Output list
  Bool_t           fMcMode;                 // monte carlo mode
  Int_t            fRecalibrator;           // custom recalibrator? 
  Double_t         fPhimin;           // Cuts. 
  Double_t         fPhimax;           // Cuts. 
  Double_t         fEtamin;           // Cuts. 
  Double_t         fEtamax;           // Cuts. 
  AliESDtrackCuts *fTrackCuts;     // Track cuts
  AliESDEvent     *fEsdEv;                  //!pointer to input esd event
  AliAODEvent     *fAodEv;                  //!pointer to input aod event
  TH1F            *h1_nClusters;  //!  # of clusters/evt
  TH1F            *h1_zvtx;  //!  # of clusters/evt
  TH1F            *h1_trigger;  //!  # of clusters/evt
  TH1F            *h1_E;        //!  energy spectrum
  TH1F            *h1_Phi;        //!  Pt spectrum

  TH2F            *h2_PiMotherID; //! first bin is primaries, second is non-primaries
  TH2F            *h2_GaMotherID; //! first bin is primaries, second is non-primaries
  TH3F            *h3_gE_RecTruth; //! E_in/E_rec vs E_rec (for 4 categories)
  TH3F            *h3_gE_RecTruth_ncellscut; //! E_in/E_rec vs E_rec (for 4 categories)
  
  TH1F            *h1_Pi0TruthPt;        //!  Pt spectrum from MC! 
  TH1F            *h1_PriPi0TruthPt;        //!  Pt spectrum from MC! 
  TH1F            *h1_Pi0TruthPtEmcal;        //!  Pt spectrum from MC! 
  TH1F            *h1_PriPi0TruthPtEmcal;        //!  Pt spectrum from MC! 

  TH2F            *h2_Pi0TruthPhiEta;    //!  etaphi spectrum from MC! 
  TH2F            *h2_PriPi0TruthPhiEta;    //!  etaphi spectrum from MC! 
  TH2F            *h2_Pi0TruthPhiEtaEmcal;    //!  etaphi spectrum from MC! 
  TH2F            *h2_PriPi0TruthPhiEtaEmcal;    //!  etaphi spectrum from MC! 

  TH1F            *h1_TruthPhotonsEmcal;        //!  Pt spectrum from MC! 
  TH2F            *h2_TruthPhotonsPhiEta;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhotonsEmcal;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhotonsNCellsCut;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhotonsTrackMatchCut;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhotonsAllCut;        //!  Pt spectrum from MC! 
  TH2F            *h2_PhotonsPhiEtaIsEmcal;        //!  Pt spectrum from MC! 

  TH1F            *h1_dR_RealMC;        //!  Pt spectrum from MC! 

  TH1F            *h1_Eta;       //!  pseudorapidity spectrum
  TH1F            *h1_Chi2;       //!  pseudorapidity spectrum
  TH1F            *h1_nTrkMatch;       //!  pseudorapidity spectrum
  TH1F            *h1_nCells;       //!  pseudorapidity spectrum
  TH1F            *h1_ClusterDisp;       //!  cluster dispersion
  TH2F            *h2_Ellipse;       //!  ellipse axis?
  TH2F            *h2_EtaPt;       //!  2d histogram Y - pseudorap spectrum
  TH2F            *h2_dphi_deta;       //!  2dimensional E vs mom
  TH2F            *h2_dphi_deta_mix;       //!  2dimensional E vs mom
  TH2F            *h2_DispRes;       //!  2dimensional E vs mom
  TH2F            *h2_cells_M02;       //!  

  std::vector<TLorentzVector> Photons[poolDepth][zvtx_bins][mult_bins];
  std::vector<Int_t> TriggerList;
  
  AliAnalysisTaskSDMGammaMC(const AliAnalysisTaskSDMGammaMC&); // not implemented
  AliAnalysisTaskSDMGammaMC& operator=(const AliAnalysisTaskSDMGammaMC&); // not implemented
    
  ClassDef(AliAnalysisTaskSDMGammaMC, 1); // example of analysis
};
#endif
