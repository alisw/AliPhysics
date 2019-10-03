#ifndef ALIANALYSISTASKEMCALMESONGGSDM_H
#define ALIANALYSISTASKEMCALMESONGGSDM_H

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
#include "AliAnalysisUtils.h"

class AliAnalysisTaskEMCALMesonGGSDM : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALMesonGGSDM();
  AliAnalysisTaskEMCALMesonGGSDM(const char *name);
  virtual ~AliAnalysisTaskEMCALMesonGGSDM();
    
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
    
  void         SetMcMode(Bool_t b)                 { fMcMode       = b;             }
  void         SetRecalScheme(Int_t kRecalibrator) { fRecalibrator = kRecalibrator; }
  void         SetMyMCType(char *kMyMCType)        { fMyMCType     = kMyMCType; }
  void         SetdRmin_ClustTrack(Double_t kdRmin_ClustTrack)    { fdRmin_ClustTrack = kdRmin_ClustTrack; }
  void         SetFidPhiMinMax(Double_t kPhimin, Double_t kPhimax){ fPhimin = kPhimin; fPhimax = kPhimax; }
  void         SetFidEtaMinMax(Double_t kEtamin, Double_t kEtamax){ fEtamin = kEtamin; fEtamax = kEtamax; }
  
 private:
  static const int zvtx_bins = 8;
  static const int mult_bins = 7;
  static const unsigned int poolDepth = 80;
  
  Int_t GetMultBin(Int_t mult);
  Int_t GetZvtxBin(Double_t vertZ);
  Int_t isGoodEsdCluster(AliESDCaloCluster* esdclust);
  Int_t isGoodAodCluster(AliAODCaloCluster* aodclust);
  Double_t getDeltaPhi(TLorentzVector p1, TLorentzVector p2);
  Double_t getDeltaEta(TLorentzVector p1, TLorentzVector p2);
  Double_t PrivateEnergyRecal(Double_t energy, Int_t iCalib);
  Double_t GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const;
  Int_t IsPhysPrimJ(AliMCEvent *mcEvent, Int_t iTrack);
  Int_t IsLongLivedOrK(Int_t MyPDGcode);
  Int_t IsMyMCHeaderType(Int_t iTrack, char *MyType, AliMCEvent *mcEvent) const;

  TList           *fOutput;        //! Output list
  Bool_t           fMcMode;                 // monte carlo mode
  char            *fMyMCType;               // monte carlo primary particles
  Int_t            fRecalibrator;           // custom recalibrator? 
  Double_t         fdRmin_ClustTrack; // Cuts. 
  Double_t         fPhimin;           // Cuts. 
  Double_t         fPhimax;           // Cuts. 
  Double_t         fEtamin;           // Cuts. 
  Double_t         fEtamax;           // Cuts. 
  AliESDtrackCuts *fTrackCuts;              // Track cuts
  AliESDEvent     *fEsdEv;                  //!pointer to input esd event
  AliAODEvent     *fAodEv;                  //!pointer to input aod event
  TH1F            *h1_nClusters;  //!  # of clusters/evt
  TH1F            *h1_zvtx;  //!  # of clusters/evt
  TH1F            *h1_trigger;  //!  # of clusters/evt
  TH1F            *h1_M;        //!  Mass spectrum
  TH1F            *h1_M_mix;    //!  Mass spectrum
  TH1F            *h1_E;        //!  energy spectrum
  TH2F            *h2_PhiEtaCluster;        //!  phi vs eta for the cluster CoG
  TH2F            *h2_PhiEtaClusterCut;     //!  phi vs eta for the cluster CoG w/ cuts
  TH2F            *h2_PhiEtaMaxCell;        //!  phi vs eta for the maximum cell
  TH2F            *h2_PhiEtaMaxCellCut;     //!  phi vs eta for the maximum cell w/ cuts
  TH1F            *h1_dR_ClustTrk;        //!  Track Matching
  TH2F            *h2_gE_RecTruth; //! gamma E, rec/truth, first bin is primaries, second is non-primaries
  TH2F            *h2_eop_E;        //!  e over p vs E. simple
  TH2F            *h2_eop_pT;        //!  e over p vs pT. simple
  TH2F            *h2_E_time;        //!  cluster energy vs time.
  
  TH1F            *h1_Pi0TruthPt;        //!  Pt spectrum from MC! 
  TH1F            *h1_K0Pi0TruthPt;        //!  Pt spectrum from MC! 
  TH1F            *h1_PriPi0TruthPt;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhysPi0TruthPt;        //!  Pt spectrum from MC! 

  TH1F            *h1_Pi0TruthPtEmcal;        //!  Pt spectrum from MC! 
  TH1F            *h1_K0Pi0TruthPtEmcal;        //!  Pt spectrum from MC! 
  TH1F            *h1_PriPi0TruthPtEmcal;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhysPi0TruthPtEmcal;        //!  Pt spectrum from MC! 

  TH1F            *h1_Pi0TruthPtPhi2piEta065;        //!  Pt spectrum from MC! 
  TH1F            *h1_K0Pi0TruthPtPhi2piEta065;        //!  Pt spectrum from MC! 
  TH1F            *h1_PriPi0TruthPtPhi2piEta065;        //!  Pt spectrum from MC! 
  TH1F            *h1_PhysPi0TruthPtPhi2piEta065;        //!  Pt spectrum from MC! 

  TH1F            *h1_Pi0TruthPtPhi2piEta1;        //!  Pt spectrum from MC!   
  TH1F            *h1_K0Pi0TruthPtPhi2piEta1;        //!  Pt spectrum from MC!   
  TH1F            *h1_PriPi0TruthPtPhi2piEta1;        //!  Pt spectrum from MC!   
  TH1F            *h1_PhysPi0TruthPtPhi2piEta1;        //!  Pt spectrum from MC!   

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
  TH2F            *h2_Mpt_Pri;       //!  2dimensional mass vs mom primary pions
  TH2F            *h2_Mpt_Sec;       //!  2dimensional mass vs mom secondary pions
  TH3F            *h3_MptR_Sec;       //!  2dimensional mass vs production radius, secondary pions
  TH3F            *h3_MptR_K0s;       //!  2dimensional mass vs production radius, K0s pions
  TH3F            *h3_MptR_Mat;       //!  2dimensional mass vs production radius, pions from material
  TH2F            *h2_PtR_MatM;       //!  2dimensional pt vs production radius, pions from material (that merged). pi mass assumed.
  TH2F            *h2_Mpt_Pri_conv;       //!  2dimensional mass vs mom primary pions
  TH2F            *h2_Mpt_Sec_conv;       //!  2dimensional mass vs mom secondary pions
  TH3F            *h3_MptR_Sec_conv;       //!  2dimensional mass vs production radius, secondary pions
  TH3F            *h3_MptR_K0s_conv;       //!  2dimensional mass vs production radius, K0s pions
  TH3F            *h3_MptR_Mat_conv;       //!  2dimensional mass vs production radius, pions from material
  TH1F            *h1_eConversionR;        //!  conversion point (radius)
  TH1F            *h1_PriPi0Mother;       //!  the parent ID of every sec pi0 mother
  TH1F            *h1_SecPi0Mother;       //!  the parent ID of every pri pi0 mother

  TH1F            *h1_Chi2;       //!  pseudorapidity spectrum
  TH1F            *h1_nTrkMatch;       //!  pseudorapidity spectrum
  TH1F            *h1_nCells;       //!  pseudorapidity spectrum
  TH1F            *h1_ClusterDisp;       //!  cluster dispersion
  TH2F            *h2_Ellipse;       //!  ellipse axis?
  TH2F            *h2_EtaPt;       //!  2d histogram Y - pseudorap spectrum
  TH3F            *h3_MptAsymm;       //!  2dimensional E vs mom
  TH3F            *h3_MptAsymm_mix;       //!  2dimensional E vs mom
  TH2F            *h2_dphi_deta;       //!  2dimensional E vs mom
  TH2F            *h2_dphi_deta_mix;       //!  2dimensional E vs mom
  TH2F            *h2_DispRes;       //!  2dimensional E vs mom
  TH2F            *h2_cells_M02;       //!  

  //AliPIDCombined *fPIDCombined; // for E/p
  //AliPIDResponse *fPIDResponse; // for E/p

  std::vector<TLorentzVector> Photons[poolDepth][zvtx_bins][mult_bins]; //!
  std::vector<Int_t> TriggerList; //!
  
  AliAnalysisUtils*   fHelperClass;           //! Vertex selection helper

  AliAnalysisTaskEMCALMesonGGSDM(const AliAnalysisTaskEMCALMesonGGSDM&); // not implemented
  AliAnalysisTaskEMCALMesonGGSDM& operator=(const AliAnalysisTaskEMCALMesonGGSDM&); // not implemented
    
  ClassDef(AliAnalysisTaskEMCALMesonGGSDM, 1); // example of analysis
};
#endif
