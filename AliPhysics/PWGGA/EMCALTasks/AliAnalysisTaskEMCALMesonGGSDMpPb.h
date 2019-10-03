#ifndef ALIANALYSISTASKEMCALMesonGGSDMPPB_H
#define ALIANALYSISTASKEMCALMesonGGSDMPPB_H

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

class AliAnalysisTaskEMCALMesonGGSDMpPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALMesonGGSDMpPb();
  AliAnalysisTaskEMCALMesonGGSDMpPb(const char *name);
  virtual ~AliAnalysisTaskEMCALMesonGGSDMpPb();
    
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
    
  void         SetMcMode(Bool_t b)                 { fMcMode       = b;             }
  void         SetRecalScheme(Int_t kRecalibrator) { fRecalibrator = kRecalibrator; }
  void         SetdRmin_ClustTrack(Double_t kdRmin_ClustTrack)    { fdRmin_ClustTrack = kdRmin_ClustTrack; }
  void         SetFidPhiMinMax(Double_t kPhimin, Double_t kPhimax){ fPhimin = kPhimin; fPhimax = kPhimax; }
  void         SetFidEtaMinMax(Double_t kEtamin, Double_t kEtamax){ fEtamin = kEtamin; fEtamax = kEtamax; }
  
 private:
  static const int zvtx_bins = 16;
  static const int mult_bins = 25;
  static const unsigned int poolDepth = 25;
  static const int cent_bins = 4;
  
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
  
  TList           *fOutput;        //! Output list
  Bool_t           fMcMode;                 // monte carlo mode
  Int_t            fRecalibrator;           // custom recalibrator? 
  Double_t         fdRmin_ClustTrack; // Cuts. 
  Double_t         fPhimin;           // Cuts. 
  Double_t         fPhimax;           // Cuts. 
  Double_t         fEtamin;           // Cuts. 
  Double_t         fEtamax;           // Cuts. 
  AliESDtrackCuts *fTrackCuts;     //! Track cuts
  AliESDEvent     *fEsdEv;                  //!pointer to input esd event
  AliAODEvent     *fAodEv;                  //!pointer to input aod event
  TH1F           * h1_nClusters[cent_bins];  //! # of clusters/evt
  TH1F            *h1_zvtx;  //! vertex distribution
  TH1F            *h1_trigger;  //! # of clusters/evt
  TH1F            *h1_centrality;  //! # of clusters/evt
  TH1F           * h1_M[cent_bins];        //! Mass spectrum
  TH1F           * h1_M_mix[cent_bins];    //! Mass spectrum
  TH1F           * h1_E[cent_bins];        //! energy spectrum
  TH2F            *h2_PhiEtaCluster;        //!  phi vs eta for the cluster CoG
  TH2F            *h2_PhiEtaClusterCut;     //!  phi vs eta for the cluster CoG w/ cuts
  TH2F            *h2_PhiEtaMaxCell;        //!  phi vs eta for the maximum cell
  TH2F            *h2_PhiEtaMaxCellCut;     //!  phi vs eta for the maximum cell w/ cuts
  TH1F           * h1_dR_ClustTrk[cent_bins];        //! Track Matching
  TH2F            *h2_gE_RecTruth; //! gamma E, rec/truth, first bin is primaries, second is non-primaries
  TH2F            *h2_eop_E;        //!  e over p vs E. simple
  TH2F            *h2_eop_pT;        //!  e over p vs pT. simple
  TH2F            *h2_E_time;        //!  cluster energy vs time.

  TH1F           * h1_Pi0TruthPt[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_PriPi0TruthPt[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_Pi0TruthPtEmcal[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_PriPi0TruthPtEmcal[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_Pi0TruthPtPhi2piEta065[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_Pi0TruthPtPhi2piEta1[cent_bins];        //! Pt spectrum from MC!   

  TH2F            *h2_Pi0TruthPhiEta;    //! etaphi spectrum from MC! 
  TH2F            *h2_PriPi0TruthPhiEta;    //! etaphi spectrum from MC! 
  TH2F            *h2_Pi0TruthPhiEtaEmcal;    //! etaphi spectrum from MC! 
  TH2F            *h2_PriPi0TruthPhiEtaEmcal;    //! etaphi spectrum from MC! 
  TH2F            *h2_Pi0TruthPhiEta_Phi2piEta065;        //! Pt spectrum from MC! 
  TH2F            *h2_Pi0TruthPhiEta_Phi2piEta1;        //! Pt spectrum from MC!   

  TH1F           * h1_TruthPhotonsEmcal[cent_bins];        //! Pt spectrum from MC! 
  TH2F            *h2_TruthPhotonsPhiEta;        //! Pt spectrum from MC! 
  TH1F           * h1_PhotonsEmcal[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_PhotonsNCellsCut[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_PhotonsTrackMatchCut[cent_bins];        //! Pt spectrum from MC! 
  TH1F           * h1_PhotonsAllCut[cent_bins];        //! Pt spectrum from MC! 
  TH2F            *h2_PhotonsPhiEtaIsEmcal;        //! Pt spectrum from MC! 

  TH1F           * h1_dR_RealMC[cent_bins];        //! Pt spectrum from MC! 

  TH1F           * h1_Chi2[cent_bins];       //! pseudorapidity spectrum
  TH1F           * h1_nTrkMatch[cent_bins];       //! pseudorapidity spectrum
  TH1F           * h1_nCells[cent_bins];       //! pseudorapidity spectrum
  TH1F           * h1_ClusterDisp[cent_bins];       //! cluster dispersion
  TH2F           * h2_Ellipse[cent_bins];       //! ellipse axis?
  TH2F           * h2_EtaPt[cent_bins];       //! 2d histogram Y - pseudorap spectrum
  TH3F           * h3_MptAsymm[cent_bins];       //!  3dimensional E vs mom
  TH3F           * h3_MptAsymm_mix[cent_bins];       //!  3dimensional E vs mom
  TH2F           * h2_dphi_deta[cent_bins];       //! 2dimensional E vs mom
  TH2F           * h2_dphi_deta_mix[cent_bins];       //! 2dimensional E vs mom
  TH2F           * h2_DispRes[cent_bins];       //! 2dimensional E vs mom
  TH2F           * h2_cells_M02[cent_bins];       //! 

  std::vector<TLorentzVector> Photons[poolDepth][zvtx_bins][mult_bins];
  std::vector<Int_t> TriggerList;
 
  AliAnalysisUtils*   fHelperClass;           //! Vertex selection helper

  AliAnalysisTaskEMCALMesonGGSDMpPb(const AliAnalysisTaskEMCALMesonGGSDMpPb&); // not implemented
  AliAnalysisTaskEMCALMesonGGSDMpPb& operator=(const AliAnalysisTaskEMCALMesonGGSDMpPb&); // not implemented
    
  ClassDef(AliAnalysisTaskEMCALMesonGGSDMpPb, 1); // example of analysis
};
#endif
