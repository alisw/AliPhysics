#ifndef AliAnalysisTaskHadronElectronAnalysis_H
#define AliAnalysisTaskHadronElectronAnalysis_H
class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliEventPoolManager;
class AliMultSelection;
class AliAnalysisUtils;
class AliGenEventHeader;
//All relevant AliPhysics includes (this list will continue to grow):
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH3D.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"


//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskHadronElectronAnalysis : public AliAnalysisTaskSE {

    
 public:
    
    enum EnhanceSigOrNot {kMB,kEnhance};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    
  AliAnalysisTaskHadronElectronAnalysis();
  AliAnalysisTaskHadronElectronAnalysis(const char *name);
  virtual ~AliAnalysisTaskHadronElectronAnalysis();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies();
  void SetEtaCut(float etaCut);
  void SetInvariantMassCut(double invariantMassCut) {fInvariantMassCut=invariantMassCut;};
  void SetPartnerElectronPtCut(double partnerElectronPtCut);
  void SetTPCNSigmaElecMin(double TPCNSigmaElecMin);
  void SetTPCNSigmaElecMax(double TPCNSigmaElecMax);
  void SetPartnerTPCNSigmaElecMin(double partnerTPCNSigmaElecMin);
  void SetPartnerTPCNSigmaElecMax(double partnerTPCNSigmaElecMax);
  void SetAssocTrkBitElec(double assocTrkBitElec);
  void SetAssocTrkBitH(double assocTrkBitH);
  void SetPartnerTrkBit(double partnerTrkBit);
  void SetTrigTrkBit(double trigTrkBit);
  void SetUseSPDKBoth(bool useSPDKBoth);
  void SetRunOnMC(bool runOnMC);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

 private:

 // Parameters that are set by setter functions
  float fEtaCut;
  double fInvariantMassCut;
  double fPartnerElectronPtCut;
  double fTPCNSigmaElecMin;
  double fTPCNSigmaElecMax;
  double fPartnerTPCNSigmaElecMin;
  double fPartnerTPCNSigmaElecMax;
  float fAssocTrkBitElec; 
  float fAssocTrkBitH;
  float fPartnerTrkBit;
  float fTrigTrkBit;
  bool fUseSPDKBoth;
  bool fRunOnMC;
  
// Old parameters (do what you will with them)
  float MULT_LOW;
  float MULT_HIGH;

  TString CENT_ESTIMATOR;
  TString EFF_FILE_PATH;
    AliVEvent* fVevent;
    TClonesArray* fMCarray;
    AliAODMCHeader* fMCHeader;
    Int_t fNTotMCpart;
    Int_t fNembMCpi0;
    Int_t fNembMCeta;
    Int_t  fNpureMC;
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  
 Bool_t fIsFrmEmbPi0;
 Bool_t fIsFrmEmbEta;
    Bool_t verbose;
    Bool_t fIsNew;
    Bool_t MC;
    Int_t ftype;//!
    Int_t fWeightPi0;
    Int_t fWeightEta;
    Int_t fWeight;
    TH1F *fRealInclsElecPt;
    TH1F *fNonHFeTrkPt;
    TH1F *fNonHFeEmbTrkPt;
    TH1F *fPi0eEmbWeightTrkPt;
    TH1F *fNonHFeEmbWeightTrkPt;
    TH1F *fRecoNonHFeTrkPt;
    TH1F *fEtaeEmbWeightTrkPt;
    TH1F *fRecoNonHFeEmbTrkPt;
    TH1F *fRecoPi0eEmbWeightTrkPt;
    TH1F *fRecoNonHFeEmbWeightTrkPt;
    TH1F *fRecoEtaeEmbWeightTrkPt;
    TF1  *fPi0Weight;//!
    TF1  *fEtaWeight;//!
    
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
    TH2D* fTriggerEffElec;
    TH2D* fAssociatedEffElec;
    TH2D* fTriggerEffH;
    TH2D* fAssociatedEffH;
    TH2D* fTriggEffVsPt;
    TH2D* fAssoEffVsPt;
    TH1D* fTrigPhiDist;
    TH1D* fAssHPhiDist;
    TH1D* fAssElecPhiDist;
   
    TH1D* fmultPercentile;
    TH1D* fmultPercentileCuts;

  TH2D* fDedx; //!>! dedx plot
  TH2D* fDedxCuts; //!>! dedx plot
  TH2D* fDedxTOF; //!>! dedx plot
  TH2D* fDedxTOFCuts; //!>! dedx plot
  TH2D* fBetaDedx;//!>! dedx plot
  TH2D* fBeta; //!>! dedx plot
  TH2D* fNSigma; //!>! dedx plot
  TH2D* fNSigmaCuts; //!>! dedx plot
    TH2D*  f2DPPtHist;
    TH2D* fNSigmaAssoCuts; //!>! dedx plot
    TH2D* fNSigmaAssoCutsPt; //!>! dedx plot
    TH2D* f2DHEpTDist;//!>! dedx plot
    TH2D* f2DHHpTDist;//!>! dedx plot
    TH2D* fHEDeltaEtaDeltaPhi;

    
    TH1D* fAllHadronsPt;
    TH1D* fAssHadronsPt;
    TH1D* fTriggerPt;
    TH1D* fElecPtBeforeCuts;
    TH1D* fElecPtAfterCuts;
    TH1D* fNumberTriggersPerEvent;
    TH1D*  fNumberElecsPerEvent;
    TH1D* fTOFHits;
    TH1D* fTPCCrossedRows;
    TH1D* fHEassNoTrig;
    TH1D* fHEnoAssTrig;
    TH1D* fHHassNoTrig;
    TH1D* fHHnoAssTrig;
    TH1D* fUnlikesignMassDist;
    TH1D* fLikesignMassDist;
  THnSparseF* fLooseDist;  //!>! single particle all hadron dist (no cuts at all)
  THnSparseF* fTriggerDist;  //!>! single particle trigger dist
  THnSparseF* fAssociatedHDist;  //!>! single particle associated hadron dist
  THnSparseF* fElectronDist;  //!>! single particle electron dist

  THnSparseF* fDphiHElec;  //!>! hadron-electron correlation hist
  THnSparseF* fDphiHUSElec;  //!>! hadron-electron correlation hist
  THnSparseF* fDphiHLSElec;  //!>! hadron-electron correlation hist
    THnSparseF* fDphiHUSElecEff;
    THnSparseF* fDphiHLSElecEff;
    THnSparseF* fDphiHElecEff;  //!>! hadron-electron correlation hist (efficiency corrected)bvg
  THnSparseF* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparseF* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparseF* fDphiTriggerTrigger;   //!>! trigger-trigger correlation hist
  THnSparseF* fDphiHElecMixed; //!>! hadron-electron mixed correlation hist
  THnSparseF* fDphiHULSElecMixed;
  THnSparseF* fDphiHLSElecMixed;
  THnSparseF* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparseF* fDphiTriggerTriggerMixed;   //!>! mixed trigger-trigger correlation hist

    THnSparseD* fSprsPi0EtaWeightCal;
  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection
    //!
  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, double multPercentile);
  void FillSingleParticleDist(std::vector<AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist, double multPercentile);
  void MakeSameHElecCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> corrElec_list, THnSparse* fDphi, double zVtx,double multPercentile, bool eff=true);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx,double multPercentile, bool eff=true);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff=true);
  void MakeMixedHElecCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> corrElec_list, THnSparse* fDphi, double zVtx,double multPercentile, bool eff=true);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, double multPercentile,bool eff=true);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCutsHadrons(AliAODTrack *track);
    bool PassAssociatedCutsElectrons(AliAODTrack *track);
    bool PassPartnerCutsElectrons(AliAODTrack *track);
  bool PassTPCCuts(AliAODTrack *track);
    void Make2DPtDistributionHistogram(std::vector<AliAODTrack*> trigger_list,std::vector<AliAODTrack*> corrElec_list, std::vector<AliAODTrack*> associated_h_list);

    std::vector<AliAODTrack*> GetUnlikeSignVector(std::vector<AliAODTrack*> corrElec_list,std::vector<AliAODTrack*> partnerElec_list);
    std::vector<AliAODTrack*> GetLikeSignVector(std::vector<AliAODTrack*> corrElec_list,std::vector<AliAODTrack*> partnerElec_list);

    bool GetNMCPartProducedpPbNew();
    bool GetNMCPartProducedpPbOld();
    void GetPi0EtaWeightpPbNew(THnSparse *SparseWeight);
    void GetPi0EtaWeightpPbOld(THnSparse *SparseWeight);
    Int_t GetPi0EtaType(AliAODMCParticle *part);
    Bool_t GetNonHFEEffiRecoTag(AliVTrack *track);
    Bool_t GetNonHFEEffiDenompPbNew(AliVTrack *track);
    Bool_t GetNonHFEEffiDenompPbOld(AliVTrack *track);
    Bool_t IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt);
  ClassDef(AliAnalysisTaskHadronElectronAnalysis, 3);

};
#endif
