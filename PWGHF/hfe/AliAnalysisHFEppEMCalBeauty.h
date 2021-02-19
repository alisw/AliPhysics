
///////////////////////////////////////////////////////////////////
//                                                               //            
// AliAnalysisHFEppEMCalBeauty.h                                 //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

#ifndef AliAnalysisHFEppEMCalBeauty_cxx
#define AliAnalysisHFEppEMCalBeauty_cxx

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;
class THnSparse;
class TClonesArray;
class TObjArray;


class AliAnalysisUtils;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEextraCuts;

class AliHFEpid;
class AliHFEpidQAmanager;

class AliEMCALTriggerPatchInfo;
class AliEMCALGeometry;
class AliMagF;

class AliAODEvent;
class AliVEvent;
class AliAODTrack;

class AliPIDResponse;
class AliAODMCHeader;
class AliSelectNonHFE;
class AliEventPoolManager;
class AliEventPool;
class AliGenEventHeader;


#include<iostream>
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliAODInputHandler.h"

#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliPID.h"
#include "AliHFEpid.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskPIDResponse.h"

#include "AliAODMCParticle.h"


class AliAnalysisHFEppEMCalBeauty : public AliAnalysisTaskSE {
 public:

  enum EnhanceSigOrNot {kMB,kEnhance};
  enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};//0,1,2,3,4,5
  enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};//0,1,2,3,4,5,6


AliAnalysisHFEppEMCalBeauty() : AliAnalysisTaskSE(),

fIsMC(kFALSE),
fIsAOD(kTRUE),
ftrigger(AliVEvent::kINT7),
 // emcal correction
fUseTender(kTRUE),
// flag for emcal dcal
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
// trigger events selection
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),

fRecalIP(kTRUE),

fEtarange(0.7),
fTPCNCrRows(70),
fRatioCrossedRowOverFindable(0.8),
fITSNclus(3),
fTPCNclusPID(60),
fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),
fDCAxyCut(1),
fDCAzCut(2),
fTPCnsigmin(-1),
fTPCnsigmax(3),
fCutEopEMin(0.9),
fCutEopEMax(1.2),
fM02Min(0.05),
fM02Max1(0.9),
fM02Max2(0.7),
fM02Max3(0.5),

fInvmassCut(0.14),
fAssoTPCCluster(60),
fAssoITSRefit(kTRUE),
fAssopTMin(0.1),
fAssoEtarange(0.9),
fAssoTPCnsig(3.0),


fTenderClusterName("caloClusters"),
fTenderTrackName("tracks"),
fTracks_tender(0),
fCaloClusters_tender(0),

 // events
fAOD(0),
fOutputList(0), 
fHistEvent(0),
fNentries(0),

fHistVx(0),
fHistVxwc(0),
fHistVy(0),
fHistVywc(0),
fHistVz(0),
fHistVzwc(0),

fHistMul(0),
fHistPt(0),
EtaPhiWoC(0),
EtaPhiWC(0),
EtaPhiAfTCATM(0),

fHistEta(0),
fHistPhi(0),
fHistEtaPhi_TPC(0),
EMCalEta_TPCpT(0),

fHistdca(0),
fHistdcaxy(0),
fHistdcaxywc(0),
fHistdcaz(0),
fHistdcazwc(0),
 //PID Cut
fPID(0),   
fPidResponse(0),
fHistBethe(0),
fnSigmaVsP_TPC(0),
fnSigmaVsP_TOF(0),


fHistClustE(0),
fEMCClsEtaPhi(0), 
fHistoNCells(0),
fHistoTimeEMC(0),

fHistPtMatch(0),
fEMCTrkMatch(0),
fEMCClsEtaPhiTrkMatch(0),
fEMCTrkMatch_Phi(0),
fEMCTrkMatch_Eta(0),

fvalueElectron(0),
fSparseElectron(0),

//MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
fPthfeGenerated(0),
fPthfe_rec(0),
fPthfe_rec_TrkSel(0),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),

fMCparticleMother(0),
fMCparticleGMother(0),
fMCparticleGGMother(0),
fMCparticleGGGMother(0),

//Used in the function FindMother
fIsHFE1(kFALSE),
fIsHFE2(kFALSE),
fIsNonHFE(kFALSE),
fIsFromD(kFALSE),
fIsFromBarionB(kFALSE),
fIsFromMesonB(kFALSE),
fIsFromBarionBD(kFALSE),
fIsFromMesonBD(kFALSE),
fIsFromPi0(kFALSE),
fIsFromEta(kFALSE),
fIsFromGamma(kFALSE),

//EID Cuts
fTrkDCA(-999.0),

fEop_AftEID(0),
fEopNL_AftEID(0),

fHadEovp_AftEID(0),
fHadEovpNL_AftEID(0),
fHadPt_AftEID(0),

fHadDCA(0),
fInclsElecPt(0),
fInclElecDCA(0),

fNElecInEvt(0),
fNEle(0),

fTPCnSigma(-999.0),

fTPCnSigmaHadMin(-10),
fTPCnSigmaHadMax(-3.5),

fInvmassULSPt(0),
fInvmassLSPt(0),
fCalculateNonHFEEffi(1),

fULSElecPt(0),
fLSElecPt(0),
fULSElecDCA(0),
fLSElecDCA(0),

//nonhfe efficiency

fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),

fPi0Weight(0),
fEtaWeight(0),

fCalculateWeight(kFALSE), 
fSprsPi0EtaWeightCal(0),
fPi0EtaSpectraSp(0),
pi0MC(0),
etaMC(0),
gammaMC(0),

fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),

fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),

fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
fNonHFeEmbInvmassLS(0),
fNonHFeEmbInvmassULS(0),
fNonHFeEmbWeightInvmassLS(0),
fNonHFeEmbWeightInvmassULS(0),
fPi0EmbInvmassLS(0),
fPi0EmbInvmassULS(0),
fPi0EmbWeightInvmassLS(0),
fPi0EmbWeightInvmassULS(0),
fEtaEmbInvmassLS(0),
fEtaEmbInvmassULS(0),
fEtaEmbWeightInvmassLS(0),
fEtaEmbWeightInvmassULS(0),

fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0),

fHadConvRadius(0),
fIncleConvRadius(0),
fNonHFeConvRadius(0),
fHFeConvRadius(0),

fNonHFeEmbTrkRConv(0),
fPi0eEmbWeightTrkRConv(0),
fNonHFeEmbWeightTrkRConv(0),
fEtaeEmbWeightTrkRConv(0),

fRecoNonHFeEmbRConv(0),
fRecoPi0eEmbWeightTrkRConv(0),
fRecoNonHFeEmbWeightTrkRConv(0),
fRecoEtaeEmbWeightTrkRConv(0),

fRVsULSElecPt(0),
fRVsLSElecPt(0)

{
fPID = new AliHFEpid("hfePid");
fvalueElectron = new Double_t[6];
//fvalueRadius = new Double_t[4];
}
  
   
  AliAnalysisHFEppEMCalBeauty(const char *name);
  virtual ~AliAnalysisHFEppEMCalBeauty(); 
  
  virtual void   Init();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

//-----------------selections cuts----------------------------------------------------------------  

  void SetMCAnalysis(Bool_t isMC){fIsMC=isMC;}
  void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};

  void SetTrigger(AliVEvent::EOfflineTriggerTypes trigger){ftrigger =trigger;}

  void SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
  void SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
  void SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};

  void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
  void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
  void SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;};
  void SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;};

  void SwitchRecalImpPar(Bool_t fSwitchRIP) {fRecalIP = fSwitchRIP;};

  //----------Setter for Track and PID cuts
 
  void SetEtaRange(Double_t Etarange){fEtarange=Etarange;}
  void SetMinTPCCluster(Int_t TPCNCrRows){fTPCNCrRows=TPCNCrRows;}
  void SetMinRatioCrossedRowOverFindable(Double_t RatioCrossedRowOverFindable){fRatioCrossedRowOverFindable=RatioCrossedRowOverFindable;}
  void SetMinITSCluster(Int_t ITSNclus){fITSNclus=ITSNclus;}
  void SetMinTPCClusterPID(Int_t TPCNclusPID){fTPCNclusPID=TPCNclusPID;}
  void SetHitsOnSPDLayers(Bool_t SPDBoth,Bool_t SPDAny, Bool_t SPDFirst){fSPDBoth=SPDBoth; if(!SPDBoth)fSPDAny=SPDAny; if(!SPDBoth && !SPDAny)fSPDFirst=SPDFirst;}
  void SetDCACut(Double_t DCAxyCut,Double_t DCAzCut){ fDCAxyCut=DCAxyCut; fDCAzCut=DCAzCut;}
  void SetTPCnsigma(Double_t TPCnsigmin,Double_t TPCnsigmax){fTPCnsigmin=TPCnsigmin;  fTPCnsigmax=TPCnsigmax;}
  void SetEopE(Double_t EopEMin,Double_t EopEMax){ fCutEopEMin=EopEMin; fCutEopEMax=EopEMax;}
  void SetShowerShapeEM02(Double_t M02Min, Double_t M02Max1, Double_t M02Max2, Double_t M02Max3) {fM02Min = M02Min; fM02Max1 = M02Max1; fM02Max2 = M02Max2; fM02Max3 = M02Max3;};
   
  //------------Setters for Photonic Electron Selection Cuts
  void SetInvMassCut(Double_t InvmassCut){fInvmassCut=InvmassCut;}
  void SetAssoTPCclus(Int_t AssoTPCCluster){fAssoTPCCluster=AssoTPCCluster;}
  void SetAssoITSrefit(Bool_t AssoITSRefit){fAssoITSRefit= AssoITSRefit;}
  void SetAssopTMin(Double_t AssopTMin){fAssopTMin = AssopTMin;}
  void SetAssoEtarange(Double_t AssoEtarange){fAssoEtarange=AssoEtarange;}
  void SetAssoTPCnsig(Double_t AssoTPCnsig){fAssoTPCnsig=AssoTPCnsig;}

//-----------------getters----------------------------------------------------------------

  Bool_t     GetTenderSwitch()    {return fUseTender;};
  Bool_t     GetEMCalTriggerEG1() { return fEMCEG1; };
  Bool_t     GetEMCalTriggerEG2() { return fEMCEG2; };
  Bool_t     GetEMCalTriggerDG1() { return fDCalDG1; };
  Bool_t     GetEMCalTriggerDG2() { return fDCalDG2; };

  Int_t ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx);                                                   // track selection cuts

  Bool_t  PassEIDCuts(AliAODTrack *track, AliAODCaloCluster *clust, Bool_t &Hadtrack);

  void GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff);   //trk clus matching     

  Bool_t  GetNonHFEEffiDenom(AliAODTrack *track);
  Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
  Int_t   GetPi0EtaType(AliAODMCParticle *part);
  Bool_t  GetNonHFEEffiRecoTag(AliAODTrack *track);
  Bool_t  GetNonHFEEffiULSLS(AliAODTrack *track, AliAODTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
  
  Int_t GetHFE(AliAODMCParticle *, TClonesArray *);
  Int_t GetElecSourceType(AliAODMCParticle *,Double_t &ptm);

  void    SetNonHFEEffi(Bool_t fSwitch) {fCalculateNonHFEEffi = fSwitch;};
  void    GetPi0EtaWeight(THnSparse *SparseWeight);
  void    SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC);		

  void TrackConvRadius(AliAODTrack* track, Double_t &R);
 
  void    RecalImpactParam(const AliAODTrack * const track, Double_t dz[2], Double_t covar[3]);
  AliAODVertex*   RemoveDaughtersFromPrimaryVtx(const AliAODTrack * const track);

  AliHFEpid *GetPID() const {return fPID;};
//______________________________________________________________________
    

 private:
  
  Bool_t        fIsMC;
  Bool_t        fIsAOD;         //!flag for AOD analysis

  AliVEvent::EOfflineTriggerTypes ftrigger;

  Bool_t fUseTender;// switch to add tender
  Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
  Bool_t fFlagClsTypeDCAL;//switch to select DCAL c

  Bool_t fEMCEG1;//EMcal Threshold EG1
  Bool_t fEMCEG2;//EMcal Threshold SetReferenceMultiplicityEG2
  Bool_t fDCalDG1;//DCal Threshold DG1
  Bool_t fDCalDG2;//DCal Threshold DG2

  Bool_t fRecalIP;// 

  //------------Track and PID cut variables--------------
  Double_t fEtarange; 
  Int_t fTPCNCrRows; 
  Double_t fRatioCrossedRowOverFindable; 
  Int_t fITSNclus;  
  Int_t fTPCNclusPID;  
  Bool_t fSPDBoth;  
  Bool_t fSPDAny;  
  Bool_t fSPDFirst;  
  Double_t fDCAxyCut;  
  Double_t fDCAzCut;    

  Double_t fTPCnSigma;//!
  
  Double_t fTPCnsigmin;  
  Double_t fTPCnsigmax;  
  Double_t fCutEopEMin;
  Double_t fCutEopEMax;
  Double_t fM02Min;//!
  Double_t fM02Max1;//!
  Double_t fM02Max2;//!
  Double_t fM02Max3;//!
  
  Double_t fInvmassCut; //    invariant mass cut value
  Int_t fAssoTPCCluster;  
  Bool_t fAssoITSRefit;  
  Double_t fAssopTMin;  
  Double_t fAssoEtarange;  
  Double_t fAssoTPCnsig;  
  

  TString fTenderClusterName;//
  TString fTenderTrackName;//
  TClonesArray* fTracks_tender;//Tender tracks
  TClonesArray* fCaloClusters_tender;//Tender cluster    



  Bool_t GetNMCPartProduced();
  Bool_t FindMother(Int_t mcIndex);
  Bool_t IsHFelectronsMC(AliAODTrack *track);

 

  TList       *fOutputList;//! Output list
  AliAODEvent *fAOD;//! AOD object

  AliHFEpid   *fPID;//! 
  AliPIDResponse   *fPidResponse;//!pid response 

  Double_t            fTPCnSigmaHadMin;//!
  Double_t            fTPCnSigmaHadMax;//!

  TH1F        *fHistEvent;//! 
  TH1F        *fNentries;//! 

  TH1F        *fHistVx;//! 
  TH1F        *fHistVy;//! 
  TH1F        *fHistVz;//!  
  
  TH1F        *fHistVxwc;//!
  TH1F        *fHistVywc;//!
  TH1F        *fHistVzwc;//! 

  TH1F        *fHistMul;//! 

  TH1F        *fHistPt;//! Pt spectrum

  TH2F        *EtaPhiWoC;//!
  TH2F        *EtaPhiWC;//!
  TH2F        *EtaPhiAfTCATM;//!

  TH1F        *fHistEta;//! 
  TH1F        *fHistPhi;//!  
  TH2F        *fHistEtaPhi_TPC;//! 
  TH1F        *EMCalEta_TPCpT;//!  

  TH2F        *fHistdca;//! 
  TH2F        *fHistdcaxy;//!      
  TH2F        *fHistdcaxywc;//! 
  TH2F        *fHistdcaz;//! 
  TH2F        *fHistdcazwc;//! 

  
  TH2F  *fHistBethe;//! 
  TH2F  *fnSigmaVsP_TPC;//!
  TH2F  *fnSigmaVsP_TOF;//!

  TH1F  *fHistClustE;//!
  TH2F  *fEMCClsEtaPhi;//!
  TH2F  *fHistoNCells;//!
  TH2F  *fHistoTimeEMC;//!


  TH1F  *fHistPtMatch;//!
  TH2F  *fEMCTrkMatch;//!
  TH2F  *fEMCClsEtaPhiTrkMatch;//!
  TH2F  *fEMCTrkMatch_Phi;//!      
  TH2F  *fEMCTrkMatch_Eta;//! 

  Double_t* fvalueElectron;//!
  THnSparse* fSparseElectron;//! Electron information

  TH2F   *fRVsULSElecPt;//!
  TH2F   *fRVsLSElecPt;//!

  TH2F   *fHadConvRadius;//!
  TH2F   *fIncleConvRadius;//!
  TH2F   *fNonHFeConvRadius;//!
  TH2F   *fHFeConvRadius;//!

  TH2F   *fNonHFeEmbTrkRConv;//!
  TH2F   *fPi0eEmbWeightTrkRConv;//!
  TH2F   *fNonHFeEmbWeightTrkRConv;//!
  TH2F   *fEtaeEmbWeightTrkRConv;//!

  TH2F   *fRecoNonHFeEmbRConv;//!
  TH2F   *fRecoPi0eEmbWeightTrkRConv;//!
  TH2F   *fRecoNonHFeEmbWeightTrkRConv;//!
  TH2F   *fRecoEtaeEmbWeightTrkRConv;//!


  Int_t     fNEle;//!
  Double_t  fTrkDCA;//!
 
  TH1F  *fHadPt_AftEID;//!

  TH2F  *fHadEovp_AftEID;//!
  TH2F  *fHadEovpNL_AftEID;//!

  TH2F  *fHadDCA;//!
  TH1F  *fInclsElecPt;//!  
  TH2F  *fInclElecDCA;//!

  TH2F  *fEop_AftEID;//!
  TH2F  *fEopNL_AftEID;//!

  TH1F  *fNElecInEvt;//!

  TH2F  *fInvmassULSPt;//!
  TH2F  *fInvmassLSPt;//!

  TH1F  *fULSElecPt;//!
  TH1F  *fLSElecPt;//!

  TH2F  *fULSElecDCA;//! 
  TH2F  *fLSElecDCA;//!  
   
 //Used in the function FindMother
  Bool_t        fIsHFE1;//!
  Bool_t        fIsHFE2;//!
  Bool_t        fIsNonHFE;//!
  Bool_t        fIsFromD;//!
  Bool_t        fIsFromBarionB;//!
  Bool_t        fIsFromMesonB;//!
  Bool_t        fIsFromBarionBD;//!
  Bool_t        fIsFromMesonBD;//!
  Bool_t        fIsFromPi0;//!
  Bool_t        fIsFromEta;//!
  Bool_t        fIsFromGamma;//!


  AliAODMCHeader    *fMCHeader;//!
  AliAODMCParticle  *fMCparticle;//!
  TClonesArray      *fMCArray;//

  AliAODMCParticle  *fMCparticleMother;//!
  AliAODMCParticle  *fMCparticleGMother;//!
  AliAODMCParticle  *fMCparticleGGMother;//!
  AliAODMCParticle  *fMCparticleGGGMother;//!


  Int_t               fNTotMCpart;
  Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
  Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
  Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator

  TH1F* fPthfeGenerated;//!
  TH1F* fPthfe_rec;//!
  TH1F* fPthfe_rec_TrkSel;//

    //non hfe

  Bool_t     fIsFrmEmbPi0;//!
  Bool_t     fIsFrmEmbEta;//!
  Int_t      ftype;//!
  Double_t   fWeight;//!

  Double_t   fWeightPi0;//!
  Double_t   fWeightEta;//!
  
  TF1       *fPi0Weight;//!
  TF1       *fEtaWeight;//!

  Bool_t    fCalculateWeight;//!
  THnSparse *fSprsPi0EtaWeightCal;//!
  THnSparseF  *fPi0EtaSpectraSp;//!
 
  TH1F        *pi0MC;//!
  TH1F        *etaMC;//!
  TH1F        *gammaMC;//!

  TH1F      *fRealInclsElecPt;//!
  TH1F      *fNonHFeTrkPt;//!
  TH1F      *fNonHFeEmbTrkPt;//!
  TH1F      *fNonHFeEmbWeightTrkPt;//!
  TH1F      *fPi0eEmbWeightTrkPt;//!
  TH1F      *fEtaeEmbWeightTrkPt;//!
  
  TH1F      *fRecoNonHFeTrkPt;//!
  TH1F      *fRecoNonHFeEmbTrkPt;//!
  TH1F      *fRecoNonHFeEmbWeightTrkPt;//!
  TH1F      *fRecoPi0eEmbWeightTrkPt;//!
  TH1F      *fRecoEtaeEmbWeightTrkPt;//!

 TH1F       *fNonHFePairInvmassLS;//!
 TH1F       *fNonHFePairInvmassULS;//!
 TH1F       *fNonHFeEmbInvmassLS;//!
 TH1F       *fNonHFeEmbInvmassULS;//!
 TH1F       *fNonHFeEmbWeightInvmassLS;//!
 TH1F       *fNonHFeEmbWeightInvmassULS;//!
 TH1F       *fPi0EmbInvmassLS;//!
 TH1F       *fPi0EmbInvmassULS;//!
 TH1F       *fPi0EmbWeightInvmassLS;//!
 TH1F       *fPi0EmbWeightInvmassULS;//!
 TH1F       *fEtaEmbInvmassLS;//!
 TH1F       *fEtaEmbInvmassULS;//!
 TH1F       *fEtaEmbWeightInvmassLS;//!
 TH1F       *fEtaEmbWeightInvmassULS;//!

 TH1F       *fRecoLSeEmbTrkPt;//!
 TH1F       *fRecoLSeEmbWeightTrkPt;//!
 TH1F       *fRecoPi0LSeEmbWeightTrkPt;//!
 TH1F       *fRecoEtaLSeEmbWeightTrkPt;//!
 TH1F       *fRecoULSeEmbTrkPt;//!
 TH1F       *fRecoULSeEmbWeightTrkPt;//!
 TH1F       *fRecoPi0ULSeEmbWeightTrkPt;//!
 TH1F       *fRecoEtaULSeEmbWeightTrkPt;//!

  Bool_t    fCalculateNonHFEEffi;//!


  AliAnalysisHFEppEMCalBeauty(const AliAnalysisHFEppEMCalBeauty&); // not implemented
  AliAnalysisHFEppEMCalBeauty& operator=(const AliAnalysisHFEppEMCalBeauty&); // not implemented
  
  ClassDef(AliAnalysisHFEppEMCalBeauty, 1); // example of analysis
};

#endif
