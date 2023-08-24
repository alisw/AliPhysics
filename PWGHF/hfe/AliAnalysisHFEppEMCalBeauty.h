
///////////////////////////////////////////////////////////////////
//                                                               //            
// AliAnalysisHFEppEMCalBeauty.h                                 //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

#ifndef AliAnalysisHFEppEMCalBeauty_H
#define AliAnalysisHFEppEMCalBeauty_H

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

#include "AliHFEextraCuts.h"
#include "AliHFEmcQA.h"

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

class AliHFEextraCuts;
class AliHFEmcQA;


class AliAnalysisHFEppEMCalBeauty : public AliAnalysisTaskSE {
 public:

  AliAnalysisHFEppEMCalBeauty();
  AliAnalysisHFEppEMCalBeauty(const char *name);
  virtual ~AliAnalysisHFEppEMCalBeauty(); 
  
  virtual void   Init();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);

  enum EnhanceSigOrNot {kMB,kEnhance};
  enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};//0,1,2,3,4,5
  enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};//0,1,2,3,4,5,6

	enum SourceType{
	  kDirectCharm=1,     // electrons from primary charmed hadrons 
	  kDirectBeauty=2,    // electrons from primary beauty hadrons
	  kBeautyCharm=3,     // electrons from charmed hadrons decaying from the beauty hadrons (B->D->e)
	  kGamma=4,           // should be obsolete -> please let me know if you see something! 
	  kPi0=5,             // electrons from pi0 Dalitz
	  kEta=6,             // electrons from eta Dalitz
	  kOmega=7,           // electrons from omega decay (Dalitz and di-electrons)
	  kPhi=8,            // electrons from phi decay (di-electron)
	  kEtaPrime=9,       // electrons from eta prime decay (Dalitz and 2charged-pions&di-electrons) 
	  kRho0=10,           // electrons from rho decay (di-electron)
	  kK0s2P=11,          // electrons from secondary pions from K0s
	  kK0l2P=12,          // electrons from secondary pions from K0l
	  kLamda2P=13,        // electrons from secondary pions from Lamda
	  kSigma2P=14,        // electrons from secondary pions from Sigma
	  kK2P=15,            // electrons from secondary pions from K  
	  kElse=16,            // all the other sources which was not in this enumeration
	  kMisID=17,           // not the electrons (hadrons)
	  kGammaPi0=18,       // electrons from photon conversion where the photon originated from pi0
	  kGammaEta=19,       // electrons from photon conversion where the photon originated from eta
	  kGammaOmega=20,     // electrons from photon conversion where the photon originated from omega
	  kGammaPhi=21,       // electrons from photon conversion where the photon originated from phi
	  kGammaEtaPrime=22,  // electrons from photon conversion where the photon originated from eta prime
	  kGammaRho0=23,      // electrons from photon conversion where the photon originated from rho
	  kGammaK0s2P=24,     // electrons from photon conversion where the photon originated from secondary pion from K0s
	  kGammaK0l2P=25,     // electrons from photon conversion where the photon originated from secondary pion from K0l
	  kGammaLamda2P=26,   // electrons from photon conversion where the photon originated from secondary pion from lamda
	  kGammaSigma2P=27,   // electrons from photon conversion where the photon originated from secondary pion from sigma
	  kGammaK2P=28,        // electrons from photon conversion where the photon originated from secondary pion from K
	  kJpsi=29,           // electrons from primary J/psi decay
	  kB2Jpsi=30,         // electrons from J/psi decay where the J/psi originated from the beauty hadrons
	  kKe3=31,            // Ke3 electrons
	  kK0L=38,            // K0L -> e +X
	  kPromptD0=39,
	  kNonPromptD=40,
	  kPromptB=41,
	  kPromptLc=42,
	  //kDirectGamma=43,
	  kDalitzGamma=44,
	  kBaryonGamma=45
	};

	// Setter for B corr (central)
	void SetBcorrCentLow(TF1* BcorrFcentL) { fBmesonCorrCentLow = BcorrFcentL; }
	void SetBcorrCentHigh(TF1* BcorrFcentH) { fBmesonCorrCentHigh = BcorrFcentH; }
	void SetBcorrMinLow(TF1* BcorrFminL) { fBmesonCorrMinLow = BcorrFminL; }
	void SetBcorrMinHigh(TF1* BcorrFminH) { fBmesonCorrMinHigh = BcorrFminH; }
	void SetBcorrMaxLow(TF1* BcorrFmaxL) { fBmesonCorrMaxLow = BcorrFmaxL; }
	void SetBcorrMaxHigh(TF1* BcorrFmaxH) { fBmesonCorrMaxHigh = BcorrFmaxH; }

	// Setter for D corr (central)
	void SetDcorrFtn(TF1* DcorrF) { fDmesonCorr = DcorrF; }
	void SetDcorrFtnVar1(TF1* DcorrFVar1) { fDmesonCorrVar1 = DcorrFVar1; }
	void SetDcorrFtnVar2(TF1* DcorrFVar2) { fDmesonCorrVar2 = DcorrFVar2; }

	// Setter for Lc corr (central)
	void SetLccorrFtn(TF1* LccorrF) {fLcCorr = LccorrF;}

//-----------------selections cuts----------------------------------------------------------------  
  void SetMCAnalysis(Bool_t isMC){fIsMC=isMC;}
  void SetAODanalysis(Bool_t IsAOD) {fIsAOD=IsAOD;}
  void SetTrigger(AliVEvent::EOfflineTriggerTypes trigger){ftrigger=trigger;}

  void SetPi0EtaWeightCalc(Bool_t IsPi0EtaWeight){fCalcPi0EtaWeight=IsPi0EtaWeight;}
  void SetNonHFEEffi(Bool_t SwitchNHFEeffi){fCalculateNonHFEEffi=SwitchNHFEeffi;}
  void SetElecRecoEffi(Bool_t SwitchEleRecoEffi){fCalculateElecRecoEffi=SwitchEleRecoEffi;}
  void SwitchMCTemplateWeightCalc(Bool_t SwitchMCTempWeight){fCalculateMCTemplWeightCalc=SwitchMCTempWeight;}
  void SwitchFillMCTemplate(Bool_t SwitchFillTemp){fFillMCTemplates=SwitchFillTemp;} 
  
  void SetTenderSwitch(Bool_t usetender){fUseTender=usetender;}
  void SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC=flagClsEMC;}
  void SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL=flagClsDCAL;}

  void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;}
  void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;}
  void SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;}
  void SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;}

  //----------Setter for Track and PID cuts
 
  void SetEtaRange(Double_t Etarange){fEtarange=Etarange;}
  void SetMinTPCCluster(Int_t TPCNCrRows){fTPCNCrRows=TPCNCrRows;}
  void SetMinRatioCrossedRowOverFindable(Double_t RatioCrossedRowOverFindable){fRatioCrossedRowOverFindable=RatioCrossedRowOverFindable;}
  void SetMinITSCluster(Int_t ITSNclus){fITSNclus=ITSNclus;}
  void SetMinTPCClusterPID(Int_t TPCNclusPID){fTPCNclusPID=TPCNclusPID;}
  void SetHitsOnSPDLayers(Bool_t SPDBoth,Bool_t SPDAny, Bool_t SPDFirst){fSPDBoth=SPDBoth; if(!SPDBoth)fSPDAny=SPDAny; if(!SPDBoth && !SPDAny)fSPDFirst=SPDFirst;}
  void SetEMCalMatching(Double_t DeltaEta,Double_t DeltaPhi){fCutDeltaEta=DeltaEta; fCutDeltaPhi=DeltaPhi;} //trk clus matching cut
  void SetDCACut(Double_t DCAxyCut,Double_t DCAzCut){ fDCAxyCut=DCAxyCut; fDCAzCut=DCAzCut;}
  void SetTPCnsigma(Double_t TPCnsigmin,Double_t TPCnsigmax){fTPCnsigmin=TPCnsigmin;  fTPCnsigmax=TPCnsigmax;}
  void SetEopE(Double_t EopEMin,Double_t EopEMax){ fCutEopEMin=EopEMin; fCutEopEMax=EopEMax;}
  void SetShowerShapeEM02(Double_t M02Min, Double_t M02Max1, Double_t M02Max2, Double_t M02Max3) {fM02Min = M02Min; fM02Max1 = M02Max1; fM02Max2 = M02Max2; fM02Max3 = M02Max3;}
   
  //------------Setters for Photonic Electron Selection Cuts
  void SetInvMassCut(Double_t InvmassCut){fInvmassCut=InvmassCut;}
  void SetAssoTPCclus(Int_t AssoTPCCluster){fAssoTPCCluster=AssoTPCCluster;}
  void SetAssoITSrefit(Bool_t AssoITSRefit){fAssoITSRefit=AssoITSRefit;}
  void SetAssopTMin(Double_t AssopTMin){fAssopTMin=AssopTMin;}
  void SetAssoEtarange(Double_t AssoEtarange){fAssoEtarange=AssoEtarange;}
  void SetAssoTPCnsig(Double_t AssoTPCnsig){fAssoTPCnsig=AssoTPCnsig;}


  void SwitchRecalImpPar(Bool_t fSwitchRIP){fRecalIP=fSwitchRIP;}

  void  SetDmesonWeightHist(TH1 *D1, TH1 *D2, TH1 *D3);
  void  SetBmesonWeightHist(TH1 *B1, TH1 *B2, TH1 *B3);
  void  SetDmesonWeightHistPbPb(TH1 *D0, TH1 *DPlus, TH1 *Ds, TH1 *Lc);

//-----------------getters----------------------------------------------------------------
  Bool_t  GetTenderSwitch()    {return fUseTender;};
  Bool_t  GetEMCalTriggerEG1() { return fEMCEG1; };
  Bool_t  GetEMCalTriggerEG2() { return fEMCEG2; };
  Bool_t  GetEMCalTriggerDG1() { return fDCalDG1; };
  Bool_t  GetEMCalTriggerDG2() { return fDCalDG2; };

  Int_t ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx); // track selection cuts

  Bool_t  PassEIDCuts(AliAODTrack *track, AliAODCaloCluster *clust, Bool_t &Hadtrack);

  void GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff);   //trk clus eta-phi difference     

  Bool_t  GetNonHFEEffiDenom(AliAODTrack *track);
  Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
  Int_t   GetPi0EtaType(AliAODMCParticle *part);
  Bool_t  GetNonHFEEffiRecoTag(AliAODTrack *track);
  Bool_t  GetNonHFEEffiULSLS(AliAODTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
  
  Int_t GetHFE(AliAODMCParticle *, TClonesArray *);
  Int_t GetElecSourceType(AliAODMCParticle *,Double_t &ptm);

  void  GetPi0EtaWeight(THnSparse *SparseWeight);
  void  SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC);		

  void TrackConvRadius(AliAODTrack* track, Double_t &R);
 
  void RecalImpactParam(const AliAODTrack * const track, Double_t dz[2], Double_t covar[3]);
  AliAODVertex*   RemoveDaughtersFromPrimaryVtx(const AliAODTrack * const track);

  AliHFEpid *GetPID() const {return fPID;};
//______________________________________________________________________
  void    GetMCTemplateWeight();
  Bool_t  GetMCDCATemplates(AliAODTrack *track, Double_t TrkDCA);
  void    GetBWeight(AliAODMCParticle *Part, Double_t &BCentWeight, Double_t &BMinWeight, Double_t &BMaxWeight);
  void    GetDWeight(AliAODMCParticle *Part, Double_t &DCentWeight, Double_t &DMinWeight, Double_t &DMaxWeight);
  void    GetDWeightPbPb(AliAODMCParticle *Part, Int_t PDG, Double_t &DCentWeight);
//______________________________________________________________________
  void    GetElectronFromStack();    
  void    GetTrackHFStatus(AliAODTrack *track, Bool_t &IsMCEle, Bool_t &IsMCPPEle, Bool_t &IsMCHFEle, Bool_t &IsMCBEle, Bool_t &IsMCDEle);
  void    GetEIDRecoEffi(AliAODTrack *track, AliAODCaloCluster *clust, Bool_t IsMCPPEle, Bool_t IsMCHFEle, Bool_t IsMCBEle, Bool_t IsMCDEle, Double_t fTPCnSigma);

	int GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg);
	int GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &hfpt, double &hfeta);

 private:
  
  Bool_t fIsMC;//
  Bool_t fIsAOD;//   

  Bool_t fCalcPi0EtaWeight;//
  Bool_t fCalculateNonHFEEffi;// cant be made kFALSE for IsMC
  Bool_t fCalculateElecRecoEffi;//
  Bool_t fCalculateMCTemplWeightCalc;//
  Bool_t fFillMCTemplates;//

  Bool_t fRecalIP;// 

  AliVEvent::EOfflineTriggerTypes ftrigger;

  Bool_t fUseTender;// switch to add tender
  Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
  Bool_t fFlagClsTypeDCAL;//switch to select DCAL c

  Bool_t fEMCEG1;//EMcal Threshold EG1
  Bool_t fEMCEG2;//EMcal Threshold SetReferenceMultiplicityEG2
  Bool_t fDCalDG1;//DCal Threshold DG1
  Bool_t fDCalDG2;//DCal Threshold DG2

  //------------Track and PID cut variables--------------
  Double_t fEtarange;// 
  Int_t fTPCNCrRows;//
  Double_t fRatioCrossedRowOverFindable;// 
  Int_t fITSNclus;//  
  Int_t fTPCNclusPID; // 
  Bool_t fSPDBoth;//  
  Bool_t fSPDAny;//  
  Bool_t fSPDFirst;//  
  Double_t fDCAxyCut;//  
  Double_t fDCAzCut;//    
  Double_t fCutDeltaEta;//
  Double_t fCutDeltaPhi;//

  Double_t fTPCnSigma;//
  Double_t fTPCnsigmin;//  
  Double_t fTPCnsigmax;//  
  Double_t fCutEopEMin;//
  Double_t fCutEopEMax;//
  Double_t fM02Min;//
  Double_t fM02Max1;//
  Double_t fM02Max2;//
  Double_t fM02Max3;//
  
  Double_t fInvmassCut;// 
  Int_t fAssoTPCCluster;//  
  Bool_t fAssoITSRefit;//  
  Double_t fAssopTMin;//  
  Double_t fAssoEtarange;//  
  Double_t fAssoTPCnsig;//  
  
  TString fTenderClusterName;//!
  TString fTenderTrackName;//!
  TClonesArray* fTracks_tender;//!Tender tracks
  TClonesArray* fCaloClusters_tender;//!Tender cluster    

  Bool_t GetNMCPartProduced();
  Bool_t FindMother(Int_t mcIndex);
  Bool_t IsHFelectronsMC(AliAODTrack *track);
  
  TList       *fOutputList;//! Output list
  AliAODEvent *fAOD;//! AOD object
  AliHFEpid   *fPID;//! 
  AliPIDResponse   *fPidResponse;//!

//______________________AliHFEClass________________________________________________

  AliHFEextraCuts *fExtraCuts;//! 
  AliHFEmcQA      *fHFEmcQA; //!

	// generated level
	TH1F	*hBhadronPt;
	TH1F	*hBhadronPtCorr;
	TH1F	*hD0Pt;
	TH1F	*hD0PtCorr;
	TH1F	*hLcPt;
  //Reco effi
	TH1F	*hGenBePt_All;
	TH1F	*hRecBePt_TrkCuts;
	TH1F	*hRecBePt_EMCMatch;
	TH1F	*hRecBePt_EovP;
	TH1F	*hRecBePt_TPCnsig;
	TH1F	*hRecBePt_SS;
	// dca
	TH2F	*dcaTrack;
	TH2F	*dcaPion;
	TH2F	*dcaBeauty;
	TH2F	*dcaBeautyCorr;
	TH2F	*dcaBeautyCorrVar1;
	TH2F	*dcaBeautyCorrVar2;
	TH2F	*dcaBzero;
	TH2F	*dcaBplus;
	TH2F	*dcaBszero;
	TH2F	*dcaLb;
	TH2F	*DelecVsDmother;
	TH2F	*dcaDmeson;
	TH2F	*dcaDmesonCorr;
	TH2F	*dcaDmesonCorrVar1;
	TH2F	*dcaDmesonCorrVar2;
	TH2F	*dcaDzero;
	TH2F	*dcaDplus;
	TH2F	*dcaDsplus;
	TH2F	*dcaLc;
	TH2F	*dcaDalitz;
	TH2F	*dcaConv;	

	// B corr
  TF1   *fBmesonCorrCentLow;
  TF1   *fBmesonCorrCentHigh;
  TF1   *fBmesonCorrMinLow;
  TF1   *fBmesonCorrMinHigh;
  TF1   *fBmesonCorrMaxLow;
  TF1   *fBmesonCorrMaxHigh;

  // D corr
  TF1   *fDmesonCorr;
  TF1   *fDmesonCorrVar1;
  TF1   *fDmesonCorrVar2;

  // Lc corr
  TF1       *fLcCorr;
  TRandom3  *fRnd;

  	//===============NonHFE========================================
	 AliSelectNonHFE *fNonHFE;//!
   TH1F   *fPte_ULS; //! ULS elec Pt
   TH1F   *fPte_LS;//! LS elec pt  
   TH1F   *fInvmassLS1; //! LS Invmass 
   TH1F   *fInvmassULS1;//! ULS Invmass
   TH2F   *fDCAULSElec; //! ULS elec Pt
   TH2F   *fDCALSElec;//! LS elec pt 


  Double_t    fTPCnSigmaHadMin;//!
  Double_t    fTPCnSigmaHadMax;//!

  TH1F  *fHistEvent;//! 
  TH1F  *fNentries;//! 

  TH1F  *fHistVz;//!  
  TH1F  *fHistVzwc;//! 
  TH1F  *fHistPt;//! 
  TH1F  *fHistEta;//! 
  TH1F  *fHistPhi;//!  
  TH2F  *fHistdca;//!  
  TH2F  *fHistBethe;//! 
  TH2F  *fnSigmaVsP_TPC;//!

  TH1F  *fHistClustE;//!
  TH2F  *fEMCClsEtaPhi;//!
  TH2F  *fHistoNCells;//!
  TH2F  *fHistoTimeEMC;//!

  TH1F  *fHistPtMatch;//!
  TH2F  *fEMCTrkMatch;//!
  TH2F  *fEMCClsEtaPhiTrkMatch;//!
  TH2F  *fEMCTrkMatch_Phi;//!      
  TH2F  *fEMCTrkMatch_Eta;//! 

  TF1   *fFuncPtDepEta;//! pT dependent Eta cut 
  TF1   *fFuncPtDepPhi;//! pT dependent phi cut

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
  TH2F  *fHadEovpNL_AftEID;//!

  TH2F  *fHadDCA;//!
  TH1F  *fInclsElecPt;//!  
  TH2F  *fInclElecDCA;//!

  TH2F  *fEopNL_AftEID;//!

  TH1F  *fNElecInEvt;//!

  TH2F  *fInvmassULSPt;//!
  TH2F  *fInvmassLSPt;//!

  TH1F  *fULSElecPt;//!
  TH1F  *fLSElecPt;//!

  TH2F  *fULSElecDCA;//! 
  TH2F  *fLSElecDCA;//!  
   
 //Used in the function FindMother
  Bool_t        fIsHFE1;//
  Bool_t        fIsHFE2;//
  Bool_t        fIsNonHFE;//
  Bool_t        fIsFromD;//
  Bool_t        fIsFromBarionB;//
  Bool_t        fIsFromMesonB;//
  Bool_t        fIsFromBarionBD;//
  Bool_t        fIsFromMesonBD;//
  Bool_t        fIsFromPi0;//
  Bool_t        fIsFromEta;//
  Bool_t        fIsFromGamma;//

  AliAODMCHeader    *fMCHeader;//!
  AliAODMCParticle  *fMCparticle;//!
  TClonesArray      *fMCArray;//!

  AliAODMCParticle  *fMCparticleMother;//!
  AliAODMCParticle  *fMCparticleGMother;//!
  AliAODMCParticle  *fMCparticleGGMother;//!
  AliAODMCParticle  *fMCparticleGGGMother;//!

  Int_t       fNTotMCpart;//!
  Int_t       fNpureMC;//! N of particles from main generator (Hijing/Pythia)
  Int_t       fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
  Int_t       fNembMCeta; //! N > fNembMCeta = particles from eta generator

  //non hfe
  Bool_t     fIsFrmEmbPi0;//
  Bool_t     fIsFrmEmbEta;//
  Int_t      ftype;//!
  Double_t   fWeight;//!

  Double_t   fWeightPi0;//!
  Double_t   fWeightEta;//!
  
  TF1       *fPi0Weight;//!
  TF1       *fEtaWeight;//!

  THnSparse *fSprsPi0EtaWeightCal;//!

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

 Int_t     fnBinsDCAHisto;//!

 TH1F       *fBHadpT;//!
 TH1F       *fBMesonpT;//!
 TH1F       *fBDHadpT;//!
 TH1F       *fDHadpT;//!
 TH1F       *fDMesonpT;//!
 TH1F       *fD0pT;//!
 TH1F       *fDPluspT;//!
 TH1F       *fDspT;//!
 TH1F       *fLambdaCpT;//!

 TH1F       *fDcent;//
 TH1F       *fDUp;//
 TH1F       *fDDown;//
 TH1F       *fBcent;//
 TH1F       *fBMin;//
 TH1F       *fBMax;//
 TH1F       *fD0;//
 TH1F       *fDPlus;//
 TH1F       *fDs;//
 TH1F       *fLc;//

 Double_t  fWeightB;//!
 Double_t  fWeightBMin;//!
 Double_t  fWeightBMax;//!
 Double_t  fWeightD;//!
 Double_t  fWeightDUp;//!
 Double_t  fWeightDDown;//!

 TH2F       *fDElecDCA;//!
 TH2F       *fBElecDCA;//!
 TH2F       *fBHadElecDCA;//!
 TH2F       *fBMesonElecDCA;//!
 TH2F       *fBBaryonElecDCA;//!
 TH2F       *fDHadElecDCA;//!
 TH2F       *fDMesonElecDCA;//!
 TH2F       *fDBaryonElecDCA;//!
 TH2F       *fLambdaCElecDCA;//!
 TH2F       *fD0ElecDCA;//!

 THnSparse  *fSprsTemplatesNoWeight;//!
 THnSparse  *fSprsTemplatesWeight;//!
 THnSparse  *fSprsTemplatesWeightVar1;//!
 THnSparse  *fSprsTemplatesWeightVar2;//!

 TH1F     *fInclElePhysPriAll;//!
 TH1F     *fHFEPhysPriAll;//!
 TH1F     *fBEPhysPriAll;//!
 TH1F     *fDEPhysPriAll;//!

 TH1F     *fInclElePhysPriTrkCuts;//!
 TH1F     *fHFEPhysPriTrkCuts;//!
 TH1F     *fBEPhysPriTrkCuts;//!
 TH1F     *fDEPhysPriTrkCuts;//!

 TH1F     *fInclElePhysPriEMCMatch;//!
 TH1F     *fHFEPhysPriEMCMatch;//!
 TH1F     *fBEPhysPriEMCMatch;//!
 TH1F     *fDEPhysPriEMCMatch;//!

 TH1F     *fInclElePhysPriEovP;//!
 TH1F     *fHFEPhysPriEovP;//!
 TH1F     *fBEPhysPriEovP;//!
 TH1F     *fDEPhysPriEovP;//!
    
 TH1F    *fInclElePhysPriTPCnsig;//!
 TH1F    *fHFEPhysPriTPCnsig;//!
 TH1F    *fBEPhysPriTPCnsig;//!
 TH1F    *fDEPhysPriTPCnsig;//!

 TH1F   *fInclElePhysPriSS;//!
 TH1F   *fHFEPhysPriSS;//!
 TH1F   *fBEPhysPriSS;//!
 TH1F   *fDEPhysPriSS;//!

  AliAnalysisHFEppEMCalBeauty(const AliAnalysisHFEppEMCalBeauty&); // not implemented
  AliAnalysisHFEppEMCalBeauty& operator=(const AliAnalysisHFEppEMCalBeauty&); // not implemented
  
  ClassDef(AliAnalysisHFEppEMCalBeauty, 1); // example of analysis
};

#endif
