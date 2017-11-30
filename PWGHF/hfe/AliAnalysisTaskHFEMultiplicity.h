
#ifndef AliAnalysisTaskHFEMultiplicity_H
#define AliAnalysisTaskHFEMultiplicity_H
class TH1F;
class TH2F;
class TLorentzVector;
class THnSparse;
class TRandom3;

class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliAODMCParticle; // sample
class AliEMCALTriggerPatchInfo;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliMultSelection;


#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "TProfile.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliVEvent.h"
#include "AliSelectNonHFE.h"




class AliAnalysisTaskHFEMultiplicity : public AliAnalysisTaskSE  
{
 public:

  enum HijingOrNot {kHijing,kElse};
  enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
  
  AliAnalysisTaskHFEMultiplicity();
  AliAnalysisTaskHFEMultiplicity(const char *name);
  virtual                 ~AliAnalysisTaskHFEMultiplicity();
  virtual void   	  Init();
  virtual void 		  LocalInit() {Init();}
  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  
	
  
  Bool_t  		PassEventSelect(AliAODEvent *fAOD);
  Bool_t		Passtrackcuts(AliAODTrack *atrack);
  Bool_t 		PassEIDCuts(AliAODTrack *track, AliAODCaloCluster *clust);
  void    		GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff );
  Int_t 		GetNcharged();
  void 			GetElectronFromStack();
  void   	        GetPi0EtaWeight(THnSparse *SparseWeight);
  Bool_t  		GetNMCPartProduced();
  Int_t                 GetPi0EtaType(AliAODMCParticle *part);
  void    		SwitchPi0EtaWeightCalc(Bool_t fSwitch) {fCalculateWeight = fSwitch;};
  Bool_t  		GetTenderSwitch() {return fUseTender;};     
  void    		SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
  void    		SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
  void    		SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
  void 			SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  
  Bool_t                GetEMCalTriggerEG1() { return fEMCEG1; };
  Bool_t                GetEMCalTriggerEG2() { return fEMCEG2; };
  void                  SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
  void                  SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
  Bool_t                GetEMCalTriggerDG1() { return fDCalDG1; };
  Bool_t                GetEMCalTriggerDG2() { return fDCalDG2; };
  void                  SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;};
  void                  SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;};
  Bool_t  		IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);  
  void    		SelectNonHFElectron(Int_t itrack, AliAODTrack *track, Int_t iMCmom, Int_t MomPDG, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS, Int_t vzeroMultCorr,Int_t correctednAcc);
  void 			SetReferenceMultiplicity(Double_t multi){fRefMult=multi;}

  void 			SetMultiProfileLHC16s(TProfile * hprof){
			 if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
        		fMultEstimatorAvg[0]=new TProfile(*hprof);
        		
    			}
  void 			SetMultiProfileLHC16r(TProfile * hprof){
        		 if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
        		 fMultEstimatorAvg[1]=new TProfile(*hprof);
    			}
  
private:
  
  AliAODEvent*          fAOD;//! input event
  AliPIDResponse*       fpidResponse;//!pid response
  
  TString		fTenderClusterName;//
  TString	        fTenderTrackName;//
  
  TClonesArray*         fTracks_tender;//Tender tracks     
  TClonesArray*         fCaloClusters_tender;//Tender cluster      
  AliAODMCParticle*	fMCparticle;
  AliAODMCHeader*	fMCHeader;//!
  TClonesArray*	        fMCArray;//! MC array
  TProfile* 		GetEstimatorHistogram(const AliAODEvent *fAOD);

  Bool_t                fUseTender;// switch to add tender
  Bool_t                fFlagClsTypeEMC;//switch to select EMC clusters
  Bool_t                fFlagClsTypeDCAL;//switch to select DCAL c
  Double_t              fTPCnSigma;//!
  Bool_t                fEMCEG1;//EMcal Threshold EG1
  Bool_t                fEMCEG2;//EMcal Threshold SetReferenceMultiplicityEG2
  Bool_t                fDCalDG1;//DCal Threshold DG1
  Bool_t                fDCalDG2;//DCal Threshold DG2
  Bool_t 		fRejectPUFromSPD;
  Bool_t                fCalculateElectronEffi;//

  Bool_t              	fCalculateWeight;//
  Int_t               	fNTotMCpart; //! N of total MC particles produced by generator
  Int_t              	fNpureMC;//! N of particles from main generator (Hijing/Pythia)
  Int_t               	fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
  Int_t               	fNembMCeta; //! N > fNembMCeta = particles from eta generator
  
  TList*                fOutputList;//! output list
  TList*		fListProfiles; // list of profile histos for z-vtx correction
 
  TH1F*			fNevents;//! no of events
  TH1F*			fClusPhi;//! Cluster Phi
  TH1F*			fClusEta;//! Cluster Eta
  TH2F*			fClusEtaPhi;//! Cluster Eta Phi
  TH1F*			fNCells;//! Number of cells in a cluster
  TH1F*			fClusE;//! Cluster Energy
  TH1F*			fClusT;//! Cluster Time
  TH1F*			fCellE;//! Cell Energy
  TH1F*			fCellT;//! Cell Time
  TH1F*			fVtxZ;//! Vertex z
  TH1F*			fVtxX;//! Vertex x
  TH1F*			fVtxY;//! Vertex y
  TH2F*			fTPCdEdx;//! dedx vs pt
  TH2F*			fTPCnsigma;//! TPC Nsigma
  TH1F*			fTrkPt;//! track pt
  TH1F*			fTrketa;//! track eta
  TH1F*			fTrkphi;//! track phi
  TH1F*			fTrkMatchTrkPt;//!
  TH1F*			fTrkMatchTrketa;//! cluster matched track eta
  TH1F*			fTrkMatchTrkphi;//! cluster matched track phi
  TH2F*			fTrkMatchClusetaphi;//! matched cluster eta phi
  TH2F*			fEMCTrkMatchcluster;//!Distance of EMC cluster to closest track in x and z


  TH1F*                 fInvmassLS;//!
  TH1F*                 fInvmassULS;//!
  TH2F*                 fInvmassLSPt;//!
  TH2F*                 fInvmassULSPt;//!
  TH1F*                 fULSElecPt;//!
  TH1F*                 fLSElecPt;//!
//MC histograms-------------------------------------------
  TH1F*                	fNonHFePairInvmassLS;//!
  TH1F*               	fNonHFePairInvmassULS;//!
  TH1F*                	fNonHFeEmbInvmassLS;//!
  TH1F*                	fNonHFeEmbInvmassULS;//!
  TH1F*                	fNonHFeEmbWeightInvmassLS;//!
  TH1F*                	fNonHFeEmbWeightInvmassULS;//!
  TH1F*                	fPi0EmbInvmassLS;//!
  TH1F*                	fPi0EmbInvmassULS;//!
  TH1F*                	fPi0EmbWeightInvmassLS;//!
  TH1F*                	fPi0EmbWeightInvmassULS;//!
  TH1F*                	fEtaEmbInvmassLS;//!
  TH1F*                	fEtaEmbInvmassULS;//!
  TH1F*                	fEtaEmbWeightInvmassLS;//!
  TH1F*                	fEtaEmbWeightInvmassULS;//!
  TH2F*                	fRecoLSeTrkPt;//!
  TH2F*                	fRecoLSeEmbTrkPt;//!
  TH2F*                	fRecoLSeEmbWeightTrkPt;//!
  TH2F*                	fRecoPi0LSeEmbWeightTrkPt;//!
  TH2F*                	fRecoEtaLSeEmbWeightTrkPt;//!
  TH2F*                	fRecoULSeTrkPt;//!
  TH2F*                	fRecoULSeEmbTrkPt;//!
  TH2F*                	fRecoULSeEmbWeightTrkPt;//!
  TH2F*                	fRecoPi0ULSeEmbWeightTrkPt;//!
  TH2F*                	fRecoEtaULSeEmbWeightTrkPt;//!

  Bool_t 		fReadMC;    //flag for access to MC
  Bool_t                fIsFrmEmbPi0;//!
  Bool_t                fIsFrmEmbEta;//!
  Int_t                 ftype;//!
  Double_t              fWeight;//!
  TH2F*                 fInclsElecPt;//!
  TH1F*                 fInclsElecPtAll;//!
  TH2F*                 fInclsElecPtReco;//!
  TH2F*                 fInclseEMCALElecPtReco;//!
  TH2F*                 fInclseDCALElecPtReco;//!
  TH2F*			fHFElecPtAll;//!
  TH2F*			fHFElecPtReco_wTPC;//!
  TH2F*			fHFElecPtReco_wCalo;//!
  TH2F*                 fMissingEmbEtaEleTrkPt;//!
  TF1*                  fPi0Weight;//!
  TF1*                  fEtaWeight;//!
  TH2F*                	fNonHFeTrkPt;//!
  TH2F*              	fNonHFeEmbAllTypeTrkPt;//!
  TH2F*                	fNonHFeEmbTrkPt;//!
  TH2F*                	fNonHFeEmbWeightTrkPt;//!
  TH2F*                	fPi0eEmbWeightTrkPt;//!
  TH2F*                	fEtaeEmbWeightTrkPt;//!
  TH2F*                	fRecoNonHFeTrkPt;//!
  TH2F*                	fRecoNonHFeEmbTrkPt;//!
  TH2F*                	fRecoNonHFeEmbWeightTrkPt;//!
  TH2F*                	fRecoPi0eEmbWeightTrkPt;//!
  TH2F*                	fRecoEtaeEmbWeightTrkPt;//!

  TH2F*                	fNonHFeEmbMomTrkPt;//!
  TH2F*                	fNonHFeEmbMomWeightTrkPt;//!
  TH2F*                	fPi0eEmbMomWeightTrkPt;//!
  TH2F*                	fEtaeEmbMomWeightTrkPt;//!
  TH2F*                	fRecoNonHFeMomEmbTrkPt;//!
  TH2F*                	fRecoNonHFeEmbMomWeightTrkPt;//!
  TH2F*                	fRecoPi0eEmbMomWeightTrkPt;//!
  TH2F*                	fRecoEtaeEmbMomWeightTrkPt;//!
  THnSparse*		fSprsPi0EtaWeightCal;//!

//--------------------------------------------------------
  TProfile*		fMultEstimatorAvg[2];
  Double_t 		fRefMult;
  TRandom3*		gRandom;//!< random number generator
 
  THnSparse* 	        fSparseElectron; //! Electron information
  THnSparse*		fSparseClusE;//cluster energy before track matching
  THnSparse* 	        fSparseLSElectron; //! Electron information
  THnSparse* 	        fSparseULSElectron; //! Electron information
  Double_t*             fvalueCluE; //!
  Double_t*             fvaluePHElectron;//!
  Double_t*             fvalueElectron; //!
 
 
  THnSparse* 	        fSparseMulti; //! Multiplicity information
  Double_t*             fvalueMulti; //!

 

  
 

  AliAnalysisTaskHFEMultiplicity(const AliAnalysisTaskHFEMultiplicity& ); // not implemented
  AliAnalysisTaskHFEMultiplicity& operator=(const AliAnalysisTaskHFEMultiplicity& ); // not implemented

  ClassDef(AliAnalysisTaskHFEMultiplicity, 1);
};



#endif
