
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
class AliAODMCHeader;
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





class AliAnalysisTaskHFEMultiplicity : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskHFEMultiplicity();
  AliAnalysisTaskHFEMultiplicity(const char *name);
  virtual                 ~AliAnalysisTaskHFEMultiplicity();
  virtual void   	  Init();
  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  
	
  
  Bool_t  		PassEventSelect(AliAODEvent *fAOD);
  Bool_t		Passtrackcuts(AliAODTrack *atrack);
  void    		GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff );
  void    		ClusterInfo(); 
  
  Bool_t  		GetTenderSwitch() {return fUseTender;};     
  void    		SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
  void    		SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
  void    		SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
  
  Bool_t                GetEMCalTriggerEG1() { return fEMCEG1; };
  Bool_t                GetEMCalTriggerEG2() { return fEMCEG2; };
  void                  SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
  void                  SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
  Bool_t                GetEMCalTriggerDG1() { return fDCalDG1; };
  Bool_t                GetEMCalTriggerDG2() { return fDCalDG2; };
  void                  SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;};
  void                  SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;};  
  void    		SelectNonHFElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS);

  void 			SetReferenceMultiplicity(Double_t multi){fRefMult=multi;}

  void 			SetMultiProfileLHC16s(TProfile * hprof){
        		fMultEstimatorAvg[0]=hprof;
    			}
  void 			SetMultiProfileLHC16r(TProfile * hprof){
        		fMultEstimatorAvg[1]=hprof;
    			}
  
private:
  
  AliAODEvent*          fAOD;//! input event
  AliPIDResponse*       fpidResponse;//!pid response
  
  TString		fTenderClusterName;//
  TString	        fTenderTrackName;//
  
  TClonesArray*         fTracks_tender;//Tender tracks     
  TClonesArray*         fCaloClusters_tender;//Tender cluster      
  AliAODMCParticle*	fMCparticle;
  TClonesArray*	fMCarray;//! MC array
  TProfile* 		GetEstimatorHistogram(const AliAODEvent *fAOD);

  Bool_t                fUseTender;// switch to add tender
  Bool_t                fFlagClsTypeEMC;//switch to select EMC clusters
  Bool_t                fFlagClsTypeDCAL;//switch to select DCAL c
  Bool_t                fEMCEG1;//EMcal Threshold EG1
  Bool_t                fEMCEG2;//EMcal Threshold EG2
  Bool_t                fDCalDG1;//DCal Threshold DG1
  Bool_t                fDCalDG2;//DCal Threshold DG2
  Bool_t 		fRejectPUFromSPD;
  
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

  TH1F*                fInvmassLS;//!
  TH1F*                fInvmassULS;//!
  TH2F*                fInvmassLSPt;//!
  TH2F*                fInvmassULSPt;//!
  TH1F*                fULSElecPt;//!
  TH1F*                fLSElecPt;//!

  Bool_t 		fReadMC;    //flag for access to MC
  TProfile*		fMultEstimatorAvg[2];
  Double_t 		fRefMult;
  TRandom3*		gRandom;//!< random number generator
 
  THnSparse* 	        fSparseElectron; //! Electron information
 
  THnSparse* 	        fSparseLSElectron; //! Electron information
  THnSparse* 	        fSparseULSElectron; //! Electron information
  Double_t*             fvaluePHElectron;//!
  Double_t*             fvalueElectron; //!
 
  THnSparse* 	        fSparseMulti; //! Multiplicity information
  Double_t*             fvalueMulti; //!

 

  AliAnalysisTaskHFEMultiplicity(const AliAnalysisTaskHFEMultiplicity& ); // not implemented
  AliAnalysisTaskHFEMultiplicity& operator=(const AliAnalysisTaskHFEMultiplicity& ); // not implemented

  ClassDef(AliAnalysisTaskHFEMultiplicity, 1);
};



#endif
