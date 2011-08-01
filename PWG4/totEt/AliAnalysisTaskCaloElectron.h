#ifndef ALIANALYSISTASKCALOELECTRON_H
#define ALIANALYSISTASKCALOELECTRON_H

//_________________________________________________________________________
//  Task for transverse energy studies
//
//*-- Author: Marcelo G. Munhoz (USP)
//_________________________________________________________________________

class TH1F;
class TH2F;
class TTree;
class TVector3;
class TRefArray;
class TString;
class AliESDEvent;
class AliESDVertex;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpid;
class AliESDCaloCluster;
class AliESDCaloCells;
class AliEMCALTrack;
class AliEMCALGeometry;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliStack;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCaloElectron : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskCaloElectron(const char *name="AliAnalysisTaskCaloElectron", Bool_t isMc=kFALSE);
	virtual ~AliAnalysisTaskCaloElectron();
	
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
	virtual void Run(); 
	Double_t RunPhotonic(AliESDtrack *track=0, TString type="", AliESDCaloCluster *caloCluster=0);
	
	virtual Bool_t IsGoodEvent();
	virtual Bool_t IsGoodTrack(AliESDtrack *track);
	virtual Bool_t IsElectronTPC(AliESDtrack *track);
	virtual Bool_t IsElectronEMCAL(AliESDCaloCluster *caloCluster, AliESDtrack *track);
	virtual Bool_t IsElectronMC(Int_t iPart=0, TString type = "Any");
	Int_t GetMCMother(TParticle *part);
	Int_t GetMCMother(Int_t iPart);

	
	Double_t CalcMass(AliESDtrack *trackIN, AliESDtrack *track);
	AliESDtrackCuts* GetTrackCuts() const;
	AliESDtrackCuts* GetTrackPhCuts() const;
	virtual Bool_t GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos);
	void InitCaloUtil();
	AliESDCaloCluster* FindMatch(AliEMCALTrack *track, TVector3 &trackPos);
	AliESDCaloCluster* FindMCMatch(Int_t iPart);

	virtual void SetTriggerSelection(Bool_t triggerFlag = kFALSE, const char *trigger = "");
	virtual Int_t CheckPhysicsSelection(Int_t runNumber);
	virtual Bool_t IsPhysicsSelected();
	
private:
    Bool_t fIsMc; // Are we analysing MC data
	Bool_t fSelectTrigger;        //Select trigger

	TString	fTrigger;		//Trigger

    Int_t fCurrentRunNum; // The current run number
	
	Float_t fgMassElectron;//Marcelo add comment
	Float_t fCutRMatch;//Marcelo add comment
	Float_t fCutEoverPMin;//Marcelo add comment
	Float_t fCutEoverPMax;//Marcelo add comment
	Float_t fCutdEdxMin;//Marcelo add comment
	Float_t fCutdEdxMax;//Marcelo add comment
	
	AliESDEvent *fESD;//Marcelo add comment    //! ESD object
	AliESDVertex* fVertex;//Marcelo add comment
	AliStack *fStack;//Marcelo add comment
	AliESDtrackCuts *fEsdTrackCuts;//Marcelo add comment
	AliESDtrackCuts *fEsdTrackPhCuts;//Marcelo add comment
	AliESDpid *fPID;//Marcelo add comment
	AliEMCALGeometry *fGeoUt;//Marcelo add comment
	AliEMCALGeometry *fGeom;//Marcelo add comment
	AliEMCALRecoUtils *fRecoUtil;//Marcelo add comment
	TRefArray *fCaloClusters;//Marcelo add comment
	AliESDCaloCells *fCells;//Marcelo add comment	


	Float_t fVtxStat, fEta, fPt, fP, fDedx, fNSigma, fDcaXY, fDcaZ, fE, fEtaClu, fPhiClu, fRes, fInvMass;//Marcelo add comment
	Float_t fEtaMCIn, fPtMCIn, fPMCIn, fEMCIn, fvXMCIn, fvYMCIn, fvZMCIn, fIsPhyPrimIn, fMomIn;//Marcelo add comment
	Float_t fEtaMC, fPtMC, fPMC, fDedxMC, fNSigmaMC, fDcaXYMC, fDcaZMC, fEMC, fEtaCluMC, fPhiCluMC, fResMC, fEMCMatch, fPhiCluMCMatch, fEtaCluMCMatch, fResMCMatch, fInvMassMC, fIsPhyPrim, fMom;//Marcelo add comment
	
	TList       *fOutputList; //! Output list
	
	TTree		*fTree;//Marcelo add comment
	TTree		*fTreeMCIn;//Marcelo add comment
	TTree		*fTreeMC;//Marcelo add comment
	
	TH2F		*fHistVertexZ;//Marcelo add comment
	TH2F		*fHistVertexZAll;//Marcelo add comment
	
	
	TH2F        *fHistEtaPt;//Marcelo add comment 
	TH2F        *fHistEtaPteTPC;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCMatch;//Marcelo add comment 
	TH2F        *fHistEtaPteEMCAL;//Marcelo add comment 
	
	TH2F        *fHistEtaPtStackeMC;//Marcelo add comment 
	TH2F        *fHistEtaPteMC;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCeMC;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCMatcheMC;//Marcelo add comment 
	TH2F        *fHistEtaPteEMCALeMC;//Marcelo add comment 
	TH2F        *fHistEtaPteMCMatch;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCeMCMatch;//Marcelo add comment 

	TH2F        *fHistEtaPtStackePhMC;//Marcelo add comment 
	TH2F        *fHistEtaPtePhMC;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCePhMC;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCMatchePhMC;//Marcelo add comment 
	TH2F        *fHistEtaPteEMCALePhMC;//Marcelo add comment 
	TH2F        *fHistEtaPtePhMCMatch;//Marcelo add comment 
	TH2F        *fHistEtaPteTPCePhMCMatch;//Marcelo add comment 
	
	
	TH1F		*fHistCluE;//Marcelo add comment
	TH1F		*fHistCluEeTPCMatch;//Marcelo add comment
	TH1F		*fHistCluEeEMCAL;//Marcelo add comment
	
	TH1F		*fHistCluEeMC;//Marcelo add comment
	TH1F		*fHistCluEeTPCMatcheMC;//Marcelo add comment
	TH1F		*fHistCluEeEMCALeMC;//Marcelo add comment
	
	TH1F		*fHistCluEePhMC;//Marcelo add comment
	TH1F		*fHistCluEeTPCMatchePhMC;//Marcelo add comment
	TH1F		*fHistCluEeEMCALePhMC;//Marcelo add comment
	
	TH2F		*fHistCluEtaPhi;//Marcelo add comment
	TH2F		*fHistCluEtaPhieTPCMatch;//Marcelo add comment
	TH2F		*fHistCluEtaPhieEMCAL;//Marcelo add comment
	
	TH2F		*fHistCluEtaPhieMC;//Marcelo add comment
	TH2F		*fHistCluEtaPhieTPCMatcheMC;//Marcelo add comment
	TH2F		*fHistCluEtaPhieEMCALeMC;//Marcelo add comment
	
	TH2F		*fHistCluEtaPhiePhMC;//Marcelo add comment
	TH2F		*fHistCluEtaPhieTPCMatchePhMC;//Marcelo add comment
	TH2F		*fHistCluEtaPhieEMCALePhMC;//Marcelo add comment
	
	
	TH2F		*fHistVtxXYeMC;//Marcelo add comment
	
	
	TH2F		*fHistdEdx;//Marcelo add comment
	TH2F		*fHistdEdxeTPC;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCMatch;//Marcelo add comment	
	TH2F		*fHistdEdxeEMCAL;//Marcelo add comment
	
	TH2F		*fHistdEdxeMC;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCeMC;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCMatcheMC;//Marcelo add comment	
	TH2F		*fHistdEdxeEMCALeMC;//Marcelo add comment
	TH2F		*fHistdEdxeMCMatch;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCeMCMatch;//Marcelo add comment	

	TH2F		*fHistdEdxePhMC;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCePhMC;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCMatchePhMC;//Marcelo add comment	
	TH2F		*fHistdEdxeEMCALePhMC;//Marcelo add comment
	TH2F		*fHistdEdxePhMCMatch;//Marcelo add comment	
	TH2F		*fHistdEdxeTPCePhMCMatch;//Marcelo add comment	
	
	
	TH2F        *fHistDeltaRZeTPCMatch;//Marcelo add comment 
	
	TH2F        *fHistDeltaRZeTPCMatcheMC;//Marcelo add comment 
	TH2F        *fHistDeltaRZeMCMatch;//Marcelo add comment 
	TH2F        *fHistDeltaRZeTPCeMCMatch;//Marcelo add comment 

	TH2F        *fHistDeltaRZeTPCMatchePhMC;//Marcelo add comment 
	TH2F        *fHistDeltaRZePhMCMatch;//Marcelo add comment 
	TH2F        *fHistDeltaRZeTPCePhMCMatch;//Marcelo add comment 
	
	
	TH2F        *fHistResPteTPCMatch;//Marcelo add comment 
	
	TH2F        *fHistResPteTPCMatcheMC;//Marcelo add comment 
	TH2F        *fHistResPteMCMatch;//Marcelo add comment 
	TH2F        *fHistResPteTPCeMCMatch;//Marcelo add comment 
	
	TH2F        *fHistResPteTPCMatchePhMC;//Marcelo add comment 
	TH2F        *fHistResPtePhMCMatch;//Marcelo add comment 
	TH2F        *fHistResPteTPCePhMCMatch;//Marcelo add comment 
	
	
	TH2F		*fHistEoverPPteTPCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteEMCAL;//Marcelo add comment
	TH2F		*fHistEoverPPeTPCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPeEMCAL;//Marcelo add comment
	TH2F		*fHistEoverPPeEMCALtrigger;//Marcelo add comment
	TH2F		*fHistEoverPEeTPCMatch;//Marcelo add comment
	TH2F		*fHistEoverPEeEMCAL;//Marcelo add comment
	
	
	TH2F		*fHistEoverPPteTPCMatcheMC;//Marcelo add comment
	TH2F		*fHistEoverPPteEMCALeMC;//Marcelo add comment
	TH2F		*fHistEoverPPteMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteEMCALeMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteTPCeMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteTPCeEMCALeMCMatch;//Marcelo add comment
	
	TH2F		*fHistEoverPPteTPCMatchePhMC;//Marcelo add comment
	TH2F		*fHistEoverPPteEMCALePhMC;//Marcelo add comment
	TH2F		*fHistEoverPPtePhMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteEMCALePhMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteTPCePhMCMatch;//Marcelo add comment
	TH2F		*fHistEoverPPteTPCeEMCALePhMCMatch;//Marcelo add comment
	
	TH2F		*fHistMassLikePteTPC;//Marcelo add comment
	TH2F		*fHistMassUnLikePteTPC;//Marcelo add comment
	TH2F		*fHistMassLikePteEMCAL;//Marcelo add comment
	TH2F		*fHistMassUnLikePteEMCAL;//Marcelo add comment
	TH2F		*fHistMassLikePeEMCAL;//Marcelo add comment
	TH2F		*fHistMassUnLikePeEMCAL;//Marcelo add comment
	TH2F		*fHistMassLikeEeEMCAL;//Marcelo add comment
	TH2F		*fHistMassUnLikeEeEMCAL;//Marcelo add comment

	TH2F		*fHistMassLikePteTPCePhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePteTPCePhMC;//Marcelo add comment
	TH2F		*fHistMassLikePteEMCALePhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePteEMCALePhMC;//Marcelo add comment
	TH2F		*fHistMassLikePeEMCALePhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePeEMCALePhMC;//Marcelo add comment
	TH2F		*fHistMassLikeEeEMCALePhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikeEeEMCALePhMC;//Marcelo add comment
	
	TH2F		*fHistMassLikePteTPCeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePteTPCeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassLikePteEMCALeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePteEMCALeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassLikePeEMCALeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikePeEMCALeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassLikeEeEMCALeNonPhMC;//Marcelo add comment
	TH2F		*fHistMassUnLikeEeEMCALeNonPhMC;//Marcelo add comment
	
	
	AliAnalysisTaskCaloElectron(const AliAnalysisTaskCaloElectron&); // not implemented
	AliAnalysisTaskCaloElectron& operator=(const AliAnalysisTaskCaloElectron&); // not implemented
	
	ClassDef(AliAnalysisTaskCaloElectron, 1); // example of analysis
};

#endif
