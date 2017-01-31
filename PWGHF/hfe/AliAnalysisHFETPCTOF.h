#ifndef AliAnalysisHFETPCTOF_cxx
#define AliAnalysisHFETPCTOF_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Heavy-flavour electron analysis in pp collisions     //
//      (and also Pb-Pb)             								  //
//																	  //
//		v1.0														  //
//                                                                    //
//	    Authors 							                          //
//		Camila de Conti (camila.de.conti@cern.ch)				      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliPIDResponse;
class AliCentrality;
class AliAODEvent;
class AliVEvent;
class AliAODMCHeader;
class AliSelectNonHFE;
class AliEventPoolManager;
class AliEventPool;
class TObjArray;

//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
//#include "AliMultSelection.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisHFETPCTOF : public AliAnalysisTaskSE 
{
//______________________________________________________________________
	public:
	AliAnalysisHFETPCTOF();
	AliAnalysisHFETPCTOF(const char *name);
	virtual ~AliAnalysisHFETPCTOF();
  
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

	//Setters
	void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
	void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
	
	void SetMCanalysis() {fIsMC = kTRUE;};
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetPPanalysis(Bool_t IsPP) {fIsPP = IsPP;};


	//Getters
	AliHFEpid *GetPID() const {return fPID;};
//______________________________________________________________________
  
//______________________________________________________________________
	private:
	
//Function to process track cuts
	Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

//Find Mothers (Finde HFE and NonHFE from MC information)
	Bool_t FindMother(Int_t mcIndex);

//Flags for specifcs analysis
    Bool_t				fIsMC;  
    Bool_t				fIsPP;  
    
//Weight to normalize the amount of pi and eta    
    Double_t CalculateWeight(Int_t pdg_particle, Double_t x);


//Used in the function FindMother
	Bool_t				fIsHFE1;
	Bool_t				fIsHFE2;
	Bool_t				fIsNonHFE;
	Bool_t				fIsFromD;
	Bool_t				fIsFromB;
	Bool_t				fIsFromPi0;
	Bool_t				fIsFromEta;
	Bool_t				fIsFromGamma;
	
//General variables
	AliESDEvent 			*fESD;
	AliAODEvent 		   	*fAOD;				/// new
	AliVEvent 		      	*fVevent;			/// new
	TList       			*fOutputList;
	AliPIDResponse 			*fPidResponse;
	AliSelectNonHFE 		*fNonHFE;
	
	Bool_t				fMassCutFlag;
	Bool_t				fAngleCutFlag;
	Bool_t				fChi2CutFlag;
	Bool_t				fDCAcutFlag;
	
    Double_t		    	fMassCut;
    Double_t			    fAngleCut;
    Double_t			    fChi2Cut;
    Double_t			    fDCAcut;
    
    Double_t			    fPtMinAsso;
    Int_t			        fTpcNclsAsso;
    
    AliESDtrackCuts         *fPartnerCuts;
    
	
//For the case of AOD analysis
	Bool_t					fIsAOD;					//flag for AOD analysis

//Vertex selection
	Float_t					fZvtx;
	
//EMCal
	//AliESDCaloCluster 		*fClus;
	AliVCluster				*fClus;
	
//Histograms for the analysis
	TH1F				*fNevent;
	TH1F				*fNevent2;	
	TH1F				*fCent;	
	TH2F				*fTPC_p1;
	TH2F				*fTPC_p2;
	TH2F				*fTPC_p3;
	TH2F				*fTPCnsigma_p1;
	TH2F				*fTPCnsigma_p2;
	TH2F				*fTPCnsigma_p3;
	TH2F				*fTPCnsigma_pt1;
	TH2F				*fTPCnsigma_pt2;
	TH2F				*fTPCnsigma_pt3;
	TH2F                *fTOFnsigma_p1;
	TH2F                *fTOFnsigma_p2;
	TH2F                *fTOFnsigma_p3;
	TH2F                *fTOFnsigma_pt1;
	TH2F                *fTOFnsigma_pt2;
	TH2F                *fTOFnsigma_pt3;
	TH2F                *fITSnsigma_p1;
	TH2F                *fITSnsigma_p2;
	TH2F                *fITSnsigma_p3;
	TH2F                *fITSnsigma_pt1;
	TH2F                *fITSnsigma_pt2;
	TH2F                *fITSnsigma_pt3;
	//TH2F				**fTPCnsigma_EoverP;
	TH2F                *fTPCnsigma_TOFnsigma1;
	TH2F                *fTPCnsigma_TOFnsigma2;
	TH2F                *fTPCnsigma_TOFnsigma3;
	//TH1F				**fVtxZ;
    //TH1F				**fNClusters;
    //TH2F                **fEtaPhiTracks;
    //TH2F                **fEtaPhiClusters;
	TH1F				*fInvMass;
	TH1F				*fInvMassBack;
	TH2F                *fTPCnsigma_p_after_tof;
	TH2F                *fTPCnsigma_pt_after_tof;
	TH2F                *fTPCnsigma_p_after_tof_its;
	TH2F                *fTPCnsigma_pt_after_tof_its;
	TH2F                *fTPCnsigma_p_after_tof_inverse;
	TH2F                *fTPCnsigma_pt_after_tof_inverse;
	
	TH1F				*fDCA;
	TH1F				*fDCABack;
	TH1F				*fAngle;
	TH1F				*fAngleBack;

    TH1F                *fPtElec_ULS;
    TH1F                *fPtElec_LS;
    TH1F                *fPtElec_ULS_MC;
    TH1F                *fPtElec_LS_MC;
    TH1F                *fPtElec_ULS_MC_from_pi0;
    TH1F                *fPtElec_LS_MC_from_pi0;
    TH1F                *fPtElec_ULS_MC_from_eta;
    TH1F                *fPtElec_LS_MC_from_eta;
    TH1F                *fPtElec_ULS_MC_from_gamma;
    TH1F                *fPtElec_LS_MC_from_gamma;
    TH1F                *fPtElec;
    TH2F                *fDCAz_pt;
    TH2F                *fDCAxy_pt;
    TH2F                *fDCAz_pt_gamma;
    TH2F                *fDCAxy_pt_gamma;
    TH2F                *fDCAz_pt_light;
    TH2F                *fDCAxy_pt_light;
    TH2F                *fDCAz_pt_B;
    TH2F                *fDCAxy_pt_B;
    TH2F                *fDCAz_pt_C;
    TH2F                *fDCAxy_pt_C;
    TH1F                *fPt_elec_phot;
    TH1F                *fPt_elec_from_pi0;
    TH1F                *fPt_elec_from_eta;
    TH1F                *fPt_elec_from_gamma;
    TH1F                *fPt_HFelec_passTOFandTPC_2;	
    TH1F                *fPt_HFelec_passTPC_2;	
    TH1F				*fPt_HFelec_pass_track;
    TH1F				*fPt_HFelec_pass_trackandTOF;
    TH1F                *fPtMCpi0;
    TH1F				*fPtMCpi0_feeddown;
    TH1F				*fPtMCpi0_nomother;
    TH1F                *fPtMCeta;
    TH1F				*fPtMCeta_feeddown;
    TH1F				*fPtMCeta_nomother;
    TH1F                *fPtMCeta_afterNorm;
    TH1F                *fPtMCpi0_afterNorm;
    TH1F                *fPtMCpi02;
    TH1F                *fPtMCeta2;
    TH1F                *fPtMCpi03;
    TH1F                *fPtMCeta3;
    TH1F                *fPtHFEMC_2;

		
//For the HFE package
	AliHFEcuts 			*fCuts;                 		// Cut Collection for HFE
	AliCFManager 		*fCFM;                  		// Correction Framework Manager
	AliHFEpid 			*fPID;                  		// PID
	AliHFEpidQAmanager 	*fPIDqa;						// PID QA manager
	
//Others
	AliStack 			*fMCstack;						//
	Bool_t              fRejectKinkMother;				//
	TParticle 			*fMCtrack;
	TParticle 			*fMCtrackMother;
	TParticle 			*fMCtrackGMother;
	TParticle 			*fMCtrackGGMother;
	TParticle 			*fMCtrackGGGMother;
	TClonesArray 		*fMCarray;
	AliAODMCHeader 		*fMCheader;
	AliAODMCParticle 	*fMCparticle;
	AliAODMCParticle 	*fMCparticleMother;
	AliAODMCParticle 	*fMCparticleGMother;
	AliAODMCParticle 	*fMCparticleGGMother;
	AliAODMCParticle 	*fMCparticleGGGMother;
	AliMCEventHandler	*fEventHandler;
	AliMCEvent			*fMCevent;
    
//______________________________________________________________________

	AliAnalysisHFETPCTOF(const AliAnalysisHFETPCTOF&); 			// not implemented
	AliAnalysisHFETPCTOF& operator=(const AliAnalysisHFETPCTOF&); 		// not implemented
  
	ClassDef(AliAnalysisHFETPCTOF, 1); 								// example of analysis
//______________________________________________________________________
};

///_________________________________________________________________________________________________
///Class copied from : $ALICE_ROOT/PWGCF/Correlations/DPhi/AliAnalysisTaskLongRangeCorrelations.h
///Author: Christoph Mayer
class EHCParticle : public TObject {
public:
  EHCParticle(Double_t eta=0, Double_t phi=0, Double_t pt=0)
    : fEta(eta), fPhi(phi), fPt(pt) {}
  virtual ~EHCParticle() {}

  Double_t Eta() const { return fEta; }
  Double_t Phi() const { return fPhi; }
  Double_t Pt() const { return fPt; }

protected:
private:
  EHCParticle(const EHCParticle&);
  EHCParticle& operator=(const EHCParticle&);

  Double_t fEta;
  Double_t fPhi;
  Double_t fPt;
  
  ClassDef(EHCParticle, 1);
} ;
///_________________________________________________________________________________________________

#endif
