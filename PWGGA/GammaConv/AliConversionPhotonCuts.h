#ifndef ALICONVERSIONPHOTONCUTS_H
#define ALICONVERSIONPHOTONCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Friederike Bock

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliPIDResponse;
class AliKFVertex;
class TH1F;
class TH2F;
class TF1;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;

using namespace std;

class AliConversionPhotonCuts : public AliAnalysisCuts {
      
	public: 
		

		enum cutIds {
			kv0FinderType,                
			ketaCut,                                     
			kRCut,                     
			ksinglePtCut,                 
			kclsTPCCut,                   
			kededxSigmaCut,               
			kpidedxSigmaCut,              
			kpiMomdedxSigmaCut,        
			kpiMaxMomdedxSigmaCut,        
			kLowPRejectionSigmaCut,       
			kTOFelectronPID,              
			kQtMaxCut,                    
			kchi2GammaCut,                
			kPsiPair, 
			kdoPhotonAsymmetryCut,
			kCosPAngle,
			kElecShare,
			kToCloseV0s,
			kDcaRPrimVtx,
			kDcaZPrimVtx,
			kInPlaneOutOfPlane,
			kNCuts
		};

		enum photonCuts {
				kPhotonIn=0,
				kOnFly,
				kNoTracks,
				kTrackCuts,
				kdEdxCuts,
				kConvPointFail,
				kPhotonCuts,
				kEventPlane,
				kPhotonOut
		};


		Bool_t SetCutIds(TString cutString); 
		Int_t fCuts[kNCuts];
		Bool_t SetCut(cutIds cutID, Int_t cut);
		Bool_t UpdateCutString();

		static const char * fgkCutNames[kNCuts];

		Double_t GetCosineOfPointingAngle(const AliConversionPhotonBase * photon, AliVEvent * event) const; 
		Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
		void FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);
		void SetPreSelectionCutFlag(Bool_t preSelFlag){fPreSelCut = preSelFlag;}   

		AliConversionPhotonCuts(const char *name="V0Cuts", const char * title="V0 Cuts");
		AliConversionPhotonCuts(const AliConversionPhotonCuts&);
		AliConversionPhotonCuts& operator=(const AliConversionPhotonCuts&);

		virtual ~AliConversionPhotonCuts();                            //virtual destructor

		static AliConversionPhotonCuts * GetStandardCuts2010PbPb();
		static AliConversionPhotonCuts * GetStandardCuts2010pp();

		void InitAODpidUtil(Int_t type);
		Bool_t InitPIDResponse();
		void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
		AliPIDResponse * GetPIDResponse() { return fPIDResponse;}

		
		virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
		virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

		TString GetCutNumber();
		
		// Cut Selection
		Bool_t PhotonIsSelected(AliConversionPhotonBase * photon, AliVEvent  * event);
		Bool_t PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack,Bool_t checkForConvertedGamma=kTRUE);
		Bool_t PhotonIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray,Bool_t checkForConvertedGamma=kTRUE);
		Bool_t ElectronIsSelectedMC(TParticle *particle,AliStack *fMCStack);
		Bool_t TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack);
		Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE);
		Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack, Bool_t bMCDaughtersInAcceptance=kFALSE);
			
		void PrintCuts();
		void PrintCutsWithValues();

		void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
		void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
		TList *GetCutHistograms(){return fHistograms;}
		void FillPhotonCutIndex(Int_t photoncut){if(hCutIndex)hCutIndex->Fill(photoncut);}
		void FillV0EtaBeforedEdxCuts(Float_t v0Eta){if(hEtaDistV0s)hEtaDistV0s->Fill(v0Eta);}
		void FillV0EtaAfterdEdxCuts(Float_t v0Eta){if(hEtaDistV0sAfterdEdxCuts)hEtaDistV0sAfterdEdxCuts->Fill(v0Eta);}
		
		static AliVTrack * GetTrack(AliVEvent * event, Int_t label);
		static AliESDtrack *GetESDTrack(AliESDEvent * event, Int_t label);
		
		///Cut functions
		Bool_t SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex);
		Bool_t SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex);
		Bool_t AcceptanceCuts(AliConversionPhotonBase *photon);
		Bool_t AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg);
		Bool_t dEdxCuts(AliVTrack * track);
		Bool_t ArmenterosQtCut(AliConversionPhotonBase *photon);
		Bool_t AsymmetryCut(AliConversionPhotonBase *photon,AliVEvent *event);
		Bool_t PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event);
		Bool_t SelectV0Finder(Bool_t onfly){
			if(onfly == fUseOnFlyV0Finder) return kTRUE;
			else return kFALSE;
		}
		Bool_t PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event);
		Bool_t CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event);
		Bool_t PsiPairCut(const AliConversionPhotonBase * photon) const;
		Bool_t CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const;
		Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
		Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);

		UChar_t DeterminePhotonQualityAOD(AliAODConversionPhoton*, AliVEvent*);
		Bool_t InPlaneOutOfPlaneCut(Double_t photonPhi, Double_t eventPlaneAngle = -100, Bool_t fill = kTRUE);
		Int_t GetInPlaneOutOfPlaneCut(){return fInPlaneOutOfPlane;}

		// Set Individual Cuts
		Bool_t SetRCut(Int_t RCut);
		Bool_t SetV0Finder(Int_t v0FinderType);
		Bool_t SetChi2GammaCut(Int_t chi2GammaCut);
		Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
		Bool_t SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut);
		Bool_t SetSinglePtCut(Int_t singlePtCut);
		Bool_t SetTPCClusterCut(Int_t clsTPCCut);
		Bool_t SetEtaCut(Int_t etaCut);
		Bool_t SetMinMomPiondEdxCut(Int_t piMinMomdedxSigmaCut);
		Bool_t SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut);
		Bool_t SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut);
		Bool_t SetQtMaxCut(Int_t QtMaxCut);
		Bool_t SetTOFElectronPIDCut(Int_t TOFelectronPID);
		Bool_t SetTRDElectronCut(Int_t TRDElectronCut);
		Bool_t SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut);
		Bool_t SetCosPAngleCut(Int_t cosCut);
		Bool_t SetPsiPairCut(Int_t psiCut);
		Bool_t SetSharedElectronCut(Int_t sharedElec);
		Bool_t SetToCloseV0sCut(Int_t toClose);
		Bool_t SetDCARPhotonPrimVtxCut(Int_t DCARPhotonPrimVtx);
		Bool_t SetDCAZPhotonPrimVtxCut(Int_t DCAZPhotonPrimVtx);
		Bool_t SetInPlaneOutOfPlane(Int_t inOutPlane);
		void SetIsHeavyIon(Int_t isHeavyIon){fIsHeavyIon=isHeavyIon;}
		Int_t GetFirstTPCRow(Double_t radius);

		// Request Flags
		Bool_t UseElecSharingCut(){return fDoSharedElecCut;}
		Bool_t UseToCloseV0sCut(){return fDoToCloseV0sCut;}
		Double_t GetEtaCut(){return fEtaCut;}
			
	protected:
		TList 			*fHistograms;							//
		AliPIDResponse 	*fPIDResponse;							//

		//cuts
		Double_t 		fMaxR; 									// r cut
		Double_t 		fMinR;									// r cut
		Double_t 		fEtaCut;								// eta cut
		Double_t 		fEtaCutMin;								// eta cut
		Double_t 		fPtCut;									// pt cut
		Double_t 		fSinglePtCut;							// pt cut for electron/positron
		Double_t 		fMaxZ;									// z cut
		Double_t 		fMinClsTPC;								// minimum clusters in the TPC
		Double_t 		fMinClsTPCToF;							// minimum clusters to findable clusters
		Double_t 		fLineCutZRSlope; 						// linecut
		Double_t 		fLineCutZValue; 						// linecut
		Double_t 		fLineCutZRSlopeMin; 					// linecut
		Double_t 		fLineCutZValueMin; 						// linecut
		Double_t 		fChi2CutConversion; 					// chi2cut
		Double_t 		fPIDProbabilityCutNegativeParticle; 	//
		Double_t 		fPIDProbabilityCutPositiveParticle; 	//
		Bool_t   		fDodEdxSigmaCut;						// flag to use the dEdxCut based on sigmas
		Bool_t   		fDoTOFsigmaCut;							// flag to use TOF pid cut RRnewTOF
		Double_t 		fPIDTRDEfficiency;						// required electron efficiency for TRD PID
		Bool_t   		fDoTRDPID;								// flag to use TRD pid
		Double_t 		fPIDnSigmaAboveElectronLine;			// sigma cut
		Double_t 		fPIDnSigmaBelowElectronLine;			// sigma cut
		Double_t 		fTofPIDnSigmaAboveElectronLine;			// sigma cut RRnewTOF
		Double_t 		fTofPIDnSigmaBelowElectronLine;			// sigma cut RRnewTOF 
		Double_t 		fPIDnSigmaAbovePionLine;     			// sigma cut
		Double_t 		fPIDnSigmaAbovePionLineHighPt;			// sigma cut
		Double_t 		fPIDMinPnSigmaAbovePionLine; 			// sigma cut
		Double_t 		fPIDMaxPnSigmaAbovePionLine; 			// sigma cut
		Double_t 		fDoKaonRejectionLowP;   				// Kaon rejection at low p
		Double_t 		fDoProtonRejectionLowP; 				// Proton rejection at low p
		Double_t 		fDoPionRejectionLowP;   				// Pion rejection at low p
		Double_t 		fPIDnSigmaAtLowPAroundKaonLine; 		// sigma cut
		Double_t 		fPIDnSigmaAtLowPAroundProtonLine;		// sigma cut
		Double_t 		fPIDnSigmaAtLowPAroundPionLine;			// sigma cut
		Double_t 		fPIDMinPKaonRejectionLowP; 				// Momentum limit to apply kaon rejection
		Double_t 		fPIDMinPProtonRejectionLowP;			// Momentum limit to apply proton rejection
		Double_t 		fPIDMinPPionRejectionLowP; 				// Momentum limit to apply proton rejection
		Bool_t   		fDoQtGammaSelection; 					// Select gammas using qtMax
		Bool_t   		fDo2DQt; 								// Select gammas using ellipse cut
		Double_t 		fQtMax; 								// Maximum Qt from Armenteros to select Gammas
		Double_t 		fNSigmaMass; 							// nsigma cut
		Bool_t 			fUseEtaMinCut; 							// flag
		Bool_t 			fUseOnFlyV0Finder; 						// flag
		Bool_t   		fDoPhotonAsymmetryCut; 					// flag to use the PhotonAsymetryCut
		Double_t 		fMinPPhotonAsymmetryCut; 				// Min Momentum for Asymmetry Cut
		Double_t 		fMinPhotonAsymmetry;  					// Asymmetry Cut
		Bool_t 			fUseCorrectedTPCClsInfo; 				// flag to use corrected tpc cl info
		Bool_t 			fUseTOFpid; 							// flag to use tof pid
		Float_t 		fOpeningAngle; 							// min opening angle for meson
		Float_t 		fPsiPairCut;							//
		Bool_t 			fDo2DPsiPairChi2;						//
		Float_t 		fCosPAngleCut;							//
		Bool_t 			fDoToCloseV0sCut; 						//
		Double_t 		fminV0Dist; 							//
		Bool_t 			fDoSharedElecCut; 						//
		Bool_t 			fDoPhotonQualitySelectionCut; 			//
		Int_t 			fPhotonQualityCut; 						//
		TRandom3 		fRandom; 								//
		Int_t 			fElectronArraySize; 					// Size of electron array
		Int_t 			*fElectronLabelArray; 					//[fElectronArraySize]
		Double_t 		fDCAZPrimVtxCut; 						// cut value for the maximum distance in Z between the photon & the primary vertex [cm]
		Double_t 		fDCARPrimVtxCut; 						// cut value for the maximum distance in R between the photon & the primary vertex [cm]
		Int_t 			fInPlaneOutOfPlane; 					// In-Plane Out-Of Plane Analysis
		Float_t 		fConversionPointXArray; 				// Array with conversion Point x
		Float_t 		fConversionPointYArray; 				// Array with conversion Point y
		Float_t 		fConversionPointZArray; 				// Array with conversion Point z
		TObjString 		*fCutString; 							// cut number used for analysis
		Int_t			fIsHeavyIon;							// flag for pp (0), PbPb (1), pPb (2)
		
		// Histograms
		TH1F			*hEtaDistV0s; 							// eta-distribution of all V0s after Finder selection
		TH1F			*hEtaDistV0sAfterdEdxCuts; 				// eta-distribution of all V0s after Finder selection after dEdx cuts
		TH1F 			*hdEdxCuts;  							// bookkeeping for dEdx cuts
		TH2F 			*hTPCdEdxbefore; 						// TPC dEdx before cuts
		TH2F 			*hTPCdEdxafter; 						// TPC dEdx after cuts
		TH2F 			*hTPCdEdxSigbefore; 					// TPC Sigma dEdx before cuts
		TH2F 			*hTPCdEdxSigafter; 						// TPC Sigm dEdx after cuts
		TH2F 			*hTOFbefore; 							// TOF before cuts
		TH2F		 	*hTOFSigbefore; 						// TOF Sigma before cuts
		TH2F 			*hTOFSigafter; 							// TOF Sigma after cuts
		TH2F 			*hPsiPairDeltaPhiafter; 				// TOF Sigma after cuts
		TH1F 			*hTrackCuts; 							// bookkeeping for track cuts
		TH1F 			*hPhotonCuts; 							// bookkeeping for photon specific cuts
		TH1F 			*hInvMassbefore; 						// e+e- inv mass distribution before cuts
		TH2F 			*hArmenterosbefore;		 				// armenteros podolanski plot before cuts
		TH1F 			*hInvMassafter; 						// e+e- inv mass distribution after cuts
		TH2F 			*hArmenterosafter;  					// armenteros podolanski plot after cuts
		TH1F 			*hAcceptanceCuts; 						// bookkeeping for acceptance cuts
		TH1F 			*hCutIndex; 							// bookkeeping for cuts
		TH1F 			*hEventPlanePhi; 						// EventPlaneAngle Minus Photon Angle
		Bool_t 			fPreSelCut; 							// Flag for preselection cut used in V0Reader

	private:
	
		ClassDef(AliConversionPhotonCuts,1)
};


inline void AliConversionPhotonCuts::InitAODpidUtil(Int_t type) {
	if (!fPIDResponse) fPIDResponse = new AliAODpidUtil();
	Double_t alephParameters[5];
	// simulation
	alephParameters[0] = 2.15898e+00/50.;
	alephParameters[1] = 1.75295e+01;
	alephParameters[2] = 3.40030e-09;
	alephParameters[3] = 1.96178e+00;
	alephParameters[4] = 3.91720e+00;
	fPIDResponse->GetTOFResponse().SetTimeResolution(80.);
	
	// data
	if (type==1){
		alephParameters[0] = 0.0283086/0.97;
		alephParameters[1] = 2.63394e+01;
		alephParameters[2] = 5.04114e-11;
		alephParameters[3] = 2.12543e+00;
		alephParameters[4] = 4.88663e+00;
		fPIDResponse->GetTOFResponse().SetTimeResolution(130.);
		fPIDResponse->GetTPCResponse().SetMip(50.);
	}
	
	fPIDResponse->GetTPCResponse().SetBetheBlochParameters(
		alephParameters[0],alephParameters[1],alephParameters[2],
		alephParameters[3],alephParameters[4]);
	
	fPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}

#endif
