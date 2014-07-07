#ifndef ALIV0READERV1_H
#define ALIV0READERV1_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliConversionPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliExternalTrackParam.h"
#include "TObject.h"
#include "AliMCEvent.h"   // for CF
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include <vector>
#include "AliESDpid.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliAnalysisManager.h"

class AliConversionPhotonBase;
class TRandom3;
class AliStack;
class TList;
class AliKFConversionPhoton;
class TString;
class TClonesArray;
class TH1F;
class TH2F;
class AliAODConversionPhoton;

using namespace std;

class AliV0ReaderV1 : public AliAnalysisTaskSE {

	public:

		AliV0ReaderV1(const char *name="V0ReaderV1");
		virtual ~AliV0ReaderV1();                            //virtual destructor
		void UserCreateOutputObjects();
		virtual Bool_t Notify();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *);
		virtual void Init();

		Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
		Bool_t IsEventSelected(){return fEventIsSelected;}

		// Return Reconstructed Gammas
		TClonesArray *GetReconstructedGammas(){return fConversionGammas;}
		Int_t GetNReconstructedGammas(){if(fConversionGammas){return fConversionGammas->GetEntriesFast();}else{return 0;}}
		
		AliConversionPhotonCuts *GetConversionCuts(){return fConversionCuts;}
		AliConvEventCuts *GetEventCuts(){return fEventCuts;}
		TList *GetCutHistograms(){if(fConversionCuts){return fConversionCuts->GetCutHistograms();}return NULL;}
		TList *GetEventCutHistograms(){if(fEventCuts){return fEventCuts->GetCutHistograms();}return NULL;}
		// Set Options

		void CountTracks();
		void SetConversionCuts(const TString cut);
		void SetConversionCuts(AliConversionPhotonCuts *cuts){fConversionCuts=cuts;}
		void SetEventCuts(const TString cut);
		void SetEventCuts(AliConvEventCuts *cuts){fEventCuts=cuts;}

		void SetUseOwnXYZCalculation(Bool_t flag){fUseOwnXYZCalculation=flag;}
		void SetUseConstructGamma(Bool_t flag){fUseConstructGamma=flag;}
		void SetUseAODConversionPhoton(Bool_t b){if(b){cout<<"Setting Outputformat to AliAODConversionPhoton "<<endl;}else{cout<<"Setting Outputformat to AliKFConversionPhoton "<<endl;};kUseAODConversionPhoton=b;}
		void SetCreateAODs(Bool_t k=kTRUE){fCreateAOD=k;}
		void SetDeltaAODFilename(TString s){fDeltaAODFilename=s;}
		void SetDeltaAODBranchName(TString string) { fDeltaAODBranchName = string;AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));}
		void RelabelAODs(Bool_t relabel=kTRUE){fRelabelAODs=relabel;}
		Bool_t AreAODsRelabeled(){return fRelabelAODs;}
		void RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate);
		TString GetPeriodName(){return fPeriodName;}
		Int_t GetNumberOfPrimaryTracks(){return fNumberOfPrimaryTracks;}
	
	protected:
		// Reconstruct Gammas
		Bool_t ProcessESDV0s();
		AliKFConversionPhoton *ReconstructV0(AliESDv0* fCurrentV0,Int_t currentV0Index);
		void FillAODOutput();
		void FindDeltaAODBranchName();
		Bool_t GetAODConversionGammas();

		// Getter Functions

		const AliExternalTrackParam *GetExternalTrackParam(AliESDv0 *fCurrentV0,Int_t &tracklabel,Int_t charge);
		const AliExternalTrackParam *GetExternalTrackParamP(AliESDv0 *fCurrentV0,Int_t &tracklabel){return GetExternalTrackParam(fCurrentV0,tracklabel,1);};
		const AliExternalTrackParam *GetExternalTrackParamN(AliESDv0 *fCurrentV0,Int_t &tracklabel){return GetExternalTrackParam(fCurrentV0,tracklabel,-1);};
		AliKFParticle *GetPositiveKFParticle(AliAODv0 *fCurrentV0,Int_t fTrackLabel[2]);
		AliKFParticle *GetNegativeKFParticle(AliAODv0 *fCurrentV0,Int_t fTrackLabel[2]);
		AliKFParticle *GetPositiveKFParticle(AliESDv0 *fCurrentV0,Int_t fTrackLabel[2]);
		AliKFParticle *GetNegativeKFParticle(AliESDv0 *fCurrentV0,Int_t fTrackLabel[2]);

		Bool_t GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3],Double_t dca[2]);
		Bool_t GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2]);
		Double_t GetPsiPair(const AliESDv0* v0, const AliExternalTrackParam *positiveparam,const AliExternalTrackParam *negativeparam) const;

		AliConversionPhotonCuts	*fConversionCuts; 			// Pointer to the ConversionCut Selection
		AliConvEventCuts 	*fEventCuts; 				// Pointer to the ConversionCut Selection
		TClonesArray 		*fConversionGammas; 		// TClonesArray holding the reconstructed photons
		Bool_t 				fUseImprovedVertex; 		// set flag to improve primary vertex estimation by adding photons
		Bool_t 				fUseOwnXYZCalculation; 		//flag that determines if we use our own calculation of xyz (markus)
		Bool_t 				fUseConstructGamma; 		//flag that determines if we use ConstructGamma method from AliKF
		Bool_t 				kUseAODConversionPhoton; 	// set flag to use AOD instead of KF output format for photons
		Bool_t 				fCreateAOD; 				// set flag for AOD creation
		TString     		fDeltaAODBranchName;		// File where Gamma Conv AOD is located, if not in default AOD
		TString 			fDeltaAODFilename; 			// set filename for delta/satellite aod
		Bool_t 				fRelabelAODs; 				//
		Bool_t 				fEventIsSelected;
		Int_t 				fNumberOfPrimaryTracks;	 	// Number of Primary Tracks in AOD or ESD
		TString 			fPeriodName;
		
	private:
		AliV0ReaderV1(AliV0ReaderV1 &original);
		AliV0ReaderV1 &operator=(const AliV0ReaderV1 &ref);


    ClassDef(AliV0ReaderV1, 4)
};

inline void AliV0ReaderV1::SetConversionCuts(const TString cut){
    if(fConversionCuts != NULL){
	delete fConversionCuts;
		fConversionCuts=NULL;
    }
    if(fConversionCuts == NULL){
		fConversionCuts=new AliConversionPhotonCuts("V0ReaderCuts","V0ReaderCuts");
		fConversionCuts->InitializeCutsFromCutString(cut.Data());
    }
}

inline void AliV0ReaderV1::SetEventCuts(const TString cut){
    if(fEventCuts != NULL){
	delete fEventCuts;
		fEventCuts=NULL;
    }
    if(fEventCuts == NULL){
		fEventCuts=new AliConvEventCuts("V0ReaderEventCuts","V0ReaderEventCuts");
		fEventCuts->InitializeCutsFromCutString(cut.Data());
    }
}

#endif
