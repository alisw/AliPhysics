/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *								          *
 * Authors: Svein Lindal, Daniel Lohner	                     		  *
 * Version 1.0								  *
 *	                                                		  *
 *                                                                        *
 * based on: on older version (see aliroot up to v5-04-42-AN)             *
 *           AliV0Reader.cxx                                              *
 *           Authors: Kathrin Koch, Kenneth Aamodt, Ana Marin             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class reconstructing conversion photons from V0s
//---------------------------------------------
////////////////////////////////////////////////

// The AliAODConversionPhotons will return the Position (== ESDTrack->GetID())of the positive (negative) ESDTrack 
// in the Array of ESDTracks stored in the ESDEvent via GetTrackLabelPositive() (GetTrackLabelNagative()).
// Since it is the Position in the Array it is always positive.
// After the conversion ESD --> AOD each AOD track will give you the former Position in the ESDArray via the
// Function AODTrack->GetID(). AODTracks are stored for different TrackParameter (e.g. HybridTracks, TPC only ...). No Standard 
// AODTracks (not copies of AliExternalTrackParam) will be stored with negative values of AODTrack->GetID() (Formula: (ESDTrack->GetID()+1)*-1)
// This leads to the possibility of more than one AOD track with the same negative ID. For that reason all AOD tracks are additionally flaged.
// In this analysis we should therefore only use AODTracks with positive values of AODTrack->GetID().

// If you want to find the AODTrack corresponding to the daugher track of a AliAODConversionPhoton you have to find the AODTrack with
// AODTrack->GetID() == GetTrackLabelPositive() (GetTrackLabelNagative()).

#include <TGeoGlobalMagField.h>

#include "AliV0ReaderV1.h"
#include "AliKFParticle.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "TRandom3.h"
#include "AliGenCocktailEventHeader.h"
#include "TList.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjArray.h"

class iostream;


using namespace std;

ClassImp(AliV0ReaderV1)

//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(const char *name) : AliAnalysisTaskSE(name),
	fConversionCuts(NULL),
	fEventCuts(NULL),
	fConversionGammas(NULL),
	fUseImprovedVertex(kTRUE),
	fUseOwnXYZCalculation(kTRUE),
	fUseConstructGamma(kFALSE),
	kUseAODConversionPhoton(kTRUE),
	fCreateAOD(kFALSE),
	fDeltaAODBranchName("GammaConv"),
	fDeltaAODFilename("AliAODGammaConversion.root"),
	fRelabelAODs(kFALSE),
	fEventIsSelected(kFALSE),
	fNumberOfPrimaryTracks(0),
	fPeriodName(""),
	fUseMassToZero(kTRUE),
	fProduceV0findingEffi(kFALSE),
	fMCPhotonLabelArray(NULL),
	fNMCRecPhotons(0),
	fHistograms(NULL),
	fHistoMCGammaPtvsR(NULL),
	fHistoMCGammaPtvsPhi(NULL),
	fHistoMCGammaPtvsEta(NULL),
	fHistoMCGammaRvsPhi(NULL),
	fHistoMCGammaRvsEta(NULL),
	fHistoMCGammaPhivsEta(NULL),
	fHistoRecMCGammaPtvsR(NULL),
	fHistoRecMCGammaPtvsPhi(NULL),
	fHistoRecMCGammaPtvsEta(NULL),
	fHistoRecMCGammaRvsPhi(NULL),
	fHistoRecMCGammaRvsEta(NULL),
	fHistoRecMCGammaPhivsEta(NULL),
	fHistoRecMCGammaMultiPt(NULL),
	fHistoRecMCGammaMultiPtvsEta(NULL),
	fHistoRecMCGammaMultiR(NULL),
	fHistoRecMCGammaMultiPhi(NULL),
	fStrFoundGammas("")
{
	// Default constructor

	DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliV0ReaderV1::~AliV0ReaderV1()
{
	// default deconstructor

	if(fConversionGammas){
		fConversionGammas->Delete();// Clear Objects
		delete fConversionGammas;
		fConversionGammas=0x0;
	}
}

/*
//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(AliV0ReaderV1 &original) : AliAnalysisTaskSE(original),
fConversionCuts(NULL),
fConversionGammas(NULL),
fUseImprovedVertex(original.fUseImprovedVertex),
fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
fUseConstructGamma(original.fUseConstructGamma),
kUseAODConversionPhoton(original.kUseAODConversionPhoton),
fCreateAOD(original.fCreateAOD),
fDeltaAODBranchName(original.fDeltaAODBranchName),
fDeltaAODFilename(original.fDeltaAODFilename),
fEventIsSelected(original.fEventIsSelected)
{
// Default constructor

DefineInput(0, TChain::Class());
}

//____________________________________________________________
AliV0ReaderV1 &AliV0ReaderV1::operator=(const AliV0ReaderV1 &ref){
//
// Assignment operator
// Only copies pointers, object is not the owner of the references
//
if(this != &ref){
AliAnalysisTaskSE::operator=(ref);
fUseImprovedVertex=ref.fUseImprovedVertex;
fUseOwnXYZCalculation=ref.fUseOwnXYZCalculation;
fUseConstructGamma=ref.fUseConstructGamma;
kUseAODConversionPhoton=ref.kUseAODConversionPhoton;
fCreateAOD=ref.fCreateAOD;
fDeltaAODBranchName=ref.fDeltaAODBranchName;
fDeltaAODFilename=ref.fDeltaAODFilename;
fEventIsSelected=ref.fEventIsSelected;
}
return *this;
}
*/

//________________________________________________________________________
void AliV0ReaderV1::Init()
{
	// Initialize function to be called once before analysis
	if(fConversionCuts==NULL){
		if(fConversionCuts==NULL)AliError("No Conversion Cut Selection initialized");
	}
	if(fEventCuts==NULL){
		if(fEventCuts==NULL)AliError("No Event Cut Selection initialized");
	}

	if(fCreateAOD){kUseAODConversionPhoton=kTRUE;}

	if(fConversionGammas != NULL){
		delete fConversionGammas;
		fConversionGammas=NULL;
	}

	if(fConversionGammas == NULL){
		if(kUseAODConversionPhoton){
			fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);}
		else{
			fConversionGammas = new TClonesArray("AliKFConversionPhoton",100);}
	}
	fConversionGammas->Delete();//Reset the TClonesArray
}

//________________________________________________________________________
void AliV0ReaderV1::UserCreateOutputObjects()
{
	// Create AODs

	if(fCreateAOD){
		if (fEventCuts){
			fDeltaAODBranchName.Append("_");
			fDeltaAODBranchName.Append(fEventCuts->GetCutNumber());	
		}	
		if(fConversionCuts){
			fDeltaAODBranchName.Append("_");
			fDeltaAODBranchName.Append(fConversionCuts->GetCutNumber());
			fDeltaAODBranchName.Append("_gamma");
		}
		fConversionGammas->SetName(fDeltaAODBranchName.Data());

		AddAODBranch("TClonesArray", &fConversionGammas, fDeltaAODFilename.Data());
		AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFilename.Data());
	}
	
	
	if (fProduceV0findingEffi){
		TH1::AddDirectory(kFALSE);
		if(fHistograms != NULL){
			delete fHistograms;
			fHistograms 				= NULL;
		}
		if(fHistograms==NULL){
			fHistograms 				= new TList();
			fHistograms->SetOwner(kTRUE);
			fHistograms->SetName(Form("V0FindingEfficiencyInput_%s",fEventCuts->GetCutNumber().Data()));
		}

		fHistoMCGammaPtvsR 				= new TH2F("MCconvGamma_Pt_R","MC converted gamma Pt vs R (|eta| < 0.9)",250,0.0,25,400,0,200);
		fHistoMCGammaPtvsR->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoMCGammaPtvsR->SetYTitle("R_{MC,conv} (cm)");
		fHistograms->Add(fHistoMCGammaPtvsR);
		fHistoMCGammaPtvsEta 			= new TH2F("MCconvGamma_Pt_Eta","MC converted gamma Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
		fHistoMCGammaPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoMCGammaPtvsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoMCGammaPtvsEta);
		fHistoMCGammaPtvsPhi 			= new TH2F("MCconvGamma_Pt_Phi","MC converted gamma Pt vs Phi (|eta| < 0.9) ",250,0.0,25,400,0,2*TMath::Pi());
		fHistoMCGammaPtvsPhi->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoMCGammaPtvsPhi->SetYTitle("#varphi_{MC} (rad)");		
		fHistograms->Add(fHistoMCGammaPtvsPhi);
		fHistoMCGammaRvsPhi 			= new TH2F("MCconvGamma_R_Phi","MC converted gamma R vs Phi (|eta| < 0.9) ",400,0,200,400,0,2*TMath::Pi());
		fHistoMCGammaRvsPhi->SetXTitle("R_{MC,conv} (cm)");
		fHistoMCGammaRvsPhi->SetYTitle("#varphi_{MC} (rad)");		
		fHistograms->Add(fHistoMCGammaRvsPhi);
		fHistoMCGammaRvsEta 			= new TH2F("MCconvGamma_R_Eta","MC converted gamma R vs Eta ",400,0,200,280,-1.4,1.4);
		fHistoMCGammaRvsEta->SetXTitle("R_{MC,conv} (cm)");
		fHistoMCGammaRvsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoMCGammaRvsEta);
		fHistoMCGammaPhivsEta 			= new TH2F("MCconvGamma_Phi_Eta","MC converted gamma Phi vs Eta ",400,0,2*TMath::Pi(),280,-1.4,1.4);
		fHistoMCGammaPhivsEta->SetXTitle("#phi_{MC} (rad)");
		fHistoMCGammaPhivsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoMCGammaPhivsEta);
		fHistoRecMCGammaPtvsR 			= new TH2F("RecMCconvGamma_Pt_R","rec MC converted gamma Pt vs R (|eta| < 0.9)",250,0.0,25,400,0,200);
		fHistoRecMCGammaPtvsR->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoRecMCGammaPtvsR->SetYTitle("R_{MC,conv} (cm)");
		fHistograms->Add(fHistoRecMCGammaPtvsR);
		fHistoRecMCGammaPtvsEta 		= new TH2F("RecMCconvGamma_Pt_Eta","rec MC converted gamma Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
		fHistoRecMCGammaPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoRecMCGammaPtvsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoRecMCGammaPtvsEta);
		fHistoRecMCGammaPtvsPhi 		= new TH2F("RecMCconvGamma_Pt_Phi","rec MC converted gamma Pt vs Phi (|eta| < 0.9) ",250,0.0,25,400,0,2*TMath::Pi());
		fHistoRecMCGammaPtvsPhi->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoRecMCGammaPtvsPhi->SetYTitle("#varphi_{MC} (rad)");		
		fHistograms->Add(fHistoRecMCGammaPtvsPhi);
		fHistoRecMCGammaRvsPhi 			= new TH2F("RecMCconvGamma_R_Phi","rec MC converted gamma R vs Phi (|eta| < 0.9) ",400,0,200,400,0,2*TMath::Pi());
		fHistoRecMCGammaRvsPhi->SetXTitle("R_{MC,conv} (cm)");
		fHistoRecMCGammaRvsPhi->SetYTitle("#varphi_{MC} (rad)");		
		fHistograms->Add(fHistoRecMCGammaRvsPhi);
		fHistoRecMCGammaRvsEta			= new TH2F("RecMCconvGamma_R_Eta","rec MC converted gamma R vs Eta ",400,0,200,280,-1.4,1.4);
		fHistoRecMCGammaRvsEta->SetXTitle("R_{MC,conv} (cm)");
		fHistoRecMCGammaRvsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoRecMCGammaRvsEta);
		fHistoRecMCGammaPhivsEta 		= new TH2F("RecMCconvGamma_Phi_Eta","rec MC converted gamma Phi vs Eta ",400,0,2*TMath::Pi(),280,-1.4,1.4);
		fHistoRecMCGammaPhivsEta->SetXTitle("#phi_{MC} (rad)");
		fHistoRecMCGammaPhivsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoRecMCGammaPhivsEta);
		
		fHistoRecMCGammaMultiPtvsEta 	= new TH2F("RecMCconvGammaMulti_Pt_Eta","rec MC converted gamma (at least double counted) Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
		fHistoRecMCGammaMultiPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
		fHistoRecMCGammaMultiPtvsEta->SetYTitle("#eta_{MC}");		
		fHistograms->Add(fHistoRecMCGammaMultiPtvsEta);

		fHistoRecMCGammaMultiPt			= new TH1F("RecMCconvGammaMulti_Pt","rec MC converted gamma (at least double counted) Pt (|eta| < 0.9)",250,0.0,25);
		fHistoRecMCGammaMultiPt->SetXTitle("p_{MC,T} (GeV/c)");
		fHistograms->Add(fHistoRecMCGammaMultiPt);
		fHistoRecMCGammaMultiR			= new TH1F("RecMCconvGammaMulti_R","rec MC converted gamma (at least double counted) R (|eta| < 0.9)",400,0,200);
		fHistoRecMCGammaMultiR->SetXTitle("R_{MC,conv} (cm)");
		fHistograms->Add(fHistoRecMCGammaMultiR);
		fHistoRecMCGammaMultiPhi		= new TH1F("RecMCconvGammaMulti_Phi","rec MC converted gamma (at least double counted) Phi (|eta| < 0.9)",400,0,2*TMath::Pi());
		fHistoRecMCGammaMultiPhi->SetXTitle("#phi_{MC} (rad)");
		fHistograms->Add(fHistoRecMCGammaMultiPhi);
		
	}	

}
//________________________________________________________________________
Bool_t AliV0ReaderV1::Notify()
{
	if (fPeriodName.CompareTo("") == 0){
		AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
		if(man) {
			AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
			if (inputHandler){
				TTree* tree = (TTree*) inputHandler->GetTree();
				TFile* file = (TFile*) tree->GetCurrentFile();
				TString fileName(file->GetName());
				TObjArray *arr = fileName.Tokenize("/");
				for (Int_t i = 0; i < arr->GetEntriesFast();i++ ){
					TObjString* testObjString = (TObjString*)arr->At(i);
					if (testObjString->GetString().Contains("LHC")){
						fPeriodName = testObjString->GetString();
						i = arr->GetEntriesFast();
					}
				}
			}
		}
		if(!fEventCuts->GetDoEtaShift()) return kTRUE; // No Eta Shift requested, continue
		if(fEventCuts->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			fEventCuts->GetCorrectEtaShiftFromPeriod(fPeriodName);
			fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
			return kTRUE;
		} else {
			printf(" Gamma Conversion Reader %s_%s :: Eta Shift Manually Set to %f \n\n",
					(fEventCuts->GetCutNumber()).Data(),(fConversionCuts->GetCutNumber()).Data(),fEventCuts->GetEtaShift());
			fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
		}
	}
	return kTRUE;
}
//________________________________________________________________________
void AliV0ReaderV1::UserExec(Option_t *option){

	AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
	if(esdEvent) {
		if (!TGeoGlobalMagField::Instance()->GetField()) esdEvent->InitMagneticField(); 
	}

	// Check if correctly initialized
	if(!fConversionGammas)Init();
	
	// User Exec
	fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);

	// AOD Output
	FillAODOutput();
	
}

//________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{

	//Reset the TClonesArray
	fConversionGammas->Delete();

	fInputEvent=inputEvent;
	fMCEvent=mcEvent;

	if(!fInputEvent){
		AliError("No Input event");
		return kFALSE;
	}
	if(!fEventCuts){AliError("No EventCuts");return kFALSE;}
	if(!fConversionCuts){AliError("No ConversionCuts");return kFALSE;}

	
	// Count Primary Tracks Event
	CountTracks();

	// Event Cuts
	if(!fEventCuts->EventIsSelected(fInputEvent,fMCEvent))return kFALSE;
	
	// Set Magnetic Field
	AliKFParticle::SetField(fInputEvent->GetMagneticField());
	
	if(fInputEvent->IsA()==AliAODEvent::Class() && fProduceV0findingEffi){
		fProduceV0findingEffi = kFALSE;
		AliWarning("V0finding effi cannot be run on AODs ");
	}	
	
	if(fProduceV0findingEffi){
		CreatePureMCHistosForV0FinderEffiESD();
		fStrFoundGammas = "";
	}	


	if(fInputEvent->IsA()==AliESDEvent::Class()){
		ProcessESDV0s();
	}
	if(fInputEvent->IsA()==AliAODEvent::Class()){
		GetAODConversionGammas();
	}

	return kTRUE;
}
///________________________________________________________________________
void AliV0ReaderV1::FillAODOutput()
{
   // Fill AOD Output with reconstructed Photons

	if(fInputEvent->IsA()==AliESDEvent::Class()){
		///Make sure delta aod is filled if standard aod is filled (for synchronization when reading aod with standard aod)
		if(fCreateAOD) {
			AliAODHandler * aodhandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
			if (aodhandler && aodhandler->GetFillAOD()) {
				AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
				//PostData(0, fConversionGammas);

			}
		}
	}
}

///________________________________________________________________________
const AliExternalTrackParam *AliV0ReaderV1::GetExternalTrackParam(AliESDv0 *fCurrentV0,Int_t &tracklabel,Int_t charge){

	// Get External Track Parameter with given charge

	if(!(charge==1||charge==-1)){AliError("Charge not defined");return 0x0;}

	// Check for sign flip
	if(fCurrentV0){
		if(!fCurrentV0->GetParamN()||!fCurrentV0->GetParamP())return 0x0;
		if(!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex())||!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))return 0x0;
		if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))->Charge()==charge){
			tracklabel=fCurrentV0->GetPindex();
			return fCurrentV0->GetParamP();}
		if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex()))->Charge()==charge){
			tracklabel=fCurrentV0->GetNindex();
			return fCurrentV0->GetParamN();}
	}
	return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessESDV0s()
{
	// Process ESD V0s for conversion photon reconstruction
	AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);

	AliKFConversionPhoton *fCurrentMotherKFCandidate=NULL;

	if(fESDEvent){

		for(Int_t currentV0Index=0;currentV0Index<fESDEvent->GetNumberOfV0s();currentV0Index++){
			AliESDv0 *fCurrentV0=(AliESDv0*)(fESDEvent->GetV0(currentV0Index));
			if(!fCurrentV0){
				printf("Requested V0 does not exist");
				continue;
			}

			fCurrentMotherKFCandidate=ReconstructV0(fCurrentV0,currentV0Index);

			if(fCurrentMotherKFCandidate){
				// Add Gamma to the TClonesArray

				if(kUseAODConversionPhoton){
					new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(fCurrentMotherKFCandidate);
					AliAODConversionPhoton * currentConversionPhoton = (AliAODConversionPhoton*)(fConversionGammas->At(fConversionGammas->GetEntriesFast()-1));
					currentConversionPhoton->SetMass(fCurrentMotherKFCandidate->M());
					if (fUseMassToZero) currentConversionPhoton->SetMassToZero();
				} else {
					new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliKFConversionPhoton(*fCurrentMotherKFCandidate);
				}

				delete fCurrentMotherKFCandidate;
				fCurrentMotherKFCandidate=NULL;
			}
		}
	}
	return kTRUE;
}

///________________________________________________________________________
AliKFConversionPhoton *AliV0ReaderV1::ReconstructV0(AliESDv0 *fCurrentV0,Int_t currentV0Index)
{
	//   cout << currentV0Index << endl;
	// Reconstruct conversion photon from ESD v0
	fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonIn);

	//checks if on the fly mode is set
	if(!fConversionCuts->SelectV0Finder(fCurrentV0->GetOnFlyStatus())){
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kOnFly);
		return 0x0;
	}

	if (fMCEvent && fProduceV0findingEffi ) FillRecMCHistosForV0FinderEffiESD(fCurrentV0);
	
	// TrackLabels
	Int_t currentTrackLabels[2]={-1,-1};

	// Get Daughter KF Particles

	const AliExternalTrackParam *fCurrentExternalTrackParamPositive=GetExternalTrackParamP(fCurrentV0,currentTrackLabels[0]);
	//    cout << fCurrentExternalTrackParamPositive << "\t" << currentTrackLabels[0] << endl;
	const AliExternalTrackParam *fCurrentExternalTrackParamNegative=GetExternalTrackParamN(fCurrentV0,currentTrackLabels[1]);
	//    cout << fCurrentExternalTrackParamNegative << "\t" << currentTrackLabels[1] << endl;
	if(!fCurrentExternalTrackParamPositive||!fCurrentExternalTrackParamNegative)return 0x0;

	// Apply some Cuts before Reconstruction

	
	
	AliVTrack * posTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[0]);
	AliVTrack * negTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[1]);
	if(!negTrack || !posTrack) {
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kNoTracks);
		return 0x0;
	}
	// Track Cuts
	if(!fConversionCuts->TracksAreSelected(negTrack, posTrack)){
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kTrackCuts);
		return 0x0;
	}
	
	fConversionCuts->FillV0EtaBeforedEdxCuts(fCurrentV0->Eta());
	if (!fConversionCuts->dEdxCuts(posTrack)) {
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kdEdxCuts);
		return 0x0;
	}
	// PID Cuts
	if(!fConversionCuts->dEdxCuts(negTrack)) {
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kdEdxCuts);
		return 0x0;
	}
	fConversionCuts->FillV0EtaAfterdEdxCuts(fCurrentV0->Eta());
	// Reconstruct Photon
	AliKFConversionPhoton *fCurrentMotherKF=NULL;
	//    fUseConstructGamma = kFALSE;
	//    cout << "construct gamma " << endl;
	AliKFParticle fCurrentNegativeKFParticle(*(fCurrentExternalTrackParamNegative),11);
	//    cout << fCurrentExternalTrackParamNegative << "\t" << endl;
	AliKFParticle fCurrentPositiveKFParticle(*(fCurrentExternalTrackParamPositive),-11);
	//    cout << fCurrentExternalTrackParamPositive << "\t"  << endl;
	//    cout << currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
	//    cout << "construct gamma " <<fUseConstructGamma << endl;
	
	// Reconstruct Gamma
	if(fUseConstructGamma){
		fCurrentMotherKF = new AliKFConversionPhoton();
		fCurrentMotherKF->ConstructGamma(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
	}else{
		fCurrentMotherKF = new AliKFConversionPhoton(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
		fCurrentMotherKF->SetMassConstraint(0,0.0001);
	}

	// Set Track Labels

	fCurrentMotherKF->SetTrackLabels(currentTrackLabels[0],currentTrackLabels[1]);

	// Set V0 index

	fCurrentMotherKF->SetV0Index(currentV0Index);

	//Set MC Label
	if(fMCEvent){

		AliStack *fMCStack= fMCEvent->Stack();

		Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelPositive())->GetLabel());
		Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelNegative())->GetLabel());

// 		cout << "rec: " <<  currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
// 		cout << "recProp: " <<  fCurrentMotherKF->GetTrackLabelPositive() << "\t" << fCurrentMotherKF->GetTrackLabelNegative() << endl;
// 		cout << "MC: " <<  labeln << "\t" << labelp << endl;
		
		TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
		TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);

		if(fPositiveMCParticle&&fNegativeMCParticle){
			fCurrentMotherKF->SetMCLabelPositive(labelp);
			fCurrentMotherKF->SetMCLabelNegative(labeln);
		}
	}

	// Update Vertex (moved for same eta compared to old)
	//      cout << currentV0Index <<" \t before: \t" << fCurrentMotherKF->GetPx() << "\t" << fCurrentMotherKF->GetPy() << "\t" << fCurrentMotherKF->GetPz()  << endl;
	if(fUseImprovedVertex == kTRUE){
		AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
		//        cout << "Prim Vtx: " << primaryVertexImproved.GetX() << "\t" << primaryVertexImproved.GetY() << "\t" << primaryVertexImproved.GetZ() << endl;
		primaryVertexImproved+=*fCurrentMotherKF;
		fCurrentMotherKF->SetProductionVertex(primaryVertexImproved);
	}
	// SetPsiPair

	Double_t PsiPair=GetPsiPair(fCurrentV0,fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative);
	fCurrentMotherKF->SetPsiPair(PsiPair);
	
	// Recalculate ConversionPoint
	Double_t dca[2]={0,0};
	if(fUseOwnXYZCalculation){
		Double_t convpos[3]={0,0,0};
		if(!GetConversionPoint(fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative,convpos,dca)){
			fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kConvPointFail);
			delete fCurrentMotherKF;
			fCurrentMotherKF=NULL;
			return 0x0;
		}

		fCurrentMotherKF->SetConversionPoint(convpos);
	}

	if(fCurrentMotherKF->GetNDF() > 0.)
		fCurrentMotherKF->SetChi2perNDF(fCurrentMotherKF->GetChi2()/fCurrentMotherKF->GetNDF());   //->Photon is created before all chi2 relevant changes are performed, set it "by hand"


	// Set Dilepton Mass (moved down for same eta compared to old)
	fCurrentMotherKF->SetMass(fCurrentMotherKF->M());

	// Apply Photon Cuts

	if(!fConversionCuts->PhotonCuts(fCurrentMotherKF,fInputEvent)){
		fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonCuts);
		delete fCurrentMotherKF;
		fCurrentMotherKF=NULL;
		return 0x0;
	}

	//    cout << currentV0Index <<" \t after: \t" <<fCurrentMotherKF->GetPx() << "\t" << fCurrentMotherKF->GetPy() << "\t" << fCurrentMotherKF->GetPz()  << endl;

	fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonOut);
	return fCurrentMotherKF;
}

///________________________________________________________________________
Double_t AliV0ReaderV1::GetPsiPair(const AliESDv0* v0, const AliExternalTrackParam *positiveparam,const AliExternalTrackParam *negativeparam) const {
	//
	// Angle between daughter momentum plane and plane
	//

	AliExternalTrackParam nt(*negativeparam);
	AliExternalTrackParam pt(*positiveparam);

	Float_t magField = fInputEvent->GetMagneticField();

	Double_t xyz[3] = {0.,0.,0.};
	v0->GetXYZ(xyz[0],xyz[1],xyz[2]);

	// Double_t pPlus[3]  = {pt.Px(),pt.Py(),pt.Pz()};
	// Double_t pMinus[3] = {nt.Px(),nt.Py(),nt.Pz()};

	// Double_t u[3] = {pPlus[0]+pMinus[0],pPlus[1]+pMinus[1],pPlus[2]+pMinus[2]};
	// Double_t normu = sqrt( (u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]) );

	// u[0] = u[0] / normu;
	// u[1] = u[1] / normu;
	// u[2] = u[2] / normu;

	// Double_t normpPlus = sqrt( (pPlus[0]*pPlus[0]) + (pPlus[1]*pPlus[1]) + (pPlus[2]*pPlus[2]) );
	// Double_t normpMinus = sqrt( (pMinus[0]*pMinus[0]) + (pMinus[1]*pMinus[1]) + (pMinus[2]*pMinus[2]) );

	// pPlus[0] = pPlus[0] / normpPlus;
	// pPlus[1] = pPlus[1] / normpPlus;
	// pPlus[2] = pPlus[2] / normpPlus;

	// pMinus[0] = pMinus[0] / normpMinus;
	// pMinus[1] = pMinus[1] / normpMinus;
	// pMinus[2] = pMinus[2] / normpMinus;

	// Double_t v[3] = {0,0,0}; // pPlus X pMinus
	// v[0] = (pPlus[1]*pMinus[2]) - (pPlus[2]*pMinus[1]);
	// v[1] = (pPlus[2]*pMinus[0]) - (pPlus[0]*pMinus[2]);
	// v[2] = (pPlus[0]*pMinus[1]) - (pPlus[1]*pMinus[0]);

	// Double_t w[3] = {0,0,0}; // u X v
	// w[0] = (u[1]*v[2]) - (u[2]*v[1]);
	// w[1] = (u[2]*v[0]) - (u[0]*v[2]);
	// w[2] = (u[0]*v[1]) - (u[1]*v[0]);

	// Double_t z[3] = {0,0,1};
	// Double_t wc[3] = {0,0,0}; // u X v
	// wc[0] = (u[1]*z[2]) - (u[2]*z[1]);
	// wc[1] = (u[2]*z[0]) - (u[0]*z[2]);
	// wc[2] = (u[0]*z[1]) - (u[1]*z[0]);

	// Double_t PhiV = TMath::ACos((w[0]*wc[0]) + (w[1]*wc[1]) + (w[2]*wc[2]));
	//return abs(PhiV);


	// TVector3 pPlus(pt.Px(),pt.Py(),pt.Pz());
	// TVector3 pMinus(nt.Px(),nt.Py(),nt.Pz());

	// TVector3 u = pMinus + pPlus;
	// u = u*(1/u.Mag());

	// TVector3 pHPlus = pPlus*(1/pPlus.Mag());
	// TVector3 pHMinus = pMinus*(1/pMinus.Mag());

	// TVector3 v = pHPlus.Cross(pHMinus);
	// TVector3 w = u.Cross(v);
	// TVector3 z(0,0,1);
	// TVector3 wc = u.Cross(z);

	// Double_t PhiV = w * wc;

	Double_t mn[3] = {0,0,0};
	Double_t mp[3] = {0,0,0};

	v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
	v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

	Double_t deltat = 1.;
	deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
	Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated

	Double_t momPosProp[3] = {0,0,0};
	Double_t momNegProp[3] = {0,0,0};

	Double_t psiPair = 4.;
	if(nt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
	if(pt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency

	pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
	nt.GetPxPyPz(momNegProp);

	Double_t pEle =
		TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter

	Double_t pPos =
		TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter

	Double_t scalarproduct =
		momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta

	Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

	//    psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));
	psiPair = TMath::ASin(deltat/chipair);
	return psiPair;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2]){

	// Get Center of the helix track parametrization

	Int_t charge=track->Charge();
	Double_t b=fInputEvent->GetMagneticField();

	Double_t	helix[6];
	track->GetHelixParameters(helix,b);

	Double_t xpos =	helix[5];
	Double_t ypos =	helix[0];
	Double_t radius = TMath::Abs(1./helix[4]);
	Double_t phi = helix[2];

	if(phi < 0){
		phi = phi + 2*TMath::Pi();
	}

	phi -= TMath::Pi()/2.;
	Double_t xpoint =	radius * TMath::Cos(phi);
	Double_t ypoint =	radius * TMath::Sin(phi);

	if(b<0){
		if(charge > 0){
			xpoint = - xpoint;
			ypoint = - ypoint;
		}
	}
	if(b>0){
		if(charge < 0){
			xpoint = - xpoint;
			ypoint = - ypoint;
		}
	}
	center[0] =	xpos + xpoint;
	center[1] =	ypos + ypoint;

	return 1;
}
///________________________________________________________________________
Bool_t AliV0ReaderV1::GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3],Double_t dca[2]){

	// Recalculate Conversion Point

	if(!pparam||!nparam)return kFALSE;

	Double_t helixcenterpos[2];
	GetHelixCenter(pparam,helixcenterpos);

	Double_t helixcenterneg[2];
	GetHelixCenter(nparam,helixcenterneg);

	Double_t helixpos[6];
	pparam->GetHelixParameters(helixpos,fInputEvent->GetMagneticField());
	Double_t posradius = TMath::Abs(1./helixpos[4]);

	Double_t helixneg[6];
	nparam->GetHelixParameters(helixneg,fInputEvent->GetMagneticField());
	Double_t negradius = TMath::Abs(1./helixneg[4]);

	// Calculate xy-position

	Double_t xpos = helixcenterpos[0];
	Double_t ypos = helixcenterpos[1];
	Double_t xneg = helixcenterneg[0];
	Double_t yneg = helixcenterneg[1];

	convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
	convpos[1] = (ypos*negradius+ yneg*posradius)/(negradius+posradius);


	// Calculate Track XY vertex position

	Double_t deltaXPos = convpos[0] -	xpos;
	Double_t deltaYPos = convpos[1] -	ypos;

	Double_t deltaXNeg = convpos[0] -	xneg;
	Double_t deltaYNeg = convpos[1] -	yneg;

	Double_t alphaPos =	TMath::Pi() + TMath::ATan2(-deltaYPos,-deltaXPos);
	Double_t alphaNeg =	TMath::Pi() + TMath::ATan2(-deltaYNeg,-deltaXNeg);

	Double_t vertexXNeg =	xneg +	TMath::Abs(negradius)*TMath::Cos(alphaNeg);
	Double_t vertexYNeg =	yneg +	TMath::Abs(negradius)*TMath::Sin(alphaNeg);

	Double_t vertexXPos =	xpos +	TMath::Abs(posradius)*TMath::Cos(alphaPos);
	Double_t vertexYPos =	ypos +	TMath::Abs(posradius)*TMath::Sin(alphaPos);

	AliExternalTrackParam p(*pparam);
	AliExternalTrackParam n(*nparam);

	TVector2 vertexPos(vertexXPos,vertexYPos);
	TVector2 vertexNeg(vertexXNeg,vertexYNeg);

	// Convert to local coordinate system
	TVector2 vertexPosRot=vertexPos.Rotate(-p.GetAlpha());
	TVector2 vertexNegRot=vertexNeg.Rotate(-n.GetAlpha());

	// Propagate Track Params to Vertex

	if(!p.PropagateTo(vertexPosRot.X(),fInputEvent->GetMagneticField()))return kFALSE;
	if(!n.PropagateTo(vertexNegRot.X(),fInputEvent->GetMagneticField()))return kFALSE;

	// Check whether propagation was sucessful

	if(TMath::Abs(vertexPos.Mod()-TMath::Sqrt(p.GetX()*p.GetX()+p.GetY()*p.GetY()))>0.01)return kFALSE;
	if(TMath::Abs(vertexNeg.Mod()-TMath::Sqrt(n.GetX()*n.GetX()+n.GetY()*n.GetY()))>0.01)return kFALSE;

	// Calculate z position

	convpos[2] = (p.GetZ()*negradius+n.GetZ()*posradius)/(negradius+posradius);

	// Calculate DCA
	TVector2 vdca=vertexPos-vertexNeg;
	dca[0]=vdca.Mod();
	dca[1]=TMath::Abs(n.GetZ()-p.GetZ());

	return kTRUE;
	}
	//________________________________________________________________________
	Bool_t AliV0ReaderV1::GetAODConversionGammas(){

	// Get reconstructed conversion photons from satellite AOD file

	AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);

	if(fAODEvent){

		if(fConversionGammas == NULL){
			fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);
		}
		fConversionGammas->Delete();//Reset the TClonesArray

		//Get Gammas from satellite AOD gamma branch

		AliAODConversionPhoton *gamma=0x0;

		TClonesArray *fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));
		
		if(!fInputGammas){
			FindDeltaAODBranchName();
			fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));}
		if(!fInputGammas){AliError("No Gamma Satellites found");return kFALSE;}
		// Apply Selection Cuts to Gammas and create local working copy
		if(fInputGammas){
			for(Int_t i=0;i<fInputGammas->GetEntriesFast();i++){
				gamma=dynamic_cast<AliAODConversionPhoton*>(fInputGammas->At(i));
				if(gamma){
				if(fRelabelAODs)RelabelAODPhotonCandidates(gamma);
				if(fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
					new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(*gamma);}
				}
			}
		}
	}

	if(fConversionGammas->GetEntries()){return kTRUE;}

	return kFALSE;
}

//________________________________________________________________________
void AliV0ReaderV1::FindDeltaAODBranchName(){

	// Find delta AOD branchname containing reconstructed photons

	TList *list=fInputEvent->GetList();
	for(Int_t ii=0;ii<list->GetEntries();ii++){
		TString name((list->At(ii))->GetName());
		if(name.BeginsWith(fDeltaAODBranchName)&&name.EndsWith("gamma")){
			fDeltaAODBranchName=name;
			AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
			return;
		}
	}
}
//________________________________________________________________________
void AliV0ReaderV1::RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate){

	// Relabeling For AOD Event
	// ESDiD -> AODiD
	// MCLabel -> AODMCLabel
	Bool_t AODLabelPos = kFALSE;
	Bool_t AODLabelNeg = kFALSE;

	for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
		AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
		if(!AODLabelPos){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
				PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
				PhotonCandidate->SetLabelPositive(i);
				AODLabelPos = kTRUE;
			}
		}
		if(!AODLabelNeg){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
				PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
				PhotonCandidate->SetLabelNegative(i);
				AODLabelNeg = kTRUE;
			}
		}
		if(AODLabelNeg && AODLabelPos){
			return;
		}
	}
	if(!AODLabelPos || !AODLabelNeg){
		AliError(Form("NO AOD Daughters Found Pos: %i %i Neg: %i %i",AODLabelPos,PhotonCandidate->GetTrackLabelPositive(),AODLabelNeg,PhotonCandidate->GetTrackLabelNegative()));
	}

}

//************************************************************************
// This function counts the number of primary tracks in the event, the 
// implementation for ESD and AOD had to be different due to the different
// filters which are already applied on AOD level the so-called filter 
// bits, we tried to replicate the conditions in both but an exact match 
// could not be reached. The cuts are similar to the ones used by the 
// jet-group
//************************************************************************
//________________________________________________________________________
void AliV0ReaderV1::CountTracks(){

	// In Case of MC count only MB tracks
	// if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
	//    fEventCuts->GetNotRejectedParticles(1,NULL,fMCEvent);   
	// }
	// if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
	//    fEventCuts->GetNotRejectedParticles(1,NULL,fInputEvent);   
	// }
	
	if(fInputEvent->IsA()==AliESDEvent::Class()){
		// Using standard function for setting Cuts
		Bool_t selectPrimaries=kTRUE;
		AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
		EsdTrackCuts->SetMaxDCAToVertexZ(2);
		EsdTrackCuts->SetEtaRange(-0.8, 0.8);
		EsdTrackCuts->SetPtRange(0.15);
		fNumberOfPrimaryTracks = 0;
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(!EsdTrackCuts->AcceptTrack(curTrack)) continue;
			//if(fMCEvent && !(fEventCuts->IsParticleFromBGEvent(abs(curTrack->GetLabel()),fMCEvent->Stack(),fInputEvent))) continue;
			fNumberOfPrimaryTracks++;
		}
		delete EsdTrackCuts;
		EsdTrackCuts=0x0;
	}
	else if(fInputEvent->IsA()==AliAODEvent::Class()){
		fNumberOfPrimaryTracks = 0;
		for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
			AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
			if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
			if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
			if(abs(curTrack->Eta())>0.8) continue;
			if(curTrack->Pt()<0.15) continue;
			//if(fMCEvent && !(fEventCuts->IsParticleFromBGEvent(abs(curTrack->GetLabel()),NULL,fInputEvent))) continue;
			//if(abs(curTrack->ZAtDCA())>2) continue; // Only Set For TPCOnly tracks
			fNumberOfPrimaryTracks++;
		}
	}

	return;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ParticleIsConvertedPhoton(AliStack *MCStack, TParticle *particle, Double_t etaMax, Double_t rMax, Double_t zMax){
	// MonteCarlo Photon Selection
	if(!MCStack)return kFALSE;
	
	if (particle->GetPdgCode() == 22){
		// check whether particle is within eta range
		if( abs(particle->Eta()) > etaMax ) return kFALSE;
		// check if particle doesn't have a photon as mother
		if(particle->GetMother(0) >-1 && MCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
			return kFALSE; // no photon as mothers!
		}
		// looking for conversion gammas (electron + positron from pairbuilding (= 5) )
		TParticle* ePos = NULL;
		TParticle* eNeg = NULL;
		if(particle->GetNDaughters() >= 2){
			for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
				TParticle *tmpDaughter = MCStack->Particle(daughterIndex);
				if(tmpDaughter->GetUniqueID() == 5){
					if(tmpDaughter->GetPdgCode() == 11){
						eNeg = tmpDaughter;
					} else if(tmpDaughter->GetPdgCode() == -11){
						ePos = tmpDaughter;
					}
				}
			}
		}
		if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
			return kFALSE;
		}
		// check if electrons are in correct eta window
		if( abs(ePos->Eta()) > etaMax ||
			abs(eNeg->Eta()) > etaMax )
			return kFALSE;

		// check if photons have converted in reconstructable range
		if(ePos->R() > rMax){
			return kFALSE; // cuts on distance from collision point
		} 
		if(abs(ePos->Vz()) > zMax){
			return kFALSE;  // outside material
		}
		if(abs(eNeg->Vz()) > zMax){
			return kFALSE;  // outside material
		}

		
		Double_t lineCutZRSlope = tan(2*atan(exp(-etaMax)));
		Double_t lineCutZValue = 7.;
		if( ePos->R() <= ((abs(ePos->Vz()) * lineCutZRSlope) - lineCutZValue)){
			return kFALSE;  // line cut to exclude regions where we do not reconstruct
		}
		if( eNeg->R() <= ((abs(eNeg->Vz()) * lineCutZRSlope) - lineCutZValue)){
			return kFALSE; // line cut to exclude regions where we do not reconstruct
		}
		if (ePos->Pt() < 0.05 || eNeg->Pt() < 0.05){
			return kFALSE;
		}	

		return kTRUE;
	}
	return kFALSE;
}

///_______________________________________________________________________
void AliV0ReaderV1::CreatePureMCHistosForV0FinderEffiESD(){
	
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
// 	cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;
	
	AliStack *fMCStack= fMCEvent->Stack();	
	// Loop over all primary MC particle	
	for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
		if (fEventCuts->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 
			// fill primary histogram
			TParticle* particle = (TParticle *)fMCStack->Particle(i);
			if (!particle) continue;
			if (ParticleIsConvertedPhoton(fMCStack, particle, 0.9, 180.,250. )){
				TParticle *tmpDaughter = fMCStack->Particle(particle->GetFirstDaughter());
				if (!tmpDaughter) continue;
				fHistoMCGammaPtvsR->Fill(particle->Pt(),tmpDaughter->R());
				fHistoMCGammaPtvsPhi->Fill(particle->Pt(),particle->Phi());
				fHistoMCGammaRvsPhi->Fill(tmpDaughter->R(),particle->Phi());
			}	
			if (ParticleIsConvertedPhoton(fMCStack, particle, 1.4, 180.,250. )){
				TParticle *tmpDaughter = fMCStack->Particle(particle->GetFirstDaughter());
				if (!tmpDaughter) continue;
				fHistoMCGammaPtvsEta->Fill(particle->Pt(),particle->Eta());
				fHistoMCGammaRvsEta->Fill(tmpDaughter->R(),particle->Eta());
				fHistoMCGammaPhivsEta->Fill(particle->Phi(),particle->Eta());
			}	
		}
	}	
}	

///_______________________________________________________________________
void AliV0ReaderV1::FillRecMCHistosForV0FinderEffiESD( AliESDv0* currentV0){
	
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
// 	cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

	Int_t tracklabelPos=currentV0->GetPindex();
	Int_t tracklabelNeg=currentV0->GetNindex();
	
	AliStack *fMCStack= fMCEvent->Stack();

	Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,tracklabelPos)->GetLabel());
	Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,tracklabelNeg)->GetLabel());

	TParticle* negPart = (TParticle *)fMCStack->Particle(labeln);
	TParticle* posPart = (TParticle *)fMCStack->Particle(labelp);
	
	if ( negPart == NULL || posPart == NULL ) return;
// 	if (!(negPart->GetPdgCode() == 11)) return;
// 	if (!(posPart->GetPdgCode() == -11)) return;
	UInt_t motherlabelNeg = negPart->GetFirstMother();
	UInt_t motherlabelPos = posPart->GetFirstMother();
	
// 	cout << "mother neg " << motherlabelNeg << " mother pos " << motherlabelPos << endl;
	if (motherlabelNeg == motherlabelPos && negPart->GetFirstMother() != -1){
		if (fEventCuts->IsConversionPrimaryESD( fMCStack, negPart->GetFirstMother(), mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 
			
			TParticle* mother =  (TParticle *)fMCStack->Particle(motherlabelNeg);
			if (mother->GetPdgCode() == 22 ){
				if (!CheckIfContainedInStringAndAppend(fStrFoundGammas,motherlabelNeg ) ){
					if (ParticleIsConvertedPhoton(fMCStack, mother, 0.9, 180.,250. )){
						fHistoRecMCGammaPtvsR->Fill(mother->Pt(),negPart->R());
						fHistoRecMCGammaPtvsPhi->Fill(mother->Pt(),mother->Phi());
						fHistoRecMCGammaRvsPhi->Fill(negPart->R(),mother->Phi());
					} 
					if (ParticleIsConvertedPhoton(fMCStack, mother, 1.4, 180.,250. )){
						fHistoRecMCGammaPtvsEta->Fill(mother->Pt(),mother->Eta());
						fHistoRecMCGammaRvsEta->Fill(negPart->R(),mother->Eta());
						fHistoRecMCGammaPhivsEta->Fill(mother->Phi(),mother->Eta());
					}	
// 					cout << "new gamma found" << endl;
				} else {
					if (ParticleIsConvertedPhoton(fMCStack, mother, 0.9, 180.,250. )){
						fHistoRecMCGammaMultiPt->Fill(mother->Pt());
						fHistoRecMCGammaMultiPhi->Fill(mother->Phi());
						fHistoRecMCGammaMultiR->Fill(negPart->R());
					} 
					if (ParticleIsConvertedPhoton(fMCStack, mother, 1.4, 180.,250. )){
						fHistoRecMCGammaMultiPtvsEta->Fill(mother->Pt(),mother->Eta());
					}	
// 					cout << "this one I had already: " << motherlabelNeg << endl << "-----------------------"  << endl;
				}	
// 				cout << "event gammas: " << fStrFoundGammas.Data() << endl;
			}	
		}	
	}	
}	

//_________________________________________________________________________________
Bool_t AliV0ReaderV1::CheckIfContainedInString(TString input, Int_t tobechecked){
	TObjArray *arr = input.Tokenize(",");
	for (Int_t i = 0; i < arr->GetEntriesFast();i++){
		TString tempStr = ((TObjString*)arr->At(i))->GetString();
		if (tempStr.Atoi() == tobechecked) return kTRUE;
	}	
	return kFALSE;
}

//_________________________________________________________________________________
Bool_t AliV0ReaderV1::CheckIfContainedInStringAndAppend(TString &input, Int_t tobechecked){
	TObjArray *arr = input.Tokenize(",");
	Bool_t isContained = kFALSE;
	for (Int_t i = 0; i < arr->GetEntriesFast();i++){
		TString tempStr = ((TObjString*)arr->At(i))->GetString();
		if (tempStr.Atoi() == tobechecked) isContained= kTRUE;
	}	
	if (!isContained)input.Append(Form("%i,",tobechecked));	
	return isContained;
}

//________________________________________________________________________
void AliV0ReaderV1::Terminate(Option_t *)
{

}
