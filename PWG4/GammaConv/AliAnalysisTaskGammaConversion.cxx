/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																																				*
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt												*
 * Version 1.1																														*
 *																																				*
 * Permission to use, copy, modify and distribute this software and its	 *
 * documentation strictly for non-commercial purposes is hereby granted	 *
 * without fee, provided that the above copyright notice appears in all	 *
 * copies and that both the copyright notice and this permission notice	 *
 * appear in the supporting documentation. The authors make no claims		 *
 * about the suitability of this software for any purpose. It is					*
 * provided "as is" without express or implied warranty.									*
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
///////////////////////////////////////////////

// root
#include <TChain.h>

// analysis
#include "AliAnalysisTaskGammaConversion.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "TNtuple.h"
//#include "AliCFManager.h"	// for CF
//#include "AliCFContainer.h"	 // for CF
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODPWG4Particle.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODConversionPhoton.h"
#include "AliGammaConversionBGHandler.h"
#include "AliESDCaloCluster.h" // for combining PHOS and GammaConv
#include "AliKFVertex.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEventHeader.h"
#include <AliMCEventHandler.h>
#include "TRandom3.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliKFConversionPhoton.h"
#include "AliKFConversionMother.h"

class AliESDTrackCuts;
class AliCFContainer;
class AliCFManager;
class AliKFVertex;
class AliAODHandler;
class AliAODEvent;
class ALiESDEvent;
class AliMCEvent;
class AliMCEventHandler;
class AliESDInputHandler;
class AliAnalysisManager;
class Riostream;
class TFile;
class TInterpreter;
class TSystem;
class TROOT;

ClassImp(AliAnalysisTaskGammaConversion)


AliAnalysisTaskGammaConversion::AliAnalysisTaskGammaConversion():
AliAnalysisTaskSE(),
	fV0Reader(NULL),
	fStack(NULL),
	fMCTruth(NULL),		// for CF
	fGCMCEvent(NULL),		// for CF
	fESDEvent(NULL),	
	fOutputContainer(NULL),
	fCFManager(0x0),	 // for CF
	fHistograms(NULL),
	fTriggerCINT1B(kFALSE),
	fDoMCTruth(kFALSE),
	fDoNeutralMeson(kFALSE),
	fDoOmegaMeson(kFALSE),
	fDoJet(kFALSE),
	fDoChic(kFALSE),
	fRecalculateV0ForGamma(kFALSE),
	fKFReconstructedGammasTClone(NULL),
	fKFReconstructedPi0sTClone(NULL),
	fKFRecalculatedGammasTClone(NULL),
	fCurrentEventPosElectronTClone(NULL),
	fCurrentEventNegElectronTClone(NULL),
	fKFReconstructedGammasCutTClone(NULL),
	fPreviousEventTLVNegElectronTClone(NULL),
	fPreviousEventTLVPosElectronTClone(NULL),	
       //  fElectronRecalculatedv1(),
       // fElectronRecalculatedv2(),
	fElectronMass(-1),
	fGammaMass(-1),
	fPi0Mass(-1),
	fEtaMass(-1),
	fGammaWidth(-1),
	fPi0Width(-1),
	fEtaWidth(-1),
	fMinOpeningAngleGhostCut(0.),
	fEsdTrackCuts(NULL),
	fCalculateBackground(kFALSE),
	fWriteNtuple(kFALSE),
	fGammaNtuple(NULL),
	fNeutralMesonNtuple(NULL),
	fTotalNumberOfAddedNtupleEntries(0),
	fChargedParticles(NULL),
	fChargedParticlesId(),
	fGammaPtHighest(0.),
	fMinPtForGammaJet(1.),
	fMinIsoConeSize(0.2),
	fMinPtIsoCone(0.7),
	fMinPtGamChargedCorr(0.5),
	fMinPtJetCone(0.5),
	fLeadingChargedIndex(-1),
	fLowPtMapping(1.),
	fHighPtMapping(3.),
	fDoCF(kFALSE),
	fAODGamma(NULL),
	//fAODPi0(NULL),
	//fAODOmega(NULL),
	fAODBranchName("GammaConv"),
	fKFCreateAOD(kTRUE),
	fKFForceAOD(kFALSE),
	fKFDeltaAODFileName(""),
	fDoNeutralMesonV0MCCheck(kFALSE),
	fUseTrackMultiplicityForBG(kTRUE),
	fMoveParticleAccordingToVertex(kFALSE),
	fApplyChi2Cut(kFALSE),
	fNRandomEventsForBG(15),
	fNDegreesPMBackground(15),
	fDoRotation(kTRUE),
	fCheckBGProbability(kTRUE),
	fRemovePileUp(kFALSE),
	fSelectV0AND(kFALSE),
	fTriggerAnalysis(NULL),
	fMultiplicity(0),
	fUseMultiplicity(0), 
	fUseMultiplicityBin(0),
	fUseHBTMultiplicity(0),
	fUseHBTMultiplicityBin(0),
	fUseCentrality(0), 
	fUseCentralityBin(0),
	fRandom(0)
{
	// Default constructor

	/*	 Kenneth: the default constructor should not have any define input/output or the call to SetESDtrackCuts
	// Common I/O in slot 0
	DefineInput (0, TChain::Class());
	DefineOutput(0, TTree::Class());
	
	// Your private output
	DefineOutput(1, TList::Class());
	
	// Define standard ESD track cuts for Gamma-hadron correlation 
	SetESDtrackCuts();
	*/
}

AliAnalysisTaskGammaConversion::AliAnalysisTaskGammaConversion(const char* name):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
	fStack(NULL),
	fMCTruth(NULL),		// for CF
	fGCMCEvent(NULL),		// for CF
	fESDEvent(NULL),	
	fOutputContainer(0x0),
	fCFManager(0x0),	 // for CF
	fHistograms(NULL),
	fTriggerCINT1B(kFALSE),
	fDoMCTruth(kFALSE),
	fDoNeutralMeson(kFALSE),
	fDoOmegaMeson(kFALSE),
	fDoJet(kFALSE),
	fDoChic(kFALSE),
	fRecalculateV0ForGamma(kFALSE),
	fKFReconstructedGammasTClone(NULL),
	fKFReconstructedPi0sTClone(NULL),
	fKFRecalculatedGammasTClone(NULL),
	fCurrentEventPosElectronTClone(NULL),
	fCurrentEventNegElectronTClone(NULL),
	fKFReconstructedGammasCutTClone(NULL),
	fPreviousEventTLVNegElectronTClone(NULL),
	fPreviousEventTLVPosElectronTClone(NULL),	
      //  fElectronRecalculatedv1(),
      //  fElectronRecalculatedv2(),
	fElectronMass(-1),
	fGammaMass(-1),
	fPi0Mass(-1),
	fEtaMass(-1),
	fGammaWidth(-1),
	fPi0Width(-1),
	fEtaWidth(-1),
	fMinOpeningAngleGhostCut(0.),
	fEsdTrackCuts(NULL),
	fCalculateBackground(kFALSE),
	fWriteNtuple(kFALSE),
	fGammaNtuple(NULL),
	fNeutralMesonNtuple(NULL),
	fTotalNumberOfAddedNtupleEntries(0),
	fChargedParticles(NULL),
	fChargedParticlesId(),
	fGammaPtHighest(0.),
	fMinPtForGammaJet(1.),
	fMinIsoConeSize(0.2),
	fMinPtIsoCone(0.7),
	fMinPtGamChargedCorr(0.5),
	fMinPtJetCone(0.5),
	fLeadingChargedIndex(-1),
	fLowPtMapping(1.),
	fHighPtMapping(3.),
	fDoCF(kFALSE),
	fAODGamma(NULL),
	//fAODPi0(NULL),
	//fAODOmega(NULL),
	fAODBranchName("GammaConv"),
	fKFCreateAOD(kTRUE),
	fKFForceAOD(kFALSE),
	fKFDeltaAODFileName(""),
	fDoNeutralMesonV0MCCheck(kFALSE),
	fUseTrackMultiplicityForBG(kTRUE),
	fMoveParticleAccordingToVertex(kFALSE),
	fApplyChi2Cut(kFALSE),
	fNRandomEventsForBG(15),
	fNDegreesPMBackground(15),
	fDoRotation(kTRUE),
	fCheckBGProbability(kTRUE),
	fRemovePileUp(kFALSE),
	fSelectV0AND(kFALSE),
	fTriggerAnalysis(NULL),
	fMultiplicity(0),
	fUseMultiplicity(0), 
	fUseMultiplicityBin(0), 
	fUseHBTMultiplicity(0),
	fUseHBTMultiplicityBin(0),
	fUseCentrality(0), 
	fUseCentralityBin(0),
	fRandom(0)
{
	// Common I/O in slot 0, don't define when inheriting from AnalysisTaskSE
	// DefineInput (0, TChain::Class());	
	// DefineOutput(0, TTree::Class()); 
	
	// Your private output
	DefineOutput(1, TList::Class());
	DefineOutput(2, AliCFContainer::Class());	// for CF
	
	
	// Define standard ESD track cuts for Gamma-hadron correlation 
	SetESDtrackCuts();

}

AliAnalysisTaskGammaConversion::~AliAnalysisTaskGammaConversion() 
{
	// Remove all pointers
	
	if(fOutputContainer){
		fOutputContainer->Clear() ; 
		delete fOutputContainer ;
	}
	if(fHistograms){
		delete fHistograms;
	}
	if(fV0Reader){
		delete fV0Reader;
	}
	
	// for CF
	if(fCFManager){
		delete fCFManager;
	}

	if(fEsdTrackCuts){
		delete fEsdTrackCuts;
	}

	//Delete AODs
	if (fAODGamma) {
		fAODGamma->Clear();
		delete fAODGamma;
	}
	fAODGamma = NULL;

	/*if (fAODPi0) {
		fAODPi0->Clear();
		delete fAODPi0;
	}
	fAODPi0 = NULL;

	if (fAODOmega) {
		fAODOmega->Clear();
		delete fAODOmega;
	}
	fAODOmega = NULL;
       */
	if(fTriggerAnalysis) {
		delete fTriggerAnalysis;
	}


}


void AliAnalysisTaskGammaConversion::Init()
{
	// Initialization
	// AliLog::SetGlobalLogLevel(AliLog::kError);
}
void AliAnalysisTaskGammaConversion::SetESDtrackCuts()
{
	// SetESDtrackCuts
	if (fEsdTrackCuts!=NULL){
		delete fEsdTrackCuts;
	}
	fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
	//standard cuts from:
	//http://aliceinfo.cern.ch/alicvs/viewvc/PWG0/dNdEta/CreateCuts.C?revision=1.4&view=markup

	// Cuts used up to 3rd of March

	//	fEsdTrackCuts->SetMinNClustersTPC(50);
	//	 fEsdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
	//	 fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
	//	 fEsdTrackCuts->SetRequireITSRefit(kTRUE);
	//	 fEsdTrackCuts->SetMaxNsigmaToVertex(3);
	//	 fEsdTrackCuts->SetRequireSigmaToVertex(kTRUE);

	//------- To be tested-----------
	// Cuts used up to 26th of Agost
	//		Int_t minNClustersTPC = 70;
	//		Double_t maxChi2PerClusterTPC = 4.0;
	//		Double_t maxDCAtoVertexXY = 2.4; // cm
	//		Double_t maxDCAtoVertexZ	= 3.2; // cm
	//		fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	//		fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
	//		fEsdTrackCuts->SetRequireITSRefit(kTRUE);
	//		//	 fEsdTrackCuts->SetRequireTPCStandAlone(kTRUE);
	//		fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	//		fEsdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
	//		fEsdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
	//		fEsdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
	//		fEsdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
	//		fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
	//		fEsdTrackCuts->SetEtaRange(-0.8, 0.8);
	//		fEsdTrackCuts->SetPtRange(0.15);

	//		fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);


	// Using standard function	for setting Cuts
	Bool_t selectPrimaries=kTRUE;
	fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
	fEsdTrackCuts->SetMaxDCAToVertexZ(2);
	fEsdTrackCuts->SetEtaRange(-0.8, 0.8);
	fEsdTrackCuts->SetPtRange(0.15);
	
	//-----	From Jacek 10.03.03 ------------------/
	//		 minNClustersTPC = 70;
	//		 maxChi2PerClusterTPC = 4.0;
	//		 maxDCAtoVertexXY = 2.4; // cm
	//		 maxDCAtoVertexZ	= 3.2; // cm

	//		 esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	//		 esdTrackCuts->SetRequireTPCRefit(kFALSE);
	//		 esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
	//		 esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	//		 esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
	//		 esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
	//		 esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
	//		 esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
	//		 esdTrackCuts->SetDCAToVertex2D(kTRUE);
	


	//	fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	//	fV0Reader->SetESDtrackCuts(fEsdTrackCuts);
}

void AliAnalysisTaskGammaConversion::UserExec(Option_t */*option*/)
{
	// Execute analysis for current event

	//	Load the esdpid from the esdhandler if exists (tender was applied) otherwise set the Bethe Bloch parameters
	Int_t eventQuality=-1;

	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliESDInputHandler *esdHandler=0x0;
	if ( (esdHandler=dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler())) && esdHandler->GetESDpid() ){
		AliV0Reader::SetESDpid(esdHandler->GetESDpid());
	} else {
		//load esd pid bethe bloch parameters depending on the existance of the MC handler
		// yes: MC parameters
		// no:	data parameters
		if (!AliV0Reader::GetESDpid()){
			if (fMCEvent ) {
				AliV0Reader::InitESDpid();
			} else {
				AliV0Reader::InitESDpid(1);
			}
		}
	} 

	if(fAODGamma) fAODGamma->Delete();
       // if(fAODPi0) fAODPi0->Delete();
       // if(fAODOmega) fAODOmega->Delete();


	//	if(fV0Reader == NULL){ // coverty does not permit this test
	// Write warning here cuts and so on are default if this ever happens
	//	}

	if (fMCEvent ) {
		// To avoid crashes due to unzip errors. Sometimes the trees are not there.

		AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
		if (!mcHandler){ 
			AliError("Could not retrive MC event handler!"); 

			eventQuality=0;
			fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
			return; 
		}
		if (!mcHandler->InitOk() ){
			eventQuality=0;
			fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
			return;
		}
		if (!mcHandler->TreeK() ){
			eventQuality=0;
			fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
			return;
		}
		if (!mcHandler->TreeTR() ) {
			eventQuality=0;
			fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
			return;
		}
	}

	fV0Reader->SetInputAndMCEvent(InputEvent(), MCEvent());

	fV0Reader->Initialize();
	fDoMCTruth = fV0Reader->GetDoMCTruth();

	
	if(fKFReconstructedGammasTClone == NULL){
		fKFReconstructedGammasTClone = new TClonesArray("AliKFConversionPhoton",0);
	}
	if(fCurrentEventPosElectronTClone == NULL){
		fCurrentEventPosElectronTClone = new TClonesArray("AliESDtrack",0);
	}
	if(fCurrentEventNegElectronTClone == NULL){
		fCurrentEventNegElectronTClone = new TClonesArray("AliESDtrack",0);
	}
	if(fKFReconstructedGammasCutTClone == NULL){
		fKFReconstructedGammasCutTClone = new TClonesArray("AliKFConversionPhoton",0);
	}
	if(fPreviousEventTLVNegElectronTClone == NULL){
		fPreviousEventTLVNegElectronTClone = new TClonesArray("TLorentzVector",0);
	}
	if(fPreviousEventTLVPosElectronTClone == NULL){
		fPreviousEventTLVPosElectronTClone	= new TClonesArray("TLorentzVector",0);
	}
	if(fChargedParticles == NULL){
		fChargedParticles = new TClonesArray("AliESDtrack",0);
	}

	if(fKFReconstructedPi0sTClone == NULL){
		fKFReconstructedPi0sTClone = new TClonesArray("AliKFConversionMother",0);
	}
 
/*	if(fKFRecalculatedGammasTClone == NULL){
		fKFRecalculatedGammasTClone = new TClonesArray("AliKFParticle",0);
	}
  */
	if(fTriggerAnalysis== NULL){
		fTriggerAnalysis = new AliTriggerAnalysis;
	}

	//clear TClones
	fKFReconstructedGammasTClone->Delete();
	fCurrentEventPosElectronTClone->Delete();
	fCurrentEventNegElectronTClone->Delete();
	fKFReconstructedGammasCutTClone->Delete();
	fPreviousEventTLVNegElectronTClone->Delete();
	fPreviousEventTLVPosElectronTClone->Delete();
	fKFReconstructedPi0sTClone->Delete();
     //   fKFRecalculatedGammasTClone->Delete();

	//clear vectors

	fChargedParticles->Delete();	

	fChargedParticlesId.clear();	


	//Clear the data in the v0Reader
	//	fV0Reader->UpdateEventByEventData();

	//Take Only events with proper trigger
	/*
		if(fTriggerCINT1B){
		if(!fV0Reader->GetESDEvent()->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")) return;
		}
	*/
	Bool_t v0A			 = fTriggerAnalysis->IsOfflineTriggerFired(fV0Reader->GetESDEvent(), AliTriggerAnalysis::kV0A);
	Bool_t v0C			 = fTriggerAnalysis->IsOfflineTriggerFired(fV0Reader->GetESDEvent(), AliTriggerAnalysis::kV0C);
	Bool_t v0AND = v0A && v0C;

	if(fSelectV0AND && !v0AND){
		eventQuality=5;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		if(fDoMCTruth){
			if(!fV0Reader->GetIsHeavyIon()){
	CheckMesonProcessTypeEventQuality(eventQuality);
			}
		}

		return;
	}

	if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
		//		cout<< "Event not taken"<< endl;
		eventQuality=1;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		if(fDoMCTruth){
			if(!fV0Reader->GetIsHeavyIon()){
				CheckMesonProcessTypeEventQuality(eventQuality);
			}
		}
		return; // aborts if the primary vertex does not have contributors.
	}

	 

	if(!fV0Reader->CheckForPrimaryVertexZ() ){
		eventQuality=2;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		if(fDoMCTruth){
			if(!fV0Reader->GetIsHeavyIon()){
				CheckMesonProcessTypeEventQuality(eventQuality);
			}
		}
		return;
	}


	if(fV0Reader->GetESDEvent()->GetPrimaryVertexTracks()->GetNContributors()>0) {
		fHistograms->FillHistogram("ESD_GlobalPrimaryVtxZ",fV0Reader->GetESDEvent()->GetPrimaryVertex()->GetZ());
	}else{
		if(fV0Reader->GetESDEvent()->GetPrimaryVertexSPD()->GetNContributors()>0) {
			fHistograms->FillHistogram("ESD_SPDPrimaryVtxZ",fV0Reader->GetESDEvent()->GetPrimaryVertex()->GetZ());
		}
	}

	if(fRemovePileUp && fV0Reader->GetESDEvent()->IsPileupFromSPD()) {
		eventQuality=4;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	

	Int_t tTracklet=0, tITSTPC=0, tITSPure=0;
	fV0Reader->GetESDEvent()->EstimateMultiplicity(tTracklet, tITSTPC, tITSPure, 0.8);
	Int_t fMultiplicityITS = tITSPure;

	fMultiplicity =	fEsdTrackCuts->CountAcceptedTracks(fV0Reader->GetESDEvent());
	Int_t fMultiplicityStandard = fMultiplicity;

	if( fUseHBTMultiplicity==1) {
		fMultiplicity = fMultiplicityITS;

	}


	fHistograms->FillHistogram("ESD_MultiplicityDeviation",fMultiplicityStandard,fMultiplicityITS);

	

	if(fUseMultiplicity!=0 && CalculateMultiplicityBin()!=fUseMultiplicityBin ){
		eventQuality=6;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}

	
	if(fUseHBTMultiplicity!=0 && CalculateMultiplicityBin()!=fUseHBTMultiplicityBin ){
		eventQuality=6;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}


	if(fV0Reader->GetIsHeavyIon()){
		if(fUseCentrality>0){
			AliCentrality *esdCentrality = fV0Reader->GetESDEvent()->GetCentrality();
			Int_t centralityC = -1;

			if(fUseCentrality==1){
	centralityC = esdCentrality->GetCentralityClass10("V0M");
	if( centralityC != fUseCentralityBin ){
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
			}

			if(fUseCentrality==2){
	centralityC = esdCentrality->GetCentralityClass10("CL1");
	if( centralityC != fUseCentralityBin ){
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
			}

			////////////////////////////////////// RRnew start /////////////////////////////////////////////////////
			if(fUseCentrality==3){
	centralityC = esdCentrality->GetCentralityClass10("V0M");
	if( (fUseCentralityBin == 0) && (centralityC!=0) ){ // 0-10%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;		
	}
	if( (fUseCentralityBin == 1) && (centralityC!=1) ){ // 10-20%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 2) && (centralityC!=2) && (centralityC!=3) ){ // 20-40%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 3) && (centralityC!=0) && (centralityC!=1) ){ // 0-20%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 4) && (centralityC!=4) && (centralityC!=5) ){ // 40-60%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 6) && (centralityC!=6) && (centralityC!=7) && (centralityC!=8) ){ // 60-90%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 7) && (centralityC!=6) && (centralityC!=7) ){ // 60-80%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 8) && (centralityC>=8) ){ // 0-80%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 9) && (centralityC>=9) ){ // 0-90%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
			}

			if(fUseCentrality==4){
	centralityC = esdCentrality->GetCentralityClass10("CL1");
	if( (fUseCentralityBin == 0) && (centralityC!=0) ){ // 0-10%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;		
	}
	if( (fUseCentralityBin == 1) && (centralityC!=1) ){ // 10-20%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 2) && (centralityC!=2) && (centralityC!=3) ){ // 20-40%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 4) && (centralityC!=4) && (centralityC!=5) ){ // 40-60%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
	if( (fUseCentralityBin == 6) && (centralityC!=6) && (centralityC!=7) && (centralityC!=8) ){ // 60-90%
		eventQuality=7;
		fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}
			}
			////////////////////////////////////// RRnew end ///////////////////////////////////////////////////////

		}
	}


	eventQuality=3;

	fHistograms->FillHistogram("ESD_EventQuality",eventQuality);



	fHistograms->FillHistogram("ESD_NumberOfGoodESDTracks",fMultiplicity);
	if (fV0Reader->GetNumberOfContributorsVtx()>=1){
		fHistograms->FillHistogram("ESD_NumberOfGoodESDTracksVtx",fMultiplicity);
	} 


	// Process the MC information
	if(fDoMCTruth){
		ProcessMCData();
	}


	
	//Process the v0 information with no cuts
	ProcessV0sNoCut();

	// Process the v0 information
	ProcessV0s();
	

	//Fill Gamma AOD
	if(fKFCreateAOD) {
		FillAODWithConversionGammas() ; 
	}


	// Process reconstructed gammas
	if(fDoNeutralMeson == kTRUE){
		ProcessGammasForNeutralMesonAnalysis();

	}
	
	if(fDoMCTruth == kTRUE){
		CheckV0Efficiency();
	}
	//Process reconstructed gammas electrons for Chi_c Analysis
	if(fDoChic == kTRUE){
		ProcessGammaElectronsForChicAnalysis();
	}
	// Process reconstructed gammas for gamma Jet/hadron correlations
	if(fDoJet == kTRUE){
		ProcessGammasForGammaJetAnalysis();
	}
	
	//calculate background if flag is set
	if(fCalculateBackground){
		CalculateBackground();
	}

	if(fDoNeutralMeson == kTRUE){
		//		 ProcessConvPHOSGammasForNeutralMesonAnalysis();
		if(fDoOmegaMeson == kTRUE){
			ProcessGammasForOmegaMesonAnalysis();
		}
	}


	//Must set fForceAOD to true for the AOD to get filled. (Unless called by other task)
	if(fKFForceAOD) {
		if (!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()) { 
			AliFatal("Cannot run ESD filter without an output event handler");
	 
		} else {
			if(fAODGamma && fAODGamma->GetEntriesFast() > 0) {
	AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
	AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);	 
			}
		}
	
	}

	///Make sure delta aod is filled if standard aod is filled (for synchronization when reading aod with standard aod)
	if(fKFCreateAOD) {
		AliAODHandler * aodhandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
		if (aodhandler && aodhandler->GetFillAOD()) {
			AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);	 
		}
	}


	//Clear the data in the v0Reader
	fV0Reader->UpdateEventByEventData();
	//if(fRecalculateV0ForGamma==kTRUE){
       // 	RecalculateV0ForGamma();
       // }
	PostData(1, fOutputContainer);
	PostData(2, fCFManager->GetParticleContainer());	// for CF
	
}

// void AliAnalysisTaskGammaConversion::ConnectInputData(Option_t *option){
//	 // see header file for documentation
//	 //	printf("	 ConnectInputData %s\n", GetName());

//	 AliAnalysisTaskSE::ConnectInputData(option);

//	 if(fV0Reader == NULL){
//		 // Write warning here cuts and so on are default if this ever happens
//	 }
//	 fV0Reader->Initialize();
//	 fDoMCTruth = fV0Reader->GetDoMCTruth();
// }

void AliAnalysisTaskGammaConversion::CheckMesonProcessTypeEventQuality(Int_t evtQ){
	// Check meson process type event quality
	fStack= MCEvent()->Stack();
	fGCMCEvent=MCEvent();

 

	for (Int_t iTracks = 0; iTracks < fStack->GetNprimary(); iTracks++) {
		TParticle* particle = (TParticle *)fStack->Particle(iTracks);
		if (!particle) {
			//print warning here
			continue;
		}
		//		 if(particle->GetPdgCode()!=111 || particle->GetPdgCode()!=221){
		//			 continue;
		//		 }
		
		Double_t rapidity;
		if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
			rapidity=8.;
		}
		else{
			rapidity = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())));
		}	
		

		if( particle->GetPdgCode()==111){
			fHistograms->FillHistogram("MC_Test_AllPi0_Pt", particle->Pt());
			if(particle->GetNDaughters()==2){
				fHistograms->FillHistogram("MC_Test_2DaughPi0_Pt", particle->Pt());
				if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
					fHistograms->FillHistogram("MC_Test_2DaughPi0_Rap_Pt", particle->Pt());
				}
			}
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_Test_Pi0_Rap_Pt", particle->Pt());
			}
		}
			
			
		if(particle->GetPdgCode()==221){
			fHistograms->FillHistogram("MC_Test_AllEta_Pt", particle->Pt());
			if(particle->GetNDaughters()==2){
				fHistograms->FillHistogram("MC_Test_2DaughEta_Pt", particle->Pt());
				if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
					fHistograms->FillHistogram("MC_Test_2DaughEta_Rap_Pt", particle->Pt());
				}
			}
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_Test_Eta_Rap_Pt", particle->Pt());
			}
			
		}


		if(particle->GetPdgCode()!=111){		 //Pi0
			continue;
		}


		if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut() ) continue; 
		if(fV0Reader->GetIsHeavyIon()) continue;

		if(evtQ==1){
			switch(GetProcessType(fGCMCEvent)){
			case	kProcSD:
				fHistograms->FillHistogram("MC_SD_EvtQ1_Pi0_Pt", particle->Pt());
				break;
			case	kProcDD:
				fHistograms->FillHistogram("MC_DD_EvtQ1_Pi0_Pt", particle->Pt());
				break;
			case	kProcND:
				fHistograms->FillHistogram("MC_ND_EvtQ1_Pi0_Pt", particle->Pt());
				break;
			default:
				AliError("Unknown Process");
			}
		}
		if(evtQ==2){
			switch(GetProcessType(fGCMCEvent)){
			case	kProcSD:
				fHistograms->FillHistogram("MC_SD_EvtQ2_Pi0_Pt", particle->Pt());
			break;
			case	kProcDD:
				fHistograms->FillHistogram("MC_DD_EvtQ2_Pi0_Pt", particle->Pt());
			break;
			case	kProcND:
				fHistograms->FillHistogram("MC_ND_EvtQ2_Pi0_Pt", particle->Pt());
			break;
			default:
				AliError("Unknown Process");
			}
		}

		if(evtQ==4){
			switch(GetProcessType(fGCMCEvent)){
			case	kProcSD:
				fHistograms->FillHistogram("MC_SD_EvtQ4_Pi0_Pt", particle->Pt());
			break;
			case	kProcDD:
				fHistograms->FillHistogram("MC_DD_EvtQ4_Pi0_Pt", particle->Pt());
			break;
			case	kProcND:
				fHistograms->FillHistogram("MC_ND_EvtQ4_Pi0_Pt", particle->Pt());
			break;
			default:
				AliError("Unknown Process");
			}
		}

		if(evtQ==5){
			switch(GetProcessType(fGCMCEvent)){
			case	kProcSD:
				fHistograms->FillHistogram("MC_SD_EvtQ5_Pi0_Pt", particle->Pt());
			break;
			case	kProcDD:
				fHistograms->FillHistogram("MC_DD_EvtQ5_Pi0_Pt", particle->Pt());
			break;
			case	kProcND:
				fHistograms->FillHistogram("MC_ND_EvtQ5_Pi0_Pt", particle->Pt());
			break;
			default:
				AliError("Unknown Process");
			}
		}

	}

}

void AliAnalysisTaskGammaConversion::ProcessMCData(){
	// see header file for documentation
	//InputEvent(), MCEvent());
	/* TestAnaMarin
		fStack = fV0Reader->GetMCStack();
		fMCTruth = fV0Reader->GetMCTruth();	// for CF
		fGCMCEvent = fV0Reader->GetMCEvent();	// for CF
	*/
	fStack= MCEvent()->Stack();
	fGCMCEvent=MCEvent();
		
	// for CF
	Double_t containerInput[3];
	if(fDoCF){
		if(!fGCMCEvent) cout << "NO MC INFO FOUND" << endl;
		fCFManager->SetEventInfo(fGCMCEvent);
	} 
	// end for CF
		
	if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
		return; // aborts if the primary vertex does not have contributors.
	}
	
	Int_t nCharged = 0;
	Int_t nCharged150MeV = 0;

	for (Int_t iTracks = 0; iTracks < fStack->GetNprimary(); iTracks++) {
		//	for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
		TParticle* particle = (TParticle *)fStack->Particle(iTracks);



		if (!particle) {
			//print warning here
			continue;
		}
				


		
		///////////////////////Begin Chic Analysis/////////////////////////////
		if(fDoChic) {
			if(particle->GetPdgCode() == 443){//Is JPsi	
				if(particle->GetNDaughters()==2){
					if(TMath::Abs(fStack->Particle(particle->GetFirstDaughter())->GetPdgCode()) == 11 &&
						TMath::Abs(fStack->Particle(particle->GetLastDaughter())->GetPdgCode()) == 11){

						TParticle* daug0 = fStack->Particle(particle->GetFirstDaughter());
						TParticle* daug1 = fStack->Particle(particle->GetLastDaughter());
						if(TMath::Abs(daug0->Eta()) < 0.9 && TMath::Abs(daug1->Eta()) < 0.9)
							fHistograms->FillTable("Table_Electrons",3);//e+ e-	from J/Psi inside acceptance
								
						if( TMath::Abs(daug0->Eta()) < 0.9){
							if(daug0->GetPdgCode() == -11)
								fHistograms->FillTable("Table_Electrons",1);//e+	from J/Psi inside acceptance
							else
								fHistograms->FillTable("Table_Electrons",2);//e-	 from J/Psi inside acceptance
									
						}
						if(TMath::Abs(daug1->Eta()) < 0.9){
							if(daug1->GetPdgCode() == -11)
								fHistograms->FillTable("Table_Electrons",1);//e+	from J/Psi inside acceptance
							else
								fHistograms->FillTable("Table_Electrons",2);//e-	 from J/Psi inside acceptance
						}
					}
				}
			}
			//							const int CHI_C0	 = 10441;
			//							const int CHI_C1	 = 20443;
			//							const int CHI_C2	 = 445
			if(particle->GetPdgCode() == 22){//gamma from JPsi
				if(particle->GetMother(0) > -1){
					if(fStack->Particle(particle->GetMother(0))->GetPdgCode() == 10441 ||
					fStack->Particle(particle->GetMother(0))->GetPdgCode() == 20443 ||
					fStack->Particle(particle->GetMother(0))->GetPdgCode() == 445){
						if(TMath::Abs(particle->Eta()) < 1.2)
							fHistograms->FillTable("Table_Electrons",17);// gamma from chic inside accptance
					}
				}
			}
			if(particle->GetPdgCode() == 10441 || particle->GetPdgCode() == 20443 || particle->GetPdgCode() == 445){
				if( particle->GetNDaughters() == 2){
					TParticle* daug0 = fStack->Particle(particle->GetFirstDaughter());
					TParticle* daug1 = fStack->Particle(particle->GetLastDaughter());
								
					if( (daug0->GetPdgCode() == 443 || daug0->GetPdgCode() == 22) && (daug1->GetPdgCode() == 443 || daug1->GetPdgCode() == 22) ){
						if( daug0->GetPdgCode() == 443){
							TParticle* daugE0 = fStack->Particle(daug0->GetFirstDaughter());
							TParticle* daugE1 = fStack->Particle(daug0->GetLastDaughter());
							if( TMath::Abs(daug1->Eta()) < 1.2 && TMath::Abs(daugE0->Eta()) < 0.9 && TMath::Abs(daugE1->Eta()) < 0.9 )
								fHistograms->FillTable("Table_Electrons",18);
								
						} else if (daug1->GetPdgCode() == 443){
							TParticle* daugE0 = fStack->Particle(daug1->GetFirstDaughter());
							TParticle* daugE1 = fStack->Particle(daug1->GetLastDaughter());
							if( TMath::Abs(daug0->Eta()) < 1.2 && TMath::Abs(daugE0->Eta()) < 0.9 && TMath::Abs(daugE1->Eta()) < 0.9 )
							fHistograms->FillTable("Table_Electrons",18);
						}//else if
					}//gamma o Jpsi
				}//GetNDaughters
			}
		}
				
		/////////////////////End Chic Analysis////////////////////////////
				
		//		if(TMath::Abs(particle->Eta())> fV0Reader->GetEtaCut() )	continue;




		if(particle->Eta() <0.8 && particle->Eta() > (-0.8) && particle->GetPDG()->Charge() != 0) {
		   nCharged++;
		   if(particle->Pt()>0.150){
		   nCharged150MeV++;}
		}
		
		
				
		if(particle->R()>fV0Reader->GetMaxRCut())	continue; // cuts on distance from collision point
				
		Double_t tmpPhi=particle->Phi();
				
		if(particle->Phi()> TMath::Pi()){
			tmpPhi = particle->Phi()-(2*TMath::Pi());
		}
				
		Double_t rapidity;
		if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
			rapidity=8.;
		} else{
			rapidity = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())));
		}	
				
		if( particle->GetPdgCode()==111){

			//AM move here, otherwise for evt we consider only pi0 to 2 g
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				if(!fV0Reader->GetIsHeavyIon()) {
				
					switch(GetProcessType(fGCMCEvent)){
					case	kProcSD:
						fHistograms->FillHistogram("MC_SD_EvtQ3_Pi0_Pt", particle->Pt());
						break;
					case	kProcDD:
						fHistograms->FillHistogram("MC_DD_EvtQ3_Pi0_Pt", particle->Pt());
						break;
					case	kProcND:
						fHistograms->FillHistogram("MC_ND_EvtQ3_Pi0_Pt", particle->Pt());
						break;
					default:
						AliError("Unknown Process");
					}	
				}	
			}
			fHistograms->FillHistogram("MC_Test_AllPi0_Pt", particle->Pt());
			if(particle->GetNDaughters()==2){
				fHistograms->FillHistogram("MC_Test_2DaughPi0_Pt", particle->Pt());
				if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
					fHistograms->FillHistogram("MC_Test_2DaughPi0_Rap_Pt", particle->Pt());
				}
			}
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_Test_Pi0_Rap_Pt", particle->Pt());
			}
		}
			
			
		if(particle->GetPdgCode()==221){
			fHistograms->FillHistogram("MC_Test_AllEta_Pt", particle->Pt());
			if(particle->GetNDaughters()==2){
				fHistograms->FillHistogram("MC_Test_2DaughEta_Pt", particle->Pt());
				if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
					fHistograms->FillHistogram("MC_Test_2DaughEta_Rap_Pt", particle->Pt());
				}
			}
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_Test_Eta_Rap_Pt", particle->Pt());
			}
		}

		if(  particle->GetPdgCode()==310 ){
			if(TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_K0S_Pt", particle->Pt());
			}
		}

		if(iTracks<=fStack->GetNprimary() ){
			if ( particle->GetPdgCode()== -211 ||	particle->GetPdgCode()== 211 ||
				particle->GetPdgCode()== 2212 ||	particle->GetPdgCode()==-2212 ||
				particle->GetPdgCode()== 321	||	particle->GetPdgCode()==-321 ){
				if(TMath::Abs(particle->Eta())> 0.8 )	continue;	 // Eta cut used in charged particle spectrum
				//if( !particle->IsPhysicalPrimary() ){
				// cout<<"not Physical primary"<< particle->IsPhysicalPrimary()<<endl;
				//}
				if(particle->GetMother(0)>-1){
					// cout<<"Mother ::"<<fStack->Particle(particle->GetMother(0))->GetPdgCode()<<endl;
					if (fStack->Particle(particle->GetMother(0))->GetPdgCode()== -211 ||fStack->Particle(particle->GetMother(0))->GetPdgCode()== -3122 ){
						//						cout<<"Mother K0, lambda::"<<fStack->Particle(particle->GetMother(0))->GetPdgCode()<<endl;
					continue;
					}
				}
				
				fHistograms->FillHistogram("MC_PhysicalPrimaryCharged_Pt", particle->Pt());

				if (particle->GetPdgCode() == 211 ) fHistograms->FillHistogram("MC_PiPlus_Pt", particle->Pt());
				if (particle->GetPdgCode() == 321 ) fHistograms->FillHistogram("MC_KaonPlus_Pt", particle->Pt());
				if (particle->GetPdgCode() == 2212 ) fHistograms->FillHistogram("MC_Proton_Pt", particle->Pt());
				if (particle->GetPdgCode() == -211 ) fHistograms->FillHistogram("MC_PiMinus_Pt", particle->Pt());
				if (particle->GetPdgCode() == -321 ) fHistograms->FillHistogram("MC_KaonMinus_Pt", particle->Pt());
				if (particle->GetPdgCode() == -2212 ) fHistograms->FillHistogram("MC_AntiProton_Pt", particle->Pt()); 
			}
			if(TMath::Abs(particle->Eta())<=0.8 ){
				if (particle->GetPdgCode() == 111 ) fHistograms->FillHistogram("MC_Pi0_Test_Pt", particle->Pt());
			}
		}


		//process the gammas
		if (particle->GetPdgCode() == 22){
			if(TMath::Abs(particle->Eta())> fV0Reader->GetEtaCut() || TMath::Abs(particle->Eta())< fV0Reader->GetEtaCutMin())	continue;			

			if(particle->GetMother(0) >-1 && fStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
				continue; // no photon as mothers!
			}
			
			if(particle->GetMother(0) >= fStack->GetNprimary()){
				continue; // the gamma has a mother, and it is not a primary particle
			}
				
			if(particle->GetMother(0) >-1){
				fHistograms->FillHistogram("MC_DecayAllGamma_Pt", particle->Pt()); // All
				switch(fStack->Particle(particle->GetMother(0))->GetPdgCode()){
					case 111: // Pi0
						fHistograms->FillHistogram("MC_DecayPi0Gamma_Pt", particle->Pt());
						break;
					case 113: // Rho0
						fHistograms->FillHistogram("MC_DecayRho0Gamma_Pt", particle->Pt());
						break;
					case 221: // Eta
						fHistograms->FillHistogram("MC_DecayEtaGamma_Pt", particle->Pt());
						break;
					case 223: // Omega
						fHistograms->FillHistogram("MC_DecayOmegaGamma_Pt", particle->Pt());
						break;
					case 310: // K_s0
						fHistograms->FillHistogram("MC_DecayK0sGamma_Pt", particle->Pt());
						break;
					case 331: // Eta'
						fHistograms->FillHistogram("MC_DecayEtapGamma_Pt", particle->Pt());
						break;
					case 333: // Phi
						fHistograms->FillHistogram("MC_DecayPhiGamma_Pt", particle->Pt());
						break;
				}
			}
			
			fHistograms->FillHistogram("MC_allGamma_Energy", particle->Energy());
			fHistograms->FillHistogram("MC_allGamma_Pt", particle->Pt());
			fHistograms->FillHistogram("MC_allGamma_Eta", particle->Eta());
			fHistograms->FillHistogram("MC_allGamma_Phi", tmpPhi);
			fHistograms->FillHistogram("MC_allGamma_Rapid", rapidity);
					
			// for CF
			if(fDoCF){
				containerInput[0] = particle->Pt();
				containerInput[1] = particle->Eta();
				if(particle->GetMother(0) >=0){
					containerInput[2] = fStack->Particle(particle->GetMother(0))->GetMass();
				} else{
					containerInput[2]=-1;
				}
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);					// generated gamma
			}
			
			if(particle->GetMother(0) < 0 || //Phojet p+p -> Direct Photons have no mother
			((particle->GetMother(0) > -1) && 
			((TMath::Abs(fStack->Particle(particle->GetMother(0))->GetPdgCode()) < 10)|| 
			(TMath::Abs(fStack->Particle(particle->GetMother(0))->GetPdgCode()) ==21) )) //Pythia p+p -> Direct Photons have quarksor gluons as mother
			){	 // direct gamma
				fHistograms->FillHistogram("MC_allDirectGamma_Energy",particle->Energy());
				fHistograms->FillHistogram("MC_allDirectGamma_Pt", particle->Pt());
				fHistograms->FillHistogram("MC_allDirectGamma_Eta", particle->Eta());
				fHistograms->FillHistogram("MC_allDirectGamma_Phi", tmpPhi);
				fHistograms->FillHistogram("MC_allDirectGamma_Rapid", rapidity);				
			}
					
			// looking for conversion (electron + positron from pairbuilding (= 5) )
			TParticle* ePos = NULL;
			TParticle* eNeg = NULL;
					
			if(particle->GetNDaughters() >= 2){
				for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
					TParticle *tmpDaughter = fStack->Particle(daughterIndex);
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
				continue;
			}
					
					
			Double_t ePosPhi = ePos->Phi();
			if(ePos->Phi()> TMath::Pi()) ePosPhi = ePos->Phi()-(2*TMath::Pi());
					
			Double_t eNegPhi = eNeg->Phi();
			if(eNeg->Phi()> TMath::Pi()) eNegPhi = eNeg->Phi()-(2*TMath::Pi());
							
			if(ePos->Pt()<fV0Reader->GetSinglePtCut() || eNeg->Pt()<fV0Reader->GetSinglePtCut()){
				continue; // no reconstruction below the Pt cut
			}
					
			if( TMath::Abs(ePos->Eta())> fV0Reader->GetEtaCut() || TMath::Abs(ePos->Eta())< fV0Reader->GetEtaCutMin()  || 
			    TMath::Abs(eNeg->Eta())> fV0Reader->GetEtaCut() || TMath::Abs(eNeg->Eta())< fV0Reader->GetEtaCutMin() ) {
				continue;
			}	
					
			if(ePos->R()>fV0Reader->GetMaxRCut()){
				continue; // cuts on distance from collision point
			}

			if(TMath::Abs(ePos->Vz()) > fV0Reader->GetMaxZCut()){
				continue;	 // outside material
			}
			if(TMath::Abs(eNeg->Vz()) > fV0Reader->GetMaxZCut()){
				continue;	 // outside material
			}

			if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue())){
				continue;							 // line cut to exclude regions where we do not reconstruct
			} else if ( fV0Reader->GetEtaCutMin() != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())	){
				continue;
			}
					
			if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue())){
				continue;							 // line cut to exclude regions where we do not reconstruct
			} else if ( fV0Reader->GetEtaCutMin() != -0.1 &&   eNeg->R() >= ((TMath::Abs(eNeg->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())	){
				continue;
			}
					
			// for CF
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructable);	// reconstructable gamma	
			}
			fHistograms->FillHistogram("MC_ConvGamma_Energy", particle->Energy());
			fHistograms->FillHistogram("MC_ConvGamma_Pt", particle->Pt());
			// Move down, in the if that mother exists
			//			if(fStack->Particle(particle->GetMother(0))->GetPdgCode() == 221)
			//			fHistograms->FillHistogram("MC_ConvEtaGamma_Pt", particle->Pt());
			fHistograms->FillHistogram("MC_ConvGamma_Eta", particle->Eta());
			fHistograms->FillHistogram("MC_ConvGamma_Phi", tmpPhi);
			fHistograms->FillHistogram("MC_ConvGamma_Rapid", rapidity);
			fHistograms->FillHistogram("MC_ConvGamma_Pt_Eta", particle->Pt(),particle->Eta());
					
			fHistograms->FillHistogram("MC_E_Energy", eNeg->Energy());
			fHistograms->FillHistogram("MC_E_Pt", eNeg->Pt());
			fHistograms->FillHistogram("MC_E_Eta", eNeg->Eta());
			fHistograms->FillHistogram("MC_E_Phi", eNegPhi);
					
			fHistograms->FillHistogram("MC_P_Energy", ePos->Energy());
			fHistograms->FillHistogram("MC_P_Pt", ePos->Pt());
			fHistograms->FillHistogram("MC_P_Eta", ePos->Eta());
			fHistograms->FillHistogram("MC_P_Phi", ePosPhi);
					
			if(!fV0Reader->GetIsHeavyIon()){
					
				// begin Mapping 
				Int_t rBin		= fHistograms->GetRBin(ePos->R());
				Int_t zBin		= fHistograms->GetZBin(ePos->Vz());
				Int_t phiBin	= fHistograms->GetPhiBin(particle->Phi());
				Double_t rFMD=30;
				Double_t rITSTPCMin=50;
				Double_t rITSTPCMax=80;
				
				TVector3 vtxPos(ePos->Vx(),ePos->Vy(),ePos->Vz());	
				
				TString nameMCMappingPhiR="";
				nameMCMappingPhiR.Form("MC_Conversion_Mapping_Phi%02d_R%02d",phiBin,rBin);
				// fHistograms->FillHistogram(nameMCMappingPhiR, ePos->Vz(), particle->Eta());
				
				TString nameMCMappingPhi="";
				nameMCMappingPhi.Form("MC_Conversion_Mapping_Phi%02d",phiBin);
				//			fHistograms->FillHistogram(nameMCMappingPhi, particle->Eta());
				//fHistograms->FillHistogram(nameMCMappingPhi, ePos->Vz(), particle->Eta());
				
				TString nameMCMappingR="";
				nameMCMappingR.Form("MC_Conversion_Mapping_R%02d",rBin);
				//			fHistograms->FillHistogram(nameMCMappingR, particle->Eta());
				//fHistograms->FillHistogram(nameMCMappingR,ePos->Vz(), particle->Eta());
				
				TString nameMCMappingPhiInR="";
				nameMCMappingPhiInR.Form("MC_Conversion_Mapping_Phi_in_R_%02d",rBin);
				//			fHistograms->FillHistogram(nameMCMappingPhiInR, tmpPhi);
				fHistograms->FillHistogram(nameMCMappingPhiInR, vtxPos.Phi());
				
				TString nameMCMappingZInR="";
				nameMCMappingZInR.Form("MC_Conversion_Mapping_Z_in_R_%02d",rBin);
				fHistograms->FillHistogram(nameMCMappingZInR,ePos->Vz() );
				
				
				TString nameMCMappingPhiInZ="";
				nameMCMappingPhiInZ.Form("MC_Conversion_Mapping_Phi_in_Z_%02d",zBin);
				//			fHistograms->FillHistogram(nameMCMappingPhiInR, tmpPhi);
				fHistograms->FillHistogram(nameMCMappingPhiInZ, vtxPos.Phi());
				
				
				if(ePos->R()<rFMD){
					TString nameMCMappingFMDPhiInZ="";
					nameMCMappingFMDPhiInZ.Form("MC_Conversion_Mapping_FMD_Phi_in_Z_%02d",zBin);
					fHistograms->FillHistogram(nameMCMappingFMDPhiInZ, vtxPos.Phi());
				}
			
				if(ePos->R()>rITSTPCMin	&& ePos->R()<rITSTPCMax){
					TString nameMCMappingITSTPCPhiInZ="";
					nameMCMappingITSTPCPhiInZ.Form("MC_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",zBin);
					fHistograms->FillHistogram(nameMCMappingITSTPCPhiInZ, vtxPos.Phi());
				}
			
				TString nameMCMappingRInZ="";
				nameMCMappingRInZ.Form("MC_Conversion_Mapping_R_in_Z_%02d",zBin);
				fHistograms->FillHistogram(nameMCMappingRInZ,ePos->R() );
				
				if(particle->Pt() > fLowPtMapping && particle->Pt()< fHighPtMapping){
					TString nameMCMappingMidPtPhiInR="";
					nameMCMappingMidPtPhiInR.Form("MC_Conversion_Mapping_MidPt_Phi_in_R_%02d",rBin);
					fHistograms->FillHistogram(nameMCMappingMidPtPhiInR, vtxPos.Phi());
					
					TString nameMCMappingMidPtZInR="";
					nameMCMappingMidPtZInR.Form("MC_Conversion_Mapping_MidPt_Z_in_R_%02d",rBin);
					fHistograms->FillHistogram(nameMCMappingMidPtZInR,ePos->Vz() );
					
					
					TString nameMCMappingMidPtPhiInZ="";
					nameMCMappingMidPtPhiInZ.Form("MC_Conversion_Mapping_MidPt_Phi_in_Z_%02d",zBin);
					fHistograms->FillHistogram(nameMCMappingMidPtPhiInZ, vtxPos.Phi());
					
					
					if(ePos->R()<rFMD){
						TString nameMCMappingMidPtFMDPhiInZ="";
						nameMCMappingMidPtFMDPhiInZ.Form("MC_Conversion_Mapping_MidPt_FMD_Phi_in_Z_%02d",zBin);
						fHistograms->FillHistogram(nameMCMappingMidPtFMDPhiInZ, vtxPos.Phi());
					}
				
					TString nameMCMappingMidPtRInZ="";
					nameMCMappingMidPtRInZ.Form("MC_Conversion_Mapping_MidPt_R_in_Z_%02d",zBin);
					fHistograms->FillHistogram(nameMCMappingMidPtRInZ,ePos->R() );
					
				}
			}

			//end mapping
					
			fHistograms->FillHistogram("MC_Conversion_R",ePos->R());
			fHistograms->FillHistogram("MC_Conversion_ZR",ePos->Vz(),ePos->R());
			fHistograms->FillHistogram("MC_Conversion_XY",ePos->Vx(),ePos->Vy());
			fHistograms->FillHistogram("MC_Conversion_OpeningAngle",GetMCOpeningAngle(ePos, eNeg));
			fHistograms->FillHistogram("MC_ConvGamma_E_AsymmetryP",particle->P(),eNeg->P()/particle->P());
			fHistograms->FillHistogram("MC_ConvGamma_P_AsymmetryP",particle->P(),ePos->P()/particle->P());


			if(particle->GetMother(0) < 0 || //Phojet p+p -> Direct Photons have no mother
				 ((particle->GetMother(0) > -1) && 
					((TMath::Abs(fStack->Particle(particle->GetMother(0))->GetPdgCode()) < 10)|| 
					 (TMath::Abs(fStack->Particle(particle->GetMother(0))->GetPdgCode()) ==21) )) //Pythia p+p -> Direct Photons have quarksor gluons as mother
				 ){	 // direct gamma, still inside converted
				fHistograms->FillHistogram("MC_ConvDirectGamma_Energy",particle->Energy());
				fHistograms->FillHistogram("MC_ConvDirectGamma_Pt", particle->Pt());
				fHistograms->FillHistogram("MC_ConvDirectGamma_Eta", particle->Eta());
				fHistograms->FillHistogram("MC_ConvDirectGamma_Phi", tmpPhi);
				fHistograms->FillHistogram("MC_ConvDirectGamma_Rapid", rapidity);
						
			} else{	 // mother exits 
				if(fStack->Particle(particle->GetMother(0))->GetPdgCode() == 221)
					fHistograms->FillHistogram("MC_ConvEtaGamma_Pt", particle->Pt());
			/*	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==10441 ||//chic0 
				fStack->Particle(particle->GetMother(0))->GetPdgCode()==20443 ||//psi2S
				fStack->Particle(particle->GetMother(0))->GetPdgCode()==445	//chic2
				){ 
				fMCGammaChic.push_back(particle);
				}
			*/
			}	// end if mother exits
		} // end if particle is a photon
				
				
				
		// process motherparticles (2 gammas as daughters)
		// the motherparticle had already to pass the R and the eta cut, but no line cut.
		// the line cut is just valid for the conversions!

		// RR primary Pi0 debug ////////////////////////////////////////////////////////////////////////////////////////////
		if (particle->GetPdgCode()==111){
			if( TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
				fHistograms->FillHistogram("MC_Pi0_Pt_vs_Rapid_allDaughters", particle->Pt(),rapidity);
			}
		}
		// end RR primary Pi0 debug ////////////////////////////////////////////////////////////////////////////////////////
				
		if(particle->GetNDaughters() == 2){
					
			TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
			TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
					
			if(daughter0->GetPdgCode() != 22 || daughter1->GetPdgCode() != 22) continue; //check for gamma gamma daughters

			if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut() ) continue; 

			// Check the acceptance for both gammas
			Bool_t gammaEtaCut = kTRUE;
			if(TMath::Abs(daughter0->Eta()) > fV0Reader->GetEtaCut() || TMath::Abs(daughter0->Eta()) < fV0Reader->GetEtaCutMin() || 
			   TMath::Abs(daughter1->Eta()) > fV0Reader->GetEtaCut() || TMath::Abs(daughter1->Eta()) < fV0Reader->GetEtaCutMin()	) gammaEtaCut = kFALSE;
			
			Bool_t gammaRCut = kTRUE;
			if(daughter0->R() > fV0Reader->GetMaxRCut() || daughter1->R() > fV0Reader->GetMaxRCut()	) gammaRCut = kFALSE;
			
			// check for conversions now -> have to pass eta, R and line cut!
			Bool_t daughter0Electron = kFALSE;
			Bool_t daughter0Positron = kFALSE;
			Bool_t daughter1Electron = kFALSE;
			Bool_t daughter1Positron = kFALSE;
					
			if(daughter0->GetNDaughters() >= 2){	// first gamma
				for(Int_t TrackIndex=daughter0->GetFirstDaughter();TrackIndex<=daughter0->GetLastDaughter();TrackIndex++){
					TParticle *tmpDaughter = fStack->Particle(TrackIndex);
					if(tmpDaughter->GetUniqueID() == 5){
						if(tmpDaughter->GetPdgCode() == 11){	
							if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
								if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
									if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1	){
										if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
											daughter0Electron = kTRUE;
										}
									}
								}			
							}
						} else if(tmpDaughter->GetPdgCode() == -11){
							if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
								if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
									if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1	){
										if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
											daughter0Positron = kTRUE;
										}
									}
								}			
							}					
						}
					}
				}
			}
					
			if(daughter1->GetNDaughters() >= 2){	 // second gamma
				for(Int_t TrackIndex=daughter1->GetFirstDaughter();TrackIndex<=daughter1->GetLastDaughter();TrackIndex++){
					TParticle *tmpDaughter = fStack->Particle(TrackIndex);
					if(tmpDaughter->GetUniqueID() == 5){
						if(tmpDaughter->GetPdgCode() == 11){
							if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
								if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
									if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1	){
										if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
											daughter1Electron = kTRUE;
										}
									}
								}			
							}
						} else if(tmpDaughter->GetPdgCode() == -11){
							if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
								if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
									if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1	){
										if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
											daughter1Positron = kTRUE;
										}
									}
								}			
							}					
						}
					}
				}
			}	
			if(particle->GetPdgCode()==111){		 //Pi0
				if( iTracks >= fStack->GetNprimary()){
					fHistograms->FillHistogram("MC_Pi0_Secondaries_Eta", particle->Eta());
					fHistograms->FillHistogram("MC_Pi0_Secondaries_Rapid", rapidity);
					fHistograms->FillHistogram("MC_Pi0_Secondaries_Phi", tmpPhi);
					fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt", particle->Pt());
					fHistograms->FillHistogram("MC_Pi0_Secondaries_Energy", particle->Energy());
					fHistograms->FillHistogram("MC_Pi0_Secondaries_R", particle->R());
					fHistograms->FillHistogram("MC_Pi0_Secondaries_ZR", particle->Vz(),particle->R());
					fHistograms->FillHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
					fHistograms->FillHistogram("MC_Pi0_Secondaries_XY", particle->Vx(),particle->Vy());//only fill from one daughter to avoid multiple filling
									
					if(gammaEtaCut && gammaRCut){	
						//if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
						fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
						fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
						if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
							fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
							fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
						}
						if(particle->GetMother(0) > -1){
							Int_t pdgPionMother=fStack->Particle(particle->GetMother(0))->GetPdgCode();
							if(pdgPionMother ==310 ){
								fHistograms->FillHistogram("MC_K0S_Pt_FromPi0", fStack->Particle(particle->GetMother(0))->Pt());
							}
						}
					}
				} else{

					fHistograms->FillHistogram("MC_Pi0_Eta", particle->Eta());	
					fHistograms->FillHistogram("MC_Pi0_Rapid", rapidity);
					fHistograms->FillHistogram("MC_Pi0_Phi", tmpPhi);
					fHistograms->FillHistogram("MC_Pi0_Pt", particle->Pt());
					//if(fGCMCEvent->IsFromBGEvent(iTracks)) fHistograms->FillHistogram("MC_Pi0_Pt_under", particle->Pt());
					fHistograms->FillHistogram("MC_Pi0_Pt_vs_Rapid", particle->Pt(),rapidity);
					fHistograms->FillHistogram("MC_Pi0_Energy", particle->Energy());
					fHistograms->FillHistogram("MC_Pi0_R", particle->R());
					fHistograms->FillHistogram("MC_Pi0_ZR", particle->Vz(),particle->R());
					fHistograms->FillHistogram("MC_Pi0_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
					fHistograms->FillHistogram("MC_Pi0_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
					if(TMath::Abs(particle->Eta())<0.9)fHistograms->FillHistogram("MC_Pi0_Pt_Fiducial", particle->Pt());
					if(particle->GetMother(0) > -1){
						Int_t pdgPionMother=fStack->Particle(particle->GetMother(0))->GetPdgCode();
						if(pdgPionMother ==221 ){
							fHistograms->FillHistogram("MC_Pi0_Pt_FromEta", particle->Pt());
						}else	if( pdgPionMother==223 ){
							fHistograms->FillHistogram("MC_Pi0_Pt_FromOmega", particle->Pt());
						}else if(TMath::Abs(pdgPionMother) ==213 ){
							fHistograms->FillHistogram("MC_Pi0_Pt_FromRhoPlusMinus", particle->Pt());
						}else{
							fHistograms->FillHistogram("MC_Pi0_Pt_FromOthers", particle->Pt());
						}
					}else{
						fHistograms->FillHistogram("MC_Pi0_Pt_Direct", particle->Pt());
					}

					if(gammaEtaCut && gammaRCut){
			//		if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
						fHistograms->FillHistogram("MC_Pi0_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
						fHistograms->FillHistogram("MC_Pi0_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
						if(TMath::Abs(particle->Eta())<0.9)fHistograms->FillHistogram("MC_Pi0_Pt_withinAcceptance_Fiducial", particle->Pt());
						if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
							fHistograms->FillHistogram("MC_Pi0_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
							fHistograms->FillHistogram("MC_Pi0_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
							fHistograms->FillHistogram("MC_Pi0_ZR_ConvGamma_withinAcceptance", particle->Vz(),particle->R());
							fHistograms->FillHistogram("MC_Pi0_ConvGamma_OpeningAngle_Pt", particle->Pt(),GetMCOpeningAngle(daughter0,daughter1));
							fHistograms->FillHistogram("MC_Pi0_ConvGamma_PtGamma_Pt", particle->Pt(),daughter0->Pt());
							fHistograms->FillHistogram("MC_Pi0_ConvGamma_PtGamma_Pt", particle->Pt(),daughter1->Pt());

							Double_t alfa=0.;
							if((daughter0->Energy()+daughter1->Energy()) > 0.){
								alfa= TMath::Abs((daughter0->Energy()-daughter1->Energy())/(daughter0->Energy()+daughter1->Energy()));
							}
							fHistograms->FillHistogram("MC_Pi0_alpha",alfa);
							if(TMath::Abs(particle->Eta())<0.9)fHistograms->FillHistogram("MC_Pi0_Pt_ConvGamma_withinAcceptance_Fiducial", particle->Pt());
						}
					}
				}
			}
					
			if(particle->GetPdgCode()==221){	 //Eta
				fHistograms->FillHistogram("MC_Eta_Eta", particle->Eta());
				fHistograms->FillHistogram("MC_Eta_Rapid", rapidity);
				fHistograms->FillHistogram("MC_Eta_Phi",tmpPhi);
				fHistograms->FillHistogram("MC_Eta_Pt", particle->Pt());
				fHistograms->FillHistogram("MC_Eta_Pt_vs_Rapid", particle->Pt(),rapidity);
				fHistograms->FillHistogram("MC_Eta_Energy", particle->Energy());
				fHistograms->FillHistogram("MC_Eta_R", particle->R());
				fHistograms->FillHistogram("MC_Eta_ZR", particle->Vz(),particle->R());
				fHistograms->FillHistogram("MC_Eta_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
				fHistograms->FillHistogram("MC_Eta_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
							
				if(gammaEtaCut && gammaRCut){	
					//	if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
					fHistograms->FillHistogram("MC_Eta_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
					fHistograms->FillHistogram("MC_Eta_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
					if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
						fHistograms->FillHistogram("MC_Eta_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
						fHistograms->FillHistogram("MC_Eta_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
						fHistograms->FillHistogram("MC_Eta_ZR_ConvGamma_withinAcceptance", particle->Vz(),particle->R());
						fHistograms->FillHistogram("MC_Eta_ConvGamma_OpeningAngle_Pt", particle->Pt(),GetMCOpeningAngle(daughter0,daughter1));
						fHistograms->FillHistogram("MC_Eta_ConvGamma_PtGamma_Pt", particle->Pt(),daughter0->Pt());
						fHistograms->FillHistogram("MC_Eta_ConvGamma_PtGamma_Pt", particle->Pt(),daughter1->Pt());

						Double_t alfa=0.;
						if((daughter0->Energy()+daughter1->Energy()) > 0.){
							alfa= TMath::Abs((daughter0->Energy()-daughter1->Energy())/(daughter0->Energy()+daughter1->Energy()));
						}
						fHistograms->FillHistogram("MC_Eta_alpha",alfa);
					}			
				}			
			}
					
			// all motherparticles with 2 gammas as daughters
			fHistograms->FillHistogram("MC_Mother_R", particle->R());
			fHistograms->FillHistogram("MC_Mother_ZR", particle->Vz(),particle->R());
			fHistograms->FillHistogram("MC_Mother_XY", particle->Vx(),particle->Vy());
			fHistograms->FillHistogram("MC_Mother_Mass", particle->GetCalcMass());
			fHistograms->FillHistogram("MC_Mother_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
			fHistograms->FillHistogram("MC_Mother_Energy", particle->Energy());
			fHistograms->FillHistogram("MC_Mother_Pt", particle->Pt());
			fHistograms->FillHistogram("MC_Mother_Eta", particle->Eta());
			fHistograms->FillHistogram("MC_Mother_Rapid", rapidity);
			fHistograms->FillHistogram("MC_Mother_Phi",tmpPhi);
			fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt",particle->GetMass(),particle->Pt());			
					
			if(gammaEtaCut && gammaRCut){	
				//			if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
				fHistograms->FillHistogram("MC_Mother_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
				fHistograms->FillHistogram("MC_Mother_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
				fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt_withinAcceptance",particle->GetMass(),particle->Pt());			
				if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
					fHistograms->FillHistogram("MC_Mother_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
					fHistograms->FillHistogram("MC_Mother_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
					fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt_ConvGamma_withinAcceptance",particle->GetMass(),particle->Pt());			
				}			
			}	// end passed R and eta cut
		} // end if(particle->GetNDaughters() == 2)
	}// end for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++)
	
	fHistograms->FillHistogram("ESD_TrueMult_vs_MeasMult", nCharged,fMultiplicity);
	fHistograms->FillHistogram("ESD_TrueMult_vs_MeasMult_pt150MeV", nCharged150MeV,fMultiplicity);
	//cout << "True   " <<  nCharged << "     Meas.     " << fMultiplicity << endl;

} // end ProcessMCData



void AliAnalysisTaskGammaConversion::FillNtuple(){
	//Fills the ntuple with the different values
	
	if(fGammaNtuple == NULL){
		return;
	}
	Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
	for(Int_t i=0;i<numberOfV0s;i++){
		Float_t values[27] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		AliESDv0 * cV0 = fV0Reader->GetV0(i);
		Double_t negPID=0;
		Double_t posPID=0;
		fV0Reader->GetPIDProbability(negPID,posPID);
		values[0]=cV0->GetOnFlyStatus();
		values[1]=fV0Reader->CheckForPrimaryVertex();
		values[2]=negPID;
		values[3]=posPID;
		values[4]=fV0Reader->GetX();
		values[5]=fV0Reader->GetY();
		values[6]=fV0Reader->GetZ();
		values[7]=fV0Reader->GetXYRadius();
		values[8]=fV0Reader->GetMotherCandidateNDF();
		values[9]=fV0Reader->GetMotherCandidateChi2();
		values[10]=fV0Reader->GetMotherCandidateEnergy();
		values[11]=fV0Reader->GetMotherCandidateEta();
		values[12]=fV0Reader->GetMotherCandidatePt();
		values[13]=fV0Reader->GetMotherCandidateMass();
		values[14]=fV0Reader->GetMotherCandidateWidth();
		//		values[15]=fV0Reader->GetMotherMCParticle()->Pt();	 MOVED TO THE END, HAS TO BE CALLED AFTER HasSameMother NB: still has the same entry in the array
		values[16]=fV0Reader->GetOpeningAngle();
		values[17]=fV0Reader->GetNegativeTrackEnergy();
		values[18]=fV0Reader->GetNegativeTrackPt();
		values[19]=fV0Reader->GetNegativeTrackEta();
		values[20]=fV0Reader->GetNegativeTrackPhi();
		values[21]=fV0Reader->GetPositiveTrackEnergy();
		values[22]=fV0Reader->GetPositiveTrackPt();
		values[23]=fV0Reader->GetPositiveTrackEta();
		values[24]=fV0Reader->GetPositiveTrackPhi();
		values[25]=fV0Reader->HasSameMCMother();
		if(values[25] != 0){
			values[26]=fV0Reader->GetMotherMCParticlePDGCode();
			values[15]=fV0Reader->GetMotherMCParticle()->Pt();
		}
		fTotalNumberOfAddedNtupleEntries++;
		fGammaNtuple->Fill(values);
	}
	fV0Reader->ResetV0IndexNumber();
	
}

void AliAnalysisTaskGammaConversion::ProcessV0sNoCut(){
	// Process all the V0's without applying any cuts to it
	
	Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
	for(Int_t i=0;i<numberOfV0s;i++){
		/*AliESDv0 * cV0 = */fV0Reader->GetV0(i);
		
		if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
			continue;
		}

		//		if( !fV0Reader->GetV0(i)->GetOnFlyStatus()){
		if( !fV0Reader->CheckV0FinderStatus(i)){
			continue;
		}


		if( !((fV0Reader->GetNegativeESDTrack())->GetStatus() & AliESDtrack::kTPCrefit) || 
				!((fV0Reader->GetPositiveESDTrack())->GetStatus() & AliESDtrack::kTPCrefit) ){
			continue;
		}

		if( fV0Reader->GetNegativeESDTrack()->GetSign()== fV0Reader->GetPositiveESDTrack()->GetSign()){
			continue;
		}

		if( fV0Reader->GetNegativeESDTrack()->GetKinkIndex(0) > 0 || 
			fV0Reader->GetPositiveESDTrack()->GetKinkIndex(0) > 0) {			
			continue;
		}
		if(TMath::Abs(fV0Reader->GetMotherCandidateEta())> fV0Reader->GetEtaCut() || TMath::Abs(fV0Reader->GetMotherCandidateEta())< fV0Reader->GetEtaCutMin()){
			continue;
		}
		if(TMath::Abs(fV0Reader->GetPositiveTrackEta())> fV0Reader->GetEtaCut() || TMath::Abs(fV0Reader->GetPositiveTrackEta())< fV0Reader->GetEtaCutMin()){
			continue;
		}
		if(TMath::Abs(fV0Reader->GetNegativeTrackEta())> fV0Reader->GetEtaCut() || TMath::Abs(fV0Reader->GetNegativeTrackEta())< fV0Reader->GetEtaCutMin()){
			continue;
		}
		if( fV0Reader->GetXYRadius() <= ((TMath::Abs(fV0Reader->GetZ())*fV0Reader->GetLineCutZRSlope())-fV0Reader->GetLineCutZValue()) ){ // cuts out regions where we do not reconstruct
			continue; 
		} else if ( fV0Reader->GetEtaCutMin() != -0.1 && fV0Reader->GetXYRadius() > ((TMath::Abs(fV0Reader->GetZ())*fV0Reader->GetLineCutZRSlopeMin())-fV0Reader->GetLineCutZValueMin()) ) {
			continue; 
		}

		fHistograms->FillHistogram("ESD_NoCutAllV0_Pt", fV0Reader->GetMotherCandidatePt());

		// RRnewTOF start ///////////////////////////////////////////////
		UInt_t statusPos = fV0Reader->GetPositiveESDTrack()->GetStatus();
		UInt_t statusNeg = fV0Reader->GetNegativeESDTrack()->GetStatus(); 

		Double_t t0pos = fV0Reader->GetESDpid()->GetTOFResponse().GetStartTime(fV0Reader->GetPositiveTrackP());
		Double_t t0neg = fV0Reader->GetESDpid()->GetTOFResponse().GetStartTime(fV0Reader->GetNegativeTrackP());

		Double_t timesPos[5];
		fV0Reader->GetPositiveESDTrack()->GetIntegratedTimes(timesPos);
		Double_t timesNeg[5];
		fV0Reader->GetNegativeESDTrack()->GetIntegratedTimes(timesNeg);

		Double_t TOFsignalPos =	fV0Reader->GetPositiveTrackTOFsignal();
		Double_t TOFsignalNeg =	fV0Reader->GetNegativeTrackTOFsignal();

		Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
		Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];

		if( (statusPos & AliESDtrack::kTOFpid) && !(statusPos & AliESDtrack::kTOFmismatch) ) fHistograms->FillHistogram("ESD_WOCutConvGamma_EandP_P_dT", fV0Reader->GetPositiveTrackP(), dTpos);
		if( (statusNeg & AliESDtrack::kTOFpid) && !(statusNeg & AliESDtrack::kTOFmismatch) ) fHistograms->FillHistogram("ESD_WOCutConvGamma_EandP_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);
		// RRnewTOF end /////////////////////////////////////////////////
		
		if(fDoMCTruth){
			
			if(fV0Reader->HasSameMCMother() == kFALSE){
				continue;
			}
			
			TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
			TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
			
			if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
				continue;
			}
			if(negativeMC->GetPdgCode() == positiveMC->GetPdgCode()){
				continue;
			}

			if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5){ // id 5 is conversion
				continue;
			}
			
			if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){				
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Eta", fV0Reader->GetMotherCandidateEta());				
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
							
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Pt_Chi2", fV0Reader->GetMotherCandidatePt(), fV0Reader->GetMotherCandidateChi2());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_Eta_Chi2", fV0Reader->GetMotherCandidateEta(), fV0Reader->GetMotherCandidateChi2());
							
				fHistograms->FillHistogram("ESD_NoCutConversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
				fHistograms->FillHistogram("ESD_NoCutConversion_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_NoCutConversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_NoCutConversion_OpeningAngle", fV0Reader->GetOpeningAngle());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_CosPointingAngle", fV0Reader->GetCosPointingAngle());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_DcaDaughters", fV0Reader->GetDcaDaughters());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_NormDcaDistDaughters", fV0Reader->GetNormDcaDistDaughters());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_LikelihoodAP", fV0Reader->GetLikelihoodAP());

				fHistograms->FillHistogram("ESD_NoCutConvGamma_E_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetNegativeTrackP()/fV0Reader->GetMotherCandidateP());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_P_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetPositiveTrackP()/fV0Reader->GetMotherCandidateP());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_E_dEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetNegativeTrackTPCdEdx());
				fHistograms->FillHistogram("ESD_NoCutConvGamma_P_dEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetPositiveTrackTPCdEdx());

				//store MCTruth properties
				fHistograms->FillHistogram("ESD_NoCutConvGamma_MC_Pt_Eta", fV0Reader->GetMotherMCParticle()->Pt(),fV0Reader->GetMotherMCParticle()->Eta());
				fHistograms->FillHistogram("ESD_NoCutConversion_MC_ZR", negativeMC->Vz(),negativeMC->R());
				fHistograms->FillHistogram("ESD_NoCutConversion_MC_XY", negativeMC->Vx(),negativeMC->Vy());
			}
		}
	}
	fV0Reader->ResetV0IndexNumber();
}

void AliAnalysisTaskGammaConversion::ProcessV0s(){
	// see header file for documentation


	if(fWriteNtuple == kTRUE){
		FillNtuple();
	}
	
	Int_t nSurvivingV0s=0;
	fV0Reader->ResetNGoodV0s();
	while(fV0Reader->NextV0()){
		nSurvivingV0s++;
		

		TVector3 vtxConv(fV0Reader->GetX(),fV0Reader->GetY(), fV0Reader->GetZ());
	
		//-------------------------- filling v0 information -------------------------------------
		fHistograms->FillHistogram("ESD_Conversion_R", fV0Reader->GetXYRadius());
		fHistograms->FillHistogram("ESD_Conversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
		fHistograms->FillHistogram("ESD_Conversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
		fHistograms->FillHistogram("ESD_Conversion_OpeningAngle", fV0Reader->GetOpeningAngle());		
		fHistograms->FillHistogram("ESD_ConversionMapping_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
		fHistograms->FillHistogram("ESD_ConversionMapping_ZPhi", fV0Reader->GetZ(),vtxConv.Phi());
		fHistograms->FillHistogram("ESD_ConversionMapping_RPhi", fV0Reader->GetXYRadius(),vtxConv.Phi());
			
		if ( fV0Reader->GetMotherCandidatePt() >= 2.){
			fHistograms->FillHistogram("ESD_Conversion_minPt_R", fV0Reader->GetXYRadius());
			fHistograms->FillHistogram("ESD_Conversion_minPt_Z", fV0Reader->GetZ());
		}

		// Specific histograms for beam pipe studies
		if( TMath::Abs(fV0Reader->GetZ()) < fV0Reader->GetLineCutZValue() && TMath::Abs(fV0Reader->GetZ()) > fV0Reader->GetLineCutZValueMin()){
			fHistograms->FillHistogram("ESD_Conversion_XY_BeamPipe", fV0Reader->GetX(),fV0Reader->GetY());
			fHistograms->FillHistogram("ESD_Conversion_RPhi_BeamPipe", vtxConv.Phi(),fV0Reader->GetXYRadius());
		}

	
		fHistograms->FillHistogram("ESD_E_Energy", fV0Reader->GetNegativeTrackEnergy());
		fHistograms->FillHistogram("ESD_E_Pt", fV0Reader->GetNegativeTrackPt());
		fHistograms->FillHistogram("ESD_E_Eta", fV0Reader->GetNegativeTrackEta());
		fHistograms->FillHistogram("ESD_E_Phi", fV0Reader->GetNegativeTrackPhi());
		fHistograms->FillHistogram("ESD_E_nTPCClusters", fV0Reader->GetNegativeTracknTPCClusters());
		fHistograms->FillHistogram("ESD_E_nITSClusters", fV0Reader->GetNegativeTracknITSClusters());
		if ( fV0Reader->GetNegativeTrackPt()> 0.150){
			fHistograms->FillHistogram("ESD_E_minPt_nTPCClusters", fV0Reader->GetNegativeTracknTPCClusters());
		}
		if ( fV0Reader->GetNegativeTracknITSClusters()==0){
			fHistograms->FillHistogram("ESD_E_onlyTPC_nTPCClusters", fV0Reader->GetNegativeTracknTPCClusters());
		}

		Double_t eClsToF= 0;
		if(!fV0Reader->GetUseCorrectedTPCClsInfo()){
			if(fV0Reader->GetNegativeTracknTPCFClusters()!=0 ){
				eClsToF=(Double_t)fV0Reader->GetNegativeTracknTPCClusters()/(Double_t)fV0Reader->GetNegativeTracknTPCFClusters();
			}
		} else {
			eClsToF= fV0Reader->GetNegativeESDTrack()->GetTPCClusterInfo(2,0,fV0Reader->GetFirstTPCRow(fV0Reader->GetXYRadius()));
		}
			fHistograms->FillHistogram("ESD_E_nTPCClustersToFP", fV0Reader->GetNegativeTrackP(),eClsToF );
			fHistograms->FillHistogram("ESD_E_nTPCClustersToFR", fV0Reader->GetXYRadius(),eClsToF );

		if(fV0Reader->GetNegativeTracknTPCClusters()!=0 ){
			fHistograms->FillHistogram("ESD_E_TPCchi2", fV0Reader->GetNegativeTrackTPCchi2()/(Double_t)fV0Reader->GetNegativeTracknTPCClusters());
		}



		fHistograms->FillHistogram("ESD_P_Energy", fV0Reader->GetPositiveTrackEnergy());
		fHistograms->FillHistogram("ESD_P_Pt", fV0Reader->GetPositiveTrackPt());
		fHistograms->FillHistogram("ESD_P_Eta", fV0Reader->GetPositiveTrackEta());
		fHistograms->FillHistogram("ESD_P_Phi", fV0Reader->GetPositiveTrackPhi());
		fHistograms->FillHistogram("ESD_P_nTPCClusters", fV0Reader->GetPositiveTracknTPCClusters());
		fHistograms->FillHistogram("ESD_P_nITSClusters", fV0Reader->GetPositiveTracknITSClusters());
		if ( fV0Reader->GetPositiveTrackPt()> 0.150){
			fHistograms->FillHistogram("ESD_P_minPt_nTPCClusters", fV0Reader->GetPositiveTracknTPCClusters());
		}
		if (fV0Reader->GetPositiveTracknITSClusters()==0){
			fHistograms->FillHistogram("ESD_P_onlyTPC_nTPCClusters", fV0Reader->GetPositiveTracknTPCClusters());
		}

		Double_t pClsToF= 0;
		if(!fV0Reader->GetUseCorrectedTPCClsInfo()){
			if(fV0Reader->GetPositiveTracknTPCFClusters()!=0){
				pClsToF = (Double_t)fV0Reader->GetPositiveTracknTPCClusters()/(Double_t)fV0Reader->GetPositiveTracknTPCFClusters();
			}
		} else {
			pClsToF= fV0Reader->GetPositiveESDTrack()->GetTPCClusterInfo(2,0,fV0Reader->GetFirstTPCRow(fV0Reader->GetXYRadius()));
		}

		fHistograms->FillHistogram("ESD_P_nTPCClustersToFP",fV0Reader->GetPositiveTrackP(), pClsToF);
		fHistograms->FillHistogram("ESD_P_nTPCClustersToFR",fV0Reader->GetXYRadius(), pClsToF);

		if(fV0Reader->GetPositiveTracknTPCClusters()!=0){
			fHistograms->FillHistogram("ESD_P_TPCchi2", fV0Reader->GetPositiveTrackTPCchi2()/(Double_t)fV0Reader->GetPositiveTracknTPCClusters());
		}



		fHistograms->FillHistogram("ESD_ConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
		fHistograms->FillHistogram("ESD_ConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
		fHistograms->FillHistogram("ESD_ConvGamma_Eta", fV0Reader->GetMotherCandidateEta());
		fHistograms->FillHistogram("ESD_ConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
		fHistograms->FillHistogram("ESD_ConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
		fHistograms->FillHistogram("ESD_ConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
		fHistograms->FillHistogram("ESD_ConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
		fHistograms->FillHistogram("ESD_ConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
		fHistograms->FillHistogram("ESD_ConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
		fHistograms->FillHistogram("ESD_ConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
		
		fHistograms->FillHistogram("ESD_ConvGamma_Pt_Chi2", fV0Reader->GetMotherCandidatePt(), fV0Reader->GetMotherCandidateChi2());
		fHistograms->FillHistogram("ESD_ConvGamma_Eta_Chi2", fV0Reader->GetMotherCandidateEta(), fV0Reader->GetMotherCandidateChi2());
		
		fHistograms->FillHistogram("ESD_ConvGamma_CosPointingAngle", fV0Reader->GetCosPointingAngle());
		fHistograms->FillHistogram("ESD_ConvGamma_DcaDaughters", fV0Reader->GetDcaDaughters());
		fHistograms->FillHistogram("ESD_ConvGamma_NormDcaDistDaughters", fV0Reader->GetNormDcaDistDaughters());
		fHistograms->FillHistogram("ESD_ConvGamma_LikelihoodAP", fV0Reader->GetLikelihoodAP());
		
		fHistograms->FillHistogram("ESD_ConvGamma_E_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetNegativeTrackP()/fV0Reader->GetMotherCandidateP());
		fHistograms->FillHistogram("ESD_ConvGamma_P_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetPositiveTrackP()/fV0Reader->GetMotherCandidateP());
		fHistograms->FillHistogram("ESD_ConvGamma_E_dEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetNegativeTrackTPCdEdx());
		fHistograms->FillHistogram("ESD_ConvGamma_P_dEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetPositiveTrackTPCdEdx());
		fHistograms->FillHistogram("ESD_ConvGamma_E_SigdEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetNegativeESDTrack(),AliPID::kElectron));
		fHistograms->FillHistogram("ESD_ConvGamma_P_SigdEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetPositiveESDTrack(),AliPID::kElectron));
		fHistograms->FillHistogram("ESD_ConvGamma_PiPl_SigdEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetNegativeESDTrack(),AliPID::kPion));
		fHistograms->FillHistogram("ESD_ConvGamma_PiMi_SigdEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetPositiveESDTrack(),AliPID::kPion));
		fHistograms->FillHistogram("ESD_ConvGamma_KPl_SigdEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetNegativeESDTrack(),AliPID::kKaon));
		fHistograms->FillHistogram("ESD_ConvGamma_KMi_SigdEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetPositiveESDTrack(),AliPID::kKaon));
		fHistograms->FillHistogram("ESD_ConvGamma_PPl_SigdEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetNegativeESDTrack(),AliPID::kProton));
		fHistograms->FillHistogram("ESD_ConvGamma_PMi_SigdEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetPositiveESDTrack(),AliPID::kProton));
		fHistograms->FillHistogram("ESD_ConvGamma_MuPl_SigdEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetNegativeESDTrack(),AliPID::kMuon));
		fHistograms->FillHistogram("ESD_ConvGamma_MuMi_SigdEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetESDpid()->NumberOfSigmasTPC(fV0Reader->GetPositiveESDTrack(),AliPID::kMuon));
		

		UInt_t statusPos = fV0Reader->GetPositiveESDTrack()->GetStatus(); //moved up here from below RRnewTOF
		UInt_t statusNeg = fV0Reader->GetNegativeESDTrack()->GetStatus(); 
		// RRnewTOF start ///////////////////////////////////////////////
		Double_t t0pos = fV0Reader->GetESDpid()->GetTOFResponse().GetStartTime(fV0Reader->GetPositiveTrackP());
		Double_t t0neg = fV0Reader->GetESDpid()->GetTOFResponse().GetStartTime(fV0Reader->GetNegativeTrackP());

		Double_t timesPos[5];
		fV0Reader->GetPositiveESDTrack()->GetIntegratedTimes(timesPos);
		Double_t timesNeg[5];
		fV0Reader->GetNegativeESDTrack()->GetIntegratedTimes(timesNeg);

		Double_t TOFsignalPos =	fV0Reader->GetPositiveTrackTOFsignal();
		Double_t TOFsignalNeg =	fV0Reader->GetNegativeTrackTOFsignal();

		Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
		Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];

		if( (statusPos & AliESDtrack::kTOFpid) && !(statusPos & AliESDtrack::kTOFmismatch) ) fHistograms->FillHistogram("ESD_ConvGamma_EandP_P_dT", fV0Reader->GetPositiveTrackP(), dTpos);
		if( (statusNeg & AliESDtrack::kTOFpid) && !(statusNeg & AliESDtrack::kTOFmismatch) ) fHistograms->FillHistogram("ESD_ConvGamma_EandP_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);
		// RRnewTOF end /////////////////////////////////////////////////

		Double_t negPID=0;
		Double_t posPID=0;
		fV0Reader->GetPIDProbability(negPID,posPID);
		fHistograms->FillHistogram("ESD_ConvGamma_E_EProbP",fV0Reader->GetNegativeTrackP(),negPID);
		fHistograms->FillHistogram("ESD_ConvGamma_P_EProbP",fV0Reader->GetPositiveTrackP(),posPID);

		Double_t negPIDmupi=0;
		Double_t posPIDmupi=0;
		fV0Reader->GetPIDProbabilityMuonPion(negPIDmupi,posPIDmupi);
		fHistograms->FillHistogram("ESD_ConvGamma_E_mupiProbP",fV0Reader->GetNegativeTrackP(),negPIDmupi);
		fHistograms->FillHistogram("ESD_ConvGamma_P_mupiProbP",fV0Reader->GetPositiveTrackP(),posPIDmupi);

		Double_t armenterosQtAlfa[2];
		fV0Reader->GetArmenterosQtAlfa(fV0Reader-> GetNegativeKFParticle(), 
					 fV0Reader-> GetPositiveKFParticle(), 
					 fV0Reader->GetMotherCandidateKFCombination(),
					 armenterosQtAlfa);
	 
		fHistograms->FillHistogram("ESD_ConvGamma_alfa_qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
		fHistograms->FillHistogram("ESD_ConvGamma_Pt_Qt",fV0Reader->GetMotherCandidatePt(),armenterosQtAlfa[0]);
	
		if(!fV0Reader->GetIsHeavyIon()){
			fHistograms->FillHistogram("3DPlots_Conversion_XYZ", fV0Reader->GetX(),fV0Reader->GetY(),fV0Reader->GetZ());		
			//fHistograms->FillHistogram("3DPlots_Conversion_ZRPhi", fV0Reader->GetZ(),fV0Reader->GetXYRadius(), vtxConv.Phi());

			// begin mapping
			Int_t rBin		= fHistograms->GetRBin(fV0Reader->GetXYRadius());
			Int_t zBin		= fHistograms->GetZBin(fV0Reader->GetZ());
			Int_t phiBin	= fHistograms->GetPhiBin(fV0Reader->GetNegativeTrackPhi());
			Double_t rFMD=25;
			Double_t rITSTPCMin=45;
			Double_t rITSTPCMax=80;
			
			//		Double_t motherCandidateEta= fV0Reader->GetMotherCandidateEta();
			
			TString nameESDMappingPhiR="";
			nameESDMappingPhiR.Form("ESD_Conversion_Mapping_Phi%02d_R%02d",phiBin,rBin);
			//fHistograms->FillHistogram(nameESDMappingPhiR, fV0Reader->GetZ(), motherCandidateEta);
			
			TString nameESDMappingPhi="";
			nameESDMappingPhi.Form("ESD_Conversion_Mapping_Phi%02d",phiBin);
			//fHistograms->FillHistogram(nameESDMappingPhi, fV0Reader->GetZ(), motherCandidateEta);
			
			TString nameESDMappingR="";
			nameESDMappingR.Form("ESD_Conversion_Mapping_R%02d",rBin);
			//fHistograms->FillHistogram(nameESDMappingR, fV0Reader->GetZ(), motherCandidateEta);	
			
			TString nameESDMappingPhiInR="";
			nameESDMappingPhiInR.Form("ESD_Conversion_Mapping_Phi_in_R_%02d",rBin);
			//		fHistograms->FillHistogram(nameESDMappingPhiInR, fV0Reader->GetMotherCandidatePhi());
			fHistograms->FillHistogram(nameESDMappingPhiInR, vtxConv.Phi());
			
			TString nameESDMappingZInR="";
			nameESDMappingZInR.Form("ESD_Conversion_Mapping_Z_in_R_%02d",rBin);
			fHistograms->FillHistogram(nameESDMappingZInR, fV0Reader->GetZ());
			
			TString nameESDMappingPhiInZ="";
			nameESDMappingPhiInZ.Form("ESD_Conversion_Mapping_Phi_in_Z_%02d",zBin);
			//		fHistograms->FillHistogram(nameESDMappingPhiInR, fV0Reader->GetMotherCandidatePhi());
			fHistograms->FillHistogram(nameESDMappingPhiInZ, vtxConv.Phi());
			
			if(fV0Reader->GetXYRadius()<rFMD){
				TString nameESDMappingFMDPhiInZ="";
				nameESDMappingFMDPhiInZ.Form("ESD_Conversion_Mapping_FMD_Phi_in_Z_%02d",zBin);
				fHistograms->FillHistogram(nameESDMappingFMDPhiInZ, vtxConv.Phi());
				fHistograms->FillHistogram("ESD_ConversionMapping_FMD_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
			}
			if(fV0Reader->GetXYRadius()>rFMD && fV0Reader->GetXYRadius()<rITSTPCMin){
				fHistograms->FillHistogram("ESD_ConversionMapping_FMD2_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
			}

			if(fV0Reader->GetXYRadius()>rITSTPCMin && fV0Reader->GetXYRadius()<rITSTPCMax){
				TString nameESDMappingITSTPCPhiInZ="";
				nameESDMappingITSTPCPhiInZ.Form("ESD_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",zBin);
				fHistograms->FillHistogram(nameESDMappingITSTPCPhiInZ, vtxConv.Phi());
				fHistograms->FillHistogram("ESD_ConversionMapping_ITSTPC_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
			}
			
			TString nameESDMappingRInZ="";
			nameESDMappingRInZ.Form("ESD_Conversion_Mapping_R_in_Z_%02d",zBin);
			fHistograms->FillHistogram(nameESDMappingRInZ, fV0Reader->GetXYRadius());
			
			if(fV0Reader->GetMotherCandidatePt() > fLowPtMapping && fV0Reader->GetMotherCandidatePt()< fHighPtMapping){
				TString nameESDMappingMidPtPhiInR="";
				nameESDMappingMidPtPhiInR.Form("ESD_Conversion_Mapping_MidPt_Phi_in_R_%02d",rBin);
				fHistograms->FillHistogram(nameESDMappingMidPtPhiInR, vtxConv.Phi());
				
				TString nameESDMappingMidPtZInR="";
				nameESDMappingMidPtZInR.Form("ESD_Conversion_Mapping_MidPt_Z_in_R_%02d",rBin);
				fHistograms->FillHistogram(nameESDMappingMidPtZInR, fV0Reader->GetZ());
				
				TString nameESDMappingMidPtPhiInZ="";
				nameESDMappingMidPtPhiInZ.Form("ESD_Conversion_Mapping_MidPt_Phi_in_Z_%02d",zBin);
				fHistograms->FillHistogram(nameESDMappingMidPtPhiInZ, vtxConv.Phi());
				if(fV0Reader->GetXYRadius()<rFMD){
					TString nameESDMappingMidPtFMDPhiInZ="";
					nameESDMappingMidPtFMDPhiInZ.Form("ESD_Conversion_Mapping_MidPt_FMD_Phi_in_Z_%02d",zBin);
					fHistograms->FillHistogram(nameESDMappingMidPtFMDPhiInZ, vtxConv.Phi());
				}
				TString nameESDMappingMidPtRInZ="";
				nameESDMappingMidPtRInZ.Form("ESD_Conversion_Mapping_MidPt_R_in_Z_%02d",zBin);
				fHistograms->FillHistogram(nameESDMappingMidPtRInZ, fV0Reader->GetXYRadius());	
			}
		}

		// end mapping
		
                   
		new((*fKFReconstructedGammasTClone)[fKFReconstructedGammasTClone->GetEntriesFast()])	AliKFConversionPhoton(fV0Reader);
	     
		 
		//----------------------------------- checking for "real" conversions (MC match) --------------------------------------
		if(fDoMCTruth){
			TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
			TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
			Double_t rFMD=25;
			Double_t rITSTPCMin=45;
			Double_t rITSTPCMax=80;
 
			if(fV0Reader->HasSameMCMother() == kFALSE){
				fHistograms->FillHistogram("ESD_TrueConvCombinatorial_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConvCombinatorial_Z", fV0Reader->GetZ());
				fHistograms->FillHistogram("ESD_TrueConvCombSelected_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				if ( fV0Reader->GetMotherCandidatePt() > 2. ) {
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialMinPt_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialMinPt_Z", fV0Reader->GetZ());			
				}
				fHistograms->FillHistogram("ESD_TrueConvCombinatorial_Pt", fV0Reader->GetMotherCandidatePt());
				fHistograms->FillHistogram("ESD_TrueConvCombinatorialDaughter_Pt", negativeMC->Pt(),positiveMC->Pt());
				if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialElec_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialElec_Pt", fV0Reader->GetMotherCandidatePt());
				}
				if(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==211){
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPi_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPi_Pt", fV0Reader->GetMotherCandidatePt());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPiDaughter_Pt", negativeMC->Pt(),positiveMC->Pt());
				}
				if((TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==2211) ||
					(TMath::Abs(negativeMC->GetPdgCode())==2212 && TMath::Abs(positiveMC->GetPdgCode())==211)){
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPiP_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPiP_Pt", fV0Reader->GetMotherCandidatePt());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialPiPDaughter_Pt", negativeMC->Pt(),positiveMC->Pt());
				}
				if((TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==211) ||
					(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==11)){
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialElecPi_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvCombinatorialElecPi_Pt", fV0Reader->GetMotherCandidatePt());
				}
				if( (statusPos & AliESDtrack::kTOFpid) && ( TMath::Abs(positiveMC->GetPdgCode()) != 11 ) && !(statusPos & AliESDtrack::kTOFmismatch) )
					fHistograms->FillHistogram("ESD_TrueConvCombinatorial_DaughtersNotElec_P_dT", fV0Reader->GetPositiveTrackP(), dTpos);//RRnewTOF
				if( (statusNeg & AliESDtrack::kTOFpid) && ( TMath::Abs(negativeMC->GetPdgCode()) != 11 ) && !(statusNeg & AliESDtrack::kTOFmismatch) )
					fHistograms->FillHistogram("ESD_TrueConvCombinatorial_DaughtersNotElec_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);//RRnewTOF
				continue;
			} 

						// Moved up to check true electron background
						//			TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
						//			TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
						
			if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
				fHistograms->FillHistogram("ESD_TrueConvHadronicBck_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConvHadronicBck_Z", fV0Reader->GetZ());
				if ( fV0Reader->GetMotherCandidatePt() > 2. ) {
					fHistograms->FillHistogram("ESD_TrueConvHadronicBckMinPt_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvHadronicBckMinPt_Z", fV0Reader->GetZ());			
				}
				fHistograms->FillHistogram("ESD_TrueConvHadronicBck_Pt", fV0Reader->GetMotherCandidatePt());
				fHistograms->FillHistogram("ESD_TrueConvHadronicBckDaughter_Pt", negativeMC->Pt(),positiveMC->Pt());
				if( (statusPos & AliESDtrack::kTOFpid) && !(statusPos & AliESDtrack::kTOFmismatch) )
					fHistograms->FillHistogram("ESD_TrueConvHadronicBck_Daughters_P_dT", fV0Reader->GetPositiveTrackP(), dTpos);//RRnewTOF
				if( (statusNeg & AliESDtrack::kTOFpid) && !(statusNeg & AliESDtrack::kTOFmismatch) )
					fHistograms->FillHistogram("ESD_TrueConvHadronicBck_Daughters_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);//RRnewTOF
				if((TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==2211) ||
					(TMath::Abs(negativeMC->GetPdgCode())==2212 && TMath::Abs(positiveMC->GetPdgCode())==211)){
					fHistograms->FillHistogram("ESD_TrueConvLambda_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvLambda_Pt", fV0Reader->GetMotherCandidatePt());
				}
				if(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==211){
					fHistograms->FillHistogram("ESD_TrueConvMeson_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvMeson_Pt", fV0Reader->GetMotherCandidatePt());
				}
				continue;
			}
			 

			if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
				continue;
			}

			//UInt_t statusPos = fV0Reader->GetPositiveESDTrack()->GetStatus(); moved higher
			//UInt_t statusNeg = fV0Reader->GetNegativeESDTrack()->GetStatus(); 
			UChar_t itsPixelPos = fV0Reader->GetPositiveESDTrack()->GetITSClusterMap();
			UChar_t itsPixelNeg = fV0Reader->GetNegativeESDTrack()->GetITSClusterMap();

			// Using the UniqueID Phojet does not get the Dalitz right
			//			if( (negativeMC->GetUniqueID() == 4 && positiveMC->GetUniqueID() ==4) ||
			//	(negativeMC->GetUniqueID() == 0 && positiveMC->GetUniqueID() ==0) ){// fill r distribution for Dalitz decays 
			if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 111){ //pi0
				fHistograms->FillHistogram("ESD_TrueDalitzContamination_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_Z", fV0Reader->GetZ());		
				if ( fV0Reader->GetMotherCandidatePt() > 2. ) {
					fHistograms->FillHistogram("ESD_TrueConvDalitzPi0MinPt_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvDalitzPi0MinPt_Z", fV0Reader->GetZ());			
				}
			
				//--------Histos for HFE 

				if(statusPos & AliESDtrack::kTOFpid){
					fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_SinglePos_R", fV0Reader->GetXYRadius());
					if( TESTBIT(itsPixelPos, 0) ){ 
						fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_SinglePos_kFirst_R", fV0Reader->GetXYRadius());
					}
				}
				if(statusNeg & AliESDtrack::kTOFpid){
					fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_SingleNeg_R", fV0Reader->GetXYRadius());
					if( TESTBIT(itsPixelNeg, 0) ){ 
						fHistograms->FillHistogram("ESD_TrueConvDalitzPi0_SingleNeg_kFirst_R", fV0Reader->GetXYRadius());
					}
				}
				//--------------------------------------------------------

			}
			if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 221){ //eta
				fHistograms->FillHistogram("ESD_TrueConvDalitzEta_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConvDalitzEta_Z", fV0Reader->GetZ());				
				if ( fV0Reader->GetMotherCandidatePt() > 2. ) {
					fHistograms->FillHistogram("ESD_TrueConvDalitzEtaMinPt_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvDalitzEtaMinPt_Z", fV0Reader->GetZ());			
				}
			}

			//}

			if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5){// check if the daughters come from a conversion
				continue;
			}
			
			if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
				if(fDoCF){
					Double_t containerInput[3];
					containerInput[0] = fV0Reader->GetMotherCandidatePt();
					containerInput[1] = fV0Reader->GetMotherCandidateEta();
					containerInput[2] = fV0Reader->GetMotherCandidateMass();
					fCFManager->GetParticleContainer()->Fill(containerInput,kStepTrueGamma); // for CF 
				}

				// RRnewTOF start ///////////////////////////////////////////////
				if( (statusPos & AliESDtrack::kTOFpid) && !(statusPos & AliESDtrack::kTOFmismatch) ) {
					fHistograms->FillHistogram("ESD_TrueConvGamma_EandP_P_dT", fV0Reader->GetPositiveTrackP(), dTpos);
				}
				if( (statusNeg & AliESDtrack::kTOFpid) && !(statusNeg & AliESDtrack::kTOFmismatch) ) {
					fHistograms->FillHistogram("ESD_TrueConvGamma_EandP_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);
				}
								// RRnewTOF end /////////////////////////////////////////////////
				if (fV0Reader->HasSameMCMother() == kTRUE){
					fHistograms->FillHistogram("ESD_TrueConvGammaSelected_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
					fHistograms->FillHistogram("ESD_TrueConvGammaSelected_Pt_Qt",fV0Reader->GetMotherCandidatePt(),armenterosQtAlfa[0]);	
				}
				// RRnewTOF end /////////////////////////////////////////////////

				fHistograms->FillHistogram("ESD_TrueConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
				if(negativeMC->GetMother(0) <= fStack->GetNprimary()){ // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
					fHistograms->FillHistogram("ESD_TrueConvPrimaryGamma_Pt", fV0Reader->GetMotherCandidatePt());
					fHistograms->FillHistogram("ESD_TrueConvPrimaryGamma_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConvPrimaryGamma_Z", fV0Reader->GetZ());				
					if(  fV0Reader->GetMotherCandidatePt() > 2. ) {
						fHistograms->FillHistogram("ESD_TrueConvPrimaryGammaMinPt_R", fV0Reader->GetXYRadius());
						fHistograms->FillHistogram("ESD_TrueConvPrimaryGammaMinPt_Z", fV0Reader->GetZ());				
					}
				}
				if(fV0Reader->GetMotherMCParticle()->GetMother(0) > -1){
					if(fStack->Particle(fV0Reader->GetMotherMCParticle()->GetMother(0))->GetPdgCode() == 221){ // Use just gamma from eta for ratio esdtruegamma / mcconvgamma
						fHistograms->FillHistogram("ESD_TrueConvEtaGamma_Pt", fV0Reader->GetMotherCandidatePt());
					}
				}
				fHistograms->FillHistogram("ESD_TrueConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Eta", fV0Reader->GetMotherCandidateEta());				
				fHistograms->FillHistogram("ESD_TrueConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
				fHistograms->FillHistogram("ESD_TrueConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
				fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLength",fV0Reader->GetNegativeNTPCClusters());
				fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLength",fV0Reader->GetPositiveNTPCClusters());
				fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass",fV0Reader->GetNegativeNTPCClusters(),fV0Reader->GetMotherCandidateMass());
				fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass",fV0Reader->GetPositiveNTPCClusters(),fV0Reader->GetMotherCandidateMass());
							
				fHistograms->FillHistogram("ESD_TrueConvGamma_Pt_Chi2", fV0Reader->GetMotherCandidatePt(), fV0Reader->GetMotherCandidateChi2());
				fHistograms->FillHistogram("ESD_TrueConvGamma_Eta_Chi2", fV0Reader->GetMotherCandidateEta(), fV0Reader->GetMotherCandidateChi2());
				if ( fV0Reader->GetMotherCandidatePt() > 2. ) {
					fHistograms->FillHistogram("ESD_TrueConversionMinPt_R", fV0Reader->GetXYRadius());
					fHistograms->FillHistogram("ESD_TrueConversionMinPt_Z", fV0Reader->GetZ());			
				}
							
				fHistograms->FillHistogram("ESD_TrueConversion_E_nTPCClustersToFR", fV0Reader->GetXYRadius(),eClsToF );
				fHistograms->FillHistogram("ESD_TrueConversion_P_nTPCClustersToFR",fV0Reader->GetXYRadius(), pClsToF);

				fHistograms->FillHistogram("ESD_TrueConversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
				fHistograms->FillHistogram("ESD_TrueConversion_R", fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConversion_Z", fV0Reader->GetZ());				
				fHistograms->FillHistogram("ESD_TrueConversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConversionMapping_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("ESD_TrueConversionMapping_ZPhi", fV0Reader->GetZ(),vtxConv.Phi());
				fHistograms->FillHistogram("ESD_TrueConversionMapping_RPhi", fV0Reader->GetXYRadius(),vtxConv.Phi());
				if(fV0Reader->GetXYRadius()<rFMD){
					fHistograms->FillHistogram("ESD_TrueConversionMapping_FMD_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
				}
				if(fV0Reader->GetXYRadius()>rFMD && fV0Reader->GetXYRadius()<rITSTPCMin){
					fHistograms->FillHistogram("ESD_TrueConversionMapping_FMD2_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
				}
				if(fV0Reader->GetXYRadius()>rITSTPCMin && fV0Reader->GetXYRadius()<rITSTPCMax){
					fHistograms->FillHistogram("ESD_TrueConversionMapping_ITSTPC_ZPhi",fV0Reader->GetZ() ,vtxConv.Phi());
				}
				fHistograms->FillHistogram("ESD_TrueConversion_OpeningAngle", fV0Reader->GetOpeningAngle());

				//----Histos for HFE-------------------------------------- 

				if(statusPos & AliESDtrack::kTOFpid){
					fHistograms->FillHistogram("ESD_TrueConversion_SinglePos_R", positiveMC->R(),fV0Reader->GetPositiveMCParticle()->Pt());
					if( TESTBIT(itsPixelPos, 0) ){ 
						fHistograms->FillHistogram("ESD_TrueConversion_SinglePos_kFirst_R", positiveMC->R(),fV0Reader->GetPositiveMCParticle()->Pt());
					}
				}
				if(statusNeg & AliESDtrack::kTOFpid){
					fHistograms->FillHistogram("ESD_TrueConversion_SingleNeg_R", negativeMC->R(),fV0Reader->GetNegativeMCParticle()->Pt());
					if( TESTBIT(itsPixelNeg, 0) ){ 
						fHistograms->FillHistogram("ESD_TrueConversion_SingleNeg_kFirst_R", negativeMC->R(),fV0Reader->GetNegativeMCParticle()->Pt());
					}
				}
				//--------------------------------------------------------

				fHistograms->FillHistogram("ESD_TrueConvGamma_CosPointingAngle", fV0Reader->GetCosPointingAngle());
				fHistograms->FillHistogram("ESD_TrueConvGamma_DcaDaughters", fV0Reader->GetDcaDaughters());
				fHistograms->FillHistogram("ESD_TrueConvGamma_NormDcaDistDaughters", fV0Reader->GetNormDcaDistDaughters());
				fHistograms->FillHistogram("ESD_TrueConvGamma_LikelihoodAP", fV0Reader->GetLikelihoodAP());
				if (fV0Reader->GetMotherCandidateP() != 0) {
					fHistograms->FillHistogram("ESD_TrueConvGamma_E_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetNegativeTrackP()/fV0Reader->GetMotherCandidateP());
					fHistograms->FillHistogram("ESD_TrueConvGamma_P_AsymmetryP",fV0Reader->GetMotherCandidateP(),fV0Reader->GetPositiveTrackP()/fV0Reader->GetMotherCandidateP());
				} else { 
					cout << "Error::fV0Reader->GetNegativeTrackP() == 0 !!!" << endl;		
				}
				
				fHistograms->FillHistogram("ESD_TrueConvGamma_E_dEdxP",fV0Reader->GetNegativeTrackP(),fV0Reader->GetNegativeTrackTPCdEdx());
				fHistograms->FillHistogram("ESD_TrueConvGamma_P_dEdxP",fV0Reader->GetPositiveTrackP(),fV0Reader->GetPositiveTrackTPCdEdx());
					
				//store MCTruth properties
				fHistograms->FillHistogram("ESD_TrueConvGamma_MC_Pt_Eta", fV0Reader->GetMotherMCParticle()->Pt(),fV0Reader->GetMotherMCParticle()->Eta());
				fHistograms->FillHistogram("ESD_TrueConversion_MC_ZR", negativeMC->Vz(),negativeMC->R());
				fHistograms->FillHistogram("ESD_TrueConversion_MC_XY", negativeMC->Vx(),negativeMC->Vy());
				
				//resolution
				Double_t mcpt	 = fV0Reader->GetMotherMCParticle()->Pt();
				Double_t esdpt	= fV0Reader->GetMotherCandidatePt();
				Double_t resdPt = 0.;
				if(mcpt > 0){
					resdPt = ((esdpt - mcpt)/mcpt)*100.;
				} else if(mcpt < 0){
					cout<<"Pt of MC particle is negative, this will cause wrong calculation of resPt"<<endl; 
				}
							
				fHistograms->FillHistogram("Resolution_Gamma_dPt_Pt", mcpt, resdPt);
				fHistograms->FillHistogram("Resolution_MCPt_ESDPt", mcpt,esdpt);
				fHistograms->FillHistogram("Resolution_Gamma_dPt_Phi", fV0Reader->GetMotherCandidatePhi(), resdPt);
				if (esdpt> 0.150){
					fHistograms->FillHistogram("Resolution_Gamma_minPt_dPt_Phi", fV0Reader->GetMotherCandidatePhi(), resdPt);
				}
						
				Double_t resdZ = 0.;
				if(fV0Reader->GetNegativeMCParticle()->Vz() != 0){
					resdZ = ((fV0Reader->GetZ() -fV0Reader->GetNegativeMCParticle()->Vz())/fV0Reader->GetNegativeMCParticle()->Vz())*100.;
				}
				Double_t resdZAbs = 0.;
				resdZAbs = (fV0Reader->GetZ() -fV0Reader->GetNegativeMCParticle()->Vz());

				fHistograms->FillHistogram("Resolution_dZAbs_VS_R", fV0Reader->GetNegativeMCParticle()->R(), resdZAbs);
				fHistograms->FillHistogram("Resolution_dZAbs_VS_Z", fV0Reader->GetNegativeMCParticle()->Vz(), resdZAbs);
				fHistograms->FillHistogram("Resolution_dZ", fV0Reader->GetNegativeMCParticle()->Vz(), resdZ);
				fHistograms->FillHistogram("Resolution_MCZ_ESDZ", fV0Reader->GetNegativeMCParticle()->Vz(),fV0Reader->GetZ());
					
				// new for dPt_Pt-histograms for Electron and Positron
				Double_t mcEpt;
				Double_t resEdPt=0;
				UInt_t kTRDoutN, statusN;
				Int_t nITSclsE;
				// new for dPt_Pt-histograms for Electron and Positron
				if (fV0Reader->GetNegativeMCParticle()->GetPdgCode() == 11) {
					mcEpt = fV0Reader->GetNegativeMCParticle()->Pt();
				} else {
					mcEpt = fV0Reader->GetPositiveMCParticle()->Pt();
				}
				if (mcEpt > 0){ 
					resEdPt = ((fV0Reader->GetNegativeTrackPt()-mcEpt)/mcEpt)*100.;
				}
				statusN = fV0Reader->GetNegativeESDTrack()->GetStatus(); 
				kTRDoutN =	(statusN & AliESDtrack::kTRDout);
				nITSclsE= fV0Reader->GetNegativeTracknITSClusters();

				// filling Resolution_Pt_dPt with respect to the Number of ITS clusters for Positrons
				switch(nITSclsE){
				case 0: // 0 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS0", mcEpt, resEdPt);
					break;
				case 1:	// 1 ITS cluster
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS1", mcEpt, resEdPt);
					break;
				case 2:	// 2 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS2", mcEpt, resEdPt);
					break;
				case 3: // 3 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS3", mcEpt, resEdPt);
					break;
				case 4: // 4 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS4", mcEpt, resEdPt);
					break;
				case 5: // 5 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS5", mcEpt, resEdPt);
					break;
				case 6: // 6 ITS clusters
					fHistograms->FillHistogram("Resolution_E_dPt_Pt_ITS6", mcEpt, resEdPt);
					break;
				}
				//Filling histograms with respect to Electron resolution
				fHistograms->FillHistogram("Resolution_E_dPt_Pt", mcEpt, resEdPt);
				fHistograms->FillHistogram("Resolution_E_dPt_Phi", fV0Reader->GetNegativeTrackPhi(), resEdPt);
				if (fV0Reader->GetNegativeTrackPt()> 0.150){
					fHistograms->FillHistogram("Resolution_E_minPt_dPt_Phi", fV0Reader->GetNegativeTrackPhi(), resEdPt);
				}

				if(kTRDoutN){
					fHistograms->FillHistogram("Resolution_E_nTRDtracklets_ESDPt", fV0Reader->GetNegativeTrackPt(), fV0Reader->GetNegativeESDTrack()->GetTRDntracklets());
					fHistograms->FillHistogram("Resolution_E_nTRDtracklets_MCPt", mcEpt, fV0Reader->GetNegativeESDTrack()->GetTRDntracklets());	
					fHistograms->FillHistogram("Resolution_E_nTRDclusters_ESDPt",fV0Reader->GetNegativeTrackPt(), fV0Reader->GetNegativeESDTrack()->GetTRDncls());
					fHistograms->FillHistogram("Resolution_E_nTRDclusters_MCPt",mcEpt, fV0Reader->GetNegativeESDTrack()->GetTRDncls());
					fHistograms->FillHistogram("Resolution_E_TRDsignal_ESDPt", fV0Reader->GetNegativeTrackPt(), fV0Reader->GetNegativeESDTrack()->GetTRDsignal());
				}

				Double_t mcPpt;
				if (fV0Reader->GetPositiveMCParticle()->GetPdgCode() == -11) {
					mcPpt = fV0Reader->GetPositiveMCParticle()->Pt();
				} else {
					mcPpt = fV0Reader->GetNegativeMCParticle()->Pt();
				}

				Double_t resPdPt = 0;
				if (mcPpt > 0){ 
					resPdPt = ((fV0Reader->GetPositiveTrackPt()-mcPpt)/mcPpt)*100.;
				}

				UInt_t statusP = fV0Reader->GetPositiveESDTrack()->GetStatus(); 
				//		 AliESDtrack * posTr= fV0Reader->GetPositiveESDTrack();
				UInt_t kTRDoutP =	(statusP & AliESDtrack::kTRDout);
				
				Int_t nITSclsP = fV0Reader->GetPositiveTracknITSClusters();
				// filling Resolution_Pt_dPt with respect to the Number of ITS clusters for Positrons
				switch(nITSclsP){
				case 0: // 0 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS0", mcPpt, resPdPt);
					break;
				case 1:	// 1 ITS cluster
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS1", mcPpt, resPdPt);
					break;
				case 2:	// 2 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS2", mcPpt, resPdPt);
					break;
				case 3: // 3 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS3", mcPpt, resPdPt);
					break;
				case 4: // 4 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS4", mcPpt, resPdPt);
					break;
				case 5: // 5 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS5", mcPpt, resPdPt);
					break;
				case 6: // 6 ITS clusters
					fHistograms->FillHistogram("Resolution_P_dPt_Pt_ITS6", mcPpt, resPdPt);
					break;
				}
				//Filling histograms with respect to Positron resolution
				fHistograms->FillHistogram("Resolution_P_dPt_Pt", mcPpt, resPdPt);
				fHistograms->FillHistogram("Resolution_P_dPt_Phi", fV0Reader->GetPositiveTrackPhi(), resPdPt);
				if (fV0Reader->GetPositiveTrackPt()> 0.150){
					fHistograms->FillHistogram("Resolution_P_minPt_dPt_Phi", fV0Reader->GetPositiveTrackPhi(), resPdPt);
				}

				if(kTRDoutP){
					fHistograms->FillHistogram("Resolution_P_nTRDtracklets_ESDPt", fV0Reader->GetPositiveTrackPt(), fV0Reader->GetPositiveESDTrack()->GetTRDntracklets());
					fHistograms->FillHistogram("Resolution_P_nTRDtracklets_MCPt", mcPpt, fV0Reader->GetPositiveESDTrack()->GetTRDntracklets());
					fHistograms->FillHistogram("Resolution_P_nTRDclusters_ESDPt",fV0Reader->GetPositiveTrackPt(), fV0Reader->GetPositiveESDTrack()->GetTRDncls());
					fHistograms->FillHistogram("Resolution_P_nTRDclusters_MCPt",mcPpt, fV0Reader->GetPositiveESDTrack()->GetTRDncls());
					fHistograms->FillHistogram("Resolution_P_TRDsignal_ESDPt", fV0Reader->GetPositiveTrackPt(), fV0Reader->GetPositiveESDTrack()->GetTRDsignal());
				}

			
				Double_t resdR = 0.;
				if(fV0Reader->GetNegativeMCParticle()->R() != 0){
					resdR = ((fV0Reader->GetXYRadius() - fV0Reader->GetNegativeMCParticle()->R())/fV0Reader->GetNegativeMCParticle()->R())*100.;
				}
				Double_t resdRAbs = 0.;
				resdRAbs = (fV0Reader->GetXYRadius() - fV0Reader->GetNegativeMCParticle()->R());

				fHistograms->FillHistogram("Resolution_dRAbs_VS_R", fV0Reader->GetNegativeMCParticle()->R(), resdRAbs);
				fHistograms->FillHistogram("Resolution_dR", fV0Reader->GetNegativeMCParticle()->R(), resdR);
				fHistograms->FillHistogram("Resolution_MCR_ESDR", fV0Reader->GetNegativeMCParticle()->R(),fV0Reader->GetXYRadius());
				fHistograms->FillHistogram("Resolution_R_dPt", fV0Reader->GetNegativeMCParticle()->R(), resdPt);
				if (esdpt> 0.150){
					fHistograms->FillHistogram("Resolution_minPt_R_dPt", fV0Reader->GetNegativeMCParticle()->R(), resdPt);
				}

				Double_t resdPhiAbs=0.;
				resdPhiAbs=0.;
				resdPhiAbs= (fV0Reader->GetMotherCandidatePhi()-fV0Reader->GetNegativeMCParticle()->Phi());
				fHistograms->FillHistogram("Resolution_MCPhi_ESDPhi", fV0Reader->GetNegativeMCParticle()->Phi(),fV0Reader->GetMotherCandidatePhi());
				fHistograms->FillHistogram("Resolution_dPhiAbs_VS_R", fV0Reader->GetNegativeMCParticle()->R(), resdPhiAbs);
				fHistograms->FillHistogram("Resolution_dPhiAbs_VS_Phi", fV0Reader->GetNegativeMCParticle()->Phi(), resdPhiAbs);
			}//if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22)
		}//if(fDoMCTruth)
	}//while(fV0Reader->NextV0)
	fHistograms->FillHistogram("ESD_NumberOfSurvivingV0s", nSurvivingV0s);
	fHistograms->FillHistogram("ESD_NumberOfV0s", fV0Reader->GetNumberOfV0s());
	fHistograms->FillHistogram("ESD_NumberOfContributorsVtx", fV0Reader->GetNumberOfContributorsVtx());
	fV0Reader->ResetV0IndexNumber();
}

 

//_____________________________________________________________________________________
void AliAnalysisTaskGammaConversion::ProcessGammasForOmegaMesonAnalysis(){
	// omega meson analysis pi0+gamma decay
	for(Int_t firstPi0Index=0;firstPi0Index<fKFReconstructedPi0sTClone->GetEntriesFast();firstPi0Index++){
		AliKFConversionMother * omegaCandidatePi0Daughter = (AliKFConversionMother *)fKFReconstructedPi0sTClone->At(firstPi0Index);
		for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){

			AliKFConversionPhoton * omegaCandidateGammaDaughter = (AliKFConversionPhoton *)fKFReconstructedGammasTClone->At(firstGammaIndex);
                           if(omegaCandidatePi0Daughter->GetGammaLabel(0)==firstGammaIndex || omegaCandidatePi0Daughter->GetGammaLabel(1)==firstGammaIndex){
				continue;
			}

			AliKFParticle omegaCandidate(*omegaCandidatePi0Daughter,*omegaCandidateGammaDaughter);
			Double_t massOmegaCandidate = 0.;
			Double_t widthOmegaCandidate = 0.;

			omegaCandidate.GetMass(massOmegaCandidate,widthOmegaCandidate);

			if ( massOmegaCandidate > 733 && massOmegaCandidate < 833 ) {
				//AddOmegaToAOD(&omegaCandidate, massOmegaCandidate, firstPi0Index, firstGammaIndex);
			}
			
			fHistograms->FillHistogram("ESD_Omega_InvMass_vs_Pt",massOmegaCandidate ,omegaCandidate.GetPt());
			fHistograms->FillHistogram("ESD_Omega_InvMass",massOmegaCandidate);
 
			//delete omegaCandidate;

		}// end of omega reconstruction in pi0+gamma channel

		if(fDoJet == kTRUE){
			AliKFParticle* negPiKF=NULL;
			AliKFParticle* posPiKF=NULL;
			Int_t piPlusMotherLabel=-1;
			Int_t piMinusMotherLabel=-1;
			
			// look at the pi+pi+pi0 channel 
			for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
				AliESDtrack* posTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
				if (posTrack->GetSign()<0) continue;
				if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kPion))>2.) continue;
				if (posPiKF) delete posPiKF; posPiKF=NULL;
				posPiKF = new AliKFParticle( *(posTrack) ,211);
				if(fDoMCTruth){
					TParticle * positiveMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[iCh])->GetLabel()));
					if(positiveMCParticle->GetMother(0)>-1){
						piPlusMotherLabel = positiveMCParticle->GetMother(0);
					}
				}
				
				for(Int_t jCh=0;jCh<fChargedParticles->GetEntriesFast();jCh++){
					AliESDtrack* negTrack = (AliESDtrack*)(fChargedParticles->At(jCh));
					if( negTrack->GetSign()>0) continue;
					if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kPion))>2.) continue;
					if (negPiKF) delete negPiKF; negPiKF=NULL;
					negPiKF = new AliKFParticle( *(negTrack) ,-211);
					AliKFParticle omegaCandidatePipPinPi0(*omegaCandidatePi0Daughter,*posPiKF,*negPiKF);
					Double_t massOmegaCandidatePipPinPi0 = 0.;
					Double_t widthOmegaCandidatePipPinPi0 = 0.;
					
					omegaCandidatePipPinPi0.GetMass(massOmegaCandidatePipPinPi0,widthOmegaCandidatePipPinPi0);

					if ( massOmegaCandidatePipPinPi0 > 733 && massOmegaCandidatePipPinPi0 < 833 ) {
						// AddOmegaToAOD(&omegaCandidatePipPinPi0, massOmegaCandidatePipPinPi0, -1, -1);
					}
					
					fHistograms->FillHistogram("ESD_OmegaPipPinPi0_InvMass_vs_Pt",massOmegaCandidatePipPinPi0 ,omegaCandidatePipPinPi0.GetPt());
					fHistograms->FillHistogram("ESD_OmegaPipPinPi0_InvMass",massOmegaCandidatePipPinPi0);


					if(fDoMCTruth){
						TParticle * negativeMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[jCh])->GetLabel()));
						if(negativeMCParticle->GetMother(0)>-1){
							piMinusMotherLabel = negativeMCParticle->GetMother(0);
							if(	piMinusMotherLabel == 	piPlusMotherLabel){
								Int_t geantCode=fStack->Particle(TMath::Abs(piPlusMotherLabel))->GetPdgCode();
								if(geantCode == 221 || geantCode == 223){
									fHistograms->FillHistogram("ESD_TrueOmegaPipPinPi0_InvMass_vs_Pt",massOmegaCandidatePipPinPi0 ,omegaCandidatePipPinPi0.GetPt());
								}
							}				

						}

					}
	
					//	delete omegaCandidatePipPinPi0;
				}
			}

			if (posPiKF) delete posPiKF; posPiKF=NULL;		 if (negPiKF) delete negPiKF; negPiKF=NULL;

		} // checking ig gammajet because in that case the chargedparticle list is created

	}
	////////////////////////

	// gamma+ pi+
	if(fDoJet == kTRUE){
		for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){
			AliKFConversionPhoton * rhoCandidateGammaDaughter = (AliKFConversionPhoton *)fKFReconstructedGammasTClone->At(firstGammaIndex);
			Int_t gamma1MotherLabel=-1;
			if(fDoMCTruth){
				Int_t indexKF1 = rhoCandidateGammaDaughter->GetV0Index();
				if(indexKF1<fV0Reader->GetNumberOfV0s()){
					fV0Reader->GetV0(indexKF1);//updates to the correct v0

					if(fV0Reader->HasSameMCMother() == kTRUE){
						//cout<<"This v0 is a real v0!!!!"<<endl;
						TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
						TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
						if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){
							if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){
								if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
									gamma1MotherLabel=fV0Reader->GetMotherMCParticle()->GetFirstMother();
								}
							}
						}
					}
				}
			}

			AliKFParticle* posPiKF=NULL;
			AliKFParticle* negPiKF=NULL;
			Int_t piPlusMotherLabel=-1;
			Int_t piMinusMotherLabel=-1;
			
 			for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
				AliESDtrack* posTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
				if (posTrack->GetSign()<0) continue;
				if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kPion))>2.) continue;
				if (posPiKF) delete posPiKF; posPiKF=NULL;
				posPiKF = new AliKFParticle( *(posTrack) ,211);
				AliKFParticle rhoPlusCandidate(*posPiKF,*rhoCandidateGammaDaughter);
				Double_t massRhoPlusCandidate = 0.;
				Double_t widthRhoPlusCandidate = 0.;

				rhoPlusCandidate.GetMass(massRhoPlusCandidate,widthRhoPlusCandidate);
				fHistograms->FillHistogram("ESD_RhoPlus_InvMass_vs_Pt",massRhoPlusCandidate ,rhoPlusCandidate.GetPt());
				fHistograms->FillHistogram("ESD_RhoPlus_InvMass",massRhoPlusCandidate);


				if(fDoMCTruth){
					TParticle * positiveMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[iCh])->GetLabel()));
					if(positiveMCParticle->GetMother(0)>-1){
						piPlusMotherLabel = positiveMCParticle->GetMother(0);
						if(piPlusMotherLabel == gamma1MotherLabel){
							//Int_t geantCode=fStack->Particle(TMath::Abs(pionMotherLabel))->GetPdgCode();
              //cout<<"RhoPlus::" << geantCode<< endl;
							fHistograms->FillHistogram("ESD_TrueRhoPlus_InvMass_vs_Pt",massRhoPlusCandidate ,rhoPlusCandidate.GetPt());
						}
					}
				}
				for(Int_t jCh=0;jCh<fChargedParticles->GetEntriesFast();jCh++){
					AliESDtrack* negTrack = (AliESDtrack*)(fChargedParticles->At(jCh));
					if (negTrack->GetSign()>0) continue;
					if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kPion))>2.) continue;
					if (negPiKF) delete negPiKF; negPiKF=NULL;
					negPiKF = new AliKFParticle( *(negTrack) ,-211);
					AliKFParticle rho0Candidate(*posPiKF,*negPiKF,*rhoCandidateGammaDaughter);
					Double_t massRho0Candidate = 0.;
					Double_t widthRho0Candidate = 0.;
					
					rho0Candidate.GetMass(massRho0Candidate,widthRho0Candidate);
					fHistograms->FillHistogram("ESD_Rho0_InvMass_vs_Pt",massRho0Candidate ,rho0Candidate.GetPt());
					fHistograms->FillHistogram("ESD_Rho0_InvMass",massRho0Candidate);

					if(fDoMCTruth){
						TParticle * negativeMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[jCh])->GetLabel()));
						if(negativeMCParticle->GetMother(0)>-1){
							piMinusMotherLabel = negativeMCParticle->GetMother(0);
							if(	piMinusMotherLabel == 	piPlusMotherLabel && piMinusMotherLabel==gamma1MotherLabel ){
								Int_t geantCode=fStack->Particle(TMath::Abs(piPlusMotherLabel))->GetPdgCode();
								if(geantCode == 221 || geantCode == 113){
									fHistograms->FillHistogram("ESD_TrueRho0_InvMass_vs_Pt",massRho0Candidate ,rho0Candidate.GetPt());
								}
							}				
						}
					}
				}

			}

			if (posPiKF) delete posPiKF; posPiKF=NULL;		 if (negPiKF) delete negPiKF; negPiKF=NULL;


 			for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
				AliESDtrack* negTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
				if (negTrack->GetSign()>0) continue;
				if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kPion))>2.) continue;
				if (negPiKF) delete negPiKF; negPiKF=NULL;
				negPiKF = new AliKFParticle( *(negTrack) ,-211);
				AliKFParticle rhoMinusCandidate(*negPiKF,*rhoCandidateGammaDaughter);
				Double_t massRhoMinusCandidate = 0.;
				Double_t widthRhoMinusCandidate = 0.;

				rhoMinusCandidate.GetMass(massRhoMinusCandidate,widthRhoMinusCandidate);
				fHistograms->FillHistogram("ESD_RhoMinus_InvMass_vs_Pt",massRhoMinusCandidate ,rhoMinusCandidate.GetPt());
				fHistograms->FillHistogram("ESD_RhoMinus_InvMass",massRhoMinusCandidate);

				if(fDoMCTruth){
					TParticle * negativeMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[iCh])->GetLabel()));
					Int_t pionMotherLabel=-1;

					if(negativeMCParticle->GetMother(0)>-1){
						pionMotherLabel = negativeMCParticle->GetMother(0);
						if(pionMotherLabel == gamma1MotherLabel){
							//Int_t geantCode=fStack->Particle(TMath::Abs(pionMotherLabel))->GetPdgCode();
              //cout<<"RhoMinus::" << geantCode<< endl;
							fHistograms->FillHistogram("ESD_TrueRhoMinus_InvMass_vs_Pt",massRhoMinusCandidate ,rhoMinusCandidate.GetPt());
						}
					}
				}
			}
			if (posPiKF) delete posPiKF; posPiKF=NULL;		 if (negPiKF) delete negPiKF; negPiKF=NULL;

			AliKFParticle* posProtKF=NULL;

			
 			for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
				AliESDtrack* posTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
				if (posTrack->GetSign()<0) continue;
				if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kProton))>2.) continue;
				if (posProtKF) delete posProtKF; posProtKF=NULL;
				posProtKF = new AliKFParticle( *(posTrack) ,2212);
				AliKFParticle deltaPlusCandidate(*posProtKF,*rhoCandidateGammaDaughter);
				Double_t massDeltaPlusCandidate = 0.;
				Double_t widthDeltaPlusCandidate = 0.;

				deltaPlusCandidate.GetMass(massDeltaPlusCandidate,widthDeltaPlusCandidate);
				fHistograms->FillHistogram("ESD_DeltaPlus_InvMass_vs_Pt",massDeltaPlusCandidate ,deltaPlusCandidate.GetPt());
				fHistograms->FillHistogram("ESD_DeltaPlus_InvMass",massDeltaPlusCandidate);

				if(fDoMCTruth){
					TParticle * positiveMCParticle = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(fChargedParticlesId[iCh])->GetLabel()));
					Int_t protonMotherLabel=-1;

					if(positiveMCParticle->GetMother(0)>-1){
						protonMotherLabel = positiveMCParticle->GetMother(0);
						if(protonMotherLabel == gamma1MotherLabel){
							//Int_t geantCode=fStack->Particle(TMath::Abs(pionMotherLabel))->GetPdgCode();
              //cout<<"RhoPlus::" << geantCode<< endl;
							fHistograms->FillHistogram("ESD_TrueDeltaPlus_InvMass_vs_Pt",massDeltaPlusCandidate ,deltaPlusCandidate.GetPt());
						}
					}
				}

			}
			if (posPiKF) delete posPiKF; posPiKF=NULL;		 if (negPiKF) delete negPiKF; negPiKF=NULL; if (posProtKF) delete posProtKF; posProtKF=NULL;
		}
	}
		
	if(fCalculateBackground){

		AliGammaConversionBGHandler * bgHandler = fV0Reader->GetBGHandler();
		
		Int_t zbin= bgHandler->GetZBinIndex(fV0Reader->GetVertexZ());
		Int_t mbin = 0;
		if(fUseTrackMultiplicityForBG == kTRUE){
			mbin = bgHandler->GetMultiplicityBinIndex(fV0Reader->CountESDTracks());
		} else {
			mbin = bgHandler->GetMultiplicityBinIndex(fV0Reader->GetNGoodV0s());
		}
		
		AliGammaConversionBGHandler::GammaConversionVertex *bgEventVertex = NULL;

		// Background calculation for the omega
		for(Int_t nEventsInBG=0;nEventsInBG <fV0Reader->GetNBGEvents();nEventsInBG++){
			AliGammaConversionKFVector * previousEventV0s = bgHandler->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			
			if(fMoveParticleAccordingToVertex == kTRUE){
				bgEventVertex = bgHandler->GetBGEventVertex(zbin,mbin,nEventsInBG);
			}
			for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
				AliKFParticle previousGoodV0 = (AliKFParticle)(*(previousEventV0s->at(iPrevious)));

				if(fMoveParticleAccordingToVertex == kTRUE){
					MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
				}

				for(Int_t firstPi0Index=0;firstPi0Index<fKFReconstructedPi0sTClone->GetEntriesFast();firstPi0Index++){
					AliKFParticle * omegaCandidatePi0Daughter = (AliKFParticle *)fKFReconstructedPi0sTClone->At(firstPi0Index);
					AliKFParticle * omegaBckCandidate = new AliKFParticle(*omegaCandidatePi0Daughter,previousGoodV0);
					Double_t massOmegaBckCandidate = 0.;
					Double_t widthOmegaBckCandidate = 0.;
					
					omegaBckCandidate->GetMass(massOmegaBckCandidate,widthOmegaBckCandidate);


					fHistograms->FillHistogram("ESD_Omega_Bck_InvMass_vs_Pt",massOmegaBckCandidate ,omegaBckCandidate->GetPt());
					fHistograms->FillHistogram("ESD_Omega_Bck_InvMass",massOmegaBckCandidate);

					delete omegaBckCandidate; 
				}

				// Bck for gamma pi+ pi-

				AliKFParticle* posPiKF=NULL;
				AliKFParticle* negPiKF=NULL;
			
				for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
					AliESDtrack* posTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
					if (posTrack->GetSign()<0) continue;
					if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kPion))>2.) continue;
					if (posPiKF) delete posPiKF; posPiKF=NULL;
					posPiKF = new AliKFParticle( *(posTrack) ,211);

					for(Int_t jCh=0;jCh<fChargedParticles->GetEntriesFast();jCh++){
						AliESDtrack* negTrack = (AliESDtrack*)(fChargedParticles->At(jCh));
						if (negTrack->GetSign()>0) continue;
						if(TMath::Abs(fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kPion))>2.) continue;
						if (negPiKF) delete negPiKF; negPiKF=NULL;
						negPiKF = new AliKFParticle( *(negTrack) ,-211);
						AliKFParticle rho0BckCandidate(*posPiKF,*negPiKF,previousGoodV0);
						Double_t massRho0BckCandidate = 0.;
						Double_t widthRho0BckCandidate = 0.;
						
						rho0BckCandidate.GetMass(massRho0BckCandidate,widthRho0BckCandidate);
						fHistograms->FillHistogram("ESD_Rho0Bck_InvMass_vs_Pt",massRho0BckCandidate ,rho0BckCandidate.GetPt());
						fHistograms->FillHistogram("ESD_Rho0Bck_InvMass",massRho0BckCandidate);
					}
					
				}

				if (posPiKF) delete posPiKF; posPiKF=NULL;		 if (negPiKF) delete negPiKF; negPiKF=NULL; 


			}
		}
	} // end of checking if background calculation is available
}


void AliAnalysisTaskGammaConversion::ProcessGammasForNeutralMesonAnalysis(){
	// see header file for documentation
	
	//	for(UInt_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammas.size();firstGammaIndex++){
	//		for(UInt_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFReconstructedGammas.size();secondGammaIndex++){

	fESDEvent = fV0Reader->GetESDEvent();

	if(fKFReconstructedGammasTClone->GetEntriesFast()>fV0Reader->GetNumberOfV0s()){
		cout<<"Warning, number of entries in the tclone is bigger than number of v0s"<<endl;
	}

	for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){
		for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();secondGammaIndex++){
			
			//			AliKFParticle * twoGammaDecayCandidateDaughter0 = &fKFReconstructedGammas[firstGammaIndex];
			//			AliKFParticle * twoGammaDecayCandidateDaughter1 = &fKFReconstructedGammas[secondGammaIndex];
			
			AliKFConversionPhoton * twoGammaDecayCandidateDaughter0 = (AliKFConversionPhoton *)fKFReconstructedGammasTClone->At(firstGammaIndex);
			AliKFConversionPhoton * twoGammaDecayCandidateDaughter1 = (AliKFConversionPhoton *)fKFReconstructedGammasTClone->At(secondGammaIndex);

                           if(twoGammaDecayCandidateDaughter0->GetTrackLabelPositive()==twoGammaDecayCandidateDaughter1->GetTrackLabelPositive()||
			   twoGammaDecayCandidateDaughter0->GetTrackLabelPositive()==twoGammaDecayCandidateDaughter1->GetTrackLabelNegative()||
			   twoGammaDecayCandidateDaughter0->GetTrackLabelNegative()==twoGammaDecayCandidateDaughter1->GetTrackLabelPositive()||
			   twoGammaDecayCandidateDaughter0->GetTrackLabelNegative()==twoGammaDecayCandidateDaughter1->GetTrackLabelNegative())continue;
			
			AliKFConversionMother *twoGammaCandidate = new AliKFConversionMother(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);
			twoGammaCandidate->SetGammaLabels(firstGammaIndex,secondGammaIndex);

			Double_t massTwoGammaCandidate = 0.;
			Double_t widthTwoGammaCandidate = 0.;
			Double_t chi2TwoGammaCandidate =10000.;	
			twoGammaCandidate->GetMass(massTwoGammaCandidate,widthTwoGammaCandidate);
			//			if(twoGammaCandidate->GetNDF()>0){
			//	chi2TwoGammaCandidate = twoGammaCandidate->GetChi2()/twoGammaCandidate->GetNDF();
			chi2TwoGammaCandidate = twoGammaCandidate->GetChi2();
				
			fHistograms->FillHistogram("ESD_Mother_Chi2",chi2TwoGammaCandidate);
			if((chi2TwoGammaCandidate>0 && chi2TwoGammaCandidate<fV0Reader->GetChi2CutMeson()) || fApplyChi2Cut == kFALSE){
					
				TVector3 momentumVectorTwoGammaCandidate(twoGammaCandidate->GetPx(),twoGammaCandidate->GetPy(),twoGammaCandidate->GetPz());
				TVector3 spaceVectorTwoGammaCandidate(twoGammaCandidate->GetX(),twoGammaCandidate->GetY(),twoGammaCandidate->GetZ());
								
                                  Double_t openingAngleTwoGammaCandidate = twoGammaCandidate->GetOpeningAngle();

				Double_t rapidity=twoGammaCandidate->GetRapidity();
			      
				if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut()){
					delete twoGammaCandidate;
					continue;	 // rapidity cut
				}

			         Double_t alfa=twoGammaCandidate->GetAlpha();

								
				if(openingAngleTwoGammaCandidate < fMinOpeningAngleGhostCut){
					delete twoGammaCandidate;
					continue;	 // minimum opening angle to avoid using ghosttracks
				}
			
				if(alfa>fV0Reader->GetAlphaMinCutMeson() && alfa<fV0Reader->GetAlphaCutMeson()){
					fHistograms->FillHistogram("ESD_Mother_GammaDaughter_OpeningAngle", openingAngleTwoGammaCandidate);
					fHistograms->FillHistogram("ESD_Mother_Energy", twoGammaCandidate->GetE());
					fHistograms->FillHistogram("ESD_Mother_Pt", momentumVectorTwoGammaCandidate.Pt());
					fHistograms->FillHistogram("ESD_Mother_Eta", momentumVectorTwoGammaCandidate.Eta());
					fHistograms->FillHistogram("ESD_Mother_Rapid", rapidity);					
					fHistograms->FillHistogram("ESD_Mother_Phi", spaceVectorTwoGammaCandidate.Phi());
					fHistograms->FillHistogram("ESD_Mother_Mass", massTwoGammaCandidate);
					fHistograms->FillHistogram("ESD_Mother_alfa", alfa);
					if( (massTwoGammaCandidate > 0.1) && (massTwoGammaCandidate < 0.15) ){
						fHistograms->FillHistogram("ESD_Mother_alfa_Pi0", alfa);
						fHistograms->FillHistogram("ESD_Mother_Pt_alpha_Pi0", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
					}
					if( (massTwoGammaCandidate > 0.5) && (massTwoGammaCandidate < 0.57) ){
						fHistograms->FillHistogram("ESD_Mother_alfa_Eta", alfa);
						fHistograms->FillHistogram("ESD_Mother_Pt_alpha_Eta", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
					}

					fHistograms->FillHistogram("ESD_Mother_R", spaceVectorTwoGammaCandidate.Pt());		// Pt in Space == R!!!
					fHistograms->FillHistogram("ESD_Mother_ZR", twoGammaCandidate->GetZ(), spaceVectorTwoGammaCandidate.Pt());
					fHistograms->FillHistogram("ESD_Mother_XY", twoGammaCandidate->GetX(), twoGammaCandidate->GetY());
					fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
					fHistograms->FillHistogram("ESD_Mother_InvMass",massTwoGammaCandidate);			
					fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt_alpha",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
				}
				if(alfa<0.1){
					fHistograms->FillHistogram("ESD_Mother_InvMass_vs_E_alpha",massTwoGammaCandidate ,twoGammaCandidate->GetE());
				}

				if(fCalculateBackground){
					/* Kenneth, just for testing*/
					AliGammaConversionBGHandler * bgHandlerTest = fV0Reader->GetBGHandler();
						
					Int_t zbin= bgHandlerTest->GetZBinIndex(fV0Reader->GetVertexZ());
					Int_t mbin=0;
					Int_t multKAA=0;
					if(fUseTrackMultiplicityForBG == kTRUE){
						multKAA=fV0Reader->CountESDTracks();
						mbin = bgHandlerTest->GetMultiplicityBinIndex(fV0Reader->CountESDTracks());
					}
					else{// means we use #v0s for multiplicity
						multKAA=fV0Reader->GetNGoodV0s();
						mbin = bgHandlerTest->GetMultiplicityBinIndex(fV0Reader->GetNGoodV0s());
					}
					//			cout<<"Filling bin number "<<zbin<<" and "<<mbin<<endl;
					//			cout<<"Corresponding to z = "<<fV0Reader->GetVertexZ()<<" and m = "<<multKAA<<endl;
					if(alfa>fV0Reader->GetAlphaMinCutMeson() && alfa<fV0Reader->GetAlphaCutMeson()){
						fHistograms->FillHistogram(Form("%d%dESD_Mother_InvMass",zbin,mbin),massTwoGammaCandidate);
						fHistograms->FillHistogram(Form("%d%dESD_Mother_InvMass_vs_Pt",zbin,mbin),massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
						/* end Kenneth, just for testing*/
						fHistograms->FillHistogram(Form("%dESD_Mother_InvMass_vs_Pt",mbin),massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
					}
				}
				/*		if(fCalculateBackground){
						AliGammaConversionBGHandler * bgHandler = fV0Reader->GetBGHandler();
						Int_t mbin= bgHandler->GetMultiplicityBinIndex(fV0Reader->CountESDTracks());
						fHistograms->FillHistogram(Form("%dESD_Mother_InvMass_vs_Pt",mbin),massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
						}*/
				//		if(fDoNeutralMesonV0MCCheck){
				if(fDoMCTruth){
					//Kenneth: Checking the eta of the gamma to check the difference between 0.9 and 1.2
				        Int_t indexKF1 = twoGammaDecayCandidateDaughter0->GetV0Index();
					if(indexKF1<fV0Reader->GetNumberOfV0s()){
						fV0Reader->GetV0(indexKF1);//updates to the correct v0
						Double_t eta1 = fV0Reader->GetMotherCandidateEta();
						Bool_t isRealPi0=kFALSE;
						Bool_t isRealEta=kFALSE;
						Int_t gamma1MotherLabel=-1;
						if(fV0Reader->HasSameMCMother() == kTRUE){
							//cout<<"This v0 is a real v0!!!!"<<endl;
							TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
							TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
							if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){
								if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){
									if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
										gamma1MotherLabel=fV0Reader->GetMotherMCParticle()->GetFirstMother();
									}
								}
								if(fV0Reader->GetMotherMCParticle()->GetPdgCode() ==111){
									gamma1MotherLabel=-111;
								}
								if(fV0Reader->GetMotherMCParticle()->GetPdgCode() ==221){
									gamma1MotherLabel=-221;
								}
							}
						}
                                                Int_t indexKF2 = twoGammaDecayCandidateDaughter1->GetV0Index();
						if(indexKF1 == indexKF2){
							cout<<"index of the two KF particles are the same.... should not happen"<<endl;
						}
						if(indexKF2<fV0Reader->GetNumberOfV0s()){
							fV0Reader->GetV0(indexKF2);
							Double_t eta2 = fV0Reader->GetMotherCandidateEta();
							Int_t gamma2MotherLabel=-1;
							if(fV0Reader->HasSameMCMother() == kTRUE){
							TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
							TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
							if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){
								if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){
									if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
										gamma2MotherLabel=fV0Reader->GetMotherMCParticle()->GetFirstMother();
									}
								}
								if(fV0Reader->GetMotherMCParticle()->GetPdgCode() ==111){
									gamma2MotherLabel=-111;
								}
								if(fV0Reader->GetMotherMCParticle()->GetPdgCode() ==221){
									gamma2MotherLabel=-221;
								}
									
								}
							}
							if(gamma1MotherLabel>=0 && gamma1MotherLabel==gamma2MotherLabel){
								if(fV0Reader->CheckIfPi0IsMother(gamma1MotherLabel)){
									isRealPi0=kTRUE;
								}
								if(fV0Reader->CheckIfEtaIsMother(gamma1MotherLabel)){
									isRealEta=kTRUE;
								}
							}

							//cout << "alpha   " << alfa << endl;
							if(isRealPi0)fHistograms->FillHistogram("ESD_TruePi0_alpha",alfa);
							if(isRealEta)fHistograms->FillHistogram("ESD_TrueEta_alpha",alfa);

							if(alfa>fV0Reader->GetAlphaMinCutMeson() && alfa<fV0Reader->GetAlphaCutMeson()){
								if(TMath::Abs(eta1)>0.9 && TMath::Abs(eta2)>0.9){
							//			fHistograms->FillHistogram("ESD_Mother_InvMass_1212",massTwoGammaCandidate);
							//			fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt1212",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
									if(isRealPi0 || isRealEta){
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_1212",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_OpeningAngle_1212",openingAngleTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt1212",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt_alpha",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										if( (isRealPi0) && (massTwoGammaCandidate > 0.1) && (massTwoGammaCandidate < 0.15) )
											fHistograms->FillHistogram("ESD_TruePi0_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
										if( (isRealEta) && (massTwoGammaCandidate > 0.5) && (massTwoGammaCandidate < 0.57) )
											fHistograms->FillHistogram("ESD_TrueEta_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha

									}

									if(!isRealPi0 && !isRealEta){
										if(gamma1MotherLabel>-1 && gamma2MotherLabel>-1){
											fHistograms->FillHistogram("ESD_TrueBckGG_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										} else {
											fHistograms->FillHistogram("ESD_TrueBckCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
										if(gamma1MotherLabel==-111 || gamma2MotherLabel==-111 || gamma1MotherLabel==-221 || gamma2MotherLabel==-221){
											fHistograms->FillHistogram("ESD_TruePi0DalitzCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
									}
								} else if(TMath::Abs(eta1)>0.9 || TMath::Abs(eta2)>0.9){
									//			fHistograms->FillHistogram("ESD_Mother_InvMass_0912",massTwoGammaCandidate);
									//			fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt0912",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
									
									if(isRealPi0 || isRealEta){
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_0912",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_OpeningAngle_0912",openingAngleTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt0912",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt_alpha",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										if( (isRealPi0) && (massTwoGammaCandidate > 0.1) && (massTwoGammaCandidate < 0.15) )
											fHistograms->FillHistogram("ESD_TruePi0_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
										if( (isRealEta) && (massTwoGammaCandidate > 0.5) && (massTwoGammaCandidate < 0.57) )
											fHistograms->FillHistogram("ESD_TrueEta_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
									}
									if(!isRealPi0 && !isRealEta){
										if(gamma1MotherLabel>-1 && gamma2MotherLabel>-1){
											fHistograms->FillHistogram("ESD_TrueBckGG_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}else{
											fHistograms->FillHistogram("ESD_TrueBckCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
										if(gamma1MotherLabel==-111 || gamma2MotherLabel==-111 || gamma1MotherLabel==-221 || gamma2MotherLabel==-221){
											fHistograms->FillHistogram("ESD_TruePi0DalitzCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
									}
								} else {
									//			fHistograms->FillHistogram("ESD_Mother_InvMass_0909",massTwoGammaCandidate);
									//			fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt0909",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
									if(isRealPi0 || isRealEta){
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_0909",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_OpeningAngle_0909",openingAngleTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt0909",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										fHistograms->FillHistogram("ESD_TruePi0_InvMass",massTwoGammaCandidate);
										fHistograms->FillHistogram("ESD_TruePi0_InvMass_vs_Pt_alpha",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
										if( (isRealPi0) && (massTwoGammaCandidate > 0.1) && (massTwoGammaCandidate < 0.15) )
											fHistograms->FillHistogram("ESD_TruePi0_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
										if( (isRealEta) && (massTwoGammaCandidate > 0.5) && (massTwoGammaCandidate < 0.57) )
											fHistograms->FillHistogram("ESD_TrueEta_Pt_alpha", momentumVectorTwoGammaCandidate.Pt(), alfa); //RR_alpha
										if(gamma1MotherLabel > fV0Reader->GetMCStack()->GetNprimary()){
											fHistograms->FillHistogram("ESD_TruePi0Sec_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
											TParticle * secPi0MC = (TParticle*)fStack->Particle(fV0Reader->GetMotherMCParticle()->GetFirstMother());
											if (secPi0MC->GetMother(0) >-1){
												if(fStack->Particle(secPi0MC->GetMother(0))->GetPdgCode()==kK0Short){
													fHistograms->FillHistogram("ESD_TruePi0SecFromK0S_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
													if(massTwoGammaCandidate>0.09 && massTwoGammaCandidate<0.145){
														fHistograms->FillHistogram("ESD_K0SFromSecPi0_Pt",fStack->Particle(secPi0MC->GetMother(0))->Pt());
													}
												}
											}
										}
									}
									if(!isRealPi0 && !isRealEta){
										if(gamma1MotherLabel>-1 && gamma2MotherLabel>-1){
											fHistograms->FillHistogram("ESD_TrueBckGG_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}else{
											fHistograms->FillHistogram("ESD_TrueBckCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
										if(gamma1MotherLabel==-111 || gamma2MotherLabel==-111 || gamma1MotherLabel==-221 || gamma2MotherLabel==-221 ){
											fHistograms->FillHistogram("ESD_TruePi0DalitzCont_InvMass_vs_Pt",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
										}
									}
								}
							}
						}
					}
				}
				if(alfa>fV0Reader->GetAlphaMinCutMeson() && alfa<fV0Reader->GetAlphaCutMeson()){
					if ( TMath::Abs(twoGammaDecayCandidateDaughter0->GetEta())<0.9 &&	TMath::Abs(twoGammaDecayCandidateDaughter1->GetEta())<0.9 ){
						fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt_Fiducial",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
						fHistograms->FillHistogram("ESD_Mother_InvMass_Fiducial",massTwoGammaCandidate);
					}
						
					if(TMath::Abs(twoGammaDecayCandidateDaughter0->GetEta())>0.9 && TMath::Abs(twoGammaDecayCandidateDaughter1->GetEta())>0.9){
						fHistograms->FillHistogram("ESD_Mother_InvMass_1212",massTwoGammaCandidate);
						fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt1212",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
					}
					else if(TMath::Abs(twoGammaDecayCandidateDaughter0->GetEta())>0.9 || TMath::Abs(twoGammaDecayCandidateDaughter1->GetEta())>0.9){
						fHistograms->FillHistogram("ESD_Mother_InvMass_0912",massTwoGammaCandidate);
						fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt0912",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
					}
					else{
						fHistograms->FillHistogram("ESD_Mother_InvMass_0909",massTwoGammaCandidate);
						fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt0909",massTwoGammaCandidate,momentumVectorTwoGammaCandidate.Pt());
					}

					Double_t lowMassPi0=0.1;
					Double_t highMassPi0=0.15;
					if ( ( massTwoGammaCandidate > lowMassPi0) && (massTwoGammaCandidate < highMassPi0) ){
						new((*fKFReconstructedPi0sTClone)[fKFReconstructedPi0sTClone->GetEntriesFast()])	AliKFConversionMother(*twoGammaCandidate);
					      }				
						
					if( fKFCreateAOD ) {
						lowMassPi0=0.08;
						highMassPi0=0.2;
						Double_t lowMassEta=0.4;
						Double_t highMassEta=0.7;

						if ( ( massTwoGammaCandidate > lowMassPi0) && (massTwoGammaCandidate < highMassPi0) ){

						    TagDaughter(twoGammaCandidate->GetGammaLabel(0));
						    TagDaughter(twoGammaCandidate->GetGammaLabel(1));
						//  AddPionToAOD(twoGammaCandidate);
						} else	if ( ( massTwoGammaCandidate > lowMassEta) && (massTwoGammaCandidate < highMassEta) ){
						    TagDaughter(twoGammaCandidate->GetGammaLabel(0));
						    TagDaughter(twoGammaCandidate->GetGammaLabel(1));
						    // AddPionToAOD(twoGammaCandidate);
						}
					} // if create aod

				}
			}
			delete twoGammaCandidate;
		}
	}
}

///__________________________________________________________________________________
void AliAnalysisTaskGammaConversion::AddGammaToAOD(AliKFConversionPhoton * kfParticle) {

    //Fill AOD with particles
    TClonesArray *branch=fAODGamma;
    if(branch){
		new((*branch)[branch->GetEntriesFast()])  AliAODConversionPhoton(kfParticle);
    } else {
		return;
	}

    //Add PID information with ESD tender (AOD implementation is not complete)

    AliAODConversionPhoton *gamma=static_cast<AliAODConversionPhoton*>(fAODGamma->At(fAODGamma->GetEntriesFast()-1));

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    AliPIDResponse *fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

    fESDEvent = fV0Reader->GetESDEvent();

    if(fESDEvent){
	Int_t labelp=((AliESDv0*)fESDEvent->GetV0(kfParticle->GetV0Index()))->GetPindex();
	Int_t labeln=((AliESDv0*)fESDEvent->GetV0(kfParticle->GetV0Index()))->GetNindex();

	AliESDtrack *trackpos=fESDEvent->GetTrack(labelp);
	AliESDtrack *trackneg=fESDEvent->GetTrack(labeln);

	if(trackpos&&trackneg&&fPIDResponse){

	    Float_t fNSigmadEdxPositive[5];
	    Float_t fNSigmadEdxNegative[5];

	    fNSigmadEdxPositive[0]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kElectron);
	    fNSigmadEdxPositive[1]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kMuon);
	    fNSigmadEdxPositive[2]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kPion);
	    fNSigmadEdxPositive[3]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kKaon);
	    fNSigmadEdxPositive[4]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kProton);

	    fNSigmadEdxNegative[0]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kElectron);
	    fNSigmadEdxNegative[1]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kMuon);
	    fNSigmadEdxNegative[2]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kPion);
	    fNSigmadEdxNegative[3]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kKaon);
	    fNSigmadEdxNegative[4]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kProton);

	    gamma->SetNSigmadEdx(fNSigmadEdxPositive,fNSigmadEdxNegative);
	}
    }
}

/*///__________________________________________________________________________________
void AliAnalysisTaskGammaConversion::AddPionToAOD(AliKFConversionMother * kfParticle) {

    //Add pions to AOD
    TClonesArray *branch=fAODPi0;

    if(branch){
	new((*branch)[branch->GetEntriesFast()])  AliAODConversionMother(kfParticle);
    }

    TagDaughter(kfParticle->GetGammaLabel(0));
    TagDaughter(kfParticle->GetGammaLabel(1));
}

///__________________________________________________________________________________
void AliAnalysisTaskGammaConversion::AddOmegaToAOD(AliKFParticle * kfParticle, Int_t daughter1, Int_t daughter2) {

    //Add omegas to AOD

    TClonesArray *branch=fAODOmega;
    // Get Daughters
   // AliAODConversionPhoton * fdaughter1 = dynamic_cast<AliAODConversionPhoton*>(fAODGamma->At(daughter1));
   // AliAODConversionPhoton * fdaughter2 = dynamic_cast<AliAODConversionPhoton*>(fAODGamma->At(daughter2));

    if(branch){
     //   new((*branch)[branch->GetEntriesFast()])  AliAODConversionMother(kfParticle);
		}

	TagDaughter(daughter1);
	TagDaughter(daughter2);


}
*/
///__________________________________________________________________________________
void AliAnalysisTaskGammaConversion::TagDaughter(Int_t gammaIndex) {
	//Set conversion tag on pion daughters
	AliAODConversionPhoton * daughter = dynamic_cast<AliAODConversionPhoton*>(fAODGamma->At(gammaIndex));
	if(daughter) {
		daughter->SetTag(kTRUE);
	} else {
		AliError("Daughter not in gamma tree!!");
	}
}

///___________________________________________________________________________________
void AliAnalysisTaskGammaConversion::FillAODWithConversionGammas(){
	// Fill AOD with reconstructed Gamma
	for(Int_t gammaIndex=0;gammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();gammaIndex++){
		AliKFConversionPhoton * gammakf = dynamic_cast<AliKFConversionPhoton*>(fKFReconstructedGammasTClone->At(gammaIndex));
		if(gammakf) {
                      AddGammaToAOD(gammakf);
		}
	}
}

/*
	void AliAnalysisTaskGammaConversion::ProcessConvPHOSGammasForNeutralMesonAnalysis(){

	// see header file for documentation
	// Analyse Pi0 with one photon from Phos and 1 photon from conversions
	


	Double_t vtx[3];
	vtx[0] = fV0Reader->GetPrimaryVertex()->GetX();
	vtx[1] = fV0Reader->GetPrimaryVertex()->GetY();
	vtx[2] = fV0Reader->GetPrimaryVertex()->GetZ();


	// Loop over all CaloClusters and consider only the PHOS ones:
	AliESDCaloCluster *clu;
	TLorentzVector pPHOS;
	TLorentzVector gammaPHOS;
	TLorentzVector gammaGammaConv;
	TLorentzVector pi0GammaConvPHOS;
	TLorentzVector gammaGammaConvBck;
	TLorentzVector pi0GammaConvPHOSBck;


	for (Int_t i=0; i<fV0Reader->GetESDEvent()->GetNumberOfCaloClusters(); i++) {
	clu = fV0Reader->GetESDEvent()->GetCaloCluster(i);
	if ( !clu->IsPHOS() || clu->E()<0.1 ) continue;
	clu ->GetMomentum(pPHOS ,vtx);
	for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){
	AliKFParticle * twoGammaDecayCandidateDaughter0 = (AliKFParticle *)fKFReconstructedGammasTClone->At(firstGammaIndex);
	gammaGammaConv.SetXYZM(twoGammaDecayCandidateDaughter0->Px(),twoGammaDecayCandidateDaughter0->Py(),twoGammaDecayCandidateDaughter0->Pz(),0.);
	gammaPHOS.SetXYZM(pPHOS.Px(),pPHOS.Py(),pPHOS.Pz(),0.);
	pi0GammaConvPHOS=gammaGammaConv+gammaPHOS;
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvPHOS",pi0GammaConvPHOS.M());
	fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvPHOS",pi0GammaConvPHOS.M(),pi0GammaConvPHOS.Pt());

	TVector3 v3D0(twoGammaDecayCandidateDaughter0->Px(),twoGammaDecayCandidateDaughter0->Py(),twoGammaDecayCandidateDaughter0->Pz());
	TVector3 v3D1(gammaPHOS.Px(),gammaPHOS.Py(),gammaPHOS.Pz());
	Double_t opanConvPHOS= v3D0.Angle(v3D1);
	if ( opanConvPHOS < 0.35){
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvPHOS_OpanLow",pi0GammaConvPHOS.M());
	}else{
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvPHOS_OpanHigh",pi0GammaConvPHOS.M());
	}

	}

	//	Now the LorentVector pPHOS is obtained and can be paired with the converted proton
	}
	//==== End of the PHOS cluster selection ============
	TLorentzVector pEMCAL;
	TLorentzVector gammaEMCAL;
	TLorentzVector pi0GammaConvEMCAL;
	TLorentzVector pi0GammaConvEMCALBck;

	for (Int_t i=0; i<fV0Reader->GetESDEvent()->GetNumberOfCaloClusters(); i++) {
	clu = fV0Reader->GetESDEvent()->GetCaloCluster(i);
	if ( !clu->IsEMCAL()	|| clu->E()<0.1 ) continue;
	if (clu->GetNCells() <= 1) continue;
	if ( clu->GetTOF()*1e9 < 550	|| clu->GetTOF()*1e9 > 750) continue;

	clu ->GetMomentum(pEMCAL ,vtx);
	for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){
	AliKFParticle * twoGammaDecayCandidateDaughter0 = (AliKFParticle *)fKFReconstructedGammasTClone->At(firstGammaIndex);
	gammaGammaConv.SetXYZM(twoGammaDecayCandidateDaughter0->Px(),
	twoGammaDecayCandidateDaughter0->Py(),
	twoGammaDecayCandidateDaughter0->Pz(),0.);
	gammaEMCAL.SetXYZM(pEMCAL.Px(),pEMCAL.Py(),pEMCAL.Pz(),0.);
	pi0GammaConvEMCAL=gammaGammaConv+gammaEMCAL;
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvEMCAL",pi0GammaConvEMCAL.M());
	fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvEMCAL",pi0GammaConvEMCAL.M(),pi0GammaConvEMCAL.Pt());
	TVector3 v3D0(twoGammaDecayCandidateDaughter0->Px(),
	twoGammaDecayCandidateDaughter0->Py(),
	twoGammaDecayCandidateDaughter0->Pz());
	TVector3 v3D1(gammaEMCAL.Px(),gammaEMCAL.Py(),gammaEMCAL.Pz());


	Double_t opanConvEMCAL= v3D0.Angle(v3D1);
	if ( opanConvEMCAL < 0.35){
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvEMCAL_OpanLow",pi0GammaConvEMCAL.M());
	}else{
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvEMCAL_OpanHigh",pi0GammaConvEMCAL.M());
	}

	}
	if(fCalculateBackground){
	for(Int_t nEventsInBG=0;nEventsInBG <fV0Reader->GetNBGEvents();nEventsInBG++){
	AliGammaConversionKFVector * previousEventV0s = fV0Reader->GetBGGoodV0s(nEventsInBG);
	for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
	AliKFParticle previousGoodV0 = (AliKFParticle)(*(previousEventV0s->at(iPrevious)));
	gammaGammaConvBck.SetXYZM(previousGoodV0.Px(),
	previousGoodV0.Py(),
	previousGoodV0.Pz(),0.);
	pi0GammaConvEMCALBck=gammaGammaConvBck+gammaEMCAL;
	fHistograms->FillHistogram("ESD_Mother_InvMass_GammaConvEMCAL_Bck",pi0GammaConvEMCALBck.M());
	fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt_GammaConvEMCAL_Bck",pi0GammaConvEMCALBck.M(),
	pi0GammaConvEMCALBck.Pt());
	}
	}
			
	//	Now the LorentVector pEMCAL is obtained and can be paired with the converted proton
	} // end of checking if background photons are available
	}
	//==== End of the PHOS cluster selection ============

	}
*/

void AliAnalysisTaskGammaConversion::MoveParticleAccordingToVertex(AliKFParticle * particle,const AliGammaConversionBGHandler::GammaConversionVertex *vertex){
	//see header file for documentation

	Double_t dx = vertex->fX - fESDEvent->GetPrimaryVertex()->GetX();
	Double_t dy = vertex->fY - fESDEvent->GetPrimaryVertex()->GetY();
	Double_t dz = vertex->fZ - fESDEvent->GetPrimaryVertex()->GetZ();
	
	//	cout<<"dx, dy, dz: ["<<dx<<","<<dy<<","<<dz<<"]"<<endl;
	particle->X() = particle->GetX() - dx;
	particle->Y() = particle->GetY() - dy;
	particle->Z() = particle->GetZ() - dz;
}

void AliAnalysisTaskGammaConversion::RotateKFParticle(AliKFParticle * kfParticle,Double_t angle){
	// Before rotate needs to be moved to position 0,0,0, ; move back after rotation
	Double_t dx = fESDEvent->GetPrimaryVertex()->GetX()-0.;
	Double_t dy = fESDEvent->GetPrimaryVertex()->GetY()-0.;
	Double_t dz = fESDEvent->GetPrimaryVertex()->GetZ()-0.;
	
	kfParticle->X() = kfParticle->GetX() - dx;
	kfParticle->Y() = kfParticle->GetY() - dy;
	kfParticle->Z() = kfParticle->GetZ() - dz;


	// Rotate the kf particle
	Double_t c = cos(angle);
	Double_t s = sin(angle);
	
	Double_t mA[8][ 8];
	for( Int_t i=0; i<8; i++ ){
		for( Int_t j=0; j<8; j++){
			mA[i][j] = 0;
		}
	}
	for( int i=0; i<8; i++ ){
		mA[i][i] = 1;
	}
	mA[0][0] =	c;	mA[0][1] = s;
	mA[1][0] = -s;	mA[1][1] = c;
	mA[3][3] =	c;	mA[3][4] = s;
	mA[4][3] = -s;	mA[4][4] = c;
	
	Double_t mAC[8][8];
	Double_t mAp[8];
	for( Int_t i=0; i<8; i++ ){
		mAp[i] = 0;
		for( Int_t k=0; k<8; k++){
			mAp[i]+=mA[i][k] * kfParticle->GetParameter(k);
		}
	}
	
	for( Int_t i=0; i<8; i++){
		kfParticle->Parameter(i) = mAp[i];
	}

	for( Int_t i=0; i<8; i++ ){
		for( Int_t j=0; j<8; j++ ){
			mAC[i][j] = 0;
			for( Int_t k=0; k<8; k++ ){
	mAC[i][j]+= mA[i][k] * kfParticle->GetCovariance(k,j);
			}
		}
	}

	for( Int_t i=0; i<8; i++ ){
		for( Int_t j=0; j<=i; j++ ){
			Double_t xx = 0;
			for( Int_t k=0; k<8; k++){
	xx+= mAC[i][k]*mA[j][k];
			}
			kfParticle->Covariance(i,j) = xx;
		}
	}

	Double_t dx1 = 0.-fESDEvent->GetPrimaryVertex()->GetX();
	Double_t dy1 = 0.-fESDEvent->GetPrimaryVertex()->GetY();
	Double_t dz1 = 0.-fESDEvent->GetPrimaryVertex()->GetZ();
	
	kfParticle->X() = kfParticle->GetX() - dx1;
	kfParticle->Y() = kfParticle->GetY() - dy1;
	kfParticle->Z() = kfParticle->GetZ() - dz1;

}


void AliAnalysisTaskGammaConversion::CalculateBackground(){
	// see header file for documentation


	TClonesArray * currentEventV0s = fV0Reader->GetCurrentEventGoodV0s();

	AliGammaConversionBGHandler * bgHandler = fV0Reader->GetBGHandler();
	
	Int_t zbin= bgHandler->GetZBinIndex(fV0Reader->GetVertexZ());
	Int_t mbin = 0;
	if(fUseTrackMultiplicityForBG == kTRUE){
		mbin = bgHandler->GetMultiplicityBinIndex(fV0Reader->CountESDTracks());
	}
	else{
		mbin = bgHandler->GetMultiplicityBinIndex(fV0Reader->GetNGoodV0s());
	}

	if(fDoRotation == kTRUE){

		for(Int_t iCurrent=0;iCurrent<currentEventV0s->GetEntriesFast();iCurrent++){
			AliKFParticle currentEventGoodV0 = *(AliKFParticle *)(currentEventV0s->At(iCurrent)); 
			for(Int_t iCurrent2=iCurrent+1;iCurrent2<currentEventV0s->GetEntriesFast();iCurrent2++){
	for(Int_t nRandom=0;nRandom<fNRandomEventsForBG;nRandom++){
	
		AliKFParticle currentEventGoodV02 = *(AliKFParticle *)(currentEventV0s->At(iCurrent2));

		if(fCheckBGProbability == kTRUE){
			Double_t massBGprob =0.;
			Double_t widthBGprob = 0.;
			AliKFParticle *backgroundCandidateProb = new AliKFParticle(currentEventGoodV0,currentEventGoodV02);
			backgroundCandidateProb->GetMass(massBGprob,widthBGprob);
			if(massBGprob>0.1 && massBGprob<0.14){
				if(fRandom.Rndm()>bgHandler->GetBGProb(zbin,mbin)){
		delete backgroundCandidateProb;
		continue;
				}
			}
			delete backgroundCandidateProb;
		}
	
		Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;

		Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
		
		RotateKFParticle(&currentEventGoodV02,rotationValue);

		AliKFParticle *backgroundCandidate = new AliKFParticle(currentEventGoodV0,currentEventGoodV02);

		Double_t massBG =0.;
		Double_t widthBG = 0.;
		Double_t chi2BG =10000.;	
		backgroundCandidate->GetMass(massBG,widthBG);
		//		if(backgroundCandidate->GetNDF()>0){
		chi2BG = backgroundCandidate->GetChi2();
		if((chi2BG>0 && chi2BG<fV0Reader->GetChi2CutMeson())	|| fApplyChi2Cut == kFALSE){
		
			TVector3 momentumVectorbackgroundCandidate(backgroundCandidate->GetPx(),backgroundCandidate->GetPy(),backgroundCandidate->GetPz());
			TVector3 spaceVectorbackgroundCandidate(backgroundCandidate->GetX(),backgroundCandidate->GetY(),backgroundCandidate->GetZ());
		
			Double_t openingAngleBG = currentEventGoodV0.GetAngle(currentEventGoodV02);
		
			Double_t rapidity;
			if(backgroundCandidate->GetE() - backgroundCandidate->GetPz() == 0 || backgroundCandidate->GetE() + backgroundCandidate->GetPz() == 0) {
				rapidity=8.;
			} else{
				rapidity = 0.5*(TMath::Log((backgroundCandidate->GetE() +backgroundCandidate->GetPz()) / (backgroundCandidate->GetE()-backgroundCandidate->GetPz())));
			}
			if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut() ){
				delete backgroundCandidate;	 
				continue;	 // rapidity cut
			}			
					
		
			Double_t alfa=0.0;
			if( (currentEventGoodV0.GetE()+currentEventGoodV02.GetE()) != 0){
				alfa=TMath::Abs((currentEventGoodV0.GetE()-currentEventGoodV02.GetE())
						/(currentEventGoodV0.GetE()+currentEventGoodV02.GetE()));
			}
		
		
			if(openingAngleBG < fMinOpeningAngleGhostCut ){
				delete backgroundCandidate;	 
				continue;	 // minimum opening angle to avoid using ghosttracks
			}			
		
			// original
			if(alfa>fV0Reader->GetAlphaMinCutMeson() &&	 alfa<fV0Reader->GetAlphaCutMeson()){
				fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle", openingAngleBG);
				fHistograms->FillHistogram("ESD_Background_Energy", backgroundCandidate->GetE());
				fHistograms->FillHistogram("ESD_Background_Pt",	momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram("ESD_Background_Eta", momentumVectorbackgroundCandidate.Eta());
				fHistograms->FillHistogram("ESD_Background_Rapidity", rapidity);
				fHistograms->FillHistogram("ESD_Background_Phi", spaceVectorbackgroundCandidate.Phi());
				fHistograms->FillHistogram("ESD_Background_Mass", massBG);
				fHistograms->FillHistogram("ESD_Background_R", spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
				fHistograms->FillHistogram("ESD_Background_ZR", backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram("ESD_Background_XY", backgroundCandidate->GetX(), backgroundCandidate->GetY());
				fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt",massBG,momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram("ESD_Background_InvMass",massBG);
				fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_alpha",massBG,momentumVectorbackgroundCandidate.Pt());

				if(massBG>0.1 && massBG<0.15){
		fHistograms->FillHistogram("ESD_Background_alfa_Pi0", alfa);
		fHistograms->FillHistogram("ESD_Background_Pt_alpha_Pi0", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
				}
				if(massBG>0.5 && massBG<0.57){
		fHistograms->FillHistogram("ESD_Background_alfa_Eta", alfa);
		fHistograms->FillHistogram("ESD_Background_Pt_alpha_Eta", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
				}

				if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(currentEventGoodV02.GetEta())<0.9 ){
		fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_Fiducial",massBG,momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram("ESD_Background_InvMass_Fiducial",massBG);
				}
				
				fHistograms->FillHistogram(Form("%d%dESD_Background_GammaDaughter_OpeningAngle",zbin,mbin), openingAngleBG);
				fHistograms->FillHistogram(Form("%d%dESD_Background_Energy",zbin,mbin), backgroundCandidate->GetE());
				fHistograms->FillHistogram(Form("%d%dESD_Background_Pt",zbin,mbin),	momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram(Form("%d%dESD_Background_Eta",zbin,mbin), momentumVectorbackgroundCandidate.Eta());
				fHistograms->FillHistogram(Form("%d%dESD_Background_Rapidity",zbin,mbin), rapidity);
				fHistograms->FillHistogram(Form("%d%dESD_Background_Phi",zbin,mbin), spaceVectorbackgroundCandidate.Phi());
				fHistograms->FillHistogram(Form("%d%dESD_Background_Mass",zbin,mbin), massBG);
				fHistograms->FillHistogram(Form("%d%dESD_Background_R",zbin,mbin), spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
				fHistograms->FillHistogram(Form("%d%dESD_Background_ZR",zbin,mbin), backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram(Form("%d%dESD_Background_XY",zbin,mbin), backgroundCandidate->GetX(), backgroundCandidate->GetY());
				fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass",zbin,mbin),massBG);
				
				if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(currentEventGoodV02.GetEta())<0.9 ){
		fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt_Fiducial",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_Fiducial",zbin,mbin),massBG);
				}
			}
			if(alfa<0.1){
				fHistograms->FillHistogram("ESD_Background_InvMass_vs_E_alpha",massBG ,backgroundCandidate->GetE());
			}

		}
		//}
		delete backgroundCandidate;			
	}
			}
		}
	}
	else{ // means no rotation
		AliGammaConversionBGHandler::GammaConversionVertex *bgEventVertex = NULL;
			
		if(fUseTrackMultiplicityForBG){
			//		cout<<"Using charged track multiplicity for background calculation"<<endl;
			for(Int_t nEventsInBG=0;nEventsInBG <fV0Reader->GetNBGEvents();nEventsInBG++){

	AliGammaConversionKFVector * previousEventV0s = bgHandler->GetBGGoodV0s(zbin,mbin,nEventsInBG);//fV0Reader->GetBGGoodV0s(nEventsInBG);
			
	if(fMoveParticleAccordingToVertex == kTRUE){
		bgEventVertex = bgHandler->GetBGEventVertex(zbin,mbin,nEventsInBG);
	}

	for(Int_t iCurrent=0;iCurrent<currentEventV0s->GetEntriesFast();iCurrent++){
		AliKFParticle currentEventGoodV0 = *(AliKFParticle *)(currentEventV0s->At(iCurrent)); 
		for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
			AliKFParticle previousGoodV0 = (AliKFParticle)(*(previousEventV0s->at(iPrevious)));
			AliKFParticle previousGoodV0test = (AliKFParticle)(*(previousEventV0s->at(iPrevious)));

			//cout<<"Primary Vertex event: ["<<fESDEvent->GetPrimaryVertex()->GetX()<<","<<fESDEvent->GetPrimaryVertex()->GetY()<<","<<fESDEvent->GetPrimaryVertex()->GetZ()<<"]"<<endl;
			//cout<<"BG prim Vertex event: ["<<bgEventVertex->fX<<","<<bgEventVertex->fY<<","<<bgEventVertex->fZ<<"]"<<endl;
		
			//cout<<"XYZ of particle before transport: ["<<previousGoodV0.X()<<","<<previousGoodV0.Y()<<","<<previousGoodV0.Z()<<"]"<<endl;
			if(fMoveParticleAccordingToVertex == kTRUE){
				MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
			}
			//cout<<"XYZ of particle after transport: ["<<previousGoodV0.X()<<","<<previousGoodV0.Y()<<","<<previousGoodV0.Z()<<"]"<<endl;

			AliKFParticle *backgroundCandidate = new AliKFParticle(currentEventGoodV0,previousGoodV0);
	
			Double_t massBG =0.;
			Double_t widthBG = 0.;
			Double_t chi2BG =10000.;	
			backgroundCandidate->GetMass(massBG,widthBG);

			//	if(backgroundCandidate->GetNDF()>0){
			//		chi2BG = backgroundCandidate->GetChi2()/backgroundCandidate->GetNDF();
			chi2BG = backgroundCandidate->GetChi2();
			if((chi2BG>0 && chi2BG<fV0Reader->GetChi2CutMeson()) || fApplyChi2Cut == kFALSE){
					
				TVector3 momentumVectorbackgroundCandidate(backgroundCandidate->GetPx(),backgroundCandidate->GetPy(),backgroundCandidate->GetPz());
				TVector3 spaceVectorbackgroundCandidate(backgroundCandidate->GetX(),backgroundCandidate->GetY(),backgroundCandidate->GetZ());
					
				Double_t openingAngleBG = currentEventGoodV0.GetAngle(previousGoodV0);
					
				Double_t rapidity;
			
				if(backgroundCandidate->GetE() - backgroundCandidate->GetPz() <= 0 || backgroundCandidate->GetE() + backgroundCandidate->GetPz() <= 0){
		cout << "Error: |Pz| > E !!!! " << endl;
		rapidity=8.;
				} else {
		rapidity = 0.5*(TMath::Log((backgroundCandidate->GetE() +backgroundCandidate->GetPz()) / (backgroundCandidate->GetE()-backgroundCandidate->GetPz())));
				}				
				if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut() ){
		delete backgroundCandidate;	 
		continue;	 // rapidity cut
				}			
							
	
				Double_t alfa=0.0;
				if( (currentEventGoodV0.GetE()+previousGoodV0.GetE()) != 0){
		alfa=TMath::Abs((currentEventGoodV0.GetE()-previousGoodV0.GetE())
				/(currentEventGoodV0.GetE()+previousGoodV0.GetE()));
				}
			
					
				if(openingAngleBG < fMinOpeningAngleGhostCut ){
		delete backgroundCandidate;	 
		continue;	 // minimum opening angle to avoid using ghosttracks
				}			

				// original
				if(alfa>fV0Reader->GetAlphaMinCutMeson() &&	 alfa<fV0Reader->GetAlphaCutMeson()){
		fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle", openingAngleBG);
		fHistograms->FillHistogram("ESD_Background_Energy", backgroundCandidate->GetE());
		fHistograms->FillHistogram("ESD_Background_Pt",	momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram("ESD_Background_Eta", momentumVectorbackgroundCandidate.Eta());
		fHistograms->FillHistogram("ESD_Background_Rapidity", rapidity);
		fHistograms->FillHistogram("ESD_Background_Phi", spaceVectorbackgroundCandidate.Phi());
		fHistograms->FillHistogram("ESD_Background_Mass", massBG);
		fHistograms->FillHistogram("ESD_Background_R", spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
		fHistograms->FillHistogram("ESD_Background_ZR", backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram("ESD_Background_XY", backgroundCandidate->GetX(), backgroundCandidate->GetY());
		fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt",massBG,momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram("ESD_Background_InvMass",massBG);
		fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_alpha",massBG,momentumVectorbackgroundCandidate.Pt());

		if(massBG>0.1 && massBG<0.15){
			fHistograms->FillHistogram("ESD_Background_alfa_Pi0", alfa);
			fHistograms->FillHistogram("ESD_Background_Pt_alpha_Pi0", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
		}
		if(massBG>0.5 && massBG<0.57){
			fHistograms->FillHistogram("ESD_Background_alfa_Eta", alfa);
			fHistograms->FillHistogram("ESD_Background_Pt_alpha_Eta", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
		}

		if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(previousGoodV0.GetEta())<0.9 ){
			fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_Fiducial",massBG,momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram("ESD_Background_InvMass_Fiducial",massBG);
		}

		// test
		fHistograms->FillHistogram(Form("%d%dESD_Background_GammaDaughter_OpeningAngle",zbin,mbin), openingAngleBG);
		fHistograms->FillHistogram(Form("%d%dESD_Background_Energy",zbin,mbin), backgroundCandidate->GetE());
		fHistograms->FillHistogram(Form("%d%dESD_Background_Pt",zbin,mbin),	momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram(Form("%d%dESD_Background_Eta",zbin,mbin), momentumVectorbackgroundCandidate.Eta());
		fHistograms->FillHistogram(Form("%d%dESD_Background_Rapidity",zbin,mbin), rapidity);
		fHistograms->FillHistogram(Form("%d%dESD_Background_Phi",zbin,mbin), spaceVectorbackgroundCandidate.Phi());
		fHistograms->FillHistogram(Form("%d%dESD_Background_Mass",zbin,mbin), massBG);
		fHistograms->FillHistogram(Form("%d%dESD_Background_R",zbin,mbin), spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
		fHistograms->FillHistogram(Form("%d%dESD_Background_ZR",zbin,mbin), backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram(Form("%d%dESD_Background_XY",zbin,mbin), backgroundCandidate->GetX(), backgroundCandidate->GetY());
		fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
		fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass",zbin,mbin),massBG);
		
		if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(previousGoodV0.GetEta())<0.9 ){
			fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt_Fiducial",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_Fiducial",zbin,mbin),massBG);
		}
		//		}
				}
				if(alfa<0.1){
		fHistograms->FillHistogram("ESD_Background_InvMass_vs_E_alpha",massBG ,backgroundCandidate->GetE());
				}

			}
			delete backgroundCandidate;			
		}
	}
			}
		}
		else{ // means using #V0s for multiplicity

			//		cout<<"Using the v0 multiplicity to calculate background"<<endl;
		
			fHistograms->FillHistogram("ESD_Background_z_m",zbin,mbin);
			fHistograms->FillHistogram("ESD_Mother_multpilicityVSv0s",fV0Reader->CountESDTracks(),fV0Reader->GetNumberOfV0s());

			for(Int_t nEventsInBG=0;nEventsInBG <fV0Reader->GetNBGEvents();nEventsInBG++){
	AliGammaConversionKFVector * previousEventV0s = bgHandler->GetBGGoodV0s(zbin,mbin,nEventsInBG);// fV0Reader->GetBGGoodV0s(nEventsInBG);
	if(previousEventV0s){
	
		if(fMoveParticleAccordingToVertex == kTRUE){
			bgEventVertex = bgHandler->GetBGEventVertex(zbin,mbin,nEventsInBG);
		}

		for(Int_t iCurrent=0;iCurrent<currentEventV0s->GetEntriesFast();iCurrent++){
			AliKFParticle currentEventGoodV0 = *(AliKFParticle *)(currentEventV0s->At(iCurrent)); 
			for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
				AliKFParticle previousGoodV0 = (AliKFParticle)(*(previousEventV0s->at(iPrevious)));

				if(fMoveParticleAccordingToVertex == kTRUE){
		MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
				}

				AliKFParticle *backgroundCandidate = new AliKFParticle(currentEventGoodV0,previousGoodV0);
				Double_t massBG =0.;
				Double_t widthBG = 0.;
				Double_t chi2BG =10000.;	
				backgroundCandidate->GetMass(massBG,widthBG);

				/*			if(backgroundCandidate->GetNDF()>0){
					chi2BG = backgroundCandidate->GetChi2()/backgroundCandidate->GetNDF();
					{//remember to remove
					TVector3 momentumVectorbackgroundCandidate(backgroundCandidate->GetPx(),backgroundCandidate->GetPy(),backgroundCandidate->GetPz());
					TVector3 spaceVectorbackgroundCandidate(backgroundCandidate->GetX(),backgroundCandidate->GetY(),backgroundCandidate->GetZ());
				
					Double_t openingAngleBG = currentEventGoodV0.GetAngle(previousGoodV0);
					fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle_nochi2", openingAngleBG);
					}
				*/
				chi2BG = backgroundCandidate->GetChi2();
				if((chi2BG>0 && chi2BG<fV0Reader->GetChi2CutMeson()) || fApplyChi2Cut == kFALSE){
		TVector3 momentumVectorbackgroundCandidate(backgroundCandidate->GetPx(),backgroundCandidate->GetPy(),backgroundCandidate->GetPz());
		TVector3 spaceVectorbackgroundCandidate(backgroundCandidate->GetX(),backgroundCandidate->GetY(),backgroundCandidate->GetZ());
					
		Double_t openingAngleBG = currentEventGoodV0.GetAngle(previousGoodV0);
					
		Double_t rapidity;
		if(backgroundCandidate->GetE() - backgroundCandidate->GetPz() == 0 || backgroundCandidate->GetE() + backgroundCandidate->GetPz() == 0){
			rapidity=8.;
		}else{
			rapidity = 0.5*(TMath::Log((backgroundCandidate->GetE() +backgroundCandidate->GetPz()) / (backgroundCandidate->GetE()-backgroundCandidate->GetPz())));
		}		
		if(TMath::Abs(rapidity) > fV0Reader->GetRapidityMesonCut() ){
			delete backgroundCandidate;	 
			continue;	 // rapidity cut
		}			
								

		Double_t alfa=0.0;
		if( (currentEventGoodV0.GetE()+previousGoodV0.GetE()) != 0){
			alfa=TMath::Abs((currentEventGoodV0.GetE()-previousGoodV0.GetE())
					/(currentEventGoodV0.GetE()+previousGoodV0.GetE()));
		}
			
					
		if(openingAngleBG < fMinOpeningAngleGhostCut ){
			delete backgroundCandidate;	 
			continue;	 // minimum opening angle to avoid using ghosttracks
		}			

		if(alfa>fV0Reader->GetAlphaMinCutMeson() &&	 alfa<fV0Reader->GetAlphaCutMeson()){
			fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle", openingAngleBG);
			fHistograms->FillHistogram("ESD_Background_Energy", backgroundCandidate->GetE());
			fHistograms->FillHistogram("ESD_Background_Pt",	momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram("ESD_Background_Eta", momentumVectorbackgroundCandidate.Eta());
			fHistograms->FillHistogram("ESD_Background_Rapidity", rapidity);
			fHistograms->FillHistogram("ESD_Background_Phi", spaceVectorbackgroundCandidate.Phi());
			fHistograms->FillHistogram("ESD_Background_Mass", massBG);
			fHistograms->FillHistogram("ESD_Background_R", spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
			fHistograms->FillHistogram("ESD_Background_ZR", backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram("ESD_Background_XY", backgroundCandidate->GetX(), backgroundCandidate->GetY());
			fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt",massBG,momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram("ESD_Background_InvMass",massBG);
			

			fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_alpha",massBG,momentumVectorbackgroundCandidate.Pt());

			if(massBG>0.1 && massBG<0.15){
				fHistograms->FillHistogram("ESD_Background_alfa_Pi0", alfa);
				fHistograms->FillHistogram("ESD_Background_Pt_alpha_Pi0", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
			}
			if(massBG>0.5 && massBG<0.57){
				fHistograms->FillHistogram("ESD_Background_alfa_Eta", alfa);
				fHistograms->FillHistogram("ESD_Background_Pt_alpha_Eta", momentumVectorbackgroundCandidate.Pt(), alfa); //RR_alpha
			}

			if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(previousGoodV0.GetEta())<0.9 ){
				fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt_Fiducial",massBG,momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram("ESD_Background_InvMass_Fiducial",massBG);
			}
			
			if(massBG>0.5 && massBG<0.6){
				fHistograms->FillHistogram("ESD_Background_alfa_pt0506",momentumVectorbackgroundCandidate.Pt(),alfa);
			}
			if(massBG>0.3 && massBG<0.4){
				fHistograms->FillHistogram("ESD_Background_alfa_pt0304",momentumVectorbackgroundCandidate.Pt(),alfa);
			}
			
			// test
			fHistograms->FillHistogram(Form("%d%dESD_Background_GammaDaughter_OpeningAngle",zbin,mbin), openingAngleBG);
			fHistograms->FillHistogram(Form("%d%dESD_Background_Energy",zbin,mbin), backgroundCandidate->GetE());
			fHistograms->FillHistogram(Form("%d%dESD_Background_Pt",zbin,mbin),	momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram(Form("%d%dESD_Background_Eta",zbin,mbin), momentumVectorbackgroundCandidate.Eta());
			fHistograms->FillHistogram(Form("%d%dESD_Background_Rapidity",zbin,mbin), rapidity);
			fHistograms->FillHistogram(Form("%d%dESD_Background_Phi",zbin,mbin), spaceVectorbackgroundCandidate.Phi());
			fHistograms->FillHistogram(Form("%d%dESD_Background_Mass",zbin,mbin), massBG);
			fHistograms->FillHistogram(Form("%d%dESD_Background_R",zbin,mbin), spaceVectorbackgroundCandidate.Pt());	// Pt in Space == R!!!!
			fHistograms->FillHistogram(Form("%d%dESD_Background_ZR",zbin,mbin), backgroundCandidate->GetZ(), spaceVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram(Form("%d%dESD_Background_XY",zbin,mbin), backgroundCandidate->GetX(), backgroundCandidate->GetY());
			fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
			fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass",zbin,mbin),massBG);
			
			if ( TMath::Abs(currentEventGoodV0.GetEta())<0.9 &&	TMath::Abs(previousGoodV0.GetEta())<0.9 ){
				fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_vs_Pt_Fiducial",zbin,mbin),massBG,momentumVectorbackgroundCandidate.Pt());
				fHistograms->FillHistogram(Form("%d%dESD_Background_InvMass_Fiducial",zbin,mbin),massBG);
			}
		}

		if(alfa<0.1){
			fHistograms->FillHistogram("ESD_Background_InvMass_vs_E_alpha",massBG ,backgroundCandidate->GetE());
		}
		//	}
				}
				delete backgroundCandidate;			
			}
		}
	}
			}
		} // end else (means use #v0s as multiplicity)
	} // end no rotation
}


void AliAnalysisTaskGammaConversion::ProcessGammasForGammaJetAnalysis(){
	//ProcessGammasForGammaJetAnalysis
	
	Double_t distIsoMin;
	
	CreateListOfChargedParticles();
	
	
	//	for(UInt_t gammaIndex=0;gammaIndex<fKFReconstructedGammas.size();gammaIndex++){
	for(Int_t gammaIndex=0;gammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();gammaIndex++){
		AliKFParticle * currentGamma = (AliKFParticle*)fKFReconstructedGammasTClone->At(gammaIndex);
		TVector3 momentumVectorCurrentGamma(currentGamma->GetPx(),currentGamma->GetPy(),currentGamma->GetPz());
		if( momentumVectorCurrentGamma.Pt()> fMinPtForGammaJet){
			distIsoMin=GetMinimumDistanceToCharge(gammaIndex);
			if (distIsoMin > fMinIsoConeSize && fLeadingChargedIndex>=0){
				CalculateJetCone(gammaIndex);
			}
		}
	}
}

//____________________________________________________________________
Bool_t AliAnalysisTaskGammaConversion::IsGoodImpPar(const AliESDtrack *const track)
{
	//
	// check whether particle has good DCAr(Pt) impact
	// parameter. Only for TPC+ITS tracks (7*sigma cut)
	// Origin: Andrea Dainese
	//

	Float_t d0z0[2],covd0z0[3];
	track->GetImpactParameters(d0z0,covd0z0);
	Float_t sigma= 0.0050+0.0060/TMath::Power(track->Pt(),0.9);
	Float_t d0max = 7.*sigma;
	if(TMath::Abs(d0z0[0]) < d0max) return kTRUE;

	return kFALSE;
}


void AliAnalysisTaskGammaConversion::CreateListOfChargedParticles(){
	// CreateListOfChargedParticles
	
	fESDEvent = fV0Reader->GetESDEvent();
	Int_t numberOfESDTracks=0;
	for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
		AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
		
		if(!curTrack){
			continue;
		}
		// Not needed if Standard function used.
		//		 if(!IsGoodImpPar(curTrack)){
		//			 continue;
		//		 }
		
		if(fEsdTrackCuts->AcceptTrack(curTrack) ){
			new((*fChargedParticles)[fChargedParticles->GetEntriesFast()])	AliESDtrack(*curTrack);
			//			fChargedParticles.push_back(curTrack);
			fChargedParticlesId.push_back(iTracks);
			numberOfESDTracks++;
		}
	}
	// Moved to UserExec using CountAcceptedTracks function. runjet is not needed by default
	//	 fHistograms->FillHistogram("ESD_NumberOfGoodESDTracks",numberOfESDTracks);
	//	 cout<<"esdtracks::"<< numberOfESDTracks<<endl;
	//	 if (fV0Reader->GetNumberOfContributorsVtx()>=1){
	//		 fHistograms->FillHistogram("ESD_NumberOfGoodESDTracksVtx",numberOfESDTracks);
	//	 } 
}

/*void AliAnalysisTaskGammaConversion::RecalculateV0ForGamma(){
	//recalculates v0 for gamma

	Double_t massE=0.00051099892;
	TLorentzVector curElecPos;
	TLorentzVector curElecNeg;
	TLorentzVector curGamma;

	TLorentzVector curGammaAt;
	TLorentzVector curElecPosAt;
	TLorentzVector curElecNegAt;
	AliKFVertex primVtxGamma(*(fESDEvent->GetPrimaryVertex()));
	AliKFVertex primVtxImprovedGamma = primVtxGamma;

	const AliESDVertex *vtxT3D=fESDEvent->GetPrimaryVertex();

	Double_t xPrimaryVertex=vtxT3D->GetXv();
	Double_t yPrimaryVertex=vtxT3D->GetYv();
	Double_t zPrimaryVertex=vtxT3D->GetZv();
	// Float_t primvertex[3]={xPrimaryVertex,yPrimaryVertex,zPrimaryVertex};

	Float_t nsigmaTPCtrackPos;
	Float_t nsigmaTPCtrackNeg;
	Float_t nsigmaTPCtrackPosToPion;
	Float_t nsigmaTPCtrackNegToPion;
	AliKFParticle* negKF=NULL;
	AliKFParticle* posKF=NULL;

	for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
		AliESDtrack* posTrack = fESDEvent->GetTrack(iTracks);
		if(!posTrack){
			continue;
		}
		if (posKF) delete posKF; posKF=NULL;
		if(posTrack->GetSign()<0) continue;
		if(!(posTrack->GetStatus() & AliESDtrack::kTPCrefit))continue;
		if(posTrack->GetKinkIndex(0)>0 ) continue;
		if(posTrack->GetNcls(1)<50)continue;
		Double_t momPos[3];
		//		posTrack->GetConstrainedPxPyPz(momPos);
		posTrack->GetPxPyPz(momPos);
		AliESDtrack *ptrk=fESDEvent->GetTrack(iTracks);
		curElecPos.SetXYZM(momPos[0],momPos[1],momPos[2],massE);
		if(TMath::Abs(curElecPos.Eta())<0.9) continue;
		posKF = new AliKFParticle( *(posTrack),-11);

		nsigmaTPCtrackPos = fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
		nsigmaTPCtrackPosToPion = fV0Reader->GetESDpid()->NumberOfSigmasTPC(posTrack,AliPID::kPion);

		if ( nsigmaTPCtrackPos>5.|| nsigmaTPCtrackPos<-2.){
			continue;
		}
	
		if(pow((momPos[0]*momPos[0]+momPos[1]*momPos[1]+momPos[2]*momPos[2]),0.5)>0.5 && nsigmaTPCtrackPosToPion<1){
			continue;
		}



		for(Int_t jTracks = 0; jTracks < fESDEvent->GetNumberOfTracks(); jTracks++){
			AliESDtrack* negTrack = fESDEvent->GetTrack(jTracks);
			if(!negTrack){
				continue;
			}
			if (negKF) delete negKF; negKF=NULL;
			if(negTrack->GetSign()>0) continue;
			if(!(negTrack->GetStatus() & AliESDtrack::kTPCrefit))continue;
			if(negTrack->GetKinkIndex(0)>0 ) continue;
			if(negTrack->GetNcls(1)<50)continue;
			Double_t momNeg[3];
			//		negTrack->GetConstrainedPxPyPz(momNeg);
			negTrack->GetPxPyPz(momNeg);

			nsigmaTPCtrackNeg = fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kElectron);		 
			nsigmaTPCtrackNegToPion = fV0Reader->GetESDpid()->NumberOfSigmasTPC(negTrack,AliPID::kPion);
			if ( nsigmaTPCtrackNeg>5. || nsigmaTPCtrackNeg<-2.){
				continue;
			}
			if(pow((momNeg[0]*momNeg[0]+momNeg[1]*momNeg[1]+momNeg[2]*momNeg[2]),0.5)>0.5 && nsigmaTPCtrackNegToPion<1){
				continue;
			}
			AliESDtrack *ntrk=fESDEvent->GetTrack(jTracks);
			curElecNeg.SetXYZM(momNeg[0],momNeg[1],momNeg[2],massE);
			if(TMath::Abs(curElecNeg.Eta())<0.9) continue;
			negKF = new AliKFParticle( *(negTrack) ,11);

			Double_t b=fESDEvent->GetMagneticField();
			Double_t xn, xp, dca=ntrk->GetDCA(ptrk,b,xn,xp);
			AliExternalTrackParam nt(*ntrk), pt(*ptrk);
			nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);


			//--- Like in ITSV0Finder
			AliExternalTrackParam ntAt0(*ntrk), ptAt0(*ptrk);
			Double_t xxP,yyP,alphaP;
			Double_t rP[3];

			//		 if (!ptAt0.GetGlobalXYZat(ptAt0->GetX(),xxP,yyP,zzP)) continue;
			if (!ptAt0.GetXYZAt(ptAt0.GetX(),b,rP)) continue;
			xxP=rP[0];
			yyP=rP[1];
			alphaP = TMath::ATan2(yyP,xxP);


			ptAt0.Propagate(alphaP,0,b);
			Float_t ptfacP	= (1.+100.*TMath::Abs(ptAt0.GetC(b)));

			//		 Double_t distP			= ptAt0.GetY();
			Double_t normP			= ptfacP*TMath::Sqrt(ptAt0.GetSigmaY2());
			Double_t normdist0P = TMath::Abs(ptAt0.GetY()/normP);
			Double_t normdist1P = TMath::Abs((ptAt0.GetZ()-zPrimaryVertex)/(ptfacP*TMath::Sqrt(ptAt0.GetSigmaZ2())));
			Double_t normdistP	= TMath::Sqrt(normdist0P*normdist0P+normdist1P*normdist1P);
	

			Double_t xxN,yyN,alphaN;
			Double_t rN[3];
			//		 if (!ntAt0.GetGlobalXYZat(ntAt0->GetX(),xxN,yyN,zzN)) continue;
			if (!ntAt0.GetXYZAt(ntAt0.GetX(),b,rN)) continue;
			xxN=rN[0];
			yyN=rN[1];
 
			alphaN = TMath::ATan2(yyN,xxN);

			ntAt0.Propagate(alphaN,0,b);

			Float_t ptfacN	= (1.+100.*TMath::Abs(ntAt0.GetC(b)));
			//		 Double_t distN			= ntAt0.GetY();
			Double_t normN			= ptfacN*TMath::Sqrt(ntAt0.GetSigmaY2());
			Double_t normdist0N = TMath::Abs(ntAt0.GetY()/normN);
			Double_t normdist1N = TMath::Abs((ntAt0.GetZ()-zPrimaryVertex)/(ptfacN*TMath::Sqrt(ntAt0.GetSigmaZ2())));
			Double_t normdistN	= TMath::Sqrt(normdist0N*normdist0N+normdist1N*normdist1N);
	
			//-----------------------------

			Double_t momNegAt[3];
			nt.GetPxPyPz(momNegAt);
			curElecNegAt.SetXYZM(momNegAt[0],momNegAt[1],momNegAt[2],massE);

			Double_t momPosAt[3];
			pt.GetPxPyPz(momPosAt);
			curElecPosAt.SetXYZM(momPosAt[0],momPosAt[1],momPosAt[2],massE);
			if(dca>1){
				continue;
			}

			//		 Double_t dneg= negTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
			//		 Double_t dpos= posTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
			AliESDv0 vertex(nt,jTracks,pt,iTracks);
		

			Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);

 

			//	cout<< "v0Rr::"<< v0Rr<<endl;
			// if (pvertex.GetRr()<0.5){
			// continue;
			//}
			//		 cout<<"vertex.GetChi2V0()"<<vertex.GetChi2V0()<<endl;
			if(cpa<0.9)continue;
			//		 if (vertex.GetChi2V0() > 30) continue;
			//		 cout<<"xp+xn::"<<xp<<" "<<xn<<endl;
			if ((xn+xp) < 0.4) continue;
			if (TMath::Abs(ntrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<0.05)
			if (TMath::Abs(ptrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<0.05) continue;

			//cout<<"pass"<<endl;

			AliKFParticle v0GammaC;
			v0GammaC+=(*negKF);
			v0GammaC+=(*posKF);
			v0GammaC.SetMassConstraint(0,0.001);
			primVtxImprovedGamma+=v0GammaC;
			v0GammaC.SetProductionVertex(primVtxImprovedGamma);


			curGamma=curElecNeg+curElecPos;
			curGammaAt=curElecNegAt+curElecPosAt;
		 
			// invariant mass versus pt of K0short
		 
			Double_t chi2V0GammaC=100000.;
			if( v0GammaC.GetNDF() != 0) {
				chi2V0GammaC = v0GammaC.GetChi2()/v0GammaC.GetNDF();
			}else{
				cout<< "ERROR::v0K0C.GetNDF()" << endl;
			}

			if(chi2V0GammaC<200 &&chi2V0GammaC>0 ){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_RecalculateV0_InvMass",v0GammaC.GetMass());
					fHistograms->FillHistogram("ESD_RecalculateV0_Pt",v0GammaC.GetPt());
					fHistograms->FillHistogram("ESD_RecalculateV0_E_dEdxP",curElecNegAt.P(),negTrack->GetTPCsignal());
					fHistograms->FillHistogram("ESD_RecalculateV0_P_dEdxP",curElecPosAt.P(),posTrack->GetTPCsignal());
					fHistograms->FillHistogram("ESD_RecalculateV0_cpa",cpa);
					fHistograms->FillHistogram("ESD_RecalculateV0_dca",dca);
					fHistograms->FillHistogram("ESD_RecalculateV0_normdistP",normdistP);
					fHistograms->FillHistogram("ESD_RecalculateV0_normdistN",normdistN);

					new((*fKFRecalculatedGammasTClone)[fKFRecalculatedGammasTClone->GetEntriesFast()])	AliKFParticle(v0GammaC);
					fElectronRecalculatedv1.push_back(iTracks);
					fElectronRecalculatedv2.push_back(jTracks);
				}
			}
		}
	}
 
	for(Int_t firstGammaIndex=0;firstGammaIndex<fKFRecalculatedGammasTClone->GetEntriesFast();firstGammaIndex++){
		for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFRecalculatedGammasTClone->GetEntriesFast();secondGammaIndex++){
			AliKFParticle * twoGammaDecayCandidateDaughter0 = (AliKFParticle *)fKFRecalculatedGammasTClone->At(firstGammaIndex);
			AliKFParticle * twoGammaDecayCandidateDaughter1 = (AliKFParticle *)fKFRecalculatedGammasTClone->At(secondGammaIndex);
			
			if(fElectronRecalculatedv1[firstGammaIndex]==fElectronRecalculatedv1[secondGammaIndex]){
				continue;
			}
			if( fElectronRecalculatedv2[firstGammaIndex]==fElectronRecalculatedv2[secondGammaIndex]){
				continue;
			}
			
			AliKFParticle twoGammaCandidate(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_RecalculateGG_InvMass",twoGammaCandidate.GetMass());		
				fHistograms->FillHistogram("ESD_RecalculateGG_InvMass_vs_Pt",twoGammaCandidate.GetMass(),twoGammaCandidate.GetPt());		
			}
		}
	}
}
        */
void AliAnalysisTaskGammaConversion::CalculateJetCone(Int_t gammaIndex){
	// CaculateJetCone
	
	Double_t cone;
	Double_t coneSize=0.3;
	Double_t ptJet=0;
	
	//	AliKFParticle * currentGamma = &fKFReconstructedGammas[gammaIndex];
	AliKFParticle * currentGamma = (AliKFParticle*)fKFReconstructedGammasTClone->At(gammaIndex);

	TVector3 momentumVectorCurrentGamma(currentGamma->GetPx(),currentGamma->GetPy(),currentGamma->GetPz());
	
	AliESDtrack* leadingCharged = (AliESDtrack*)(fChargedParticles->At(fLeadingChargedIndex));

	Double_t momLeadingCharged[3];
	leadingCharged->GetConstrainedPxPyPz(momLeadingCharged);
	
	TVector3 momentumVectorLeadingCharged(momLeadingCharged[0],momLeadingCharged[1],momLeadingCharged[2]);
	
	Double_t phi1=momentumVectorLeadingCharged.Phi();
	Double_t eta1=momentumVectorLeadingCharged.Eta();
	Double_t phi3=momentumVectorCurrentGamma.Phi();
	
	for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
		AliESDtrack* curTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
		Int_t chId = fChargedParticlesId[iCh];
		if(fLeadingChargedIndex==chId || fLeadingChargedIndex==chId) continue;
		Double_t mom[3];
		curTrack->GetConstrainedPxPyPz(mom);
		TVector3 momentumVectorChargedParticle(mom[0],mom[1],mom[2]);
		Double_t phi2=momentumVectorChargedParticle.Phi();
		Double_t eta2=momentumVectorChargedParticle.Eta();
		
		
		cone=100.;
		if( TMath::Abs(phi2 - phi1) <= ( TMath::TwoPi()-coneSize) ){
			cone = TMath::Sqrt(	TMath::Power((eta2-eta1),2)+ TMath::Power((phi2-phi1),2) );
		}else{
			if( (phi2 - phi1)> TMath::TwoPi()-coneSize ){
				cone = TMath::Sqrt(	TMath::Power((eta2-eta1),2)+ TMath::Power((phi2-TMath::TwoPi()-phi1),2) );
			}
			if( (phi2 - phi1)< -(TMath::TwoPi()-coneSize) ){
				cone = TMath::Sqrt(	TMath::Power((eta2-eta1),2)+ TMath::Power((phi2+TMath::TwoPi()-phi1),2) );
			}
		}
		
		if(cone <coneSize&& momentumVectorChargedParticle.Pt()>fMinPtJetCone ){
			ptJet+= momentumVectorChargedParticle.Pt();
			Double_t ffzHdrGam = momentumVectorChargedParticle.Pt()/momentumVectorCurrentGamma.Pt();
			Double_t imbalanceHdrGam=-momentumVectorChargedParticle.Dot(momentumVectorCurrentGamma)/momentumVectorCurrentGamma.Mag2();
			fHistograms->FillHistogram("ESD_FFzHdrGam",ffzHdrGam);
			fHistograms->FillHistogram("ESD_ImbalanceHdrGam",imbalanceHdrGam);
			
		}
		
		Double_t dphiHdrGam=phi3-phi2;
		if ( dphiHdrGam < (-TMath::PiOver2())){
			dphiHdrGam+=(TMath::TwoPi());
		}
		
		if ( dphiHdrGam > (3.*TMath::PiOver2()) ){
			dphiHdrGam-=(TMath::TwoPi());
		}
		
		if (momentumVectorChargedParticle.Pt()>fMinPtGamChargedCorr){
			fHistograms->FillHistogram("ESD_dphiHdrGamIsolated",dphiHdrGam);
		}
	}//track loop
	
	
	}



Double_t AliAnalysisTaskGammaConversion::GetMinimumDistanceToCharge(Int_t indexHighestPtGamma){
	// GetMinimumDistanceToCharge
	
	Double_t fIsoMin=100.;
	Double_t ptLeadingCharged=-1.;

	fLeadingChargedIndex=-1;
	
	AliKFConversionPhoton * gammaHighestPt = (AliKFConversionPhoton*)fKFReconstructedGammasTClone->At(indexHighestPtGamma);
	TVector3 momentumVectorgammaHighestPt(gammaHighestPt->GetPx(),gammaHighestPt->GetPy(),gammaHighestPt->GetPz());
	
	Double_t phi1=momentumVectorgammaHighestPt.Phi();
	Double_t eta1=momentumVectorgammaHighestPt.Eta();
	
	for(Int_t iCh=0;iCh<fChargedParticles->GetEntriesFast();iCh++){
		AliESDtrack* curTrack = (AliESDtrack*)(fChargedParticles->At(iCh));
		Int_t chId = fChargedParticlesId[iCh];
		if(gammaHighestPt->GetTrackLabelPositive()==chId || gammaHighestPt->GetTrackLabelNegative()==chId) continue;
		Double_t mom[3];
		curTrack->GetConstrainedPxPyPz(mom);
		TVector3 momentumVectorChargedParticle(mom[0],mom[1],mom[2]);
		Double_t phi2=momentumVectorChargedParticle.Phi();
		Double_t eta2=momentumVectorChargedParticle.Eta();
		Double_t iso=pow(	(pow( (eta1-eta2),2)+ pow((phi1-phi2),2)),0.5 );
		
		if(momentumVectorChargedParticle.Pt()>fMinPtIsoCone ){
			if (iso<fIsoMin){
	fIsoMin=iso;
			}
		}
		
		Double_t dphiHdrGam=phi1-phi2;
		if ( dphiHdrGam < (-TMath::PiOver2())){
			dphiHdrGam+=(TMath::TwoPi());
		}
		
		if ( dphiHdrGam > (3.*TMath::PiOver2()) ){
			dphiHdrGam-=(TMath::TwoPi());
		}
		if (momentumVectorChargedParticle.Pt()>fMinPtGamChargedCorr){
			fHistograms->FillHistogram("ESD_dphiHdrGam",dphiHdrGam);
		}
		
		if (dphiHdrGam>0.9*TMath::Pi() && dphiHdrGam<1.1*TMath::Pi()){
			if (momentumVectorChargedParticle.Pt()> ptLeadingCharged && momentumVectorChargedParticle.Pt()>0.1*momentumVectorgammaHighestPt.Pt()){
				ptLeadingCharged=momentumVectorChargedParticle.Pt();
				fLeadingChargedIndex=iCh;
			}
		}
		
	}//track loop
	fHistograms->FillHistogram("ESD_MinimumIsoDistance",fIsoMin);
	return fIsoMin;
	
}

Int_t	AliAnalysisTaskGammaConversion::GetIndexHighestPtGamma(){
	//GetIndexHighestPtGamma
	
	Int_t indexHighestPtGamma=-1;
	//Double_t 
	fGammaPtHighest = -100.;
	
	for(Int_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammasTClone->GetEntriesFast();firstGammaIndex++){
		AliKFParticle * gammaHighestPtCandidate = (AliKFParticle*)fKFReconstructedGammasTClone->At(firstGammaIndex);
		TVector3 momentumVectorgammaHighestPtCandidate(gammaHighestPtCandidate->GetPx(),gammaHighestPtCandidate->GetPy(),gammaHighestPtCandidate->GetPz());
		if (momentumVectorgammaHighestPtCandidate.Pt() > fGammaPtHighest){
			fGammaPtHighest=momentumVectorgammaHighestPtCandidate.Pt();
			//gammaHighestPt = gammaHighestPtCandidate;
			indexHighestPtGamma=firstGammaIndex;
		}
	}
	
	return indexHighestPtGamma;
	
}


void AliAnalysisTaskGammaConversion::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	AliDebug(1,"Do nothing in Terminate");
}

void AliAnalysisTaskGammaConversion::UserCreateOutputObjects()
{
	
	if(fKFCreateAOD) {

		//AOD
		if(!fAODGamma) fAODGamma = new TClonesArray("AliAODConversionPhoton", 0);
		else fAODGamma->Delete();
		fAODGamma->SetName(Form("%s_gamma", fAODBranchName.Data()));
		
	       /* if(!fAODPi0) fAODPi0 = new TClonesArray("AliAODConversionPhoton", 0);
		else fAODPi0->Delete();
		fAODPi0->SetName(Form("%s_Pi0", fAODBranchName.Data()));
		
		if(!fAODOmega) fAODOmega = new TClonesArray("AliAODConversionPhoton", 0);
		else fAODOmega->Delete();
		fAODOmega->SetName(Form("%s_Omega", fAODBranchName.Data()));
	       */
		//If delta AOD file name set, add in separate file. Else add in standard aod file. 
		if(GetDeltaAODFileName().Length() > 0) {
			AddAODBranch("TClonesArray", &fAODGamma, GetDeltaAODFileName().Data());
		      //  AddAODBranch("TClonesArray", &fAODPi0, GetDeltaAODFileName().Data());
		      //  AddAODBranch("TClonesArray", &fAODOmega, GetDeltaAODFileName().Data());
			AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(GetDeltaAODFileName().Data());
		} else	{
			AddAODBranch("TClonesArray", &fAODGamma);
		      //  AddAODBranch("TClonesArray", &fAODPi0);
		      //  AddAODBranch("TClonesArray", &fAODOmega);
		}
	}

	// Create the output container
	if(fOutputContainer != NULL){
		delete fOutputContainer;
		fOutputContainer = NULL;
	}
	if(fOutputContainer == NULL){
		fOutputContainer = new TList();
		fOutputContainer->SetOwner(kTRUE);
	}
	
	//Adding the histograms to the output container
	fHistograms->GetOutputContainer(fOutputContainer);
	
	
	if(fWriteNtuple){
		if(fGammaNtuple == NULL){
			fGammaNtuple = new TNtuple("V0ntuple","V0ntuple","OnTheFly:HasVertex:NegPIDProb:PosPIDProb:X:Y:Z:R:MotherCandidateNDF:MotherCandidateChi2:MotherCandidateEnergy:MotherCandidateEta:MotherCandidatePt:MotherCandidateMass:MotherCandidateWidth:MCMotherCandidatePT:EPOpeningAngle:ElectronEnergy:ElectronPt:ElectronEta:ElectronPhi:PositronEnergy:PositronPt:PositronEta:PositronPhi:HasSameMCMother:MotherMCParticlePIDCode",50000);
		}
		if(fNeutralMesonNtuple == NULL){
			fNeutralMesonNtuple = new TNtuple("NeutralMesonNtuple","NeutralMesonNtuple","test");
		}
		TList * ntupleTList = new TList();
		ntupleTList->SetOwner(kTRUE);
		ntupleTList->SetName("Ntuple");
		ntupleTList->Add((TNtuple*)fGammaNtuple);
		fOutputContainer->Add(ntupleTList);
	}
	
	fOutputContainer->SetName(GetName());

	PostData(1, fOutputContainer);
	PostData(2, fCFManager->GetParticleContainer());	// for CF

}

Double_t AliAnalysisTaskGammaConversion::GetMCOpeningAngle(const TParticle* const daughter0, const TParticle* const daughter1) const{
	//helper function
	TVector3 v3D0(daughter0->Px(),daughter0->Py(),daughter0->Pz());
	TVector3 v3D1(daughter1->Px(),daughter1->Py(),daughter1->Pz());
	return v3D0.Angle(v3D1);
}

void AliAnalysisTaskGammaConversion::CheckV0Efficiency(){
	// see header file for documentation

	vector<Int_t> indexOfGammaParticle;
	
	fStack = fV0Reader->GetMCStack();
	
	if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
		return; // aborts if the primary vertex does not have contributors.
	}
	
	for (Int_t iTracks = 0; iTracks < fStack->GetNprimary(); iTracks++) {
		TParticle* particle = (TParticle *)fStack->Particle(iTracks);
		if(particle->GetPdgCode()==22){		 //Gamma
			if(particle->GetNDaughters() >= 2){
				TParticle* electron=NULL;
				TParticle* positron=NULL; 
				for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
					TParticle *tmpDaughter = fStack->Particle(daughterIndex);
					if(tmpDaughter->GetUniqueID() == 5){
						if(tmpDaughter->GetPdgCode() == 11){
							electron = tmpDaughter;
						}
						else if(tmpDaughter->GetPdgCode() == -11){
							positron = tmpDaughter;
						}
					}
				}
				if(electron!=NULL && positron!=0){
					if(electron->R()<160){
						indexOfGammaParticle.push_back(iTracks);
					}
				}
			}
		}
	}
	
	Int_t nFoundGammas=0;
	Int_t nNotFoundGammas=0;
	
	Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
	for(Int_t i=0;i<numberOfV0s;i++){
		fV0Reader->GetV0(i);
		
		if(fV0Reader->HasSameMCMother() == kFALSE){
			continue;
		}
		
		TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
		TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
		
		if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
			continue;
		}
		if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
			continue;
		}
		
		if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
			//TParticle * v0Gamma = fV0Reader->GetMotherMCParticle();
			for(UInt_t mcIndex=0;mcIndex<indexOfGammaParticle.size();mcIndex++){
				if(negativeMC->GetFirstMother()==indexOfGammaParticle[mcIndex]){
					nFoundGammas++;
				}
				else{
					nNotFoundGammas++;
				}
			}
		}
	}
}


void AliAnalysisTaskGammaConversion::ProcessGammaElectronsForChicAnalysis(){
	// see header file for documantation
	
	fESDEvent = fV0Reader->GetESDEvent();
	
	
	TClonesArray * vESDeNegTemp = new TClonesArray("AliESDtrack",0);
	TClonesArray * vESDePosTemp = new TClonesArray("AliESDtrack",0);
	TClonesArray * vESDxNegTemp = new TClonesArray("AliESDtrack",0);
	TClonesArray * vESDxPosTemp = new TClonesArray("AliESDtrack",0);
	TClonesArray * vESDeNegNoJPsi = new TClonesArray("AliESDtrack",0);
	TClonesArray * vESDePosNoJPsi = new TClonesArray("AliESDtrack",0);
	
	/*
		vector <AliESDtrack*> vESDeNegTemp(0);
		vector <AliESDtrack*> vESDePosTemp(0);
		vector <AliESDtrack*> vESDxNegTemp(0);
		vector <AliESDtrack*> vESDxPosTemp(0);
		vector <AliESDtrack*> vESDeNegNoJPsi(0);
		vector <AliESDtrack*> vESDePosNoJPsi(0); 
	*/
	
	
	fHistograms->FillTable("Table_Electrons",0);//Count number of Events
	
	for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
		AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
		
		if(!curTrack){
			//print warning here
			continue;
		}
		
		double p[3];if(!curTrack->GetConstrainedPxPyPz(p))continue;
		double r[3];curTrack->GetConstrainedXYZ(r);
		
		TVector3 rXYZ(r);
		
		fHistograms->FillTable("Table_Electrons",4);//Count number of ESD tracks
		
		Bool_t flagKink			 =	kTRUE;
		Bool_t flagTPCrefit	 =	kTRUE;
		Bool_t flagTRDrefit	 =	kTRUE;
		Bool_t flagITSrefit	 =	kTRUE;
		Bool_t flagTRDout		 =	kTRUE;
		Bool_t flagVertex		 =	kTRUE;
		
		
		//Cuts ---------------------------------------------------------------
		
		if(curTrack->GetKinkIndex(0) > 0){
			fHistograms->FillHistogram("Table_Electrons",5);//Count kink
			flagKink = kFALSE;
		}
		
		ULong_t trkStatus = curTrack->GetStatus();
		
		ULong_t tpcRefit = (trkStatus & AliESDtrack::kTPCrefit);
		
		if(!tpcRefit){
			fHistograms->FillHistogram("Table_Electrons",9);//Count not TPCrefit
			flagTPCrefit = kFALSE;
		}
		
		ULong_t itsRefit = (trkStatus & AliESDtrack::kITSrefit);
		if(!itsRefit){
			fHistograms->FillHistogram("Table_Electrons",10);//Count not ITSrefit
			flagITSrefit = kFALSE;
		}
		
		ULong_t trdRefit = (trkStatus & AliESDtrack::kTRDrefit);
		
		if(!trdRefit){
			fHistograms->FillHistogram("Table_Electrons",8); //Count not TRDrefit
			flagTRDrefit = kFALSE;
		}
		
		ULong_t trdOut = (trkStatus & AliESDtrack::kTRDout);
		
		if(!trdOut) {
			fHistograms->FillHistogram("Table_Electrons",7); //Count not TRDout
			flagTRDout = kFALSE;
		}
		
		double nsigmaToVxt = GetSigmaToVertex(curTrack);
		
		if(nsigmaToVxt > 3){
			fHistograms->FillHistogram("Table_Electrons",6); //Count Tracks with number of sigmas > 3
			flagVertex = kFALSE;
		}
		
		if(! (flagKink && flagTPCrefit && flagITSrefit && flagTRDrefit && flagTRDout && flagVertex ) ) continue;
		fHistograms->FillHistogram("Table_Electrons",11);//Count Tracks passed Cuts
		
		
		Stat_t pid, weight;
		GetPID(curTrack, pid, weight);
		
		if(pid!=0){
			fHistograms->FillHistogram("Table_Electrons",12); //Count Tracks with pid != 0
		}
		
		if(pid == 0){
			fHistograms->FillHistogram("Table_Electrons",13); //Count Tracks with pid != 0
		}
		
		
		
		
		
		
		TLorentzVector curElec;
		curElec.SetXYZM(p[0],p[1],p[2],fElectronMass);
		
		
		if(fDoMCTruth){		
			Int_t labelMC = TMath::Abs(curTrack->GetLabel());
			TParticle* curParticle = fStack->Particle(labelMC);
			if(curTrack->GetSign() > 0){
				if( pid == 0){
					fHistograms->FillHistogram("MC_ElectronPosNegPt",curParticle->Pt());
					fHistograms->FillHistogram("MC_ElectronPosNegEta",curParticle->Eta());
				} else {
					fHistograms->FillHistogram("MC_ElectronPosNegPt",curParticle->Pt());
					fHistograms->FillHistogram("MC_ElectronPosNegEta",curParticle->Eta());
				}
			}
		}
		
		
		if(curTrack->GetSign() > 0){
			
			//		 vESDxPosTemp.push_back(curTrack);
			new((*vESDxPosTemp)[vESDxPosTemp->GetEntriesFast()])	AliESDtrack(*curTrack);
			
			if( pid == 0){
				fHistograms->FillHistogram("ESD_ElectronPosNegPt",curElec.Pt());
				fHistograms->FillHistogram("ESD_ElectronPosPt",curElec.Pt());
				//	fHistograms->FillHistogram("MC_ElectronPosNegPt",curParticle->Pt());
				fHistograms->FillHistogram("ESD_ElectronPosNegEta",curElec.Eta());
				//	fHistograms->FillHistogram("MC_ElectronPosNegEta",curParticle->Eta());
				//	vESDePosTemp.push_back(curTrack);
				new((*vESDePosTemp)[vESDePosTemp->GetEntriesFast()])	AliESDtrack(*curTrack);
			}
			
		} else {

			new((*vESDxNegTemp)[vESDxNegTemp->GetEntriesFast()])	AliESDtrack(*curTrack);
			
			if( pid == 0){
					
				fHistograms->FillHistogram("ESD_ElectronPosNegPt",curElec.Pt());
				fHistograms->FillHistogram("ESD_ElectronNegPt",curElec.Pt());
				fHistograms->FillHistogram("ESD_ElectronPosNegEta",curElec.Eta());
				new((*vESDeNegTemp)[vESDeNegTemp->GetEntriesFast()])	AliESDtrack(*curTrack);
					
			}
			
		}
		
	}
	
	
	Bool_t ePosJPsi = kFALSE;
	Bool_t eNegJPsi = kFALSE;		
	Bool_t ePosPi0	= kFALSE;
	Bool_t eNegPi0	= kFALSE;
	
	UInt_t iePosJPsi=0,ieNegJPsi=0,iePosPi0=0,ieNegPi0=0;
	
	for(Int_t iNeg=0; iNeg < vESDeNegTemp->GetEntriesFast(); iNeg++){
		if(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDeNegTemp->At(iNeg)))->GetLabel()))->GetPdgCode() == 11)
		if(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDeNegTemp->At(iNeg)))->GetLabel()))->GetMother(0) > -1){
		Int_t labelMother = fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDeNegTemp->At(iNeg)))->GetLabel()))->GetMother(0);
			TParticle* partMother = fStack ->Particle(labelMother);
			if (partMother->GetPdgCode() == 111){
				ieNegPi0 = iNeg;
				eNegPi0 = kTRUE;
			}
			if(partMother->GetPdgCode() == 443){ //Mother JPsi
				fHistograms->FillTable("Table_Electrons",14);
				ieNegJPsi = iNeg;
				eNegJPsi = kTRUE;
			} else {	
				//		vESDeNegNoJPsi.push_back(vESDeNegTemp[iNeg]);
				new((*vESDeNegNoJPsi)[vESDeNegNoJPsi->GetEntriesFast()])	AliESDtrack(*(AliESDtrack*)(vESDeNegTemp->At(iNeg)));
				//		cout<<"ESD No Positivo JPsi "<<endl;
			}
		}
	}	
	
	for(Int_t iPos=0; iPos < vESDePosTemp->GetEntriesFast(); iPos++){
		if(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDePosTemp->At(iPos)))->GetLabel()))->GetPdgCode() == -11)
		if(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDePosTemp->At(iPos)))->GetLabel()))->GetMother(0) > -1){
			Int_t labelMother = fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDePosTemp->At(iPos)))->GetLabel()))->GetMother(0);
			TParticle* partMother = fStack ->Particle(labelMother);
			if (partMother->GetPdgCode() == 111){
				iePosPi0 = iPos;
				ePosPi0 = kTRUE;
			}
			if(partMother->GetPdgCode() == 443){ //Mother JPsi
				fHistograms->FillTable("Table_Electrons",15);
				iePosJPsi = iPos;
				ePosJPsi = kTRUE;
			}
			else{
				//		vESDePosNoJPsi.push_back(vESDePosTemp[iPos]);
				new((*vESDePosNoJPsi)[vESDePosNoJPsi->GetEntriesFast()])	AliESDtrack(*(AliESDtrack*)(vESDePosTemp->At(iPos)));		
				//		cout<<"ESD No Negativo JPsi "<<endl;
			}				
		}
	}
	
	if( eNegJPsi && ePosJPsi ){
		TVector3 tempeNegV,tempePosV;
		tempeNegV.SetXYZ(((AliESDtrack*)(vESDeNegTemp->At(ieNegJPsi)))->Px(),((AliESDtrack*)(vESDeNegTemp->At(ieNegJPsi)))->Py(),((AliESDtrack*)(vESDeNegTemp->At(ieNegJPsi)))->Pz());			
		tempePosV.SetXYZ(((AliESDtrack*)(vESDePosTemp->At(iePosJPsi)))->Px(),((AliESDtrack*)(vESDePosTemp->At(iePosJPsi)))->Py(),((AliESDtrack*)(vESDePosTemp->At(iePosJPsi)))->Pz());
		fHistograms->FillTable("Table_Electrons",16);
		fHistograms->FillHistogram("ESD_ElectronPosNegJPsiAngle",tempeNegV.Angle(tempePosV));	
		fHistograms->FillHistogram("MC_ElectronPosNegJPsiAngle",GetMCOpeningAngle(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDeNegTemp->At(ieNegJPsi)))->GetLabel())),
												fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDePosTemp->At(iePosJPsi)))->GetLabel()))));	
	}
	
	if( eNegPi0 && ePosPi0 ){
		TVector3 tempeNegV,tempePosV;
		tempeNegV.SetXYZ(((AliESDtrack*)(vESDeNegTemp->At(ieNegPi0)))->Px(),((AliESDtrack*)(vESDeNegTemp->At(ieNegPi0)))->Py(),((AliESDtrack*)(vESDeNegTemp->At(ieNegPi0)))->Pz());
		tempePosV.SetXYZ(((AliESDtrack*)(vESDePosTemp->At(iePosPi0)))->Px(),((AliESDtrack*)(vESDePosTemp->At(iePosPi0)))->Py(),((AliESDtrack*)(vESDePosTemp->At(iePosPi0)))->Pz());
		fHistograms->FillHistogram("ESD_ElectronPosNegPi0Angle",tempeNegV.Angle(tempePosV));
		fHistograms->FillHistogram("MC_ElectronPosNegPi0Angle",GetMCOpeningAngle(fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDeNegTemp->At(ieNegPi0)))->GetLabel())),
											 fStack->Particle(TMath::Abs(((AliESDtrack*)(vESDePosTemp->At(iePosPi0)))->GetLabel()))));	 
	}
	
	
	FillAngle("ESD_eNegePosAngleBeforeCut",GetTLorentzVector(vESDeNegTemp),GetTLorentzVector(vESDePosTemp));
	
	CleanWithAngleCuts(*vESDeNegTemp,*vESDePosTemp,*fKFReconstructedGammasTClone);
	
	//	vector <TLorentzVector> vCurrentTLVeNeg = GetTLorentzVector(fCurrentEventNegElectron);
	//	vector <TLorentzVector> vCurrentTLVePos = GetTLorentzVector(fCurrentEventPosElectron);
	
	TClonesArray vCurrentTLVeNeg = GetTLorentzVector(fCurrentEventNegElectronTClone);
	TClonesArray vCurrentTLVePos = GetTLorentzVector(fCurrentEventPosElectronTClone);
	
	
	FillAngle("ESD_eNegePosAngleAfterCut",vCurrentTLVeNeg,vCurrentTLVePos);
	
	
	
	
	//FillAngle("ESD_eNegePosAngleAfterCut",CurrentTLVeNeg,CurrentTLVePos);
	
	
	FillElectronInvMass("ESD_InvMass_ePluseMinus",vCurrentTLVeNeg,vCurrentTLVePos);
	FillElectronInvMass("ESD_InvMass_xPlusxMinus",GetTLorentzVector(vESDxNegTemp),GetTLorentzVector(vESDxPosTemp));
	
	
	
	FillGammaElectronInvMass("ESD_InvMass_GammaePluseMinusChiC","ESD_InvMass_GammaePluseMinusChiCDiff",*fKFReconstructedGammasCutTClone,vCurrentTLVeNeg,vCurrentTLVePos);
	
	FillGammaElectronInvMass("ESD_InvMass_GammaePluseMinusPi0","ESD_InvMass_GammaePluseMinusPi0Diff",
				 *fKFReconstructedGammasCutTClone,vCurrentTLVeNeg,vCurrentTLVePos);
	
	//BackGround
	
	//Like Sign e+e-
	ElectronBackground("ESD_ENegBackground",vCurrentTLVeNeg);
	ElectronBackground("ESD_EPosBackground",vCurrentTLVePos);
	ElectronBackground("ESD_EPosENegBackground",vCurrentTLVeNeg);
	ElectronBackground("ESD_EPosENegBackground",vCurrentTLVePos);
	
	//				Like Sign e+e- no JPsi
	ElectronBackground("ESD_EPosENegNoJPsiBG",GetTLorentzVector(vESDeNegNoJPsi));
	ElectronBackground("ESD_EPosENegNoJPsiBG",GetTLorentzVector(vESDePosNoJPsi));
	
	//Mixed Event
	
	if( fCurrentEventPosElectronTClone->GetEntriesFast() > 0 && fCurrentEventNegElectronTClone->GetEntriesFast() > 0 && fKFReconstructedGammasCutTClone->GetEntriesFast() > 0 ){
		FillGammaElectronInvMass("ESD_EPosENegGammaBackgroundMX","ESD_EPosENegGammaBackgroundMXDiff",
					 *fKFReconstructedGammasCutTClone,*fPreviousEventTLVNegElectronTClone,*fPreviousEventTLVPosElectronTClone);
		*fPreviousEventTLVNegElectronTClone = vCurrentTLVeNeg;
		*fPreviousEventTLVPosElectronTClone = vCurrentTLVePos;
		
	}
	
	/*
	//Photons P
	Double_t vtx[3];
	vtx[0]=0;vtx[1]=0;vtx[2]=0;
	for(UInt_t i=0;i<fKFReconstructedGammasChic.size();i++){
	 
	//			if(fMCGammaChicTempCut[i]->GetMother(0) < 0) continue;
	 
	 
	 
	Int_t tempLabel = fStack->Particle(fMCGammaChicTempCut[i]->GetMother(0))->GetPdgCode();
	//			cout<<" Label Pedro Gonzalez " <<tempLabel <<endl;
	 
	//			cout<<" Label Distance"<<fKFReconstructedGammasChic[i].GetDistanceFromVertex(vtx)<<endl;
	 
	if( tempLabel == 10441 || tempLabel == 20443 || tempLabel == 445 )
	 
	fHistograms->FillHistogram("ESD_PhotonsMomentum",fKFReconstructedGammasChic[i].GetMomentum());
	 
	 
	}
	 
	 
	*/


	vESDeNegTemp->Delete();
	vESDePosTemp->Delete();
	vESDxNegTemp->Delete();
	vESDxPosTemp->Delete();
	vESDeNegNoJPsi->Delete();
	vESDePosNoJPsi->Delete();

	delete vESDeNegTemp;
	delete vESDePosTemp;
	delete vESDxNegTemp;
	delete vESDxPosTemp;
	delete vESDeNegNoJPsi;
	delete vESDePosNoJPsi;	
}

/*
	void AliAnalysisTaskGammaConversion::FillAngle(TString histoName,vector <TLorentzVector> tlVeNeg, vector <TLorentzVector> tlVePos){
	//see header file for documentation
	for( UInt_t iNeg=0; iNeg < tlVeNeg.size(); iNeg++){
	for (UInt_t iPos=0; iPos < tlVePos.size(); iPos++){
	fHistograms->FillHistogram(histoName.Data(),tlVeNeg[iNeg].Vect().Angle(tlVePos[iPos].Vect()));
	}
	}
	}
*/
void AliAnalysisTaskGammaConversion::FillAngle(TString histoName,TClonesArray const tlVeNeg, TClonesArray const tlVePos){
	//see header file for documentation
	for( Int_t iNeg=0; iNeg < tlVeNeg.GetEntriesFast(); iNeg++){
		for (Int_t iPos=0; iPos < tlVePos.GetEntriesFast(); iPos++){
			fHistograms->FillHistogram(histoName.Data(),((TLorentzVector*)(tlVeNeg.At(iNeg)))->Vect().Angle(((TLorentzVector*)(tlVePos.At(iPos)))->Vect()));
		}
	}
}

void AliAnalysisTaskGammaConversion::FillElectronInvMass(TString histoName, TClonesArray const eNeg, TClonesArray const ePos){
	//see header file for documentation
	for( Int_t n=0; n < eNeg.GetEntriesFast(); n++){
		TLorentzVector en = (*(TLorentzVector*)(eNeg.At(n)));
		for (Int_t p=0; p < ePos.GetEntriesFast(); p++){
			TLorentzVector ep = (*(TLorentzVector*)(ePos.At(p)));
			TLorentzVector np = ep + en;
			fHistograms->FillHistogram(histoName.Data(),np.M());
		}
	}
}

void AliAnalysisTaskGammaConversion::FillGammaElectronInvMass(TString histoMass,TString histoDiff,TClonesArray const fKFGammas,
										TClonesArray const tlVeNeg,TClonesArray const tlVePos)
{
	//see header file for documentation
	
	for( Int_t iNeg=0; iNeg < tlVeNeg.GetEntriesFast(); iNeg++ ){
		
		for (Int_t iPos=0; iPos < tlVePos.GetEntriesFast(); iPos++){
			
			TLorentzVector xy = *((TLorentzVector *)(tlVePos.At(iPos))) + *((TLorentzVector *)(tlVeNeg.At(iNeg)));
			
			for (Int_t iGam=0; iGam < fKFGammas.GetEntriesFast(); iGam++){
				
	//	AliKFParticle * gammaCandidate = &fKFGammas[iGam];
	AliKFParticle * gammaCandidate = (AliKFParticle *)(fKFGammas.At(iGam));
	TLorentzVector g;
				
	g.SetXYZM(gammaCandidate->GetPx(),gammaCandidate->GetPy(),gammaCandidate->GetPz(),fGammaMass);
	TLorentzVector xyg = xy + g;
	fHistograms->FillHistogram(histoMass.Data(),xyg.M());
	fHistograms->FillHistogram(histoDiff.Data(),(xyg.M()-xy.M()));
			}
		}
	}
	
}
void AliAnalysisTaskGammaConversion::ElectronBackground(TString hBg, TClonesArray e)
{
	// see header file for documentation
	for(Int_t i=0; i < e.GetEntriesFast(); i++) {
		for (Int_t j=i+1; j < e.GetEntriesFast(); j++) {
			TLorentzVector ee = (*(TLorentzVector*)(e.At(i))) + (*(TLorentzVector*)(e.At(j)));
			fHistograms->FillHistogram(hBg.Data(),ee.M());
		}
	}
}


void AliAnalysisTaskGammaConversion::CleanWithAngleCuts(TClonesArray const negativeElectrons,
							TClonesArray const positiveElectrons, 
							TClonesArray const gammas){
	// see header file for documentation
	
	UInt_t	sizeN = negativeElectrons.GetEntriesFast();
	UInt_t	sizeP = positiveElectrons.GetEntriesFast();
	UInt_t	sizeG = gammas.GetEntriesFast();
	
	
	
	vector <Bool_t> xNegBand(sizeN);
	vector <Bool_t> xPosBand(sizeP);
	vector <Bool_t> gammaBand(sizeG);
	
	
	for(UInt_t iNeg=0; iNeg < sizeN; iNeg++) xNegBand[iNeg]=kTRUE;
	for(UInt_t iPos=0; iPos < sizeP; iPos++) xPosBand[iPos]=kTRUE;
	for(UInt_t iGam=0; iGam < sizeG; iGam++) gammaBand[iGam]=kTRUE;
	
	
	for(UInt_t iPos=0; iPos < sizeP; iPos++){
		
		Double_t aP[3]; 
		((AliESDtrack*)(positiveElectrons.At(iPos)))->GetConstrainedPxPyPz(aP); 
		
		TVector3 ePosV(aP[0],aP[1],aP[2]);
		
		for(UInt_t iNeg=0; iNeg < sizeN; iNeg++){
			
			Double_t aN[3]; 
			((AliESDtrack*)(negativeElectrons.At(iNeg)))->GetConstrainedPxPyPz(aN); 
			TVector3 eNegV(aN[0],aN[1],aN[2]);
			
			if(ePosV.Angle(eNegV) < 0.05){ //e+e- from gamma
				xPosBand[iPos]=kFALSE;
				xNegBand[iNeg]=kFALSE;
			}
			
			for(UInt_t iGam=0; iGam < sizeG; iGam++){
				AliKFParticle* gammaCandidate = (AliKFParticle*)gammas.At(iGam);
				TVector3 gammaCandidateVector(gammaCandidate->Px(),gammaCandidate->Py(),gammaCandidate->Pz());
				if(ePosV.Angle(gammaCandidateVector) < 0.05 || eNegV.Angle(gammaCandidateVector) < 0.05)
					gammaBand[iGam]=kFALSE;
			}
		}
	}
	
	for(UInt_t iPos=0; iPos < sizeP; iPos++){
		if(xPosBand[iPos]){
			new((*fCurrentEventPosElectronTClone)[fCurrentEventPosElectronTClone->GetEntriesFast()]) AliESDtrack((*(AliESDtrack*)(positiveElectrons.At(iPos))));
			//			fCurrentEventPosElectron.push_back(positiveElectrons[iPos]);
		}
	}
	for(UInt_t iNeg=0;iNeg < sizeN; iNeg++){
		if(xNegBand[iNeg]){
			new((*fCurrentEventNegElectronTClone)[fCurrentEventNegElectronTClone->GetEntriesFast()]) AliESDtrack((*(AliESDtrack*)(negativeElectrons.At(iNeg))));
			//			fCurrentEventNegElectron.push_back(negativeElectrons[iNeg]);
		}
	}
	for(UInt_t iGam=0; iGam < sizeG; iGam++){
		if(gammaBand[iGam]){
			new((*fKFReconstructedGammasCutTClone)[fKFReconstructedGammasCutTClone->GetEntriesFast()]) AliKFParticle((*(AliKFParticle*)(gammas.At(iGam))));
			//fKFReconstructedGammasCut.push_back(*(AliKFParticle*)gammas->At(iGam));
		}
	}
}


void	AliAnalysisTaskGammaConversion::GetPID(const AliESDtrack *track, Stat_t &pid, Stat_t &weight)
{
	// see header file for documentation
	pid = -1;
	weight = -1;
	
	double wpart[5];
	double wpartbayes[5];
	
	//get probability of the diffenrent particle types
	track->GetESDpid(wpart);
	
	// Tentative particle type "concentrations"
	double c[5]={0.01, 0.01, 0.85, 0.10, 0.05};
	
	//Bayes' formula
	double rcc = 0.;
	for (int i = 0; i < 5; i++){
		rcc+=(c[i] * wpart[i]);
	}
	
	
	
	for (int i=0; i<5; i++) {
		if( rcc>0 || rcc<0){//Kenneth: not sure if the rcc<0 should be there, this is from fixing a coding violation where we are not allowed to say: rcc!=0 (RC19)	
			wpartbayes[i] = c[i] * wpart[i] / rcc;
		}
	}
	
	Float_t max=0.;
	int ipid=-1;
	//find most probable particle in ESD pid
	//0:Electron - 1:Muon - 2:Pion - 3:Kaon - 4:Proton
	for (int i = 0; i < 5; i++){
			if (wpartbayes[i] > max){
				ipid = i;
				max = wpartbayes[i];
			}
		}
	
	pid = ipid;
	weight = max;
}

double AliAnalysisTaskGammaConversion::GetSigmaToVertex(const AliESDtrack* t)
{
	// Calculates the number of sigma to the vertex.
	
	Float_t b[2];
	Float_t bRes[2];
	Float_t bCov[3];
	t->GetImpactParameters(b,bCov);
	if (bCov[0]<=0 || bCov[2]<=0) {
		AliDebug(1, "Estimated b resolution lower or equal zero!");
		bCov[0]=0; bCov[2]=0;
	}
	bRes[0] = TMath::Sqrt(bCov[0]);
	bRes[1] = TMath::Sqrt(bCov[2]);
	
	// -----------------------------------
	// How to get to a n-sigma cut?
	//
	// The accumulated statistics from 0 to d is
	//
	// ->	Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
	// ->	1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
	//
	// It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
	// Can this be expressed in a different way?
	
	if (bRes[0] == 0 || bRes[1] ==0)
		return -1;
	
	double d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
	
	// stupid rounding problem screws up everything:
	// if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
	if (TMath::Exp(-d * d / 2) < 1e-10)
		return 1000;
	
	
	d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
	return d;
}

//vector <TLorentzVector> AliAnalysisTaskGammaConversion::GetTLorentzVector(vector <AliESDtrack*> esdTrack){
TClonesArray AliAnalysisTaskGammaConversion::GetTLorentzVector(TClonesArray *const esdTrack){
	//Return TLoresntz vector of track?
	//	vector <TLorentzVector> tlVtrack(0);
	TClonesArray array("TLorentzVector",0); 
	
	for(Int_t itrack=0; itrack < esdTrack->GetEntriesFast(); itrack++){
		double p[3]; 
		//esdTrack[itrack]->GetConstrainedPxPyPz(p);
		((AliESDtrack*)(esdTrack->At(itrack)))->GetConstrainedPxPyPz(p);
		TLorentzVector currentTrack;
		currentTrack.SetXYZM(p[0],p[1],p[2],fElectronMass);
		new((array)[array.GetEntriesFast()])	TLorentzVector(currentTrack);
		//		tlVtrack.push_back(currentTrack);
	}
	
	return array;
	
	//	return tlVtrack;
}
Int_t AliAnalysisTaskGammaConversion::GetProcessType(const AliMCEvent * mcEvt) {

	// Determine if the event was generated with pythia or phojet and return the process type

	// Check if mcEvt is fine
	if (!mcEvt) { // coverty does not allow this, the check is done elsewhere
		AliFatal("NULL mc event");
		return -1;
	} 

	// Determine if it was a pythia or phojet header, and return the correct process type
	AliGenPythiaEventHeader * headPy	= 0;
	AliGenDPMjetEventHeader * headPho = 0;
	AliGenEventHeader * htmp = mcEvt->GenEventHeader();
	if(!htmp) {
		AliFatal("Cannot Get MC Header!!");
		return -1;
	}
	if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
		headPy =	(AliGenPythiaEventHeader*) htmp;
	} else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
		headPho = (AliGenDPMjetEventHeader*) htmp;
	} else {
		AliError("Unknown header");
	}

	// Determine process type
	if(headPy)	 {
		if(headPy->ProcessType() == 92 || headPy->ProcessType() == 93) {
			// single difractive
			return kProcSD;
		} else if (headPy->ProcessType() == 94) {
			// double diffractive
			return kProcDD;
		}
		else if(headPy->ProcessType() != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94) {		
			// non difractive
			return kProcND; 
		}
	} else if (headPho) {
		if(headPho->ProcessType() == 5 || headPho->ProcessType() == 6 ) {
			// single difractive
			return kProcSD;
		} else if (headPho->ProcessType() == 7) { 
			// double diffractive
			return kProcDD;			
		} else if(headPho->ProcessType() != 5 && headPho->ProcessType() != 6	&& headPho->ProcessType() != 7 ) {
			// non difractive
			return kProcND; 
		}			 
	}
	

	// no process type found?
	AliError(Form("Unknown header: %s", htmp->IsA()->GetName()));
	return kProcUnknown;
}


Int_t AliAnalysisTaskGammaConversion::CalculateMultiplicityBin(){
	// Get Centrality bin

	Int_t multiplicity = 0;

	if ( fUseMultiplicity == 1 ) {

		if (fMultiplicity>= 0 && fMultiplicity<= 5) multiplicity=1;
		if (fMultiplicity>= 6 && fMultiplicity<= 9) multiplicity=2;
		if (fMultiplicity>=10 && fMultiplicity<=14) multiplicity=3;
		if (fMultiplicity>=15 && fMultiplicity<=22) multiplicity=4;
		if (fMultiplicity>=23 )										 multiplicity=5;

	}

	if ( fUseHBTMultiplicity == 1 ) {
		
			if (fMultiplicity>= 0 && fMultiplicity<=11)	multiplicity=1;
			if (fMultiplicity>=12 && fMultiplicity<=16)	multiplicity=2;
			if (fMultiplicity>=17 && fMultiplicity<=22)	multiplicity=3;
			if (fMultiplicity>=23 && fMultiplicity<=29)	multiplicity=4;
			if (fMultiplicity>=30 && fMultiplicity<=36)	multiplicity=5;
			if (fMultiplicity>=37 && fMultiplicity<=44)	multiplicity=6;
			if (fMultiplicity>=45 && fMultiplicity<=57)	multiplicity=7;
			if (fMultiplicity>=58 && fMultiplicity<=149) multiplicity=8;
		
	        /*
		if (fMultiplicity>= 0 && fMultiplicity<=5) multiplicity=1;
		if (fMultiplicity>=6 && fMultiplicity<=11) multiplicity=2;
		if (fMultiplicity>=12 && fMultiplicity<=16) multiplicity=3;
		if (fMultiplicity>=17 && fMultiplicity<=22) multiplicity=4;
		if (fMultiplicity>=23  ) 		 multiplicity=5;
		*/

	}

	return multiplicity;
}
