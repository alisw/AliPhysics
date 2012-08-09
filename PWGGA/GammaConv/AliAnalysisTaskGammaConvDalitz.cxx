/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Pedro González, Pedro Ladrón de Guevara, Ernesto López Torres, *
 *         Eulogio Serradilla                                             *
 * Version 2                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include <vector>

#include "TParticle.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"

#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliGammaConversionHistograms.h"
#include "AliV0Reader.h"
#include "AliKFVertex.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAnalysisTaskGammaConvDalitz.h"
#include "TH1.h"

ClassImp( AliAnalysisTaskGammaConvDalitz )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitz::AliAnalysisTaskGammaConvDalitz():
	AliAnalysisTaskSE(),
	fStack(0),
	fGCMCEvent(0),
	fESDEvent(0),
	fEposCandidateIndex(),
	fEnegCandidateIndex(),
	fGammaCandidatePosIndex(),
	fGammaCandidateNegIndex(),
	fGammaCandidates(0),
	fGammaPool(0),
	fPoolMaxSize(10),
	fGamPoolPos(0),
	fBGEventHandler(0),
	fOutputContainer(0),
	fMCTruth(0),
	fV0Reader(0),
	fESDpid(0),
	fESDtrackCuts(0),
	fITSsaTrackCuts(0),
	fESDpidCuts(0),
	fRemovePileUp(kFALSE),
	fSelectV0AND(kFALSE),                 // Select V0AND
	fTriggerAnalysis(0),
	fMultiplicity(0),
	fUseMultiplicity(0),
	fUseMultiplicityBin(0),
	fUseHBTMultiplicity(0),
	fUseHBTMultiplicityBin(0),
	fUseCentrality(0),
	fUseCentralityBin(0),   

	fHistograms(0),
	fStandalone(kFALSE),
	fDoMC(kFALSE),
	fComputeBkg(kTRUE),
	fUseBayesPID(kFALSE),
	fUseTrackIndexCut(kTRUE),
	fUsePsiPairCut(kTRUE),
	fUseMassCut(kFALSE),
	fUseGammaCut(kFALSE),
	fUseCorrectedTPCClsInfo(kFALSE),
	fReadMagFieldSign(kTRUE),
	fUseAliKF(kFALSE),
	fBackgroundType(2),
	fMinClsTPCToF(0.),
	fMagFieldSign(1),
	fkElectronMass(0.00051099891),
	fPsiPairCut(0.45),
	fDeltaPhiCutMin(0.),
	fDeltaPhiCutMax(0.12),
	fMassCutMin(0.),
	fMassCutMax(0.1),
	/*fNSigmaBelowElecTPCbethe(-2.),
	fNSigmaAboveElecTPCbethe(3.),
	fNSigmaAbovePionTPCbethe(3.),
	fNSigmaAboveKaonTPCbethe(3.),
	fNSigmaAboveProtonTPCbethe(3.),*/
	fDoTOFsigmaCut(kFALSE), 
	fPIDnSigmaAboveElectronLine(100), // sigma cut
	fPIDnSigmaBelowElectronLine(-100), // sigma cut
	fTofPIDnSigmaAboveElectronLine(100), // sigma cut RRnewTOF
	fTofPIDnSigmaBelowElectronLine(-100), // sigma cut RRnewTOF 
	fPIDnSigmaAbovePionLine(-100),     // sigma cut
	fPIDnSigmaAbovePionLineHighPt(-100),     // sigma cut
	fPIDMinPnSigmaAbovePionLine(100), // sigma cut
	fPIDMaxPnSigmaAbovePionLine(100), // sigma cut
	fDoKaonRejectionLowP(kFALSE),   // Kaon rejection at low p
	fDoProtonRejectionLowP(kFALSE), // Proton rejection at low p
	fDoPionRejectionLowP(kFALSE),   // Pion rejection at low p*/
	fPIDnSigmaAtLowPAroundKaonLine(0), // sigma cut
	fPIDnSigmaAtLowPAroundProtonLine(0), // sigma cut
	fPIDnSigmaAtLowPAroundPionLine(0), // sigma cut
	fPIDMinPKaonRejectionLowP(0), // Momentum limit to apply kaon rejection
	fPIDMinPProtonRejectionLowP(0), // Momentum limit to apply proton rejection
	fPIDMinPPionRejectionLowP(0), // Momentum limit to apply proton rejection
	fTrkSelectionCriteria(kGlobalTrack)
{
//
// Default constructor
//
	AdoptITSsaTrackCuts();
	AdoptESDtrackCuts();
	AdoptESDpidCuts();
	
	fGammaPool = new TClonesArray("AliKFParticle", fPoolMaxSize);
	fGammaPool->SetOwner(kTRUE);
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitz::AliAnalysisTaskGammaConvDalitz( const char* name ):
	AliAnalysisTaskSE( name ),
	fStack(0),
	fGCMCEvent(0),
	fESDEvent(0),
	fEposCandidateIndex(),
	fEnegCandidateIndex(),
	fGammaCandidatePosIndex(),
	fGammaCandidateNegIndex(),
	fGammaCandidates(0),
	fGammaPool(0),
	fPoolMaxSize(10),
	fGamPoolPos(0),
	fBGEventHandler(0),
	fOutputContainer(0),
	fMCTruth(0),
	fV0Reader(0),
	fESDpid(0),
	fESDtrackCuts(0),
	fITSsaTrackCuts(0),
	fESDpidCuts(0),
	fRemovePileUp(kFALSE),
	fSelectV0AND(kFALSE),                 // Select V0AND
	fTriggerAnalysis(0),
	fMultiplicity(0),
	fUseMultiplicity(0),
	fUseMultiplicityBin(0),
	fUseHBTMultiplicity(0),
	fUseHBTMultiplicityBin(0),
	fUseCentrality(0),
	fUseCentralityBin(0),
	fHistograms(0),
	fStandalone(kFALSE),
	fDoMC(kFALSE),
	fComputeBkg(kTRUE),
	fUseBayesPID(kFALSE),
	fUseTrackIndexCut(kTRUE),
	fUsePsiPairCut(kTRUE),
	fUseMassCut(kFALSE),
	fUseGammaCut(kFALSE),
	fUseCorrectedTPCClsInfo(kFALSE),
	fReadMagFieldSign(kTRUE),
	fUseAliKF(kFALSE),
	fBackgroundType(2),
	fMinClsTPCToF(0.),
	fMagFieldSign(1),
	fkElectronMass(0.00051099891),
	fPsiPairCut(0.45),
	fDeltaPhiCutMin(0.),
	fDeltaPhiCutMax(0.12),
	fMassCutMin(0.),
	fMassCutMax(0.1),
	/*fNSigmaBelowElecTPCbethe(-2.),
	fNSigmaAboveElecTPCbethe(3.),
	fNSigmaAbovePionTPCbethe(3.),
	fNSigmaAboveKaonTPCbethe(3.),
	fNSigmaAboveProtonTPCbethe(3.),*/
	fDoTOFsigmaCut(kFALSE), 
	fPIDnSigmaAboveElectronLine(100), // sigma cut
	fPIDnSigmaBelowElectronLine(-100), // sigma cut
	fTofPIDnSigmaAboveElectronLine(100), // sigma cut RRnewTOF
	fTofPIDnSigmaBelowElectronLine(-100), // sigma cut RRnewTOF 
	fPIDnSigmaAbovePionLine(-100),     // sigma cut
	fPIDnSigmaAbovePionLineHighPt(-100),     // sigma cut
	fPIDMinPnSigmaAbovePionLine(100), // sigma cut
	fPIDMaxPnSigmaAbovePionLine(100), // sigma cut
	fDoKaonRejectionLowP(kFALSE),   // Kaon rejection at low p
	fDoProtonRejectionLowP(kFALSE), // Proton rejection at low p
	fDoPionRejectionLowP(kFALSE),   // Pion rejection at low p*/
	fPIDnSigmaAtLowPAroundKaonLine(0), // sigma cut
	fPIDnSigmaAtLowPAroundProtonLine(0), // sigma cut
	fPIDnSigmaAtLowPAroundPionLine(0), // sigma cut
	fPIDMinPKaonRejectionLowP(0), // Momentum limit to apply kaon rejection
	fPIDMinPProtonRejectionLowP(0), // Momentum limit to apply proton rejection
	fPIDMinPPionRejectionLowP(0), // Momentum limit to apply proton rejection
	fTrkSelectionCriteria(kGlobalTrack)
{
	// Common I/O in slot 0
	DefineInput (0, TChain::Class());

	// Your private output
	DefineOutput(1, TList::Class());
//  DefineOutput(2, AliCFContainer::Class());  // for CF

	AdoptITSsaTrackCuts();
	AdoptESDtrackCuts();
	AdoptESDpidCuts();
// fkElectronMass = TDatabasePDG::Instance()->GetParticle( ::kElectron )->Mass(); //

	fGammaPool = new TClonesArray("AliKFParticle", fPoolMaxSize);
	fGammaPool->SetOwner(kTRUE);
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitz::~AliAnalysisTaskGammaConvDalitz()
{
//
// virtual destructor
//

	if( fOutputContainer )          delete fOutputContainer;
	if( fHistograms )               delete fHistograms;
	if( fStandalone && fV0Reader )  delete fV0Reader;
	if( fITSsaTrackCuts )           delete fITSsaTrackCuts;
	if( fESDtrackCuts )             delete fESDtrackCuts;
	if( fESDpidCuts )               delete fESDpidCuts;
	if( fGammaCandidates)           delete fGammaCandidates;
	if( fGammaPool )                delete fGammaPool;
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::ConnectInputData(Option_t *option)
{
//
// Connect Input Data
//
	if( fDebug ) AliInfo("=> ConnectInputData");

	AliAnalysisTaskSE::ConnectInputData(option);

	if( fV0Reader == 0 )
	{
		AliFatal("There is not pointer to AliV0Reader object!!!");
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::UserCreateOutputObjects()
{
//
// Create ouput objects
//
	if( fDebug ) AliInfo("=> UserCreateOutputObjects");

	// Create the output container
	if( fOutputContainer != 0 )
	{
		delete fOutputContainer;
	}

	fOutputContainer = new TList();

	// Add the histograms to the output container
	fHistograms->GetOutputContainer( fOutputContainer );
	fOutputContainer->SetOwner(kTRUE);

	PostData( 1, fOutputContainer );
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::UserExec(Option_t */*option*/)
{
//
// Execute analysis for current event
//        
	if( fDebug ) AliInfo("=> UserExec");

	if( fV0Reader == 0 )	{
		AliFatal("no pointer to AliV0Reader");
		return;
	}

	if(fTriggerAnalysis == 0){
		fTriggerAnalysis = new AliTriggerAnalysis;
	}

	// Create list of gamma candidates in standalone mode
	// otherwise use the created ones by AliAnalysisTaskGammaConversion
	if( fStandalone ) {
		AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
		AliESDInputHandler *esdHandler=0;
		if ( (esdHandler=dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler())) && esdHandler->GetESDpid() ){
			AliV0Reader::SetESDpid(esdHandler->GetESDpid());
		} else {
			//load esd pid bethe bloch parameters depending on the existance of the MC handler
			// yes: MC parameters
			// no:  data parameters
			if (!AliV0Reader::GetESDpid()){
				if (MCEvent() ) {
					AliV0Reader::InitESDpid();
				} else {
					AliV0Reader::InitESDpid(1);
				}
			}
		}
	
		if (MCEvent() ) {
		// To avoid crashes due to unzip errors. Sometimes the trees are not there.
			AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

			if (!mcHandler){ 
				AliError("Could not retrive MC event handler!"); 
				return; 
			}

			if (!mcHandler->InitOk() ){
				return;
			}
			if (!mcHandler->TreeK() ){
				return;
			}
			if (!mcHandler->TreeTR() ) {
				return;
			}
		}

		fV0Reader->SetInputAndMCEvent( InputEvent(), MCEvent() );
		fV0Reader->Initialize();
	}

	Bool_t v0A                   = fTriggerAnalysis->IsOfflineTriggerFired(fV0Reader->GetESDEvent(), AliTriggerAnalysis::kV0A);
	Bool_t v0C                   = fTriggerAnalysis->IsOfflineTriggerFired(fV0Reader->GetESDEvent(), AliTriggerAnalysis::kV0C);
	Bool_t v0AND = v0A && v0C;

	if( fTriggerAnalysis ){
		delete fTriggerAnalysis;
		fTriggerAnalysis=0;
	}


	if(fSelectV0AND && !v0AND){
		return;
	}

	if( fV0Reader->CheckForPrimaryVertex() == kFALSE )
	{
		if( fDebug ) AliInfo("no contributors to primary vertex");
		return;
	}

	if( fV0Reader->CheckForPrimaryVertexZ() == kFALSE  )
	{		
		if( fDebug ) AliInfo("z vertex out of range");
		return;
	}

	if(fRemovePileUp && fV0Reader->GetESDEvent()->IsPileupFromSPD()) {
		// eventQuality=4;	
		return;
	}	

	Int_t tTracklet=0, tITSTPC=0, tITSPure=0;
	fV0Reader->GetESDEvent()->EstimateMultiplicity(tTracklet, tITSTPC, tITSPure, 0.8);
	Int_t fMultiplicityITS = tITSPure;

	AliESDtrackCuts* kEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);


	kEsdTrackCuts->SetMaxDCAToVertexZ(2);
	kEsdTrackCuts->SetEtaRange(-0.8, 0.8);
	kEsdTrackCuts->SetPtRange(0.15);

	fMultiplicity = kEsdTrackCuts->CountAcceptedTracks(fV0Reader->GetESDEvent());
		
	if(  kEsdTrackCuts ) delete kEsdTrackCuts;

	if( fUseHBTMultiplicity==1) {
		fMultiplicity = fMultiplicityITS;
	}

	if(fUseMultiplicity!=0 && CalculateMultiplicityBin()!=fUseMultiplicityBin ){
		//eventQuality=6;
		//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}

			
	if(fUseHBTMultiplicity!=0 && CalculateMultiplicityBin()!=fUseHBTMultiplicityBin ){
		//eventQuality=6;
		//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
		return;
	}


	if(fV0Reader->GetIsHeavyIon()){
		if(fUseCentrality>0){
			AliCentrality *esdCentrality = fV0Reader->GetESDEvent()->GetCentrality();
			Int_t centralityC = -1;

			if(fUseCentrality==1){
				centralityC = esdCentrality->GetCentralityClass10("V0M");
				if( centralityC != fUseCentralityBin ){
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
			}

			if(fUseCentrality==2){
				centralityC = esdCentrality->GetCentralityClass10("CL1");
				if( centralityC != fUseCentralityBin ){
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
			}

			////////////////////////////////////// RRnew start /////////////////////////////////////////////////////
			if(fUseCentrality==3){
				centralityC = esdCentrality->GetCentralityClass10("V0M");
				if( (fUseCentralityBin == 0) && (centralityC!=0) ){ // 0-10%
						//eventQuality=7;
						//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
						return;         
				}
				if( (fUseCentralityBin == 1) && (centralityC!=1) ){ // 10-20%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 2) && (centralityC!=2) && (centralityC!=3) ){ // 20-40%
					// eventQuality=7;
				//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 3) && (centralityC!=0) && (centralityC!=1) ){ // 0-20%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 4) && (centralityC!=4) && (centralityC!=5) ){ // 40-60%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 6) && (centralityC!=6) && (centralityC!=7) && (centralityC!=8) ){ // 60-90%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 7) && (centralityC!=6) && (centralityC!=7) ){ // 60-80%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 8) && (centralityC>=8) ){ // 0-80%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 9) && (centralityC>=9) ){ // 0-90%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
			}

			if(fUseCentrality==4){
				centralityC = esdCentrality->GetCentralityClass10("CL1");
				if( (fUseCentralityBin == 0) && (centralityC!=0) ){ // 0-10%
				//eventQuality=7;
				//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
				return;         
				}
				if( (fUseCentralityBin == 1) && (centralityC!=1) ){ // 10-20%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 2) && (centralityC!=2) && (centralityC!=3) ){ // 20-40%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 4) && (centralityC!=4) && (centralityC!=5) ){ // 40-60%
					//eventQuality=7;
					//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
					return;
				}
				if( (fUseCentralityBin == 6) && (centralityC!=6) && (centralityC!=7) && (centralityC!=8) ){ // 60-90%
						//eventQuality=7;
						//fHistograms->FillHistogram("ESD_EventQuality",eventQuality);
						return;
				}
			}
		}
	}


	// Get Pointers
	fBGEventHandler = fV0Reader->GetBGHandler();
	fESDpid = fV0Reader->GetESDpid();
	fESDEvent = fV0Reader->GetESDEvent();
	if(fDoMC && MCEvent()){
		fStack= MCEvent()->Stack();
		fGCMCEvent=MCEvent();
	}

	// Read the magnetic field sign from ESD
	if ( fReadMagFieldSign == kTRUE ){
		fMagFieldSign = (fESDEvent->GetMagneticField() < 0) ? 1 : -1;
	}

	// Process MC information
	if(fDoMC){
		ProcessMCData();
	}

	if( fStandalone ){
		while(fV0Reader->NextV0()){}; //SelectGammas
		fV0Reader->ResetV0IndexNumber();
	}

	CreateListOfDalitzPairCandidates();
	ProcessGammaElectronsForDalitzAnalysis();
	
	if ( fStandalone ){
		fV0Reader->UpdateEventByEventData();
	}

	PostData( 1, fOutputContainer );
}


void AliAnalysisTaskGammaConvDalitz::Terminate(Option_t */*option*/)
{
//
    if( fDebug ) AliInfo("Not to do anything in Terminate");
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::AdoptITSsaTrackCuts( AliESDtrackCuts* esdCuts )
{
//
// set user ITSsa track cuts
//
	if( fITSsaTrackCuts ) delete fITSsaTrackCuts;

	if( esdCuts ){
		fITSsaTrackCuts = esdCuts;
	} else {
		// default cuts
		fITSsaTrackCuts = new AliESDtrackCuts("Default ITSsa track cuts for Pi0 Dalitz decay");	
		fITSsaTrackCuts->SetEtaRange( -0.9, 0.9 );
		fITSsaTrackCuts->SetAcceptKinkDaughters(kFALSE);

		fITSsaTrackCuts->SetMinNClustersITS(2);
		fITSsaTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
		fITSsaTrackCuts->SetRequireITSRefit(kTRUE);
		
		fITSsaTrackCuts->SetRequireSigmaToVertex(kTRUE);
		fITSsaTrackCuts->SetMaxNsigmaToVertex(3);
		
		fITSsaTrackCuts->SetRequireITSStandAlone(kTRUE);
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::AdoptESDtrackCuts( AliESDtrackCuts* esdCuts )
{
//
// set user global track cuts
//
	if( fESDtrackCuts ) delete fESDtrackCuts;

	if( esdCuts ){
		fESDtrackCuts = esdCuts;
	} else {
		//default cuts
		fESDtrackCuts = new AliESDtrackCuts("Default global track cuts for Pi0 Dalitz decay");

		fESDtrackCuts->SetEtaRange( -0.9, 0.9 );
		fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);

		fESDtrackCuts->SetMinNClustersITS(2);
		fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
		fESDtrackCuts->SetRequireITSRefit(kTRUE);

		fESDtrackCuts->SetRequireSigmaToVertex(kTRUE);
		fESDtrackCuts->SetMaxNsigmaToVertex(3);

		fESDtrackCuts->SetMinNClustersTPC(80);
		fESDtrackCuts->SetMaxChi2PerClusterTPC(4.);
		fESDtrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
		fESDtrackCuts->SetRequireTPCRefit(kTRUE);
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::AdoptESDpidCuts( AliESDpidCuts* esdPIDCuts )
{
//
// set user pid cuts
//
	if( fESDpidCuts ) delete fESDpidCuts;
	if( esdPIDCuts ) {
		fESDpidCuts = esdPIDCuts;
	} else { // default cuts
		fESDpidCuts = new AliESDpidCuts("Electrons", "Electron PID cuts");
	//   fESDpidCuts->SetTPCnSigmaCut(AliPID::kElectron, 3.);
		fESDpidCuts->SetTPCnSigmaCut(AliPID::kElectron, -4., 6.);
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::ProcessMCData()
{
//
// Process generation
//
	if( fDebug ) AliInfo("=> ProcessMCData");

	fHistograms->FillTable("Table_Generation", 0);  //number of events

	for ( Int_t i = 0; i < fStack->GetNprimary(); i++ ) {
		TParticle* iParticle = fStack->Particle( i );
		if( !iParticle ) continue;

		Float_t etaMin = 0.;
		Float_t etaMax = 0.;
		fESDtrackCuts->GetEtaRange(etaMin,etaMax);
		Double_t rapidity = Rapidity(iParticle);

		//_____________________________________JPsi ___________________________________________________
		if(iParticle->GetPdgCode() == 443){
			if(iParticle->GetNDaughters()==2){
				if(TMath::Abs(fStack->Particle(iParticle->GetFirstDaughter())->GetPdgCode()) == 11 &&
					TMath::Abs(fStack->Particle(iParticle->GetLastDaughter())->GetPdgCode()) == 11){
					if( TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ) {
						fHistograms->FillHistogram("MC_Jpsi_EE_Pt",iParticle->Pt());    
						
						TParticle* daug0 = fStack->Particle(iParticle->GetFirstDaughter());
						TParticle* daug1 = fStack->Particle(iParticle->GetLastDaughter());
						
						if( ( daug0->Eta() > etaMin && daug0->Eta() < etaMax ) && 
								( daug1->Eta() > etaMin && daug1->Eta() < etaMax ) ){
							fHistograms->FillHistogram("MC_Jpsi_EE_Acceptance_Pt_Vs_Y",iParticle->Pt(),rapidity);
						}
					}
				}
			}
		}

		//_____________________________________Chic_____________________________________________________
		if(iParticle->GetPdgCode() == 10441 || iParticle->GetPdgCode() == 20443 || iParticle->GetPdgCode() == 445){
			if( iParticle->GetNDaughters() == 2){
				TParticle* daug0 = fStack->Particle(iParticle->GetFirstDaughter());
				TParticle* daug1 = fStack->Particle(iParticle->GetLastDaughter());
												
				if( ( daug0->GetPdgCode() == 443 || daug0->GetPdgCode() == 22) && (daug1->GetPdgCode() == 443 || daug1->GetPdgCode() == 22) ){				
					if( TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ) {
						if( iParticle->GetPdgCode() == 10441 ) fHistograms->FillHistogram("MC_Chic0_All_Pt",iParticle->Pt() );
						if( iParticle->GetPdgCode() == 20443 ) fHistograms->FillHistogram("MC_Chic1_All_Pt",iParticle->Pt() );
						if( iParticle->GetPdgCode() ==   445 ) fHistograms->FillHistogram("MC_Chic2_All_Pt",iParticle->Pt() );
							
						TParticle* daugEPos;
						TParticle* daugENeg;
						TParticle* Gamma;
						
						Bool_t daugGammaEle = kFALSE;
						Bool_t daugGammaPos = kFALSE;
					
						if( daug0->GetPdgCode() == 443){
								daugEPos = fStack->Particle(daug0->GetFirstDaughter());
								daugENeg = fStack->Particle(daug0->GetLastDaughter());
								Gamma = daug1;				
						} else if (daug1->GetPdgCode() == 443){		
								daugEPos = fStack->Particle(daug1->GetFirstDaughter());
								daugENeg = fStack->Particle(daug1->GetLastDaughter());
								Gamma = daug0;		
						}
			
						if( Gamma->GetNDaughters() >= 2 ){
							for( Int_t tIndex=Gamma->GetFirstDaughter(); tIndex<=Gamma->GetLastDaughter(); ++tIndex ){
								TParticle* tmpDaughter = fStack->Particle(tIndex);
								if( tmpDaughter->GetUniqueID() != kPPair ) continue; // check if the daughters come from a conversion
								if( tmpDaughter->GetPdgCode() == ::kElectron ){ // e+
										// daugGammaElectronAll = kTRUE;
									if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && 
											TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
										if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
											if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1   ){
												if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
													daugGammaEle = kTRUE;
												}
											}
										}
									}
								} else if( tmpDaughter->GetPdgCode() == ::kPositron ) {
									//daugGammaPositronAll = kTRUE;
									if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
										if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
											if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1   ){
												if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
													daugGammaPos = kTRUE;
												}
											}
										}
									}
								}
							}
						}        


						Bool_t gammaChicEtaCut = kTRUE;
						if( TMath::Abs(Gamma->Eta()) > fV0Reader->GetEtaCut() || TMath::Abs(Gamma->Eta()) < fV0Reader->GetEtaCutMin() ) {
							gammaChicEtaCut = kFALSE;
						}

											
						if ( gammaChicEtaCut && Gamma->R() < fV0Reader->GetMaxRCut()  &&
							TMath::Abs( daugEPos->Eta() ) < etaMax  && TMath::Abs( daugENeg->Eta() ) < etaMax ){
							if( iParticle->GetPdgCode() == 10441 ) fHistograms->FillHistogram("MC_Chic0_Acceptance_Pt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
							if( iParticle->GetPdgCode() == 20443 ) fHistograms->FillHistogram("MC_Chic1_Acceptance_Pt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
							if( iParticle->GetPdgCode() ==   445 ) fHistograms->FillHistogram("MC_Chic2_Acceptance_Pt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
				
							if( daugGammaEle && daugGammaPos ){
									if( iParticle->GetPdgCode() == 10441 ) fHistograms->FillHistogram("MC_Chic0_Acceptance_ConvGammaPt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
									if( iParticle->GetPdgCode() == 20443 ) fHistograms->FillHistogram("MC_Chic1_Acceptance_ConvGammaPt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
									if( iParticle->GetPdgCode() ==   445 ) fHistograms->FillHistogram("MC_Chic2_Acceptance_ConvGammaPt_Vs_Y", iParticle->Pt(), Rapidity(iParticle) );
							}    
						}
					}
				}
			}
		}

		//_____________________________Pi0 or Eta______________________________________________________________
		if ( iParticle->GetPdgCode() != ::kPi0 && iParticle->GetPdgCode() != 221) continue; 
									
		if( iParticle->GetNDaughters() == 2 &&
			fStack->Particle(iParticle->GetFirstDaughter())->GetPdgCode() == ::kGamma &&
			fStack->Particle(iParticle->GetLastDaughter())->GetPdgCode() == ::kGamma ) {

			if( iParticle->GetPdgCode() == ::kPi0 ){
				if( TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ){
					fHistograms->FillTable("Table_Generation", 1);  // pi0 -> gg
				}
			}
		}
		if ( iParticle->GetNDaughters() != 3 ) continue;    // Num == 3 (e+,e-,gamma)
	
		// Check for Pi0 Dalitz decay
		TParticle* eposPi0 = 0;
		TParticle* enegPi0 = 0;
		TParticle* gammaPi0 = 0;
		
		for( Int_t idxPi0 = iParticle->GetFirstDaughter(); idxPi0 <= iParticle->GetLastDaughter(); idxPi0++ ){
			switch(fStack->Particle(idxPi0)->GetPdgCode()) {
				case ::kPositron:
					eposPi0 = fStack->Particle(idxPi0);
					break;
				case ::kElectron:
					enegPi0 = fStack->Particle(idxPi0);
					break;
				case ::kGamma:
					gammaPi0 = fStack->Particle(idxPi0);
					break;
			}
		}
	
		if (eposPi0==0 || enegPi0==0 || gammaPi0==0) continue;

		// check for gamma conversion
		Bool_t daugGammaElectron    = kFALSE;
		Bool_t daugGammaPositron    = kFALSE;  // acceptance
		Bool_t daugGammaElectronAll = kFALSE;
		Bool_t daugGammaPositronAll = kFALSE;


		Bool_t gammaEtaCut = kTRUE;
		if( TMath::Abs(gammaPi0->Eta()) > fV0Reader->GetEtaCut() || TMath::Abs(gammaPi0->Eta()) < fV0Reader->GetEtaCutMin() )
		gammaEtaCut = kFALSE;

		Bool_t gammaRCut = kTRUE;
		if( gammaPi0->R() > fV0Reader->GetMaxRCut()  ) 
		gammaRCut = kFALSE;

		// is the gamma converted? -> has 2 daughter e+e-
		// are e+ e- from gamma in the acceptance for the V0s

		if( gammaPi0->GetNDaughters() >= 2 ){
			for( Int_t tIndex=gammaPi0->GetFirstDaughter(); tIndex<=gammaPi0->GetLastDaughter(); ++tIndex ){
				TParticle* tmpDaughter = fStack->Particle(tIndex);
				if( tmpDaughter->GetUniqueID() != kPPair ) continue; // check if the daughters come from a conversion
			
				if( tmpDaughter->GetPdgCode() == ::kElectron ) { // e+
					daugGammaElectronAll = kTRUE;

					if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
						if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() ) ){
							if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1   ){
								if(tmpDaughter->R()< fV0Reader->GetMaxRCut()){
									daugGammaElectron = kTRUE;
								}
							}
						}
					}
				}
				else if( tmpDaughter->GetPdgCode() == ::kPositron ){
					daugGammaPositronAll = kTRUE;

					if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(tmpDaughter->Eta()) >= fV0Reader->GetEtaCutMin() ){
						if(  tmpDaughter->R() > ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue() )) {
							if ( (fV0Reader->GetEtaCutMin() != -0.1 &&  tmpDaughter->R() < ((TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlopeMin()) - fV0Reader->GetLineCutZValueMin())) || fV0Reader->GetEtaCutMin() == -0.1   ){
								if(tmpDaughter->R()< fV0Reader->GetMaxRCut()) {
										daugGammaPositron = kTRUE;
								}
							} 
						}
					}
				}
			}
		}


		if( iParticle->GetPdgCode() == ::kPi0 ) {					// found a Pi0 Dalitz decay
			if( TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ) {
				fHistograms->FillHistogram("MC_Pi0_All_Pt",iParticle->Pt());
				fHistograms->FillTable("Table_Generation", 2);
				fHistograms->FillHistogram("MC_Pi0Dalitz_P", iParticle->P());
				fHistograms->FillHistogram("MC_Pi0Dalitz_Pt", iParticle->Pt());

				fHistograms->FillHistogram("MC_Pi0Dalitz_Eta", iParticle->Eta());
				fHistograms->FillHistogram("MC_Pi0Dalitz_Pt_vs_Y", iParticle->Pt(),Rapidity(iParticle));
				fHistograms->FillHistogram("MC_EposDalitz_Pt", eposPi0->Pt());
				fHistograms->FillHistogram("MC_EposDalitz_Eta", eposPi0->Eta());
				fHistograms->FillHistogram("MC_EnegDalitz_Pt", enegPi0->Pt());
				fHistograms->FillHistogram("MC_EnegDalitz_Eta", enegPi0->Eta());
				fHistograms->FillHistogram("MC_GammaPi0Dalitz_Pt", gammaPi0->Pt());
				fHistograms->FillHistogram("MC_GammaPi0Dalitz_Eta", gammaPi0->Eta());
				
				// Angle between the gamma and the plane e+e-
				TVector3 ePosMom( eposPi0->Px(), eposPi0->Py(), eposPi0->Pz() );
				TVector3 eNegMom( enegPi0->Px(), enegPi0->Py(), enegPi0->Pz() );
				TVector3 gamMom( gammaPi0->Px(), gammaPi0->Py() , gammaPi0->Pz() );
				TVector3 planeEposEneg =  eNegMom.Cross( ePosMom );
				Double_t anglePlaneGamma = planeEposEneg.Angle(gamMom);
				
				fHistograms->FillHistogram("MC_EposEnegDalitz_Angle", ePosMom.Angle(eNegMom) );
				
				fHistograms->FillHistogram("MC_EposEnegDalitz_GammaPi0_Angle", anglePlaneGamma);
				fHistograms->FillHistogram("MC_EposEnegDalitz_GammaPi0_Angle_vs_P", anglePlaneGamma, gammaPi0->P());
				fHistograms->FillHistogram("MC_EposEnegDalitz_GammaPi0_Angle_vs_Pt", anglePlaneGamma, gammaPi0->Pt());
			
				if(  daugGammaElectronAll && daugGammaPositronAll ) {
					TParticle* tmpDaughter = fStack->Particle( gammaPi0->GetFirstDaughter() );
					fHistograms->FillHistogram("MC_GC_GammaPi0Dalitz_All_Z_vs_R",tmpDaughter->Vz(),tmpDaughter->R() );
				}

	
				if ( TMath::Abs( eposPi0->Eta() ) < etaMax  && TMath::Abs( enegPi0->Eta() ) < etaMax ) {
					fHistograms->FillHistogram("MC_Acceptance_EposDalitz_Pt", eposPi0->Pt());
					fHistograms->FillHistogram("MC_Acceptance_EposDalitz_Eta", eposPi0->Eta());
					fHistograms->FillHistogram("MC_Acceptance_EnegDalitz_Pt", enegPi0->Pt());
					fHistograms->FillHistogram("MC_Acceptance_EnegDalitz_Eta", enegPi0->Eta());
					fHistograms->FillHistogram("MC_Acceptance_DalitzPair_EposPt_vs_EnegPt", eposPi0->Pt(), enegPi0->Pt());
				}

			
				if ( gammaEtaCut == kTRUE && gammaRCut == kTRUE &&
					TMath::Abs(  eposPi0->Eta() ) < etaMax  && TMath::Abs( enegPi0->Eta() ) < etaMax ){

					fHistograms->FillTable("Table_Generation",3);  //

					fHistograms->FillHistogram("MC_Acceptance_Pi0Dalitz_Pt",iParticle->Pt());
					fHistograms->FillHistogram("MC_Acceptance_Pi0Dalitz_Eta",iParticle->Eta());
					fHistograms->FillHistogram("MC_Acceptance_Pi0Dalitz_Radius",iParticle->R());
					fHistograms->FillHistogram("MC_Acceptance_GammaPi0Dalitz_Pt",gammaPi0->Pt());
					fHistograms->FillHistogram("MC_Acceptance_GammaPi0Dalitz_Eta",gammaPi0->Eta());
					fHistograms->FillHistogram("MC_Acceptance_EposPi0Dalitz_Pt",eposPi0->Pt());
					fHistograms->FillHistogram("MC_Acceptance_EposPi0Dalitz_Eta",eposPi0->Eta());
					fHistograms->FillHistogram("MC_Acceptance_EnegPi0Dalitz_Pt",enegPi0->Pt());
					fHistograms->FillHistogram("MC_Acceptance_EnegPi0Dalitz_Eta",enegPi0->Eta());
					fHistograms->FillHistogram("MC_Acceptance_DalitzPair_OpeningAngle", ePosMom.Angle(eNegMom) );
					fHistograms->FillHistogram("MC_Pi0_Acceptance_Pt_vs_Y",iParticle->Pt(),Rapidity(iParticle));


					// Pi0 within acceptance with gamma converted

					if ( daugGammaElectron && daugGammaPositron ){
						fHistograms->FillTable("Table_Generation",4); //

						fHistograms->FillHistogram("MC_Acceptance_GC_Pi0Dalitz_Pt",iParticle->Pt());
						fHistograms->FillHistogram("MC_Acceptance_GC_Pi0Dalitz_Eta",iParticle->Eta());
						fHistograms->FillHistogram("MC_Acceptance_GC_EposPi0Dalitz_Pt",eposPi0->Pt());
						fHistograms->FillHistogram("MC_Acceptance_GC_EposPi0Dalitz_Eta",eposPi0->Eta());
						fHistograms->FillHistogram("MC_Acceptance_GC_EnegPi0Dalitz_Pt",enegPi0->Pt());
						fHistograms->FillHistogram("MC_Acceptance_GC_EnegPi0Dalitz_Eta",enegPi0->Eta());
						fHistograms->FillHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Pt",gammaPi0->Pt());
						fHistograms->FillHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Eta",gammaPi0->Eta());
						//fHistograms->FillHistogram("MC_Acceptance_GC_Gamma_Angle",anglePlaneGamma);
						//fHistograms->FillHistogram("MC_Acceptance_GC_Gamma_Angle_vs_Pt",anglePlaneGamma,gammaPi0->Pt());
						TParticle* tmpDaughter = fStack->Particle( gammaPi0->GetFirstDaughter() );
						fHistograms->FillHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Z_vs_R",tmpDaughter->Vz(),tmpDaughter->R() );
						fHistograms->FillHistogram("MC_Pi0_Acceptance_GammaConv_Pt_vs_Y",iParticle->Pt(), Rapidity(iParticle));
					}
				}
			}
		} else if( iParticle->GetPdgCode() == 221  ){ //	found Eta Dalitz decay

			if(  TMath::Abs(rapidity) < fV0Reader->GetRapidityMesonCut() ) {
			fHistograms->FillHistogram("MC_Eta_All_Pt",iParticle->Pt());

				if (  gammaEtaCut && gammaRCut &&
						TMath::Abs( eposPi0->Eta() ) < etaMax  && TMath::Abs( enegPi0->Eta() ) < etaMax  ){
					fHistograms->FillHistogram("MC_Eta_Acceptance_Pt_vs_Y",iParticle->Pt(),rapidity);

					if ( daugGammaElectron && daugGammaPositron ){
						fHistograms->FillHistogram("MC_Eta_Acceptance_GammaConv_Pt_vs_Y",iParticle->Pt(),rapidity);

					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::CreateListOfDalitzPairCandidates()
{
//
// Dalitz pair candidates
//
	if( fDebug ) AliInfo("=> CreateListOfDalitzPairCandidates");
	
	fEposCandidateIndex.clear();
	fEnegCandidateIndex.clear();

	fHistograms->FillTable("Table_Cuts", 0);
	
	for( Int_t i = 0; i < fESDEvent->GetNumberOfTracks(); ++i ){
		AliESDtrack* iTrack = fESDEvent->GetTrack(i);
		if ( !iTrack ) continue;

		Double_t p[3];
	
		if ( !iTrack->GetConstrainedPxPyPz(p) ) continue;

		TVector3 iMom(p[0],p[1],p[2]);

		//
		// Check track cuts and find track type
		//

		Bool_t isTrackAccepted = 0;
		Int_t trackType = -1;
		Double_t kClsToF = 0.;

		switch(fTrkSelectionCriteria){
			case kITSsaTrack:
				isTrackAccepted = fITSsaTrackCuts->AcceptTrack( iTrack );
				trackType = kITSsaTrack;
				break;
			
			case kGlobalTrack:
				isTrackAccepted = fESDtrackCuts->AcceptTrack( iTrack );
				trackType = kGlobalTrack;
				if( fUseCorrectedTPCClsInfo ){
					if( iTrack->GetTPCNclsF()!=0){
						kClsToF = (Double_t) ( iTrack->GetNcls(1)/iTrack->GetTPCNclsF() );
					}
				} else {
					kClsToF= iTrack->GetTPCClusterInfo(2,0,fV0Reader->GetFirstTPCRow(fV0Reader->GetXYRadius()));
				}
					
				if( kClsToF < fMinClsTPCToF ) {

								isTrackAccepted = kFALSE;
				}

				break;
			case kITSsaGlobalTrack:
				if(fITSsaTrackCuts->AcceptTrack( iTrack ) || fESDtrackCuts->AcceptTrack( iTrack ))
				{
					isTrackAccepted = kTRUE;
					if(fITSsaTrackCuts->AcceptTrack( iTrack )) trackType = kITSsaTrack;
					else trackType = kGlobalTrack;
				}
				break;
		}
		
		if(!isTrackAccepted) continue;
		
		//
		// PID
		//
		
		Int_t pid=-1;
		Int_t pidMC=-1;

		if(fUseBayesPID){
			pid = GetBayesPid(iTrack,trackType);
		} else {
			pid = GetNSigmaPid(iTrack,trackType);
		}
		
		if( fDoMC ){
			pidMC = GetMonteCarloPid(iTrack);
			// pid table
			Int_t iLabel = TMath::Abs(iTrack->GetLabel());
			TParticle* iParticle = fStack->Particle(iLabel);
			FillPidTable(iParticle, pid);
		}
		
		// ITS standalone tracks
		if( trackType == kITSsaTrack){
			Double_t mom = iTrack->GetP();
			Double_t signal = iTrack->GetITSsignal();
			
			fHistograms->FillHistogram( "ESD_ITSsa_dEdx_vs_P", mom, signal );
			
			if( pid == AliPID::kElectron ) {
				fHistograms->FillHistogram( "ESD_ITSsa_PidCut_dEdx_vs_P", mom, signal );
				if(fDoMC && pid == pidMC) {
					fHistograms->FillHistogram( "MC_ESD_ITSsa_PidCut_dEdx_vs_P", mom, signal );
				}
			}
			
			if( fDoMC && pidMC == AliPID::kElectron){
				fHistograms->FillHistogram( "MC_ESD_ITSsa_Electron_dEdx_vs_P", mom, signal );
			}
		}
		
		else { // global tracks
			const AliExternalTrackParam *in = iTrack->GetInnerParam();
			Double_t mom = in->GetP();
			Double_t signal = iTrack->GetTPCsignal();
			
			fHistograms->FillHistogram( "ESD_TPC_dEdx_vs_P", mom, signal );

			UInt_t statusTrack = iTrack->GetStatus(); //moved up here from below RRnewTOF
	
			Double_t timesTrack[5];
			iTrack->GetIntegratedTimes(timesTrack);
	
			Double_t t0pos = fV0Reader->GetESDpid()->GetTOFResponse().GetStartTime(iTrack->GetP());
	
		
			
			Double_t TOFsignaTrack = iTrack->GetTOFsignal();
			

			Double_t dTtrack = TOFsignaTrack - t0pos - timesTrack[0];
			

			if( (statusTrack & AliESDtrack::kTOFpid) && !(statusTrack & AliESDtrack::kTOFmismatch) )
				fHistograms->FillHistogram("ESD_TOF_E_P_dT", iTrack->GetP(), dTtrack); //if( (statusNeg & AliESDtrack::kTOFpid) && !(statusNeg & AliESDtrack::kTOFmismatch) ) fHistograms->FillHistogram("ESD_ConvGamma_EandP_P_dT", fV0Reader->GetNegativeTrackP(), dTneg);

			if( fDoMC && pidMC == AliPID::kElectron ){
				fHistograms->FillHistogram( "MC_ESD_TPC_Electron_dEdx_vs_P", mom, signal );
				if( (statusTrack & AliESDtrack::kTOFpid) && !(statusTrack & AliESDtrack::kTOFmismatch) ){
					fHistograms->FillHistogram("MC_ESD_TOF_E_P_dT", iTrack->GetP(), dTtrack);
				}
			}
			
			if( pid == AliPID::kElectron ){
				fHistograms->FillHistogram("ESD_TPC_PidCut_dEdx_vs_P", mom, signal );
				if( (statusTrack & AliESDtrack::kTOFpid) && !(statusTrack & AliESDtrack::kTOFmismatch) ){
					fHistograms->FillHistogram("ESD_TOF_PidCut_E_P_dT", iTrack->GetP(), dTtrack); 
				}
				if(fDoMC && pid == pidMC) {
					fHistograms->FillHistogram( "MC_ESD_TPC_PidCut_dEdx_vs_P", mom, signal );
					if( (statusTrack & AliESDtrack::kTOFpid) && !(statusTrack & AliESDtrack::kTOFmismatch) ){
						fHistograms->FillHistogram("MC_ESD_TOF_PidCut_E_P_dT", iTrack->GetP(), dTtrack); 	
					}
				}
			}
		}
		
		if( AliPID::kElectron != pid) continue;
		
		// electron track candidates from here
		
		if( iTrack->GetSign() > 0 ){
			fEposCandidateIndex.push_back(i);
		} else {
			fEnegCandidateIndex.push_back(i);
		}
	}
	
	// gamma candidates
	GetGammaCandidates(fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex);
	
	if(fDoMC) {
		TClonesArray* pi0Dalitz = FindParticleDalitz(fEposCandidateIndex, fEnegCandidateIndex, fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex,1);
		((TH1*)fHistograms->GetValue("Table_Cuts"))->Fill(1.,(Double_t)pi0Dalitz->GetEntriesFast());
		delete pi0Dalitz;
	}
	
	if(fUseTrackIndexCut) {// remove repeated tracks
		ESDtrackIndexCut(fEposCandidateIndex,fEnegCandidateIndex, fGammaCandidates);
		if(fDoMC) {
			TClonesArray* pi0Dalitz = FindParticleDalitz(fEposCandidateIndex, fEnegCandidateIndex, fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex,1);
			((TH1*)fHistograms->GetValue("Table_Cuts"))->Fill(2.,(Double_t)pi0Dalitz->GetEntriesFast());
			delete pi0Dalitz;
		}
	}
	
	if(fUsePsiPairCut) {// remove electrons from gamma conversions
		PsiPairCut(fEposCandidateIndex,fEnegCandidateIndex);
		if(fDoMC) {
			TClonesArray* pi0Dalitz = FindParticleDalitz(fEposCandidateIndex, fEnegCandidateIndex, fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex,1);
			((TH1*)fHistograms->GetValue("Table_Cuts"))->Fill(3.,(Double_t)pi0Dalitz->GetEntriesFast());
			delete pi0Dalitz;
		}
	}
	
	if( fUseMassCut ) {
		MassCut(fEposCandidateIndex, fEnegCandidateIndex);
		if(fDoMC) {
			TClonesArray* pi0Dalitz = FindParticleDalitz(fEposCandidateIndex, fEnegCandidateIndex, fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex,1);
			((TH1*)fHistograms->GetValue("Table_Cuts"))->Fill(4.,(Double_t)pi0Dalitz->GetEntriesFast());
			delete pi0Dalitz;
		}
	}
	
	if(fUseGammaCut) {
		AngleEposEnegGammaCut(fEposCandidateIndex,fEnegCandidateIndex,fV0Reader->GetCurrentEventGoodV0s(), fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex);
		
		if(fDoMC) {
			TClonesArray* pi0Dalitz = FindParticleDalitz(fEposCandidateIndex, fEnegCandidateIndex, fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex,1);
			((TH1*)fHistograms->GetValue("Table_Cuts"))->Fill(5.,(Double_t)pi0Dalitz->GetEntriesFast());
			delete pi0Dalitz;
		}
	}
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::ProcessGammaElectronsForDalitzAnalysis()
{
//
// Process gamma and electrons for pi0 Dalitz decay
//
	if( fDebug ) AliInfo("=> ProcessGammaElectronsForDalitzAnalysis");
	
	fHistograms->FillTable( "Table_Reconstruction", 0); // number of events

	TClonesArray* ePosCandidates = IndexToAliKFParticle(fEposCandidateIndex, ::kPositron);
	
	for(Int_t i=0; i < ePosCandidates->GetEntriesFast(); ++i) {
		AliKFParticle* epos = (AliKFParticle*) ePosCandidates->At(i);
		fHistograms->FillHistogram("ESD_EposCandidates_Pt", epos->GetPt());
		fHistograms->FillHistogram("ESD_EposCandidates_Eta", epos->GetEta());
		fHistograms->FillTable( "Table_Reconstruction", 1);
	}
	
	TClonesArray* eNegCandidates = IndexToAliKFParticle(fEnegCandidateIndex, ::kElectron);
	
	for(Int_t i=0; i < eNegCandidates->GetEntriesFast(); ++i) {
		AliKFParticle* eneg = (AliKFParticle*) eNegCandidates->At(i);
		fHistograms->FillHistogram("ESD_EnegCandidates_Pt", eneg->GetPt());
		fHistograms->FillHistogram("ESD_EnegCandidates_Eta", eneg->GetEta());
		fHistograms->FillTable( "Table_Reconstruction", 2);
	}
	
	TClonesArray* dalitzPairCandidates = FindDalitzPair(ePosCandidates, eNegCandidates);
	for(Int_t i=0; i < dalitzPairCandidates->GetEntriesFast(); ++i) {
		TLorentzVector* dalitz = (TLorentzVector*)dalitzPairCandidates->At(i);
		fHistograms->FillHistogram("ESD_DalitzPairCandidates_Pt", dalitz->Pt());
		fHistograms->FillHistogram("ESD_DalitzPairCandidates_InvMass", dalitz->M());
      fHistograms->FillHistogram("ESD_DalitzPairCandidates_InvMass_vs_Pt",dalitz->M(),dalitz->Pt());
	}
	
	// gamma candidates
	for(Int_t i=0; i < fGammaCandidates->GetEntriesFast(); ++i) {
		AliKFParticle* gamma = (AliKFParticle*) fGammaCandidates->At(i);
		fHistograms->FillHistogram("ESD_GammaCandidates_Pt", gamma->GetPt());
		fHistograms->FillHistogram("ESD_GammaCandidates_Eta", gamma->GetEta());
	}
	
	// psi pair for all candidates
	//if(fUsePsiPairCut)
	FillPsiPair(ePosCandidates,eNegCandidates,"ESD_EposEneg_PsiPair_vs_DPhi");
	
	// Angle epos,eneg gamma
	FillAngle(ePosCandidates, fGammaCandidates, "ESD_EposEneg_GammaCandidates_Angle");
	FillAngle(eNegCandidates, fGammaCandidates, "ESD_EposEneg_GammaCandidates_Angle");
	
	
	for(Int_t iPos=0; iPos < ePosCandidates->GetEntriesFast(); ++iPos) {
		AliKFParticle* lPosKF = (AliKFParticle*)ePosCandidates->At(iPos);

		for(Int_t iNeg=0; iNeg < eNegCandidates->GetEntriesFast(); ++iNeg) {
			AliKFParticle* lNegKF = (AliKFParticle*)eNegCandidates->At(iNeg);
			AliKFParticle lPosNeg(*lPosKF,*lNegKF );
			
			for(Int_t iGam=0; iGam < fGammaCandidates->GetEntriesFast(); ++iGam) {
				AliKFParticle* lGamKF = (AliKFParticle*)fGammaCandidates->At(iGam);
				
				AliKFParticle lPosNegGam( *lPosKF, *lNegKF, *lGamKF );

				Double_t lDiffMass = lPosNegGam.GetMass() - lPosNeg.GetMass();

				Double_t RapidityParticle     = RapidityKF(lPosNegGam); 
				if (	TMath::Abs(RapidityParticle) < fV0Reader->GetRapidityMesonCut()){ 
					fHistograms->FillHistogram("ESD_EposEnegGamma_P", lPosNegGam.GetP());
					fHistograms->FillHistogram("ESD_EposEnegGamma_Eta",lPosNegGam.GetEta());
					fHistograms->FillHistogram("ESD_EposEnegGamma_Phi",lPosNegGam.GetPhi()); 
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMass",lPosNegGam.GetMass());
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMass_vs_Eta",lPosNegGam.GetMass(),lPosNegGam.GetEta());
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMass_vs_Y",lPosNegGam.GetMass(),RapidityParticle);
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMass_vs_Pt",lPosNegGam.GetMass(),lPosNegGam.GetPt());
					fHistograms->FillHistogram("ESD_EposEnegGamma_Pt_vs_Y",lPosNegGam.GetPt(),RapidityParticle);
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMassDiff",lDiffMass );
					fHistograms->FillHistogram("ESD_EposEnegGamma_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
				}
			}
		}
	}

	delete dalitzPairCandidates;
	//delete pi0Candidates;
	
	if(fComputeBkg) {

		// 1) e+e- dalitz
		for(Int_t i=0; i < ePosCandidates->GetEntriesFast(); ++i) {
			AliKFParticle* epos1 = (AliKFParticle*) ePosCandidates->At(i);

			for(Int_t j=i+1; j < ePosCandidates->GetEntriesFast(); ++j) {
				AliKFParticle* epos2 = (AliKFParticle*) ePosCandidates->At(j);
				AliKFParticle ePosePos( *epos1,*epos2 );
				fHistograms->FillHistogram("ESD_BKG_LikeSign_InvMass",ePosePos.GetMass());
				fHistograms->FillHistogram("ESD_BKG_LikeSign_InvMass_vs_Pt",ePosePos.GetMass(),ePosePos.GetPt());
			}
		}
		for(Int_t i=0; i < eNegCandidates->GetEntriesFast(); ++i) {
			AliKFParticle* eneg1 = (AliKFParticle*) eNegCandidates->At(i);

			for(Int_t j=i+1; j < eNegCandidates->GetEntriesFast(); ++j) {
				AliKFParticle* eneg2 = (AliKFParticle*) eNegCandidates->At(j);
				AliKFParticle eNegeNeg( *eneg1,*eneg2 );
				fHistograms->FillHistogram("ESD_BKG_LikeSign_InvMass",eNegeNeg.GetMass());
				fHistograms->FillHistogram("ESD_BKG_LikeSign_InvMass_vs_Pt",eNegeNeg.GetMass(),eNegeNeg.GetPt());
			}
		}

		TClonesArray* pi0Bkg;
		TClonesArray* gammaBGHandler;
		TClonesArray* elecBGHandler;

		if( fBackgroundType == 1 ) {

			// 1) e+e- with with gammas used in the signal
			pi0Bkg = FindParticleDalitz(ePosCandidates, eNegCandidates, fGammaPool,1);
		
			for(Int_t i=0; i < pi0Bkg->GetEntriesFast(); ++i) {
				TLorentzVector* pi0 = (TLorentzVector*)pi0Bkg->At(i);
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass", pi0->M());
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass_vs_Pt",pi0->M(),pi0->Pt());
			}
			
			if(ePosCandidates->GetEntriesFast() > 0 &&
				eNegCandidates->GetEntriesFast() > 0 &&
				fGammaCandidates->GetEntriesFast() > 0) {
				
				UpdateGammaPool(fGammaCandidates);
			}
			
			delete pi0Bkg;
		} else if( fBackgroundType == 2) { // 2) e+e- with gammas from a pool of events
			gammaBGHandler = GammasFromBGHandler();
			pi0Bkg = FindParticleDalitz(ePosCandidates, eNegCandidates, gammaBGHandler,2);
			
			for(Int_t i=0; i < pi0Bkg->GetEntriesFast(); ++i) {
				TLorentzVector* pi0 = (TLorentzVector*)pi0Bkg->At(i);
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass", pi0->M());
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass_vs_Pt",pi0->M(), pi0->Pt());
			}
		
			delete pi0Bkg;
         delete gammaBGHandler;
		}
		
		// 3) e+ with e-, gamma from a pool of events
		else if( fBackgroundType == 3 ){

			elecBGHandler = ElectronFromBGHandler();
			gammaBGHandler = GammasFromBGHandler();
			pi0Bkg = FindParticleDalitz(ePosCandidates, elecBGHandler, gammaBGHandler,3);
		
			for(Int_t i=0; i < pi0Bkg->GetEntriesFast(); ++i){
				TLorentzVector* pi0 = (TLorentzVector*)pi0Bkg->At(i);
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass", pi0->M());
				fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMass_vs_Pt",pi0->M(), pi0->Pt());
			}
		
			if(eNegCandidates->GetEntriesFast() > 0) {
				UpdateElectronPool(eNegCandidates);
			}
		
			delete gammaBGHandler;
			delete elecBGHandler;
			delete pi0Bkg;
		}
	}
	
	delete ePosCandidates;
	delete eNegCandidates;
	
	if(fDoMC) {
		TClonesArray* ePosPi0Dalitz = FindElectronFromPi0Dalitz(fEposCandidateIndex, ::kPositron);
		for(Int_t i=0; i < ePosPi0Dalitz->GetEntriesFast(); ++i) {
			AliKFParticle* epos = (AliKFParticle*) ePosPi0Dalitz->At(i);

			fHistograms->FillHistogram("MC_ESD_Pi0_EposDalitz_Pt", epos->GetPt());
			fHistograms->FillHistogram("MC_ESD_Pi0_EposDalitz_Eta", epos->GetEta());

			fHistograms->FillTable( "Table_Reconstruction", 3);
		}
		
		TClonesArray* eNegPi0Dalitz = FindElectronFromPi0Dalitz(fEnegCandidateIndex, ::kElectron);
		for(Int_t i=0; i < eNegPi0Dalitz->GetEntriesFast(); ++i){
			AliKFParticle* eneg = (AliKFParticle*) eNegPi0Dalitz->At(i);
			fHistograms->FillHistogram("MC_ESD_Pi0_EnegDalitz_Pt", eneg->GetPt());
			fHistograms->FillHistogram("MC_ESD_Pi0_EnegDalitz_Eta", eneg->GetEta());
			fHistograms->FillTable( "Table_Reconstruction", 4);
		}
		
		TClonesArray* dalitzPairPi0 = FindDalitzPair(fEposCandidateIndex, fEnegCandidateIndex,1);
		for(Int_t i=0; i < dalitzPairPi0->GetEntriesFast(); ++i) {
			TLorentzVector* dalitz = (TLorentzVector*) dalitzPairPi0->At(i);
			fHistograms->FillHistogram("MC_ESD_Pi0_DalitzPair_Pt", dalitz->Pt());
			fHistograms->FillHistogram("MC_ESD_Pi0_DalitzPair_Mass", dalitz->M());
			fHistograms->FillHistogram( "Table_Reconstruction", 5 );
		}
		
		TClonesArray* dalitzPairEta = FindDalitzPair(fEposCandidateIndex, fEnegCandidateIndex,2);
		for(Int_t i=0; i < dalitzPairEta->GetEntriesFast(); ++i) {
			TLorentzVector* dalitz = (TLorentzVector*) dalitzPairEta->At(i);
			fHistograms->FillHistogram("MC_ESD_Eta0_DalitzPair_Pt", dalitz->Pt());
			fHistograms->FillHistogram("MC_ESD_Eta0_DalitzPair_InvMass", dalitz->M());
		}

		TClonesArray* lJpsiAll = FindJpsi(fEposCandidateIndex, fEnegCandidateIndex,-1);

		for(Int_t i=0; i < lJpsiAll->GetEntriesFast(); ++i){
				TLorentzVector* jpsi = (TLorentzVector*) lJpsiAll->At(i);
				fHistograms->FillHistogram("MC_ESD_Jpsi_Pt",jpsi->Pt());
				fHistograms->FillHistogram("MC_ESD_Jpsi_InvMass",jpsi->M());
				fHistograms->FillHistogram("MC_ESD_Jpsi_InvMass_vs_Pt",jpsi->M(),jpsi->Pt());
		}


		TClonesArray* lJpsiChic0 = FindJpsi(fEposCandidateIndex, fEnegCandidateIndex,0);

		for(Int_t i=0; i < lJpsiChic0->GetEntriesFast(); ++i) {
				TLorentzVector* jpsi = (TLorentzVector*) lJpsiChic0->At(i);
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic0_Pt",jpsi->Pt());
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic0_InvMass",jpsi->M());
		}
		TClonesArray* lJpsiChic1 = FindJpsi(fEposCandidateIndex, fEnegCandidateIndex,1);

		for(Int_t i=0; i < lJpsiChic1->GetEntriesFast(); ++i) {
				TLorentzVector* jpsi = (TLorentzVector*) lJpsiChic1->At(i);
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic1_Pt",jpsi->Pt());
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic1_InvMass",jpsi->M());
		}
		TClonesArray* lJpsiChic2 = FindJpsi(fEposCandidateIndex, fEnegCandidateIndex,2);

		for(Int_t i=0; i < lJpsiChic2->GetEntriesFast(); ++i) {
				TLorentzVector* jpsi = (TLorentzVector*) lJpsiChic2->At(i);
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic2_Pt",jpsi->Pt());
				fHistograms->FillHistogram("MC_ESD_Jpsi_Chic2_InvMass",jpsi->M());
		}

		// psi pair for dalitz pairs
		//if(fUsePsiPairCut)
		FillPsiPair(ePosPi0Dalitz,eNegPi0Dalitz,"MC_ESD_Pi0_DalitzPair_PsiPair_vs_DPhi");
		
		delete ePosPi0Dalitz;
		delete eNegPi0Dalitz;
		delete dalitzPairPi0;
		delete dalitzPairEta;
		delete lJpsiAll;
		delete lJpsiChic0;
		delete lJpsiChic1;
		delete lJpsiChic2;
		// all gammas
		TClonesArray* gamma = FindGamma(fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex);
		for(Int_t i=0; i < gamma->GetEntriesFast(); ++i){
			AliKFParticle* iGamma = (AliKFParticle*) gamma->At(i);
			fHistograms->FillHistogram("MC_ESD_Gamma_Pt", iGamma->GetPt());
			fHistograms->FillHistogram("MC_ESD_Gamma_Eta", iGamma->GetEta());
		}
		
		delete gamma;
		
		// gamma from pi0 dalitz
		TClonesArray* gammaPi0Dalitz = FindGammaFromPi0Dalitz(fGammaCandidates, fGammaCandidatePosIndex, fGammaCandidateNegIndex);
		for(Int_t i=0; i < gammaPi0Dalitz->GetEntriesFast(); ++i)
		{
			AliKFParticle* iGamma = (AliKFParticle*) gammaPi0Dalitz->At(i);
			fHistograms->FillHistogram("MC_ESD_GammaPi0Dalitz_Pt", iGamma->GetPt());
			fHistograms->FillHistogram("MC_ESD_GammaPi0Dalitz_Eta", iGamma->GetEta());
			fHistograms->FillTable( "Table_Reconstruction", 6);
		}
		
		delete gammaPi0Dalitz;
		
		for( UInt_t iPos=0; iPos < fEposCandidateIndex.size(); ++iPos ) {
			AliESDtrack* posTrack = fESDEvent->GetTrack(fEposCandidateIndex[iPos]);
			Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
			Int_t posLabel = TMath::Abs(posTrack->GetLabel());
							
			AliKFParticle posKF( *posTrack->GetConstrainedParam(), ::kPositron );
		
			for( UInt_t iNeg=0; iNeg < fEnegCandidateIndex.size(); ++iNeg ) {
				AliESDtrack* negTrack = fESDEvent->GetTrack(fEnegCandidateIndex[iNeg]);
				Int_t negLabel = TMath::Abs(negTrack->GetLabel());

				Int_t jpsiLabel = fStack->Particle(negLabel)->GetMother(0);

				if( !HaveSameMother(posLabel,negLabel) ) continue; //Check if both particles have same mother

				Bool_t kPi0Mother  = IsDalitzPair(posLabel,negLabel,1); //comes from Pi0
				Bool_t kEta0Mother = IsDalitzPair(posLabel,negLabel,2); //comes from Eta0
				Bool_t kJpsiMother  = kFALSE;
				Bool_t kGammaMother = kFALSE;


				if ( fStack->Particle(jpsiLabel)->GetPdgCode() == 443  ) {
					kJpsiMother = kTRUE;
				} else if( fStack->Particle(jpsiLabel)->GetPdgCode() ==  ::kGamma   ) {
					kGammaMother = kTRUE;
				}

				
				if( !kPi0Mother  &&  //1: Pi0->DalitzDecay 
						!kEta0Mother &&  //2: Eta0->DalitzDecay
						!kJpsiMother &&  //3: Jpsi
						!kGammaMother   ) continue;
				
				Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
				AliKFParticle negKF( *negTrack->GetConstrainedParam(), ::kElectron );
				AliKFParticle lPosNeg(posKF,negKF );
				
				for(Int_t iGam=0; iGam < fGammaCandidates->GetEntriesFast(); ++iGam) {
					AliKFParticle* gamKF = (AliKFParticle*)fGammaCandidates->At(iGam);
					
					Int_t labelv1 = TMath::Abs((fESDEvent->GetTrack(fGammaCandidatePosIndex[iGam]))->GetLabel());
					Int_t labelv2 = TMath::Abs((fESDEvent->GetTrack(fGammaCandidateNegIndex[iGam]))->GetLabel());
					
					if( !HaveSameMother(labelv1,labelv2) ) continue;
					
					Int_t labelGamma = TMath::Abs(fStack->Particle(labelv1)->GetMother(0));
					
					if( fStack->Particle(labelGamma)->GetPdgCode() != ::kGamma ) continue;

					AliKFParticle lPosNegGam( posKF, negKF, *gamKF );

					Double_t lDiffMass = lPosNegGam.GetMass() - lPosNeg.GetMass();
					Double_t RapidityParticle     = RapidityKF(lPosNegGam); 

					if (	TMath::Abs(RapidityParticle) < fV0Reader->GetRapidityMesonCut()){ 
						if( kPi0Mother == kTRUE || kEta0Mother == kTRUE ){
							if( !HaveSameMother(labelGamma, posLabel) ) continue;
						
							if( IsDalitzPair(posLabel,negLabel,1) == kTRUE ){
								fHistograms->FillHistogram("MC_ESD_Pi0_P", lPosNegGam.GetP());
								fHistograms->FillHistogram("MC_ESD_Pi0_Phi", lPosNegGam.GetPhi());
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMass", lPosNegGam.GetMass());
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMass_vs_Pt",lPosNegGam.GetMass(),lPosNegGam.GetPt());
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMass_vs_Y", lPosNegGam.GetMass(), RapidityParticle);
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMass_vs_Eta", lPosNegGam.GetMass(),lPosNegGam.GetEta());
								fHistograms->FillHistogram("MC_ESD_Pi0_Pt_vs_Y",lPosNegGam.GetPt(), RapidityParticle);
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMassDiff",lDiffMass);
								fHistograms->FillHistogram("MC_ESD_Pi0_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
							} else if( IsDalitzPair(posLabel,negLabel,2) == kTRUE ) {
								fHistograms->FillHistogram("MC_ESD_Eta0_P", lPosNegGam.GetP());
								fHistograms->FillHistogram("MC_ESD_Eta0_Phi", lPosNegGam.GetPhi());
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMass",lPosNegGam.GetMass());
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMass_vs_Pt", lPosNegGam.GetMass(),lPosNegGam.GetPt());
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMass_vs_Y",lPosNegGam.GetMass(),RapidityParticle);
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMass_vs_Eta",lPosNegGam.GetMass(),lPosNegGam.GetEta());
								fHistograms->FillHistogram("MC_ESD_Eta0_Pt_vs_Y", lPosNegGam.GetPt(),RapidityParticle);
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMassDiff",lDiffMass);
								fHistograms->FillHistogram("MC_ESD_Eta0_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
							}
							fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMassDiff",lDiffMass);
							fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
							fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMass",lPosNegGam.GetMass());
							fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMass_vs_Pt",lPosNegGam.GetMass(),lPosNegGam.GetPt());
						} else if( kGammaMother == kTRUE ) {
							if( !HaveSameMother(labelGamma, jpsiLabel ) ) continue;
							if( fStack->Particle(labelGamma)->GetMother(0) < 0 ) continue;

							Int_t MotherGammaGammalabel = fStack->Particle(labelGamma)->GetMother(0);
							TParticle* parMotherGammaGamma = fStack->Particle(MotherGammaGammalabel);
							

							if( parMotherGammaGamma->GetPdgCode() == ::kPi0 || parMotherGammaGamma->GetPdgCode() == 221 ){
								fHistograms->FillHistogram("MC_ESD_GammaGamma_InvMassDiff",lDiffMass);
								fHistograms->FillHistogram("MC_ESD_GammaGamma_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
								fHistograms->FillHistogram("MC_ESD_GammaGamma_InvMass",lPosNegGam.GetMass());
								fHistograms->FillHistogram("MC_ESD_GammaGamma_InvMass_vs_Pt",lPosNegGam.GetMass(),lPosNegGam.GetPt());
							}
						} else if( kJpsiMother == kTRUE ){
							if( !HaveSameMother(labelGamma, jpsiLabel ) ) continue;


							TParticle *jpsiParticle = fStack->Particle(jpsiLabel);

							if ( jpsiParticle->GetMother(0) < 0 ) continue; 


							Int_t chicLabel = jpsiParticle->GetMother(0);

							Int_t MotherpdgCode = fStack->Particle(chicLabel)->GetPdgCode();


							if( MotherpdgCode == 10441 || MotherpdgCode == 20443 || MotherpdgCode == 445 ) {
								fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMassDiff",lDiffMass);
								fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMassDiff_vs_Pt",lDiffMass,lPosNegGam.GetPt());
								fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMass",lPosNegGam.GetMass());
								fHistograms->FillHistogram("MC_ESD_EposEnegGamma_InvMass_vs_Pt",lPosNegGam.GetMass(),lPosNegGam.GetPt());
							}
						}
					}
				} 
			}
		}
		
		// psi pair for electrons from gamma conversions assuming they came from main vertex
		// if(fUsePsiPairCut)
		for(UInt_t i=0; i < fEposCandidateIndex.size(); ++i) {
			AliESDtrack* posTrack = fESDEvent->GetTrack(fEposCandidateIndex[i]);
			Int_t posLabel = TMath::Abs(posTrack->GetLabel());
			
			for(UInt_t j=0; j < fEnegCandidateIndex.size(); ++j) {
				AliESDtrack* negTrack = fESDEvent->GetTrack(fEnegCandidateIndex[j]);
				Int_t negLabel = TMath::Abs(negTrack->GetLabel());
				
				if(!IsFromGammaConversion(posLabel,negLabel)) continue;
				
				Double_t psiPair = GetPsiPair(posTrack, negTrack);
				Double_t deltaPhi = fMagFieldSign * TVector2::Phi_mpi_pi( negTrack->GetConstrainedParam()->Phi()-posTrack->GetConstrainedParam()->Phi());
				
				fHistograms->FillHistogram("MC_ESD_EposEnegGamma_PsiPair_vs_DPhi", deltaPhi, psiPair);
			}
		}
		// FIXME: eta -> e+e-gamma
	}
}

//--------------------------------------------------------------------------
Double_t AliAnalysisTaskGammaConvDalitz::Rapidity(const TParticle* p) const
{
//
// Get rapidity
//
	const double kEPSILON=1.e-16;
	
	if(p->Energy() - TMath::Abs(p->Pz()) < kEPSILON ) {
		return 1.e10;
	}
	return 0.5*TMath::Log( (p->Energy()+p->Pz()) / (p->Energy()-p->Pz()) );
}

Double_t AliAnalysisTaskGammaConvDalitz::RapidityKF(const AliKFParticle p) const
{
//
// Get rapidity
//
	const double kEPSILON=1.e-16;
	
	if((p.GetE()-TMath::Abs(p.GetPz())) < kEPSILON ) {
		return 1.e10;
	}
	return 0.5*TMath::Log( (p.GetE()+p.GetPz()) / (p.GetE()-p.GetPz()) );
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::FillPsiPair(const TClonesArray* pos, const TClonesArray* neg, const TString& hName)
{
//
// Fill histogram with psipair(pos,neg)
//
	for(Int_t i=0; i < pos->GetEntriesFast(); ++i ){
		AliKFParticle* posKF = (AliKFParticle*) pos->At(i);
		for( Int_t j=0; j < neg->GetEntriesFast(); ++j ){
			AliKFParticle* negKF = (AliKFParticle*) neg->At(j);
			Double_t psiPair = GetPsiPair(posKF, negKF);
			Double_t deltaPhi = fMagFieldSign * TVector2::Phi_mpi_pi( negKF->GetPhi() - posKF->GetPhi());
			fHistograms->FillHistogram(hName, deltaPhi, psiPair);
		}
	}
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::FillAngle(const TClonesArray* x, const TClonesArray* y, const TString& hName)
{
//
// Fill histogram with angle(x,y)
//
	for(Int_t i=0; i < x->GetEntriesFast(); ++i ){
		AliKFParticle* xKF = (AliKFParticle*) x->At(i);
		TVector3 xMom(xKF->Px(),xKF->Py(),xKF->Pz());
		for( Int_t j=0; j < y->GetEntriesFast(); ++j ){
			AliKFParticle* yKF = (AliKFParticle*) y->At(j);
			TVector3 yMom(yKF->Px(),yKF->Py(),yKF->Pz());
			fHistograms->FillHistogram(hName, xMom.Angle(yMom));
		}
	}
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::FillPidTable(const TParticle* p, Int_t pid)
{
//
// Fill table with pid info
//
	Int_t iGen=-1;
	switch(TMath::Abs(p->GetPdgCode())){
		case ::kElectron:   iGen=0; break;
		case ::kMuonMinus:  iGen=1; break;
		case ::kPiPlus:     iGen=2; break;
		case ::kKPlus:      iGen=3; break;
		case ::kProton:     iGen=4; break;
		default: iGen=-1;
	}
	
	int jRec=-1;
	if(pid > -1 && pid < 5) jRec = pid;
	
	if ((iGen > -1) && (jRec > -1)){
		fHistograms->FillTable("Table_PID", iGen, jRec);
		// sum
		fHistograms->FillTable("Table_PID", iGen, 5);
		fHistograms->FillTable("Table_PID", 5, jRec);
	}
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::GetGammaCandidates(TClonesArray*& gamma, vector<Int_t>& posIndex, vector<Int_t>& negIndex)
{
//
// Make a copy of gamma candidates from V0reader
//
	posIndex.clear();
	negIndex.clear();
	
	if(gamma) delete gamma;
	
	TClonesArray* gammaV0 = fV0Reader->GetCurrentEventGoodV0s();
	
	gamma = new TClonesArray("AliKFParticle", gammaV0->GetEntriesFast());
	gamma->SetOwner(kTRUE);
	
	// make a copy
	for(Int_t i=0; i < gammaV0->GetEntriesFast(); ++i){
		AliKFParticle* gamKF = (AliKFParticle*)gammaV0->At(i);
		new ((*gamma)[i]) AliKFParticle(*gamKF);
		posIndex.push_back(fV0Reader->GetPindex(i));
		negIndex.push_back(fV0Reader->GetNindex(i));
	}
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::IndexToAliKFParticle(const vector<Int_t>& index, Int_t PDG)
{
//
// Convert track index vector to AliKFParticle array
//
	TClonesArray* indexKF = new TClonesArray("AliKFParticle",index.size());
	indexKF->SetOwner(kTRUE);
	
	for(UInt_t i = 0; i < index.size(); ++i){
		AliESDtrack* t = fESDEvent->GetTrack(index[i]);
		new((*indexKF)[i]) AliKFParticle(*t->GetConstrainedParam(), PDG);
	}
	
	return indexKF;
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindElectronFromPi0Dalitz(const vector<Int_t>& candidates, Int_t PDG)
{
//
// Find true electrons from pi0 Dalitz decay candidates with MC
//
	TClonesArray* elec = new TClonesArray("AliKFParticle");
	elec->SetOwner(kTRUE);
	
	for(UInt_t i=0, j=0; i < candidates.size(); ++i){
		AliESDtrack* track = fESDEvent->GetTrack(candidates[i]);
		Int_t trackLabel = TMath::Abs(track->GetLabel());
		
		if( fStack->Particle(trackLabel)->GetPdgCode() != PDG ) continue;
		if( !IsPi0DalitzDaughter(trackLabel) ) continue;
		
		new ((*elec)[j++]) AliKFParticle(*track->GetConstrainedParam(), PDG);
	}
	
	return elec;
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindGammaFromPi0Dalitz(const TClonesArray* gamma, const vector<Int_t>& posIdx, const vector<Int_t>& negIdx)
{
//
// Find true gammas from pi0 Dalitz decay candidates with MC
//
	TClonesArray* gammaPi0 = new TClonesArray("AliKFParticle");
	gammaPi0->SetOwner(kTRUE);
	
	for(Int_t i=0, j=0; i < gamma->GetEntriesFast(); ++i) {
		AliKFParticle* gamKF = (AliKFParticle*)gamma->At(i);
		
		Int_t labelv1 = TMath::Abs((fESDEvent->GetTrack(posIdx[i]))->GetLabel());
		Int_t labelv2 = TMath::Abs((fESDEvent->GetTrack(negIdx[i]))->GetLabel());
		
		if( !HaveSameMother(labelv1,labelv2) ) continue;
		
		Int_t labelGamma = TMath::Abs(fStack->Particle(labelv1)->GetMother(0));
		
		if( fStack->Particle(labelGamma)->GetPdgCode() != ::kGamma ) continue;
		
		if( !IsPi0DalitzDaughter( labelGamma) ) continue;
		
		new ((*gammaPi0)[j++]) AliKFParticle(*gamKF);
	}
	
	return gammaPi0;
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindGamma(const TClonesArray* gamma, const vector<Int_t>& posIdx, const vector<Int_t>& negIdx)
{
//
// Find true gammas from gamma candidates with MC
//
	TClonesArray* gammaConv = new TClonesArray("AliKFParticle");
	gammaConv->SetOwner(kTRUE);
	
	for(Int_t i=0, j=0; i < gamma->GetEntriesFast(); ++i){
		AliKFParticle* gamKF = (AliKFParticle*)gamma->At(i);
		
		Int_t labelv1 = TMath::Abs((fESDEvent->GetTrack(posIdx[i]))->GetLabel());
		Int_t labelv2 = TMath::Abs((fESDEvent->GetTrack(negIdx[i]))->GetLabel());
		
		if( !HaveSameMother(labelv1,labelv2) ) continue;
		
		Int_t labelGamma = TMath::Abs(fStack->Particle(labelv1)->GetMother(0));
		
		if( fStack->Particle(labelGamma)->GetPdgCode() != ::kGamma ) continue;
		
		new ((*gammaConv)[j++]) AliKFParticle(*gamKF);
	}
	
	return gammaConv;
}

//--------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::ESDtrackIndexCut(vector<Int_t>& pos, vector<Int_t>& neg, const TClonesArray* gamma)
{
//
// Remove repeated electron candidate tracks
// according to the gamma candidate array
//
	vector<Bool_t> posTag(pos.size(),kTRUE);
	vector<Bool_t> negTag(neg.size(),kTRUE);
	
	for(Int_t i=0; i < gamma->GetEntriesFast(); ++i){
		Int_t gamPosIndex = fGammaCandidatePosIndex[i];
		Int_t gamNegIndex = fGammaCandidateNegIndex[i];
		
		for( UInt_t j=0; j < pos.size(); ++j ){
			if(pos[j] == gamPosIndex || pos[j] == gamNegIndex) posTag[j] = kFALSE;
		}
		
		for( UInt_t j=0; j < neg.size(); ++j )	{
			if(neg[j] == gamPosIndex || neg[j] == gamNegIndex) negTag[j] = kFALSE;
		}
	}
	
	CleanArray(pos, posTag);
	CleanArray(neg, negTag);
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::PsiPairCut(vector<Int_t>& pos, vector<Int_t>& neg)
{
//
// Remove electron candidates from gamma conversions
// according to the Psi pair angle
//
	vector<Bool_t> posTag(pos.size(), kTRUE);
	vector<Bool_t> negTag(neg.size(), kTRUE);

	for( UInt_t i=0; i < pos.size(); ++i ){
		AliESDtrack* posTrack = fESDEvent->GetTrack(pos[i]);

		for( UInt_t j=0; j < neg.size(); ++j ){
			AliESDtrack* negTrack = fESDEvent->GetTrack(neg[j]);

			Double_t psiPair = GetPsiPair(posTrack, negTrack);
			Double_t deltaPhi = fMagFieldSign * TVector2::Phi_mpi_pi( negTrack->GetConstrainedParam()->Phi()-posTrack->GetConstrainedParam()->Phi());

			if(IsFromGammaConversion( psiPair, deltaPhi )) {
					posTag[i] = kFALSE;
					negTag[j] = kFALSE;
			}
		}
	}

	CleanArray(pos, posTag);
	CleanArray(neg, negTag);
}

//-----------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::MassCut(vector<Int_t>& pos, vector<Int_t>& neg)
{
//
// Remove electron candidates pairs 
// which have mass not in the range (fMassCutMin,fMassCutMax)
//
	vector<Bool_t> posTag(pos.size(), kTRUE);
	vector<Bool_t> negTag(neg.size(), kTRUE);

	for( UInt_t i=0; i < pos.size(); ++i ) {
		AliESDtrack* posTrack = fESDEvent->GetTrack(pos[i]);

		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		TLorentzVector posLV;
		posLV.SetXYZM(posMom[0],posMom[1],posMom[2],fkElectronMass);

		for( UInt_t j=0; j < neg.size(); ++j ){
			AliESDtrack* negTrack = fESDEvent->GetTrack(neg[j]);

			Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
			TLorentzVector negLV;
			negLV.SetXYZM(negMom[0],negMom[1],negMom[2],fkElectronMass);

			TLorentzVector posnegLV = posLV + negLV;

			if( (posnegLV.M() < fMassCutMin) || (posnegLV.M() > fMassCutMax) ){
					posTag[i] = kFALSE;
					negTag[j] = kFALSE;
			}
		}
	}

	CleanArray(pos, posTag);
	CleanArray(neg, negTag);
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::CleanArray(vector<Int_t>& x, const vector<Bool_t>& tag)
{
//
// Clean the x array according to the tag parameter
//
	vector<Int_t> tmp;
	for(UInt_t i=0; i< x.size(); ++i){
		if(tag[i]) tmp.push_back(x[i]);
	}
	
	x = tmp;
}

//--------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::AngleEposEnegGammaCut( const vector<Int_t>& posIdx, const vector<Int_t>& negIdx, const TClonesArray* candidates, TClonesArray*& gamma, vector<Int_t>& posGamIdx, vector<Int_t>& negGamIdx)
{
//
// Remove gamma candidates according to
// the angle between the plane e+,e- and the gamma
//
	vector<Bool_t> gammaTag(candidates->GetEntriesFast(), kTRUE);
	
	for( UInt_t iPos=0; iPos < posIdx.size(); ++iPos ) {
		AliESDtrack* posTrack = fESDEvent->GetTrack(posIdx[iPos]);
		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		TVector3 xMom(posMom[0],posMom[1],posMom[2]);
		
		for( UInt_t iNeg=0; iNeg < negIdx.size(); ++iNeg ) {
			AliESDtrack* negTrack = fESDEvent->GetTrack(negIdx[iNeg]);
			Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
			TVector3 yMom(negMom[0],negMom[1],negMom[2]);
			
			// normal vector to x+y- plane
			TVector3 planePosNeg = xMom.Cross(yMom);
			for(Int_t i=0; i < candidates->GetEntriesFast(); ++i) {
				AliKFParticle* gamKF = (AliKFParticle*)candidates->At(i);
				TVector3 gamMom(gamKF->Px(),gamKF->Py(),gamKF->Pz());
				if (planePosNeg.Angle(gamMom) < 1. || planePosNeg.Angle(gamMom) > 2.) {
					gammaTag[i] = kFALSE;
				}
			}
		}
	}
	
	// Rebuild gamma candidates array
	
	if(gamma) delete gamma;
	gamma = new TClonesArray("AliKFParticle");
	gamma->SetOwner(kTRUE);
	
	posGamIdx.clear();
	negGamIdx.clear();
	
	for(Int_t i=0, j=0; i < candidates->GetEntriesFast(); ++i) {
		AliKFParticle* iGamma = (AliKFParticle*)candidates->At(i);
		if(gammaTag[i]) {
			new ((*gamma)[j++]) AliKFParticle(*iGamma);
			posGamIdx.push_back(fV0Reader->GetPindex(i));
			negGamIdx.push_back(fV0Reader->GetNindex(i));
		}
	}
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindDalitzPair(const TClonesArray* pos, const TClonesArray* neg)
{
//
// Find Dalitz pair candidates
//
	TClonesArray* dalitz = new TClonesArray("TLorentzVector");
	dalitz->SetOwner(kTRUE);
	
	for( Int_t iPos=0, j=0; iPos < pos->GetEntriesFast(); ++iPos ) {
		AliKFParticle* posKF = (AliKFParticle*)pos->At(iPos);
		
		TLorentzVector posLV;
		posLV.SetXYZM(posKF->Px(),posKF->Py(),posKF->Pz(),fkElectronMass);
		
		for( Int_t iNeg=0; iNeg < neg->GetEntriesFast(); ++iNeg ) {
			AliKFParticle* negKF = (AliKFParticle*)neg->At(iNeg);
			
			TLorentzVector negLV;
			negLV.SetXYZM(negKF->Px(),negKF->Py(),negKF->Pz(),fkElectronMass);
			
			if(fUseAliKF) {
				AliKFParticle posNeg( *posKF, *negKF);
				
				TLorentzVector posNegLV;
				posNegLV.SetXYZM(posNeg.Px(), posNeg.Py(), posNeg.Pz(), posNeg.GetMass());
				new ((*dalitz)[j++]) TLorentzVector(posNegLV);
			} else {
				new ((*dalitz)[j++]) TLorentzVector(posLV + negLV);
			}
		}
	}
	
	return dalitz;
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindParticleDalitz(const TClonesArray* pos, const TClonesArray* neg, const TClonesArray* gamma,Int_t opc)
{
//
// Find pi0 Dalitz decay candidates
//
	TClonesArray* pi0 = new TClonesArray("TLorentzVector");
	pi0->SetOwner(kTRUE);
	
	for( Int_t iPos=0, j=0; iPos < pos->GetEntriesFast(); ++iPos ){
		AliKFParticle* posKF = (AliKFParticle*)pos->At(iPos);
		
		TLorentzVector posLV;
		posLV.SetXYZM(posKF->Px(),posKF->Py(),posKF->Pz(),fkElectronMass);
		
		for( Int_t iNeg=0; iNeg < neg->GetEntriesFast(); ++iNeg ) {
			AliKFParticle* negKF = (AliKFParticle*)neg->At(iNeg);
			
			TLorentzVector negLV;
			negLV.SetXYZM(negKF->Px(),negKF->Py(),negKF->Pz(),fkElectronMass);
			AliKFParticle posNegKF(*posKF,*negKF);

			for(Int_t iGam=0; iGam < gamma->GetEntriesFast(); ++iGam){
				AliKFParticle* gamKF = (AliKFParticle*)gamma->At(iGam);
				AliKFParticle posNegGam( *posKF, *negKF, *gamKF );
				
				Double_t lDiffMass = posNegGam.GetMass() - posNegKF.GetMass(); 

				if( opc == 1 ) {
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff",lDiffMass );
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff_vs_Pt",lDiffMass,posNegGam.GetPt());
				} else if ( opc == 2 ) {
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff",lDiffMass );
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff_vs_Pt",lDiffMass,posNegGam.GetPt());
				} else if ( opc == 3 ) {
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff",lDiffMass );
					fHistograms->FillHistogram("ESD_BKG_BGHandler_InvMassDiff_vs_Pt",lDiffMass,posNegGam.GetPt());
				}

				if(fUseAliKF) {
					TLorentzVector posNegGamLV;
					posNegGamLV.SetXYZM(posNegGam.Px(),posNegGam.Py(),posNegGam.Pz(),posNegGam.GetMass());
					new ((*pi0)[j++]) TLorentzVector(posNegGamLV);
				} else {
					TLorentzVector gamLV;
					gamLV.SetXYZM(gamKF->Px(),gamKF->Py(),gamKF->Pz(),0);
					new ((*pi0)[j++]) TLorentzVector(posLV + negLV + gamLV);
				}
			}
		}
	}
	
	return pi0;
}

//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindDalitzPair(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx,Int_t motherOpc)
{
//
// Find true Dalitz pairs from Dalitz pair candidats with MC
//
	TClonesArray* dalitz = new TClonesArray("TLorentzVector");
	dalitz->SetOwner(kTRUE);
	
	for( UInt_t iPos=0, j=0; iPos < posIdx.size(); ++iPos ) {
		AliESDtrack* posTrack = fESDEvent->GetTrack(posIdx[iPos]);
		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		Int_t posLabel = TMath::Abs(posTrack->GetLabel());
		
		TLorentzVector posLV;
		posLV.SetXYZM(posMom[0],posMom[1],posMom[2],fkElectronMass);
		
		AliKFParticle posKF( *posTrack->GetConstrainedParam(), ::kPositron );
		
		for( UInt_t iNeg=0; iNeg < negIdx.size(); ++iNeg ) {
			AliESDtrack* negTrack = fESDEvent->GetTrack(negIdx[iNeg]);
			Int_t negLabel = TMath::Abs(negTrack->GetLabel());
			
			if(!IsDalitzPair(posLabel,negLabel,motherOpc)) continue;
			
			if(fUseAliKF) {
				AliKFParticle negKF( *negTrack->GetConstrainedParam(), ::kElectron );
				AliKFParticle posNeg( posKF, negKF);
				
				TLorentzVector posNegLV;
				posNegLV.SetXYZM(posNeg.Px(),posNeg.Py(),posNeg.Pz(),posNeg.GetMass());
				
				new ((*dalitz)[j++]) TLorentzVector(posNegLV);
			} else {// TLorentzVector
			
				Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
				
				TLorentzVector negLV;
				negLV.SetXYZM(negMom[0],negMom[1],negMom[2],fkElectronMass);
				
				new ((*dalitz)[j++]) TLorentzVector(posLV + negLV);
			}
		}
	}
	return dalitz;
}

TClonesArray* AliAnalysisTaskGammaConvDalitz::FindJpsi(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx,Int_t motherOpc)
{
//
// Find true Jpsi's
// If mother
// -1: from the all sources
// 0: from the Chic_0
// 1: from the Chic_1
// 2: from the Chic_2

	TClonesArray* jPsi = new TClonesArray("TLorentzVector");
	jPsi->SetOwner(kTRUE);
	
	for( UInt_t iPos=0, j=0; iPos < posIdx.size(); ++iPos )
	{
		AliESDtrack* posTrack = fESDEvent->GetTrack(posIdx[iPos]);
		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		Int_t posLabel = TMath::Abs(posTrack->GetLabel());

		if( fStack->Particle(posLabel)->GetPdgCode() != ::kPositron ) continue;
		
		TLorentzVector posLV;
		posLV.SetXYZM(posMom[0],posMom[1],posMom[2],fkElectronMass);
		
		AliKFParticle posKF( *posTrack->GetConstrainedParam(), ::kPositron );
		
		for( UInt_t iNeg=0; iNeg < negIdx.size(); ++iNeg ) {
			AliESDtrack* negTrack = fESDEvent->GetTrack(negIdx[iNeg]);
			Int_t negLabel = TMath::Abs(negTrack->GetLabel());

			if( fStack->Particle(negLabel)->GetPdgCode() != ::kElectron ) continue;
			
			if( !HaveSameMother(posLabel,negLabel) ) continue;
			
			Int_t motherLabel = fStack->Particle(negLabel)->GetMother(0);
			TParticle *motherJpsi = fStack->Particle(motherLabel);

			if( motherJpsi->GetPdgCode() != 443 ){
					continue;
			}
			if( motherOpc > -1 ) { 
				if( motherJpsi->GetMother(0) < 0 ) continue;
				
				TParticle *gmotherChic = fStack->Particle(motherJpsi->GetMother(0));
				Int_t pdgCode = gmotherChic->GetPdgCode();

				Bool_t lson = kTRUE;
				
				switch(motherOpc){

						case 0:     if ( pdgCode != 10441 )
										lson = kFALSE;
										break;
						case 1:     if ( pdgCode != 20443 )
										lson = kFALSE;
										break;
						case 2:     if ( pdgCode != 445   )
										lson = kFALSE;
										break;
				}
				if( lson == kFALSE ) continue;
			}
			if(fUseAliKF) {
				AliKFParticle negKF( *negTrack->GetConstrainedParam(), ::kElectron );
				AliKFParticle posNeg( posKF, negKF);
				
				TLorentzVector posNegLV;
				posNegLV.SetXYZM(posNeg.Px(),posNeg.Py(),posNeg.Pz(),posNeg.GetMass());
				
				new ((*jPsi)[j++]) TLorentzVector(posNegLV);
			} else {// TLorentzVector
			
						Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
						
						TLorentzVector negLV;
						negLV.SetXYZM(negMom[0],negMom[1],negMom[2],fkElectronMass);
						
						new ((*jPsi)[j++]) TLorentzVector(posLV + negLV);
			}
		}
	}
	return jPsi;
}



//--------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindParticleDalitz(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx, const TClonesArray* gamma, const vector<Int_t>& posGam, const vector<Int_t>& negGam,Int_t motherOpc)
{
//
// Find true pi0 Dalitz decay from pi0 candidates with MC
//
	TClonesArray* pi0 = new TClonesArray("TLorentzVector");
	pi0->SetOwner(kTRUE);
	
	for( UInt_t iPos=0, j=0; iPos < posIdx.size(); ++iPos ) {
		AliESDtrack* posTrack = fESDEvent->GetTrack(posIdx[iPos]);
		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		Int_t posLabel = TMath::Abs(posTrack->GetLabel());
		
		TLorentzVector posLV;
		posLV.SetXYZM(posMom[0],posMom[1],posMom[2],fkElectronMass);
		
		AliKFParticle posKF( *posTrack->GetConstrainedParam(), ::kPositron );
		
		for( UInt_t iNeg=0; iNeg < negIdx.size(); ++iNeg ) {
			AliESDtrack* negTrack = fESDEvent->GetTrack(negIdx[iNeg]);
			Int_t negLabel = TMath::Abs(negTrack->GetLabel());
		
			if( !HaveSameMother(posLabel,negLabel) ) continue; //Check if both particles have same mother
			if(!IsDalitzPair(posLabel,negLabel,motherOpc)) continue; //check if mohter is eta0 or pi0
		
			Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
			
			TLorentzVector negLV;
			negLV.SetXYZM(negMom[0],negMom[1],negMom[2],fkElectronMass);
			
			AliKFParticle negKF( *negTrack->GetConstrainedParam(), ::kElectron );
			
			for(Int_t iGam=0; iGam < gamma->GetEntriesFast(); ++iGam) {
				AliKFParticle* gamKF = (AliKFParticle*)gamma->At(iGam);
				
				Int_t labelv1 = TMath::Abs((fESDEvent->GetTrack(posGam[iGam]))->GetLabel());
				Int_t labelv2 = TMath::Abs((fESDEvent->GetTrack(negGam[iGam]))->GetLabel());
				
				if( !HaveSameMother(labelv1,labelv2) ) continue;
				
				Int_t labelGamma = TMath::Abs(fStack->Particle(labelv1)->GetMother(0));
				
				if( fStack->Particle(labelGamma)->GetPdgCode() != ::kGamma ) continue;

				if( !HaveSameMother(labelGamma, posLabel) ) continue;
				
				if(fUseAliKF) {
					AliKFParticle posNegGam( posKF, negKF, *gamKF );
					TLorentzVector posNegGamLV;
					posNegGamLV.SetXYZM(posNegGam.Px(),posNegGam.Py(),posNegGam.Pz(),posNegGam.GetMass());
					new ((*pi0)[j++]) TLorentzVector(posNegGamLV);
				} else { // TLorentzVector
					TLorentzVector gamLV;
					gamLV.SetXYZM(gamKF->Px(),gamKF->Py(),gamKF->Pz(),0);
					
					new ((*pi0)[j++]) TLorentzVector(posLV + negLV + gamLV);
				}
			}
		}
	}
	
	return pi0;
}
TClonesArray* AliAnalysisTaskGammaConvDalitz::FindParticleChic(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx, const TClonesArray* gamma, const vector<Int_t>& posGam, const vector<Int_t>& negGam,Int_t motherOpc)
{
//
// Find true pi0 Dalitz decay from pi0 candidates with MC
//
	TClonesArray* chic = new TClonesArray("TLorentzVector");
	chic->SetOwner(kTRUE);
	
	for( UInt_t iPos=0, j=0; iPos < posIdx.size(); ++iPos ) {
		AliESDtrack* posTrack = fESDEvent->GetTrack(posIdx[iPos]);
		Double_t posMom[3]; posTrack->GetConstrainedPxPyPz(posMom);
		Int_t posLabel = TMath::Abs(posTrack->GetLabel());

		if( fStack->Particle(posLabel)->GetPdgCode() != ::kPositron ) continue;
		
		TLorentzVector posLV;
		posLV.SetXYZM(posMom[0],posMom[1],posMom[2],fkElectronMass);
		
		AliKFParticle posKF( *posTrack->GetConstrainedParam(), ::kPositron );
		
		for( UInt_t iNeg=0; iNeg < negIdx.size(); ++iNeg ) {
			AliESDtrack* negTrack = fESDEvent->GetTrack(negIdx[iNeg]);
			Int_t negLabel = TMath::Abs(negTrack->GetLabel());

			if( fStack->Particle(negLabel)->GetPdgCode() != ::kElectron ) continue;
			if( !HaveSameMother(posLabel,negLabel) ) continue;

			Int_t jpsiLabel = fStack->Particle(negLabel)->GetMother(0);

			if( fStack->Particle(jpsiLabel)->GetPdgCode() != 443 ) continue;
			
			TParticle *jpsiParticle = fStack->Particle(jpsiLabel);

			if ( jpsiParticle->GetMother(0) < 0 ) continue; 
			Int_t chicLabel = jpsiParticle->GetMother(0);

			Int_t pdgCode = fStack->Particle(chicLabel)->GetPdgCode();
			Bool_t lSon = kTRUE;
					
			switch(motherOpc){
				case 0:     if ( pdgCode != 10441 )
								lSon = kFALSE;
								break;
				case 1:     if ( pdgCode != 20443 )
								lSon = kFALSE;
								break;
				case 2:     if ( pdgCode != 445   )
								lSon = kFALSE;
								break;
			}


			if( lSon == kFALSE ) continue;

			Double_t negMom[3]; negTrack->GetConstrainedPxPyPz(negMom);
			
			TLorentzVector negLV;
			negLV.SetXYZM(negMom[0],negMom[1],negMom[2],fkElectronMass);
			
			AliKFParticle negKF( *negTrack->GetConstrainedParam(), ::kElectron );
			
			for(Int_t iGam=0; iGam < gamma->GetEntriesFast(); ++iGam) {
				AliKFParticle* gamKF = (AliKFParticle*)gamma->At(iGam);
				
				Int_t labelv1 = TMath::Abs((fESDEvent->GetTrack(posGam[iGam]))->GetLabel());
				Int_t labelv2 = TMath::Abs((fESDEvent->GetTrack(negGam[iGam]))->GetLabel());
				
				if( !HaveSameMother(labelv1,labelv2) ) continue;
				
				Int_t labelGamma = TMath::Abs(fStack->Particle(labelv1)->GetMother(0));
				
				if( fStack->Particle(labelGamma)->GetPdgCode() != ::kGamma ) continue;

				if( !HaveSameMother(labelGamma, jpsiLabel) ) continue;

				if(fUseAliKF) {
					AliKFParticle posNegGam( posKF, negKF, *gamKF );
					TLorentzVector posNegGamLV;
					posNegGamLV.SetXYZM(posNegGam.Px(),posNegGam.Py(),posNegGam.Pz(),posNegGam.GetMass());
					new ((*chic)[j++]) TLorentzVector(posNegGamLV);
				} else {// TLorentzVector
		
					TLorentzVector gamLV;
					gamLV.SetXYZM(gamKF->Px(),gamKF->Py(),gamKF->Pz(),0);
					
					new ((*chic)[j++]) TLorentzVector(posLV + negLV + gamLV);
				}
			}
		}
	}
	
	return chic;
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskGammaConvDalitz::UpdateGammaPool(const TClonesArray* gamma)
{
//
// Update gamma event pool for background computation
//
	if( fDebug ) AliInfo("=> UpdateGammaPool");
	
	// cycle
	for(Int_t j=0; j< gamma->GetEntriesFast(); ++j) {
		if((AliKFParticle*)fGammaPool->At(fGamPoolPos)) delete (AliKFParticle*)fGammaPool->RemoveAt(fGamPoolPos);
		new ((*fGammaPool)[fGamPoolPos]) AliKFParticle( *((AliKFParticle*)gamma->At(j)));
		++fGamPoolPos;
		if(fGamPoolPos == fPoolMaxSize) {
			fGamPoolPos = 0;
		}
	}
}

void AliAnalysisTaskGammaConvDalitz::UpdateElectronPool(TClonesArray* elec) // FIXME: const
{
//
// Update electron event pool for background computation
//
	Int_t multiplicity = fV0Reader->CountESDTracks();
	fBGEventHandler->AddElectronEvent(elec,fESDEvent->GetPrimaryVertex()->GetZ(),multiplicity);
}

//-----------------------------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::GammasFromBGHandler() const
{
//
// Gamma copy from events with same multiplicity and Z
//
	if( fDebug ) AliInfo("=> GammasFromBGHandler");

	Int_t zbin = fBGEventHandler->GetZBinIndex(fV0Reader->GetVertexZ());
	Int_t mbin = fBGEventHandler->GetMultiplicityBinIndex(fV0Reader->CountESDTracks());

	TClonesArray* gammaPool = new TClonesArray("AliKFParticle");
	gammaPool->SetOwner(kTRUE);
	
	for( Int_t iEventBG=0; iEventBG < fV0Reader->GetNBGEvents(); ++iEventBG ) {
		AliGammaConversionKFVector* gammaV0s = fBGEventHandler->GetBGGoodV0s(zbin,mbin,iEventBG);
		for( UInt_t i = 0; i < gammaV0s->size(); ++i) {
			new ((*gammaPool)[i]) AliKFParticle( *((AliKFParticle*)gammaV0s->at(i)) );
		}
	}
	return gammaPool;
}

//-----------------------------------------------------------------------------------------------
TClonesArray* AliAnalysisTaskGammaConvDalitz::ElectronFromBGHandler() const
{
//
// Electron copy from events with same multiplicity and Z
//
	if( fDebug ) AliInfo("=> ElectronFromBGHandler");
	
	TClonesArray* electronPool = new TClonesArray("AliKFParticle");
	electronPool->SetOwner(kTRUE);

	Int_t multiplicity = fV0Reader->CountESDTracks();
	for( Int_t iEventBG=0; iEventBG < fV0Reader->GetNBGEvents(); ++iEventBG ) {
		AliGammaConversionKFVector* electronNeg =  fBGEventHandler->GetBGGoodENeg(iEventBG,fESDEvent->GetPrimaryVertex()->GetZ(),multiplicity);
		for (UInt_t i = 0; i < electronNeg->size(); ++i  ) {
			new ((*electronPool)[i]) AliKFParticle( *((AliKFParticle*)electronNeg->at(i)) );
		}
	}
	
	return electronPool;
}

//-----------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskGammaConvDalitz::GetMonteCarloPid(const AliESDtrack* t) const
{
//
// Get track pid according to MC
//
    Int_t label   = TMath::Abs(t->GetLabel());
    Int_t pdgCode = TMath::Abs(fStack->Particle(label)->GetPdgCode());

    switch(pdgCode) {
		case ::kElectron:  return AliPID::kElectron;
		case ::kMuonMinus: return AliPID::kMuon;
		case ::kPiPlus:    return AliPID::kPion;
		case ::kKPlus:     return AliPID::kKaon;
		case ::kProton:    return AliPID::kProton;
    }
    return -1;
}

//-----------------------------------------------------------------------------------------------
//FIXME PID ITS
// NOTE prior should be estimated from data
// NOTE: move to config

Int_t AliAnalysisTaskGammaConvDalitz::GetBayesPid(const AliESDtrack* t, Int_t trackType ) const
{
//
// Get track pid according to Bayes' formula
//
	double priors[AliPID::kSPECIES] = {0.009, 0.01, 0.82, 0.10, 0.05};
	Double_t detectoProb[AliPID::kSPECIES];

	if( trackType == kITSsaTrack ) { // ITS standalone pid
		t->GetITSpid( detectoProb );
	} else { // global track
		t->GetESDpid( detectoProb );
	}

	AliPID bayesPID( detectoProb );
	return bayesPID.GetMostProbable( priors );
}




Int_t AliAnalysisTaskGammaConvDalitz::GetNSigmaPid(const AliESDtrack* track, Int_t trackType ) const
{
//
// Get track pid according to a n-sigma cut around ITS and/or TPC signals
//

	if( trackType == kITSsaTrack) {  // ITS standalone tracks
		Double_t mom = track->GetP();

		// ITS Number of sigmas (FIXME: add new fESDpidCuts)
		// NOTE there is not AliESDpidCuts::SetITSnSigmaCut yet
		Double_t nElecSigma = fESDpid->NumberOfSigmasITS(track, AliPID::kElectron );
		Double_t nPionSigma = fESDpid->NumberOfSigmasITS(track, AliPID::kPion );

		if( nElecSigma < 4. && nElecSigma > -3. && mom < .2  &&  nPionSigma < -1.5 ) {
			return AliPID::kElectron;
		}
	} else { // global track
		fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_E_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron));

		if( fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
				fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron)>fPIDnSigmaAboveElectronLine){
					return -1;
		}
		fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_E_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron));
		fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_Pi_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
		fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_Mu_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kMuon));

		if( track->P()>fPIDMinPnSigmaAbovePionLine && track->P()<fPIDMaxPnSigmaAbovePionLine ){
			if( fESDpid->NumberOfSigmasTPC(track,AliPID::kPion)<fPIDnSigmaAbovePionLine){
				return -1;
			}
		}
		
		// High Pt
		
		if( track->P()>fPIDMaxPnSigmaAbovePionLine ){
			if( fESDpid->NumberOfSigmasTPC(track,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){
				return -1;
			}
		}

		fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_Pi_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
		fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_Mu_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kMuon));

		if(fDoKaonRejectionLowP == kTRUE){
			
			fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_K_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));
			if( track->P()<fPIDMinPKaonRejectionLowP ){      
				if( TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
					return -1;
				}
			}

			fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_K_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));
		}

		if(fDoProtonRejectionLowP == kTRUE){
			fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_P_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kProton));
			if( track->P()<fPIDMinPProtonRejectionLowP ){ 	
				if( TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
					return -1;
				}
			}
			fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_P_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kProton));
		}

		if(fDoPionRejectionLowP == kTRUE){
			fHistograms->FillHistogram("ESD_DalitzBeforeCorresCut_Low_Pi_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
	
			if( track->P()<fPIDMinPPionRejectionLowP ){
				if( TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){
					return -1;
				}
			}
			fHistograms->FillHistogram("ESD_DalitzAfterCorresCut_Low_Pi_SigdEdxP",track->P(),fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
		}

		if( fDoTOFsigmaCut == kTRUE ){ 
			if( (track->GetStatus() & AliESDtrack::kTOFpid) && !(track->GetStatus() & AliESDtrack::kTOFmismatch) ){

				Double_t t0pos = fESDpid->GetTOFResponse().GetStartTime(track->P());
				Double_t nSigma = fESDpid->NumberOfSigmasTOF(track, AliPID::kElectron, t0pos);
				fHistograms->FillHistogram("ESD_TOF_BeforeCorresCutE_P_SigdT", track->GetP(), nSigma); //if( (statusNeg & AliESD 
				if( (nSigma>fTofPIDnSigmaAboveElectronLine) || (nSigma<fTofPIDnSigmaBelowElectronLine) ){
					return -1;
				}    
				fHistograms->FillHistogram("ESD_TOF_AfterCorresCutE_P_SigdT", track->GetP(), nSigma); //if( (statusNeg & AliESD 
			}
		}
		return AliPID::kElectron;

	// NOTE: add other particle types
	}
	return -1;
}


//-----------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskGammaConvDalitz::IsDalitzPair( Int_t posLabel, Int_t negLabel,Int_t motherOpc ) const
{
//
// Returns true if the two particles is a Dalitz pair
//
//motherOpc 1:  for Pi0Dalitz
//motherOpc 2:  for EtaDalitz 

	if(!HaveSameMother(posLabel, negLabel)) return kFALSE;

	TParticle* pos = fStack->Particle( posLabel );
	TParticle* neg = fStack->Particle( negLabel );

	if( pos->GetPdgCode() != ::kPositron ) return kFALSE;
	if( neg->GetPdgCode() != ::kElectron ) return kFALSE;
	
	//if( pos->GetUniqueID() != ::kPDecay ) return kFALSE;
	//if( neg->GetUniqueID() != ::kPDecay ) return kFALSE;
	
	Int_t motherLabel = pos->GetMother(0);
	if( motherLabel < 0 ) return kFALSE;

	TParticle* mother = fStack->Particle( motherLabel );
		
	if( mother->GetNDaughters() != 3)    return kFALSE;

	if( motherOpc == 1 ){ //Pi0Dalitz 
		if( mother->GetPdgCode() != ::kPi0 ) return kFALSE;
	} else if(motherOpc == 2) {
		if( mother->GetPdgCode() != 221 ) return kFALSE; 
	} else {
		return kFALSE;
	}
	// NOTE: one of them must be a photon

	return kTRUE;
}

//-----------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskGammaConvDalitz::IsPi0DalitzDaughter( Int_t label ) const
{
//
// Returns true if the particle comes from Pi0 -> e+ e- gamma
//
	Bool_t ePlusFlag  = kFALSE;
	Bool_t eMinusFlag = kFALSE;
	Bool_t gammaFlag  = kFALSE;
	
	Int_t motherLabel = fStack->Particle( label )->GetMother(0);
	
	if( motherLabel < 0 ) return kFALSE;
	
	TParticle* mother = fStack->Particle( motherLabel );
	
	if ( mother->GetPdgCode() != ::kPi0 ) return kFALSE;
	
	if ( mother->GetNDaughters() != 3 ) return kFALSE;
	
	for( Int_t idx = mother->GetFirstDaughter(); idx <= mother->GetLastDaughter(); ++idx ) {
		switch( fStack->Particle(idx)->GetPdgCode()) {
			case ::kPositron:
				ePlusFlag  = kTRUE;
				break;
			case ::kElectron:
				eMinusFlag = kTRUE;
				break;
			case ::kGamma:
				gammaFlag  = kTRUE;
				break;
		}
	}
	return ( ePlusFlag && eMinusFlag && gammaFlag );
}

//--------------------------------------------------------------------------
Bool_t AliAnalysisTaskGammaConvDalitz::IsFromGammaConversion( Double_t psiPair, Double_t deltaPhi ) const
{
//
// Returns true if it is a gamma conversion according to psi pair value
//
	return ( (deltaPhi > fDeltaPhiCutMin  &&  deltaPhi < fDeltaPhiCutMax) &&
	TMath::Abs(psiPair) < ( fPsiPairCut - fPsiPairCut/fDeltaPhiCutMax * deltaPhi ) );
}

//--------------------------------------------------------------------------
Bool_t AliAnalysisTaskGammaConvDalitz::IsFromGammaConversion( Int_t posLabel, Int_t negLabel ) const
{
//
// Returns true if it is a gamma conversion according to MC
//
	if( !HaveSameMother(posLabel,negLabel) ) return kFALSE;
	
	TParticle* pos = fStack->Particle( posLabel );
	TParticle* neg = fStack->Particle( negLabel );
	
	if( pos->GetPdgCode() != ::kPositron ) return kFALSE;
	if( neg->GetPdgCode() != ::kElectron ) return kFALSE;

	if( pos->GetUniqueID() != kPPair ) return kFALSE;

	Int_t motherLabel = pos->GetMother(0);
	if( motherLabel < 0 ) return kFALSE;

	TParticle* mother = fStack->Particle( motherLabel );

	if( mother->GetPdgCode() != ::kGamma ) return kFALSE;
	
	return kTRUE;
}

//-----------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskGammaConvDalitz::HaveSameMother( Int_t label1, Int_t label2 ) const
{
//
// Returns true if the two particle have the same mother
//
	if(fStack->Particle( label1 )->GetMother(0) < 0 ) return kFALSE;
	return (fStack->Particle( label1 )->GetMother(0) == fStack->Particle( label2 )->GetMother(0));
}

//-----------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskGammaConvDalitz::GetPsiPair( const AliESDtrack* trackPos, const AliESDtrack* trackNeg ) const
{
//
// This angle is a measure for the contribution of the opening in polar
// direction Δ0 to the opening angle ξ Pair
//
// Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
//      Master Thesis. Thorsten Dahms. 2005
// https://twiki.cern.ch/twiki/pub/ALICE/GammaPhysicsPublications/tdahms_thesis.pdf
//
	Double_t momPos[3];
	Double_t momNeg[3];
	if( trackPos->GetConstrainedPxPyPz(momPos) == 0 ) trackPos->GetPxPyPz( momPos );
	if( trackNeg->GetConstrainedPxPyPz(momNeg) == 0 ) trackNeg->GetPxPyPz( momNeg );

	TVector3 posDaughter;
	TVector3 negDaughter;

	posDaughter.SetXYZ( momPos[0], momPos[1], momPos[2] );
	negDaughter.SetXYZ( momNeg[0], momNeg[1], momNeg[2] );

	Double_t deltaTheta = negDaughter.Theta() - posDaughter.Theta();
	Double_t openingAngle =  posDaughter.Angle( negDaughter );  //TMath::ACos( posDaughter.Dot(negDaughter)/(negDaughter.Mag()*posDaughter.Mag()) );
	if( openingAngle < 1e-20 ) return 0.;
	Double_t psiAngle = TMath::ASin( deltaTheta/openingAngle );

	return psiAngle;
}

//-----------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskGammaConvDalitz::GetPsiPair(const AliKFParticle* xPos, const AliKFParticle* yNeg ) const
{
//
// Get psi pair value
//
	TVector3 pos(xPos->GetPx(), xPos->GetPy(), xPos->GetPz());
	TVector3 neg(yNeg->GetPx(), yNeg->GetPy(), yNeg->GetPz());

	Double_t deltaTheta = neg.Theta() - pos.Theta();
	Double_t openingAngle = pos.Angle( neg );

	if( openingAngle < 1e-20 ) return 0.;

	return TMath::ASin( deltaTheta/openingAngle );
}
//-----------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskGammaConvDalitz::GetPsiPair(const TLorentzVector* xPos, const TLorentzVector* yNeg ) const
{
//
// Get psi pair value
//
	Double_t deltaTheta = yNeg->Theta() - xPos->Theta();
	Double_t openingAngle = xPos->Angle( yNeg->Vect() );
	
	if( openingAngle < 1e-20 ) return 0.;
	
	return TMath::ASin( deltaTheta/openingAngle );;
}
//------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskGammaConvDalitz::CalculateMultiplicityBin(){
	// Get Centrality bin

	Int_t multiplicity = 0;

	if ( fUseMultiplicity == 1 ) {
		if (fMultiplicity>= 0 && fMultiplicity<= 5) multiplicity=1;
		if (fMultiplicity>= 6 && fMultiplicity<= 9) multiplicity=2;
		if (fMultiplicity>=10 && fMultiplicity<=14) multiplicity=3;
		if (fMultiplicity>=15 && fMultiplicity<=22) multiplicity=4;
		if (fMultiplicity>=23 )multiplicity=5;
	}

	if ( fUseHBTMultiplicity == 1 ) {
		if (fMultiplicity>= 0 && fMultiplicity<=11)     multiplicity=1;
		if (fMultiplicity>=12 && fMultiplicity<=16)     multiplicity=2;
		if (fMultiplicity>=17 && fMultiplicity<=22)     multiplicity=3;
		if (fMultiplicity>=23 && fMultiplicity<=29)     multiplicity=4;
		if (fMultiplicity>=30 && fMultiplicity<=36)     multiplicity=5;
		if (fMultiplicity>=37 && fMultiplicity<=44)     multiplicity=6;
		if (fMultiplicity>=45 && fMultiplicity<=57)     multiplicity=7;
		if (fMultiplicity>=58 && fMultiplicity<=149) multiplicity=8;
	}

	return multiplicity;
}

