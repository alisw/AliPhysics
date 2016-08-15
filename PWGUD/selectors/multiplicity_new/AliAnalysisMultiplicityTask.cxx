#include "AliAnalysisMultiplicityTask.h"

#include <AliAnalysisManager.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliStack.h>
#include <TH1.h>
#include <TList.h>
#include "AliMultiplicityEventSelection.h"
#include "AliMultiplicityAnalysisSelection.h"

using std::map;

ClassImp ( AliAnalysisMultiplicityTask )
// Bool_t AliAnalysisMultiplicityTask::GetUseMultInV0() {
// 	return useMultInV0;
// }
//
// void AliAnalysisMultiplicityTask::SetUseMultInV0 ( Bool_t use ) {
// 	useMultInV0 = use;
// }


//_____________________________________________________________________________________________________________________________________________
AliAnalysisMultiplicityTask::AliAnalysisMultiplicityTask ( const char* name, const Int_t nOutput ) : AliAnalysisTaskSE ( name ),
	fExclusiveSelectionList ( 0 ),
	fAnalysisSelectionList ( 0 ),
	fOutput ( 0 ),
	fTrackCuts ( 0 ),
	fCollectedEvents ( 0 ),
	fUseMC ( kFALSE ),
	fNOutput ( nOutput ),
	fOutputMap(),
	fIsCocktail ( kFALSE ),
	nEvents ( 0 ),
// 	useMultInV0 ( kFALSE ),
	fEventClassification ( 0 ),
	fConsistencyCut ( 0.5 ),
	fSPDPileupThreshold ( 0.8 ),
	hVertex ( 0 ) {
	if ( nOutput < 1 ) {
		fNOutput = 0;
		AliFatalF ( "Insufficient number of outputs: %d < 1.", nOutput );
		return;
	}

	if ( nOutput > 100 ) {
		AliWarningF ( "Do you really need %d outputs?", nOutput );
	}

	for ( Int_t i = 1; i <= nOutput; i++ ) {
		DefineOutput ( i, TList::Class() );
	}
}
//_____________________________________________________________________________________________________________________________________________
AliAnalysisMultiplicityTask::~AliAnalysisMultiplicityTask() {
	//Clean-up
	if ( fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
		delete fOutput;
	}
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::SetUseMC ( Bool_t use ) {
	fUseMC = use;
}
//_____________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisMultiplicityTask::UseMC() {
	return fUseMC;
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::AddExclusiveSelection ( AliMultiplicityEventSelection* selection, const Int_t outputSlot ) {
	AddSelection ( kExclusiveSelection, selection, outputSlot );
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::AddAnalysisSelection ( AliMultiplicityEventSelection* selection,  const Int_t outputSlot ) {
	AddSelection ( kAnalysisSelection, selection, outputSlot );
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::AddSelection ( ESelection selType, AliMultiplicityEventSelection* selection, const Int_t outputSlot ) {
	if ( !selection ) {
		AliError ( "The selection pointer is NULL" );
		return;
	}

	TList* selList = 0;

	switch ( selType ) {
		case kExclusiveSelection:
			if ( !fExclusiveSelectionList ) {
				AliDebug ( 1, "Creating new Exclusive selections list." );
				fExclusiveSelectionList = new TList();
				fExclusiveSelectionList->SetOwner();
			}

			selList = fExclusiveSelectionList;
			break;

		case kAnalysisSelection:
			if ( !fAnalysisSelectionList ) {
				AliDebug ( 1, "Creating new Analysis selections list." );
				fAnalysisSelectionList = new TList();
				fAnalysisSelectionList->SetOwner();
			}

			selList = fAnalysisSelectionList;
			break;

		default
				:
			AliError ( "Unknown selection type" );
			return;
	}

	Int_t slot = outputSlot;

	if ( slot < 0 ) {
		AliErrorF ( "Invalid output slot: %d, setting to default 1", outputSlot );
		slot = 1;
	}

	if ( slot > fNOutput ) {
		AliErrorF ( "Output slot out of range: %d, setting to default 1", outputSlot );
		slot = 1;
	}

	selList->Add ( selection );
	fOutputMap[selection->GetName()] = slot;
	AliDebugF ( 1, "Selection %s to slot %d.", selection->GetName(), slot );
}
//_____________________________________________________________________________________________________________________________________________
AliESDEvent* AliAnalysisMultiplicityTask::ESDEvent() {
	AliESDEvent* event = dynamic_cast<AliESDEvent*> ( InputEvent() );

	if ( !event ) {
		AliError ( "Could not retrieve ESD event." );
		return 0;
	}

	return event;
}

//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::AddTrackCuts ( AliESDtrackCuts* externalTrackCuts ) {
	if ( !externalTrackCuts ) {
		AliWarning ( "No track cuts supplied, creating default." );
		CreateDefaultTrackCuts();
	} else {
		if ( !fTrackCuts ) {
			AliWarning ( "Creating new track cuts list." );
			fTrackCuts = new TList();
			fTrackCuts->SetOwner();
		}

		fTrackCuts->Add ( externalTrackCuts );
	}
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::AddTrackCuts ( TList* cutsList ) {
	if ( !cutsList || cutsList->IsEmpty() ) {
		AliWarning ( "Null or empty list!" );
		return;
	}

	fTrackCuts = cutsList;
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::CreateDefaultTrackCuts ( Int_t year ) {
	if ( fTrackCuts ) {
		AliWarning ( "Track cuts already exist!" );
// 		delete fTrackCuts;
		return;
	}

	fTrackCuts = new TList();
	fTrackCuts->SetOwner();

	AliESDtrackCuts* defaultcuts = 0;

	switch ( year ) {
		case 2010:
			defaultcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010 ( kTRUE, 1 ); // use number of crossed rows cut
			break;

		case 2011:
			defaultcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011 ( kTRUE, 1 ); // use number of crossed rows cut
			break;
	}

	TString ptdep ( defaultcuts->GetMaxDCAToVertexXYPtDep() );
	TString ptdepwo = TString::Format ( "1.5*(%s)", ptdep.Data() );

	AliESDtrackCuts* defaultcutsITSpure = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010 ( kTRUE, kFALSE ); //not for PID

	TString ptdepITS ( defaultcutsITSpure->GetMaxDCAToVertexXYPtDep() );
	TString ptdepITSwo = TString::Format ( "1.5*(%s)", ptdepITS.Data() );

	// quality cut on ITS+TPC tracks
	AliESDtrackCuts* GlobalTrackQuality = new AliESDtrackCuts ( "GlobalTrackQuality" );
	GlobalTrackQuality->SetMaxChi2PerClusterTPC ( 4 );

	GlobalTrackQuality->SetMinNCrossedRowsTPC ( 70 );
	GlobalTrackQuality->SetMinRatioCrossedRowsOverFindableClustersTPC ( 0.8 );

	GlobalTrackQuality->SetAcceptKinkDaughters ( kFALSE );
	GlobalTrackQuality->SetRequireTPCRefit ( kTRUE );
	GlobalTrackQuality->SetRequireITSRefit ( kTRUE );

	GlobalTrackQuality->SetDCAToVertex2D ( kFALSE );
	GlobalTrackQuality->SetRequireSigmaToVertex ( kFALSE );

	GlobalTrackQuality->SetMaxChi2PerClusterITS ( 36 );

	fTrackCuts->Add ( GlobalTrackQuality );
	//TPC only track quality
	AliESDtrackCuts* TPConlyQuality = new AliESDtrackCuts ( "TPConlyQuality" );
	TPConlyQuality->SetEtaRange ( -TPCetaLimit, TPCetaLimit );

	fTrackCuts->Add ( TPConlyQuality );

	// quality cut on pure ITS_SA tracks (complementary to ITS+TPC)
	AliESDtrackCuts* ITSpureQuality = new AliESDtrackCuts ( "ITSpureQuality" );

	ITSpureQuality->SetRequireITSRefit ( kTRUE );
	ITSpureQuality->SetMinNClustersITS ( 4 );
	ITSpureQuality->SetMaxChi2PerClusterITS ( 2.5 );
	ITSpureQuality->SetEtaRange ( -ITSetaLimit, ITSetaLimit );

	fTrackCuts->Add ( ITSpureQuality );

	// primary selection for global tracks with SPD hits
	AliESDtrackCuts* TrackWSPDhit = new AliESDtrackCuts ( "TrackWSPDhit" );
	TrackWSPDhit->SetClusterRequirementITS ( AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny );
	TrackWSPDhit->SetMaxDCAToVertexXYPtDep ( ptdep.Data() );
	TrackWSPDhit->SetMaxDCAToVertexZ ( 0.5 );
	TrackWSPDhit->SetMaxChi2TPCConstrainedGlobal ( 36 );

	fTrackCuts->Add ( TrackWSPDhit );

	// primary selection for global tracks w/o SPD hits
	AliESDtrackCuts* TrackWOSPDhit = new AliESDtrackCuts ( "TrackWOSPDhit" );
	TrackWOSPDhit->SetClusterRequirementITS ( AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone );
	TrackWOSPDhit->SetMaxDCAToVertexXYPtDep ( ptdepwo.Data() );
	TrackWOSPDhit->SetMaxDCAToVertexZ ( 0.5 );
	TrackWOSPDhit->SetMaxChi2TPCConstrainedGlobal ( 36 );

	fTrackCuts->Add ( TrackWOSPDhit );

	// primary selection for ITS pure tracks with SPD hits
	AliESDtrackCuts* TrackWSPDhitITSpure = new AliESDtrackCuts ( "TrackWSPDhitITSpure" );
	TrackWSPDhit->SetClusterRequirementITS ( AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny );
	TrackWSPDhit->SetMaxDCAToVertexXYPtDep ( ptdepITS.Data() );
	TrackWSPDhit->SetMaxDCAToVertexZ ( 0.5 );

	fTrackCuts->Add ( TrackWSPDhitITSpure );

	// primary selection for ITS pure tracks w/o SPD hits
	AliESDtrackCuts* TrackWOSPDhitITSpure = new AliESDtrackCuts ( "TrackWOSPDhitITSpure" );
	TrackWOSPDhit->SetClusterRequirementITS ( AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone );
	TrackWOSPDhit->SetMaxDCAToVertexXYPtDep ( ptdepITSwo.Data() );
	TrackWOSPDhit->SetMaxDCAToVertexZ ( 0.5 );

	fTrackCuts->Add ( TrackWOSPDhitITSpure );

	delete defaultcuts;
	delete defaultcutsITSpure;

}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::CreateOutput() {
	if ( fOutput ) {
		AliError ( "fOutput already initialised." );
		return;
	}

	fOutput = new TList();
	fOutput->SetOwner();

	for ( Int_t i = 1; i <= fNOutput; i++ ) {
		TList* list = new TList();
		list->SetOwner();
		fOutput->Add ( list );
	}

// 	if ( !fCollectedEvents ) {
// 		fCollectedEvents = new TList();
// 		fCollectedEvents->SetOwner();
// 		fCollectedEvents->SetName ( "CollectedEvents" );
// 	}
// 
// 	( dynamic_cast<TList*> ( fOutput->At ( 0 ) ) )->Add ( fCollectedEvents );

	TH1::SetDefaultSumw2 ( kTRUE );
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::PostDataLists() {
	for ( Int_t i = 1; i <= fNOutput; i++ ) {
		PostData ( i, ( ( TList* ) fOutput )->At ( i - 1 ) );
	}
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::UserCreateOutputObjects() {
	CreateOutput();

	if ( fExclusiveSelectionList ) {
		if ( !fExclusiveSelectionList->IsEmpty() ) {
			TIter nextE ( fExclusiveSelectionList );

			while ( TObject* selection = nextE() ) {
				if ( AliMultiplicityEventSelection* eventSelection = dynamic_cast<AliMultiplicityEventSelection*> ( selection ) ) {
					eventSelection->CreateHistograms();
					Int_t index = fOutputMap[eventSelection->GetName()] - 1;
					( ( TList* ) fOutput->At ( index ) )->Add ( eventSelection );
				}
			}
		}
	} else {
		AliWarning ( "Empty exclusive selections list" );
	}

	if ( fAnalysisSelectionList ) {
		if ( !fAnalysisSelectionList->IsEmpty() ) {
			TIter nextA ( fAnalysisSelectionList );

			while ( TObject* selection = nextA() ) {
				if ( AliMultiplicityEventSelection* eventSelection = dynamic_cast<AliMultiplicityEventSelection*> ( selection ) ) {
					eventSelection->CreateHistograms();
					Int_t index = fOutputMap[eventSelection->GetName()] - 1;
					( ( TList* ) fOutput->At ( index ) )->Add ( eventSelection );
				}
			}
		}
	} else {
		AliWarning ( "Empty analysis selections list" );
	}

	if ( !fTrackCuts ) {
		AliWarning ( "No track cut defined, creating default" );
		CreateDefaultTrackCuts();
	}

	nEvents = new TH1D ( "nEvents", "Event statistics", 6, -0.5, 5.5 );
	nEvents->GetXaxis()->SetBinLabel ( 1, "Input" );
	nEvents->GetXaxis()->SetBinLabel ( 2, "With vertex" );
	nEvents->GetXaxis()->SetBinLabel ( 3, "Not excluded" );
	nEvents->GetXaxis()->SetBinLabel ( 4, "Not excluded + with vertex" );
	nEvents->GetXaxis()->SetBinLabel ( 5, "MBOR" );
	nEvents->GetXaxis()->SetBinLabel ( 6, "V0AND" );
	( ( TList* ) fOutput->At ( 0 ) )->Add ( nEvents );

	hVertex = new TH2D ( "hVertex", "Vertices", 5, -0.5, 4.5, 2, -0.5, 1.5 );
	hVertex->GetXaxis()->SetBinLabel ( 1, "No vertex" );
	hVertex->GetXaxis()->SetBinLabel ( 2, "SPD only" );
	hVertex->GetXaxis()->SetBinLabel ( 3, "global only" );
	hVertex->GetXaxis()->SetBinLabel ( 4, "inconsistent" );
	hVertex->GetXaxis()->SetBinLabel ( 5, "consistent" );
	
	hVertex->GetYaxis()->SetBinLabel ( 1, "no pileup" );
	hVertex->GetYaxis()->SetBinLabel ( 2, "tagged pileup" );
	( ( TList* ) fOutput->At ( 0 ) )->Add ( hVertex );

	PostDataLists();
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::UserExec ( Option_t* ) {
	if ( UseMC() && !MCEvent() ) {
		AliError ( "Cannot load MC event even though it is requested, stopping." );
		return;
	}

// 	MCEvent()->GetPrimaryVertex()->Print();

	if ( !ESDEvent() ) {
		AliError ( "Cannot load ESD event, stopping." );
		return;
	}

	nEvents->Fill ( "Input", 1 );

	SetEventBits ( AliMultiplicityHelper::Classify ( ESDEvent(), ( UseMC() ? MCEvent() : 0 ), GetConsistencyCut() ) );

	if ( TestEventBit ( AliMultiplicityHelper::kMBOR ) ) {
		nEvents->Fill ( "MBOR", 1 );
// 		AliInfo ( "MBOR event" );
	}

	if ( TestEventBit ( AliMultiplicityHelper::kV0AND ) ) {
		nEvents->Fill ( "V0AND", 1 );
// 		AliInfo ( "V0AND event" );
	}

	Bool_t excludeESDEvent = kFALSE;
	Bool_t hasVertex = TestEventBit ( AliMultiplicityHelper::kVSPD ) || TestEventBit ( AliMultiplicityHelper::kVGlobal );
	Bool_t excludeMCEvent = !UseMC();

	if ( !hasVertex ) {
		hVertex->Fill ( "No vertex", TestEventBit( AliMultiplicityHelper::kSPDPileUp ) ? "tagged pileup" : "no pileup" , 1 );
// 		if ( TestEventBit( AliMultiplicityHelper::kVNoSPDDC ) ) hVertex->Fill( "SPD dispersion cut", 1 );
	}

	if ( TestEventBit ( AliMultiplicityHelper::kVSPD ) && !TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) {
// 		if ( TestEventBit ( AliMultiplicityHelper::kVSPD1D ) ) 
		hVertex->Fill ( "SPD only", TestEventBit( AliMultiplicityHelper::kSPDPileUp ) ? "tagged pileup" : "no pileup", 1 );

// 		if ( TestEventBit ( AliMultiplicityHelper::kVSPD3D ) ) hVertex->Fill ( "SPD3D", 1 );
	}

	if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) && !TestEventBit ( AliMultiplicityHelper::kVSPD ) ) {
		hVertex->Fill ( "global only", TestEventBit( AliMultiplicityHelper::kSPDPileUp ) ? "tagged pileup" : "no pileup", 1 );
	}

	if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) && TestEventBit ( AliMultiplicityHelper::kVSPD ) ) {
		if ( TestEventBit ( AliMultiplicityHelper::kVInconsistent ) ) {
			hVertex->Fill ( "inconsistent", TestEventBit( AliMultiplicityHelper::kSPDPileUp ) ? "tagged pileup" : "no pileup", 1 );
	// 		AliInfoF ( "inconsistent: %f", hVertex->GetBinContent(5) );
		}

		if ( TestEventBit ( AliMultiplicityHelper::kVConsistent ) ) {
			hVertex->Fill ( "consistent", TestEventBit( AliMultiplicityHelper::kSPDPileUp ) ? "tagged pileup" : "no pileup", 1 );
	// 		AliInfoF ( "consistent: %f", hVertex->GetBinContent(6) );
		}
	}

	if ( fExclusiveSelectionList ) {
		TIter nextE ( fExclusiveSelectionList );

		while ( TObject* selection = nextE() ) {
			if ( AliMultiplicityEventSelection* eventSelection = dynamic_cast<AliMultiplicityEventSelection*> ( selection ) ) {
				eventSelection->ResetAccepted();

				if ( !eventSelection->AcceptESDEvent ( ESDEvent(), GetEventBits() ) )
					excludeESDEvent = kTRUE;

				if ( UseMC() && eventSelection->GetNeedsMC() && !eventSelection->AcceptMCEvent ( MCEvent(), GetEventBits() ) )
					excludeMCEvent = kTRUE;

				eventSelection->UpdateStats();
// 				eventSelection->SetCollectEvent ( kFALSE );
			}
		}
	} else {
		AliDebug ( 1, "No exclusive selections found!" );
	}

	if ( !excludeESDEvent )
		nEvents->Fill ( "Not excluded", 1 );

	if ( hasVertex ) {
		nEvents->Fill ( "With vertex", 1 );

		if ( !excludeESDEvent )
			nEvents->Fill ( "Not excluded + with vertex", 1 );
	}

	if ( !fAnalysisSelectionList || fAnalysisSelectionList->IsEmpty() ) {
		AliWarning ( "Analysis selection list is empty or undefined." );
		PostDataLists();
		return;
	}

	excludeESDEvent = excludeESDEvent || !hasVertex;

	Bool_t skipESDEvent = kTRUE;
	Bool_t skipMCEvent = kTRUE;

	if ( fAnalysisSelectionList ) {
		TIter nextA ( fAnalysisSelectionList );

		while ( TObject* selection = nextA() ) {
			if ( AliMultiplicityEventSelection* eventSelection = dynamic_cast<AliMultiplicityEventSelection*> ( selection ) ) {
// 				eventSelection->SetCollectEvent ( kFALSE );
				eventSelection->ResetAccepted();

				if ( eventSelection->GetCorrelate() ) {
					if ( !excludeMCEvent && UseMC() && eventSelection->GetNeedsMC() ) {
						if ( eventSelection->AcceptESDandMCEvent ( ESDEvent(), MCEvent(), excludeESDEvent, GetEventBits() ) ) {
							skipMCEvent = kFALSE;
							skipESDEvent = excludeESDEvent;
						}
					}
				} else {
					if ( !excludeESDEvent ) {
						Bool_t accept = eventSelection->AcceptESDEvent ( ESDEvent(), GetEventBits() );

						if ( accept && skipESDEvent )
							skipESDEvent = kFALSE;
					}

					if ( UseMC() && !excludeMCEvent && eventSelection->GetNeedsMC() ) {
						Bool_t accept = eventSelection->AcceptMCEvent ( MCEvent(), GetEventBits() );

						if ( accept && skipMCEvent )
							skipMCEvent = kFALSE;
					}
				}

				eventSelection->UpdateStats();
			}
		}
	}

	if ( ( excludeESDEvent && excludeMCEvent ) || ( skipESDEvent && skipMCEvent ) ) {
		AliWarning ( "All analysis/exclusive selections reject the event, skipping iteration." );
		PostDataLists();
		return;
	}

	Multiplicity ( excludeESDEvent || skipESDEvent, excludeMCEvent || skipMCEvent );

// 	Bool_t collected = kFALSE;
	
	if ( fAnalysisSelectionList ) {
		TIter nextA ( fAnalysisSelectionList );

		while ( TObject* selection = nextA() ) {
			if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
				if ( ( eventSelection->IsCurrentESDEventAccepted() || eventSelection->IsPreselected() ) || ( UseMC() && eventSelection->GetNeedsMC() && eventSelection->IsCurrentMCEventAccepted() ) ) {
// 					AliInfoF ( "%s", eventSelection->GetName() );
					eventSelection->Conclude();

// 					if ( eventSelection->CollectEvent() ) {
// // 						if ( !collected ) {
// // 							CollectEvent ( eventSelection->GetName() );
// 							break;
// // 							collected
// // 						}
// 					}
				}
			}
		}
	}

	PostDataLists();
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::Multiplicity ( Bool_t excludeESDEvent, Bool_t excludeMCEvent ) {
// 	AliInfo ( "Multiplicity loop." );
	if ( !fAnalysisSelectionList ) {
		AliError ( "Multiplicity function called with no analysis selection list, aborting." );
		return;
	}

	//************************************************************************************************************
	//*************************************************************************************************************
	//generator level multiplicity (if available)
	if ( UseMC() && !excludeMCEvent ) {
		AliStack* stack = MCEvent()->Stack();

		if ( !stack ) {
			AliError ( "No stack can be loaded!" );
			return;
		} else {
			Int_t nPrimaries  = stack->GetNprimary(); //get number of particles

			for ( Int_t iParticle = 0; iParticle < nPrimaries; ++iParticle ) { //loop over particles
				if ( stack->IsPhysicalPrimary ( iParticle ) ) { // if it is a primary
					TParticle* particle = stack->Particle ( iParticle );

					if ( !particle ) {
						AliError ( "Cannot get particle info." );
						continue;
					}

					TIter nextA ( fAnalysisSelectionList );

					while ( TObject* selection = nextA() ) {
						if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
// 							if ( eventSelection->GetNeedsMC() && eventSelection->IsCurrentMCEventAccepted() )
								eventSelection->AcceptParticle ( particle, iParticle );
						}
					}
				}
			}
		}
	}

	//*************************************************************************************************************
	//estimated multiplicity


	if ( !excludeESDEvent ) {
		// bit mask for ESD tracks, to check if multiple tracklets are associated to it
		const Int_t nESDTracks = ESDEvent()->GetNumberOfTracks();
// 		Int_t highestID = 0;
//
// 		for ( Int_t iTrack = 0; iTrack < ESDEvent()->GetNumberOfTracks(); iTrack++ ) {
// 			if ( ESDEvent()->GetTrack ( iTrack )->GetLabel() > highestID ) highestID = ESDEvent()->GetTrack ( iTrack )->GetLabel();
// 		}
		const Int_t maxid = nESDTracks + 1; // used to define Bool_t array for check multiple associations of tracklets to one track. array starts at 0.
// 		Bool_t globalBits[maxid], pureITSBits[maxid];
		TBits globalBits ( maxid ), pureITSBits ( maxid );
		// flags for secondary and rejected tracks
		const Int_t kRejBit = BIT ( 15 ); // set this bit in ESD tracks if it is rejected by a cut
		const Int_t kSecBit = BIT ( 16 ); // set this bit in ESD tracks if it is secondary according to a cut

		for ( Int_t itracks = 0; itracks < nESDTracks; itracks++ ) {
			ESDEvent()->GetTrack ( itracks )->ResetBit ( kSecBit | kRejBit ); //reset bits used for flagging secondaries and rejected tracks in case they were changed before this analysis
		}

// 		if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) {
			//aux vars
// 			for ( Int_t i = 0; i < maxid; i++ ) { // set all bools to kFALSE
// 				globalBits[i] = kFALSE;
// 				pureITSBits[i] = kFALSE;
// 			}
			//The Counting
			// tracks

			//*******************************************************************************************************
			// get multiplicity from ESD tracks
			for ( Int_t iTrack = 0; iTrack < nESDTracks; iTrack++ ) { // flag the tracks
				AliESDtrack* track = ESDEvent()->GetTrack ( iTrack );
				// if track is a secondary from a V0, flag as a secondary
// 				if ( GetUseMultInV0() && track->IsOn ( AliESDtrack::kMultInV0 ) ) {
// 					track->SetBit ( kSecBit ); // secondary from V0
// 					continue;
// 				}

				// check tracks with ITS part
				//*******************************************************************************************************
				if ( track->IsOn ( AliESDtrack::kITSin ) && !track->IsOn ( AliESDtrack::kITSpureSA ) ) { // track has ITS part but is not an ITS_SA
					//*******************************************************************************************************
					// TPC+ITS
					if ( track->IsOn ( AliESDtrack::kTPCin ) ) {  // Global track, has ITS and TPC contributions
						if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "GlobalTrackQuality" ) )->AcceptTrack ( track ) ) {  // good ITSTPC track
							if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) { //have global vertex - check DCA
								if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWSPDhit" ) )->AcceptTrack ( track )
									|| dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWOSPDhit" ) )->AcceptTrack ( track ) ) {
									TIter nextA ( fAnalysisSelectionList );

									while ( TObject* selection = nextA() ) {
										if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
											if ( eventSelection->IsCurrentESDEventAccepted() ) {
												eventSelection->AcceptESDTrack ( track, AliMultiplicityEventSelection::kITSTPC, kFALSE );

												if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
													eventSelection->AcceptESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC, kFALSE );
												}
											}
										}
									}
									globalBits.SetBitNumber ( iTrack );
								} else
									track->SetBit ( kSecBit ); // large DCA -> secondary, don't count either track not associated tracklet
							} else {
								TIter nextA ( fAnalysisSelectionList );
								
								while ( TObject* selection = nextA() ) {
									if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
										if ( eventSelection->IsCurrentESDEventAccepted() ) {
											eventSelection->AcceptSecondaryESDTrack ( track, AliMultiplicityEventSelection::kITSTPC, kFALSE );
											
											if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
												eventSelection->AcceptSecondaryESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC, kFALSE );
											}
										}
									}
								}
// 								globalBits.SetBitNumber ( iTrack );
							}
						} else
							track->SetBit ( kRejBit ); // bad quality, don't count the track, but may count tracklet if associated
					}
					//*******************************************************************************************************
					// ITS complementary
					else if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "ITSpureQuality" ) )->AcceptTrack ( track ) ) {  // good ITS complementary track
						if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) { //have global vertex - check DCA
							if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWSPDhitITSpure" ) )->AcceptTrack ( track )
								|| dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWOSPDhitITSpure" ) )->AcceptTrack ( track ) ) {
								TIter nextA ( fAnalysisSelectionList );

								while ( TObject* selection = nextA() ) {
									if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
										if ( eventSelection->IsCurrentESDEventAccepted() ) {
											eventSelection->AcceptESDTrack ( track, AliMultiplicityEventSelection::kITSTPC, kTRUE );

											if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
												eventSelection->AcceptESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC, kTRUE );
											}
										}
									}
								}

	// 							globalBits[iTrack] = kTRUE;
								globalBits.SetBitNumber ( iTrack );
							} else
								track->SetBit ( kSecBit ); // large DCA -> secondary, don't count either track not associated tracklet
						} else {
							TIter nextA ( fAnalysisSelectionList );
							
							while ( TObject* selection = nextA() ) {
								if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
									if ( eventSelection->IsCurrentESDEventAccepted() ) {
										eventSelection->AcceptSecondaryESDTrack ( track, AliMultiplicityEventSelection::kITSTPC, kTRUE );
										
										if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
											eventSelection->AcceptSecondaryESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC, kTRUE );
										}
									}
								}
							}
// 							globalBits.SetBitNumber ( iTrack );
						}
					} else
						track->SetBit ( kRejBit ); // bad quality, don't count the track, but may count tracklet if associated
				}

				//*******************************************************************************************************
				// check tracks from ITS_SA_PURE
				if ( track->IsOn ( AliESDtrack::kITSin ) && track->IsOn ( AliESDtrack::kITSpureSA ) ) {
					if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "ITSpureQuality" ) )->AcceptTrack ( track ) ) { // good ITSSA track
						if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) { //have global vertex - check DCA
							if ( dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWSPDhitITSpure" ) )->AcceptTrack ( track )
								|| dynamic_cast<AliESDtrackCuts*> ( fTrackCuts->FindObject ( "TrackWOSPDhitITSpure" ) )->AcceptTrack ( track ) ) {
								TIter nextA ( fAnalysisSelectionList );

								while ( TObject* selection = nextA() ) {
									if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
										if ( eventSelection->IsCurrentESDEventAccepted() ) {
											eventSelection->AcceptESDTrack ( track, AliMultiplicityEventSelection::kITSSA, kFALSE );

											if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
												eventSelection->AcceptESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSSA, kFALSE );
											}
										}
									}
								}

	// 							pureITSBits[iTrack] = kTRUE;
								pureITSBits.SetBitNumber ( iTrack );
							} else
								track->SetBit ( kSecBit );  //WTF FIXME
						} else {
							TIter nextA ( fAnalysisSelectionList );
							
							while ( TObject* selection = nextA() ) {
								if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
									if ( eventSelection->IsCurrentESDEventAccepted() ) {
										eventSelection->AcceptSecondaryESDTrack ( track, AliMultiplicityEventSelection::kITSSA, kFALSE );
										
										if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
											eventSelection->AcceptSecondaryESDTrackMC ( track, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSSA, kFALSE );
										}
									}
								}
							}
// 							pureITSBits.SetBitNumber ( iTrack );
						}
					} else
						track->SetBit ( kRejBit );
				}
			}//ESD tracks counted
// 		}

		//tracklets
		if ( TestEventBit ( AliMultiplicityHelper::kVSPD ) ) {
			const AliMultiplicity* spdmult = ESDEvent()->GetMultiplicity();    // spd multiplicity object

			for ( Int_t iTracklet = 0; iTracklet < spdmult->GetNumberOfTracklets(); ++iTracklet ) {
				TIter nextAe ( fAnalysisSelectionList );

				while ( TObject* selection = nextAe() ) {
					if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
// 							AliInfoF ( "%s selection called", selection->GetTitle());
						if ( eventSelection->IsCurrentESDEventAccepted() ) {
							eventSelection->AcceptTracklet ( spdmult, iTracklet, AliMultiplicityEventSelection::kSPD );

							if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
								eventSelection->AcceptTrackletMC ( spdmult, iTracklet, MCEvent()->Stack(), AliMultiplicityEventSelection::kSPD );
							}
						}
					}
				}

				if ( TestEventBit ( AliMultiplicityHelper::kVGlobal ) ) {
					Int_t id1, id2, id3, id4;
					spdmult->GetTrackletTrackIDs ( iTracklet, 0, id1, id2 ); // references for eventual Global/ITS_SA tracks
					spdmult->GetTrackletTrackIDs ( iTracklet, 1, id3, id4 ); // references for eventual ITS_SA_pure tracks

					// are both clusters from the same tracks? If not, skip the tracklet (shouldn't change things much)
					if ( ( id1 != id2 && id1 >= 0 && id2 >= 0 ) || ( id3 != id4 && id3 >= 0 && id4 >= 0 ) )
						continue;

					//referenced track
					//at this point we either have id1 = id2 (id3 = id4) or only one of the ids pair is -1
					// id1>=0, id2>=0, id1=id2	: tracklet has associated track
					// id1>=0, id2 = -1		: 1st layer cluster has associated track
					// id1=-1, id2>=0		: 2nd layer cluster has associated track
					// id1=-1, id2=-1		: tracklet has no associated track

					Int_t bUsedInGlobal ( -1 );

					if ( id1 != -1 )
						bUsedInGlobal = globalBits.TestBitNumber ( id1 ) ? id1 : -1;
					else if ( id2 != -1 )
						bUsedInGlobal = globalBits.TestBitNumber ( id2 ) ? id2 : -1; // has associated global track been associated to a previous tracklet?

					Int_t bUsedInPureITS ( -1 );

					if ( id3 != -1 )
						bUsedInPureITS = pureITSBits.TestBitNumber ( id3 ) ? id3 : -1;
					else if ( id4 != -1 )
						bUsedInPureITS = pureITSBits.TestBitNumber ( id4 ) ? id4 : -1; // has associated pure ITS track been associated to a previous tracklet?

					AliESDtrack* tr_global = bUsedInGlobal >= 0 ? ESDEvent()->GetTrack ( bUsedInGlobal ) : 0;
					AliESDtrack* tr_itssa = bUsedInPureITS >= 0 ? ESDEvent()->GetTrack ( bUsedInPureITS ) : 0;

					//*******************************************************************************************************
					// count tracklets towards global+complimentary tracks
					if ( ( tr_global && !tr_global->TestBit ( kSecBit ) ) && ( tr_global &&  tr_global->TestBit ( kRejBit ) ) ) {  // count tracklet as bad quality track
						TIter nextA ( fAnalysisSelectionList );

						while ( TObject* selection = nextA() ) {
							if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
								if ( eventSelection->IsCurrentESDEventAccepted() ) {
									eventSelection->AcceptESDTrack ( tr_global, AliMultiplicityEventSelection::kITSTPC, kTRUE );

									if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
										eventSelection->AcceptESDTrackMC ( tr_global, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC, kTRUE );
									}
								}
							}
						}

						globalBits.SetBitNumber ( bUsedInGlobal ); // mark global track linked to this tracklet as "associated"
					};

					if ( bUsedInGlobal < 0 ) {//no associated track
						TIter nextA ( fAnalysisSelectionList );

						while ( TObject* selection = nextA() ) {
							if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
								if ( eventSelection->IsCurrentESDEventAccepted() ) {
									eventSelection->AcceptTracklet ( spdmult, iTracklet, AliMultiplicityEventSelection::kITSTPC );

									if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
										eventSelection->AcceptTrackletMC ( spdmult, iTracklet, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSTPC );
									}
								}
							}
						}
					}

					//*******************************************************************************************************
					// count tracklets towards ITS_SA_pure tracks
					if ( ( tr_itssa && !tr_itssa->TestBit ( kSecBit ) ) && ( tr_itssa &&  tr_itssa->TestBit ( kRejBit ) ) ) {  // count tracklet as bad quality track
						TIter nextA ( fAnalysisSelectionList );

						while ( TObject* selection = nextA() ) {
							if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
								if ( eventSelection->IsCurrentESDEventAccepted() ) {
									eventSelection->AcceptESDTrack ( tr_itssa, AliMultiplicityEventSelection::kITSSA, kTRUE );

									if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
										eventSelection->AcceptESDTrackMC ( tr_itssa, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSSA, kTRUE );
									}
								}
							}
						}

						pureITSBits.SetBitNumber ( bUsedInPureITS );
					};

					if ( bUsedInPureITS < 0 ) { //no associated track
						TIter nextA ( fAnalysisSelectionList );

						while ( TObject* selection = nextA() ) {
							if ( AliMultiplicityAnalysisSelection* eventSelection = dynamic_cast<AliMultiplicityAnalysisSelection*> ( selection ) ) {
								if ( eventSelection->IsCurrentESDEventAccepted() ) {
									eventSelection->AcceptTracklet ( spdmult, iTracklet, AliMultiplicityEventSelection::kITSSA );

									if ( UseMC() && ( !excludeMCEvent && !excludeESDEvent ) && eventSelection->GetCorrelateTracks() ) {
										eventSelection->AcceptTrackletMC ( spdmult, iTracklet, MCEvent()->Stack(), AliMultiplicityEventSelection::kITSSA );
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::Terminate ( Option_t* ) {
	if ( fOutput ) {
		AliWarning ( "Deleting output list to load a new one." );
		delete fOutput;
	}

	fOutput = new TList();

	for ( Int_t i = 1; i <= fNOutput; i++ ) {
		TList* list = dynamic_cast<TList*> ( GetOutputData ( i ) );

		if ( !list ) {
			AliErrorF ( "Unable to get output from slot %d.", i );
			fOutput->Add ( new TList() );
		} else {
// 			AliInfoF ( "Output: %d", i );
			fOutput->Add ( list );
		}
	}

	TIter nextSlot ( fOutput );

	while ( TList* output = dynamic_cast<TList*> ( nextSlot() ) ) {
		if ( output->GetEntries() > 0 ) {
// 			AliDebugF ( 2, "Reading from slot %d.", slot );
			TIter nextS ( output );

			while ( TObject* selection = nextS() ) {
				if ( AliMultiplicityEventSelection* eventSelection = dynamic_cast<AliMultiplicityEventSelection*> ( selection ) ) {
// 					AliDebugF ( 1, "Got selection from slot %d: %s.", slot, eventSelection->GetName() );

					if ( eventSelection->GetSaveHistograms() ) {
						eventSelection->SaveHistograms();
					}

					if ( !eventSelection->GetBatchMode() )
						eventSelection->Result();
				}
			}
		}
	}	
}
//_____________________________________________________________________________________________________________________________________________
void AliAnalysisMultiplicityTask::DumpLists() {
	if ( fExclusiveSelectionList )
		fExclusiveSelectionList->Dump();

	if ( fAnalysisSelectionList )
		fAnalysisSelectionList->Dump();
}

void AliAnalysisMultiplicityTask::CollectEvent ( const char* origin ) {
	const char* filename = CurrentFileName();
	TObjString* s = new TObjString ( TString::Format ( "From: %s; file = %s; event# = %d", origin, filename, ESDEvent()->GetEventNumberInFile() ).Data() );
	AliInfoF ( "COLLECT: %s", s->GetString().Data() );
	fCollectedEvents->Add ( s );
}
//_____________________________________________________________________________________________________________________________________________
// Bool_t AliAnalysisMultiplicityTask::IsCocktail() {
// 	AliGenEventHeader* h = AliMultiplicityHelper::GetGenEventHeader(MCEvent());
// 	AliGenCocktailEventHeader* hh = dynamic_cast<AliGenCocktailEventHeader*>(h);
// 	if ( hh ) {
// 		fIsCocktail = kTRUE;
// 	}
// 	return fIsCocktail;
// }
