#include "AliMultiplicityEfficiencySelection.h"
#include "TLegend.h"
#include <TAttMarker.h>
#include <TMathBase.h>
#include <TPad.h>
#include "AliMultiplicityHelper.h"
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

ClassImp ( AliMultiplicityEfficiencySelection )

void AliMultiplicityEfficiencySelection::CreateSelectionHistograms() {
	TString names[3] = {"!SD && !DD", "SD", "DD"};
	TString methods[3] = {"ITSTPC", "ITSSA", "SPD"};

	for ( Int_t i = 0; i < 3; i++ ) {
		hEventsGenerated[i] = new TH1D ( Form ( "hEventsGenerated_%d", i ), Form ( "Generator-level multiplicity (%s)", names[i].Data() ), 161, -0.5, 160.5 );
		hEventsSelected[i] = new TH1D ( Form ( "hEventsSelected_%d", i ), Form ( "Selected events per true multiplicity bin (%s)", names[i].Data() ), 161, -0.5, 160.5 );

// 		hEventsGeneratedN[i] = new TH2D ( Form ( "hEventsGeneratedN_%d", i ), Form ( "Generator-level multiplicity (%s) vs. Z", names[i].Data() ), 161, -0.5, 160.5, 2001, -20.02, 20.02 );
	}

	AliMultiplicitySelection::CreateSelectionHistograms();

	for ( Int_t i = 0; i < kSPD; ++i ) {
		for ( Int_t k = 0; k < 2; ++k ) {
			pStrangeSurvival[i * 2 + k] = new TProfile ( Form ( "pStrangeSurvival_%d_%d", i, k ), Form ( "Average strangeness survival rate per event vs. multiplicity (%s, %d)", methods[i].Data(), k ), 161, -0.5, 160.5, 0, 5, "" );
		}

		hSurvivedEta[i] = new TH1D ( Form ( "hSurvivedEta_%d", i ), Form ( "Strangeness survival vs. #eta (%s)", methods[i].Data() ), 100, -2.5, 2.5 );
	}

// 	for (Int_t i = 0; i < 2; ++i) {
// 		hMZZ[i] = new TH1D(Form("hMZZ_%d",i),"",161,-0.5,160.5);
// 		hNZc[i] = new TH1D(Form("hNZc_%d",i),"",161,-0.5,160.5);
// 	}
//
// 	h0bin = new TH1D("h0bin","",161,-0.5,160.5);
// 	h0binZ = new TH1D("h0binZ","",2001,-20.02,20.02);
// 	h0binE = new TH1D("h0binE","",201,-2.02,2.02);

	hGeneratedEta = new TH1D ( "hGeneratedEta", "Strange particles generated vs. #eta", 100, -2.5, 2.5 );
}

Long64_t AliMultiplicityEfficiencySelection::Merge ( const TCollection* mergelist ) {
	if ( !mergelist ) {
		AliError ( "NULL merge list." );
		return 0;
	}

	if ( mergelist->IsEmpty() ) {
		AliWarning ( "Empty merge list." );
		return 1;
	}

	AliMultiplicitySelection::Merge ( mergelist );

	TList ghistograms[3];
	TList shistograms[3];
	TList psurv[3 * 2];
	TList etasurv[3];
	TList etagen;
// 	TList hnz[3];
// 	TList lmzz[2];
// 	TList l0;
// 	TList l0z;
// 	TList l0e;
// 	TList lhnzc[2];

	Int_t count = 0;

	TIter nextE ( mergelist );

	while ( AliMultiplicityEfficiencySelection* selection = dynamic_cast<AliMultiplicityEfficiencySelection*> ( nextE() ) ) {
		if ( !selection ) continue;

		for ( Int_t i = 0; i < 3; i++ ) {
			ghistograms[i].Add ( selection->hEventsGenerated[i] );
			shistograms[i].Add ( selection->hEventsSelected[i] );
			nGenerated[i] += selection->nGenerated[i];
			nSelected[i] += selection->nSelected[i];
// 			hnz[i].Add( selection->hEventsGeneratedN[i] );
		}

		for ( Int_t i = 0; i < kSPD; ++i ) {
			for ( Int_t k = 0; k < 2; ++k ) {
				psurv[i * 2 + k].Add ( selection->pStrangeSurvival[i * 2 + k] );
			}

			etasurv[i].Add ( selection->hSurvivedEta[i] );
		}

// 		for (Int_t i = 0; i < 2; ++i) {
// 			lmzz[i].Add( selection->hMZZ[i] );
// 			lhnzc[i].Add( selection->hNZc[i] );
// 		}

// 		l0.Add( selection->h0bin );
// 		l0z.Add( selection->h0binZ );
// 		l0e.Add( selection->h0binE );

		etagen.Add ( selection->hGeneratedEta );
		count++;
	}

	for ( Int_t i = 0; i < 3; i++ ) {
		hEventsGenerated[i]->Merge ( &ghistograms[i] );
		hEventsSelected[i]->Merge ( &shistograms[i] );
// 		hEventsGeneratedN[i]->Merge ( &hnz[i] );
	}

	for ( Int_t i = 0; i < kSPD; ++i ) {
		for ( Int_t k = 0; k < 2; ++k ) {
			pStrangeSurvival[i * 2 + k]->Merge ( &psurv[i * 2 + k] );
		}

		hSurvivedEta[i]->Merge ( &etasurv[i] );
	}

// 	for (Int_t i = 0; i < 2; ++i) {
// 		hMZZ[i]->Merge( &lmzz[i] );
// 		hNZc[i]->Merge( &lhnzc[i] );
// 	}
// 	h0bin->Merge(&l0);
// 	h0binZ->Merge(&l0z);
// 	h0binE->Merge(&l0e);

	hGeneratedEta->Merge ( &etagen );

	return ( Long64_t ) count + 1;
}

void AliMultiplicityEfficiencySelection::Conclude( ) {
	if ( fCurrentESDEventAccepted ) AliMultiplicitySelection::Conclude();

	if ( ( currentProcessType != -1 ) ) {
		if ( fCurrentMCEventAccepted ) hEventsGenerated[currentProcessType]->Fill ( nMCparticles );

		if ( fCurrentESDEventAccepted ) hEventsSelected[currentProcessType]->Fill ( nMCparticles );

// 		hEventsGeneratedN[currentProcessType]->Fill( nMCparticles, fZcacheGen );
	}

// 	if ( fCurrentESDEventAccepted /*IsPreselected() && fMCVInCuts*/ ) {
// 		if ( nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements == 0 ) {
// 			h0bin->Fill( nMCparticles );
// 			h0binZ->Fill( fZcacheGen );
// 			TIter next ( particles );
// 			TObject* o = 0;
// 			while ( ( o = next() ) ) {
// 				TParticle* p = dynamic_cast<TParticle*>(o);
// 				if ( p ) {
// 					h0binE->Fill(p->Eta());
// 				}
// 			}
// 		}
// 		if ( fMCVInCuts && !(fGlobalInCuts || fSPDInCuts) ) {
// 			hMZZ[0]->Fill(nMCparticles);
// 		}
//
// 		if ( !fMCVInCuts && (fGlobalInCuts || fSPDInCuts) ) {
// 			hMZZ[1]->Fill(nMCparticles);
// 		}
// 	}
//
// 	if ( fPreSelected ){
// 		if ( fMCVInCuts ) {
// 			hNZc[0]->Fill(nMCparticles);
// 		}
// 		if ( fGlobalInCuts || fSPDInCuts ) {
// 			hNZc[1]->Fill(nMCparticles);
// 		}
// 	}

	if ( ( nGeneratedS[0] == 0 ) && ( nGeneratedS[1] == 0 ) ) {
		return;
	}

	for ( Int_t i = 0; i < kSPD; ++i ) {
		for ( Int_t k = 0; k < 2; ++k ) {
			if ( fCurrentESDEventAccepted && ( nGeneratedS[k] != 0 ) ) pStrangeSurvival[i * 2 + k]->Fill ( nMCparticles, nSurvived[i * 2 + k] / nGeneratedS[k], 1 );
		}
	}
}

Bool_t AliMultiplicityEfficiencySelection::AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask ) {
	AliMultiplicitySelection::AcceptESDEvent ( esd, mask );
	AliMultiplicitySelection::AcceptMCEvent ( mc, mask );

	UInt_t processes = AliMultiplicityHelper::kNonDiffractive | AliMultiplicityHelper::kSingleDiffractive | AliMultiplicityHelper::kDoubleDiffractive;
	UInt_t check = mask & processes;

	switch ( check ) {
		case AliMultiplicityHelper::kNonDiffractive:
			AliInfo ( "ND" );
			currentProcessType = 0;
			break;

		case AliMultiplicityHelper::kSingleDiffractive:
			AliInfo ( "SD" );
			currentProcessType = 1;
			break;

		case AliMultiplicityHelper::kDoubleDiffractive:
			AliInfo ( "DD" );
			currentProcessType = 2;
			break;

		default:
			AliWarning ( "Cannot determine MC process type" );
	}

	Int_t originalProcessType = currentProcessType;

	if ( bWeightDiffraction && IsCurrentMCEventAccepted() ) {
		Double_t energy = esd->GetESDRun()->GetBeamEnergy();

		if ( esd->GetESDRun()->IsBeamEnergyIsSqrtSHalfGeV() ) energy *= 2;
		
		if ( fDiffMassCut < 0 ) {
			fCurrentMCEventAccepted = AliMultiplicityHelper::CheckDiff ( mc->Stack(), originalProcessType, currentProcessType, energy, hwSD, hwNDDD );
		} else {
			fCurrentMCEventAccepted = AliMultiplicityHelper::CheckDiff ( mc->Stack(), originalProcessType, currentProcessType, energy, fDiffMassCut );
		}
		
		AliInfoF ( "process %d to %d", originalProcessType, currentProcessType );
	}

	if ( fFollowMCselection ) fCurrentESDEventAccepted = !excludeESD && IsPreselected() && fCurrentMCEventAccepted;

	AliInfoF ( "%d - %d", fCurrentMCEventAccepted, fCurrentESDEventAccepted );

	if ( ( currentProcessType != -1 ) && IsCurrentMCEventAccepted() ) {
		nGenerated[currentProcessType] += 1;

		if ( IsPreselected() ) nSelected[currentProcessType] += 1;
	}

	//count strange decay products in given eta range
	AliStack* stack = mc->Stack();
	Int_t nPrim = stack->GetNprimary();
	Int_t nParticles = stack->GetNtrack();

	TParticle* mother = 0;
	Int_t nGeneratedSA = 0;

	//count all charged and relatively stable secondary particles that have strange mother
	for ( Int_t i = nPrim; i < nParticles; ++i ) {

		Int_t mlabel = AliMultiplicityHelper::FindPrimaryMother ( stack, i );

		mother = AliMultiplicityHelper::GetMother ( stack, mlabel );

		if ( mother && ( ( TMath::Abs ( mother->GetPdgCode() ) == 3122 ) || ( mother->GetPdgCode() == 310 ) ) ) {
			if ( ( TMath::Abs ( stack->Particle ( i )->GetPDG()->Charge() ) > 0 ) && ( stack->IsStable ( stack->Particle ( i )->GetPdgCode() ) ) ) {
				if ( TMath::Abs ( stack->Particle ( i )->Eta() ) < GetEtaCut() ) {
					if ( TMath::Abs ( mother->GetPdgCode() ) == 3122 ) nGeneratedS[0] += 1;

					if ( TMath::Abs ( mother->GetPdgCode() ) == 310 ) nGeneratedS[1] += 1;
				}

				hGeneratedEta->Fill ( stack->Particle ( i )->Eta() );
				nGeneratedSA++;
			}
		}
	}

	bIsGenerated = ! ( nGeneratedSA == 0 );

	return fCurrentMCEventAccepted || fCurrentESDEventAccepted;
}

AliMultiplicityEfficiencySelection::AliMultiplicityEfficiencySelection ( const char* name, const char* title ) : AliMultiplicitySelection ( name, title ),
	currentProcessType ( -1 ),
	hGeneratedEta ( 0 ),
	bIsGenerated ( kFALSE ),
	bWeightDiffraction ( kFALSE ),
	hwSD ( 0 ),
	hwNDDD ( 0 )
/*	h0bin(0),
	h0binZ(0),
	h0binE(0)*/ {
	for ( Int_t i = 0; i < 3; i++ ) {
		hEventsGenerated[i] = 0;
// 		hEventsGeneratedN[i] = 0;
		hEventsSelected[i] = 0;
		nSelected[i] = 0;
		nGenerated[i] = 0;
	}

	nGeneratedS[0] = 0;
	nGeneratedS[1] = 0;

	for ( Int_t i = 0; i < kSPD; ++i ) {
		for ( Int_t k = 0; k < 2; ++k ) {
			nSurvived[i * 2 + k] = 0;
			pStrangeSurvival[i * 2 + k] = 0;
		}

		hSurvivedEta[i] = 0;
	}

// 	for (Int_t i = 0; i < 2; ++i) {
// 		hMZZ[i] = 0;
// 		hNZc[i] = 0;
// 	}

	SetNeedsMC ( kTRUE );
	SetCorrelate ( kTRUE );
	SetCollectEvent ( kFALSE );
	SetCorrelateTracks ( kTRUE );
	SetQueryPhysSelection ( kTRUE );
	SetFollowMCselection ( kTRUE );
}

AliMultiplicityEfficiencySelection::~AliMultiplicityEfficiencySelection() {

}

void AliMultiplicityEfficiencySelection::SaveSelectionHistograms() {
	AliMultiplicitySelection::SaveSelectionHistograms();

	for ( Int_t i = 0; i < 3; i++ ) {
		if ( hEventsGenerated[i] ) hEventsGenerated[i]->Write();

		if ( hEventsSelected[i] ) hEventsSelected[i]->Write();
	}


	TH1D* e1 = dynamic_cast<TH1D*> ( GetEfficiency ( 0 ) );

	if ( e1 ) {
		e1->SetName ( "effINEL" );
		e1->Write();
	}

	TH1D* e2 = dynamic_cast<TH1D*> ( GetEfficiency ( 1 ) );

	if ( e2 ) {
		e2->SetName ( "effSD" );
		e2->Write();
	}

	TH1D* e3 = dynamic_cast<TH1D*> ( GetEfficiency ( 2 ) );

	if ( e3 ) {
		e3->SetName ( "effNSD" );
		e3->Write();
	}

	for ( Int_t i = 0; i < kSPD; ++i ) {
		if ( pStrangeSurvival[i] ) pStrangeSurvival[i]->Write();
	}

	TH1D* etarate[kSPD];

	for ( Int_t i = 0; i < kSPD; ++i ) {
		if ( hSurvivedEta[i] ) {
			hSurvivedEta[i]->Write();
			etarate[i] = new TH1D ( *hSurvivedEta[i] );
			etarate[i]->Divide ( hSurvivedEta[i], hGeneratedEta, 1, 1, "b" );
			etarate[i]->SetName ( Form ( "hSurvivedEta_div_%d", i ) );
			etarate[i]->Write();
		}
	}

	if ( hGeneratedEta ) hGeneratedEta->Write();

	TH1D* halpha = GetAlphaSD();

	if ( halpha ) halpha->Write();
}

void AliMultiplicityEfficiencySelection::Result() {
	TCanvas* cnv = GetCanvas();
	cnv->cd()->SetLogy ( 0 );

	TLegend* l = new TLegend ( 0.8, 0.2, 0.9, 0.3 );

	TH1D* numeratorINEL = dynamic_cast<TH1D*> ( GetEfficiency ( 0 ) );
	numeratorINEL->Draw ( "P" );
	l->AddEntry ( numeratorINEL, "INEL", "lep" );

	TH1D* numeratorSD = dynamic_cast<TH1D*> ( GetEfficiency ( 1 ) );
	numeratorSD->Draw ( "P same" );
	l->AddEntry ( numeratorSD, "SD", "lep" );

	TH1D* numeratorDD = dynamic_cast<TH1D*> ( GetEfficiency ( 2 ) );
	numeratorDD->Draw ( "P same" );
	l->AddEntry ( numeratorDD, "DD", "lep" );

	TH1D* numeratorNSD = dynamic_cast<TH1D*> ( GetEfficiency ( 3 ) );
	numeratorNSD->Draw ( "P same" );
	l->AddEntry ( numeratorNSD, "NSD", "lep" );
	l->SetFillColor ( kWhite );

	l->Draw();
}

// void AliMultiplicityEfficiencySelection::ResetAccepted() {
// 	AliMultiplicityEventSelection::ResetAccepted();
// }

TH1D* AliMultiplicityEfficiencySelection::GetGeneratedMultiplicity ( Int_t type ) {
	TH1D* output = 0;

	switch ( type ) {
		case 0:
			output = new TH1D ( *hEventsGenerated[0] );
			output->Add ( hEventsGenerated[1] );
			output->Add ( hEventsGenerated[2] );
			return output;

		case 1:
			output = new TH1D ( *hEventsGenerated[1] );
			return output;

		case 2:
			output = new TH1D ( *hEventsGenerated[2] );
			return output;

		case 3:
			output = new TH1D ( *hEventsGenerated[0] );
			output->Add ( hEventsGenerated[2] );
			return output;
			break;

		default:
			AliWarning ( "Unknown process type" );
			return 0;
			break;
	}
}

TH1D* AliMultiplicityEfficiencySelection::GetSelectedMultiplicity ( Int_t type ) {
	TH1D* output = 0;
	
	switch ( type ) {
		case 0:
			output = new TH1D ( *hEventsSelected[0] );
			output->Add ( hEventsSelected[1] );
			output->Add ( hEventsSelected[2] );
			return output;
			
		case 1:
			output = new TH1D ( *hEventsSelected[1] );
			return output;
			
		case 2:
			output = new TH1D ( *hEventsSelected[2] );
			return output;
			
		case 3:
			output = new TH1D ( *hEventsSelected[0] );
			output->Add ( hEventsSelected[2] );
			return output;
			break;
			
		default:
			AliWarning ( "Unknown process type" );
			return 0;
			break;
	}
}


TH1D* AliMultiplicityEfficiencySelection::GetEfficiency ( Int_t type ) {
	if ( !hEventsGenerated[0] || !hEventsSelected[0] ) {
		AliWarning ( "Histograms are undefined!" );
		return 0;
	}

	TH1D* numeratorINEL = 0;
	TH1D* numeratorSD = 0;
	TH1D* numeratorDD = 0;
	TH1D* numeratorNSD = 0;
	TH1D* numeratorND = 0;
	TH1D* denominatorINEL = 0;
	TH1D* denominatorNSD = 0;

	switch ( type ) {
		case 0:
			denominatorINEL = new TH1D ( *hEventsGenerated[0] );
			denominatorINEL->Add ( hEventsGenerated[1] );
			denominatorINEL->Add ( hEventsGenerated[2] );

			numeratorINEL = new TH1D ( *hEventsSelected[0] );
			numeratorINEL->Add ( hEventsSelected[1] );
			numeratorINEL->Add ( hEventsSelected[2] );

			numeratorINEL->Divide ( numeratorINEL, denominatorINEL, 1, 1, "b" );
			numeratorINEL->SetTitle ( "" );
			numeratorINEL->SetMarkerStyle ( kOpenCircle );
			numeratorINEL->SetMarkerColor ( kRed );
			numeratorINEL->SetAxisRange ( 0, 1.1, "Y" );
			numeratorINEL->SetAxisRange ( -0.5, 12.5, "X" );
			numeratorINEL->SetTitle ( "" );
			return numeratorINEL;

		case 1:
			numeratorSD = new TH1D ( *hEventsSelected[1] );
			numeratorSD->Divide ( numeratorSD, hEventsGenerated[1], 1, 1, "b" );
			numeratorSD->SetTitle ( "" );
			numeratorSD->SetMarkerStyle ( kOpenDiamond );
			numeratorSD->SetMarkerColor ( kBlue );
			numeratorSD->SetAxisRange ( 0, 1.1, "Y" );
			numeratorSD->SetAxisRange ( -0.5, 12.5, "X" );
			numeratorSD->SetTitle ( "" );
			return numeratorSD;

		case 2:
			numeratorDD = new TH1D ( *hEventsSelected[2] );
			numeratorDD->Divide ( numeratorDD, hEventsGenerated[2], 1, 1, "b" );
			numeratorDD->SetTitle ( "" );
			numeratorDD->SetMarkerStyle ( kOpenSquare );
			numeratorDD->SetMarkerColor ( kGreen );
			numeratorDD->SetAxisRange ( 0, 1.1, "Y" );
			numeratorDD->SetAxisRange ( -0.5, 12.5, "X" );
			numeratorDD->SetTitle ( "" );
			return numeratorDD;

		case 3:
			numeratorNSD = new TH1D ( *hEventsSelected[0] );
			numeratorNSD->Add ( hEventsSelected[2] );

			denominatorNSD = new TH1D ( *hEventsGenerated[0] );
			denominatorNSD->Add ( hEventsGenerated[2] );

			numeratorNSD->Divide ( numeratorNSD, denominatorNSD, 1, 1, "b" );
			numeratorNSD->SetTitle ( "" );
			numeratorNSD->SetMarkerStyle ( 24 );
			numeratorNSD->SetMarkerColor ( kCyan );
			numeratorNSD->SetAxisRange ( 0, 1.1, "Y" );
			numeratorNSD->SetAxisRange ( -0.5, 12.5, "X" );
			
			return numeratorNSD;
		case 4:
			numeratorND = new TH1D ( *hEventsSelected[0] );
			numeratorND->Divide ( numeratorND, hEventsGenerated[0], 1, 1, "b" );
			numeratorND->SetTitle("");
			numeratorND->SetMarkerStyle(kFullCircle);
			numeratorND->SetMarkerColor( kViolet );
			numeratorND->SetAxisRange( 0, 1.1, "Y" );
			numeratorND->SetAxisRange( -0.5, 12.5, "X" );
			return numeratorND;

		default:
			AliError ( "unimplemented type" );
			return 0;
	}
}

Double_t AliMultiplicityEfficiencySelection::GetGlobalTriggeringEficiency ( Int_t type ) {
	Double_t n = 0, d = 0;

	switch ( type ) {
		case 0:
			n = ( Double_t ) nSelected[0] + ( Double_t ) nSelected[1] + ( Double_t ) nSelected[2];
			d = ( Double_t ) nGenerated[0] + ( Double_t ) nGenerated[1] + ( Double_t ) nGenerated[2];
			return n / d;
			break;

		case 3:
			n = ( Double_t ) nSelected[0] + ( Double_t ) nSelected[2];
			d = ( Double_t ) nGenerated[0] + ( Double_t ) nGenerated[2];
			return n / d;

		default:
			AliWarningF ( "Processtype %d not eligible!", type );
			return 0;
			break;
	}

	return 0;
}

Bool_t AliMultiplicityEfficiencySelection::AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t, Bool_t ) {
	Int_t label = track->GetLabel();

	if ( label < 0 ) return kFALSE;

	TParticle* mother = AliMultiplicityHelper::GetMother ( stack, label );

	if ( mother ) {
		if ( TMath::Abs ( mother->GetPdgCode() ) == 3122 || mother->GetPdgCode() == 310 ) {

			if ( TMath::Abs ( mother->GetPdgCode() ) == 3122 ) nSurvived[ ( tType - 1 ) * 2 + 0] += 1;

			if ( TMath::Abs ( mother->GetPdgCode() ) == 310 ) nSurvived[ ( tType - 1 ) * 2 + 1] += 1;

			if ( bIsGenerated ) hSurvivedEta[tType - 1]->Fill ( stack->Particle ( track->GetLabel() )->Eta() );

			return kTRUE;
		}
	}

	return kFALSE;
}

Bool_t AliMultiplicityEfficiencySelection::AcceptTrackletMC ( const AliMultiplicity* m, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t ) {
	Int_t label0 = m->GetLabel ( iTracklet, 0 );
	Int_t label1 = m->GetLabel ( iTracklet, 1 );

	if ( label0 < 0 && label1 < 0 ) return kFALSE;

	Int_t label = label0; if ( label < 0 ) label = label1;

	TParticle* mother0 = AliMultiplicityHelper::GetMother ( stack, label0 );
	TParticle* mother1 = AliMultiplicityHelper::GetMother ( stack, label1 );

	Bool_t accept = kFALSE;

	if ( mother0 || mother1 ) {
		if ( mother0 && ( ( TMath::Abs ( mother0->GetPdgCode() ) == 3122 ) || ( mother0->GetPdgCode() == 310 ) ) ) {
			if ( TMath::Abs ( mother0->GetPdgCode() ) == 3122 ) nSurvived[ ( tType - 1 ) * 2 + 0] += 1;

			if ( TMath::Abs ( mother0->GetPdgCode() ) == 310 ) nSurvived[ ( tType - 1 ) * 2 + 1] += 1;

			if ( bIsGenerated ) hSurvivedEta[tType - 1]->Fill ( stack->Particle ( label0 )->Eta() );

			accept = kTRUE;
		}

		if ( mother1 && ( label0 != label1 ) && ( ( TMath::Abs ( mother1->GetPdgCode() ) == 3122 ) || ( mother1->GetPdgCode() == 310 ) ) ) {
			if ( TMath::Abs ( mother1->GetPdgCode() ) == 3122 ) nSurvived[ ( tType - 1 ) * 2 + 0] += 1;

			if ( TMath::Abs ( mother1->GetPdgCode() ) == 310 ) nSurvived[ ( tType - 1 ) * 2 + 1] += 1;

			if ( bIsGenerated ) hSurvivedEta[tType - 1]->Fill ( stack->Particle ( label1 )->Eta() );

			accept = kTRUE;
		}

		return accept;
	}

	return kFALSE;
}

void AliMultiplicityEfficiencySelection::ResetCounters() {
	AliMultiplicitySelection::ResetCounters();

	nGeneratedS[0] = 0;
	nGeneratedS[1] = 0;

	for ( Int_t i = 0; i < 3; ++i ) {
		nSelected[i] = 0;
		nGenerated[i] = 0;

		for ( Int_t j = 0; j < 2; ++j ) {
			nSurvived[i * 2 + j] = 0;
		}
	}

	currentProcessType = -1;
}

void AliMultiplicityEfficiencySelection::SetWeightDiffraction ( TH2D* hwSDs, TH2D* hwNDDDs ) {
	if ( !hwSDs ) {
		AliError ( "Single-diffractive weights unavailable." );
	}

	if ( !hwNDDDs ) {
		AliError ( "Non- and double-diffractive weights unavailable." );
	}

	if ( !hwSDs || !hwNDDDs ) return;

	hwSD = new TH2D ( *hwSDs );
	hwNDDD = new TH2D ( *hwNDDDs );

	bWeightDiffraction = kTRUE;
}

TH1D* AliMultiplicityEfficiencySelection::GetAlphaSD() {
	TH1D* halpha = dynamic_cast<TH1D*> ( hEventsGenerated[1]->Clone() );
	halpha->SetName ( "alphaSD" );

	TH1D* denominatorINEL = dynamic_cast<TH1D*> ( hEventsGenerated[0]->Clone() );
	denominatorINEL->Add ( hEventsGenerated[1] );
	denominatorINEL->Add ( hEventsGenerated[2] );

	halpha->Divide ( denominatorINEL );
	delete denominatorINEL;
	return halpha;
}

TH1D* AliMultiplicityEfficiencySelection::GetAlphaDD() {
	TH1D* halpha = dynamic_cast<TH1D*> ( hEventsGenerated[2]->Clone() );
	halpha->SetName ( "alphaDD" );
	
	TH1D* denominatorINEL = dynamic_cast<TH1D*> ( hEventsGenerated[0]->Clone() );
	denominatorINEL->Add ( hEventsGenerated[1] );
	denominatorINEL->Add ( hEventsGenerated[2] );
	
	halpha->Divide ( denominatorINEL );
	delete denominatorINEL;
	return halpha;
}

TCanvas* AliMultiplicityEfficiencySelection::DrawC() {
	TCanvas* c = new TCanvas();
// 	hMZZ[0]->SetLineWidth(2);
// 	hMZZ[0]->SetLineColor(kRed);
// 	hMZZ[0]->DrawNormalized();
// 	hMZZ[1]->SetLineWidth(2);
// 	hMZZ[1]->DrawNormalized("same");
	gPad->SetLogy();
	gPad->SetGrid ( 1, 1 );
	return c;
}
