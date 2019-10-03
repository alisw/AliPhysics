#include "TLegend.h"
#include <TAttMarker.h>
#include "AliMultiplicityHelper.h"
#include "AliMultiplicityPileup.h"
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TCollection.h>
#include <TMathBase.h>
#include <TLine.h>
#include <TAttLine.h>

ClassImp ( AliMultiplicityPileup )

AliMultiplicityPileup::AliMultiplicityPileup ( const char* name, const char* title ) : AliMultiplicitySelection ( name, title ),
	hGenerated ( 0 ),
	hNHeaders ( 0 ),
	nParticlesPrimary ( 0 ),
	nParticlesSecondary ( 0 ),
	nHeaders ( 0 ),
	lHeaders ( 0 ),
	nCoincidence ( -1 ),
	nClosest ( -1 ),
	ZESD ( 0 ),
	dZESD ( 0 ),
	dZ ( 0 )

{
	for ( Int_t i = 0; i < 3; ++i ) {
		hdZvsExtra[i] = 0;
		hdZvsExtraNC[i] = 0;
		hResponse[i] = 0;
		hResponse_true[i] = 0;

// 			hSPDPileupEfficiency[i] = 0;

		for ( Int_t j = 0; j < 2; ++j ) {
			hSelected[3 * j + i] = 0;
		}
	}

// 	for (Int_t i = 0; i < 4; ++i) {
// 		hSampleBias[i] = 0;
// 	}

	for ( Int_t i = 0; i < 2; ++i ) {
		hVertexCorrelation[i] = 0;
	}

	SetNeedsMC ( kTRUE );
	SetCorrelate ( kTRUE );
	SetCollectEvent ( kFALSE );
	SetCorrelateTracks ( kTRUE );
	SetQueryPhysSelection ( kTRUE );
	ResetCounters();
}

AliMultiplicityPileup::~AliMultiplicityPileup() {

}

Bool_t AliMultiplicityPileup::AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask ) {
	fCurrentESDEventAccepted = AliMultiplicitySelection::AcceptESDEvent ( esd, mask ) && !excludeESD;
	AliMultiplicitySelection::AcceptMCEvent ( mc, mask );

	if ( fFollowMCselection ) {
		if ( fPreSelected && fCurrentMCEventAccepted ) fCurrentESDEventAccepted = kTRUE;
	}

	if ( fCurrentESDEventAccepted ) {
		nCoincidence = -1;
		nClosest = -1;

		ZESD = esd->GetPrimaryVertex()->GetZ();
		dZESD = esd->GetPrimaryVertex()->GetZRes();
		AliGenCocktailEventHeader* cheader = dynamic_cast<AliGenCocktailEventHeader*> ( AliMultiplicityHelper::GetGenEventHeader ( mc ) );
		lHeaders = cheader->GetHeaders();
		nHeaders = lHeaders->GetSize();
		hNHeaders->Fill ( nHeaders );

		//loop through headers and find the closest vertex
		Double_t deltaZ = 100; //too much
		TIter headIt ( lHeaders );
		AliGenEventHeader* header = 0;
		Int_t count = -1;
		TArrayF adZ ( lHeaders->GetSize() );
		TArrayF aZ ( lHeaders->GetSize() );

		do { // Go trhough all headers and look for the correct one
			header = dynamic_cast<AliGenEventHeader*> ( headIt() );

			if ( header ) {
				count++;
				TArrayF v;
				header->PrimaryVertex ( v );
				Double_t deltaZn = TMath::Abs ( ZESD - v[2] );

				adZ.AddAt ( deltaZn, count );
				aZ.AddAt ( v[2], count );

				if ( deltaZn < deltaZ ) {
					nCoincidence = count;
					deltaZ = deltaZn;
				}
			}
		} while ( header );

		Double_t deltaZk = 100;

		for ( Int_t i = 0; i < adZ.GetSize(); ++i ) {
			if ( ( adZ[i] < deltaZk ) && ( adZ[i] > deltaZ ) ) {
				deltaZk = adZ[i];
				nClosest = i;
			}
		}

		dZ = aZ[nClosest] - aZ[nCoincidence];
	}

	return fCurrentESDEventAccepted && fCurrentMCEventAccepted;
}

Bool_t AliMultiplicityPileup::AcceptESDTrackMC ( const AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement ) {
	Bool_t trackAccepted = kFALSE;
// 	AliMultiplicitySelection::AcceptESDTrack(track,tType,complement);

	Int_t label = TMath::Abs ( track->GetLabel() );
	Bool_t isPrimary = stack->IsPhysicalPrimary ( label );
	Int_t mlabel = 0;
	Bool_t isPrimaryMother = isPrimary;

	if ( !isPrimary ) {
		mlabel = AliMultiplicityHelper::FindPrimaryMother ( stack, label );
		isPrimaryMother = stack->IsPhysicalPrimary ( mlabel );
	}

	if ( isPrimary || isPrimaryMother ) {
		Int_t V = 0;
		Int_t nEvent = FindHeader ( isPrimary ? label : mlabel );

		if ( nEvent != nCoincidence ) {
			if ( nEvent > 0 ) {
				V = 1;
			} else {
				V = 2;
			}
		}

		if ( TMath::Abs ( track->Eta() ) < GetEtaCut() ) { trackAccepted = kTRUE; }

		if ( trackAccepted ) {
			switch ( tType ) {
				case kITSSA:
					if ( !complement ) {
						++nSelectedITSSA[V];
					}

					break;

				case kITSTPC:
					if ( complement ) {
						++nSelectedGlobalITSSAcomplements[V];
					} else {
						++nSelectedGlobal[V];
					}

					break;

				default:
					AliError ( "Wrong estimator type!" );
					break;
			}
		}
	}

	return trackAccepted;
}

Bool_t AliMultiplicityPileup::AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType ) {
	Bool_t trackAccepted = kFALSE;
// 	AliMultiplicitySelection::AcceptTracklet(m,iTracklet,tType);

	if ( UseTrackletChi2Cut() ) {
		Double_t chi2 = AliMultiplicityHelper::GetTrackletChi2 ( spdm, iTracklet );

		if ( chi2 > GetTrackletChi2Cut() ) {
			AliWarningF ( "Tracklet chi2 too high: %.3f", chi2 );
			return trackAccepted;
		}
	}

	Int_t label[2];
	Bool_t isPrimary[2];
	Int_t mlabel[2];
	Bool_t isPrimaryMother[2];

	for ( Int_t i = 0; i < 2; ++i ) {
		label[i] = TMath::Abs ( spdm->GetLabel ( iTracklet, i ) );
		isPrimary[i] = stack->IsPhysicalPrimary ( label[i] );

		if ( !isPrimary[i] ) {
			mlabel[i] = AliMultiplicityHelper::FindPrimaryMother ( stack, label[i] );
			isPrimaryMother[i] = stack->IsPhysicalPrimary ( mlabel[i] );
		}
	}

	Bool_t isTPrimary = isPrimary[0] || isPrimary[1];
	Bool_t isTPrimaryMother = isPrimaryMother[0] || isPrimaryMother[1];

	if ( isTPrimary || isTPrimaryMother ) {
		Int_t V = 0;
		Int_t id = 0;

		if ( isPrimary[0] ) { id = label[0]; }
		else if ( isPrimaryMother[0] ) { id = mlabel[0]; }
		else if ( isPrimary[1] ) { id = label[1]; }
		else if ( isPrimaryMother[1] ) { id = mlabel[1]; }

		Int_t nEvent = FindHeader ( id );

		if ( nEvent != nCoincidence ) { V = 1; }

		if ( TMath::Abs ( spdm->GetEta ( iTracklet ) ) < GetEtaCut() ) { trackAccepted = kTRUE; }

		if ( trackAccepted ) {
			switch ( tType ) {
				case kSPD:
					++nSelectedTracklets[V];
					trackAccepted = kTRUE;
					break;

				case kITSSA:
					++nSelectedITSSATrackletComplements[V];
					trackAccepted = kTRUE;
					break;

				case kITSTPC:
					++nSelectedGlobalTrackletcomplements[V];
					trackAccepted = kTRUE;
					break;

				default:
					AliError ( "Wrong estimator type!" );
					break;
			}
		}
	}

	return trackAccepted;
}

Bool_t AliMultiplicityPileup::AcceptParticle ( const TParticle* particle, const Int_t iParticle ) {
	TParticlePDG* pdg = particle->GetPDG();

	if ( TMath::Abs ( pdg->Charge() ) > 0 ) {
		if ( TMath::Abs ( particle->Eta() ) < GetEtaCut() ) {
			if ( FindHeader ( iParticle ) == nCoincidence ) {
				nParticlesPrimary++;
			}

			if ( FindHeader ( iParticle ) == nClosest ) {
				nParticlesSecondary++;
			}
		}
	}

	return AliMultiplicitySelection::AcceptParticle ( particle, iParticle );
}


void AliMultiplicityPileup::CreateSelectionHistograms() {
	AliMultiplicitySelection::CreateSelectionHistograms();
	hGenerated = new TH1D ( "hGenerated", "Generated multiplicity of coinciding event", 161, -0.5, 160.5 );
// 	lHeaders = new TList();

	for ( Int_t i = 0; i < 3; ++i ) {
		hdZvsExtra[i] = new TH2D ( Form ( "hdZvsExtra_%d", i ), "", 601, -30.1, 30.1, 161, -0.5, 160.5 );
		hdZvsExtraNC[i] = new TH2D ( Form ( "hdZvsExtraNC_%d", i ), "", 601, -30.1, 30.1, 161, -0.5, 160.5 );
		hResponse[i] = new TH2D ( Form ( "hResponseA_%d", i ), "", 161, -0.5, 160.5, 161, -0.5, 160.5 );
		hResponse_true[i] = new TH2D ( Form ( "hResponseT_%d", i ), "", 161, -0.5, 160.5, 161, -0.5, 160.5 );

// 			hSPDPileupEfficiency[i] = new TH1D ( Form("hSPDPileupEfficiency_%d",i), "", 161, -0.5, 160.5 );

		for ( Int_t j = 0; j < 2; ++j ) {
			hSelected[3 * j + i] = new TH1D ( Form ( "hSelected_%d_%s", i, j == 0 ? "prima" : "extra" ), "", 161, -0.5, 160.5 );
		}
	}

// 	for (Int_t i = 0; i < 4; ++i) {
// 		hSampleBias[i] = new TH1D ( Form("hSampleBias_%d",i), "", 161, -0.5, 160.5 );
// 	}

	for ( Int_t i = 0; i < 2; ++i ) {
		hVertexCorrelation[i] = new TH2D ( Form ( "hVertexCorrelation_%d", i ), "", 601, -30.1, 30.1, 601, -30.1, 30.1 );
	}

// 	hFakeRate = new TH1D ( "hFakeRate", "", 161, -0.5, 160.5 );

	hNHeaders = new TH1D ( "hNHeaders", "Number of MC headers", 10, -0.5, 9.5 );
}

void AliMultiplicityPileup::SaveSelectionHistograms() {
	AliMultiplicitySelection::SaveSelectionHistograms();

	if ( hGenerated ) { hGenerated->Write(); }

	if ( hNHeaders ) { hNHeaders->Write(); }

	for ( Int_t i = 0; i < 3; ++i ) {
		if ( hdZvsExtra[i] ) { hdZvsExtra[i]->Write(); }

		if ( hdZvsExtraNC[i] ) { hdZvsExtraNC[i]->Write(); }

		if ( hResponse[i] ) { hResponse[i]->Write(); }

		if ( hResponse_true[i] ) { hResponse_true[i]->Write(); }

		for ( Int_t j = 0; j < 2; ++j ) {
			if ( hSelected[3 * j + i] ) { hSelected[3 * j + i]->Write(); }
		}
	}
}

void AliMultiplicityPileup::ResetCounters() {
	AliMultiplicitySelection::ResetCounters();
	nParticlesPrimary = 0;
	nParticlesSecondary = 0;

	for ( Int_t i = 0; i < 3; ++i ) {
		nSelectedGlobal[i] = 0;
		nSelectedITSSA[i] = 0;
		nSelectedTracklets[i] = 0;
		nSelectedGlobalITSSAcomplements[i] = 0;
		nSelectedGlobalTrackletcomplements[i] = 0;
		nSelectedITSSATrackletComplements[i] = 0;
	}
}

void AliMultiplicityPileup::Conclude() {
	AliMultiplicitySelection::Conclude();
	hGenerated->Fill ( nParticlesPrimary );
// 	hSampleBias[0]->Fill( nParticlesPrimary );
// 	hSampleBias[1]->Fill( nParticlesSecondary );

	if ( fCurrentESDEventAccepted ) {
		hdZvsExtra[0]->Fill ( dZ, nSelectedGlobal[1] + nSelectedGlobalITSSAcomplements[1] + nSelectedGlobalTrackletcomplements[1] );
		hdZvsExtra[1]->Fill ( dZ, nSelectedITSSA[1] + nSelectedITSSATrackletComplements[1] );
		hdZvsExtra[2]->Fill ( dZ, nSelectedTracklets[1] );

		hdZvsExtraNC[0]->Fill ( dZ, nSelectedGlobal[2] + nSelectedGlobalITSSAcomplements[2] + nSelectedGlobalTrackletcomplements[2] );
		hdZvsExtraNC[1]->Fill ( dZ, nSelectedITSSA[2] + nSelectedITSSATrackletComplements[2] );
		hdZvsExtraNC[2]->Fill ( dZ, nSelectedTracklets[2] );

		for ( Int_t i = 0; i < 2; ++i ) {
// 					if ( fCurrentESDEventAccepted ) {
			hSelected[3 * i + 0]->Fill ( nSelectedGlobal[i] + nSelectedGlobalITSSAcomplements[i] + nSelectedGlobalTrackletcomplements[i] );
			hSelected[3 * i + 1]->Fill ( nSelectedITSSA[i] + nSelectedITSSATrackletComplements[i] );
			hSelected[3 * i + 2]->Fill ( nSelectedTracklets[i] );
// 						}
		}
	}

	if ( fCurrentESDEventAccepted ) {
		hResponse[0]->Fill ( nParticlesPrimary, nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements );
		hResponse[1]->Fill ( nParticlesPrimary, nITSTracks + nITSTrackletComplements );
		hResponse_true[0]->Fill ( nParticlesPrimary, nSelectedGlobal[0] + nSelectedGlobalITSSAcomplements[0] + nSelectedGlobalTrackletcomplements[0] );
		hResponse_true[1]->Fill ( nParticlesPrimary, nSelectedITSSA[0] + nSelectedITSSATrackletComplements[0] );
		hResponse[2]->Fill ( nParticlesPrimary, nTracklets );
		hResponse_true[2]->Fill ( nParticlesPrimary, nSelectedTracklets[0] );
	}
}

Long64_t AliMultiplicityPileup::Merge ( const TCollection* mergelist ) {
	AliMultiplicitySelection::Merge ( mergelist );

	if ( !mergelist ) {
		AliDebug ( 1, "NULL merge list." );
		return 0;
	}

	if ( mergelist->IsEmpty() ) {
		AliDebug ( 1, "Empty merge list." );
		return 1;
	}

	TList hg;
	TList hs[3 * 2];
	TList hr[3];
	TList hrt[3];
	TList hdz[3];
	TList hdzc[3];
	TList hh;

// 	TH1D* hSampleBias[4];	//!
// 	TH2D* hVertexCorrelation[2];		//!
//
// 	TH1D* hSPDPileupEfficiency[3];		//!
//
// 	TH2D* hFakeRate;					//!

// 	TList hsb[4];
	TList hvc[2];
// 	TList hple[3];
// 	TList hf;

	Int_t count = 0;

	TIter nextE ( mergelist );

	while ( AliMultiplicityPileup* selection = dynamic_cast<AliMultiplicityPileup*> ( nextE() ) ) {
		if ( !selection ) { continue; }

		hg.Add ( selection->hGenerated );
// 			hf.Add( selection->hFakeRate );

		for ( Int_t i = 0; i < 3; ++i ) {
// 				hple[i].Add( selection->hSPDPileupEfficiency[i] );
			for ( Int_t j = 0; j < 2; ++j ) {
				hs[3 * j + i].Add ( selection->hSelected[3 * j + i] );
			}

			hdz[i].Add ( selection->hdZvsExtra[i] );
			hdzc[i].Add ( selection->hdZvsExtraNC[i] );
			hr[i].Add ( selection->hResponse[i] );
			hrt[i].Add ( selection->hResponse_true[i] );
		}

// 			for (Int_t i = 0; i < 4; ++i) {
// 				hsb[i].Add( selection->hSampleBias[i] );
// 			}
		for ( Int_t i = 0; i < 2; ++i ) {
			hvc[i].Add ( selection->hVertexCorrelation[i] );
		}

		hh.Add ( selection->hNHeaders );
		count++;
	}

	hGenerated->Merge ( &hg );
	hNHeaders->Merge ( &hh );
// 	hFakeRate->Merge( &hf );

	for ( Int_t i = 0; i < 3; ++i ) {
// 		hSPDPileupEfficiency[i]->Merge( &hple[i] );
		for ( Int_t j = 0; j < 2; ++j ) {
			hSelected[3 * j + i]->Merge ( &hs[3 * j + i] );
		}

		hdZvsExtra[i]->Merge ( &hdz[i] );
		hdZvsExtraNC[i]->Merge ( &hdzc[i] );
		hResponse[i]->Merge ( &hr[i] );
		hResponse_true[i]->Merge ( &hrt[i] );
	}

// 	for (Int_t i = 0; i < 4; ++i) {
// 		hSampleBias[i]->Merge( &hsb[i] );
// 	}
	for ( Int_t i = 0; i < 2; ++i ) {
		hVertexCorrelation[i]->Merge ( &hvc[i] );
	}

	return ( Long64_t ) count + 1;
}

void AliMultiplicityPileup::Result() {
	AliMultiplicitySelection::Result();
	TCanvas* c = new TCanvas ( "c", "c", 900, 600 );
	c->Divide ( 3, 2 );

	for ( Int_t i = 0; i < 3; ++i ) {
		c->cd ( i + 1 );
		hdZvsExtra[i]->SetAxisRange ( -0.5, 40.5, "Y" );
		hdZvsExtra[i]->Draw ( "colz" );
		gPad->SetLogz();

		c->cd ( 3 + i + 1 );
		hdZvsExtraNC[i]->SetAxisRange ( -0.5, 40.5, "Y" );
		hdZvsExtraNC[i]->Draw ( "colz" );
		gPad->SetLogz();
	}

	new TCanvas ( "d", "d", 800, 600 );
	hNHeaders->Draw();

	TCanvas* e = new TCanvas ( "e", "e", 900, 600 );
	e->Divide ( 3, 2 );

	for ( Int_t i = 0; i < 3; ++i ) {
		e->cd ( i + 1 );
		hResponse[i]->SetAxisRange ( -0.5, 70.5, "X" );
		hResponse[i]->SetAxisRange ( -0.5, 70.5, "Y" );
		hResponse[i]->Draw ( "colz" );
		gPad->SetGrid ( 1, 1 );
		gPad->SetLogz();
		TLine* l = new TLine ( -0.5, -0.5, hResponse[i]->GetBinLowEdge ( hResponse[i]->GetNbinsX() ) + 1, hResponse[i]->GetBinLowEdge ( hResponse[i]->GetNbinsY() ) + 1 );
		l->SetLineStyle ( kDashed );
		l->SetLineWidth ( 2 );
		l->Draw();
	}

	for ( Int_t i = 0; i < 3; ++i ) {
		e->cd ( 3 + i + 1 );
		hResponse_true[i]->SetAxisRange ( -0.5, 70.5, "X" );
		hResponse_true[i]->SetAxisRange ( -0.5, 70.5, "Y" );
		hResponse_true[i]->Draw ( "colz" );
		gPad->SetGrid ( 1, 1 );
		gPad->SetLogz();
		TLine* l = new TLine ( -0.5, -0.5, hResponse_true[i]->GetBinLowEdge ( hResponse_true[i]->GetNbinsX() ) + 1, hResponse_true[i]->GetBinLowEdge ( hResponse_true[i]->GetNbinsY() ) + 1 );
		l->SetLineStyle ( kDashed );
		l->SetLineWidth ( 2 );
		l->Draw();
	}
}

Int_t AliMultiplicityPileup::FindHeader ( const Int_t iParticle ) {
	TIter headIt ( lHeaders );
	Int_t nproduced = 0;
	AliGenEventHeader* header = 0;
	Int_t count = -1;

	do { // Go trhough all headers and look for the correct one
		header = dynamic_cast<AliGenEventHeader*> ( headIt() );

		if ( header ) {
			nproduced += header->NProduced();
			count++;
		}
	} while ( header && iParticle >= nproduced );

	return count;
}

TH1D* AliMultiplicityPileup::GetGenerated1() {
	return hGenerated;
}

TH2D* AliMultiplicityPileup::GetExtraTvsZ ( Int_t iMethod ) {
	if ( ( iMethod >= 0 ) && ( iMethod < 3 ) ) return hdZvsExtra[iMethod];

	AliWarningF ( "Out of bounds: %d (need -1 < i < 3)", iMethod );
	return 0;
}

TH2D* AliMultiplicityPileup::GetExtraTvsZwo0 ( Int_t iMethod ) {
	if ( ( iMethod < 0 ) && ( iMethod > 2 ) ) {
		AliWarningF ( "Out of bounds: %d (need -1 < i < 3)", iMethod );
		return 0;
	}

	TH2D* sh = hdZvsExtra[iMethod];
	TH2D* hdZvsExtrawo0 = new TH2D ( Form ( "%s_wo0", sh->GetName() ), "", sh->GetNbinsX(), sh->GetXaxis()->GetBinLowEdge ( 1 ), sh->GetXaxis()->GetBinLowEdge ( sh->GetNbinsX() ) + 1, 160, 0.5, 160.5 );

	for ( Int_t i = 0; i < sh->GetNbinsX(); ++i ) {
		for ( Int_t j = 1; j < sh->GetNbinsY(); ++j ) {
			hdZvsExtrawo0->SetBinContent ( i + 1, j, sh->GetBinContent ( i + 1, j + 1 ) );
		}
	}

	return hdZvsExtrawo0;
}

TH2D* AliMultiplicityPileup::GetResponse ( Int_t iMethod, Bool_t bTrue ) {
	if ( ( iMethod < 0 ) && ( iMethod > 2 ) ) {
		AliWarningF ( "Out of bounds: %d (need -1 < i < 3)", iMethod );
		return 0;
	}

	if ( bTrue ) return hResponse_true[iMethod];

	return hResponse[iMethod];
}

Double_t AliMultiplicityPileup::GetAffectedFraction ( Int_t iMethod ) {
	if ( ( iMethod < 0 ) && ( iMethod > 2 ) ) {
		AliWarningF ( "Out of bounds: %d (need -1 < i < 3)", iMethod );
		return -1.0;
	}

	return hdZvsExtra[iMethod]->Integral ( 1, hdZvsExtra[iMethod]->GetNbinsX(), 2, hdZvsExtra[iMethod]->GetNbinsY() ) / hdZvsExtra[iMethod]->GetEntries();
}
