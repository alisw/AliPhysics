#include "AliMultiplicityResponseSelection.h"
#include <TF1.h>
#include <TLine.h>
#include <TMathBase.h>
#include <TRandom.h>
#include <TPad.h>

ClassImp ( AliMultiplicityResponseSelection )

void AliMultiplicityResponseSelection::Conclude(  ) {
	if ( fCurrentESDEventAccepted ) AliMultiplicitySelection::Conclude();
	if ( fNoVertexESD ) {
		hMinus1binRedist->Fill ( nMCparticles );
	} else {
		if ( fCurrentESDEventAccepted && fCurrentMCEventAccepted ) {
			response[kITSTPC - 1]->Fill ( nMCparticles, nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements );
			response[kITSSA - 1]->Fill ( nMCparticles, nITSTracks + nITSTrackletComplements );
			response[kSPD - 1]->Fill ( nMCparticles, nTracklets );
		}
	}
}

void AliMultiplicityResponseSelection::SaveSelectionHistograms() {
	AliMultiplicitySelection::SaveSelectionHistograms();

	for ( Int_t i = 0; i < kSPD; i++ ) {
		if ( response[i] ) response[i]->Write();
	}

// 	if ( hMultiplicityMC ) hMultiplicityMC->Write();
}

void AliMultiplicityResponseSelection::Result() {
	const Double_t sliceat = 40.0;

	TH2D* RITSTPC = NormalizeResponse ( response[kITSTPC - 1] );
	TH2D* RITSSA = NormalizeResponse ( response[kITSSA - 1] );
	TH2D* RSPD = NormalizeResponse ( response[kSPD - 1] );

	TCanvas* cnv  = GetCanvas();
	cnv->cd();

	TPad* pTl = new TPad ( "tl", "tl", 0, 0.3, 0.33, 1 );
	TPad* pTm = new TPad ( "tm", "tm", 0.33, 0.3, 0.66, 1 );
	TPad* pTr = new TPad ( "tr", "tr", 0.66, 0.3, 1, 1 );
	TPad* pBl = new TPad ( "bl", "bl", 0, 0, 0.33, 0.3 );
	TPad* pBm = new TPad ( "bm", "bm", 0.33, 0, 0.66, 0.3 );
	TPad* pBr = new TPad ( "br", "br", 0.66, 0, 1, 0.3 );
	pTl->Draw();
	pTm->Draw();
	pTr->Draw();
	pBl->Draw();
	pBm->Draw();
	pBr->Draw();


	pTl->cd()->SetLogz();
	RSPD->SetTitle ( "Response for SPD tracklets method" );
	RSPD->GetXaxis()->SetTitle ( TString::Format ( "generated charged multiplicity in |#eta|<%.1f, N_{gen}", GetEtaCut() ) );
	RSPD->GetYaxis()->SetTitle ( "measured charged multiplicity, N_{m}" );

	RSPD->Draw ( "COLZ" );
	TLine* slice_lineT = new TLine ( sliceat, 0, sliceat, RSPD->GetNbinsY() );
	slice_lineT->SetLineWidth ( 3 );
	slice_lineT->Draw();
	TLine* diagT = new TLine ( 0, 0, RSPD->GetNbinsX(), RSPD->GetNbinsY() );
	diagT->SetLineWidth ( 3 );
	diagT->SetLineStyle ( 2 );
	diagT->Draw();

	pTm->cd()->SetLogz();
	RITSSA->SetTitle ( "Response for ITSSA+ method" );
	RITSSA->GetXaxis()->SetTitle ( TString::Format ( "generated charged multiplicity in |#eta|<%.1f, N_{gen}", GetEtaCut() ) );
	RITSSA->GetYaxis()->SetTitle ( "measured charged multiplicity, N_{m}" );

	RITSSA->Draw ( "COLZ" );
	TLine* slice_lineI = new TLine ( sliceat, 0, sliceat, RITSSA->GetNbinsY() );
	slice_lineI->SetLineWidth ( 3 );
	slice_lineI->Draw();
	TLine* diagI = new TLine ( 0, 0, RITSSA->GetNbinsX(), RITSSA->GetNbinsY() );
	diagI->SetLineWidth ( 3 );
	diagI->SetLineStyle ( 2 );
	diagI->Draw();

	pTr->cd()->SetLogz();
	RITSTPC->SetTitle ( "Response for ITSTPC+ method" );
	RITSTPC->GetXaxis()->SetTitle ( TString::Format ( "generated charged multiplicity in |#eta|<%.1f, N_{gen}", GetEtaCut() ) );
	RITSTPC->GetYaxis()->SetTitle ( "measured charged multiplicity, N_{m}" );

	RITSTPC->Draw ( "COLZ" );
	TLine* slice_lineG = new TLine ( sliceat, 0, sliceat, RITSTPC->GetNbinsY() );
	slice_lineG->SetLineWidth ( 3 );
	slice_lineG->Draw();
	TLine* diagG = new TLine ( 0, 0, RITSTPC->GetNbinsX(), RITSTPC->GetNbinsY() );
	diagG->SetLineWidth ( 3 );
	diagG->SetLineStyle ( 2 );
	diagG->Draw();

	pBl->cd()->SetLogy();
	TH1D* RSPDp = RSPD->ProjectionY ( "T", sliceat, sliceat );
	RSPDp->Fit ( "gaus", "ql" );
	RSPDp->SetStats ( kFALSE );
	RSPDp->SetTitle ( "" );
	TF1* fitT = RSPDp->GetFunction ( "gaus" );
	fitT->SetLineColor ( kRed );
	RSPDp->Draw ( "P" );

	pBm->cd()->SetLogy();
	TH1D* RITSSAp = RITSSA->ProjectionY ( "I", sliceat, sliceat );
	RITSSAp->Fit ( "gaus", "ql" );
	RITSSAp->SetStats ( kFALSE );
	RITSSAp->SetTitle ( "" );
	TF1* fitI = RITSSAp->GetFunction ( "gaus" );
	fitI->SetLineColor ( kRed );
	RITSSAp->Draw ( "P" );

	pBr->cd()->SetLogy();
	TH1D* RITSTPCp = RITSTPC->ProjectionY ( "G", sliceat, sliceat );
	RITSTPCp->Fit ( "gaus", "ql" );
	RITSTPCp->SetStats ( kFALSE );
	RITSTPCp->SetTitle ( "" );
	TF1* fitG = RITSTPCp->GetFunction ( "gaus" );
	fitG->SetLineColor ( kRed );
	RITSTPCp->Draw ( "P" );
}

void AliMultiplicityResponseSelection::CreateSelectionHistograms() {
	AliMultiplicitySelection::CreateSelectionHistograms();
	TString names[3] = {"ITSTPC", "ITSSA", "SPD"};

	for ( Int_t i = 0; i < kSPD; i++ ) {
		response[i] = new TH2D ( TString::Format ( "response_%s", names[i].Data() ), TString::Format ( "Response %s", names[i].Data() ), 161, -0.5, 160.5, 161, -0.5, 160.5 );
	}

// 	hMultiplicityMC = new TH1D ( "hMultiplicityMC", "MC generator-level multiplicity", 161, -0.5, 160.5 );
	hMinus1binRedist = new TH1D ( "hMinus1binRedist", "-1 bin redistribution histogram", 161, -0.5, 160.5 ); //numbers of events with no vertex per generated multiplicity bin
}

Long64_t AliMultiplicityResponseSelection::Merge ( const TCollection* mergelist ) {
	AliMultiplicitySelection::Merge ( mergelist );

	if ( !mergelist ) {
		AliError ( "NULL merge list." );
		return 0;
	}

	if ( mergelist->IsEmpty() ) {
		AliWarning ( "Empty merge list." );
		return 1;
	}

	TList histograms[4];
// 	TList histm;
	Int_t count = 0;
	
// 	TList lresponse_fail;

	TIter nextE ( mergelist );

	while ( AliMultiplicityResponseSelection* selection = dynamic_cast<AliMultiplicityResponseSelection*> ( nextE() ) ) {
		if ( !selection ) continue;

// 		histm.Add ( selection->hMultiplicityMC );

		for ( Int_t i = 0; i < kSPD; i++ ) {
			histograms[i].Add ( selection->response[i] );
		}

		histograms[3].Add ( selection->hMinus1binRedist );

// 		if (fFollowESDselection) lresponse_fail.Add( selection->response_fail );
		count++;
	}

// 	hMultiplicityMC->Merge ( &histm );

	for ( Int_t i = 0; i < kSPD; i++ ) {
		response[i]->Merge ( &histograms[i] );
	}

	hMinus1binRedist->Merge ( &histograms[3] );
	
// 	if (fFollowESDselection) response_fail->Merge( &lresponse_fail );

	return ( Long64_t ) count + 1;
}

Bool_t AliMultiplicityResponseSelection::AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask ) {
	SetCollectEvent ( kFALSE );
	fNoVertexESD = excludeESD;

	fCurrentESDEventAccepted = AliMultiplicitySelection::AcceptESDEvent ( esd, mask ) && !excludeESD;
	AliMultiplicitySelection::AcceptMCEvent ( mc, mask );
	
	if ( fFollowMCselection ) {
		if ( fPreSelected && fCurrentMCEventAccepted ) fCurrentESDEventAccepted = kTRUE;
	}

	return fCurrentESDEventAccepted && fCurrentMCEventAccepted;
}

AliMultiplicityResponseSelection::AliMultiplicityResponseSelection ( const char* name, const char* title ) : AliMultiplicitySelection ( name, title ),
// 	hMultiplicityMC ( 0 ),
	hMinus1binRedist ( 0 ),
	hPtWeights ( 0 ),
	fRWStrangeness ( kFALSE )  {
	for ( Int_t i = 0; i < kSPD; i++ ) {
		response[i] = 0;
	}
	
	SetNeedsMC ( kTRUE );
	SetCorrelate ( kTRUE );
	SetCollectEvent ( kFALSE );
	SetCorrelateTracks ( kTRUE );
	SetFollowESDselection( kTRUE );
}

void AliMultiplicityResponseSelection::SetPtWeights ( TH1D* hweights, Bool_t RWSec ) {
	if ( hweights ) {
		hPtWeights = hweights; 
		fRWStrangeness = RWSec;
	} else {
		AliWarning ( "Weights histogram undefined" );
	}
}


AliMultiplicityResponseSelection::~AliMultiplicityResponseSelection() {

}

TH2D* AliMultiplicityResponseSelection::NormalizeResponse ( TH2D* response, TH1D* efficiency )  {
	response->Sumw2();
	TH2D* out = dynamic_cast<TH2D*> ( response->Clone() );
	out->SetName ( TString::Format ( efficiency ? "Normalized_eff_%s" : "Normalized_%s" , response->GetName() ) );

	for ( Int_t i = 1; i <= response->GetNbinsX(); i++ ) {
		Double_t slice = response->Integral ( i, i, 1, response->GetNbinsY() );

		if ( slice < 1e-10 ) slice = 1;

		if ( efficiency && ( efficiency->GetBinContent ( i ) > 0 ) ) slice /= efficiency->GetBinContent ( i );

		for ( Int_t j = 1; j <= response->GetNbinsY(); j++ ) {
			Double_t val    = response->GetBinContent ( i, j ) / slice;
			Double_t vale   = response->GetBinError ( i, j ) / slice;
			out->SetBinContent ( i, j, TMath::Finite ( val ) ? val : 0 );
			out->SetBinError ( i, j, TMath::Finite ( val ) ? vale : 0 );
		}
	}

	return out;
}

TH2D* AliMultiplicityResponseSelection::GetResponse ( AliMultiplicityEventSelection::mEstimators estimator ) {
	return response[estimator - 1];
}

TH2D* AliMultiplicityResponseSelection::GetResponseNormalized ( AliMultiplicityEventSelection::mEstimators estimator ) {
	return AliMultiplicityResponseSelection::NormalizeResponse ( GetResponse ( estimator ) );
}

TH1D* AliMultiplicityResponseSelection::GetSeedDistribution() {
	return hMultiplicityMC;
}

Int_t AliMultiplicityResponseSelection::GetTrackWeight ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators , Bool_t ) {
	Int_t factor = 1;
	Int_t label = TMath::Abs ( track->GetLabel() );

	if ( hPtWeights )  {
		Double_t check = 0;
		Bool_t isPrimary = stack->IsPhysicalPrimary ( label );

		Int_t mlabel = AliMultiplicityHelper::FindPrimaryMother ( stack, label );

		TParticle* mother = AliMultiplicityHelper::GetMother ( stack, mlabel );

// 		if (  ) {
		// it cannot be just multiplied because we cannot count "half of a particle"
		// instead a random generator decides if the particle is counted twice (if value > 1)
		// or not (if value < 0)
		Int_t nbins = hPtWeights->GetNbinsX();
		Int_t ibin = hPtWeights->FindBin ( stack->Particle ( label )->Pt() );
		Double_t content = hPtWeights->GetBinContent ( ibin );
		if ( fRWStrangeness ) {
			if ( !isPrimary && mother && ( ( TMath::Abs ( mother->GetPdgCode() ) == 3122 ) || ( mother->GetPdgCode() == 310 ) ) ) {
				if ( stack->Particle ( label )->Pt() <= 1.5 ) {
					check = content;
				} else {
					check = hPtWeights->GetBinContent ( nbins );
				}
			}
		} else {
			check = content;
		}

		if ( check > 0 ) {
			Float_t random = gRandom->Uniform();

			if ( check > 1 && random < check - 1 ) {
				factor = 2;
			} else if ( check < 1 && random < 1 - check ) {
				factor = 0;
			}
		}
// 		}
	}

	return factor;
}

Int_t AliMultiplicityResponseSelection::GetTrackWeight ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators ) {
	Int_t factor = 1;
	Int_t label0 = spdmultiplicity->GetLabel ( iTracklet, 0 );
	Int_t label1 = spdmultiplicity->GetLabel ( iTracklet, 1 );

	if ( label0 != label1 ) return factor;

	if ( hPtWeights )  {
		Double_t check = 0;
		Bool_t isPrimary = stack->IsPhysicalPrimary ( TMath::Abs ( label0 ) );
		Int_t mlabel = AliMultiplicityHelper::FindPrimaryMother ( stack, label0 );

		TParticle* mother = AliMultiplicityHelper::GetMother ( stack, mlabel );

		// it cannot be just multiplied because we cannot count "half of a particle"
		// instead a random generator decides if the particle is counted twice (if value > 1)
		// or not (if value < 0)
		Int_t nbins = hPtWeights->GetNbinsX();
		Int_t ibin = hPtWeights->FindBin ( stack->Particle ( label0 )->Pt() );
		Double_t content = hPtWeights->GetBinContent ( ibin );
		if ( fRWStrangeness ) {
			if ( !isPrimary && mother && ( ( TMath::Abs ( mother->GetPdgCode() ) == 3122 ) || ( mother->GetPdgCode() == 310 ) ) ) {
				if ( stack->Particle ( label0 )->Pt() <= 1.5 ) {
					check = content;
				} else {
					check = hPtWeights->GetBinContent ( nbins );
				}
			}
		} else {
			check = content;
		}
		
		if ( check > 0 ) {
			Float_t random = gRandom->Uniform();
			
			if ( check > 1 && random < check - 1 ) {
				factor = 2;
			} else if ( check < 1 && random < 1 - check ) {
				factor = 0;
			}
		}
	}

	return factor;
}


Bool_t AliMultiplicityResponseSelection::AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement, Bool_t count ) {
	Int_t factor = GetTrackWeight ( track, stack, tType, complement );

	if ( factor == 0 ) {
// 		AliInfo ( "Not counting track" );
		return kFALSE;
	}

	if ( factor == 2 ) {
		AliMultiplicitySelection::AcceptESDTrack ( track, tType, complement, !count );
// 		AliInfo ( "Counting double" );
	}

	Bool_t accepted = AliMultiplicitySelection::AcceptESDTrack ( track, tType, complement );

	return accepted;
}

Bool_t AliMultiplicityResponseSelection::AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t ) {
	Int_t factor = GetTrackWeight ( spdm, iTracklet, stack, tType );

	if ( factor == 0 ) {
// 		AliInfo ( "Not counting track" );
		return kFALSE;
	}

	if ( factor == 2 ) {
		AliMultiplicitySelection::AcceptTracklet ( spdm, iTracklet, tType );
// 		AliInfo ( "Counting double" );
	}

	Bool_t accepted = AliMultiplicitySelection::AcceptTracklet ( spdm, iTracklet, tType );
	
	return accepted;
}

Bool_t AliMultiplicityResponseSelection::AcceptParticle ( TParticle* particle, const Int_t, Bool_t ) {
	Int_t factor = GetParticleWeight( particle, -1 );
	
	if ( factor == 0 ) {
// 		AliInfo ( "Not counting particle" );
		return kFALSE;
	}
		
	if ( factor == 2 ) {
		AliMultiplicitySelection::AcceptParticle ( particle, -1, kFALSE );
// 		AliInfo ( "Counting particle double" );
	}
	
	Bool_t accepted = AliMultiplicitySelection::AcceptParticle ( particle, -1 );

	return accepted;
}

Int_t AliMultiplicityResponseSelection::GetParticleWeight ( TParticle* particle, const Int_t ) {
	Int_t factor = 1;
	if ( hPtWeights && !fRWStrangeness)  {
		Int_t ibin = hPtWeights->FindBin ( particle->Pt() );
		Double_t check = 0;
		check = hPtWeights->GetBinContent ( ibin );
		
		if ( check > 0 ) {
			Float_t random = gRandom->Uniform();
			
			if ( check > 1 && random < check - 1 ) {
				factor = 2;
			} else if ( check < 1 && random < 1 - check ) {
				factor = 0;
			}
		}
	}
	return factor;
}

Bool_t AliMultiplicityResponseSelection::CollectCondition() {
	Bool_t coll = !fCurrentMCEventAccepted;
	
// 	if ( TMath::Abs( fZcacheGen ) > GetZvCut() + 0.5 )  {
// 		coll = kTRUE;
// 	}
	
	return coll;
}
