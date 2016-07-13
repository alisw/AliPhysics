#include "AliMultiplicitySelection.h"
#include "AliMultiplicityHelper.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include <TCanvas.h>
#include <TMathBase.h>
#include <TList.h>
#include <TAttMarker.h>
#include <iostream>

ClassImp ( AliMultiplicitySelection )

Bool_t AliMultiplicitySelection::AcceptTracklet ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliMultiplicityEventSelection::mEstimators tType ) {
	Bool_t trackAccepted = kFALSE;

	if ( UseTrackletChi2Cut() ) {
		Double_t chi2 = AliMultiplicityHelper::GetTrackletChi2 ( spdmultiplicity, iTracklet );

		if ( chi2 > GetTrackletChi2Cut() ) {
			AliWarningF ( "Tracklet chi2 too high: %.3f", chi2 );
			return trackAccepted;
		}
	}

	if ( TMath::Abs ( spdmultiplicity->GetEta ( iTracklet ) ) < GetEtaCut() ) trackAccepted = kTRUE;
	hEta[tType-1]->Fill(spdmultiplicity->GetEta(iTracklet));
	if ( trackAccepted ) {
		switch ( tType ) {
			case kSPD:
				nTracklets++;
// 				trackAccepted = kTRUE;
				hPhi[2]->Fill(spdmultiplicity->GetPhi(iTracklet) * 180.0 / TMath::Pi());
// 				hTrackletChi2->Fill ( AliMultiplicityHelper::GetTrackletChi2 ( spdmultiplicity, iTracklet ) );
// 				hTrackletChi2[1]->Fill ( AliMultiplicityHelper::GetTrackletChi2 ( spdmultiplicity, iTracklet, kTRUE ) );
				break;

			case kITSSA:
				nITSTrackletComplements++;
// 				trackAccepted = kTRUE;
				hPhi[1]->Fill(spdmultiplicity->GetPhi(iTracklet) * 180.0 / TMath::Pi());
				break;

			case kITSTPC:
				nGlobalTrackletComplements++;
// 				trackAccepted = kTRUE;
				hPhi[0]->Fill(spdmultiplicity->GetPhi(iTracklet) * 180.0 / TMath::Pi());
				break;

			default:
				AliError ( "Wrong estimator type!" );
				break;
		}
	}
	return trackAccepted;
}

Bool_t AliMultiplicitySelection::AcceptESDTrack ( AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement, Bool_t count ) {
	Bool_t trackAccepted = kFALSE;

	if ( TMath::Abs ( track->Eta() ) < GetEtaCut() ) trackAccepted = kTRUE;
	hEta[tType-1]->Fill(track->Eta());
	if ( trackAccepted ) {
// 		AliESDtrack* t = dynamic_cast<AliESDtrack*> ( track->Clone() );
		

		switch ( tType ) {
			case kITSSA:
				if ( !complement ) {
					nITSTracks++;
					hPhi[1]->Fill(track->Phi() * 180.0 / TMath::Pi());
				}
				break;

			case kITSTPC:
				hPhi[0]->Fill(track->Phi() * 180.0 / TMath::Pi());
				if ( complement ) {
					nGlobalITScomplements++;
// 					if ( count && wITSComplements ) tracks->Add ( track );
				} else {
					nGlobalTracks++;
// 					if ( count ) tracks->Add ( track );
				}
				break;
			default:
				AliError ( "Wrong estimator type!" );
				break;
		}
	}

	return trackAccepted;
}

Bool_t AliMultiplicitySelection::AcceptESDEvent ( const AliESDEvent* esd, UInt_t mask ) {
	currentMask = mask;
	
	Bool_t vertexInCuts = IsVertexInCuts ( esd, mask );
		
	UInt_t processes = ~( AliMultiplicityHelper::kNonDiffractive | AliMultiplicityHelper::kSingleDiffractive | AliMultiplicityHelper::kDoubleDiffractive );
	UInt_t sample = sampleDef & processes;
	UInt_t check = currentMask & processes;
	
	fPreSelected = (Bool_t) ( sample == (check & sample) );
	
	isPileup = ( Bool_t ) ( ( mask & AliMultiplicityHelper::kSPDPileUp ) != 0 );
		
	if ( fSPDInCuts ) {
		hZSPD->Fill(AliMultiplicityHelper::GetVertex(AliMultiplicityHelper::kVertexSPD, esd)->GetZ(), 0);
		if ( (Bool_t) ( ( mask & AliMultiplicityHelper::kSPDPileUp ) != 0 ) ) {
			hZSPD->Fill(AliMultiplicityHelper::GetVertex(AliMultiplicityHelper::kVertexSPD, esd)->GetZ(), 1);
		}
	}
	
	fCurrentESDEventAccepted = fPreSelected && vertexInCuts;

// 	if ( isPileup ) {
// 		const AliVertex* vSPDp = AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexSPD, esd );
// 		const AliESDVertex* const vSPDplus = esd->GetPileupVertexSPD ( 0 );
// 		const Double_t distance = vSPDp->GetZ() - vSPDplus->GetZ();
// 		Int_t bin = hPileUpDist[0]->FindBin ( distance );
// 		hPileUpDist[0]->Fill ( distance, 1.0 / hPileUpDist[0]->GetBinWidth ( bin ) );
// 
// 		if ( fCurrentESDEventAccepted ) hPileUpDist[1]->Fill ( distance, 1.0 / hPileUpDist[1]->GetBinWidth ( bin ) );
// 	}

	return fCurrentESDEventAccepted;
}

// void AliMultiplicitySelection::Add0bin ( const AliESDEvent* ) {
// 	if ( GetQueryPhysSelection() ? AliMultiplicityHelper::IsSelected ( GetSelectionTrigger() ) : kTRUE ) {
// 		nNoVertexEvents++;
// 	}
// }
// 
// Long64_t AliMultiplicitySelection::GetNoVertexEvents() {
// 	return nNoVertexEvents;
// }

Long64_t AliMultiplicitySelection::Merge ( const TCollection* mergelist ) {
	if ( !mergelist ) {
		AliDebug ( 1, "NULL merge list." );
		return 0;
	}

	if ( mergelist->IsEmpty() ) {
		AliDebug ( 1, "Empty merge list." );
		return 1;
	}

	MergeStats ( mergelist );

	TList histograms[kSPD];
// 	TList hpile[kSPD];
// 	TList htpc;
	TList hmc;
// 	TList hTchi2[2];
// 	TList hpt;
	TList hptnch;
	TList hphil[kSPD];
// 	TList hptnchgen;
// 	TList hptnchgens[kSPD];
// 	TList hpd[2];
// 	TList hetal[kSPD];
// 	TList hTTl;
	TList hzs;
// 	TList lzrg[3];
// 	TList ldzrg[3];
	Int_t count = 0;

	TIter nextE ( mergelist );

	while ( AliMultiplicitySelection* selection = dynamic_cast<AliMultiplicitySelection*> ( nextE() ) ) {
		if ( !selection ) continue;

		if ( ( GetEtaCut() == selection->GetEtaCut() ) && ( GetZvCut() == selection->GetZvCut() ) ) {
			for ( Int_t i = 0; i < kSPD; i++ ) {
				histograms[i].Add ( selection->hMultiplicity[i] );
				hphil[i].Add( selection->hPhi[i] );
// 				hpile[i].Add ( selection->hPileupvsMult[i] );
// 				hetal[i].Add ( selection->hEta[i] );
			}

// 			htpc.Add ( selection->hMultiplicityTracks );
// 			hTTl.Add ( selection->hTT );
			hzs.Add ( selection->hZSPD );

			if ( GetNeedsMC() ) {
				hmc.Add ( selection->hMultiplicityMC );
// 				hptnchgen.Add ( selection->fhPtvsNchGen );

// 				for ( Int_t i = 0; i < kSPD; ++i ) {
// // 					hptnchgens[i].Add ( selection->fhPtvsNchGenS[i] );
// 				}
				
// 				for (Int_t i = 0; i < 3; ++i) {
// // 					lzrg[i].Add( selection->hZrg[i] );
// 					ldzrg[i].Add( selection->hDZg[i] );
// 				}
			}

// 			hTchi2[0].Add ( selection->hTrackletChi2[0] );
// 			hTchi2[1].Add ( selection->hTrackletChi2[1] );

			hptnch.Add ( selection->fhPtvsNch );
// 			hpd[0].Add ( selection->hPileUpDist[0] );
// 			hpd[1].Add ( selection->hPileUpDist[1] );
			count++;
		} else {
			AliWarning ( "Not merging selections with different parameters." );
		}
	}

	for ( Int_t i = 0; i < kSPD; i++ ) {
		hMultiplicity[i]->Merge ( &histograms[i] );
		hPhi[i]->Merge( &hphil[i] );
// 		hPileupvsMult[i]->Merge ( &hpile[i] );
// 		hEta[i]->Merge ( &hetal[i] );
	}

// 	hTT->Merge ( &hTTl );
	
// 	hMultiplicityTracks->Merge ( &htpc );

	if ( GetNeedsMC() ) {
		hMultiplicityMC->Merge ( &hmc );
// 		fhPtvsNchGen->Merge ( &hptnchgen );

// 		for ( Int_t i = 0; i < kSPD; ++i ) {
// 			fhPtvsNchGenS[i]->Merge ( &hptnchgens[i] );
// 		}
		
// 		for (Int_t i = 0; i < 3; ++i) {
// 			hZrg[i]->Merge(&lzrg[i]);
// 			hDZg[i]->Merge(&ldzrg[i]);
// 		}
	}

// 	hTrackletChi2[0]->Merge ( &hTchi2[0] );
// 	hTrackletChi2[1]->Merge ( &hTchi2[1] );

// 	hTrackPt->Merge ( &hpt );

	fhPtvsNch->Merge ( &hptnch );

// 	hPileUpDist[0]->Merge ( &hpd[0] );
// 	hPileUpDist[1]->Merge ( &hpd[1] );
	
	hZSPD->Merge ( &hzs );

	return ( Long64_t ) count + 1;
}

void AliMultiplicitySelection::CreateSelectionHistograms() {
	TString names[kSPD] = {"ITSTPC", "ITSSA", "SPD"};
	const Int_t nEtaBins = 21;
	Double_t etabins[nEtaBins+1];
	for (Int_t i = 0; i <= nEtaBins; ++i) {
		etabins[i] = -2.2 + 0.2 * i; 
	}

	for ( Int_t i = 0; i < kSPD; i++ ) {
		hMultiplicity[i] = new TH1D ( names[i] + " raw", names[i] + " Raw multiplicity", 161, -0.5, 160.5 );
		hPileupvsMult[i] = new TH1D ( names[i] + " pileup", names[i] + " pileup count", 161, -0.5, 160.5 );
		hEta[i] = new TH1D ( names[i] + "_eta", names[i] + " #eta", nEtaBins, etabins );
		hPhi[i] = new TH1D ( Form("hPhi_%d", i), "", 360, -0.5, 359.5 );
	}

	hTT = new TH2D ( "hTT", "", 161, -0.5, 160.5, 161, -0.5, 160.5 );
	
	hMultiplicityTracks = new TH1D ( "hTPConlyMult", "TPC only multiplicity", 161, -0.5, 160.5 );

	if ( GetNeedsMC() ) {
		hMultiplicityMC = new TH1D ( "hMultiplicityMC", "Generated multiplicity", 161, -0.5, 160.5 );
		fhPtvsNchGen = new TH2D ( "pPtvsNchGen", "<p_t> vs. N_{CH}", 161, -0.5, 160.5, 3001, -0.01, 30.01 );

		for ( Int_t i = 0; i < kSPD; ++i ) {
			fhPtvsNchGenS[i] = new TH2D ( Form ( "pPtvsNchGenS_%s", names[i].Data() ), "", 161, -0.5, 160.5, 3001, -0.01, 30.01 );
// 			particlesS[i] = new TList();
// 			particlesS[i]->SetOwner();
		}

// 		particles = new TList();
// 		particles->SetOwner();
	}

// 	hTrackletChi2 = new TH1D ( "hTrackletChi2", "Tracklet Chi2 distribution (without sin() scaling)", 501, -0.01, 5.01 );
// 	hTrackletChi2[1] = new TH1D ( "hTrackletChi2Sin", "Tracklet Chi2 distribution (with sin() scaling)", 501, -0.01, 5.01 );

// 	hTrackPt = new TH2D ( "hTrackPt", "p_t distribution from tracks",  501, -0.1, 50.1, 3, 0.5, 3.5 );
// 	hTrackPt->GetYaxis()->SetBinLabel ( 1, "global" );
// 	hTrackPt->GetYaxis()->SetBinLabel ( 2, "ITS complementary" );
// 	hTrackPt->GetYaxis()->SetBinLabel ( 3, "ITS" );

	fhPtvsNch = new TH2D ( "pPtvsNch", "<p_t> vs. N_{CH}", 161, -0.5, 160.5, 301, -0.1, 30.1 ); // only the tracks have pt info, but the Nch is for the full estimator
// 	tracks = new TList();
// 	tracks->SetOwner();
	
// 	secondaryTracks = new TList();
// 	secondaryTracks->SetOwner();
	
	hZSPD = new TH2D ( "hZ", "", 201, -10.1, 10.1, 2, -0.5, 1.5 );

// 	hPileUpDist[0] = new TH1D ( "hPileUpDist_0", "Distance between primary and secondary SPD vertices", 100, -20.0, 20.0 /*nbins, xbins*/ );
// 	hPileUpDist[1] = new TH1D ( "hPileUpDist_1", "Distance between primary and secondary SPD vertices", 100, -20.0, 20.0/*nbins, xbins*/ );
	
// 	if ( GetNeedsMC() ) {
// 		for (Int_t i = 0; i < 3; ++i) {
// 			hEtarg[i] = new TH1D ( Form("hEtarg_%d",i), "", 25, -2.5, 2.5 );
// 			hZetarg[i] = new TH2D ( Form("hZetarg_%d",i), "", 301, -30.2, 30.2, 25, -2.5, 2.5 );
// 			hZrg[i] = new TH2D ( Form("hZrg_%d", i), "", 2001, -20.02, 20.02, 2001, -20.02, 20.02 );
// 			hDZg[i] = new TH2D ( Form("hDZg_%d", i), "", 2001, -20.02, 20.02, 2001, -20.02, 20.02 );
// 		}
// 	}
// 	hZrg[2] = new TH2D ( "hZrg_2", "", 1001, -10.02, 10.02, 1001, -10.02, 10.02 );
}

AliMultiplicitySelection::AliMultiplicitySelection ( const char* name, const char* title ) : AliMultiplicityAnalysisSelection ( name, title ),
	nMCparticles ( 0 ),
	hMultiplicityTracks ( 0 ),
	wITSComplements(kFALSE),
	hMultiplicityMC ( 0 ),
// 	hTrackletChi2 ( 0 ),
// 	nNoVertexEvents ( 0 ),
	useTrackletChi2Cut ( kTRUE ),
// 	hTrackPt ( 0 ),
	fhPtvsNch ( 0 ),
	tracks ( 0 ),
// 	fPileThreshold ( 0.8 ),
	fhPtvsNchGen ( 0 ),
	particles ( 0 ),
	nGlobalSecondaryTracks ( 0 ),
	nGlobalSecondaryITScomplements ( 0 ),
	nITSSecondaryTracks ( 0 ),
	secondaryTracks ( 0 ),
	hTT ( 0 ),
	hZSPD ( 0 ) {
	for ( Int_t i = 0; i < kSPD; i++ ) {
		hMultiplicity[i] = 0;
		hPileupvsMult[i] = 0;
		fhPtvsNchGenS[i] = 0;
		particlesS[i] = 0;
		hEta[i] = 0;
		hPhi[i] = 0;
	}
	
// 	for (Int_t i = 0; i < 2; ++i) {
// // 		hEtarg[i] = 0;
// // 		hZetarg[i] = 0;
// 		hZrg[i] = 0;
// 		hDZg[i] = 0;
// // 		fZcache[i] = 0;
// 	}
// 	fZcache[2] = 0;
// 	hZrg[2] = 0;
// 	hDZg[2] = 0;

// 	hTrackletChi2[0] = 0;
// 	hTrackletChi2[1] = 0;
	hPileUpDist[0] = 0;
	hPileUpDist[1] = 0;
	SetNeedsMC ( kFALSE );
// 	SetCheckPileUp ( kFALSE );
	ResetCounters();
	SetTrackletChi2Cut ( 1.6 );
// 	SetForceConclude ( kFALSE );
}

// void AliMultiplicitySelection::SetPileupThreshold ( Double_t threshold ) {
// 	if ( ( threshold < 0 ) || threshold > 2.0 ) {
// 		AliWarningF ( "Treshold value incorrect: %.2f", threshold );
// 		return;
// 	}
// 
// // 	fPileThreshold = threshold;
// }

void AliMultiplicitySelection::Conclude() {
	Bool_t countIt = kTRUE;
	if ( ( Bool_t ) ( ( sampleDef & AliMultiplicityHelper::kSPDPileUp ) != 0 ) ) {
		if ( isPileup ) countIt = kFALSE;
	}; 
	
/*	if ( CollectCondition() ) {
		SetCollectEvent(kTRUE);
	}*/	
	
	if ( countIt ) {
		if ( fCurrentESDEventAccepted ) {
			hMultiplicity[0]->Fill ( nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements );
			hMultiplicity[1]->Fill ( nITSTracks + nITSTrackletComplements );
			hMultiplicityTracks->Fill ( wITSComplements ? nGlobalTracks + nGlobalITScomplements : nGlobalTracks );
			hMultiplicity[2]->Fill ( nTracklets );

// 			Int_t nTrk = tracks->GetEntries();

// 			for ( Int_t i = 0; i < nTrk; ++i ) {
// 				AliESDtrack* t = dynamic_cast<AliESDtrack*> ( tracks->At ( i ) );
// 				fhPtvsNch->Fill ( wITSComplements ? nGlobalTracks + nGlobalITScomplements : nGlobalTracks, t->Pt() );
// 			}
// 			Int_t nTracks = 0;
// 			if ( secondaryTracks && !secondaryTracks->IsEmpty() ) {
// 				nTracks = wITSComplements ? nGlobalSecondaryTracks : nGlobalSecondaryTracks + nGlobalSecondaryITScomplements;
// 			} else {
// 				nTracks = wITSComplements ? nGlobalTracks + nGlobalITScomplements : nGlobalTracks;
// 			}
// 			hTT->Fill(nTracklets, nTracks);
		}
	}

	if ( isPileup ) {
		if ( fCurrentESDEventAccepted ) {
			hPileupvsMult[0]->Fill ( nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements );
			hPileupvsMult[1]->Fill ( nITSTracks + nITSTrackletComplements );
			hPileupvsMult[2]->Fill ( nTracklets );
		}
	}
	if ( GetNeedsMC() ) {
		hMultiplicityMC->Fill ( nMCparticles );
// 		for ( Int_t i = 0; i < particles->GetEntries(); i++ ) {
// 			TParticle* p = dynamic_cast<TParticle*> ( particles->At ( i ) );
// 			if ( p ) fhPtvsNchGen->Fill ( nMCparticles, p->Pt() );
// 		}
		
/*		if (fPreSelected) {			
			Double_t dz[3];
			for (Int_t i = 0; i < 3; ++i) {
				dz[i] = fZcache[i] - fZcacheGen;
			}
			
// 			for (Int_t i = 0; i < 3; ++i) {
// 				hZrg[i]->Fill( fZcache[i], fZcacheGen );
// 				hDZg[i]->Fill( fZcacheGen, dz[i] );
// 			}
		}*/	
	
// 		if ( fCurrentMCEventAccepted ) {
// 			if ( fCurrentESDEventAccepted ) {
// 				for (Int_t i = 0; i < kSPD; ++i) {
// 					for (Int_t j = 0; j < particlesS[i]->GetEntries(); ++j ) {
// 						TParticle* p = dynamic_cast<TParticle*> ( particlesS[i]->At(j) );
// 						if ( p ) fhPtvsNchGenS[i]->Fill( nMCparticles, p->Pt() );
// 					}
// 				}
// 			}
// 		}
	}
}

void AliMultiplicitySelection::ResetCounters() {
	nGlobalTracks = 0;
	nGlobalITScomplements = 0;
	nGlobalTrackletComplements = 0;
	nITSTrackletComplements = 0;
	nITSTracks = 0;
	nTracklets = 0;
	nMCparticles = 0;
	
	nGlobalSecondaryTracks = 0;
	nGlobalSecondaryITScomplements = 0;
	nITSSecondaryTracks = 0;
	
// 	if ( secondaryTracks && !secondaryTracks->IsEmpty() ) {
// 		delete secondaryTracks;                                   //->Clear();
// 		secondaryTracks = new TList();
// 	}
// 
// 	if ( tracks && !tracks->IsEmpty() ) {
// 		delete tracks;                                            //->Clear();
// 		tracks = new TList();
// 	}
// 
// 	if ( particles && !particles->IsEmpty() ) {
// 		delete particles;                                         //->Clear();
// 		particles = new TList();
// 	}
// 	
// 	for (Int_t i = 0; i < kSPD; ++i) {
// 		if ( particlesS[i] && !(particlesS[i]->IsEmpty()) ) {
// 			delete particlesS[i];                                    //->Clear();
// 			particlesS[i] = new TList();
// 		}
// 	}
}


void AliMultiplicitySelection::Result() {
	TCanvas* cnv = GetCanvas();
	cnv->cd()->SetLogy ( 1 );

	for ( Int_t i = 0; i < kSPD; i++ ) {
		hMultiplicity[i]->SetMarkerStyle ( 21 + i );
		hMultiplicity[i]->SetMarkerColor ( kBlue + i * 2 );
		hMultiplicity[i]->Draw ( i == 0 ? "P" : "P SAME" );
	}

	if ( GetNeedsMC() ) {
		hMultiplicityMC->SetLineColor ( kRed );
		hMultiplicityMC->Draw ( "hist same" );
	}
}

void AliMultiplicitySelection::SaveSelectionHistograms() {
	for ( Int_t i = 0; i < kSPD; i++ ) {
		if ( hMultiplicity[i] ) hMultiplicity[i]->Write();

		if ( hPileupvsMult[i] ) hPileupvsMult[i]->Write();
		
		if ( hEta[i] ) {
			hEta[i]->Scale(1.0/hStats->GetBinContent(1),"width");
			hEta[i]->Write();
		}
	}

	if ( hMultiplicityTracks ) hMultiplicityTracks->Write();

	if ( GetNeedsMC() && hMultiplicityMC ) {
		hMultiplicityMC->Write();
	}

// 	if ( hTrackletChi2 ) hTrackletChi2->Write();

// 	if ( hTrackletChi2[1] ) hTrackletChi2[1]->Write();

// 	if ( hTrackPt ) hTrackPt->Write();

	if ( hPileUpDist[0] ) hPileUpDist[0]->Write();

	if ( hPileUpDist[1] ) hPileUpDist[1]->Write();
	
	if ( hTT ) hTT->Write();
	
}

TH1D* AliMultiplicitySelection::GetMultiplicityHistogram ( AliMultiplicityEventSelection::mEstimators estimator ) {
	return hMultiplicity[estimator - 1];
}

Bool_t AliMultiplicitySelection::AcceptParticle ( TParticle* particle, const Int_t,  Bool_t count ) {
	TParticlePDG* pdg = particle->GetPDG();

	if ( TMath::Abs ( pdg->Charge() ) > 0 ) {
		if ( TMath::Abs ( particle->Eta() ) < GetEtaCut() ) {
// 			TParticle* p = dynamic_cast<TParticle*> ( particle->Clone() );
			nMCparticles++;
// 			if ( count ) particles->Add ( particle );
			return kTRUE;
		}
	}

	return kFALSE;
}

Bool_t AliMultiplicitySelection::AcceptMCEvent ( const AliMCEvent* mc, UInt_t ) {
	fCurrentMCEventAccepted = IsVertexInCuts( mc );
	
	if ( fFollowESDselection ) {
		fCurrentMCEventAccepted = fCurrentESDEventAccepted;
	}
	
	return fCurrentMCEventAccepted;
}

TH1D* AliMultiplicitySelection::GetTrackMultiplicity() {
	return hMultiplicityTracks;
}

Double_t AliMultiplicitySelection::GetTrackletChi2Cut() {
	return trackletChi2Cut;
}

void AliMultiplicitySelection::SetTrackletChi2Cut ( Double_t cut ) {
	trackletChi2Cut = cut;
	useTrackletChi2Cut = kTRUE;
}


Double_t AliMultiplicitySelection::GetPileupFraction ( AliMultiplicityEventSelection::mEstimators estimator ) {

	Double_t nTotal = hMultiplicity[estimator - 1]->GetEntries(); if ( nTotal < 1 ) { return 0; };

	Double_t nPileup = hPileupvsMult[estimator - 1]->GetEntries();

	Double_t fraction = nPileup / nTotal;

	return fraction;
}

TProfile* AliMultiplicitySelection::GetMeanPtvsMult() {
	if ( !fhPtvsNch ) return 0;

	TProfile* p = fhPtvsNch->ProfileX();
	return p;
}

TH1D* AliMultiplicitySelection::GetTrackPt ( Int_t Ngen ) {
	if ( ( Ngen > fhPtvsNchGen->GetNbinsX() ) ) return 0;
	if ( Ngen < 0 ) return fhPtvsNchGen->ProjectionY ( );

	return fhPtvsNchGen->ProjectionY ( "_py", fhPtvsNchGen->GetXaxis()->FindFixBin ( Ngen ), fhPtvsNchGen->GetXaxis()->FindFixBin ( Ngen ) );
}

Bool_t AliMultiplicitySelection::AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t,  Bool_t count ) {
	if ( TMath::Abs ( track->Eta() ) < GetEtaCut() ) {
		Int_t label = TMath::Abs ( track->GetLabel() );

		if ( stack->IsPhysicalPrimary ( label ) ) {
			TParticle* p = dynamic_cast<TParticle*>(stack->Particle ( label )->Clone());
			if ( count ) particlesS[tType - 1]->Add ( p );
			return kTRUE;
		}
	}

	return kFALSE;
}

Bool_t AliMultiplicitySelection::AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t count ) {
	if ( TMath::Abs ( spdm->GetEta ( iTracklet ) ) < GetEtaCut() ) {
		Int_t label[2];

		for ( Int_t i = 0; i < 2; ++i ) {
			label[i] = spdm->GetLabel ( iTracklet, i );
		}

		if ( label[0] != label[1] ) return kFALSE;

		if ( stack->IsPhysicalPrimary ( label[0] ) ) {
// 			TParticle* p = dynamic_cast<TParticle*>(stack->Particle ( label[0] )->Clone());
			if ( count ) particlesS[tType - 1]->Add ( stack->Particle ( label[0] ) );
			return kTRUE;
		}
	}

	return kFALSE;
}

Bool_t AliMultiplicitySelection::AcceptSecondaryESDTrack(AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement, Bool_t count )
{
	Bool_t trackAccepted = kFALSE;
	
	if ( TMath::Abs ( track->Eta() ) < GetEtaCut() ) trackAccepted = kTRUE;
	
	if ( trackAccepted ) {
// 		AliESDtrack* t = dynamic_cast<AliESDtrack*> ( track->Clone() );
		
		switch ( tType ) {
			case kITSSA:
				if ( !complement ) {
					nITSSecondaryTracks++;
				}
				break;
				
			case kITSTPC:
				if ( complement ) {
					nGlobalSecondaryITScomplements++;
					
// 					if ( count && wITSComplements ) secondaryTracks->Add ( track );
				} else {
					nGlobalSecondaryTracks++;
// 					if ( count ) secondaryTracks->Add ( track );
				}
				break;
			default:
				AliError ( "Wrong estimator type!" );
				break;
		}
	}
	
	return trackAccepted;
}

Bool_t AliMultiplicitySelection::CollectCondition ( ) {
	Bool_t coll = kFALSE;
// 	if ( (Bool_t) ( ( currentMask & AliMultiplicityHelper::kVNoGlobal ) != 0 ) && (nTracklets > 20) ) coll = kTRUE;
// 	if ( (Bool_t) ( ( currentMask & AliMultiplicityHelper::kVNoGlobal ) != 0 ) && (Bool_t) ( ( currentMask & AliMultiplicityHelper::kVNoSPD ) != 0 ) && (nMCparticles > 30) ) coll = kTRUE;
// 	if ( ( nGlobalTracks + nGlobalITScomplements + nGlobalTrackletComplements ) == 0 ) coll = kTRUE;
	
	return coll;
}
