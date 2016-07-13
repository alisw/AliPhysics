#include "AliMultiplicityEtaPhiSelection.h"
#include "AliMultiplicityHelper.h"
#include <TCanvas.h>
#include <TMathBase.h>

ClassImp ( AliMultiplicityEtaPhiSelection )

AliMultiplicityEtaPhiSelection::AliMultiplicityEtaPhiSelection ( const char* name, const char* title ) : AliMultiplicityAnalysisSelection ( name, title ),
	fVz ( 0 ),
	nClass ( 0 ) {
	for ( Int_t i = 0; i < kSPD; i++ ) {
			hEtaPhi[i] = 0;
			for ( Int_t j = 0; j < 3; ++j ) {
					hVzEta[i + 3 * j] = 0;
					hVzNch[i + 3 * j] = 0;
				}
			nTotal[i] = 0;
		}
	hEtaChi2[0] = 0;
	hEtaChi2[1] = 0;
	hBkg[0] = 0;
	hBkg[1] = 0;
	SetNeedsMC ( kFALSE );
	SetCollectEvent ( kFALSE );
}

AliMultiplicityEtaPhiSelection::~AliMultiplicityEtaPhiSelection() {

}

Bool_t AliMultiplicityEtaPhiSelection::AcceptESDEvent ( const AliESDEvent* esd ) {
	for ( Int_t i = 0; i < kSPD; ++i ) {
			nTotal[i] = 0;
			isClass[i] = kFALSE;
		}
	fVz = AliMultiplicityHelper::GetVertex ( kPrimaryVertex, esd )->GetZ();
	TString labels[4] = {"All", "V0AND", "MBOR", "MBOR > 0"};
	isClass[0] = AliMultiplicityHelper::IsSelected ( AliVEvent::kINT7 );
	isClass[1] = AliMultiplicityHelper::IsSelected ( AliVEvent::kMB );
	nClass->Fill ( labels[0], 1 );

	if ( isClass[1] ) {
			const AliMultiplicity* spdmult = esd->GetMultiplicity();
			/*if(!isClass[2])*/
			for ( Int_t iTracklet = 0; iTracklet < spdmult->GetNumberOfTracklets(); ++iTracklet ) {
					if ( TMath::Abs ( spdmult->GetEta ( iTracklet ) ) < 1.0 ) {
							isClass[2] = kTRUE;
							break;
						}
				}
		}


	return AliMultiplicityEventSelection::AcceptESDEvent ( esd );
}


Bool_t AliMultiplicityEtaPhiSelection::AcceptESDTrack ( const AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t ) {
	switch ( tType ) {
		case kSPD:
			break;
		case kITSSA:
			hEtaPhi[kITSSA - 1]->Fill ( track->Phi(), track->Eta() );
			for ( Int_t i = 0; i < 3; ++i ) {
					if ( isClass[i] ) hVzEta[kITSSA - 1 + 3 * i]->Fill ( fVz, track->Eta() );
				}
			nTotal[kITSSA - 1] += 1;
			break;
		case kITSTPC:
			hEtaPhi[kITSTPC - 1 ]->Fill ( track->Phi(), track->Eta() );
			for ( Int_t i = 0; i < 3; ++i ) {
					if ( isClass[i] ) hVzEta[kITSTPC - 1 + 3 * i]->Fill ( fVz, track->Eta() );
				}
			nTotal[kITSTPC - 1] += 1;
			break;
		default:
			AliError ( "Wrong estimator type!" );
			break;
		}

	return kTRUE;
}

Bool_t AliMultiplicityEtaPhiSelection::AcceptTracklet ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliMultiplicityEventSelection::mEstimators tType ) {
	switch ( tType ) {
		case kSPD:
			hEtaPhi[kSPD - 1]->Fill ( spdmultiplicity->GetPhi ( iTracklet ), spdmultiplicity->GetEta ( iTracklet ) );
			for ( Int_t i = 0; i < 3; ++i ) {
					if ( isClass[i] ) hVzEta[kSPD - 1 + 3 * i]->Fill ( fVz, spdmultiplicity->GetEta ( iTracklet ) );
				}
			nTotal[kSPD - 1] += 1;
			hEtaChi2[0]->Fill ( AliMultiplicityHelper::GetTrackletChi2 ( spdmultiplicity, iTracklet ), spdmultiplicity->GetEta ( iTracklet ) );
			if ( fPhiCut ) {
					if ( fPhiCut->CheckCut( spdmultiplicity->GetPhi(iTracklet) ) ) {
						hEtaChi2[1]->Fill ( AliMultiplicityHelper::GetTrackletChi2 ( spdmultiplicity, iTracklet ), spdmultiplicity->GetEta ( iTracklet ) );
					}
				}
			break;
		case kITSSA:
			hEtaPhi[kITSSA - 1]->Fill ( spdmultiplicity->GetPhi ( iTracklet ), spdmultiplicity->GetEta ( iTracklet ) );
			for ( Int_t i = 0; i < 3; ++i ) {
					if ( isClass[i] ) hVzEta[kITSSA - 1 + 3 * i]->Fill ( fVz, spdmultiplicity->GetEta ( iTracklet ) );
				}
			nTotal[kITSSA - 1] += 1;
			break;
		case kITSTPC:
			hEtaPhi[kITSTPC - 1]->Fill ( spdmultiplicity->GetPhi ( iTracklet ), spdmultiplicity->GetEta ( iTracklet ) );
			for ( Int_t i = 0; i < 3; ++i ) {
					if ( isClass[i] ) hVzEta[kITSTPC - 1 + 3 * i]->Fill ( fVz, spdmultiplicity->GetEta ( iTracklet ) );
				}
			nTotal[kITSTPC - 1] += 1;
			break;
		default:
			AliError ( "Wrong estimator type!" );
			break;
		}

	return kTRUE;
}

void AliMultiplicityEtaPhiSelection::Conclude() {
	TString labels[4] = {"All", "V0AND", "MBOR", "MBOR > 0"};
	for ( Int_t i = 0; i < 3; ++i ) { //class loop
			if ( isClass[i] ) nClass->Fill ( labels[i + 1], 1 );
			for ( Int_t j = 0; j < kSPD; ++j ) { //estimator loop
					if ( isClass[i] ) hVzNch[j + 3 * i]->Fill ( fVz, nTotal[j] );
				}
		}
}

void AliMultiplicityEtaPhiSelection::CreateSelectionHistograms() {
	TString names[kSPD] = {"ITSTPC", "ITSSA", "SPD"};
	TString names2[3] = {"V0AND", "MBOR", "MBORgt0"};
	for ( Int_t i = 0; i < kSPD; i++ ) {
			hEtaPhi[i] = new TH2D ( "", names[i].Data(), 400, 0, 2 * TMath::Pi(), 100, -2, 2 );
			hEtaPhi[i]->GetXaxis()->SetTitle ( "#phi" );
			hEtaPhi[i]->GetYaxis()->SetTitle ( "#eta" );
			for ( Int_t j = 0; j < 3; ++j ) {
					hVzEta[i + 3 * j] = new TH2D ( names2[j] + "_ZvsEta_" + names[i], names2[j] + " vertex Z vs #eta " + names[i], 403, -20.1, 20.1, 100, -2, 2 );
					hVzNch[i + 3 * j] = new TH2D ( names2[j] + "_ZvsNch" + names[i], names2[j] + " vertex Z vs N_{CH}" + names[i], 403, -20.1, 20.1, 161, -0.5, 160.5 );
				}
		}
	TString labels[4] = {"All", "V0AND", "MBOR", "MBOR > 0"};
	nClass = new TH1D ( "nClass", "Event counts", 4, -0.5, 3.5 );
	for ( Int_t i = 0; i < 4; ++i ) {
			nClass->GetXaxis()->SetBinLabel ( i + 1, labels[i] );
		}

	hEtaChi2[0] = new TH2D ( "hEtaChi2", "#chi^{2} vs. #eta", 100, 0, 10, 501, -2.5, 2.5 );
	hEtaChi2[0]->GetXaxis()->SetTitle ( "#chi^{2}" );
	hEtaChi2[0]->GetYaxis()->SetTitle ( "#eta" );
	hEtaChi2[0]->SetDrawOption ( "colz" );

	hEtaChi2[1] = new TH2D ( "hEtaChi2_phicut", "#chi^{2} vs. #eta (#phi cut)", 100, 0, 10, 501, -2.5, 2.5 );
	hEtaChi2[1]->GetXaxis()->SetTitle ( "#chi^{2}" );
	hEtaChi2[1]->GetYaxis()->SetTitle ( "#eta" );
	hEtaChi2[1]->SetDrawOption ( "colz" );
	
	hBkg[0] = new TH2D ( "hBkg", "Background tracklets", 100, 0, 10, 501, -2.5, 2.5 );
	hBkg[1] = new TH2D ( "hBkg_phicut", "Background tracklets", 100, 0, 10, 501, -2.5, 2.5 );
}

Long64_t AliMultiplicityEtaPhiSelection::Merge ( const TCollection* mergelist ) {
	if ( !mergelist ) {
			AliDebug ( 1, "NULL merge list." );
			return 0;
		}

	if ( mergelist->IsEmpty() ) {
			AliDebug ( 1, "Empty merge list." );
			return 1;
		}

	MergeStats ( mergelist );

	TList histogramsEP[kSPD];
	TList histogramsZE[9];
	TList histogramsZN[9];
	TList histogramsEC[2];
	TList nclasses;
	TList bkg[2];
	Int_t count = 0;

	TIter nextE ( mergelist );
	while ( AliMultiplicityEtaPhiSelection* selection = dynamic_cast<AliMultiplicityEtaPhiSelection*> ( nextE() ) ) {
			if ( !selection ) continue;
			for ( Int_t i = 0; i < kSPD; i++ ) {
					histogramsEP[i].Add ( selection->hEtaPhi[i] );
					for ( Int_t j = 0; j < 3; ++j ) {
							histogramsZE[i + 3 * j].Add ( selection->hVzEta[i + 3 * j] );
							histogramsZN[i + 3 * j].Add ( selection->hVzNch[i + 3 * j] );
						}
			}
			for (Int_t i = 0; i < 2; ++i) {
				histogramsEC[i].Add ( selection->hEtaChi2[i] );
			}
			
			nclasses.Add ( selection->nClass );
			bkg[0].Add ( selection->hBkg[0] );
			bkg[1].Add ( selection->hBkg[1] );
			count++;
		}

	for ( Int_t i = 0; i < kSPD; i++ ) {
			hEtaPhi[i]->Merge ( &histogramsEP[i] );
			for ( Int_t j = 0; j < 3; ++j ) {
					hVzEta[i + 3 * j]->Merge ( &histogramsZE[i + 3 * j] );
					hVzNch[i + 3 * j]->Merge ( &histogramsZN[i + 3 * j] );
				}
	}
	for (Int_t i = 0; i < 2; ++i) {
		hEtaChi2[i]->Merge ( &histogramsEC[i] );
	}	

	nClass->Merge ( &nclasses );
	hBkg[0]->Merge ( &bkg[0] );
	hBkg[1]->Merge ( &bkg[1] );
	return ( Long64_t ) count + 1;
}

void AliMultiplicityEtaPhiSelection::Result() {
	TCanvas* cnv = GetCanvas();
	cnv->Divide ( 1, 3 );
	for ( Int_t i = 0; i < kSPD; i++ ) {
			cnv->cd ( i + 1 );
			hEtaPhi[i]->Scale ( 1.0 / hStats->GetBinContent ( 1 ), "w" );
			hEtaPhi[i]->Draw ( "colz" );
		}

	cnv->SaveAs ( TString::Format ( "%s.png", GetName() ) );
}

void AliMultiplicityEtaPhiSelection::SaveSelectionHistograms() {
	for ( Int_t i = 0; i < kSPD; i++ ) {
			if ( hEtaPhi[i] ) hEtaPhi[i]->Write();
			for ( Int_t j = 0; j < 3; ++j ) {
					if ( hVzEta[i + 3 * j] ) hVzEta[i + 3 * j]->Write();
				}
			if ( hVzNch[i] ) hVzNch[i]->Write();
		}
}

TH2D* AliMultiplicityEtaPhiSelection::GetHistogram ( AliMultiplicityEventSelection::mEstimators type ) {
	return hEtaPhi[type - 1];
}

Bool_t AliMultiplicityEtaPhiSelection::AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent*, Bool_t ) {
	AcceptESDEvent( esd );
	fCurrentMCEventAccepted = fCurrentESDEventAccepted;
	return fCurrentMCEventAccepted;
}

Bool_t AliMultiplicityEtaPhiSelection::AcceptTrackletMC ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliStack*, AliMultiplicityEventSelection::mEstimators tType ) {
	if ( tType == kSPD ) {
		Double_t chi2 = AliMultiplicityHelper::GetTrackletChi2(spdmultiplicity,iTracklet,kTRUE);
		Int_t label0 = spdmultiplicity->GetLabel ( iTracklet, 0 );
		Int_t label1 = spdmultiplicity->GetLabel ( iTracklet, 1 );
// 		Bool_t isPrimary0 = stack->IsPhysicalPrimary( TMath::Abs ( label0 ) );
// 		Bool_t isPrimary1 = stack->IsPhysicalPrimary( TMath::Abs ( label1 ) );
		if ( label0 != label1 ) {
			hBkg[0]->Fill ( chi2,  spdmultiplicity->GetEta(iTracklet) );
			if ( fPhiCut && fPhiCut->CheckCut( spdmultiplicity->GetPhi(iTracklet) ) ) {
				hBkg[1]->Fill( chi2, spdmultiplicity->GetEta(iTracklet) );
			}
			return kTRUE;
		}
	}
	return kFALSE;
}

TH1D* AliMultiplicityEtaPhiSelection::GetBkgFraction ( Double_t chi2cut, Bool_t phiCut, Bool_t rebinned ) {
	Int_t i = 0; if ( phiCut ) i = 1;
	Int_t chi2bin = hEtaChi2[i]->GetXaxis()->FindBin(chi2cut);
	TH1D* hBkgf = hBkg[i]->ProjectionY("_py1",0,chi2bin);
	TH1D* hAll = hEtaChi2[i]->ProjectionY("_py2",0,hEtaChi2[i]->GetNbinsX());
	
	if ( rebinned ) {
		const Int_t netabins = 21;
		Double_t etabinsize = 0.2;
		Double_t* etabins = new Double_t[netabins+1];
		for (Int_t j = 0; j <= netabins; ++j ) {
			etabins[j] = -2.3 + ( j + 1 ) * etabinsize;
		}
		TH1D* hBkgf_r = (TH1D*) hBkgf->Rebin(netabins,"hBkgf_r",etabins);
		TH1D* hAll_r = (TH1D*) hAll->Rebin(netabins,"hAll_r",etabins);
		hBkgf_r->Divide(hAll_r);
		delete hBkgf;
		hBkgf = hBkgf_r;
	} else {
		hBkgf->Divide(hAll);
	}
	return hBkgf;
}
