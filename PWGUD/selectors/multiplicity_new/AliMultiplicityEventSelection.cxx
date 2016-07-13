#include "AliMultiplicityEventSelection.h"

ClassImp ( AliMultiplicityEventSelection )

void AliMultiplicityEventSelection::CreateHistograms() {
	CreateStats();
	CreateSelectionHistograms();
}


AliMultiplicityEventSelection::AliMultiplicityEventSelection ( const char* name, const char* title ) : TNamed ( name, title ),
	hStats ( 0 ),
	sampleDef ( 0 ),
	currentMask ( 0 ),
	fFollowESDselection(kFALSE),
	fFollowMCselection(kFALSE),
	fZcacheGen ( 0 )
{
	for (Int_t i = 0; i < 3; ++i) {
		fZcache[i] = 0;
	}
	SetCollectEvent(kFALSE);
}

void AliMultiplicityEventSelection::SetNeedsMC ( Bool_t needs ) {
	fNeedsMC = needs;
}

Bool_t AliMultiplicityEventSelection::GetNeedsMC() {
	return fNeedsMC;
}

Bool_t AliMultiplicityEventSelection::IsCurrentESDEventAccepted() {
	return fCurrentESDEventAccepted;
}

Bool_t AliMultiplicityEventSelection::IsCurrentMCEventAccepted() {
	return fCurrentMCEventAccepted;
}

void AliMultiplicityEventSelection::ResetAccepted() {
	SetCollectEvent(kFALSE);
	fCurrentESDEventAccepted = kFALSE;
	fCurrentMCEventAccepted = kFALSE;
}

void AliMultiplicityEventSelection::CreateStats() {
	hStats = new TH1D ( "hStats", "Statistics", GetNeedsMC() ? 4 : 2, 0.5, GetNeedsMC() ? 2.5 : 4.5 );
	hStats->GetXaxis()->SetBinLabel ( 1, "Accepted" );
	hStats->GetXaxis()->SetBinLabel ( 2, "Rejected" );
	if ( GetNeedsMC() ) {
		hStats->GetXaxis()->SetBinLabel ( 3, "Accepted MC" );
		hStats->GetXaxis()->SetBinLabel ( 4, "Rejected MC" );
	}
	hStats->GetYaxis()->SetTitle ( "Count" );
}

void AliMultiplicityEventSelection::UpdateStats() {
	if ( IsCurrentESDEventAccepted() ) {
		hStats->Fill ( "Accepted", 1 );
	} else {
		hStats->Fill ( "Rejected", 1 );

	}
	if ( GetNeedsMC() && IsCurrentMCEventAccepted() ) {
		hStats->Fill ( "Accepted MC", 1 );
	} else {
		hStats->Fill ( "Rejected MC", 1 );
	}
}

Long64_t AliMultiplicityEventSelection::MergeStats ( const TCollection* mergelist ) {
// 	Printf("STATS MERGING CALLED.");
	if ( !mergelist ) {
		AliDebug ( 1, "NULL merge list." );
		return 0;
	}

	if ( mergelist->IsEmpty() ) {
		AliDebug ( 1, "Empty merge list." );
		return 1;
	}

	TList histograms;
// 	histograms.SetOwner();
	Int_t count = 0;

	TIter nextE ( mergelist );
	while ( AliMultiplicityEventSelection* selection = dynamic_cast<AliMultiplicityEventSelection*> ( nextE() ) ) {
		if ( !selection ) continue;
		histograms.Add ( selection->hStats );
		count++;
	}

	hStats->Merge ( &histograms );
	return ( Long64_t ) count + 1;
}

void AliMultiplicityEventSelection::SaveHistograms() {
	TFile* outfile = TFile::Open ( TString::Format ( "%s.root", GetName() ), "RECREATE" );
	AliInfoF ( "Saving histograms of %s", GetTitle() );
	hStats->Write();
	SaveSelectionHistograms();
	outfile->Close();
}

Bool_t AliMultiplicityEventSelection::GetSaveHistograms() {
	return fSaveHistograms;
}

void AliMultiplicityEventSelection::SetSaveHistograms ( Bool_t save ) {
	fSaveHistograms = save;
}

TCanvas* AliMultiplicityEventSelection::GetCanvas() {
	TCanvas* cnv = new TCanvas ( TString::Format ( "%s_canvas", GetName() ), TString::Format ( "%s canvas", GetTitle() ), 1280, 800 );
	return cnv;
}

Bool_t AliMultiplicityEventSelection::GetBatchMode() {
	return fBatchMode;
}

void AliMultiplicityEventSelection::SetBatchMode ( Bool_t batch ) {
	fBatchMode = batch;
}


Bool_t AliMultiplicityEventSelection::GetCorrelate() {
	return fCorrelate;
}

void AliMultiplicityEventSelection::SetCorrelate ( Bool_t correlate ) {
	fCorrelate = correlate;
}

Bool_t AliMultiplicityEventSelection::CollectEvent() {
	return fCollectEvent;
}

void AliMultiplicityEventSelection::SetCollectEvent ( Bool_t collect ) {
	fCollectEvent = collect;
}
