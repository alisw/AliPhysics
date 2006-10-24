void AliKStarAnalysis
(
	const char *path = "/home/pulvir/resonances/aliroot-v4-04-Rev-08/pythia/selections",
	Double_t ptMin = 0.0, 
	Double_t ptMax = 0.0
)
{
	gSystem->Load("libPWG2.so");
	
	TChain *tree = new TChain("selection");
	
	// Open the working directory
	void *dirp = gSystem->OpenDirectory(path);
	const char *name = 0x0;
	
	// Add all files matching *pattern* to the chain
	while((name = gSystem->GetDirEntry(dirp))) {
		if (strstr(name, ".root")) {
			cout << "Adding " << name << endl;
			tree->Add(Form("%s/%s", path, name));
		}
	}
	gSystem->FreeDirectory(dirp);
	
	// assign working parameters
	AliRsnAnalysis *analysis = new AliRsnAnalysis;
	analysis->SetEventsTree(tree);
	if (ptMin < ptMax && ptMax != 0.0) {
		analysis->SetPtBin(ptMin, ptMax);
	}
	
	// set histogram bins
	analysis->SetBins(300, 0.5, 2.0);  // 300 bins of 5 MeV each
	analysis->SetTrueMotherPDG(313);   // PDG code of K*
		
	// import all useful combinations
	analysis->AddPairDef(AliPID::kPion, '+', AliPID::kKaon, '-');
	analysis->AddPairDef(AliPID::kPion, '-', AliPID::kKaon, '+');
	analysis->AddPairDef(AliPID::kPion, '+', AliPID::kKaon, '-', kTRUE);
	analysis->AddPairDef(AliPID::kPion, '-', AliPID::kKaon, '+', kTRUE);
	
	analysis->AddMixPairDef(AliPID::kPion, '+', AliPID::kKaon, '-');
	analysis->AddMixPairDef(AliPID::kPion, '-', AliPID::kKaon, '+');
	analysis->AddMixPairDef(AliPID::kKaon, '-', AliPID::kPion, '+');
	analysis->AddMixPairDef(AliPID::kKaon, '+', AliPID::kPion, '-');
	
	// process data for signal event in same event
	analysis->Process();
	analysis->EventMix(5, 5, 0.02, kFALSE);
	
	// open output file
	TFile *outFile = 0;
	if (ptMin == ptMax && ptMin == 0.) {
		outFile = TFile::Open("kstar.invmass.root", "RECREATE");
	}
	else {
		outFile = TFile::Open(Form("kstar.invmass.pt%3.1f-%3.1f.root", ptMin, ptMax), "RECREATE");
	}
	
	// write histograms
	analysis->WriteHistograms();
	
	outFile->Close();
}
