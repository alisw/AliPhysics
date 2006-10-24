void AliLambdaStarAnalysis
(
	Int_t nevents = 650,
	const char *path_in = "myselection",
	const char *file_out = "lambda.invmass.root"
)
{
	gSystem->Load("libANALYSIS.so");
	
	TChain *tree = new TChain("selection");
	
	if (nevents <= 0) {
		// Open the working directory
		void *dirp = gSystem->OpenDirectory(Form("$HOME/lambda/%s", path_in));
		const char *name = 0x0;
		while((name = gSystem->GetDirEntry(dirp))) {
			if (strstr(name, ".root")) {
				tree->Add(Form("$HOME/lambda/%s/%s", path_in, name));
			}
		}
		gSystem->FreeDirectory(dirp);
	}
	else {
		for (Int_t iev = 0; iev <= nevents; iev++) {
			tree->Add(Form("$HOME/lambda/%s/%d.root", path_in, iev));
		}
	}
	cout << "# events: " << tree->GetEntries() << endl;
	
	//AliRsnDaughterCutPtSingle *cutproton = new AliRsnDaughterCutPtSingle(0, 1);
	//AliRsnDaughterCutPtSingle *cutkaon = new AliRsnDaughterCutPtSingle(0, 1);
	
	// assign working parameters
	AliRsnAnalysis *analysis = new AliRsnAnalysis;
	analysis->SetEventsTree(tree);
	
	//analysis->AddCutSingle(AliPID::kProton, cutproton);
	//analysis->AddCutSingle(AliPID::kKaon, cutkaon);
	
	// set histogram bins
	analysis->SetBins(800, 1.3, 2.1);   // 700 bins of 1 MeV each
	analysis->SetTrueMotherPDG(3124);   // PDG code of Lambda*
		
	// import all useful combinations
	
	analysis->AddPairDef(AliPID::kProton, '+', AliPID::kKaon, '-');
	analysis->AddPairDef(AliPID::kProton, '-', AliPID::kKaon, '+');
	analysis->AddPairDef(AliPID::kProton, '+', AliPID::kKaon, '-', kTRUE);
	analysis->AddPairDef(AliPID::kProton, '-', AliPID::kKaon, '+', kTRUE);
	
	analysis->AddMixPairDef(AliPID::kProton, '+', AliPID::kKaon, '-');
	analysis->AddMixPairDef(AliPID::kProton, '-', AliPID::kKaon, '+');
	analysis->AddMixPairDef(AliPID::kKaon, '-', AliPID::kProton, '+');
	analysis->AddMixPairDef(AliPID::kKaon, '+', AliPID::kProton, '-');
	
	

	// process data for signal event in same event
	analysis->Process();
	analysis->EventMix(5, 5, 0.02, kFALSE);
	// open output file
	TFile *outFile = 0;
	outFile = TFile::Open(file_out, "RECREATE");
	
	// write histograms
	analysis->WriteHistograms();
	
	outFile->Close();
}
