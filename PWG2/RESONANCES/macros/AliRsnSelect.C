void AliRsnSelect
(
	Int_t     event            = 1,
	Option_t *inRootPath       = "/home/pulvir/resonances/aliroot-v4-04-Rev-08/pythia/events",
	Option_t *outRootPath      = "/home/pulvir/resonances/aliroot-v4-04-Rev-08/pythia/selections",
	Double_t  ptLimit4PID      = 10000.0
)
{
	// load ANALYSIS library for output objects
	gSystem->Load("libPWG2.so");
	
	// instantiate reading manager
	AliRsnReader *reader = new AliRsnReader;
	
	// Define prior probabilities (only for ESD PID)
	reader->SetPriorProbability(AliPID::kElectron, 0.0339947); //0.0339947															
	reader->SetPriorProbability(AliPID::kMuon,     0.0192307); //0.0192307
	reader->SetPriorProbability(AliPID::kPion,     0.822957);  //0.822957 
	reader->SetPriorProbability(AliPID::kKaon,     0.0751355); //0.0751355
	reader->SetPriorProbability(AliPID::kProton,   0.0486821); //0.0486821
	reader->SetProbabilityThreshold(0.5);
	
	// Define PID method
	reader->SetPIDMethod(AliRsnReader::kPerfectPID);

	// create input/output names
	Text_t inPath[200], outFileName[200];
	sprintf(inPath, "%s/%d", inRootPath, event);
	sprintf(outFileName, "%s/%d.root", outRootPath, event);
				
	// do event reading
	cout << "Reading data in " << inPath << endl;
	TTree *events = reader->ReadTracksAndParticles(inPath, "R");
	if (!events) return;
	cout << endl;
	
	// open output file
	TFile *fileOut = TFile::Open(outFileName, "RECREATE");
	cout << "Saving  data in " << outFileName << endl;
	events->Write(events->GetName(), TObject::kOverwrite);
	fileOut->Close();
}
