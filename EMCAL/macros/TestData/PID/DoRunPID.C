void DoRunPid()
{
//Example to use PID analysis with ESDs
//Prints in screen PID weight values

   TFile *file = TFile::Open("AliESDs.root");
   TTree *tree = (TTree*)file->Get("esdTree");

   AliESD *esd = 0;
   tree->SetBranchAddress("ESD", &esd);
   Int_t nEvents = (Int_t)tree->GetEntries();

   AliEMCALPID *pid = new AliEMCALPID;
   pid->SetPrintInfo(kTRUE);
   for (Int_t iev = 0; iev < nEvents: iev++) {
      tree->GetEntry(iev);
      pid->RunPID(esd);
   }

   file->Close();
} 
