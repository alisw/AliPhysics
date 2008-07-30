//
// This macro prints all details of each track for each AliRsnEvent in a file
// which is assumed to be saved in the same directory where this is executed.
//

void AliRsnPrintTracks(const char *fileName = "AliRsnEventsESD.root")
{
    // load libraries
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libPWG2resonances.so");
    
    TFile *file = TFile::Open(fileName);
    TTree *tree = (TTree*)file->Get("aodTree");
    
    Int_t i, last, nEvents = (Int_t)tree->GetEntries();
    
    AliRsnEvent *event = 0;
    tree->SetBranchAddress("rsnEvents", &event);
    
    for (i = 0; i < nEvents; i++) {
        cout << "*** Event " << i << " ***" << endl;
        tree->GetEntry(i);
        event->Print("P");
        last = event->GetLastFastTrack(1.0);
        cout << endl;
        if (last > 0) {
            cout << "...Tracks from 0 to " << last << " have Pt > 1.0 GeV" << endl << endl;
        }
        else {
            cout << "...NO tracks with Pt > 1.0 GeV" << endl << endl;
        }
    }
}
