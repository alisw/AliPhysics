void AliITSMerge(Int_t Nfiles=1,const char* file0="galice.root",
		 const char* file1="galice_bg.root"){
    //
 
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(file0);
    if (file) {file->Close(); delete file;}
    cout << "AliITSMerge" << endl;
    file = new TFile(file0,"UPDATE");
    if (!file->IsOpen()) {
        cerr<<"Can't open "<<file0<<" !" << endl;
        return;
    } // end if !file
    file->ls();
 
    // Get AliRun object from file or return if not on file
    if (gAlice) delete gAlice;
    gAlice = (AliRun*)file->Get("gAlice");
    if (!gAlice) {
        cerr << "AliITSMerge.C : AliRun object not found on file" << endl;
        return;
    } // end if !gAlice

    if(Nfiles>2) Nfiles = 2;
    if(Nfiles<1) Nfiles = 1;
    AliRunDigitizer *manager = new AliRunDigitizer(Nfiles,1); 
    manager->SetInputStream(0,file0);
    if(Nfiles>1) manager->SetInputStream(1,file1);
    manager->SetOutputFile(file0);
    AliITSDigitizer *dITS = new AliITSDigitizer(manager);
    manager->Exec("");
}
