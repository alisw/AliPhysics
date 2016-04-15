Bool_t ProcessOutputCheb(TString filesToProcess, Int_t startRun, Int_t endRun, const char* ocdbStorage) {

  // macro that process a list of files (xml or txt) to produce then the
  // OCDB entry for the TPC SP Distortion calibration; inspired by
  // AliAnalysisAlien::MergeOutput for simplicity

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateCorrMapObj.C");
  //gROOT->LoadMacro("CreateCorrMapObj.C");
    
  Bool_t isGrid = kTRUE;
  TObjArray *listoffiles = new TObjArray();
  
  if (filesToProcess.Contains(".xml")) {
    // Merge files pointed by the xml 
    TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\");", filesToProcess.Data()));
    if (!coll) {
      ::Error("ProcessOutput", "Input XML collection empty.");
      return kFALSE;
    }
    coll->Print();
    // Iterate grid collection
    while (coll->Next()) {
      TString fname = coll->GetTURL();
      Printf("fname = %s", fname.Data());
      listoffiles->Add(new TNamed(fname.Data(),""));
    }   
  }
  
  else if (filesToProcess.Contains(".txt")) {
    TString line;
    ifstream in;
    in.open(filesToProcess);
    if (in.fail()) {
      ::Error("ProcessOutput", "File %s cannot be opened. Processing stopped." ,filesToProcess.Data());
      return kTRUE;
    }
    Int_t nfiles = 0;
    while (in.good()) {
      in >> line;
      if (line.IsNull()) continue;
      if (!line.Contains("alien:")) isGrid = kFALSE;
      nfiles++;
      listoffiles->Add(new TNamed(line.Data(),""));
    }
    in.close();
    if (!nfiles) {
      ::Error("ProcessOutput","Input file %s contains no files to be processed\n", filesToProcess.Data());
      delete listoffiles;
      return kFALSE;
    }
  }
  
  if (!listoffiles->GetEntries()) {
    ::Error("ProcessOutput","No files to process\n");
    delete listoffiles;
    return kFALSE;
  }
  
  if (startRun != endRun) {
    Printf("The processing now is only run-level, please check again!");
    return kFALSE;
  }
  
  Int_t run = startRun;
  TObject *nextfile;
  TString snextfile;
  TIter next(listoffiles);   
  TObjArray* a = new TObjArray();
  a->SetOwner(kTRUE);
  while (nextfile=next()) {
    snextfile = nextfile->GetName();
    Printf("opening file %s", snextfile.Data());
    if (isGrid) TGrid::Connect("alien://");
    TFile* ftmp = TFile::Open(Form("%s", snextfile.Data()));
    TList* l = (TList*)ftmp->GetListOfKeys();
    for (Int_t ikey = 0; ikey < l->GetEntries(); ikey++){
      TKey* k = l->At(ikey);
      TString kStr = k->GetName();
      if (kStr.Contains(Form("run%d", run))){
	AliTPCDcalibRes* dcalibRes = (AliTPCDcalibRes*)ftmp->Get(Form("%s", kStr.Data()));
	AliTPCChebCorr* c = dcalibRes->GetChebCorrObject();
	a->Add(c);
      }
    }
  }
  
  a->Print();
  CreateCorrMapObjTime(a, startRun, endRun, ocdbStorage);
  delete a;
  delete listoffiles;
  
  return kTRUE;
  
}
