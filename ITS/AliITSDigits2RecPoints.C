TFile* AccessFile(TString inFile="galice.root", TString acctype="R");

void AliITSDigits2RecPoints(TString inFile="galice.root", TString outFile="galice.root"){

  TFile *file;
  if(outFile.Data() == inFile.Data()){
    file = AccessFile(inFile,"U");
  }
  else {
    file = AccessFile(inFile);
  }
  
  TStopwatch timer;

  cout << "Creating reconstructed points from digits for the ITS..." << endl;
  const char *nulptr=0;
  AliITSreconstruction *itsr = new AliITSreconstruction(nulptr);
  if(outFile.Data() != inFile.Data())itsr->SetOutputFile(outFile);
  timer.Start();
  itsr->Init();
  itsr->Exec(); 
  timer.Stop(); 
  timer.Print();    
  delete itsr;
}

//-------------------------------------------------------------------
TFile * AccessFile(TString FileName, TString acctype){

  // Function used to open the input file and fetch the AliRun object

  if (gAlice) {delete gAlice; gAlice = 0;}
  TFile *retfil = 0;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
  if (file) {file->Close(); delete file; file = 0;}
  if(acctype.Contains("U")){
    file = new TFile(FileName,"update");
  }
  if(acctype.Contains("N") && !file){
    file = new TFile(FileName,"recreate");
  }
  if(!file) file = new TFile(FileName);   // default readonly
  if (!file->IsOpen()) {
	cerr<<"Can't open "<<FileName<<" !" << endl;
	return retfil;
  } 

  // Get AliRun object from file or return if not on file
  //  if (gAlice) {delete gAlice; gAlice = 0;}  
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
	cerr << "AliRun object not found on file"<< endl;
	return retfil;
  } 
  return file;
}
