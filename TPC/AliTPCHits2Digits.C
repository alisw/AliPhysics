void TPCHits2Digits()
{
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }  
  gROOT->LoadMacro("SetTPCParam.C");
  AliTPCParam *par=SetTPCParam();


  // Connect the Root Galice file containing Geometry, Kine and Hits
  const char * inFile = "galice.root";  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if (file) file->Close();
  file = new TFile(inFile,"UPDATE");
  // Get AliRun object from file or create it if not on file

  gAlice = (AliRun*)file->Get("gAlice");
  if (gAlice) printf("AliRun object found on file\n");
  if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");

  gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      
    
  cerr<<"Hits2Digits\n";
  //setup TPCDigitsArray 
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(par);
  TPC->SetParam(par);
  arr->MakeTree();

  TPC->SetDigitsArray(arr);
  TPC->Hits2DigitsSector(1);             
  TPC->Hits2DigitsSector(2);             
  TPC->Hits2DigitsSector(3);             
  TPC->Hits2DigitsSector(1+18);             
  TPC->Hits2DigitsSector(2+18);             
  TPC->Hits2DigitsSector(3+18);             

  TPC->Hits2DigitsSector(36+1);             
  TPC->Hits2DigitsSector(36+2);             
  TPC->Hits2DigitsSector(36+3);             
  TPC->Hits2DigitsSector(36+1+18);             
  TPC->Hits2DigitsSector(36+2+18);             
  TPC->Hits2DigitsSector(36+3+18);             
  //write results
  char treeName[100];
  sprintf(treeName,"TreeD_%s",par->GetTitle());
  TPC->GetDigitsArray()->GetTree()->Write(treeName);
  par->Write(par->GetTitle());
  file->Close();
};

