Int_t AliTPCHits2Digits()
{

  // new version by J.Belikov

  // Connect the Root Galice file containing Geometry, Kine and Hits

  const char * inFile = "galice.root";  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if (file) {file->Close(); delete file;}
  file = new TFile(inFile,"UPDATE");
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<inFile<<" !\n";
    return 1;
  }

  // Get AliRun object from file or create it if not on file
  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }

  gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

//Set response functions

  AliTPCParamSR *param=(AliTPCParamSR*)gDirectory->Get("75x40_100x60");
  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter   = new AliTPCPRF2D;
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  rf->SetGauss(param->GetZSigma(),param->GetZWidth(),1.);
  rf->SetOffset(3*param->GetZSigma());
  rf->Update();

  TDirectory *savedir=gDirectory;
  TFile *f=TFile::Open("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  if (!f->IsOpen()) { 
     cerr<<"Can't open $ALICE_ROOT/TPC/AliTPCprf2d.root !\n"
     return 3;
  }
  prfinner->Read("prf_07504_Gati_056068_d02");
  prfouter->Read("prf_10006_Gati_047051_d03");
  f->Close();
  savedir->cd();

  param->SetInnerPRF(prfinner);
  param->SetOuterPRF(prfouter); 
  param->SetTimeRF(rf);
  TPC->SetParam(param);
   
  cerr<<"Digitizing TPC...\n";

  //setup TPCDigitsArray 
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(param);
  TPC->SetParam(param);
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
  sprintf(treeName,"TreeD_%s",param->GetTitle());
  TPC->GetDigitsArray()->GetTree()->Write(treeName);

  delete gAlice; gAlice=0;
  file->Close(); delete file;
  return 0;
};

