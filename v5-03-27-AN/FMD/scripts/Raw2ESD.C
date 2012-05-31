void
Raw2ESD(const char* file="")
{
  AliCDBManager::Instance()->SetRun(0);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliGeomManager::LoadGeometry("geometry.root");

  AliRawReader* reader = 0;
  TString rawFile(file);
  if (rawFile.IsNull() && rawFile.EndsWith(".root")) 
      reader = new AliRawReaderRoot(rawFile.Data());
  else if (!rawFile.IsNull() && rawFile.EndsWith(".raw"))
    reader = new AliRawReaderDate(rawFile.Data());
  else
    reader = new AliRawReaderFile(-1);
  
  AliFMDReconstructor* reco = new AliFMDReconstructor();
  reco->Init();

  Int_t        event       = 0;
  TFile*       digitFile   = TFile::Open("reco_digits.rot", "RECREATE");
  TTree*       digitTree   = new TTree("digit", "FMD digits");

  TFile*       clusterFile   = TFile::Open("FMD.RecPoints.root", "RECREATE");
  TTree*       clusterTree = new TTree("cluster", "FMD digits");

  TFile*       esdFile     = TFile::Open("AliESDs.root", "RECREATE");
  TTree*       esdTree     = new TTree("esdTree", "ESD Treee");
  AliESDEvent* esd         = new AliESDEvent();
  esd->CreateStdContent();
  esd->WriteToTree(esdTree);
  while ((reader && reader->NextEvent())) {
    reco->ConvertDigits(reader, digitTree);
    reco->Reconstruct(digitTree, clusterTree);

    esd->SetRunNumber(0);
    esd->SetEventNumberInFile(event);
    reco->FillESD((TTree*)0, (TTree*)0, esd);
    esdTree->Fill();
    esd->Reset();
    
    event++;    
  }
  esdFile->Write();
}

  
