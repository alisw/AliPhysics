void centua1(char* bin)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libJETAN");
  
  char dir[100];
  //sprintf(dir, "/data1/jgcn/partonic_events/unquenched/%s", bin);
  sprintf(dir, "/home/guest/alice/data/cent1_nq/%s", bin);

  AliJetESDReaderHeader *krh = new AliJetESDReaderHeader(); // hijing
  //AliJetKineReaderHeader *krh = new AliJetKineReaderHeader(); // pythia
  krh->SetComment("");
  krh->SetDirectory(dir);
  krh->SetPattern("miniesd"); // hijing
  krh->SetFirstEvent(0);
  krh->SetLastEvent(100000);
  krh->SetPtCut(1.0);
  //krh->SetFastSimTPC(kFALSE); // pythia
  //krh->SetFastSimEMCAL(kFALSE); // pythia
 
  // define reader and set its header
  AliJetESDReader *kr = new AliJetESDReader(); // hijing
  //AliJetKineReader *kr = new AliJetKineReader(); // pythia
  kr->SetReaderHeader(krh);
 
  // define jet header
  AliUA1JetHeader *jh=new AliUA1JetHeader();
  jh->SetComment("UA1 jet code with radius 1");
  jh->SetMode(1);
  jh->SetRadius(0.4);
  jh->SetMinCellEt(0.);
  jh->SetEtSeed(4.);
  jh->SetLegoNbinPhi(420.);
  jh->SetLegoNbinEta(120.);
  jh->SetLegoEtaMin(-0.9);
  jh->SetLegoEtaMax(+0.9);
    
  // define jet finder. Set its header and reader
  AliUA1JetFinder *jf = new AliUA1JetFinder();
  jf->SetJetHeader(jh);
  jf->SetJetReader(kr);
  jf->SetPlotMode(kTRUE);
  jf->SetOutputFile("jets.root");
  // do the job
  jf->Run();
}
 
