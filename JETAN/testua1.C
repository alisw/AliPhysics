void testua1()
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libJETAN");

  // define reader header
  AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
  jrh->SetComment("testing");
  //jrh->SetDirectory("rfio:///castor/cern.ch/user/p/phristov/cent1_nq/104-125GeV");
  jrh->SetDirectory("/home/morsch/jets_104_125");
  jrh->SetPattern("00165");
  jrh->SetFirstEvent(0);
  jrh->SetLastEvent(1000);
  jrh->SetPtCut(2.);
  jrh->SetReadSignalOnly(kFALSE);

  AliJetKineReaderHeader *krh = new AliJetKineReaderHeader();
  krh->SetComment("testing");
  krh->SetDirectory("/home/morsch/AliRoot/newio/FASTSIM/ctest/100/uq");
  krh->SetFirstEvent(0);
  krh->SetLastEvent(1000);
  krh->SetPtCut(2.);
  krh->SetFastSimTPC(kTRUE);
  

  // define reader and set its header
  AliJetESDReader *jr = new AliJetESDReader();
  jr->SetReaderHeader(jrh);

  // define reader and set its header
  AliJetKineReader *kr = new AliJetKineReader();
  kr->SetReaderHeader(krh);

  // define jet header
  AliUA1JetHeader *jh=new AliUA1JetHeader();
  jh->SetComment("UA1 jet code with default parameters");
  jh->SetMode(0);
  jh->SetRadius(0.4);
  jh->SetMinCellEt(0.);
  jh->SetEtSeed(4.);
  jh->SetNbinPhi(420.);
  jh->SetNbinEta(120.);
  jh->SetEtaMin(-0.9);
  jh->SetEtaMax(+0.9);  
  

  // define jet finder. Set its header and reader
  AliUA1JetFinder *jf = new AliUA1JetFinder();
  jf->SetJetHeader(jh);
  jf->SetJetReader(jr);
  jf->SetPlotMode(kTRUE);
  jf->SetOutputFile("jets.root");
  // do the job
  jf->Run();
}
