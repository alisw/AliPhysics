void testpxcone()
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libJETAN");

  // define reader header
  AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
  jrh->SetComment("testing");
  jrh->SetDirectory("/home/work/alice/data/hij_central");
  jrh->SetPattern("run00");
  jrh->SetFirstEvent(0);
  jrh->SetLastEvent(1);
  jrh->SetPtCut(2);

  // define reader and set its header
  AliJetESDReader *jr = new AliJetESDReader();
  jr->SetReaderHeader(jrh);

  // define jet header, take the defaults
  AliPxconeJetHeader *jh=new AliPxconeJetHeader();
  jh->SetComment("Testing jet code");

  // define jet finder. Set its header and reader
  AliPxconeJetFinder *jf = new AliPxconeJetFinder();
  jf->SetJetHeader(jh);
  jf->SetJetReader(jr);
  jf->SetPlotMode(kTRUE);
  jf->SetOutputFile("pxcone.root");
  // do the job
  jf->Run();
}
