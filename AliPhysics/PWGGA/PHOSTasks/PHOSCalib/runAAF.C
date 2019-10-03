void runAAF()
{
  // TProof::Open("polishch:default@alice-caf.cern.ch", "masteronly"); // только мастер
  //TProof::Open("prsnko:default@ccapl0012","workers=2");             //  2 workers
//  TProof::Open("prsnko:default@alice-caf.cern.ch","workers=32");                                        //  all workers (48)
//  TProof::Open("prsnko:default@alice-caf.cern.ch","workers=29");                                        //  all workers (48)
  TProof::Open("prsnko:default@alice-caf.cern.ch");                                        //  all workers (48)
  //TProof::Open("prsnko:default@ccapl0012", "masteronly");                                        //  all workers (48)

  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice"));
  //list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice:PWGGAPHOSTasks"));

  gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-25-AN",list);
} 
