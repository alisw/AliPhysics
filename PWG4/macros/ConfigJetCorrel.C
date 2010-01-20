AliJetCorrelSelector* JetCorrelSelector(){

  ///////////////////////////////////
  // set correlation input parameters
  ///////////////////////////////////
  // set generic selections:
  UInt_t PoolDepth = 10;
  UInt_t CorrelTypes[] = {0};
  Float_t TriggBins[] = {5.,7.,10.,15.,25.};
  Float_t AssocBins[] = {0.3,0.5,1.,2.,5.,7.};
  Float_t CentrBins[] = {0.,50.,200.,500.};
  Float_t ZVertBins[] = {-30.,-15.,-5.,-1.,1.,5.,15.,30.};
  // set track selections:
  Bool_t ITSRefit = kTRUE;
  Bool_t TPCRefit = kTRUE;
  Bool_t TRDRefit = kTRUE;          // used only for electron tracks
  UInt_t MinNClusITS = 1;
  UInt_t MinNClusTPC = 50; 
  Float_t MaxITSChi2 = 3.5;         // max track Chi2 per ITS cluster
  Float_t MaxTPCChi2 = 3.5;         // max track Chi2 per TPC cluster
  Float_t MaxNsigVtx = 3.5;         // max dist to primary vertex
  Bool_t RejectKinkChild = kTRUE;   // reject track comming from a kink

  //////////////////////////////////
  // load them into selector object:
  //////////////////////////////////
  AliJetCorrelSelector* Selector = new AliJetCorrelSelector();
  Selector->SetPoolDepth(PoolDepth);
  Selector->SetCorrelTypes(sizeof(CorrelTypes)/sizeof(Int_t),CorrelTypes);
  Selector->SetBinningTrigg(sizeof(TriggBins)/sizeof(Float_t),TriggBins);
  Selector->SetBinningAssoc(sizeof(AssocBins)/sizeof(Float_t),AssocBins);
  Selector->SetBinningCentr(sizeof(CentrBins)/sizeof(Float_t),CentrBins);
  Selector->SetBinningZvert(sizeof(ZVertBins)/sizeof(Float_t),ZVertBins);
  Selector->SetITSRefit(ITSRefit); Selector->SetTPCRefit(TPCRefit);
  Selector->SetTRDRefit(TRDRefit);
  Selector->SetMinNClusITS(MinNClusITS); Selector->SetMinNClusTPC(MinNClusTPC);
  Selector->SetMaxITSChi2(MaxITSChi2); Selector->SetMaxTPCChi2(MaxTPCChi2);
  Selector->SetMaxNsigmaVtx(MaxNsigVtx);
  Selector->SetRejectKinkChild(RejectKinkChild);
  Selector->Print();

  return Selector;
}
