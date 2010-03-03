AliJetCorrelSelector* ConfigJetCorrel(){

  ///////////////////////////////////
  // set correlation input parameters
  ///////////////////////////////////
  // set generic selections:
  Bool_t kGenQA = kTRUE;       // generate QA histos
  UInt_t kDPhiNumBins = 60;    // number of bins in DeltaPhi histos
  UInt_t kDEtaNumBins = 40;    // number of bins in DeltaEta histos
  Bool_t kUseAliKF = kFALSE;   // use AliKF or TLorentzVector for parent reconstruction
  UInt_t poolDepth = 10;
  UInt_t correlTypes[] = {0};  // 0=dihadron, 1=pi0-hadron, 2=photon-hadron
  Float_t centrBins[] = {1,300};
  Float_t zVertBins[] = {-10,-7,-5,5,7,10};
  Float_t bwTriggPt = 1;   Float_t minTriggPt = 2;   Float_t maxTriggPt = 10;
  Float_t bwAssocPt = 0.5; Float_t minAssocPt = 0.5; Float_t maxAssocPt = 4;
  //TString sTrigg[] = {"ALL"}; // selects events where one of the strings is matched; "ALL"=no cut
  TString sTrigg[] = {"CINT1B-"};
  // set track selections:
  Bool_t itsRefit = kTRUE;
  Bool_t tpcRefit = kTRUE;
  Bool_t trdRefit = kTRUE;         // used only for electron tracks
  Float_t maxEta = 0.8;
  UInt_t minNClusTPC = 70;
  Float_t maxTPCChi2 = 4.0;        // max track Chi2 per TPC cluster
  Bool_t rejectKinkChild = kTRUE;  // reject track comming from a kink
  Float_t trkPairCut = 0.;         // track pair proximity cut (dist at TPC entrance)
  // code that applies next 3 cuts (NClusITS,ITSChi2,NsigVtx) currently commented out
  UInt_t minNClusITS = 0;
  Float_t maxITSChi2 = 35;        // max track Chi2 per ITS cluster
  Float_t maxNsigVtx = 35;        // max dist to primary vertex (sigma)
  Float_t maxTrkVtx  = 2.4;       // max dist to primary vertex (absolute) - temporarily instead of sigma

  //////////////////////////////////
  // load them into selector object:
  //////////////////////////////////
  AliJetCorrelSelector* selector = new AliJetCorrelSelector();
  selector->SetQA(kGenQA);
  selector->SetUseAliKF(kUseAliKF);
  selector->SetDPhiNumBins(kDPhiNumBins);
  selector->SetDEtaNumBins(kDEtaNumBins);
  selector->SetPoolDepth(poolDepth);
  selector->SetCorrelTypes(sizeof(correlTypes)/sizeof(UInt_t),correlTypes);
  selector->SetBinningCentr(sizeof(centrBins)/sizeof(Float_t),centrBins);
  selector->SetBinningZvert(sizeof(zVertBins)/sizeof(Float_t),zVertBins);
  selector->SetBinningTrigg(minTriggPt,maxTriggPt,bwTriggPt);
  selector->SetBinningAssoc(minAssocPt,maxAssocPt,bwAssocPt);
  selector->SetTriggers(sizeof(sTrigg)/sizeof(TString),sTrigg);
  selector->SetITSRefit(itsRefit);
  selector->SetTPCRefit(tpcRefit);
  selector->SetTRDRefit(trdRefit);
  selector->SetMaxEta(maxEta);
  selector->SetMinNClusITS(minNClusITS);
  selector->SetMinNClusTPC(minNClusTPC);
  selector->SetMaxITSChi2(maxITSChi2);
  selector->SetMaxTPCChi2(maxTPCChi2);
  selector->SetMaxNsigmaVtx(maxNsigVtx);
  selector->SetMaxTrkVtx(maxTrkVtx);
  selector->SetRejectKinkChild(rejectKinkChild);
  selector->SetTrkProximityCut(trkPairCut);
  selector->Show();

  return selector;
}
