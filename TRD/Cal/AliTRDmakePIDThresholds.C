//_______________________________________________________
void AliTRDmakePIDThresholds(TString filename){
  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libTRDqaRec.so")<0) return;

  const TString histnames[4] = {"fHistThreshLQ", 
                                "fHistThreshNN", 
                                "fHistPionEffLQ", 
                                "fHistPionEffNN"};
  const TString histtitles[4] = {"PID Thresholds for the 2D Likelihood Method", 
                                 "PID Thresholds for the Neural Network Method", 
                                 "Pion Efficiency for the 2D Likelihood Method", 
                                 "Pion Efficiency for the Neural Network Method"};
  const Int_t pos[4] = {3 + AliTRDpidChecker::kLQ, 
                        3 + AliTRDpidChecker::kNN, 
                        AliTRDpidChecker::kLQ,
                        AliTRDpidChecker::kNN};   
  
  AliTRDpidChecker pidchecker;
  pidchecker.SetDebugLevel(2);
  pidchecker.Load(filename.Data());
  pidchecker.PostProcess();
  TObjArray *fGraph = pidchecker.GetGraphs();

  // Save the thresholds
  TObjArray *histos = new TObjArray;
  TH1F *histo_thresh = 0x0;
  for(Int_t ihist = 0; ihist < 4; ihist++){
    g = (TGraphErrors*)fGraph->At(pos[ihist]);
    histo_thresh = CreateHistogram(histnames[ihist], histtitles[ihist]);
    CovertHisto(g, histo_thresh);
    histos->AddAt(histo_thresh, ihist);
  }

  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexander Wilk");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-21-01"); //root version
  metaData->SetComment("TRD PID thresholds based on 90\% electron efficiency");
  
  AliCDBId id("TRD/Calib/PIDThresholds", 0, AliCDBRunRange::Infinity()); 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT");
  if (!gStorLoc) {
    return;
  }
  gStorLoc->Put(histos, id, metaData); 

  return;
}

//_______________________________________________________
TH1F *CreateHistogram(TString &histname, TString &histtitle){
  const Float_t fTrackMomentum[AliTRDCalPID::kNMom + 1] = { 0.6,  0.8,  1.0,  1.5,  2.0, 3.0, 4.0,  5.0,  6.0,  8.0, 10.0, 11.0};
  return new TH1F(histname.Data(), histtitle.Data(), AliTRDCalPID::kNMom, fTrackMomentum);
}

//_______________________________________________________
void CovertHisto(TGraphErrors *g, TH1F *h){
  Double_t x, y, dy;
  for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++){
    g->GetPoint(ibin - 1, x, y);
    dy = g->GetErrorY(ibin -1);
    h->SetBinContent(ibin, y);
    h->SetBinError(ibin, dy);
  }
}

