void makeTrendingDB(const Char_t *fl)
{
// Make trending of variable list "tl" from trending file list "fl"
// The trending value list should be formated "var1:var2:var3"
// The trending file from the list should be found on a path formated "your path"/runId/TRD.PerformanceTrend.root 
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGmuon.so");

  const Int_t nt(45);
  const Char_t *tvn[nt][2] = {
    {"TRDcheckDET_NTracksEvent", "<N_{track}/Event>"},
    {"TRDcheckDET_NTracksEventRMS", "RMS(N_{track}/Event)"},
    {"TRDcheckDET_NTracksSector", "<N_{track}/Sector>"},
    {"TRDcheckDET_NClustersTrack", "<N_{cls}/Track>"},
    {"TRDcheckDET_NClustersTrackRMS", "RMS(N_{cls}/Track)"},
    {"TRDcheckDET_NClustersTracklet", "<N_{cls}/Tracklet>"},
    {"TRDcheckDET_NClustersTrackletRMS", "RMS(N_{cls}/Tracklet)"},
    {"TRDcheckDET_NTrackletsTrack", "<N_{tracklet}/Track>"},
    {"TRDcheckDET_NTrackletsTrackRMS", "RMS(N_{tracklet}/Track)"},
    {"TRDcheckDET_ChargeTracklet", "<dQdl>"},
    {"TRDcheckDET_ChargeTrackletRMS", "RMS(dQdl)"},
    {"TRDcheckDET_PHplateau", "Plateau(<PH>)"},
    {"TRDcheckDET_PHslope", "Slope(<PH>)"},
    {"TRDcheckDET_PHamplificationPeak", "Peak(<PH>)"},
//=======================================================
    {"TRDresolution_TrkInY", "#Deltay (TrkIn)"},
    {"TRDresolution_TrkInYn", "#Deltay^{-} (TrkIn)"},
    {"TRDresolution_TrkInYp", "#Deltay^{+} (TrkIn)"},
    {"TRDresolution_TrkInRCZ", "#Deltaz (TrkIn)"},
    {"TRDresolution_TrkInPhn", "#Delta#phi^{-} (TrkIn)"},
    {"TRDresolution_TrkInPhp", "#Delta#phi^{+} (TrkIn)"},
    {"TRDresolution_TrkInQp0", "MPV(dQdl^{e+}) (TrkIn)"},
    {"TRDresolution_TrkInQp1", "MPV(dQdl^{#mu+}) (TrkIn)"},
    {"TRDresolution_TrkInQp2", "MPV(dQdl^{#pi+}) (TrkIn)"},
    {"TRDresolution_TrkInQp3", "MPV(dQdl^{K+}) (TrkIn)"},
    {"TRDresolution_TrkInQp4", "MPV(dQdl^{p+}) (TrkIn)"},
    {"TRDresolution_TrkInQn0", "MPV(dQdl^{e-}) (TrkIn)"},
    {"TRDresolution_TrkInQn1", "MPV(dQdl^{#mu-}) (TrkIn)"},
    {"TRDresolution_TrkInQn2", "MPV(dQdl^{#pi-}) (TrkIn)"},
    {"TRDresolution_TrkInQn3", "MPV(dQdl^{K-}) (TrkIn)"},
    {"TRDresolution_TrkInQn4", "MPV(dQdl^{p-}) (TrkIn)"},
    {"TRDresolution_TrkInQ0", "MPV(dQdl^{e}) (TrkIn)"},
    {"TRDresolution_TrkInQ1", "MPV(dQdl^{#mu}) (TrkIn)"},
    {"TRDresolution_TrkInQ2", "MPV(dQdl^{#pi}) (TrkIn)"},
    {"TRDresolution_TrkInQ3", "MPV(dQdl^{K}) (TrkIn)"},
    {"TRDresolution_TrkInQ4", "MPV(dQdl^{p}) (TrkIn)"},
    {"TRDresolution_TrkInQSp0", "<dQdl^{e+}> (TrkIn)"},
    {"TRDresolution_TrkInQSp1", "<dQdl^{#mu+}> (TrkIn)"},
    {"TRDresolution_TrkInQSp2", "<dQdl^{pi+}> (TrkIn)"},
    {"TRDresolution_TrkInQSp3", "<dQdl^{K+}> (TrkIn)"},
    {"TRDresolution_TrkInQSp4", "<dQdl^{p+}> (TrkIn)"},
    {"TRDresolution_TrkInQSn0", "<dQdl^{e-}> (TrkIn)"},
    {"TRDresolution_TrkInQSn1", "<dQdl^{#mu-}> (TrkIn)"},
    {"TRDresolution_TrkInQSn2", "<dQdl^{#pi-}> (TrkIn)"},
    {"TRDresolution_TrkInQSn3", "<dQdl^{K-}> (TrkIn)"},
    {"TRDresolution_TrkInQSn4", "<dQdl^{p-}> (TrkIn)"}
  };
  const char *resName[] = {"Markus Fasel", "Alexandru Bercuci"},
             *resMail[] = {"M.Fasel@gsi.de", "A.Bercuci@gsi.de"};
  const char *notName[] = {"Julian Book", "Hans Beck", "Ionut Arsene", "Raphaelle Bailache", "Christoph Blume"},
             *notMail[] = {"jbook@ikf.uni-frankfurt.de", "hbeck@ikf.uni-frankfurt.de", "I.C.Arsene@gsi.de", "R.Bailhache@gsi.de", "blume@ikf.uni-frankfurt.de"};

  TFile *fDB = TFile::Open("TRD.TrendDB.root", "RECREATE");
  TTree *tDB = new TTree("trend", "Reference Trend Values");
  Double_t val[nt];
  for(Int_t it(0); it<nt; it++) tDB->Branch(tvn[it][0], &val[it], Form("%s/D", tvn[it][0]));
  gROOT->cd();

  AliTRDtrendValue *tv(NULL);
  FILE *fp = fopen(fl, "rt");
  TString sfp;
  while(sfp.Gets(fp)){
    if(!TFile::Open(sfp.Data())) continue;
    for(Int_t it(0); it<nt; it++){
      val[it] = -999;
      if(!(tv = (AliTRDtrendValue*)gFile->Get(tvn[it][0]))) {
        Warning("makeTrendingDB()", "Missing %s", tvn[it][0]);
        continue;
      }
      val[it] = tv->GetVal();
    }
    gFile->Close();
    tDB->Fill();
  }


//   TFile *fDB = TFile::Open("TRD.TrendDB.root");
//   TTree *tDB = (TTree*)gFile->Get("trend");
//   Double_t val[nt];
//   for(Int_t it(0); it<nt; it++) tDB->SetBranchAddress(tvn[it][0], &val[it]);
//   gROOT->cd();

  TString res[] = {Form("%s/%s", resName[0], resMail[0]), Form("%s/%s", resName[1], resMail[1])};
  TString notifiable;
  for(Int_t inot(0); inot<5; inot++){
    notifiable+=notName[inot];
    notifiable+="/";
    notifiable+=notMail[inot];
    if(inot<4) notifiable+=",";
  }
  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  TCanvas *c = new TCanvas("c", "Trend Distrib.", 10, 10, 500, 500);
  for(Int_t it(0); it<nt; it++){
    tDB->Draw(tvn[it][0]);
    TH1 *h = (TH1*)gROOT->FindObject("htemp");
    h->Fit("gaus", "WQ");
    c->Modified(); c->Update(); c->SaveAs(Form("%s.gif", tvn[it][0]));

    // write trending value to manager
    TF1 *f = h->GetFunction("gaus");
    printf("%s %f[%f] %f[%f]\n", tvn[it][0], h->GetMean(), f->GetParameter(1), h->GetRMS(), f->GetParameter(2));
    tm->AddValue(tvn[it][0], h->GetMean()/*f->GetParameter(1)*/, h->GetRMS()/*f->GetParameter(2)*/,
      tvn[it][1], res[it>13], notifiable);
  }
  tm->Terminate();

  fDB->cd();
  tDB->Write();
  fDB->Close();
}
