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
  TF1 f("f", "gaus", -100, 100);
  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  TCanvas *c = new TCanvas("c", "Trend Distrib.", 10, 10, 500, 500);
  Int_t ntr=tDB->GetEntries();
  for(Int_t it(0); it<nt; it++){
    tDB->Draw(tvn[it][0], "", "goff");
    Double_t *v = tDB->GetV1(), xmin(100.), xmax(-100);
    for(Int_t ir=0; ir<ntr; ir++){
      if(v[ir]<-100) continue;
      if(v[ir]<xmin) xmin = v[ir];
      if(v[ir]>xmax) xmax = v[ir];
    }
    TH1 *h = new TH1F("h", Form(";%s;entries", tvn[it][0]), 10, 0.5*(3*xmin-xmax), 0.5*(3*xmax - xmin));
    tDB->Draw(Form("%s>>h", tvn[it][0]), Form("%s>-100", tvn[it][0]));
    if(h->Integral() < 1) continue;
    f.SetParameter(0, h->Integral());
    f.SetParameter(1, h->GetMean());
    f.SetParameter(2, h->GetRMS());
    h->Fit(&f, "WQ");
    c->Modified(); c->Update(); c->SaveAs(Form("Fit_%s.gif", tvn[it][0]));

    // write trending value to manager
    Info("makeTrendingDB", "%s [%f - %f] %f[%f]", tvn[it][0], xmin, xmax, f.GetParameter(1), f.GetParameter(2));
    tm->AddValue(tvn[it][0], f.GetParameter(1), f.GetParameter(2),
      tvn[it][1], res[it>13], notifiable);
    delete h;
  }
  tm->Terminate();

  fDB->cd();
  tDB->Write();
  fDB->Close();
}
