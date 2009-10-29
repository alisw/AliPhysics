/////////////////////////////////////////
// A set of plotting utilities for the
// electron chapter of the PPR
//
// J.L.Klay (Cal Poly)
// 28-Oct-2009
////////////////////////////////////////

TH1F *alltte,    *alltrk,    *allemc;
TH1F *sumtte,    *sumtrk,    *sumemc;
TH1F *sumHFemc,  *sumHFmc;
TH1F *btte,      *btrk,      *bptemc;
TH1F *ctte,      *ctrk,      *cptemc;
TH1F *cbtte,     *cbtrk,     *cbemc;
TH1F *convtte,   *convtrk,   *convemc;
TH1F *daltte,    *daltrk,    *dalemc;
TH1F *wztte,     *wzpttrk,   *wzemc;
TH1F *othtte,    *othtrk,    *othemc;
TH1F *htte,      *htrk,      *hemc;

TH1F *allmc;
TH1F *sigemc, *bkgemc, *wallemc, *hijemc;
TH1F* belemc, *celemc, *candbmc;
TH1F *convmc, *dalmc, *wzmc, *othermc;
TH1F* mchad;

void makeData(char* hijfname = "data/scaled25Oct09/histosLHC08d6.root",
              char* jjfname = "data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
              char* bfname = "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
              char* wfname = "data/scaled25Oct09/histosWboson.root") {


  //For HIJING need to divide by the number of events, which we
  //can get from the file and do when we perform scaling
  double hijscale = 0.05*(1.E6)*0.5*7700; //0-5% * seconds*lumi*PbPb
					  //x-section                      
  //For bjet and jet-jet events
  double pyscale = (1.E6)*0.5*208*208*100/360; //seconds*lumi*Pb*Pb*acceptance
  double bscale = pyscale; //Do we need to scale by Branching ratio
			   //for forced                      
  //semi-leptonic decays?                                             
  double wscale = pyscale;

  TFile* hijfile = new TFile(hijfname);
  if(!hijfile) { printf("NO HIJING FILE\n"); return; }
  TList* hijlist = (TList*)hijfile->Get("histos");
  TH2F* hijtte = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTTE");
  TH2F* hijemc2d = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleEMCAL");
  TH2F* hijtrk = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTPCTRD");
  TH1F* hijmult = (TH1F*)histos->FindObject("AnaElectron_hRefMult");
  Int_t nEvt = hijmult->GetEntries();
  if(nEvt == 0) { printf("NO HIJING EVENTS\n"); return; }
  hijtte->Scale(hijscale/nEvt);
  hijemc2d->Scale(hijscale/nEvt);
  hijtrk->Scale(hijscale/nEvt);
  TH2F* hijmcele = (TH2F*)histos->FindObject("AnaElectron_hPtMCElectron");
  TH1F* hijmchad = (TH1F*)histos->FindObject("AnaElectron_hPtMCHadron");
  hijmcele->Scale(hijscale/nEvt);
  hijmchad->Scale(hijscale/nEvt);

  TFile* jjfile = new TFile(jjfname);
  if(!jjfile) { printf("NO JET-JET FILE\n"); return; }
  TH2F* jjtte = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* jjemc = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleEMCALScaled");
  TH2F* jjtrk = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleTPCTRDScaled");
  TH1F* jjmult = (TH1F*)jjfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtJJ = jjmult->GetEntries();
  jjtte->Scale(pyscale);
  jjemc->Scale(pyscale);
  jjtrk->Scale(pyscale);
  TH2F* jjmcele = (TH2F*)jjfile->Get("AnaElectron_hPtMCElectronScaled");
  TH1F* jjmchad = (TH1F*)jjfile->Get("AnaElectron_hPtMCHadronScaled");
  jjmcele->Scale(pyscale);
  jjmchad->Scale(pyscale);

  TFile* bfile = new TFile(bfname);
  if(!bfile) { printf("NO B-JET FILE\n"); return; }
  TH2F* bjtte = (TH2F*)bfile->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* bjemc = (TH2F*)bfile->Get("AnaElectron_hPtNPEleEMCALScaled");
  TH2F* bjtrk = (TH2F*)bfile->Get("AnaElectron_hPtNPEleTPCTRDScaled");
  TH1F* bjmult = (TH1F*)bfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtB = bjmult->GetEntries();
  bjtte->Scale(bscale);
  bjemc->Scale(bscale);
  bjtrk->Scale(bscale);
  TH2F* bjmcele = (TH2F*)bfile->Get("AnaElectron_hPtMCElectronScaled");
  TH1F* bjmchad = (TH1F*)bfile->Get("AnaElectron_hPtMCHadronScaled");
  bjmcele->Scale(bscale);
  bjmchad->Scale(bscale);

  TFile* wfile = new TFile(wfname);
  if(!wfile) { printf("NO W-BOSON FILE\n"); return; }
  TH2F* wjtte = (TH2F*)wfile->Get("AnaElectron_hPtNPEleTTE");
  TH2F* wjemc = (TH2F*)wfile->Get("AnaElectron_hPtNPEleEMCAL");
  TH2F* wjtrk = (TH2F*)wfile->Get("AnaElectron_hPtNPEleTPCTRD");
  TH1F* wjmult = (TH1F*)wfile->Get("AnaElectron_hRefMult");
  Int_t nEvtW = wjmult->GetEntries();
  wjtte->Scale(wscale); //already scaled by evts
  wjemc->Scale(wscale); //already scaled by evts
  wjtrk->Scale(wscale); //already scaled by evts
  TH2F* wjmcele = (TH2F*)wfile->Get("AnaElectron_hPtMCElectron");
  TH1F* wjmchad = (TH1F*)wfile->Get("AnaElectron_hPtMCHadron");
  wjmcele->Scale(wscale);
  wjmchad->Scale(wscale);

  printf("Event statistics: %d (HIJING)  %d (JET-JET)  %d (B-JET)  %d (W-Boson)\n",nEvt,nEvtJJ,nEvtB,nEvtW);

  TH2F* combined = (TH2F*)hijmcele->Clone();
  combined->Add(jjmcele);
  combined->Add(bjmcele);  
  combined->Add(wjmcele);
  combined->SetTitle("MC electrons in Pb+Pb in EMCAL acceptance");
  combined->SetName("CombinedMCEle");
  combined->SetXTitle("p_T (GeV/c)");

  mchad = (TH1F*)hijmchad->Clone();
  mchad->Add(jjmchad);
  mchad->Add(bjmchad);
  mchad->Add(wjmchad);
  mchad->SetTitle("MC hadrons in Pb+Pb in EMCAL acceptance");
  mchad->SetName("CombinedMCHad");
  mchad->SetXTitle("p_T (GeV/c)");

  allmc = (TH1F*)combined->ProjectionX("allmc",1,1);
  belemc = (TH1F*)combined->ProjectionX("bmc",2,2);
  celemc = (TH1F*)combined->ProjectionX("cmc",3,3);
  candbmc = (TH1F*)combined->ProjectionX("candbmc",4,4);
  convmc = (TH1F*)combined->ProjectionX("convmc",5,5);
  dalmc = (TH1F*)combined->ProjectionX("dalmc",6,6);
  wzmc = (TH1F*)combined->ProjectionX("wzmc",7,7);
  othermc = (TH1F*)combined->ProjectionX("othermc",8,8);

  //For comparing contributions                        
  wallemc = (TH1F*)wjmcele->ProjectionX("wallemc",7,7);
  sigemc = (TH1F*)bjmcele->ProjectionX("sigemc",1,1);
  bkgemc = (TH1F*)jjmcele->ProjectionX("bkgemc",1,1);
  hijemc = (TH1F*)hijmcele->ProjectionX("hijemc",1,1);

  double myscale = 1.; //we already scaled them       
  ScaleAndConfigure(allmc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(belemc,myscale,kRed,kFALSE);
  ScaleAndConfigure(sigemc,myscale,kRed,kFALSE);
  ScaleAndConfigure(celemc,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candbmc,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convmc,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(bkgemc,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalmc,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(wzmc,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(wallemc,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(mchad,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(hijemc,myscale,kGreen+2,kFALSE);

  TH2F* combTTE = (TH2F*)hijtte->Clone();
  combTTE->Add(jjtte);
  combTTE->Add(bjtte);
  combTTE->Add(wjtte);
  combTTE->SetTitle("Identified non-phot. electrons (TPC+TRD+EMCAL)");
  combTTE->SetName("CombinedEleTTE");
  combTTE->SetXTitle("p_T (GeV/c)");

  alltte = (TH1F*)combTTE->ProjectionX("alltte",1,1);
  btte = (TH1F*)combTTE->ProjectionX("btte",2,2);
  ctte = (TH1F*)combTTE->ProjectionX("ctte",3,3);
  cbtte = (TH1F*)combTTE->ProjectionX("cbtte",4,4);
  convtte = (TH1F*)combTTE->ProjectionX("convtte",5,5);
  daltte = (TH1F*)combTTE->ProjectionX("daltte",6,6);
  wztte = (TH1F*)combTTE->ProjectionX("wztte",7,7);
  othtte = (TH1F*)combTTE->ProjectionX("othtte",8,8);
  sumtte = (TH1F*)btte->Clone(); sumtte->SetName("sumtte");
  sumtte->Add(ctte); sumtte->Add(cbtte); sumtte->Add(convtte);
  sumtte->Add(daltte); sumtte->Add(wztte); //sumtte->Add(othtte);  
  htte = (TH1F*)combTTE->ProjectionX("htte",9,9);

  ScaleAndConfigure(alltte,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumtte,myscale,kBlack,kFALSE);
  ScaleAndConfigure(btte,myscale,kRed,kFALSE);
  ScaleAndConfigure(ctte,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbtte,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convtte,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(daltte,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(htte,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wztte,myscale,kOrange-7,kFALSE);

  TH2F* combEMC = (TH2F*)hijemc2d->Clone();
  combEMC->Add(jjemc);
  combEMC->Add(bjemc);
  combEMC->Add(wjemc);
  combEMC->SetTitle("Identified non-phot. electrons (EMCAL)");
  combEMC->SetName("CombinedEleEMC");
  combEMC->SetXTitle("p_T (GeV/c)");

  allemc = (TH1F*)combEMC->ProjectionX("allemc",1,1);
  bemc = (TH1F*)combEMC->ProjectionX("bemc",2,2);
  cemc = (TH1F*)combEMC->ProjectionX("cemc",3,3);
  cbemc = (TH1F*)combEMC->ProjectionX("cbemc",4,4);
  convemc = (TH1F*)combEMC->ProjectionX("convemc",5,5);
  dalemc = (TH1F*)combEMC->ProjectionX("dalemc",6,6);
  wzemc = (TH1F*)combEMC->ProjectionX("wzemc",7,7);
  othemc = (TH1F*)combEMC->ProjectionX("othemc",8,8);
  hemc = (TH1F*)combEMC->ProjectionX("hemc",9,9);

  sumemc = (TH1F*)bemc->Clone(); sumemc->SetName("sumemc");
  sumemc->Add(cemc); sumemc->Add(cbemc); sumemc->Add(convemc);
  sumemc->Add(dalemc); //sumemc->Add(othemc);
  sumemc->Add(wzemc);
  sumHFemc = (TH1F*)bemc->Clone(); sumHFemc->SetName("sumHFemc");
  sumHFemc->Add(cemc); sumHFemc->Add(cbemc); sumHFemc->Add(wzemc);

  ScaleAndConfigure(allemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumHFemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bemc,myscale,kRed,kFALSE);
  ScaleAndConfigure(cemc,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbemc,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convemc,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalemc,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(hemc,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wzemc,myscale,kOrange-7,kFALSE);

  TH2F* combTRK = (TH2F*)hijtrk->Clone();
  combTRK->Add(jjtrk);
  combTRK->Add(bjtrk);
  combTRK->Add(wjtrk);
  combTRK->SetTitle("Identified non-phot. electrons (TPC+TRD)");
  combTRK->SetName("CombinedEleTRK");
  combTRK->SetXTitle("p_T (GeV/c)");

  alltrk = (TH1F*)combTRK->ProjectionX("alltrk",1,1);
  btrk = (TH1F*)combTRK->ProjectionX("btrk",2,2);
  ctrk = (TH1F*)combTRK->ProjectionX("ctrk",3,3);
  cbtrk = (TH1F*)combTRK->ProjectionX("cbtrk",4,4);
  convtrk = (TH1F*)combTRK->ProjectionX("convtrk",5,5);
  daltrk = (TH1F*)combTRK->ProjectionX("daltrk",6,6);
  wztrk = (TH1F*)combTRK->ProjectionX("wztrk",7,7);
  othtrk = (TH1F*)combTRK->ProjectionX("othtrk",8,8);
  sumtrk = (TH1F*)btrk->Clone(); sumtrk->SetName("sumtrk");
  sumtrk->Add(ctrk); sumtrk->Add(cbtrk); sumtrk->Add(convtrk);
  sumtrk->Add(daltrk); //sumtrk->Add(othtrk);  sumtrk->Add(wztrk);
  htrk = (TH1F*)combTRK->ProjectionX("htrk",9,9);

  for(Int_t i = 1; i <= alltrk->GetNbinsX(); i++) {
    Double_t myall = alltrk->GetBinContent(i);
    Double_t partsum = btrk->GetBinContent(i) +
      ctrk->GetBinContent(i) + cbtrk->GetBinContent(i) +
      convtrk->GetBinContent(i) + daltrk->GetBinContent(i) +
      wztrk->GetBinContent(i) + othtrk->GetBinContent(i);
    Double_t mysum = partsum + htrk->GetBinContent(i);
    printf("<%d> Compare bins All: %d, Sum: %d   Had: %d partSum: %d\n",i,myall,mysum,htrk->GetBinContent(i),partsum);
  }

  ScaleAndConfigure(alltrk,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumtrk,myscale,kBlack,kFALSE);
  ScaleAndConfigure(btrk,myscale,kRed,kFALSE);
  ScaleAndConfigure(ctrk,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbtrk,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convtrk,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(daltrk,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(htrk,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wztrk,myscale,kOrange-7,kFALSE);

  TFile* fout = new TFile("CombinedCocktailHistograms.root","RECREATE");
  fout->cd();
  combined->Write();
  combTTE->Write();
  combTRK->Write();
  combEMC->Write();
  alltte->Write();
  alltrk->Write();
  allemc->Write();
  sumtte->Write();
  sumtrk->Write();
  sumemc->Write();
  btte->Write();
  btrk->Write();
  bemc->Write();
  ctte->Write();
  ctrk->Write();
  cemc->Write();
  cbtte->Write();
  cbtrk->Write();
  cbemc->Write();
  convtte->Write();
  convtrk->Write();
  convemc->Write();
  daltte->Write();
  daltrk->Write();
  dalemc->Write();
  othtte->Write();
  othtrk->Write();
  othemc->Write();
  wztte->Write();
  wztrk->Write();
  wzemc->Write();
  htte->Write();
  htrk->Write();
  hemc->Write();
  allmc->Write();
  belemc->Write();
  celemc->Write();
  candbmc->Write();
  convmc->Write();
  dalmc->Write();
  wzmc->Write();
  othermc->Write();
  mchad->Write();
  sigemc->Write();
  bkgemc->Write();
  wallemc->Write();
  hijemc->Write();
  fout->Close();

}


void ScaleAndConfigure(TH1F* hist,Double_t scale, Int_t color,Bool_t keepErr)
{
  hist->Scale(scale);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  if(keepErr == kFALSE) {
    //remove the error bars - useful for MC rates
    for(Int_t i = 1; i <= hist->GetNbinsX(); i++) {
      if(hist->GetBinContent(i) > 0.) {
        if(hist->GetBinError(i)/hist->GetBinContent(i) > 0.5) {
          Double_t avg = 0.;
          if(i > 1 && i < hist->GetNbinsX())
            avg = (hist->GetBinContent(i-1) + hist->GetBinContent(i+1))/2.;
          hist->SetBinContent(i,avg);
        }
      }
      hist->SetBinError(i,0.);
    }
  } else {
    //Set the error bars to statistics of the bin content
    for(Int_t i = 1; i <= hist->GetNbinsX(); i++) {
      if(hist->GetBinContent(i) > 0.) {
        if(hist->GetBinError(i)/hist->GetBinContent(i) > 0.5) {
          Double_t avg = 0;
          if(i > 1 && i < hist->GetNbinsX())
            avg = (hist->GetBinContent(i-1) + hist->GetBinContent(i+1))/2.;
          hist->SetBinContent(i,avg);
        }
      }
      hist->SetBinError(i,TMath::Sqrt(hist->GetBinContent(i)));
    }
  }
}

TH1F* GetPtCutHisto(TH1F* input)
{
  //Given a rate histogram vs pt, return the histogram with yield
  //above a given pTcut                                          

  TH1F* result = (TH1F*)input->Clone();
  char name[100];
  sprintf(name,"%s_ptCut",result->GetName());
  result->SetNameTitle(name,name);
  for(Int_t i = 1; i <= result->GetNbinsX(); i++) {
    Double_t val = input->Integral(i,result->GetNbinsX());
    result->SetBinContent(i,val);
    result->SetBinError(i,0.);
  }

  return result;

}
