/////////////////////////////////////////
// A set of plotting utilities for the
// electron chapter of the PPR
//
// J.L.Klay (Cal Poly)
// 28-Oct-2009
////////////////////////////////////////

TH1F *alltte,    *alltrk,    *allemc;
TH1F *sumtte,    *sumtrk,    *sumemc;  //all but misid'ed hadrons
TH1F *btte,      *btrk,      *bemc;
TH1F *ctte,      *ctrk,      *cemc;
TH1F *cbtte,     *cbtrk,     *cbemc;
TH1F *convtte,   *convtrk,   *convemc;
TH1F *daltte,    *daltrk,    *dalemc;
TH1F *wztte,     *wztrk,     *wzemc;
TH1F *htte,      *htrk,      *hemc;

TH1F *allMC;
TH1F *bMC, *cMC, *cbMC;
TH1F *convMC, *dalMC, *wzMC;
TH1F *mchad;
TH1F *allheratio, *behratio;

TF1* fpow;

void makeData(char* jjfname = "data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
              char* bfname = "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
              char* wfname = "data/scaled25Oct09/histosWboson.root") {

  //TO AVOID DOUBLE-COUNTING ELECTRONS:
  //NOTE:  Jet-Jet events are "minimum bias" which means that there
  //are B-Jets included and for the higher pThard bins the B-jets will
  //play a more important role.  So to AVOID DOUBLE-COUNTING
  //B-electrons, we will suppress the electrons in the Jet-Jet events
  //that come from B->e or B->C->e when we combine the productions
  //into the total histograms

  //TO AVOID DOUBLE-COUNTING HADRONS:
  //For the hadron yields, we need to use only the yields from the
  //Jet-Jet events + W events and NOT USE the hadrons from the B-Jet events

  //For bjet and jet-jet events
  double pyscale = (1.E6)*0.5*208*208*100/360; //seconds*lumi*Pb*Pb*acceptance
  double bscale = pyscale; //Do we need to scale by Branching ratio
			   //for forced semi-leptonic decays?
  double wscale = pyscale;

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

  printf("Event statistics: %d (JET-JET)  %d (B-JET)  %d (W-Boson)\n",nEvtJJ,nEvtB,nEvtW);

  makeMCElectrons(jjmcele,bjmcele,wjmcele,jjmchad,wjmchad);
  makeTTEElectrons(jjtte,bjtte,wjtte);
  makeEMCElectrons(jjemc,bjemc,wjemc);
  makeTRKElectrons(jjtrk,bjtrk,wjtrk);

  TFile* fout = new TFile("CombinedCocktailHistograms.root","RECREATE");
  fout->cd();
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

  wztte->Write();
  wztrk->Write();
  wzemc->Write();

  htte->Write();
  htrk->Write();
  hemc->Write();

  allMC->Write();
  bMC->Write();
  cMC->Write();
  cbMC->Write();
  convMC->Write();
  dalMC->Write();
  wzMC->Write();
  mchad->Write();

  allheratio->Write();
  behratio->Write();

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
  //  printf("####### %s #########\n",name);
  result->SetNameTitle(name,name);
  for(Int_t i = 1; i <= result->GetNbinsX(); i++) {
    Double_t val = input->Integral(i,result->GetNbinsX());
    //printf("<Bin %d> value %.2f integral above %.2f\n",i,input->GetBinContent(i),val);
    result->SetBinContent(i,val);
    result->SetBinError(i,0.);
  }
  //printf("################\n");
  return result;

}

void makeMCElectrons(TH2F* jjmcele, TH2F* bjmcele, TH2F* wjmcele, TH1F* jjmchad, TH1F* wjmchad) {

  //Jet-Jet MC electrons
  TH1F* jjallmc = (TH1F*)jjmcele->ProjectionX("jjallmc",1,1);
  TH1F* jjbmc = (TH1F*)jjmcele->ProjectionX("jjbmc",2,2);
  TH1F* jjcmc = (TH1F*)jjmcele->ProjectionX("jjcmc",3,3);
  TH1F* jjcandbmc = (TH1F*)jjmcele->ProjectionX("jjcandbmc",4,4);
  TH1F* jjconvmc = (TH1F*)jjmcele->ProjectionX("jjconvmc",5,5);
  TH1F* jjdalmc = (TH1F*)jjmcele->ProjectionX("jjdalmc",6,6);

  //Bottom-Jet MC electrons
  TH1F* bjallmc = (TH1F*)bjmcele->ProjectionX("bjallmc",1,1);
  TH1F* bjbmc = (TH1F*)bjmcele->ProjectionX("bjbmc",2,2);
  TH1F* bjcmc = (TH1F*)bjmcele->ProjectionX("bjcmc",3,3);
  TH1F* bjcandbmc = (TH1F*)bjmcele->ProjectionX("bjcandbmc",4,4);
  TH1F* bjconvmc = (TH1F*)bjmcele->ProjectionX("bjconvmc",5,5);
  TH1F* bjdalmc = (TH1F*)bjmcele->ProjectionX("bjdalmc",6,6);

  //W-Jet MC electrons
  TH1F* wjallmc = (TH1F*)wjmcele->ProjectionX("wjallmc",1,1);
  TH1F* wjbmc = (TH1F*)wjmcele->ProjectionX("wjbmc",2,2);
  TH1F* wjcmc = (TH1F*)wjmcele->ProjectionX("wjcmc",3,3);
  TH1F* wjcandbmc = (TH1F*)wjmcele->ProjectionX("wjcandbmc",4,4);
  TH1F* wjconvmc = (TH1F*)wjmcele->ProjectionX("wjconvmc",5,5);
  TH1F* wjdalmc = (TH1F*)wjmcele->ProjectionX("wjdalmc",6,6);
  TH1F* wjwzmc = (TH1F*)wjmcele->ProjectionX("wjwzmc",7,7);

  //MC Hadrons (from jj events only)
  TCanvas *ctemp = new TCanvas("ctemp");
  ctemp->Divide(2,3);
  ctemp->cd(1); gPad->SetLogy(); jjmchad->Draw();
  mchad = (TH1F*)jjmchad->Clone(); mchad->SetName("mchad");
  TH1F* temphad = (TH1F*)mchad->Clone(); temphad->SetName("temphad");
  smoothWithFit(temphad,1e5,-3,10,60,100);  
  for(Int_t i = 10; i<= mchad->GetNbinsX(); i++) {
    Double_t pt = mchad->GetBinCenter(i);
    mchad->SetBinContent(i,temphad->GetFunction("fpow")->Eval(pt));
  }
  jjmchad->Draw();
  temphad->Draw("same");
  mchad->Draw("same");
  mchad->SetTitle("MC hadrons in Pb+Pb in EMCAL acceptance");
  mchad->SetName("mchad");
  mchad->SetXTitle("p_T (GeV/c)");

  bMC = (TH1F*)bjbmc->Clone(); bMC->SetName("bMC");  //B-Jet + W-jet
  bMC->Add(wjbmc);
  ctemp->cd(2); gPad->SetLogy(); bMC->Draw();
  TH1F* foob = (TH1F*)bMC->Clone(); foob->SetName("foob");
  smoothWithFit(foob,1e6,-3,15,50,100);  
  for(Int_t i = 10; i<= bMC->GetNbinsX(); i++) bMC->SetBinContent(i,foob->GetBinContent(i));
  bMC->Draw();
  foob->Draw("same");

  cMC = (TH1F*)jjcmc->Clone(); cMC->SetName("cMC"); //Jet-Jet + W-jet
  cMC->Add(wjcmc);
  ctemp->cd(3); gPad->SetLogy(); cMC->Draw();
  TH1F* fooc = (TH1F*)cMC->Clone(); fooc->SetName("fooc");
  smoothWithFit(fooc,1e6,-3,5,30,100);  
  for(Int_t i = 5; i<= cMC->GetNbinsX(); i++) cMC->SetBinContent(i,fooc->GetBinContent(i));
  cMC->Draw();
  fooc->Draw("same");

  cbMC = (TH1F*)bjcandbmc->Clone(); cbMC->SetName("cbMC"); //B-Jet + W-jet
  cbMC->Add(wjcandbmc);
  ctemp->cd(4); gPad->SetLogy(); cbMC->Draw();
  TH1F* foocb = (TH1F*)cbMC->Clone(); foocb->SetName("foocb");
  smoothWithFit(foocb,1e6,-3,8,35,100);  
  for(Int_t i = 8; i<= cbMC->GetNbinsX(); i++) cbMC->SetBinContent(i,foocb->GetBinContent(i));
  cbMC->Draw();
  foocb->Draw("same");

  convMC = (TH1F*)jjconvmc->Clone(); convMC->SetName("convMC"); //Jet-Jet + W-jet
  convMC->Add(wjconvmc);
  ctemp->cd(5); gPad->SetLogy(); convMC->Draw();
  TH1F* fooconv = (TH1F*)convMC->Clone(); fooconv->SetName("fooconv");
  smoothWithFit(fooconv,1e6,-3,6,40,100);  
  //  for(Int_t i = 6; i<= convMC->GetNbinsX(); i++) convMC->SetBinContent(i,fooconv->GetBinContent(i));
  convMC->Draw();
  fooconv->Draw("same");

  dalMC = (TH1F*)jjdalmc->Clone(); dalMC->SetName("dalMC"); //Jet-Jet + W-jet
  dalMC->Add(wjdalmc);
  ctemp->cd(6); gPad->SetLogy(); dalMC->Draw();
  TH1F* foodal = (TH1F*)dalMC->Clone(); foodal->SetName("foodal");
  smoothWithFit(foodal,1e6,-3,3,20,100);  
  for(Int_t i = 3; i<= dalMC->GetNbinsX(); i++) dalMC->SetBinContent(i,foodal->GetBinContent(i));
  dalMC->Draw();
  foodal->Draw("same");

  wzMC = (TH1F*)wjwzmc->Clone(); wzMC->SetName("wzMC"); //W-jet only
  TCanvas* cw= new TCanvas("cw");
  cw->cd(); gPad->SetLogy(); wzMC->GetYaxis()->SetRangeUser(1,1e3); wzMC->Draw();
  TH1F* foowz = (TH1F*)wzMC->Clone(); foowz->SetName("foowz");
  TF1* fws = new TF1("fws","[0]*(1+exp((x-[1])/[2]))^-1",39,100);
  fws->SetParameters(100,30,5);
  foowz->Fit(fws,"R");
  TF1* fwzexp = new TF1("fwzexp","[0]+[1]*log(x/[2])^2",4,40);
  fwzexp->SetParameters(100,10,3);
  foowz->Fit(fwzexp,"R");
  for(Int_t i = 8; i<= wzMC->GetNbinsX(); i++) {
    Double_t pt = wzMC->GetBinCenter(i);
    if(pt < 40) wzMC->SetBinContent(i,fwzexp->Eval(pt));
    if(pt > 40) wzMC->SetBinContent(i,fws->Eval(pt));
  }
  wzMC->GetYaxis()->SetRangeUser(1,1e3);
  wzMC->Draw();
  fws->Draw("same");
  foowz->Draw("same");

  //All mc electrons is the sum of 
  //Jet-Jet: conversions + direct charm + dalitz
  //Bottom-Jet: direct bottom + indirect bottom
  //W-Jet: all (because these events are exclusive of the others)
  allMC = (TH1F*)wzMC->Clone(); allMC->SetNameTitle("allMC","All MC Electrons");
  allMC->Add(convMC); allMC->Add(cMC); allMC->Add(dalMC);
  allMC->Add(bMC); allMC->Add(cbMC);

  //Hadron/electron ratios
  allheratio = (TH1F*)allMC->Clone(); allheratio->SetName("allheratio");
  behratio = (TH1F*)bMC->Clone(); behratio->SetName("behratio");
  allheratio->SetTitle("MC hadrons and electrons in Pb+Pb, 5.5 TeV");
  allheratio->SetXTitle("p_{T} (GeV/c)");
  allheratio->SetYTitle("Hadrons/Electrons");
  for(Int_t i = 1; i <= allheratio->GetNbinsX(); i++) {
    Double_t vale = allMC->GetBinContent(i);
    Double_t valb = bMC->GetBinContent(i);
    Double_t valh = mchad->GetBinContent(i);
    //printf("pT %.2f, Hadron %.1f, Electron %.1f, B-electron
    //%.1f\n",all->GetBinCenter(i),valh,vale,valb);             
    if(vale>0) allheratio->SetBinContent(i,valh/vale);
    else allheratio->SetBinContent(i,0.);
    
    if(valb>0) behratio->SetBinContent(i,valh/valb);
    else behratio->SetBinContent(i,0.);
    
    allheratio->SetBinError(i,0.);
    behratio->SetBinError(i,0.);
  }

  double myscale = 1.; //we already scaled them       
  ScaleAndConfigure(allMC,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bMC,myscale,kRed,kFALSE);
  ScaleAndConfigure(cMC,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbMC,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convMC,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalMC,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(wzMC,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(mchad,myscale,kGreen+2,kFALSE);

  return;
}

void makeTTEElectrons(TH2F* jjtte, TH2F* bjtte, TH2F* wjtte) {
  
  //Jet-Jet TTE Electrons
  TH1F* jjalltte = (TH1F*)jjtte->ProjectionX("jjalltte",1,1);
  TH1F* jjbtte = (TH1F*)jjtte->ProjectionX("jjbtte",2,2);
  TH1F* jjctte = (TH1F*)jjtte->ProjectionX("jjctte",3,3);
  TH1F* jjcbtte = (TH1F*)jjtte->ProjectionX("jjcbtte",4,4);
  TH1F* jjconvtte = (TH1F*)jjtte->ProjectionX("jjconvtte",5,5);
  TH1F* jjdaltte = (TH1F*)jjtte->ProjectionX("jjdaltte",6,6);
  TH1F* jjwztte = (TH1F*)jjtte->ProjectionX("jjwztte",7,7);
  TH1F* jjhtte = (TH1F*)jjtte->ProjectionX("jjhtte",9,9);

  //B-Jet TTE Electrons
  TH1F* bjalltte = (TH1F*)bjtte->ProjectionX("bjalltte",1,1);
  TH1F* bjbtte = (TH1F*)bjtte->ProjectionX("bjbtte",2,2);
  TH1F* bjctte = (TH1F*)bjtte->ProjectionX("bjctte",3,3);
  TH1F* bjcbtte = (TH1F*)bjtte->ProjectionX("bjcbtte",4,4);
  TH1F* bjconvtte = (TH1F*)bjtte->ProjectionX("bjconvtte",5,5);
  TH1F* bjdaltte = (TH1F*)bjtte->ProjectionX("bjdaltte",6,6);
  TH1F* bjwztte = (TH1F*)bjtte->ProjectionX("bjwztte",7,7);
  TH1F* bjhtte = (TH1F*)bjtte->ProjectionX("bjhtte",9,9);

  //W-Jet TTE Electrons
  TH1F* wjalltte = (TH1F*)wjtte->ProjectionX("wjalltte",1,1);
  TH1F* wjbtte = (TH1F*)wjtte->ProjectionX("wjbtte",2,2);
  TH1F* wjctte = (TH1F*)wjtte->ProjectionX("wjctte",3,3);
  TH1F* wjcbtte = (TH1F*)wjtte->ProjectionX("wjcbtte",4,4);
  TH1F* wjconvtte = (TH1F*)wjtte->ProjectionX("wjconvtte",5,5);
  TH1F* wjdaltte = (TH1F*)wjtte->ProjectionX("wjdaltte",6,6);
  TH1F* wjwztte = (TH1F*)wjtte->ProjectionX("wjwztte",7,7);
  TH1F* wjhtte = (TH1F*)wjtte->ProjectionX("wjhtte",9,9);
  
  btte = (TH1F*)bjbtte->Clone(); btte->SetName("btte");  //B-Jet + W-jet
  btte->Add(wjbtte)
  TCanvas * ctemptte = new TCanvas("ctemptte");
  ctemptte->Divide(2,3);
  ctemptte->cd(1); gPad->SetLogy(); btte->Draw();
  TH1F* foobtte = (TH1F*)btte->Clone(); foobtte->SetName("foobtte");
  smoothWithFit(foobtte,1e5,-3,8,40,100);
  btte->Draw(); foobtte->Draw("same");

  ctte = (TH1F*)jjctte->Clone(); ctte->SetName("ctte"); //Jet-Jet + W-jet
  ctte->Add(wjctte);
  ctemptte->cd(2); gPad->SetLogy(); ctte->Draw();
  TH1F* fooctte = (TH1F*)ctte->Clone(); fooctte->SetName("fooctte");
  smoothWithFit(fooctte,1e5,-3,3,12,100);
  ctte->Draw(); fooctte->Draw("same");

  cbtte = (TH1F*)bjcbtte->Clone(); cbtte->SetName("cbtte"); //B-Jet + W-jet
  cbtte->Add(wjcbtte);
  ctemptte->cd(3); gPad->SetLogy(); cbtte->Draw();
  TH1F* foocbtte = (TH1F*)cbtte->Clone(); foocbtte->SetName("foocbtte");
  smoothWithFit(foocbtte,1e5,-3,5,20,100);
  cbtte->Draw(); foocbtte->Draw("same");

  convtte = (TH1F*)jjconvtte->Clone(); convtte->SetName("convtte"); //Jet-Jet + W-jet
  convtte->Add(wjconvtte);
  ctemptte->cd(4); gPad->SetLogy(); convtte->Draw();
  TH1F* fooconvtte = (TH1F*)convtte->Clone(); fooconvtte->SetName("fooconvtte");
  smoothWithFit(fooconvtte,1e5,-3,6,13,100);
  convtte->Draw(); fooconvtte->Draw("same");

  daltte = (TH1F*)jjdaltte->Clone(); daltte->SetName("daltte"); //Jet-Jet + W-jet
  daltte->Add(wjdaltte);
  ctemptte->cd(5); gPad->SetLogy(); daltte->Draw();
  TH1F* foodaltte = (TH1F*)daltte->Clone(); foodaltte->SetName("foodaltte");
  smoothWithFit(foodaltte,1e2,-5,2,10,100);
  daltte->Draw(); foodaltte->Draw("same");

  htte = (TH1F*)jjhtte->Clone(); htte->SetName("htte");
  htte->Add(wjhtte);
  ctemptte->cd(6); gPad->SetLogy(); htte->Draw();
  TH1F* foohtte = (TH1F*)htte->Clone(); foohtte->SetName("foohtte");
  smoothWithFit(foohtte,1e5,-3,10,40,100);
  htte->Draw(); foohtte->Draw("same");

  wztte = (TH1F*)wjwztte->Clone(); wztte->SetName("wztte"); //W-jet only
  TCanvas* ctempwztte = new TCanvas("ctempwztte");
  ctempwztte->cd(); gPad->SetLogy(); wztte->Draw();
  TH1F* foowztte = (TH1F*)wztte->Clone(); foowztte->SetName("foowztte");
  smoothWithFit(foowztte,1e5,-3,10,40,100);
  wztte->Draw(); foowztte->Draw("same");

  //All TTE electrons is the sum of 
  //Jet-Jet: conversions + direct charm + dalitz + misid
  //Bottom-Jet: direct bottom + indirect bottom
  //W-Jet: all (because these events are exclusive of the others)
  alltte = (TH1F*)wjalltte->Clone(); alltte->SetName("alltte");
  alltte->Add(jjconvtte); alltte->Add(jjctte); alltte->Add(jjdaltte); alltte->Add(jjhtte);
  alltte->Add(bjbtte); alltte->Add(bjcbtte);
  sumtte = (TH1F*)wjalltte->Clone(); sumtte->SetName("sumtte");
  sumtte->Add(jjconvtte); sumtte->Add(jjctte); sumtte->Add(jjdaltte);
  sumtte->Add(bjbtte); sumtte->Add(bjcbtte);

  double myscale = 1.; //we already scaled them       
  ScaleAndConfigure(alltte,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumtte,myscale,kBlack,kFALSE);
  ScaleAndConfigure(btte,myscale,kRed,kFALSE);
  ScaleAndConfigure(ctte,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbtte,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convtte,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(daltte,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(htte,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wztte,myscale,kOrange-7,kFALSE);

  return;
}

void makeEMCElectrons(TH2F* jjemc, TH2F* bjemc, TH2F* wjemc) {
  
  //Jet-Jet EMC Electrons
  TH1F* jjallemc = (TH1F*)jjemc->ProjectionX("jjallemc",1,1);
  TH1F* jjbemc = (TH1F*)jjemc->ProjectionX("jjbemc",2,2);
  TH1F* jjcemc = (TH1F*)jjemc->ProjectionX("jjcemc",3,3);
  TH1F* jjcbemc = (TH1F*)jjemc->ProjectionX("jjcbemc",4,4);
  TH1F* jjconvemc = (TH1F*)jjemc->ProjectionX("jjconvemc",5,5);
  TH1F* jjdalemc = (TH1F*)jjemc->ProjectionX("jjdalemc",6,6);
  TH1F* jjwzemc = (TH1F*)jjemc->ProjectionX("jjwzemc",7,7);
  TH1F* jjhemc = (TH1F*)jjemc->ProjectionX("jjhemc",9,9);

  //B-Jet EMC Electrons
  TH1F* bjallemc = (TH1F*)bjemc->ProjectionX("bjallemc",1,1);
  TH1F* bjbemc = (TH1F*)bjemc->ProjectionX("bjbemc",2,2);
  TH1F* bjcemc = (TH1F*)bjemc->ProjectionX("bjcemc",3,3);
  TH1F* bjcbemc = (TH1F*)bjemc->ProjectionX("bjcbemc",4,4);
  TH1F* bjconvemc = (TH1F*)bjemc->ProjectionX("bjconvemc",5,5);
  TH1F* bjdalemc = (TH1F*)bjemc->ProjectionX("bjdalemc",6,6);
  TH1F* bjwzemc = (TH1F*)bjemc->ProjectionX("bjwzemc",7,7);
  TH1F* bjhemc = (TH1F*)bjemc->ProjectionX("bjhemc",9,9);

  //W-Jet EMC Electrons
  TH1F* wjallemc = (TH1F*)wjemc->ProjectionX("wjallemc",1,1);
  TH1F* wjbemc = (TH1F*)wjemc->ProjectionX("wjbemc",2,2);
  TH1F* wjcemc = (TH1F*)wjemc->ProjectionX("wjcemc",3,3);
  TH1F* wjcbemc = (TH1F*)wjemc->ProjectionX("wjcbemc",4,4);
  TH1F* wjconvemc = (TH1F*)wjemc->ProjectionX("wjconvemc",5,5);
  TH1F* wjdalemc = (TH1F*)wjemc->ProjectionX("wjdalemc",6,6);
  TH1F* wjwzemc = (TH1F*)wjemc->ProjectionX("wjwzemc",7,7);
  TH1F* wjhemc = (TH1F*)wjemc->ProjectionX("wjhemc",9,9);
  
  bemc = (TH1F*)bjbemc->Clone(); bemc->SetName("bemc");  //B-Jet + W-jet
  bemc->Add(wjbemc);
  TCanvas * ctempemc = new TCanvas("ctempemc","",0,0,800,800);
  ctempemc->Divide(2,3);
  ctempemc->cd(1); gPad->SetLogy(); bemc->Draw();
  TH1F* foobemc = (TH1F*)bemc->Clone(); foobemc->SetName("foobemc");
  smoothWithFit(foobemc,1e5,-3,8,40,100);
  TRandom rand;
  for(Int_t i = 8; i<= bemc->GetNbinsX(); i++) {
    Double_t dither = rand.Gaus(0.,foobemc->GetBinContent(i)/2.);
    if(dither + foobemc->GetBinContent(i) < foobemc->GetBinContent(i)/100.) dither = 0.;   
    bemc->SetBinContent(i,foobemc->GetBinContent(i)+dither);    
    bemc->SetBinError(i,sqrt(foobemc->GetBinContent(i)));    
  }
  bemc->Draw(); foobemc->Draw("same");

  cemc = (TH1F*)jjcemc->Clone(); cemc->SetName("cemc"); //Jet-Jet + W-jet
  cemc->Add(wjcemc);
  ctempemc->cd(2); gPad->SetLogy(); cemc->Draw();
  TH1F* foocemc = (TH1F*)cemc->Clone(); foocemc->SetName("foocemc");
  smoothWithFit(foocemc,1e5,-3,5,14,100);
  for(Int_t i = 5; i<= cemc->GetNbinsX(); i++) {
    Double_t dither = rand.Gaus(0.,foocemc->GetBinContent(i)/2.);
    if(dither + foocemc->GetBinContent(i) < foocemc->GetBinContent(i)/100.) dither = 0.;   
    cemc->SetBinContent(i,foocemc->GetBinContent(i)+dither);    
    cemc->SetBinError(i,sqrt(foocemc->GetBinContent(i)));    
  }
  cemc->Draw(); foocemc->Draw("same");

  cbemc = (TH1F*)bjcbemc->Clone(); cbemc->SetName("cbemc"); //B-Jet + W-jet
  cbemc->Add(wjcbemc);
  ctempemc->cd(3); gPad->SetLogy(); cbemc->Draw();
  TH1F* foocbemc = (TH1F*)cbemc->Clone(); foocbemc->SetName("foocbemc");
  smoothWithFit(foocbemc,1e5,-3,8,20,100);
  for(Int_t i = 8; i<= cbemc->GetNbinsX(); i++) {
    Double_t dither = rand.Gaus(0.,foocbemc->GetBinContent(i)/2.);
    if(dither + foocbemc->GetBinContent(i) < foocbemc->GetBinContent(i)/100.) dither = 0.;   
    cbemc->SetBinContent(i,foocbemc->GetBinContent(i)+dither);    
    cbemc->SetBinError(i,sqrt(foocbemc->GetBinContent(i)));    
  }
  cbemc->Draw(); foocbemc->Draw("same");

  convemc = (TH1F*)jjconvemc->Clone(); convemc->SetName("convemc"); //Jet-Jet + W-jet
  convemc->Add(wjconvemc);
  ctempemc->cd(4); gPad->SetLogy(); convemc->Draw();
  TH1F* fooconvemc = (TH1F*)convemc->Clone(); fooconvemc->SetName("fooconvemc");
  smoothWithFit(fooconvemc,1e6,-3,5,15,100);
  for(Int_t i = 5; i<= convemc->GetNbinsX(); i++) {
    Double_t dither = rand.Gaus(0.,fooconvemc->GetBinContent(i)/2.);
    if(dither + fooconvemc->GetBinContent(i) < fooconvemc->GetBinContent(i)/100.) dither = 0.;   
    convemc->SetBinContent(i,fooconvemc->GetBinContent(i)+dither);    
    convemc->SetBinError(i,sqrt(fooconvemc->GetBinContent(i)));    
  }
  convemc->Draw(); fooconvemc->Draw("same");

  dalemc = (TH1F*)jjdalemc->Clone(); dalemc->SetName("dalemc"); //Jet-Jet + W-jet
  dalemc->Add(wjdalemc);
  ctempemc->cd(5); gPad->SetLogy(); dalemc->Draw();
  TH1F* foodalemc = (TH1F*)dalemc->Clone(); foodalemc->SetName("foodalemc");
  for(Int_t i = 18; i <= dalemc->GetNbinsX(); i++) {
    dalemc->SetBinContent(i,0);
    dalemc->SetBinError(i,0);
    foodalemc->SetBinContent(i,0);
    foodalemc->SetBinError(i,0);
  }
  //  smoothWithFit(foodalemc,1e6,-3,8.,12.,100);
  dalemc->Draw(); //foodalemc->Draw("same");

  hemc = (TH1F*)jjhemc->Clone(); hemc->SetName("hemc");
  hemc->Add(wjhemc);
  ctempemc->cd(6); gPad->SetLogy(); hemc->Draw();
  TH1F* foohemc = (TH1F*)hemc->Clone(); foohemc->SetName("foohemc");
  smoothWithFit(foohemc,1e5,-3,10,40,100);
  for(Int_t i = 10; i<= hemc->GetNbinsX(); i++) {
    Double_t dither = rand.Gaus(0.,foohemc->GetBinContent(i)/2.);
    if(dither + foohemc->GetBinContent(i) < foohemc->GetBinContent(i)/100.) dither = 0.;   
    hemc->SetBinContent(i,foohemc->GetBinContent(i)+dither);    
    hemc->SetBinError(i,sqrt(foohemc->GetBinContent(i)));    
  }
  hemc->Draw(); foohemc->Draw("same");

  wzemc = (TH1F*)wjwzemc->Clone(); wzemc->SetName("wzemc"); //W-jet only
  TCanvas* ctempwzemc = new TCanvas("ctempwzemc");
  ctempwzemc->cd(); gPad->SetLogy(); wzemc->Draw();
  TH1F* foowzemc = (TH1F*)wzemc->Clone(); foowzemc->SetName("foowzemc");
  TF1* fwsemc = new TF1("fwsemc","[0]*(1+exp((x-[1])/[2]))^-1",30,50);
  fwsemc->SetParameters(10,30,3);
  foowzemc->Fit(fwsemc,"R");
  TF1* fwzexpemc = new TF1("fwzexpemc","[0]+[1]*log(x/[2])^2",5,20);
  fwzexpemc->SetParameters(10,10,3);
  //  foowzemc->Fit(fwzexpemc,"R");
  /*
  for(Int_t i = 8; i<= wzMC->GetNbinsX(); i++) {
    Double_t pt = wzMC->GetBinCenter(i);
    if(pt < 40) wzMC->SetBinContent(i,fwzexp->Eval(pt));
    if(pt > 40) wzMC->SetBinContent(i,fws->Eval(pt));
  }
  */
  wzemc->Draw(); foowzemc->Draw("same");
  fwzexpemc->Draw("same");

  //All EMC electrons is the sum of 
  //Jet-Jet: conversions + direct charm + dalitz + misid
  //Bottom-Jet: direct bottom + indirect bottom
  //W-Jet: all (because these events are exclusive of the others)
  allemc = (TH1F*)wjallemc->Clone(); allemc->SetName("allemc");
  allemc->Add(jjconvemc); allemc->Add(jjcemc); allemc->Add(jjdalemc); allemc->Add(jjhemc);
  allemc->Add(bjbemc); allemc->Add(bjcbemc);
  sumemc = (TH1F*)wjallemc->Clone(); sumemc->SetName("sumemc");
  sumemc->Add(jjconvemc); sumemc->Add(jjcemc); sumemc->Add(jjdalemc);
  sumemc->Add(bjbemc); sumemc->Add(bjcbemc);

  double myscale = 1.; //we already scaled them       
  ScaleAndConfigure(allemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bemc,myscale,kRed,kFALSE);
  ScaleAndConfigure(cemc,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbemc,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convemc,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalemc,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(hemc,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wzemc,myscale,kOrange-7,kFALSE);

  return;
}

void makeTRKElectrons(TH2F* jjtrk, TH2F* bjtrk, TH2F* wjtrk) {
  
  //Jet-Jet TRK Electrons
  TH1F* jjalltrk = (TH1F*)jjtrk->ProjectionX("jjalltrk",1,1);
  TH1F* jjbtrk = (TH1F*)jjtrk->ProjectionX("jjbtrk",2,2);
  TH1F* jjctrk = (TH1F*)jjtrk->ProjectionX("jjctrk",3,3);
  TH1F* jjcbtrk = (TH1F*)jjtrk->ProjectionX("jjcbtrk",4,4);
  TH1F* jjconvtrk = (TH1F*)jjtrk->ProjectionX("jjconvtrk",5,5);
  TH1F* jjdaltrk = (TH1F*)jjtrk->ProjectionX("jjdaltrk",6,6);
  TH1F* jjwztrk = (TH1F*)jjtrk->ProjectionX("jjwztrk",7,7);
  TH1F* jjhtrk = (TH1F*)jjtrk->ProjectionX("jjhtrk",9,9);

  //B-Jet TRK Electrons
  TH1F* bjalltrk = (TH1F*)bjtrk->ProjectionX("bjalltrk",1,1);
  TH1F* bjbtrk = (TH1F*)bjtrk->ProjectionX("bjbtrk",2,2);
  TH1F* bjctrk = (TH1F*)bjtrk->ProjectionX("bjctrk",3,3);
  TH1F* bjcbtrk = (TH1F*)bjtrk->ProjectionX("bjcbtrk",4,4);
  TH1F* bjconvtrk = (TH1F*)bjtrk->ProjectionX("bjconvtrk",5,5);
  TH1F* bjdaltrk = (TH1F*)bjtrk->ProjectionX("bjdaltrk",6,6);
  TH1F* bjwztrk = (TH1F*)bjtrk->ProjectionX("bjwztrk",7,7);
  TH1F* bjhtrk = (TH1F*)bjtrk->ProjectionX("bjhtrk",9,9);

  //W-Jet TRK Electrons
  TH1F* wjalltrk = (TH1F*)wjtrk->ProjectionX("wjalltrk",1,1);
  TH1F* wjbtrk = (TH1F*)wjtrk->ProjectionX("wjbtrk",2,2);
  TH1F* wjctrk = (TH1F*)wjtrk->ProjectionX("wjctrk",3,3);
  TH1F* wjcbtrk = (TH1F*)wjtrk->ProjectionX("wjcbtrk",4,4);
  TH1F* wjconvtrk = (TH1F*)wjtrk->ProjectionX("wjconvtrk",5,5);
  TH1F* wjdaltrk = (TH1F*)wjtrk->ProjectionX("wjdaltrk",6,6);
  TH1F* wjwztrk = (TH1F*)wjtrk->ProjectionX("wjwztrk",7,7);
  TH1F* wjhtrk = (TH1F*)wjtrk->ProjectionX("wjhtrk",9,9);
  
  btrk = (TH1F*)bjbtrk->Clone(); btrk->SetName("btrk");  //B-Jet + W-jet
  btrk->Add(wjbtrk);
  TCanvas * ctemptrk = new TCanvas("ctemptrk","",0,0,800,800);
  ctemptrk->Divide(2,3);
  ctemptrk->cd(1); gPad->SetLogy(); btrk->Draw();
  TH1F* foobtrk = (TH1F*)btrk->Clone(); foobtrk->SetName("foobtrk");
  smoothWithFit(foobtrk,1e5,-3,8,40,100);
  TRandom rand2;
  for(Int_t i = 8; i<= btrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,foobtrk->GetBinContent(i)/2.);
    if(dither + foobtrk->GetBinContent(i) < foobtrk->GetBinContent(i)/100.) dither = 0.;   
    btrk->SetBinContent(i,foobtrk->GetBinContent(i)+dither);    
    btrk->SetBinError(i,sqrt(foobtrk->GetBinContent(i)));    
  }
  btrk->Draw(); foobtrk->Draw("same");

  ctrk = (TH1F*)jjctrk->Clone(); ctrk->SetName("ctrk"); //Jet-Jet + W-jet
  ctrk->Add(wjctrk);
  ctemptrk->cd(2); gPad->SetLogy(); ctrk->Draw();
  TH1F* fooctrk = (TH1F*)ctrk->Clone(); fooctrk->SetName("fooctrk");
  smoothWithFit(fooctrk,1e5,-3,5,20,100);
  for(Int_t i = 5; i<= ctrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,fooctrk->GetBinContent(i)/2.);
    if(dither + fooctrk->GetBinContent(i) < fooctrk->GetBinContent(i)/100.) dither = 0.;   
    ctrk->SetBinContent(i,fooctrk->GetBinContent(i)+dither);    
    ctrk->SetBinError(i,sqrt(fooctrk->GetBinContent(i)));    
  }
  ctrk->Draw(); fooctrk->Draw("same");

  cbtrk = (TH1F*)bjcbtrk->Clone(); cbtrk->SetName("cbtrk"); //B-Jet + W-jet
  cbtrk->Add(wjcbtrk);
  ctemptrk->cd(3); gPad->SetLogy(); cbtrk->Draw();
  TH1F* foocbtrk = (TH1F*)cbtrk->Clone(); foocbtrk->SetName("foocbtrk");
  smoothWithFit(foocbtrk,1e5,-3,5,20,100);
  for(Int_t i = 5; i<= cbtrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,foocbtrk->GetBinContent(i)/2.);
    if(dither + foocbtrk->GetBinContent(i) < foocbtrk->GetBinContent(i)/100.) dither = 0.;   
    cbtrk->SetBinContent(i,foocbtrk->GetBinContent(i)+dither);    
    cbtrk->SetBinError(i,sqrt(foocbtrk->GetBinContent(i)));    
  }
  cbtrk->Draw(); foocbtrk->Draw("same");

  convtrk = (TH1F*)jjconvtrk->Clone(); convtrk->SetName("convtrk"); //Jet-Jet + W-jet
  convtrk->Add(wjconvtrk);
  ctemptrk->cd(4); gPad->SetLogy(); convtrk->Draw();
  TH1F* fooconvtrk = (TH1F*)convtrk->Clone(); fooconvtrk->SetName("fooconvtrk");
  smoothWithFit(fooconvtrk,1e6,-3,5,20,100);
  for(Int_t i = 5; i<= convtrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,fooconvtrk->GetBinContent(i)/2.);
    if(dither + fooconvtrk->GetBinContent(i) < fooconvtrk->GetBinContent(i)/100.) dither = 0.;   
    //    convtrk->SetBinContent(i,fooconvtrk->GetBinContent(i)+dither);    
    //convtrk->SetBinError(i,sqrt(fooconvtrk->GetBinContent(i)));    
  }
  convtrk->Draw(); fooconvtrk->Draw("same");

  daltrk = (TH1F*)jjdaltrk->Clone(); daltrk->SetName("daltrk"); //Jet-Jet + W-jet
  daltrk->Add(wjdaltrk);
  ctemptrk->cd(5); gPad->SetLogy(); daltrk->Draw();
  TH1F* foodaltrk = (TH1F*)daltrk->Clone(); foodaltrk->SetName("foodaltrk");
  smoothWithFit(foodaltrk,1e6,-3,5.,40.,100);
  for(Int_t i = 5; i<= daltrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,foodaltrk->GetBinContent(i)/2.);
    if(dither + foodaltrk->GetBinContent(i) < foodaltrk->GetBinContent(i)/100.) dither = 0.;   
    daltrk->SetBinContent(i,foodaltrk->GetBinContent(i)+dither);    
    daltrk->SetBinError(i,sqrt(foodaltrk->GetBinContent(i)));    
  }
  daltrk->Draw(); foodaltrk->Draw("same");

  htrk = (TH1F*)jjhtrk->Clone(); htrk->SetName("htrk");
  htrk->Add(wjhtrk);
  ctemptrk->cd(6); gPad->SetLogy(); htrk->Draw();
  TH1F* foohtrk = (TH1F*)htrk->Clone(); foohtrk->SetName("foohtrk");
  smoothWithFit(foohtrk,1e5,-3,5,80,100);
  for(Int_t i = 5; i<= htrk->GetNbinsX(); i++) {
    Double_t dither = rand2.Gaus(0.,foohtrk->GetBinContent(i)/2.);
    if(dither + foohtrk->GetBinContent(i) < foohtrk->GetBinContent(i)/100.) dither = 0.;   
    htrk->SetBinContent(i,foohtrk->GetBinContent(i)+dither);    
    htrk->SetBinError(i,sqrt(foohtrk->GetBinContent(i)));    
  }
  htrk->Draw(); foohtrk->Draw("same");

  wztrk = (TH1F*)wjwztrk->Clone(); wztrk->SetName("wztrk"); //W-jet only
  TCanvas* ctempwztrk = new TCanvas("ctempwztrk");
  ctempwztrk->cd(); gPad->SetLogy(); wztrk->Draw();
  TH1F* foowztrk = (TH1F*)wztrk->Clone(); foowztrk->SetName("foowztrk");
  TF1* fwstrk = new TF1("fwstrk","[0]*(1+exp((x-[1])/[2]))^-1",30,50);
  fwstrk->SetParameters(10,30,3);
  foowztrk->Fit(fwstrk,"R");
  TF1* fwzexptrk = new TF1("fwzexptrk","[0]+[1]*log(x/[2])^2",8,40);
  fwzexptrk->SetParameters(10,10,3);
  //  foowztrk->Fit(fwzexptrk,"R");
  /*
  for(Int_t i = 8; i<= wzMC->GetNbinsX(); i++) {
    Double_t pt = wzMC->GetBinCenter(i);
    if(pt < 40) wzMC->SetBinContent(i,fwzexp->Eval(pt));
    if(pt > 40) wzMC->SetBinContent(i,fws->Eval(pt));
  }
  */
  wztrk->Draw(); foowztrk->Draw("same");

  //All TRK electrons is the sum of 
  //Jet-Jet: conversions + direct charm + dalitz + misid
  //Bottom-Jet: direct bottom + indirect bottom
  //W-Jet: all (because these events are exclusive of the others)
  alltrk = (TH1F*)wjalltrk->Clone(); alltrk->SetName("alltrk");
  alltrk->Add(jjconvtrk); alltrk->Add(jjctrk); alltrk->Add(jjdaltrk); alltrk->Add(jjhtrk);
  alltrk->Add(bjbtrk); alltrk->Add(bjcbtrk);
  sumtrk = (TH1F*)wjalltrk->Clone(); sumtrk->SetName("sumtrk");
  sumtrk->Add(jjconvtrk); sumtrk->Add(jjctrk); sumtrk->Add(jjdaltrk);
  sumtrk->Add(bjbtrk); sumtrk->Add(bjcbtrk);

  double myscale = 1.; //we already scaled them       
  ScaleAndConfigure(alltrk,myscale,kBlack,kFALSE);
  ScaleAndConfigure(sumtrk,myscale,kBlack,kFALSE);
  ScaleAndConfigure(btrk,myscale,kRed,kFALSE);
  ScaleAndConfigure(ctrk,myscale,kBlue,kFALSE);
  ScaleAndConfigure(cbtrk,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convtrk,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(daltrk,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(htrk,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wztrk,myscale,kOrange-7,kFALSE);

  return;
}

void smoothWithFit(TH1F* hist, Double_t p0, Double_t p1, Double_t min, Double_t max,Double_t remax) {

  fpow = new TF1("fpow","[0]*pow(x,[1])",min,max);
  fpow->SetParameters(p0,p1);
  hist->Fit(fpow,"R");  
  for(Int_t i = (Int_t)min; i <= (Int_t)remax; i++) {
    Double_t pt = hist->GetBinCenter(i);
    Double_t val = fpow->Eval(pt);
    hist->SetBinContent(i,val);
    hist->SetBinError(i,sqrt(val));
  }

  return;
}


