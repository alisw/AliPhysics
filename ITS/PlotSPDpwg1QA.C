#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <TLine.h>
#include <TString.h>
#include <TPaveText.h>
#endif

void ratiomodules();
void ratiochips();
void mapsinner(Bool_t isShowMaps=kFALSE);
void phiTracklet();
void phiTrackletsZ();
void foEfficiency(); 


TFile *fData=0x0;
TFile *fMc=0x0;

TList *fListData = 0x0;
TList  *fListMc  = 0x0;

TString fTitleData = "";
TString fTitleMc = "";


void PlotSPDpwg1QA(TString data, TString mc, TString titleData = "[Data]", TString titleMc = "[MC]", Bool_t isGeneralTrain = kFALSE){

  fTitleData=titleData;
  fTitleMc=titleMc;

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

 fData = TFile::Open(data.Data());
  if(!isGeneralTrain) fListData = (TList*)fData->Get("chist");
  else {
    TDirectoryFile *spddata = (TDirectoryFile*)fData->Get("SPD_Performance");
    spddata->cd();
    fListData = (TList*)spddata->Get("coutput1");
  }

  Double_t nevtsData = ((TH1I*)(fListData->FindObject("hEventsProcessed")))->GetEntries();
  printf("   #events in %s : %f \n",fTitleData.Data(),nevtsData);

  fMc = TFile::Open(mc.Data());
  if(!isGeneralTrain) fListMc = (TList*)fMc->Get("chist");
  else {
    TDirectoryFile *spdmc = (TDirectoryFile*)fMc->Get("SPD_Performance");
    spdmc->cd();
    fListMc = (TList*)spdmc->Get("coutput1");
  }
  Double_t nevtsMc = ((TH1I*)(fListMc->FindObject("hEventsProcessed")))->GetEntries();
  // phi projection

  printf("   #events in %s : %f \n",fTitleMc.Data(),nevtsMc);
  printf("Available functions : \n - ratiomodules() \n - ratiochips() \n - mapsinner(isShowMaps) \n - phiTracklet() \n - phiTrackletsZ() \n - foEfficiency() \n");
  phiTracklet();
}

void ratiomodules(){

  TH1F *data_module = ((TH1F*)fListData->FindObject("hClusterModYield"));
  data_module->Sumw2();

  TH1F *mc_module = ((TH1F*)fListMc->FindObject("hClusterModYield"));
  mc_module->Sumw2();

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  data_module->SetTitle(Form("cluster yield in %s",fTitleData.Data()));
  data_module->DrawCopy();
  c1->cd(2);
  mc_module->SetTitle(Form("cluster yield in %s",fTitleMc.Data()));
  mc_module->DrawCopy();


  TH1D *ratiomodule  = new TH1D("ratiomod",Form("Module cluster ratio %s / %s  - scaling : Integral",fTitleData.Data(),fTitleMc.Data()),mc_module->GetNbinsX(),mc_module->GetXaxis()->GetXmin(),mc_module->GetXaxis()->GetXmax());
  ratiomodule->GetXaxis()->SetTitle("module number");
  printf("data_module %f     -  mc_module %f    \n",data_module->GetEntries(),mc_module->GetEntries());
  ratiomodule->Divide(data_module,mc_module,mc_module->Integral(),data_module->Integral(),"e");

  TCanvas *ratiofit = new TCanvas("ratiofit","ratiofit",1200,600);
  ratiofit->cd();
 
  TF1 *f1inner = new TF1("f1inner", "pol0", 0, 79);
  f1inner->SetLineColor(kBlue);
  TF1 *f1outer = new TF1("f1outer", "pol0", 80, 223);
  f1outer->SetLineColor(kGreen);
  
  printf("------------- Fitting standard normalized ratio  ----------------\n");
  ratiomodule->Fit("f1inner", "R+"); 
  ratiomodule->Fit("f1outer", "R+");
  ratiomodule->GetListOfFunctions()->Print();
   
  f1inner->Draw("same");
  f1outer->Draw("same");
  
  TPaveText *parsFitEntries = new TPaveText(0.3,0.23,0.6,0.55,"NDC"); 
  parsFitEntries->AddText("Inner ");
  parsFitEntries->AddText(Form("#chi2 / ndf    %4.3f / %2.1i",f1inner->GetChisquare(),f1inner->GetNDF()));
  parsFitEntries->AddText(Form("Fit :    %4.4f +- %4.4f",f1inner->GetParameter(0),f1inner->GetParError(0)));
  parsFitEntries->AddText("");
  parsFitEntries->AddText("Outer ");
  parsFitEntries->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",f1outer->GetChisquare(),f1outer->GetNDF()));
  parsFitEntries->AddText(Form("Fit :    %4.4f +- %4.4f",f1outer->GetParameter(0),f1outer->GetParError(0)));
  ratiomodule->GetListOfFunctions()->Add(parsFitEntries);




  TCanvas *ratiofitEvents = new TCanvas("ratiofitEvents","ratiofitEvents",1200,600);

  ratiofitEvents->cd();
  TH1D *ratiomoduleEvents  = new TH1D("ratiomodEvents",Form("Module cluster ratio %s / %s  - scaling : # events",fTitleData.Data(),fTitleMc.Data()),mc_module->GetNbinsX(),mc_module->GetXaxis()->GetXmin(),mc_module->GetXaxis()->GetXmax());
  ratiomoduleEvents->GetXaxis()->SetTitle("module number");
  Double_t nEvData = ((TH1I*)(fListData->FindObject("hEventsProcessed")))->GetEntries();
  Double_t nEvMc= ((TH1I*)(fListMc->FindObject("hEventsProcessed")))->GetEntries();
  ratiomoduleEvents->Divide(data_module,mc_module,nEvMc,nEvData,"e");
  ratiomoduleEvents->Draw();

  TF1 *fInner = new TF1("fInner", "pol0", 0, 79);
  fInner->SetLineColor(kBlue);
  TF1 *fOuter = new TF1("fOuter", "pol0", 80, 223);
  fOuter->SetLineColor(kGreen);
  printf("------------- Fitting #evts normalized ratio  ----------------\n");
  ratiomoduleEvents->Fit("fInner", "R");
  ratiomoduleEvents->Fit("fOuter", "R+");
  fInner->Draw("same");
  fOuter->Draw("same");
  
  
  
  TPaveText *parsFitEvents = new TPaveText(0.3,0.23,0.6,0.55,"NDC"); 
  parsFitEvents->AddText("Inner ");
  parsFitEvents->AddText(Form("#chi2 / ndf    %4.3f / %2.1i",fInner->GetChisquare(),fInner->GetNDF()));
  parsFitEvents->AddText(Form("Fit :    %4.4f +- %4.4f",fInner->GetParameter(0),fInner->GetParError(0)));
  parsFitEvents->AddText("");
  parsFitEvents->AddText("Outer ");
  parsFitEvents->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",fOuter->GetChisquare(),fOuter->GetNDF()));
  parsFitEvents->AddText(Form("Fit :    %4.4f +- %4.4f",fOuter->GetParameter(0),fOuter->GetParError(0)));
  ratiomoduleEvents->GetListOfFunctions()->Add(parsFitEvents);
  

}
void mapsinner(Bool_t isShowMaps){
  TH2D *data_mapL1 = (TH2D*)fListData->FindObject("hLocalMapL1");
  // ------------ phi projection ---------------------
  TH1D *data_projyL1 = data_mapL1->ProjectionY();

  TString titleDataL1  = data_projyL1->GetTitle();
  titleDataL1+=fTitleData.Data();
  data_projyL1->SetTitle(titleDataL1.Data());
  //data_projyL1->Rebin(10);
  data_projyL1->SetYTitle(Form("entries / %1.2f cm",data_projyL1->GetBinWidth(0)));

  TH2D *data_mapL2 = (TH2D*)fListData->FindObject("hLocalMapL2");
  TH1D *data_projyL2 = data_mapL2->ProjectionY();
  TString titleDataL2  = data_projyL2->GetTitle();
  titleDataL2+=fTitleData.Data();
  data_projyL2->SetTitle(titleDataL2.Data());
  //data_projyL2 ->Rebin(10);
  data_projyL2->SetYTitle(Form("entries / %1.2f cm",data_projyL2->GetBinWidth(0)));

  // ------- z projection Data -----------
  TH1D *data_projyL1z = data_mapL1->ProjectionX();
  TString titleDataL1z  = data_projyL1z->GetTitle();
  titleDataL1z+=fTitleData.Data();
  data_projyL1z->SetTitle(titleDataL1z.Data());
  //data_projyL1z->Rebin(10);
  data_projyL1z->SetYTitle(Form("entries / %1.2f cm",data_projyL1z->GetBinWidth(0)));
  TH1D *data_projyL2z = data_mapL2->ProjectionX();
  TString titleDataL2z  = data_projyL2z->GetTitle();
  titleDataL2z+=fTitleData.Data();
  data_projyL2z->SetTitle(titleDataL2.Data());
  //data_projyL2z ->Rebin(10);
  data_projyL2z->SetYTitle(Form("entries / %1.2f cm",data_projyL2z->GetBinWidth(0)));

  // ------------ phi projection MC ---------------------
  TH2D *mc_mapL1 = (TH2D*)fListMc->FindObject("hLocalMapL1");
  TH1D *mc_projyL1 = mc_mapL1->ProjectionY();
  TString titleMCL1  = mc_projyL1->GetTitle();
  titleMCL1+=fTitleMc.Data();
  mc_projyL1->SetTitle(titleMCL1.Data());
  //mc_projyL1->Rebin(10);
  mc_projyL1->SetYTitle(Form("entries / %1.2f cm",mc_projyL1->GetBinWidth(0)));

  TH2D *mc_mapL2 = (TH2D*)fListMc->FindObject("hLocalMapL2");
  TH1D *mc_projyL2 = mc_mapL2->ProjectionY();
  TString titleMCL2  = mc_projyL2->GetTitle();
  titleMCL2+=fTitleMc.Data();
  mc_projyL2->SetTitle(titleMCL2.Data());
  //mc_projyL2->Rebin(10);
  mc_projyL2->SetYTitle(Form("entries / %1.2f cm",mc_projyL2->GetBinWidth(0)));

  // ------- z projection MC -----------

  TH1D *mc_projyL1z = mc_mapL1->ProjectionX();
  TString titleMCL1z  = mc_projyL1z->GetTitle();
  titleMCL1z+=fTitleMc.Data();
  mc_projyL1z->SetTitle(titleMCL1z.Data());
  //mc_projyL1z->Rebin(10);
  mc_projyL1z->SetYTitle(Form("entries / %1.2f cm",mc_projyL1z->GetBinWidth(0)));

  TH1D *mc_projyL2z = mc_mapL2->ProjectionX();
  TString titleMCL2z  = mc_projyL2z->GetTitle();
  titleMCL2z+="  [ MC ]";
  mc_projyL2z->SetTitle(titleMCL2z.Data());
  //mc_projyL2z->Rebin(10);
  mc_projyL2z->SetYTitle(Form("entries / %1.2f cm",mc_projyL2z->GetBinWidth(0)));

  if(isShowMaps) {


    TCanvas *cmapData = new TCanvas("cmapData","cmapData",1200,600);
    cmapData->Divide(2,1);
    cmapData->cd(1);
    TString titledata1 = data_mapL1->GetTitle();
    titledata1+=fTitleData.Data();
    data_mapL1->SetTitle(titledata1.Data());
    data_mapL1->Draw("colz");
    cmapData->cd(2);
    TString titledata2 = data_mapL2->GetTitle();
    titledata2+=fTitleData.Data();
    data_mapL2->SetTitle(titledata2.Data());
    data_mapL2->Draw("colz");

    TCanvas *cmapMc = new TCanvas("cmapMc","cmapMc",1200,600);
    cmapMc->Divide(2,1);

    cmapMc->cd(1);
    TString titlemc1 = mc_mapL1->GetTitle();
    titlemc1+=fTitleMc.Data();
    mc_mapL1->SetTitle(titlemc1.Data());
    mc_mapL1->Draw("colz");
    cmapMc->cd(2);
    TString titlemc2 = mc_mapL2->GetTitle();
    titlemc2+=fTitleMc.Data();
    mc_mapL2->SetTitle(titlemc2.Data());
    mc_mapL2->Draw("colz");
  }
  // booking ratios 
  // projection phi

  TH1D *ratioL1 = new TH1D("ratioL1","Data / MC - Layer 1",mc_projyL1->GetNbinsX(),mc_projyL1->GetXaxis()->GetXmin(),mc_projyL1->GetXaxis()->GetXmax());
  ratioL1->SetXTitle("(NOT GLOBAL)       direction #varphi [cm]");
  TH1D *ratioL2 = new TH1D("ratioL2","Data / MC - Layer 2",mc_projyL2->GetNbinsX(),mc_projyL2->GetXaxis()->GetXmin(),mc_projyL2->GetXaxis()->GetXmax());
  ratioL2->SetXTitle("(NOT GLOBAL)       direction #varphi [cm]");
  // projection z 
 
  TH1D *ratioL1z = new TH1D("ratioL1z","Data / MC - Layer 1",mc_projyL1z->GetNbinsX(),mc_projyL1z->GetXaxis()->GetXmin(),mc_projyL1z->GetXaxis()->GetXmax());
  ratioL1z->SetXTitle("(NOT GLOBAL)       direction z [cm]");
  TH1D *ratioL2z = new TH1D("ratioL2z","Data / MC - Layer 2",mc_projyL2z->GetNbinsX(),mc_projyL2z->GetXaxis()->GetXmin(),mc_projyL2z->GetXaxis()->GetXmax());
  ratioL2->SetXTitle("(NOT GLOBAL)       direction z [cm]");

  // making the ratios
  ratioL1->Divide(data_projyL1,mc_projyL1,mc_projyL1->GetEntries(),data_projyL1->GetEntries(),"e");
  ratioL2->Divide(data_projyL2,mc_projyL2,mc_projyL2->GetEntries(),data_projyL2->GetEntries(),"e");

  ratioL1z->Divide(data_projyL1z,mc_projyL1z,mc_projyL1z->GetEntries(),data_projyL1z->GetEntries(),"e");
  ratioL2z->Divide(data_projyL2z,mc_projyL2z,mc_projyL2z->GetEntries(),data_projyL2z->GetEntries(),"e");

  TCanvas *c = new TCanvas("cRatioPhi","cRatioPhi",1200,600);
  c->Divide(2,1);
  c->cd(1);
  ratioL1->SetYTitle(Form("ratio / %2.2f cm",ratioL1->GetBinWidth(0)));
  ratioL1->Draw();
  c->cd(2);
  ratioL2->SetYTitle(Form("ratio / %2.2f cm",ratioL2->GetBinWidth(0)));
  ratioL2->Draw();

  //TLine *l1z = new TLine(-16.5,1,16.5,1);
  //TLine *l2z = new TLine(-16.5,1,16.5,1);
  TCanvas *cz = new TCanvas("cRatioZ","cRatioZ",1200,600);
  cz->Divide(2,1);
  cz->cd(1);
  ratioL1z->SetYTitle(Form("ratio / %2.2f cm",ratioL1z->GetBinWidth(0)));
  ratioL1z->Draw();
  cz->cd(2);
  ratioL2z->SetYTitle(Form("ratio / %2.2f cm",ratioL2z->GetBinWidth(0)));
  ratioL2z->Draw();

}
void phiTracklet(){

  TCanvas *rawDist = new TCanvas("rawDist"," raw distributions ",1200,800);
  rawDist->Divide(2,2);

  TH2F *trackData = (TH2F*)fListData->FindObject("hSPDphivsSPDeta");
  trackData->SetTitle(Form("%s %s",trackData->GetTitle(),fTitleData.Data()));
  TH1D *trackDataPhi = trackData->ProjectionY();
  if(!trackDataPhi) printf("NO 1 \n");
  //trackDataPhi->SetTitle(Form("%s %s",trackDataPhi->GetTitle(),fTitleData.Data()));
  rawDist->cd(1);
  trackDataPhi->SetLineColor(kRed);
  trackDataPhi->DrawCopy();
  TH1D *trackDataEta = trackData->ProjectionX();
  if(!trackDataEta) printf("NO 2 \n");
  //trackDataEta->SetTitle(Form("%s %s",trackDataEta->GetTitle(),fTitleData.Data()));
  rawDist->cd(2);
  trackDataEta->SetLineColor(kRed);
  trackDataEta->DrawCopy();

  TH1F etaData, phiData;
  trackDataEta->Copy(etaData);
  trackDataPhi->Copy(phiData);

  TH1F etaFrac, phiFrac, mcEta, mcPhi;
  trackDataEta->Copy(etaFrac);
  trackDataPhi->Copy(phiFrac);

  TH2F *trackMc = (TH2F*)fListMc->FindObject("hSPDphivsSPDeta");
  trackMc->SetTitle(Form("%s %s",trackMc->GetTitle(),fTitleMc.Data()));

  TCanvas *tracklets = new TCanvas("tracklets","tracklets",1200,600);
  tracklets->Divide(2,1);
  tracklets->cd(1);
  tracklets->cd(1)->SetRightMargin(0.15);
  //trackData->SetTitle(Form("%s %s",trackData->GetTitle(),fTitleData.Data()));
  trackData->DrawCopy("colz");
  tracklets->cd(2);
  tracklets->cd(2)->SetRightMargin(0.15);
  //trackMc->SetTitle(Form("%s %s",trackMc->GetTitle(),fTitleMc.Data()));
  TH1D *h = (TH1D*)trackMc->DrawCopy("colz");
  fTitleData.ReplaceAll(" ","");
  fTitleMc.ReplaceAll(" ","");
  tracklets->SaveAs(Form("trackletsPhiEtaMaps_%s_%s.png",fTitleData.Data(),fTitleMc.Data()));

  TH1D *trackMcPhi = trackMc->ProjectionY();
  trackMcPhi->SetTitle(Form("%s",h->GetTitle()));
  rawDist->cd(3);
  trackMcPhi->DrawCopy();
  TH1D *trackMcEta = trackMc->ProjectionX();
  trackMcEta->SetTitle(Form("%s",h->GetTitle()));
  rawDist->cd(4);
  trackMcEta->DrawCopy();

  rawDist->SaveAs(Form("trackletsPhiEtaRaw_%s_%s.png",fTitleData.Data(),fTitleMc.Data()));

  TH1F etaMc, phiMc;
  trackMcEta->Copy(etaMc);
  trackMcPhi->Copy(phiMc);

  trackMcEta->Copy(mcEta);
  trackMcPhi->Copy(mcPhi);

  etaFrac.Scale(1./etaFrac.GetEntries());
  mcEta.Scale(1./mcEta.GetEntries());
  etaFrac.Add(&mcEta,-1);
  etaFrac.Divide(&mcEta);

  phiFrac.Scale(1./phiFrac.GetEntries());
  mcPhi.Scale(1./mcPhi.GetEntries());
  phiFrac.Add(&mcPhi,-1);
  phiFrac.Divide(&mcPhi);


  TCanvas *track = new TCanvas("track","track",1200,600);
  track->Divide(2,1);
  track->cd(1);
  phiData.SetLineColor(kRed);
  phiData.SetLineWidth(2);
  phiData.Scale(1./phiData.GetEntries());
  phiData.DrawCopy();
  phiMc.Scale(1./phiMc.GetEntries());
  phiMc.DrawCopy("same");
  track->cd(2);
  etaData.SetLineColor(kRed);
  etaData.SetLineWidth(2);
  etaData.Scale(1./etaData.GetEntries());
  etaData.DrawCopy();
  etaMc.Scale(1./etaMc.GetEntries());
  etaMc.DrawCopy("same");
  track->SaveAs(Form("trackletsPhiEtaNorm_%s_%s.png",fTitleData.Data(),fTitleMc.Data()));

  TCanvas *frac = new TCanvas("frac","fractions",1200,600);
  frac->Divide(2,1);
  frac->cd(1);
  etaFrac.SetTitle(Form(" #Delta#eta/#eta_{%s}   %s - %s ",fTitleMc.Data(),fTitleData.Data(),fTitleMc.Data()));
  etaFrac.SetLineColor(1);
  etaFrac.DrawCopy();
  frac->cd(2);
  phiFrac.SetTitle(Form(" #Delta#varphi/#varphi_{%s}   %s - %s ",fTitleMc.Data(),fTitleData.Data(),fTitleMc.Data()));
  phiFrac.SetLineColor(1);
  phiFrac.DrawCopy();

  frac->SaveAs(Form("relativeRatios_%s_%s.png",fTitleData.Data(),fTitleMc.Data()));

}

void foEfficiency(){
  TH2F *firedFoData = (TH2F*)fListData->FindObject("hFOgoodPerBCmod4");
  TH2F *firedChipsData = (TH2F*)fListData->FindObject("hFiredChipsPerBCmod4");
  TH2F mapBCmod;
  firedFoData->Copy(mapBCmod); 
  mapBCmod.Divide(firedChipsData);
  mapBCmod.SetTitle(Form("FO eff per BCmod4 in %s ",fTitleData.Data()));
  mapBCmod.GetYaxis()->SetNdivisions(4,kFALSE);
  TCanvas *c = new TCanvas("mapFo"," FO eff map",800,800);
  c->cd();
  c->cd()->SetGridy();
  mapBCmod.DrawCopy("colz");

  TH1F bcmod[4];

  TH1F *hbc = new TH1F("bc","bc",firedFoData->GetNbinsX(),(firedFoData->GetXaxis())->GetXmin(),(firedFoData->GetXaxis())->GetXmax());

  for(Int_t bc=0; bc<4; bc++){
    hbc->Clear();
    hbc->Reset();
    for(Int_t iBin=0; iBin<firedFoData->GetNbinsX(); iBin++){
      hbc->SetBinContent(iBin+1,mapBCmod.GetBinContent(iBin+1,bc+1)); 
    }
    hbc->Copy(bcmod[bc]);
    bcmod[bc].SetLineColor(bc+1); 
  }


  TH1F *h[4];
  TCanvas *bceff = new TCanvas("bceff");
  bceff->Divide(2,2);
  for(Int_t iPad=0; iPad<4; iPad++){
   TVirtualPad *pad = bceff->cd(iPad+1);
    if(iPad<2){
      pad->Divide(2,1);
      Int_t idx = -1;
      if(iPad==0) idx=0;
      if(iPad==1) idx=2;  
      pad->cd(1); bcmod[idx].DrawCopy();
      pad->cd(2); bcmod[idx+1].DrawCopy();
    }

    if(iPad==2){
      h[0]= (TH1F*)(bcmod[0].DrawCopy());
      for(Int_t bb=1; bb<4; bb++){
	h[bb] = (TH1F*)(bcmod[bb].DrawCopy("same"));
      }
    }

    TH1F *hdiff[3];

    if(iPad==3){
      hdiff[0]=(TH1F*)(h[1]->Clone());   
      hdiff[0]->Add(h[0],-1);
      hdiff[0]->DrawCopy();
      hdiff[1]=(TH1F*)(h[2]->Clone());   
      hdiff[1]->Add(h[0],-1);
      hdiff[1]->DrawCopy("same");
      hdiff[2]=(TH1F*)(h[3]->Clone());   
      hdiff[2]->Add(h[0],-1);
      hdiff[2]->DrawCopy("same");
    }
  }
}

void phiTrackletsZ(){

  TH1F* phiZposData = (TH1F*)fListData->FindObject("hSPDphiZpos");
  phiZposData->SetLineColor(kRed);
  phiZposData->SetLineWidth(2);
  TH1F* phiZnegData = (TH1F*)fListData->FindObject("hSPDphiZneg");
  phiZnegData->SetLineColor(kRed);
  phiZnegData->SetLineWidth(2);
  TCanvas *cZ = new TCanvas("cZ","cZ",1000,600);

  cZ->Divide(2,1);
  cZ->cd(1);
  phiZposData->Scale(1./phiZposData->GetEntries());
  phiZposData->DrawCopy();
  cZ->cd(2);
  phiZnegData->Scale(1./phiZnegData->GetEntries());
  phiZnegData->DrawCopy();

  TH1F* phiZposMc = (TH1F*)fListMc->FindObject("hSPDphiZpos");
  TH1F* phiZnegMc = (TH1F*)fListMc->FindObject("hSPDphiZneg");
  cZ->cd(1);
  phiZposMc->Scale(1./phiZposMc->GetEntries());
  phiZposMc->DrawCopy("same");
  cZ->cd(2);
  phiZnegMc->Scale(1./phiZnegMc->GetEntries());
  phiZnegMc->DrawCopy("same");

} 

void ratiochips(){

  TH1F *data_chip = ((TH1F*)fListData->FindObject("hClusterYield"));
  data_chip->Sumw2();

  TH1F *mc_chip = ((TH1F*)fListMc->FindObject("hClusterYield"));
  mc_chip->Sumw2();

  TCanvas *c1chip = new TCanvas("c1chip","c1chip",1200,600);
  c1chip->Divide(2,1);
  c1chip->cd(1);
  data_chip->SetTitle(Form("chip cluster yield in %s",fTitleData.Data()));
  data_chip->DrawCopy();
  c1chip->cd(2);
  mc_chip->SetTitle(Form("chip cluster yield in %s",fTitleMc.Data()));
  mc_chip->DrawCopy();


  TH1D *ratiochip  = new TH1D("ratiochip",Form("Chip cluster ratio %s / %s  - scaling : Integral",fTitleData.Data(),fTitleMc.Data()),mc_chip->GetNbinsX(),mc_chip->GetXaxis()->GetXmin(),mc_chip->GetXaxis()->GetXmax());
  ratiochip->GetXaxis()->SetTitle("chip number");
  printf("data_chip %f     -  mc_chip %f    \n",data_chip->GetEntries(),mc_chip->GetEntries());
  ratiochip->Divide(data_chip,mc_chip,mc_chip->Integral(),data_chip->Integral(),"e");
  ratiochip->SetMarkerStyle(20);
  ratiochip->SetMarkerSize(0.7);
  ratiochip->SetMarkerColor(2);

  TCanvas *ratiofitchip = new TCanvas("ratiofitchip","ratiofit",1200,600);
  ratiofitchip->cd();
  TF1 *f1innerChip = new TF1("f1innerChip", "pol0", 0, 399);
  f1innerChip->SetLineColor(kBlue);
  TF1 *f1outerChip = new TF1("f1outerChip", "pol0", 400, 1199);
  f1outerChip->SetLineColor(kGreen);
  
  printf("------------- Fitting standard normalized ratio  ----------------\n");
  ratiochip->Fit("f1innerChip", "R");
  ratiochip->Fit("f1outerChip", "R+");
  //f1inner->GetParameters(&par[0]);
  //f1outer->GetParameters(&par[1]);
  //f1->SetParameters(par);
  //ratiomodule->Fit(f1, "R+");
  f1innerChip->Draw("same");
  f1outerChip->Draw("same");

  TPaveText *parsFitChip= new TPaveText(0.3,0.23,0.6,0.55,"NDC"); 
  parsFitChip->AddText("Inner ");
  parsFitChip->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",f1innerChip->GetChisquare(),f1innerChip->GetNDF()));
  parsFitChip->AddText(Form("Fit :    %4.4f +- %4.4f",f1innerChip->GetParameter(0),f1innerChip->GetParError(0)));
  parsFitChip->AddText("");
  parsFitChip->AddText("Outer ");
  parsFitChip->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",f1outerChip->GetChisquare(),f1outerChip->GetNDF()));
  parsFitChip->AddText(Form("Fit :    %4.4f +- %4.4f",f1outerChip->GetParameter(0),f1outerChip->GetParError(0)));
  ratiochip->GetListOfFunctions()->Add(parsFitChip);



  TCanvas *ratiofitChipEvents = new TCanvas("ratiofitChipEvents","ratiofitEvents",1200,600);

  ratiofitChipEvents->cd();
  TH1D *ratiochipEvents  = new TH1D("ratiochipEvents",Form("Chip cluster ratio %s / %s  - scaling : # events",fTitleData.Data(),fTitleMc.Data()),mc_chip->GetNbinsX(),mc_chip->GetXaxis()->GetXmin(),mc_chip->GetXaxis()->GetXmax());
  ratiochipEvents->GetXaxis()->SetTitle("chip number");
  Double_t nEvData = ((TH1I*)(fListData->FindObject("hEventsProcessed")))->GetEntries();
  Double_t nEvMc= ((TH1I*)(fListMc->FindObject("hEventsProcessed")))->GetEntries();
  ratiochipEvents->Divide(data_chip,mc_chip,nEvMc,nEvData,"e");
  ratiochipEvents->Draw();
  ratiochipEvents->SetMarkerStyle(20);
  ratiochipEvents->SetMarkerSize(0.7);
  ratiochipEvents->SetMarkerColor(2);

  TF1 *fInnerChip = new TF1("fInnerChip", "pol0", 0, 399);
  fInnerChip->SetLineColor(kBlue);
  TF1 *fOuterChip = new TF1("fOuterChip", "pol0", 400, 1199);
  fOuterChip->SetLineColor(kGreen);
  ratiochipEvents->Fit("fInnerChip", "R");
  ratiochipEvents->Fit("fOuterChip", "R+");
  fInnerChip->Draw("same");
  fOuterChip->Draw("same");
  
  
  TPaveText *parsFit= new TPaveText(0.3,0.23,0.6,0.55,"NDC"); 
  parsFit->AddText("Inner ");
  parsFit->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",fInnerChip->GetChisquare(),fInnerChip->GetNDF()));
  parsFit->AddText(Form("Fit :    %4.4f +- %4.4f",fInnerChip->GetParameter(0),fInnerChip->GetParError(0)));
  parsFit->AddText("");
  parsFit->AddText("Outer ");
  parsFit->AddText(Form("#chi2 / ndf    %4.3f / %3.1i",fOuterChip->GetChisquare(),fOuterChip->GetNDF()));
  parsFit->AddText(Form("Fit :    %4.4f +- %4.4f",fOuterChip->GetParameter(0),fOuterChip->GetParError(0)));
  ratiochipEvents ->GetListOfFunctions()->Add(parsFit);
  
  
  //------------------------------------------
  
  TH1D * diffsClus[2];
  diffsClus[0]= new TH1D("diffsL1Clus"," ",80,-0.2,0.2);
  diffsClus[1]= new TH1D("diffsL2Clus"," ",80,-0.2,0.2);
  for(Int_t ibin=0; ibin<1200; ibin++){

    if(ibin<400){
      if(ratiochip->GetBinContent(ibin+1)>0) diffsClus[0]->Fill(ratiochip->GetBinContent(ibin+1)-f1innerChip->GetParameter(0));
    }else {

      if(ratiochip->GetBinContent(ibin+1)>0) diffsClus[1]->Fill(ratiochip->GetBinContent(ibin+1)-f1outerChip->GetParameter(0));

    }


  }


  TCanvas * pullsClus = new TCanvas("pullsClus","pullsClus",1200,600);
  pullsClus->Divide(2,1);
  pullsClus->cd(1);
  diffsClus[0]->SetTitle("dispersion Layer 1");
  diffsClus[0]->Rebin(2);
  diffsClus[0]->Fit("gaus","","",-0.2,0.2);
  if(diffsClus[0]->GetFunction("gaus")) diffsClus[0]->GetFunction("gaus")->SetLineColor(kBlue);
  diffsClus[0]->Draw();
  pullsClus->cd(2);
  diffsClus[1]->SetTitle("dispersion Layer 2");
  diffsClus[1]->Rebin(2);
  diffsClus[1]->Fit("gaus","","",-0.2,0.2);
  if(diffsClus[1]->GetFunction("gaus")) diffsClus[1]->GetFunction("gaus")->SetLineColor(kBlue);
  diffsClus[1]->Draw();

}
