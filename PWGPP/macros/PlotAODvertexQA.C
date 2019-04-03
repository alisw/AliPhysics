#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TSystem.h>
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>
#endif

void Draw2D(TH2F* h, Double_t maxx);
TGraphErrors* GetResolVsMult(TH2F* h, TString axis, Int_t uConv=1);

void PlotAODvertexQA(TString filename="QAresults_AOD.root", TString suffix="QA", Int_t runNumber=-1){
  
  const Int_t totTrending=16;
  Float_t vecForTrend[totTrending];
  TString varForTrending[totTrending]={"fTrackV","fSPD3D","fSPDz","fTPC","fInvalid",
				       "meanx","meany","meanz","emeanx","emeany","emeanz","meancont",
				       "fracpilSPD","fracselpilSPD","fracpilMV","fracselpilMV"};
  
  TTree* trtree=new TTree("trendingVert","tree of trending variables");
  trtree->Branch("nrun",&runNumber,"nrun/I");
  for(Int_t j=0; j<totTrending; j++){
    trtree->Branch(varForTrending[j].Data(),&vecForTrend[j],Form("%s/F",varForTrending[j].Data()));
    vecForTrend[j]=-99.;
  }

  TFile* f=new TFile(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)f->Get("CheckVertexAOD");
  if(!df){
    printf("Directory CheckVertexAOD not found in file %s\n",filename.Data());
    return;
  }
  TList* l=(TList*)df->Get(Form("clistCheckVertexAOD%s",suffix.Data()));
  if(!l){
    printf("TList clistCheckVertexAOD%s not found in file %s\n",suffix.Data(),filename.Data());
    return;    
  }

  gStyle->SetOptTitle(0);

  TH1F*	hvt=(TH1F*)l->FindObject("hAllVtxType");
  TH1F*	hpvt=(TH1F*)l->FindObject("hPrimVtxType");
  hvt->SetStats(0);
  hpvt->SetStats(0);
  hvt->GetXaxis()->SetTitleOffset(1.7);
  hpvt->GetXaxis()->SetTitleOffset(1.7);
  hvt->GetYaxis()->SetTitleOffset(1.5);
  hpvt->GetYaxis()->SetTitleOffset(1.5);

  TCanvas* cty=new TCanvas("cty","VertexType",1400,600);
  cty->Divide(2,1);
  cty->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.07);
  hvt->Draw();
  cty->cd(2);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.07);
  hpvt->Draw();
  for(Int_t ib=1; ib<=hpvt->GetNbinsX(); ib++){
    Double_t fr=hpvt->GetBinContent(ib)/hpvt->GetEntries()*100.;
    if(fr>0.09){
      TLatex* t=new TLatex(hpvt->GetBinCenter(ib)-0.3*hpvt->GetBinWidth(ib),hpvt->GetBinContent(ib)+0.015*(hpvt->GetMaximum()),Form("%.1f%%",fr));
      t->SetTextColor(kGray+1);
      t->SetTextFont(43);
      t->SetTextSize(18);
      t->Draw();
    }
  }
  cty->SaveAs("VertexAOD-Type.png");

  vecForTrend[0]=hpvt->GetBinContent(3)/hpvt->GetEntries();
  vecForTrend[1]=hpvt->GetBinContent(4)/hpvt->GetEntries();
  vecForTrend[2]=hpvt->GetBinContent(5)/hpvt->GetEntries();
  vecForTrend[3]=hpvt->GetBinContent(6)/hpvt->GetEntries();
  vecForTrend[4]=hpvt->GetBinContent(1)/hpvt->GetEntries();

  Int_t kcm2um=10000;
  Int_t kcm2mm=10;

  TH2F* hxspd=(TH2F*)l->FindObject("hXspdVsMult");
  if(!hxspd) hxspd=(TH2F*)l->FindObject("hXspdVsContrib");
  TH2F* hyspd=(TH2F*)l->FindObject("hYspdVsMult");
  if(!hyspd) hyspd=(TH2F*)l->FindObject("hYspdVsContrib");
  TH2F* hzspd=(TH2F*)l->FindObject("hZspdVsMult");
  if(!hzspd) hzspd=(TH2F*)l->FindObject("hZspdVsContrib");
  TH2F* hxtrk=(TH2F*)l->FindObject("hXtrkVsMult");
  if(!hxtrk) hxtrk=(TH2F*)l->FindObject("hXtrkVsContrib");
  TH2F* hytrk=(TH2F*)l->FindObject("hYtrkVsMult");
  if(!hytrk) hytrk=(TH2F*)l->FindObject("hYtrkVsContrib");
  TH2F* hztrk=(TH2F*)l->FindObject("hZtrkVsMult");
  if(!hztrk) hztrk=(TH2F*)l->FindObject("hZtrkVsContrib");

  Double_t maxContrib=20;
  TH1D* hContrib=hztrk->ProjectionX();
  for(Int_t jb=1; jb<=hContrib->GetNbinsX(); jb++){
    if(hContrib->GetBinContent(jb)>0) maxContrib=hContrib->GetBinLowEdge(jb+1);
  }

  TCanvas* cvs=new TCanvas("cvs","SPDVertex",1400,800);
  cvs->Divide(3,2);
  cvs->cd(1);
  Draw2D(hxspd,maxContrib);
  cvs->cd(2);
  Draw2D(hyspd,maxContrib);
  TLatex* tspd=new TLatex(0.16,0.93,"SPD Vertexer");
  tspd->SetNDC();
  tspd->SetTextFont(43);
  tspd->SetTextSize(26);
  tspd->Draw();
  cvs->cd(3);
  Draw2D(hzspd,maxContrib);
  cvs->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gxs=GetResolVsMult(hxspd,"x",kcm2um);
  gxs->GetXaxis()->SetLimits(0.,maxContrib);
  gxs->SetMinimum(0);
  gxs->SetMaximum(600.);
  gxs->Draw("APZ");
  cvs->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gys=GetResolVsMult(hyspd,"y",kcm2um);
  gys->GetXaxis()->SetLimits(0.,maxContrib);
  gys->SetMinimum(0);
  gys->SetMaximum(600.);
  gys->Draw("APZ");
  cvs->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gzs=GetResolVsMult(hzspd,"z");
  gzs->GetXaxis()->SetLimits(0.,maxContrib);
  gzs->SetMinimum(0);
  gzs->SetMaximum(10.);
  gzs->Draw("APZ");
  cvs->SaveAs("VertexAOD-SPD.png");

  //-------------------


  vecForTrend[5]=hxtrk->GetMean(2);
  vecForTrend[6]=hytrk->GetMean(2);
  vecForTrend[7]=hztrk->GetMean(2);
  vecForTrend[8]=hxtrk->GetMeanError(2);
  vecForTrend[9]=hytrk->GetMeanError(2);
  vecForTrend[10]=hztrk->GetMeanError(2);
  vecForTrend[11]=hztrk->GetMean(1);

  TCanvas* cvt=new TCanvas("cvt","TrackVertex",1400,800);
  cvt->Divide(3,2);
  cvt->cd(1);
  Draw2D(hxtrk,maxContrib);
  cvt->cd(2);
  Draw2D(hytrk,maxContrib);
  TLatex* ttrk=new TLatex(0.16,0.93,"Track Vertexer");
  ttrk->SetNDC();
  ttrk->SetTextFont(43);
  ttrk->SetTextSize(26);
  ttrk->Draw();
  cvt->cd(3);
  Draw2D(hztrk,maxContrib);
  cvt->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gxt=GetResolVsMult(hxtrk,"x",kcm2um);
  gxt->GetXaxis()->SetLimits(0.,maxContrib);
  gxt->SetMinimum(0);
  gxt->SetMaximum(600.);
  gxt->Draw("APZ");
  cvt->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gyt=GetResolVsMult(hytrk,"y",kcm2um);
  gyt->GetXaxis()->SetLimits(0.,maxContrib);
  gyt->SetMinimum(0);
  gyt->SetMaximum(600.);
  gyt->Draw("APZ");
  cvt->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gzt=GetResolVsMult(hztrk,"z");
  gzt->GetXaxis()->SetLimits(0.,maxContrib);
  gzt->SetMinimum(0);
  gzt->SetMaximum(10.);
  gzt->Draw("APZ");
  cvt->SaveAs("VertexAOD-Tracks.png");
 
  //----------

  TH2F* hxtpc=(TH2F*)l->FindObject("hXtpcVsMult");
  if(!hxtpc) hxtpc=(TH2F*)l->FindObject("hXtpcVsContrib");
  TH2F* hytpc=(TH2F*)l->FindObject("hYtpcVsMult");
  if(!hytpc) hytpc=(TH2F*)l->FindObject("hYtpcVsContrib");
  TH2F* hztpc=(TH2F*)l->FindObject("hZtpcVsMult");
  if(!hztpc) hztpc=(TH2F*)l->FindObject("hZtpcVsContrib");

  TCanvas* cvp=new TCanvas("cvp","TPCVertex",1400,800);
  cvp->Divide(3,2);
  cvp->cd(1);
  Draw2D(hxtpc,maxContrib);
  cvp->cd(2);
  Draw2D(hytpc,maxContrib);
  TLatex* ttpc=new TLatex(0.16,0.93,"TPC Vertexer");
  ttpc->SetNDC();
  ttpc->SetTextFont(43);
  ttpc->SetTextSize(26);
  ttpc->Draw();
  cvp->cd(3);
  Draw2D(hztpc,maxContrib);
  cvp->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gxp=GetResolVsMult(hxtpc,"x");
  gxp->GetXaxis()->SetLimits(0.,maxContrib);
  gxp->SetMinimum(0);
  gxp->SetMaximum(1.);
  gxp->Draw("APZ");
  cvp->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gyp=GetResolVsMult(hytpc,"y");
  gyp->GetXaxis()->SetLimits(0.,maxContrib);
  gyp->SetMinimum(0);
  gyp->SetMaximum(1.);
  gyp->Draw("APZ");
  cvp->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  TGraphErrors* gzp=GetResolVsMult(hztpc,"z");
  gzp->GetXaxis()->SetLimits(0.,maxContrib);
  gzp->SetMinimum(0);
  gzp->SetMaximum(10.);
  gzp->Draw("APZ");
  cvp->SaveAs("VertexAOD-TPC.png");

  //-----

  TH1F*	hpils=(TH1F*)l->FindObject("hNOfPileupVertSPD");
  TH1F*	hselpils=(TH1F*)l->FindObject("hNOfSelPileupVertSPD");
  TH1F*	hpilm=(TH1F*)l->FindObject("hNOfPileupVertMV");
  TH1F*	hselpilm=(TH1F*)l->FindObject("hNOfSelPileupVertMV");
  hpils->SetLineWidth(2);
  hpils->SetLineColor(1);
  hselpils->SetLineColor(2);
  hselpils->SetFillStyle(1001);
  hselpils->SetFillColor(kRed-9);
  hpilm->SetLineWidth(2);
  hpilm->SetLineColor(1);
  hselpilm->SetLineColor(2);
  hselpilm->SetFillStyle(1001);
  hselpilm->SetFillColor(kRed-9);
  hselpils->SetMinimum(0.2);
  hselpilm->SetMinimum(0.2);
  hselpils->GetXaxis()->SetTitle("Number of pileup vertices");
  hselpilm->GetXaxis()->SetTitle("Number of pileup vertices");
  hselpils->GetYaxis()->SetTitle("Entries");
  hselpilm->GetYaxis()->SetTitle("Entries");
  hselpils->GetYaxis()->SetTitleOffset(1.1);
  hselpilm->GetYaxis()->SetTitleOffset(1.1);

  TCanvas* cpil=new TCanvas("cpil","Pileup",1400,600);
  cpil->Divide(2,1);
  cpil->cd(1);
  gPad->SetLogy();
  hselpils->Draw();
  gPad->Update();
  TPaveStats* st1=(TPaveStats*)hselpils->GetListOfFunctions()->FindObject("stats");
  st1->SetTextColor(2);
  st1->SetY1NDC(0.51);
  st1->SetY2NDC(0.7);
  hpils->Draw("SAMES");
  gPad->Update();
  TPaveStats* st2=(TPaveStats*)hpils->GetListOfFunctions()->FindObject("stats");
  st2->SetTextColor(1);
  st2->SetY1NDC(0.73);
  st2->SetY2NDC(0.92);
  gPad->Modified();
  TLatex* t1=new TLatex(0.3,0.84,"SPD-vertex based tagging");
  t1->SetNDC();
  t1->SetTextFont(43);
  t1->SetTextSize(24);
  t1->Draw();
  TLegend* leg=new TLegend(0.3,0.68,0.68,0.8);
  leg->AddEntry(hpils,"Reconstructed","F");
  leg->AddEntry(hselpils,"Selected","F");
  leg->Draw();
  cpil->cd(2);
  gPad->SetLogy();
  hselpilm->Draw();
  gPad->Update();
  TPaveStats* st3=(TPaveStats*)hselpilm->GetListOfFunctions()->FindObject("stats");
  st3->SetTextColor(2);
  st3->SetY1NDC(0.51);
  st3->SetY2NDC(0.7);
  hpilm->Draw("SAMES");
  gPad->Update();
  TPaveStats* st4=(TPaveStats*)hpilm->GetListOfFunctions()->FindObject("stats");
  st4->SetTextColor(1);
  st4->SetY1NDC(0.73);
  st4->SetY2NDC(0.92);
  gPad->Modified();
  TLatex* t2=new TLatex(0.3,0.84,"Track-vertex based tagging");
  t2->SetNDC();
  t2->SetTextFont(43);
  t2->SetTextSize(24);
  t2->Draw();
  leg->Draw();
  cpil->SaveAs("VertexAOD-Pileup.png");

  vecForTrend[12]=hpils->GetMean();
  vecForTrend[13]=hselpils->GetMean();
  vecForTrend[14]=hpilm->GetMean();
  vecForTrend[15]=hselpilm->GetMean();

  trtree->Fill();

  if(runNumber>0){
    TFile* foutfile=new TFile("trendingAODvertex.root","recreate");
    trtree->Write();
    TDirectory* outdir=foutfile->mkdir(df->GetName());
    outdir->cd();
    l->Write(l->GetName(),1);
    foutfile->Close();
    delete foutfile;
  }
}

void Draw2D(TH2F* h, Double_t maxx){
  h->GetXaxis()->SetRangeUser(0.,maxx);
  TProfile* p=h->ProfileX(Form("%s_prof",h->GetName()));
  p->SetMarkerStyle(20);
  h->SetStats(0);
  h->GetYaxis()->SetTitleOffset(1.2);
  gPad->SetTickx();
  gPad->SetTicky();
  h->Draw("col");
  p->Draw("psame");
}


TGraphErrors* GetResolVsMult(TH2F* h, TString axis, Int_t uConv){
  TGraphErrors* gr=new TGraphErrors(0);
  Int_t nPoint=0;
  for(Int_t jb=1; jb<=h->GetNbinsX(); jb++){
    TH1D* htmp=h->ProjectionY("htmp",jb,jb);
    if(htmp->GetEntries()>30){
      htmp->Fit("gaus","Q0");
      TF1* fg=(TF1*)htmp->GetListOfFunctions()->FindObject("gaus");
      gr->SetPoint(nPoint,h->GetXaxis()->GetBinCenter(jb),fg->GetParameter(2)*uConv);
      gr->SetPointError(nPoint,0.5*h->GetXaxis()->GetBinWidth(jb),fg->GetParError(2)*uConv);
      nPoint++;
    }
    delete htmp;
  }
  gr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  if(uConv==10000) gr->GetYaxis()->SetTitle(Form("#sigma_{%s} (#mum)",axis.Data()));
  else if(uConv==10) gr->GetYaxis()->SetTitle(Form("#sigma_{%s} (mm)",axis.Data()));
  else if(uConv==1) gr->GetYaxis()->SetTitle(Form("#sigma_{%s} (cm)",axis.Data()));
  gr->GetYaxis()->SetTitleOffset(1.2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  return gr;
}
