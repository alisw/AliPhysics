//############################################################################################
//
//   Macro to manage the average of D0,D+,D*+ meson - charged particle azimuthal correlations
//   andrea.rossi@cern.ch 
//
//
//############################################################################################

TString codeDir=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");
void SetFitPlotMacroPath(TString strdir){
  codeDir=strdir;
}
void DrawPlotWithSystUnc(TH1D *h){
  AliHFDhadronCorrSystUnc *oUnc=new AliHFDhadronCorrSystUnc();
  oUnc->InitStandardUncertainties2010(0,10.,0.3);
  oUnc->BuildSystUncertaintyPlotVsDeltaPhi(h,1);
  
}

void MakeAverage(Double_t ptmin,Double_t ptmax,Double_t ptassocmin,Double_t ptassocmax,Int_t system=0/*0=pp, 1= pPb*/,Int_t year=2010,Int_t reflected=0,Double_t maxH=7./*just for graphics*/,Bool_t doArithmeticAv=kFALSE,Int_t baselineOpt=5,TString strDzero="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDzero/pp/Singlev2Envelope",TString strDplus="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDplus/pp",TString strDstar="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2014Nov9paperProp/InputPlotsDStar/pp/CorrectedPlots",TString v2suffix="_v2D0.00_v2had0.00"){

  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/LoadLibraries.C");
  gROOT->LoadMacro(Form("%s/FitPlots.C",codeDir.Data()));
  
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);

  TString systemStr;
  if(system==0)systemStr="pp";
  else if(system==1)systemStr="pPb";
  else if(system==2)systemStr="pPb";
  else {
    Printf("Make Average: WRONG SYSTEM INPUT");
  }
  TString histname="hDataCorrectedTempl0CentrFprompt";
  if(reflected==1)histname="hDataCorrectedTempl0CentrFpromptReflected";
  // Dzero
  TFile *f=0x0;
  f=TFile::Open(Form("%s/CanvaAndVariedHisto%sDzeroPt%.0fto%.0fassocPt%.1fto%.1f%s.root",strDzero.Data(),systemStr.Data(),ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()),"READ");
  if(!f)return;
  TH1D *hDzero=f->Get(histname.Data());
  hDzero->SetName(Form("Dzero%.0fto%.0fassoc%.0fto%.0f",ptmin,ptmax,ptassocmin*10,ptassocmax*10,v2suffix.Data()));
    
  AliHFDhadronCorrSystUnc *systDzero=(AliHFDhadronCorrSystUnc*)f->Get("SystematicUncertainty");
  systDzero->SetName("systDzero");
    
  // DSTAR
  f=0x0;
  f=TFile::Open(Form("%s/CanvaAndVariedHisto%sDstarPt%.0fto%.0fassocPt%.1fto%.1f%s.root",strDstar.Data(),systemStr.Data(),ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()),"READ");

  if(!f) return;
  TH1D *hDstar=f->Get(histname.Data());
  hDstar->SetName(Form("Dstar%.0fto%.0fassoc%.1fto%.1f",ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()));

  AliHFDhadronCorrSystUnc *systDstar=(AliHFDhadronCorrSystUnc*)f->Get("SystematicUncertainty");
  systDstar->SetName("systDstar");

  // DPLUS
  f=0x0;
  f=TFile::Open(Form("%s/CanvaAndVariedHisto%sDplusPt%.0fto%.0fassocPt%.1fto%.1f%s.root",strDplus.Data(),systemStr.Data(),ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()),"READ");
  if(!f)return;
  TH1D *hDplus=f->Get(histname.Data());
  hDplus->SetName(Form("Dplus%.0fto%.0fassoc%.1fto%.1f",ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()));
  AliHFDhadronCorrSystUnc *systDplus=(AliHFDhadronCorrSystUnc*)f->Get("SystematicUncertainty");
  systDplus->SetName("systDplus");

  TCanvas *cInput=new TCanvas("cInput","cInput",800,800);
  cInput->cd();
  hDzero->Draw();
  hDplus->Draw("Same");
  hDstar->Draw("same");


  AliHFDmesonCorrAverage *av=new AliHFDmesonCorrAverage();
  av->SetSystem(system,year);
  av->SetMethod(10);
  av->SetArithmeticAverage(doArithmeticAv);
  av->SetMomentumRanges(ptmin,ptmax,ptassocmin,ptassocmax);
  av->SetIncludeDzero(kTRUE);
  av->SetIncludeDstar(kTRUE);
  av->SetIncludeDplus(kTRUE);
  av->SetDzeroHisto(hDzero);
  av->SetDstarHisto(hDstar);
  av->SetDplusHisto(hDplus);
  av->SetDzeroSystUnc(systDzero);
  av->SetDstarSystUnc(systDstar);
  av->SetDplusSystUnc(systDplus);
  av->SetSystAreAlreadySet(kTRUE);
  av->CalculateAverage();	
  TH1D *hWeightDzero=av->GetWeightsUsedDzero();
  TH1D *hWeightDplus=av->GetWeightsUsedDplus();
  TH1D *hWeightDstar=av->GetWeightsUsedDstar();

  TCanvas *cUsedWeights=new TCanvas("cUsedWeights","cUsedWeights",800,800);
  cUsedWeights->cd();
  hWeightDzero->SetLineColor(kBlack);
  hWeightDzero->Draw();
  hWeightDstar->Draw("same");
  hWeightDstar->SetLineColor(kBlue);
  hWeightDplus->Draw("same");
  hWeightDplus->SetLineColor(kRed);

  AliHFDhadronCorrSystUnc *avSyst=av->GetAverageSystUncertainty();
  avSyst->SetName("AverageSystematicUncertainty");

  TGraphAsymmErrors *grTotSyst=avSyst->GetTotUncGraph();
  //BuildSystUncertaintyPlotVsDeltaPhi(h);

  TCanvas *c=new TCanvas(Form("c%.0fto%.0f",ptmin,ptmax),Form("c%.0fto%.0f",ptmin,ptmax),700,700);
  c->SetTicks();
  c->cd();
  TH1D *hAv=av->GetAverageHisto();
  hAv->SetTitle("");
  hAv->SetMaximum(maxH);
  hAv->Draw();
//   hAv->SetLineColor(kBlack);
//   hAv->SetLineWidth(1);
//   hAv->SetLineStyle(1);
//   hAv->SetMarkerSize(1);
//   hAv->SetMarkerStyle(20);
//   hAv->SetMarkerColor(kBlack);

  hDzero->SetTitle(Form("Dzero%.0fto%.0f",ptmin,ptmax));
  hDzero->SetLineColor(kRed);
  hDzero->SetLineWidth(2);
  hDzero->SetLineStyle(2);
  hDzero->SetMarkerSize(1);
  hDzero->SetMarkerStyle(28);
  hDzero->SetMarkerColor(kRed);
  hDzero->Draw("same");

  hDstar->SetTitle(Form("Dstar%.0fto%.0f",ptmin,ptmax));
  hDstar->SetLineColor(kBlue);
  hDstar->SetLineWidth(2);
  hDstar->SetLineStyle(2);
  hDstar->SetMarkerSize(1);
  hDstar->SetMarkerStyle(25);
  hDstar->SetMarkerColor(kBlue);
  hDstar->Draw("same");

  hDplus->SetTitle(Form("Dplus%.0fto%.0f",ptmin,ptmax));
  hDplus->SetLineColor(kBlue);
  hDplus->SetLineWidth(2);
  hDplus->SetLineStyle(2);
  hDplus->SetMarkerSize(1);
  hDplus->SetMarkerStyle(20);
  hDplus->SetMarkerColor(kGreen);
  hDplus->Draw("same");
  
  c->BuildLegend();
  
  TCanvas *cRef=new TCanvas(Form("cRef%.0fto%.0f",ptmin,ptmax),Form("cRef%.0fto%.0f",ptmin,ptmax),700,700);
  cRef->cd();
  cRef->SetTicks();
  TH1D *hAvRefl=av->ReflectHisto(hAv);
  hAvRefl->Draw();


  avSyst->BuildSystUncertaintyPlotVsDeltaPhi(hAv,1);

  TCanvas *cFit=new TCanvas(Form("cFit%.0fto%.0f",ptmin,ptmax),Form("cFit%.0fto%.0f",ptmin,ptmax),700,700);
  cFit->cd();
  cFit->SetTicks();
  TH1D *hAvFit=(TH1D*)hAv->Clone(Form("%sCp",hAv->GetName()));
  TF1 *fit=FitPlotsShort(hAvFit,1,baselineOpt,3,1);
  fit->SetName(Form("fit%.0fto%.0f",ptmin,ptmax));
  DrawLegendWithParameters(fit,cFit);
  AddTextInfo(cFit,ptmin,ptmax,ptassocmin,ptassocmax);
  TGraphAsymmErrors *grTotSyst=avSyst->GetTotUncGraph();
  grTotSyst->SetName("grTotSyst");
  grTotSyst->Draw("E2");



  TFile *fout;
  if(doArithmeticAv)fout=new TFile(Form("ArithmeticAverage%sDzeroDstarDplus%.0fto%.0f_assoc%.1fto%.1f%s.root",systemStr.Data(),ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()),"RECREATE");
  else fout=new TFile(Form("WeightedAverage%sDzeroDstarDplus%.0fto%.0f_assoc%.1fto%.1f%s.root",systemStr.Data(),ptmin,ptmax,ptassocmin,ptassocmax,v2suffix.Data()),"RECREATE");
  grTotSyst->Write();
  avSyst->Write();
  hAv->Write();
  cRef->Write();
  cFit->Write();
  c->Write();
  cUsedWeights->Write();
  fout->Close();
}




void OpenOutputFileAndDraw(TString strfile,Double_t ptminD,Double_t ptmaxD,TString strMeson="D",Int_t system=0,Double_t ptminAss=0.3,Double_t ptmaxAss=99.0,Double_t deltaeta=1,TString avString="Weighted",Double_t maxH=7./*just for graphics*/){
  gStyle->SetOptStat(0000);

  TFile *f=TFile::Open(strfile.Data(),"READ");
  AliHFDhadronCorrSystUnc *syst=(AliHFDhadronCorrSystUnc*)f->Get("AverageSystematicUncertainty");
  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();
  TH1D *hFDsub=(TH1D*)f->Get("fhDaverage");
  TGraphAsymmErrors *gr=syst->GetTotNonFlatUncGraph();
  gStyle->SetOptStat(0000);
  
  TCanvas *cDraw=new TCanvas("cDraw","cDraw",700,700);
  cDraw->cd();
  cDraw->SetLeftMargin(0.15);
  cDraw->SetRightMargin(0.05);
  cDraw->SetTicks();

  hFDsub->SetLineColor(kBlack);
  hFDsub->SetMarkerColor(kBlack);
  hFDsub->SetXTitle("#Delta#phi (rad)");
  hFDsub->SetYTitle(Form("#frac{1}{N_{%s}}#frac{dN^{assoc}}{d#Delta#phi} (rad^{-1})",strMeson.Data()));
  hFDsub->GetYaxis()->SetTitleOffset(1.3);
  hFDsub->GetYaxis()->SetRangeUser(0,maxH);
  hFDsub->SetTitle("");
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);
  hFDsub->Draw();
  gr->Draw("E2");
  TLatex *tSystem=new TLatex(0.18,0.80,"#bf{pp, #sqrt{#it{s}} = 7 TeV, L_{int} = 5 nb^{-1}}");
  if(system==1 || system==2) tSystem->SetTitle("#bf{p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, L_{int} = 50 #mub^{-1}}");
  tSystem->SetNDC();
  tSystem->SetTextSize(0.03);
  tSystem->Draw();
  
  TLatex *tptD=new TLatex(0.6,0.83,Form("#bf{%.1f<#it{p}_{T}^{%s}<%.1f GeV/#it{c}}",ptminD,ptmaxD,strMeson.Data()));
  tptD->SetNDC();
  tptD->SetTextSize(0.04);
  tptD->Draw();

  TLatex *tptAssoc;
  if(ptmaxAss!=99.0){
    tptAssoc=new TLatex(0.6,0.78,Form("#bf{%.1f<#it{p}_{T}^{assoc}<%.1f GeV/#it{c}}",ptminAss,ptmaxAss));
  }
  else tptAssoc=new TLatex(0.6,0.78,Form("#bf{#it{p}_{T}^{assoc}>%.1f GeV/#it{c}}",ptminAss));
  tptAssoc->SetNDC();
  tptAssoc->SetTextSize(0.04);
  tptAssoc->Draw();
  TLatex *tDEta=new TLatex(0.6,0.73,Form("#bf{|#Delta#eta|<%.1f}",deltaeta));
  tDEta->SetNDC();
  tDEta->SetTextSize(0.04);
  tDEta->Draw();

printf("RANGE MAX = %d\n",maxH);
  TLatex *tUncertainty;
  if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("%.0f#% correlated uncertainty",hUncCorrMin->GetBinContent(1)*100.));
  else tUncertainty=new TLatex(0.18,0.68,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.04);
  tUncertainty->Draw();

  
  TCanvas *cVaryHisto=new TCanvas("cVaryHisto","cVaryHisto",700,700);
  cVaryHisto->cd();
  hFDsub->DrawCopy();
  TH1D *hVaryUp=syst->GetVariedHisto(hFDsub,gr,1);
  TH1D *hVaryDown=syst->GetVariedHisto(hFDsub,gr,0);
  hVaryUp->Draw("same");
  hVaryDown->Draw("same");

  TString strfileout="CanvaAndVariedHisto";
  strfileout.Append(Form("%sAverage",avString.Data()));
  if(strfile.Contains("Dzero"))strfileout.Append("Dzero");
  if(strfile.Contains("Dstar"))strfileout.Append("Dstar");
  if(strfile.Contains("Dplus"))strfileout.Append("Dplus");
  if(system==0)strfileout.Append("_pp_");
  else if(system==1)strfileout.Append("_pPb_");
  else if(system==2)strfileout.Append("_pPb_");
  else strfileout.Append("WrongCollSyst");
  strfileout.Append(Form("Pt%.0fto%.0fassocPt%.1fto",ptminD,ptmaxD,ptminAss));
  if(ptmaxAss>0)strfileout.Append(Form("%.1f.root",ptmaxAss));
  else strfileout.Append("99.0.root");

  cDraw->Update();
  strfileout.ReplaceAll(".root",".png");
  cDraw->Print(strfileout.Data());
  strfileout.ReplaceAll(".png",".root");
  TFile *fout=new TFile(strfileout.Data(),"RECREATE");
  fout->cd();
  cDraw->Write();
  hVaryUp->Write();
  hVaryDown->Write();
  fout->Close();
}

void SetHistoStyle(TH1* h,Int_t color,Int_t markerStyle=20,Int_t lineStyle=1,Int_t lineWidth=1,Int_t fillStyle=-1){

  h->SetLineColor((Color_t)color);

  h->SetMarkerColor(color);
  if(fillStyle>0)h->SetFillColor(color);

 h->SetMarkerStyle(markerStyle);

  h->SetLineStyle(lineStyle);
  h->SetLineWidth(lineWidth);
  
  if(fillStyle>0)h->SetFillStyle(fillStyle);

}

void SetGraphStyle(TGraph *gr,Int_t color,Int_t markerStyle=20,Int_t lineStyle=1,Int_t lineWidth=1,Int_t fillStyle=-1){

  gr->SetLineColor((Color_t)color);

  gr->SetMarkerColor(color);
  if(fillStyle>0)gr->SetFillColor(color);

 gr->SetMarkerStyle(markerStyle);

  gr->SetLineStyle(lineStyle);
  gr->SetLineWidth(lineWidth);
  
  if(fillStyle>0)gr->SetFillStyle(fillStyle);

}

void CompareHistosFinalOutput(TString strFile1,TString strFile2,TString strTitle1="Arithmetic Av",TString strTitle2="Weighted Av"){
  TFile *f1=TFile::Open(strFile1.Data(),"READ");
  TH1D *hAv1=(TH1D*)f1->Get("fhDaverage");
  SetHistoStyle(hAv1,kBlack,20);
  hAv1->SetTitle(strTitle1.Data());
  TH1D *hAbsStat1=(TH1D*)hAv1->Clone("hAbsStat1");
  TH1D *hRelStat1=(TH1D*)hAv1->Clone("hRelStat1");
  for(Int_t j=1;j<=hAbsStat1->GetNbinsX();j++){
    hAbsStat1->SetBinContent(j,0);
    hRelStat1->SetBinError(j,hRelStat1->GetBinError(j)/hRelStat1->GetBinContent(j));
    hRelStat1->SetBinContent(j,0);
  }
  TGraphAsymmErrors *gr1=(TGraphAsymmErrors*)f1->Get("grTotSyst");
  gr1->SetTitle(strTitle1.Data());
  SetGraphStyle(gr1,kBlack,20);
  TGraphAsymmErrors *grRelSyst1=(TGraphAsymmErrors*)gr1->Clone("grRelSyst1");  
  //  grRelSyst1->SetTitle(strTitle1.Data());
  Int_t npoints=gr1->GetN();
  Double_t x,y;
  for(Int_t i=0;i<npoints;i++){
    gr1->GetPoint(i,x,y);
    gr1->SetPoint(i,x,0);
    grRelSyst1->SetPoint(i,x,0);
    grRelSyst1->SetPointEYlow(i,grRelSyst1->GetErrorYlow(i)/y);
    grRelSyst1->SetPointEYhigh(i,grRelSyst1->GetErrorYhigh(i)/y);
  }
  Printf("First file: histos loaded");

  TFile *f2=TFile::Open(strFile2.Data(),"READ");
  TH1D *hAv2=(TH1D*)f2->Get("fhDaverage");
  hAv2->SetTitle(strTitle2.Data());
  SetHistoStyle(hAv2,kRed,21,2,2);
  TH1D *hAbsStat2=(TH1D*)hAv2->Clone("hAbsStat2");
  TH1D *hRelStat2=(TH1D*)hAv2->Clone("hRelStat2");
  for(Int_t j=1;j<=hAbsStat2->GetNbinsX();j++){
    hAbsStat2->SetBinContent(j,0);
    hRelStat2->SetBinError(j,hRelStat2->GetBinError(j)/hRelStat2->GetBinContent(j));
    hRelStat2->SetBinContent(j,0);
  }
  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)f2->Get("grTotSyst");
  gr2->SetTitle(strTitle2.Data());
  SetGraphStyle(gr2,kRed,21,2,2);
  TGraphAsymmErrors *grRelSyst2=(TGraphAsymmErrors*)gr2->Clone("grRelSyst2");  
  npoints=gr2->GetN();
  for(Int_t i=0;i<npoints;i++){
    gr2->GetPoint(i,x,y);
    gr2->SetPoint(i,x,0);
    grRelSyst2->SetPoint(i,x,0);
    grRelSyst2->SetPointEYlow(i,grRelSyst2->GetErrorYlow(i)/y);
    grRelSyst2->SetPointEYhigh(i,grRelSyst2->GetErrorYhigh(i)/y);
  }

  TCanvas *cValues=new TCanvas("cValues","cValues",800,800);
  cValues->cd();
  hAv1->Draw();
  hAv2->Draw("same");

  TCanvas *cAbsStat=new TCanvas("cAbsStat","cAbsStat",800,800);
  cAbsStat->cd();
  hAbsStat1->Draw();
  hAbsStat2->Draw("same");

  TCanvas *cRelStat=new TCanvas("cRelStat","cRelStat",800,800);
  cRelStat->cd();
  hRelStat1->Draw();
  hRelStat2->Draw("same");

  TCanvas *cAbsSyst=new TCanvas("cAbsSyst","cAbsSyst",800,800);
  cAbsSyst->cd();
  gr1->Draw("ap");
  gr2->Draw("p");

  TCanvas *cRelSyst=new TCanvas("cRelSyst","cRelSyst",800,800);
  cRelSyst->cd();
  grRelSyst1->Draw("ap");
  grRelSyst2->Draw("p");



}
