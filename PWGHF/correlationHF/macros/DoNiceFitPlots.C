/* 
 Macro for fit of azimuthal correlation
 in final style
 by: Fabio (fabio.colamaria@cern.ch)
*/

TString inputdirectory = "";

void SetInputDirectory(TString strdir){
  inputdirectory=strdir;
}

void DoNiceFitPlots() {
  TFile fInPP(Form("%s/FitOutputs_pp/FitSystematics_pp_WeightedAverage_1.0_99.0.root",inputdirectory.Data()));
  TFile fInPPb(Form("%s/FitOutputs_pPb/FitSystematics_pPb_WeightedAverage_1.0_99.0.root",inputdirectory.Data()));

  TCanvas *cInPP = (TCanvas*)fInPP.Get("cFitting_0");
  TCanvas *cInPPb = (TCanvas*)fInPPb.Get("cFitting_1");

  Int_t indexpp = 2;
  Int_t indexpPb = 2;

  TCanvas *cOut1 = ExtractPad(cInPP,2);
  ModStyle(cOut1,0);
  cOut1->SetName("cNew_pp");
  cOut1->Draw();
  TCanvas *cOut2 = ExtractPad(cInPPb,2);
  ModStyle(cOut2,1);
  cOut2->SetName("cNew_pPb");
  cOut2->Draw();

  TCanvas *cUnc_pp = (TCanvas*)fInPP.Get(Form("cFinalCorrelation%d",indexpp-1));
  TCanvas *cUnc_pPb = (TCanvas*)fInPPb.Get(Form("cFinalCorrelation%d",indexpPb-1));

  TGraphAsymmErrors *graph_pp = (TGraphAsymmErrors*)cUnc_pp->FindObject(Form("grapherror%d",indexpp-1));
  TGraphAsymmErrors *graph_pPb = (TGraphAsymmErrors*)cUnc_pPb->FindObject(Form("grapherror%d",indexpPb-1));

  TPad* pad1 = (TPad*)cOut1->FindObject("pad");
  pad1->cd();
  graph_pp->SetLineColor(kBlack);
  graph_pp->SetMarkerColor(kBlack);
  graph_pp->SetFillStyle(0);
  graph_pp->Draw("E2");

  TPad* pad2 = (TPad*)cOut2->FindObject("pad");
  pad2->cd();
  graph_pPb->SetLineColor(kBlack);
  graph_pPb->SetMarkerColor(kBlack);
  graph_pPb->SetFillStyle(0);
  graph_pPb->Draw("E2");

  SaveCanvas(cOut1,Form("%s/NiceStylePlots",inputdirectory.Data()),"cFitOutput_NiceStyle_pp_WeightedAverage_1.0_99.0");
  SaveCanvas(cOut2,Form("%s/NiceStylePlots",inputdirectory.Data()),"cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0");

}

TCanvas* ExtractPad(TCanvas *c, Int_t padnum) {
  TPad *pd = (TPad*)c->GetPad(padnum);
  TCanvas *cNew = new TCanvas("cNew",c->GetTitle(),800,800);
  cNew->cd();
  pd->SetName("pad");
  pd->SetPad(0,0,1,1);
  pd->SetLeftMargin(0.13);
  pd->DrawClone();
  return cNew;
}

void ModStyle(TCanvas *c, Int_t system) {
  TH1D *h = (TH1D*)c->FindObject("fHist");
  h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()+0.3);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.12);

  Double_t max = h->GetBinContent(h->GetMaximumBin());
  h->SetMaximum(TMath::Floor(((max*2)+1)));

  TPad* pad = (TPad*)c->FindObject("pad");
  pad->SetTickx();
  pad->SetTicky();
  pad->SetMargin(0.17,0.05,0.12,0.05);

  TList *lc=pad->GetListOfPrimitives();
  Int_t entries=lc->GetEntries();
  Int_t syst=1;
  Int_t nextMeson=0;

  for(Int_t jl=0;jl<entries;jl++){
    TObject *obj=(TObject*)lc->At(jl);
    TString strName=obj->ClassName();
    if(strName.Contains("TLatex")) obj->Clear();
    if(strName.Contains("TPaveText")) {
      TPaveText *tl=(TPaveText*)obj;
      TString name=tl->GetName();
      if(!name.Contains("TPave")) obj->Clear(); 
      if(name.Contains("TPave")) {
        TPaveText *tl=(TPaveText*)lc->At(jl);
        
        TText *t1 = tl->GetLine(0);
        TLatex *l1 = new TLatex(0.65,0.74,Form("#bf{%s}",t1->GetTitle()));
        l1->SetTextSize(0.032);
        l1->Draw();
        TText *t2 = tl->GetLine(1);
        TLatex *l2 = new TLatex(0.65,0.69,Form("#bf{%s}",t2->GetTitle()));
        l2->SetTextSize(0.032);
        l2->Draw();
        TText *t3 = tl->GetLine(2);
        TLatex *l3 = new TLatex(0.65,0.64,Form("#bf{%s}",t3->GetTitle()));
        l3->SetTextSize(0.032);
        l3->Draw();
        TText *t4 = tl->GetLine(3);
        TLatex *l4 = new TLatex(0.65,0.59,Form("#bf{%s}",t4->GetTitle()));
        l4->SetTextSize(0.032);
        l4->Draw();
        TText *t5 = tl->GetLine(4);
        TLatex *l5 = new TLatex(0.65,0.54,Form("#bf{%s}",t5->GetTitle()));
        l5->SetTextSize(0.032);
        l5->Draw();
        TText *t6 = tl->GetLine(5);
        TLatex *l6 = new TLatex(0.65,0.49,Form("#bf{%s}",t6->GetTitle()));
        l6->SetTextSize(0.032);
        l6->Draw();                             

	    obj->Clear();
      }
    }
  }

  pad->cd();

  TLatex *tl1=new TLatex(0.215,0.60,Form("#bf{Average D^{0}, D^{+}, D^{*+}}"));
  tl1->SetNDC();
  tl1->SetTextSize(0.038);
  tl1->Draw("same");

  if(system==0) {
    TLatex *tl2=new TLatex(0.215,0.88,Form("#bf{ALICE, pp, #sqrt{#it{s}}=7 TeV}"));
    tl2->SetNDC();
    tl2->SetTextSize(0.038);
    tl2->Draw("same");
  } else {
    TLatex *tl2=new TLatex(0.215,0.88,Form("#bf{ALICE, p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV}"));
    tl2->SetNDC();
    tl2->SetTextSize(0.038);
    tl2->Draw("same");
  }

  if(system==0) {
    TLatex *tl3=new TLatex(0.215,0.81,Form("#bf{5<#it{p}_{T}^{D}<8 GeV/c, |#it{y}^{D}|<0.5}"));
    tl3->SetNDC();
    tl3->SetTextSize(0.038);
    tl3->Draw("same");
  } else {
    TLatex *tl3=new TLatex(0.215,0.81,Form("#bf{8<#it{p}_{T}^{D}<16 GeV/c, -0.96<#it{y}^{D}_{cms}<0.04}"));
    tl3->SetNDC();
    tl3->SetTextSize(0.038);
    tl3->Draw("same");
  }

  TLatex *tl4=new TLatex(0.215,0.74,Form("#bf{#it{p}_{T}^{assoc}>1 GeV/c, |#Delta#eta|<1.0}"));
  tl4->SetNDC();
  tl4->SetTextSize(0.038);
  tl4->Draw("same");
/*
  TLatex *tlAlice=new TLatex(0.68,0.83,Form("#bf{ALICE}"));
  tlAlice->SetNDC();
  tlAlice->Draw();
  tlAlice->SetTextSize(0.038);
*/
  if(system==0) {
    TLatex *tlUnc=new TLatex(0.20,0.67,Form("#bf{{}^{+13%}_{-10%} scale uncertainty}"));
    tlUnc->SetNDC();
    tlUnc->SetTextSize(0.042);
    tlUnc->Draw("same");
  } else {
    TLatex *tlUnc=new TLatex(0.20,0.67,Form("#bf{{}^{+10%}_{-10%} scale uncertainty}"));
    tlUnc->SetNDC();
    tlUnc->SetTextSize(0.042);
    tlUnc->Draw("same");
  }

//UNCOMMENT TO WRITE FORMULA AND BASELINE APPROACH INFO
/* 
  TLatex *tl5=new TLatex(0.18,0.6,Form("#bf{Fit function: f(#Delta#varphi) = #it{b} + #frac{Y_{NS}}{#sqrt{2#pi}#sigma_{NS}} exp#left(- #frac{(#Delta#varphi)^{2}}{2#sigma_{NS}^{2}}#right)+#frac{Y_{AS}}{#sqrt{2#pi}#sigma_{AS}} exp#left(- #frac{(#Delta#varphi-#pi)^{2}}{2#sigma_{AS}^{2}}#right)}"));
  tl5->SetNDC();
  tl5->SetTextSize(0.024);
  tl5->Draw("same");

  TLatex *tl6=new TLatex(0.18,0.52,Form("#bf{Baseline (#it{b}) fixed to the weighted average of points in #pi/4 < #Delta#varphi < #pi/2}"));
  tl6->SetNDC();
  tl6->SetTextSize(0.024);
  tl6->Draw("same");
*/
  return;
}

//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString nameoutput){
   
  c->SaveAs(Form("%s/%s.root",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.png",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.eps",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.pdf",directory.Data(),nameoutput.Data()));
  
}


