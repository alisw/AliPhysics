/* 
 Macro for fit of azimuthal correlation
 in final style
 by: Fabio (fabio.colamaria@cern.ch)
*/

TString inputdirectory = "";
Int_t pPbyear=2013;

//gStyle->SetLineStyleString(9,"80 20");

void SetInputDirectory(TString strdir){
  inputdirectory=strdir;
}

void DoNiceFitPlots() {
  TFile fInPP(Form("%s/FitOutputs_pp/FitSystematics_pp_WeightedAverage_0.3_99.0.root",inputdirectory.Data()));
  TFile fInpp2(Form("%s/FitOutputs_pp/FitSystematics_pp_WeightedAverage_1.0_2.0.root",inputdirectory.Data()));

  TCanvas *cInPP = (TCanvas*)fInPP.Get("cFitting_0");
  TCanvas *cInpp2 = (TCanvas*)fInpp2.Get("cFitting_0");
   gStyle->SetOptStat(0000);
  gStyle->SetOptFit(000);

  Int_t indexpp = 3;
  Int_t indexpp2 = 4;

  TCanvas *cOut1 = ExtractPad(cInPP,3);
  cOut1->SetName("cNew_pp");
  TCanvas *cOut2 = ExtractPad(cInpp2,4);
  cOut2->SetName("cNew_pp2");

  TCanvas *cUnc_pp = (TCanvas*)fInPP.Get(Form("cFinalCorrelation%d",indexpp-1));
  TCanvas *cUnc_pp2 = (TCanvas*)fInpp2.Get(Form("cFinalCorrelation%d",indexpp2-1));

  TGraphAsymmErrors *graph_pp = (TGraphAsymmErrors*)cUnc_pp->FindObject(Form("grapherror%d",indexpp-1));
  TGraphAsymmErrors *graph_pp2 = (TGraphAsymmErrors*)cUnc_pp2->FindObject(Form("grapherror%d",indexpp2-1));

  TPad* pad1 = (TPad*)cOut1->FindObject("pad");
  pad1->cd();
  graph_pp->SetLineColor(kBlack);
  graph_pp->SetMarkerColor(kBlack);
  graph_pp->SetFillStyle(0);
  graph_pp->Draw("E2");

  TPad* pad2 = (TPad*)cOut2->FindObject("pad");
  pad2->cd();
  graph_pp2->SetLineColor(kBlack);
  graph_pp2->SetMarkerColor(kBlack);
  graph_pp2->SetFillStyle(0);
  graph_pp2->Draw("E2");

  ModStyle(cOut1,0,0);
     TLatex *chi = new TLatex(0.18,4.3,"#chi^{2}/ndf = 1.38");
     chi->SetTextFont(43);
     chi->SetTextSize(32);
     chi->Draw();
  cOut1->Draw();
  ModStyle(cOut2,0,1);
  cOut2->Draw();

  SaveCanvas(cOut1,Form("%s/NiceStylePlots",inputdirectory.Data()),"cFitOutput_NiceStyle_pp_WeightedAverage_0.3_99.0");
  SaveCanvas(cOut2,Form("%s/NiceStylePlots",inputdirectory.Data()),"cFitOutput_NiceStyle_pp_WeightedAverage_1.0_2.0");

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

void ModStyle(TCanvas *c, Int_t system, Int_t opt) {
printf("***c list***\n");
c->ls();
  TF1 *fun = (TF1*)c->FindObject("fGausASper");
  fun->SetLineWidth(4);
  fun->SetLineColor(kGreen+3);
  fun->SetLineStyle(5);

  TF1 *fun2 = (TF1*)c->FindObject("fModGausNSperFixBeta");
  fun2->SetLineWidth(4);
  fun2->SetLineStyle(9);

  TF1 *fun3 = (TF1*)c->FindObject("fPed");
  fun3->SetLineWidth(4);
  fun3->SetLineStyle(2);

  TH1D* hSuperimp;

  TH1D *h = (TH1D*)c->FindObject("fHist");
  h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()+0.3);
  h->GetYaxis()->CenterTitle(kTRUE);
  h->GetYaxis()->SetTitleSize(0.046);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{#scale[1.25]{-1}})");
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitleSize(0.046);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.12);
  if(system==1) {
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
    h->SetMarkerStyle(21);
    hSuperimp = (TH1D*)h->Clone();
    hSuperimp->SetMarkerStyle(25);
    hSuperimp->SetMarkerColor(kRed+1);
  }
printf("***h list***\n");
h->GetListOfFunctions()->ls();
  TF1 *funfit = (TF1*)(h->GetListOfFunctions()->FindObject("kModifNSGausPeriodicityFixBeta"));
  funfit->SetLineWidth(4); 
  funfit->SetLineColor(kRed+1); 

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

  TLegend * legend = new TLegend(0.47,0.46,0.73,0.65);
  legend->SetFillColor(0);
  legend->SetMargin(0.33);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->AddEntry(funfit,"Total fit","l");
  legend->AddEntry(fun2,"Near side","l");
  legend->AddEntry(fun,"Away side","l");
  legend->AddEntry(fun3,"Baseline","l");
  legend->Draw("same");

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
        
   /*     TText *t1 = tl->GetLine(0);
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
        l6->Draw(); */                            

	    obj->Clear();
      }
    }
  }

  pad->cd();

  TLatex *tl1=new TLatex(0.215,0.87,Form("Average D^{0}, D^{+}, D^{*+}"));
  tl1->SetNDC();
  tl1->SetTextSize(0.042);
  tl1->SetTextFont(42);
  tl1->Draw("same");

  TLegend * legend2 = new TLegend(0.20,0.795,0.5,0.85);
  legend2->SetFillColor(0);
  legend2->SetMargin(0.3);
  legend2->SetTextSize(0.042);
  legend2->SetBorderSize(0);
  if(system==0) legend2->AddEntry(h,"pp, #sqrt{#it{s}} = 5.02 TeV","lp");
  else legend2->AddEntry(h,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lp");
  legend2->Draw("same");
 
  if(system==1) {
    TLegend * legendSuperimp = new TLegend(0.20,0.795,0.5,0.85);
    legendSuperimp->SetFillStyle(0);
    legendSuperimp->SetMargin(0.3);
    legendSuperimp->SetTextSize(0.042);
    legendSuperimp->SetBorderSize(0);
    legendSuperimp->AddEntry(hSuperimp,"","lp");
    legendSuperimp->Draw("same");
  }


  TLatex *tl2b=new TLatex(0.79,0.87,Form("ALICE"));
  tl2b->SetNDC();
  tl2b->SetTextSize(0.042);
    tl2b->SetTextFont(42);
  tl2b->Draw("same");

  if(system==0) {
    TLatex *tl3=new TLatex(0.215,0.75,Form("|#it{y}^{D}_{cms}| < 0.5, |#Delta#it{#eta}| < 1"));
    tl3->SetNDC();
    tl3->SetTextSize(0.042);
        tl3->SetTextFont(42);
    tl3->Draw("same");
  } else {
    TLatex *tl3=new TLatex(0.215,0.75,Form("-0.96 < #it{y}^{D}_{cms} < 0.04, |#Delta#it{#eta}| < 1"));
    tl3->SetNDC();
    tl3->SetTextSize(0.042);
        tl3->SetTextFont(42);
    tl3->Draw("same");
  }

  if(opt==0) {
    TLatex *tl4=new TLatex(0.215,0.69,Form("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > 0.3 GeV/#it{c}"));
    tl4->SetNDC();
    tl4->SetTextSize(0.042);
        tl4->SetTextFont(42);
    tl4->Draw("same");
  } else {
    TLatex *tl4;
    tl4=new TLatex(0.215,0.69,Form("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}"));
    tl4->SetNDC();
    tl4->SetTextSize(0.042);
        tl4->SetTextFont(42);
    tl4->Draw("same");
  }
/*
  TLatex *tlAlice=new TLatex(0.68,0.83,Form("#bf{ALICE}"));
  tlAlice->SetNDC();
  tlAlice->Draw();
  tlAlice->SetTextSize(0.038);
*/
  if(opt==0) {
    TLatex *tlUnc=new TLatex(0.35,0.18,Form("#pm4%s scale uncertainty","%"));
    tlUnc->SetNDC();
    tlUnc->SetTextSize(0.042);
        tlUnc->SetTextFont(42);
    tlUnc->Draw("same");
  } else {
    TLatex *tlUnc;
    tlUnc=new TLatex(0.35,0.18,Form("#pm4%s scale uncertainty","%"));
    tlUnc->SetNDC();
    tlUnc->SetTextSize(0.042);
       tlUnc->SetTextFont(42);
    tlUnc->Draw("same");
  }

//UNCOMMENT TO WRITE FORMULA AND BASELINE APPROACH INFO
 /*
  TLatex *tl5=new TLatex(0.21,0.58,Form("#bf{Fit function: f(#Delta#varphi) = #it{b} + #frac{#it{A}_{NS}}{#sqrt{2#pi}#sigma_{NS}} exp#left(- #frac{(#Delta#varphi)^{2}}{2#sigma_{NS}^{2}}#right)+#frac{#it{A}_{AS}}{#sqrt{2#pi}#sigma_{AS}} exp#left(- #frac{(#Delta#varphi-#pi)^{2}}{2#sigma_{AS}^{2}}#right)}"));
  tl5->SetNDC();
  tl5->SetTextSize(0.0255);
  tl5->Draw("same");
*/
 /* TLatex *tl6=new TLatex(0.18,0.52,Form("#bf{Baseline (#it{b}) fixed to the weighted average of points in #pi/4 < #Delta#varphi < #pi/2}"));
  tl6->SetNDC();
  tl6->SetTextSize(0.024);
  tl6->Draw("same");*/
  if(system==1) {
    hSuperimp->Draw("same");
    hSuperimp->GetFunction("kModifNSGausPeriodicityFixBeta")->SetBit(TF1::kNotDraw);
  }

  fun->Draw("same");
  fun2->Draw("same");
  fun3->Draw("same");
  funfit->Draw("same");

    h->Draw("same");
    h->GetFunction("kModifNSGausPeriodicityFixBeta")->SetBit(TF1::kNotDraw);

  return;
}

//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString nameoutput){
   
  c->SaveAs(Form("%s/%s.root",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.png",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.eps",directory.Data(),nameoutput.Data()));
  c->SaveAs(Form("%s/%s.pdf",directory.Data(),nameoutput.Data()));
  
}

