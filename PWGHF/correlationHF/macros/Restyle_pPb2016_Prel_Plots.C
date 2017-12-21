void Restyle_pPb2016_Prel_Plots() {

  gSystem->Exec("mkdir PrelPlots_SQM2016");

  /*********************************/
  // pPb/MC fit observables (NS 1)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/ComparisonPPbtoMC/ComparePPbtoMCFitResults.root");
  TCanvas *c = (TCanvas*)fInput.Get("cPPbvsMCFitResultsFinalPaperStyle");

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  TLatex *obj0=(TLatex*)lc->At(21); //shift the pT>0.3 label
  TLatex *obj0a=(TLatex*)lc->At(20); //ALICe preliminary
  obj0a->SetY(obj0a->GetY()+0.001);
  obj0a->SetTextSize(28);
  obj0->SetX(obj0->GetX()-0.16);
  TLatex *obj0a=(TLatex*)lc->At(19); //ALICE preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetY(obj0a->GetY()-0.003); 
  obj0a->SetTextSize(obj0a->GetTextSize()-3);
  obj0a->SetTitle("ALICE Preliminary");

  TPad *p2 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_2");
  TList *lc=p2->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(21);  //remove fake data legend
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(20);  //remove fake data legend
  obj2->Clear();
  obj2->SetLineColor(kWhite);

  TPad *p6 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_6");
  TList *lc=p6->GetListOfPrimitives();
  lc->ls();

  TH1D *hData = (TH1D*)lc->At(3);
  TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(15);
  TGraphAsymmErrors *hV2 = (TGraphAsymmErrors*)lc->At(17);

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *leg = new TLegend(0.25,0.54,0.55,0.68); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  leg->AddEntry(hV2,"Total data syst unc","f");
  leg->AddEntry(hBox,"Syst unc from v_{2#Delta} hypothesis","f");
  p1->cd();
  leg->DrawClone();
  c->cd();

  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults.root");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults.eps");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults.pdf");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults.png");
  
 
  /*********************************/
  // pPb/MC fit observables (NS 2)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/ComparisonPPbtoMC/ComparePPbtoMCFitResults_2.root");
  TCanvas *c = (TCanvas*)fInput.Get("cPPbvsMCFitResultsFinalPaperStyle");

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  TLatex *obj0=(TLatex*)lc->At(21); //shift the pT>0.3 label
  TLatex *obj0a=(TLatex*)lc->At(20); //ALICe preliminary
  obj0a->SetY(obj0a->GetY()+0.001);
  obj0a->SetTextSize(28);
  TLatex *obj0a=(TLatex*)lc->At(19); //ALICe preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetY(obj0a->GetY()-0.003); 
  obj0a->SetTextSize(obj0a->GetTextSize()-3);
  obj0a->SetTitle("ALICE Preliminary");

  TPad *p2 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_2");
  TList *lc=p2->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(21);  //remove fake data legend
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(20);  //remove fake data legend
  obj2->Clear();
  obj2->SetLineColor(kWhite);

  TPad *p3 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_3");
  TList *lc=p3->GetListOfPrimitives();
  lc->ls();

  TH1D *hData = (TH1D*)lc->At(2);
  TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors *hV2 = (TGraphAsymmErrors*)lc->At(18);
  hData->SetBinContent(hData->GetXaxis()->FindBin(4),0.);
  hData->SetBinError(hData->GetXaxis()->FindBin(4),0.);
  hBox->RemovePoint(0);
  hV2->RemovePoint(0);
  TH1D *hData = (TH1D*)lc->At(17);  
  hData->SetBinContent(hData->GetXaxis()->FindBin(4),0.);
  hData->SetBinError(hData->GetXaxis()->FindBin(4),0.);  
  TH1D *hData = (TH1D*)lc->At(19);  
  hData->SetBinContent(hData->GetXaxis()->FindBin(4),0.);
  hData->SetBinError(hData->GetXaxis()->FindBin(4),0.);

  TPad *p6 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_6");
  TList *lc=p6->GetListOfPrimitives();
  lc->ls();

  TH1D *hData = (TH1D*)lc->At(2);
  TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(15);
  TGraphAsymmErrors *hV2 = (TGraphAsymmErrors*)lc->At(17);
  hData->SetBinContent(hData->GetXaxis()->FindBin(4),0.);
  hData->SetBinError(hData->GetXaxis()->FindBin(4),0.);
  hBox->RemovePoint(0);
  hV2->RemovePoint(0);

  TFrame *fr = (TFrame*)lc->At(0);
  fr->SetX2(28);
  hData->GetXaxis()->SetRangeUser(0,28);

  TH1D *hData2 = (TH1D*)lc->At(18);
  hData2->SetBinContent(hData2->GetXaxis()->FindBin(4),0.);
  hData2->SetBinError(hData2->GetXaxis()->FindBin(4),0.);

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *leg = new TLegend(0.25,0.54,0.55,0.68); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  leg->AddEntry(hV2,"Total data syst unc","f");
  leg->AddEntry(hBox,"Syst unc from v_{2#Delta} hypothesis","f");
  p1->cd();
  leg->DrawClone();

  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults_2.root");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults_2.eps");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults_2.pdf");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResults_2.png");


  /*********************************/
  // pPb/MC fit observables (AS)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/ComparisonPPbtoMC/ComparePPbtoMCFitResultsAS.root");
  TCanvas *c = (TCanvas*)fInput.Get("cPPbvsMCFitResultsFinalPaperStyleAS");

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_1");
  TList *lc=p1->GetListOfPrimitives();

  TLatex *obj0=(TLatex*)lc->At(21); //shift the pT>0.3 label
  TLatex *obj0a=(TLatex*)lc->At(20); //ALICe preliminary
  obj0a->SetY(obj0a->GetY()+0.001);  
  obj0a->SetTextSize(28);
  TLatex *obj0a=(TLatex*)lc->At(19); //ALICe preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetY(obj0a->GetY()-0.003); 
  obj0a->SetTextSize(obj0a->GetTextSize()-3);
  obj0a->SetTitle("ALICE Preliminary");

  TPad *p2 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_2");
  TList *lc=p2->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(21);  //remove fake data legend
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(20);  //remove fake data legend
  obj2->Clear();
  obj2->SetLineColor(kWhite);

  TPad *p6 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_6");
  TList *lc=p6->GetListOfPrimitives();
  lc->ls();

  TFrame *fr = (TFrame*)lc->At(0);
  fr->SetX2(28);
  hData->GetXaxis()->SetRangeUser(0,28);

  TH1D *hData2 = (TH1D*)lc->At(18);
  hData2->SetBinContent(hData2->GetXaxis()->FindBin(4),0.);
  hData2->SetBinError(hData2->GetXaxis()->FindBin(4),0.);

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *leg = new TLegend(0.25,0.54,0.55,0.68); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  leg->AddEntry(hV2,"Total data syst unc","f");
  leg->AddEntry(hBox,"Syst unc from v_{2#Delta} hypothesis","f");
  p1->cd();
  leg->DrawClone();

  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS.root");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS.eps");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS.pdf");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS.png");


  /*********************************/
  // pPb/MC fit observables (AS 2)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/ComparisonPPbtoMC/ComparePPbtoMCFitResultsAS_2.root");
  TCanvas *c = (TCanvas*)fInput.Get("cPPbvsMCFitResultsFinalPaperStyleAS");

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_1");
  TList *lc=p1->GetListOfPrimitives();
  TLatex *obj0=(TLatex*)lc->At(21); //shift the pT>0.3 label
  TLatex *obj0a=(TLatex*)lc->At(20); //ALICe preliminary
  obj0a->SetY(obj0a->GetY()+0.001);  
  obj0a->SetTextSize(28);
  TLatex *obj0a=(TLatex*)lc->At(19); //ALICe preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetY(obj0a->GetY()-0.003); 
  obj0a->SetTextSize(obj0a->GetTextSize()-3);
  obj0a->SetTitle("ALICE Preliminary");

  TPad *p2 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_2");
  TList *lc=p2->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(21);  //remove fake data legend
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(20);  //remove fake data legend
  obj2->Clear();
  obj2->SetLineColor(kWhite);

  TPad *p6 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_6");
  TList *lc=p6->GetListOfPrimitives();
  lc->ls();

  TFrame *fr = (TFrame*)lc->At(0);
  fr->SetX2(28);
  hData->GetXaxis()->SetRangeUser(0,28);

  TH1D *hData2 = (TH1D*)lc->At(18);
  hData2->SetBinContent(hData2->GetXaxis()->FindBin(4),0.);
  hData2->SetBinError(hData2->GetXaxis()->FindBin(4),0.);

  TPad *p1 = (TPad*)c->FindObject("cPPbvsMCFitResultsFinalPaperStyleAS_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *leg = new TLegend(0.25,0.54,0.55,0.68); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  leg->AddEntry(hV2,"Total data syst unc","f");
  leg->AddEntry(hBox,"Syst unc from v_{2#Delta} hypothesis","f");
  p1->cd();
  leg->DrawClone();

  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS_2.root");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS_2.eps");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS_2.pdf");
  c->SaveAs("PrelPlots_SQM2016/ComparePPbtoMCFitResultsAS_2.png");


  /*********************************/
  // pPb/pp deltaPhi correls  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_UniqueCanvas_Style1.root");
  TCanvas *c = (TCanvas*)fInput.Get("cFinalPaperStyle");
 
  for(Int_t i=1; i<=9; i++){
    TPad *p = (TPad*)c->FindObject(Form("cFinalPaperStyle_%d",i));
    TList *lc=p->GetListOfPrimitives();

  if(i<7) {
    TH1D *hData = (TH1D*)lc->At(1);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.8);
    TH1D *hData2 = (TH1D*)lc->At(5);
    hData2->SetLineWidth(1);
    hData2->SetMarkerSize(0.8);
    TH1D *hData3 = (TH1D*)lc->At(7);
    hData3->SetLineWidth(1);
    hData3->SetMarkerSize(0.8);
    TH1D *hData4 = (TH1D*)lc->At(10);
    hData4->SetLineWidth(1);
    hData4->SetMarkerSize(1);

    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(2);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(3);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(6);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(9);
    hBox->SetLineWidth(1);
  } else {
    TH1D *hData = (TH1D*)lc->At(0);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.8);
    TH1D *hData2 = (TH1D*)lc->At(4);
    hData2->SetLineWidth(1);
    hData2->SetMarkerSize(0.8);
    TH1D *hData3 = (TH1D*)lc->At(6);
    hData3->SetLineWidth(1);
    hData3->SetMarkerSize(0.8);
    TH1D *hData4 = (TH1D*)lc->At(9);
    hData4->SetLineWidth(1);
    hData4->SetMarkerSize(1);

    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(1);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(2);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(5);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(8);
    hBox->SetLineWidth(1);
  }
  }


  TPad *p1 = (TPad*)c->FindObject("cFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TPaveText *obj0b=(TPaveText*)lc->At(17); //move the average label
  obj0b->SetX1NDC(0.40);
  obj0b->SetY1NDC(obj0b->GetY1NDC()-0.10); //0.03
  TPaveText *obj1a=(TPaveText*)lc->At(16); //ALICE Preliminary
  obj1a->SetX1NDC(0.25);
  obj1a->SetY1NDC(obj1a->GetY1NDC()-0.12);
  obj1a->Clear();
  //obj1a->AddText("ALICE Preliminary"); //SHYAM
  TLegend *lg = (TLegend*)lc->At(12);   //move legend
  lg->SetY1NDC(lg->GetY1NDC()-0.02);
  lg->SetY2NDC(lg->GetY2NDC()-0.04);
  lg->SetX2NDC(lg->GetX2NDC()-0.14);
  lg->SetX1NDC(lg->GetX1NDC()-0.005);

    TH1D *hData1 = (TH1D*)lc->At(1);
    TH1D *hData2 = (TH1D*)lc->At(5);
    TH1D *hData3 = (TH1D*)lc->At(7);
    TH1D *hData4 = (TH1D*)lc->At(10);
    lg->Clear();
    lg->AddEntry(hData4,"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
    lg->AddEntry((TObject*)0,"EPJC 77 (2017) 245","");
    lg->AddEntry(hData1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96 < #it{y}^{D}_{cms} < 0.04","lep");

  TLegend *lg2 = (TLegend*)lc->At(13);
  lg2->SetY1NDC(lg2->GetY1NDC()-0.02);
  lg2->SetY2NDC(lg2->GetY2NDC()-0.04);
  lg2->SetX2NDC(lg2->GetX2NDC()-0.14);
  lg2->SetX1NDC(lg2->GetX1NDC()-0.005);

    lg2->Clear();
    lg2->AddEntry(hData4,"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
    lg2->AddEntry((TObject*)0,"EPJC 77 (2017) 245","");
    lg2->AddEntry(hData3,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96 < #it{y}^{D}_{cms} < 0.04","lep");

  TPaveText *obj0c=(TPaveText*)lc->At(14); //move pt label
  obj0c->SetX1NDC(obj0c->GetX1NDC()-0.01);
  obj0c->SetY2NDC(obj0c->GetY2NDC()-0.03);
  TPaveText *obj0c=(TPaveText*)lc->At(15); //move deltaeta label
  obj0c->SetX1NDC(obj0c->GetX1NDC()-0.14);
  obj0c->SetY2NDC(obj0c->GetY2NDC()-0.02);
  
 //SHYAM -->
  //===Title Size on X and Y axis======================= 
  TH1D* h1 = (TH1D*)lc->At(1);
  h1->GetYaxis()->SetTitleSize(18.5);
  h1->GetYaxis()->SetTitleOffset(3.5);
  TPad *p = (TPad*)c->GetPad(4);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h2 = (TH1D*)lc->At(1);
  h2->GetYaxis()->SetTitleSize(18.5);
  h2->GetYaxis()->SetTitleOffset(3.5);
  TPad *p = (TPad*)c->GetPad(7);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h3 = (TH1D*)lc->At(0);
  h3->GetYaxis()->SetTitleSize(18.5);
  h3->GetYaxis()->SetTitleOffset(3.5);
  h3->GetXaxis()->SetTitleSize(22.0);
  h3->GetXaxis()->SetTitleOffset(2.57);
  TPad *p = (TPad*)c->GetPad(8);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h4 = (TH1D*)lc->At(0);
  h4->GetXaxis()->SetTitleSize(22.0);
  h4->GetXaxis()->SetTitleOffset(2.57);
  TPad *p = (TPad*)c->GetPad(9);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h5 = (TH1D*)lc->At(0);
  h5->GetXaxis()->SetTitleSize(22.0);
  h5->GetXaxis()->SetTitleOffset(2.57);
  
  //===========PaveText ALICE Preliminary on Pad3===================
  TPad * p1 = (TPad*)c->GetPad(3);
  TPaveText* t1=new TPaveText(0.22,0.62,0.60,0.70,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.068);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"ALICE Preliminary");
  t1->SetFillColor(kBlue);
  p1->cd();
  t1->Draw();
  //================D meson - charged hadron correlation on Pad1===============
  TPad * p1 = (TPad*)c->GetPad(1);
  TPaveText* t1=new TPaveText(0.27,0.790,0.60,0.88,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.052);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"D meson-charged hadron correlation");
  t1->SetFillColor(kBlue);
  p1->cd();
  t1->Draw();
//<--SHYAM


  TPad *p3 = (TPad*)c->FindObject("cFinalPaperStyle_3");
  TList *lc=p3->GetListOfPrimitives();
  TLatex *obj0a=(TLatex*)lc->At(8); //ALICE preliminary
  obj0a->SetY(obj0a->GetY()-0.18); 
  TLatex *obj0a=(TLatex*)lc->At(11); //ALICE preliminary
  obj0a->SetY(obj0a->GetY()-0.18); 

  c->SaveAs("PrelPlots_SQM2016/plotComparison_WeightedAverage_pp_pPb_UniqueCanvas_Style1.root");
  c->SaveAs("PrelPlots_SQM2016/plotComparison_WeightedAverage_pp_pPb_UniqueCanvas_Style1.eps");
  c->SaveAs("PrelPlots_SQM2016/plotComparison_WeightedAverage_pp_pPb_UniqueCanvas_Style1.pdf");
  c->SaveAs("PrelPlots_SQM2016/plotComparison_WeightedAverage_pp_pPb_UniqueCanvas_Style1.png");


  /*********************************/
  // pPb/pp Fit Observables
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/ComparisonPPtoPPb/ComparePPtoPPbFitResults.root");
  TCanvas *c = (TCanvas*)fInput.Get("cPPvsPPbFitResultsFinalPaperStyle");
 
  TPad *p1 = (TPad*)c->FindObject("cPPvsPPbFitResultsFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();
  TLatex *obj0a=(TLatex*)lc->At(8); //ALICE preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetY(obj0a->GetY()-0.003); 
  obj0a->SetTextSize(obj0a->GetTextSize()-2);
  obj0a->SetTitle("ALICE Preliminary");

  TPad *p2 = (TPad*)c->FindObject("cPPvsPPbFitResultsFinalPaperStyle_2");
  TList *lc=p2->GetListOfPrimitives();
  lc->ls();
  TH1D *hData = (TH1D*)lc->At(2);
  TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(3);
  TGraphAsymmErrors *hBox2 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors *hV2 = (TGraphAsymmErrors*)lc->At(5);
  TLegend *lg = (TLegend*)lc->At(8);   //move legend
  lg->AddEntry(hBox,"Total data syst unc (pp)","f");
  lg->AddEntry(hBox2,"Total data syst unc (p-Pb)","f");
  lg->AddEntry(hV2,"Syst unc from v_{2#Delta} hypothesis","f");
  lg->SetY1NDC(0.32);
  TLegend *lg2 = (TLegend*)lc->At(9);
  lg2->AddEntry(hBox,"Total data syst unc (pp)","f");
  lg2->AddEntry(hBox2,"Total data syst unc (p-Pb)","f");
  lg2->AddEntry(hV2,"Syst unc from v_{2#Delta} hypothesis","f");
  lg2->SetY1NDC(0.32);

  TPad *p1 = (TPad*)c->FindObject("cPPvsPPbFitResultsFinalPaperStyle_1"); //Add DeltaEta<1
  p1->cd();
  TList *lc=p1->GetListOfPrimitives();
  TLatex *obj0y=(TLatex*)lc->At(10); 
  obj0y->SetX(obj0y->GetX()-0.03);

  TLatex *obj0z=((TLatex*)lc->At(10))->Clone("DeltaEta"); 
  obj0z->SetY(obj0z->GetY()-0.08);
  obj0z->SetX(obj0z->GetX()-0.01);
  obj0z->SetTitle("|#Delta#eta| < 1");
  obj0z->Draw();

  c->SaveAs("PrelPlots_SQM2016/ComparePPtoPPbFitResults.root");
  c->SaveAs("PrelPlots_SQM2016/ComparePPtoPPbFitResults.eps");
  c->SaveAs("PrelPlots_SQM2016/ComparePPtoPPbFitResults.pdf");
  c->SaveAs("PrelPlots_SQM2016/ComparePPtoPPbFitResults.png");


  /*********************************/
  // pPb/MC deltaPhi correls (1)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonToModels/CorrelationppMC4x6_1New.root");
  TCanvas *c = (TCanvas*)fInput.Get("cFinalPaperStyle");
  TPaveText *Pt[12], *Pt1[12];

  for(Int_t i=1; i<=12; i++){
    TPad *p = (TPad*)c->FindObject(Form("cFinalPaperStyle_%d",i));
    TList *lc=p->GetListOfPrimitives();

//-->SHYAM
    if (i==1){ // Shyam Start Changes
      TPaveText *obj1a=(TPaveText*)lc->At(15); //ALICE Preliminary
      obj1a->Clear();
    }
    if (i==1){
      Pt[i]=(TPaveText*)lc->At(16);
      Pt[i]->SetX1(Pt[i]->GetX1()+0.26);
      Pt[i]->SetY1NDC(Pt[i]->GetY1NDC()-0.02);
      Pt[i]->SetTextSize(Pt[i]->GetTextSize()+1.0);
    }
    else
    {
      Pt[i]=(TPaveText*)lc->At(15);
      Pt[i]->SetX1(Pt[i]->GetX1()+0.03);
      Pt[i]->SetY1NDC(Pt[i]->GetY1NDC()-0.02);
      Pt[i]->SetTextSize(Pt[i]->GetTextSize()+1.0);
      if (i==8) Pt[i]->SetX1(Pt[i]->GetX1()-0.03);
      if (i==5) Pt[i]->SetX1(Pt[i]->GetX1()+0.24);
      else if (i==9) Pt[i]->SetX1(Pt[i]->GetX1()+0.26);
    }
//<--SHYAM

    TH1D *hData = (TH1D*)lc->At(1);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TH1D *hData = (TH1D*)lc->At(4);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TH1D *hData = (TH1D*)lc->At(12);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(2);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(5);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(13);
    hBox->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(7);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(8);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(9);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(10);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(11);
    hMC->SetLineWidth(1);
  }

  //Moving legends of simulations (arossi)...
  TPad *p2 = (TPad*)c->FindObject("cFinalPaperStyle_2");
  TList *lc=p2->GetListOfPrimitives();
  lc->ls();  

  TLegend *lgNew = (TLegend*)lc->At(16)->Clone("newLegMC_1");   //move legend
  TLegend *lg = (TLegend*)lc->At(16)->Clear();

  TPad *p3 = (TPad*)c->FindObject("cFinalPaperStyle_3");
  TList *lc=p3->GetListOfPrimitives();
  lc->ls();  

  TLegend *lgNew2 = (TLegend*)lc->At(16)->Clone("newLegMC_2");   //move legend
  TLegend *lg = (TLegend*)lc->At(16)->Clear();  

  p3->cd();
  lgNew->SetX1NDC(lgNew->GetX1NDC()+0.06);
  lgNew->Draw();
  TPad *p4 = (TPad*)c->FindObject("cFinalPaperStyle_4");
  p4->cd();
  lgNew2->SetX1NDC(lgNew2->GetX1NDC()+0.06);
  lgNew2->Draw();
  //...done

  TPad *p1 = (TPad*)c->FindObject("cFinalPaperStyle_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *lg = (TLegend*)lc->At(20);   //move legend
  lg->SetX1NDC(0.336);

//-->SHYAM
  TPaveText *obj2u=(TPaveText*)(lc->At(18)); //Kine info shifted
  obj2u->SetX1(0.32); 
  obj2u->SetTextSize(0.063);
  obj2u->SetTextFont(42);  
  obj2u->Clear();
  obj2u->AddText("-0.96 < #it{y}^{D}_{cms} < 0.04, |#Delta#eta| < 1");

  TPad * p1 = (TPad*)c->GetPad(2);
  TPaveText* t1=new TPaveText(0.27,0.62,0.60,0.80,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.07);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"ALICE Preliminary");
  t1->SetFillColor(kBlue);
  p1->cd();
  t1->Draw();
  //======Changing Y axis TitleSize===========================
  TPad *p = (TPad*)c->GetPad(1);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h1 = (TH1D*)lc->At(1);
  h1->GetYaxis()->SetTitleSize(18.5);
  h1->GetYaxis()->SetTitleOffset(4.77);
  TPad *p = (TPad*)c->GetPad(5);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h2 = (TH1D*)lc->At(1);
  h2->GetYaxis()->SetTitleSize(18.5);
  h2->GetYaxis()->SetTitleOffset(4.77);
  TPad *p = (TPad*)c->GetPad(9);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h3 = (TH1D*)lc->At(1);
  h3->GetYaxis()->SetTitleSize(18.5);
  h3->GetYaxis()->SetTitleOffset(4.77);
 //====Removing p-Pb and Energy on pad 1 (modification) and set Title "D meson-charged hadron correlation"
  TPad *p = (TPad*)c->GetPad(1);
  TList *lc=p->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(20);  
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(14);  
  obj2->Clear();
  obj2->SetLineColor(kWhite);
  TPaveText *obj0b=(TPaveText*)lc->At(17); 
  obj0b->SetX1NDC(0.40);
  obj0b->SetY1NDC(obj0b->GetY1NDC()-0.20); 
  TPaveText *obj0c=(TPaveText*)lc->At(18); 
  obj0c->SetX1NDC(0.35);
  obj0c->SetY1NDC(obj0c->GetY1NDC()-0.25);
  
  TPaveText* t1=new TPaveText(0.33,0.68,0.65,0.78,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.052);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"D meson-charged hadron correlation");
  t1->SetFillColor(kBlue);
  p->cd();
  t1->Draw(); 
//====Adding p-Pb and Energy on pad 5 and 3 Baseline Legend"
  TH1D *hData = (TH1D*)lc->At(4);
  TPad *p = (TPad*)c->GetPad(2);
  TLegend *leg = new TLegend(0.13,0.53,0.70,0.60); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.053);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  p->cd();
  leg->DrawClone();
  
  TH1D *hData = (TH1D*)lc->At(2);
  TPad *p = (TPad*)c->GetPad(2);
  TLegend *leg = new TLegend(0.13,0.45,0.60,0.52); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.053);
  leg->AddEntry(hData,"baseline-subtraction uncertainty","f");
  p->cd();
  leg->Draw();
//<--SHYAM

  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_1New.root");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_1New.eps");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_1New.pdf");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_1New.png");


  /*********************************/
  // pPb/MC deltaPhi correls (2)  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonToModels/CorrelationppMC4x6_2New.root");
  TCanvas *c = (TCanvas*)fInput.Get("cFinalPaperStyle2");
  TPaveText *Pt[12], *Pt1[12];

  for(Int_t i=1; i<=12; i++){
    TPad *p = (TPad*)c->FindObject(Form("cFinalPaperStyle2_%d",i));
    TList *lc=p->GetListOfPrimitives();

//-->SHYAM
    if (i==1){ //Added by shyam
      TPaveText *obj2a=(TPaveText*)lc->At(16); 
      obj2a->Clear();
      TPaveText *obj3a=(TPaveText*)lc->At(19); 
      obj3a->SetTextSize(obj3a->GetTextSize()+3.0);
    }  
    if (i==1){
      Pt1[i]=(TPaveText*)lc->At(14);
      Pt1[i]->SetX1(Pt1[i]->GetX1()+0.27);
      Pt1[i]->SetY1NDC(Pt1[i]->GetY1NDC()-0.02);
      Pt1[i]->SetTextSize(Pt1[i]->GetTextSize()+1.0);
    }
    else
    {
      Pt[i]=(TPaveText*)lc->At(14);
      Pt[i]->SetX1(Pt[i]->GetX1()+0.03);
      Pt[i]->SetY1NDC(Pt[i]->GetY1NDC()-0.02);
      Pt[i]->SetTextSize(Pt[i]->GetTextSize()+1.0);
      if (i==5) Pt[i]->SetX1(Pt[i]->GetX1()+0.25);
      else if (i==9) Pt[i]->SetX1(Pt[i]->GetX1()+0.27);
    }
//<--SHYAM

    TH1D *hData = (TH1D*)lc->At(1);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TH1D *hData = (TH1D*)lc->At(4);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TH1D *hData = (TH1D*)lc->At(12);
    hData->SetLineWidth(1);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerStyle(20);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(2);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(5);
    hBox->SetLineWidth(1);
    TGraphAsymmErrors *hBox = (TGraphAsymmErrors*)lc->At(13);
    hBox->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(7);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(8);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(9);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(10);
    hMC->SetLineWidth(1);
    TH1D *hMC = (TH1D*)lc->At(11);
    hMC->SetLineWidth(1);
  }

  //Moving legends of simulations (arossi)...
  TPad *p2 = (TPad*)c->FindObject("cFinalPaperStyle2_2");
  TList *lc=p2->GetListOfPrimitives();
  lc->ls();  

  TLegend *lgNew = (TLegend*)lc->At(16)->Clone("newLegMC_1");   //move legend
  TLegend *lg = (TLegend*)lc->At(16)->Clear();

  TPad *p3 = (TPad*)c->FindObject("cFinalPaperStyle2_3");
  TList *lc=p3->GetListOfPrimitives();
  lc->ls();  

  TLegend *lgNew2 = (TLegend*)lc->At(16)->Clone("newLegMC_2");   //move legend
  TLegend *lg = (TLegend*)lc->At(16)->Clear();  

  p3->cd();
  lgNew->SetX1NDC(lgNew->GetX1NDC()+0.06);
  lgNew->Draw();
  TPad *p4 = (TPad*)c->FindObject("cFinalPaperStyle2_4");
  p4->cd();
  lgNew2->SetX1NDC(lgNew2->GetX1NDC()+0.06);
  lgNew2->Draw();
  TList *lc=p4->GetListOfPrimitives();
  TLatex *objx=(TLatex*)(lc->At(6)); //Kine info shifted
  objx->SetY(objx->GetY()-0.08); 
  objx->SetY(objx->GetY()-0.08); 
  //...done

  TPad *p1 = (TPad*)c->FindObject("cFinalPaperStyle2_1");
  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLegend *lg = (TLegend*)lc->At(20);   //move legend
  lg->SetX1NDC(0.336);

//-->SHYAM
  TPaveText *obj2u=(TPaveText*)(lc->At(19)); //Kine info shifted
  obj2u->SetTextSize(0.063);
  obj2u->SetTextFont(42);  
  obj2u->Clear();
  obj2u->AddText("-0.96 < #it{y}^{D}_{cms} < 0.04, |#Delta#eta| < 1");

  TPad * p1 = (TPad*)c->GetPad(2);
  
  TPaveText* t1=new TPaveText(0.27,0.62,0.60,0.80,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.07);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"ALICE Preliminary");
  t1->SetFillColor(kBlue);
  p1->cd();
  t1->Draw();
  //======Changing Y axis TitleSize===========================
  TPad *p = (TPad*)c->GetPad(1);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h1 = (TH1D*)lc->At(1);
  h1->GetYaxis()->SetTitleSize(18.5);
  h1->GetYaxis()->SetTitleOffset(4.77);
  TPad *p = (TPad*)c->GetPad(5);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h2 = (TH1D*)lc->At(1);
  h2->GetYaxis()->SetTitleSize(18.5);
  h2->GetYaxis()->SetTitleOffset(4.77);
  TPad *p = (TPad*)c->GetPad(9);
  TList *lc=p->GetListOfPrimitives();
  TH1D* h3 = (TH1D*)lc->At(1);
  h3->GetYaxis()->SetTitleSize(18.5);
  h3->GetYaxis()->SetTitleOffset(4.77);
   //====Removing p-Pb and Energy on pad 1 (modification) and set Title "D meson-charged hadron correlation"
  TPad *p = (TPad*)c->GetPad(1);
  TList *lc=p->GetListOfPrimitives();
  TLegend *obj1=(TLegend*)lc->At(20);  
  obj1->Clear();
  obj1->SetLineColor(kWhite);
  TLegend *obj2=(TLegend*)lc->At(15);  
  obj2->Clear();
  obj2->SetLineColor(kWhite);
  TPaveText *obj0b=(TPaveText*)lc->At(18); 
  obj0b->SetX1NDC(0.40);
  obj0b->SetY1NDC(obj0b->GetY1NDC()-0.20); 
  TPaveText *obj0c=(TPaveText*)lc->At(19); 
  obj0c->SetX1NDC(0.35);
  obj0c->SetY1NDC(obj0c->GetY1NDC()-0.25);
  
  TPaveText* t1=new TPaveText(0.33,0.68,0.65,0.78,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(0.052);
 //t1->SetTextColor(kMagenta);
  t1->SetTextFont(42);
  t1->SetTextAlign(13);
  t1->AddText(0.,0.,"D meson-charged hadron correlation");
  t1->SetFillColor(kBlue);
  p->cd();
  t1->Draw();
//====Adding p-Pb and Energy on pad 5 and 3 Baseline Legend"
  TH1D *hData = (TH1D*)lc->At(4);
  TPad *p = (TPad*)c->GetPad(2);
  TLegend *leg = new TLegend(0.13,0.53,0.70,0.60); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.053);
  leg->AddEntry(hData,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
  p->cd();
  leg->DrawClone();
  
  TH1D *hData = (TH1D*)lc->At(2);
  TPad *p = (TPad*)c->GetPad(2);
  TLegend *leg = new TLegend(0.13,0.45,0.60,0.52); //new data points legend
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.053);
  leg->AddEntry(hData,"baseline-subtraction uncertainty","f");
  p->cd();
  leg->Draw();
//<--SHYAM

  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_2New.root");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_2New.eps");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_2New.pdf");
  c->SaveAs("PrelPlots_SQM2016/CorrelationppMC4x6_2New.png");


  /*********************************/
  // Fit example in nice style  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults/NiceStylePlots/cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0.root");
  TCanvas *c = (TCanvas*)fInput.Get("cNew_pPb");
  
  TPad *p1 = (TPad*)c->FindObject("pad");

  TList *lc=p1->GetListOfPrimitives();
  lc->ls();

  TLatex *obj0a=(TLatex*)(lc->At(13)); //ALICE preliminary
  obj0a->SetX(obj0a->GetX()-0.2); 
  obj0a->SetTitle("ALICE Preliminary");

  c->SaveAs("PrelPlots_SQM2016/cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0.root");
  c->SaveAs("PrelPlots_SQM2016/cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0.eps");
  c->SaveAs("PrelPlots_SQM2016/cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0.pdf");
  c->SaveAs("PrelPlots_SQM2016/cFitOutput_NiceStyle_pPb_WeightedAverage_1.0_99.0.png");


  /*********************************/
  // Average example in nice style  
  /*********************************/

  TFile fInput("$PWD/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0to99.0.root");
  TCanvas *c = (TCanvas*)fInput.Get("cDraw");

  TList *lc=c->GetListOfPrimitives();
  lc->ls();

  TLatex *objL=new TLatex(0.18,0.835,"ALICE Preliminary");
  objL->SetNDC();
  objL->SetTextFont(42);
  objL->SetTextSize(0.040);
  objL->Draw();
  TLatex *obj2aa=(TLatex*)(lc->At(3)); //ALICE preliminary
  obj2aa->SetY(obj2aa->GetY()-0.03); 
  obj2aa->SetTextSize(0.035);
  obj2aa->SetTextFont(42);  
  obj2aa->SetTitle("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  TLatex *objL1=new TLatex(0.18,0.71,"#it{L}_{int} = 300 #mub^{-1}");
  objL1->SetNDC();
  objL1->SetTextFont(42);
  objL1->SetTextSize(0.035);
  objL1->Draw();
  TLatex *obj2bb=(TLatex*)(lc->At(7)); //ALICE preliminary
  obj2bb->SetY(obj2bb->GetY()-0.03); 
  obj2bb->SetTextSize(0.035);
  TLatex *obj2c=(TLatex*)(lc->At(4));
  obj2c->SetTitle("#bf{5 < #it{p}_{T}^{D} < 8 GeV/#it{c}}");
  TLatex *obj2d=(TLatex*)(lc->At(5));
  obj2d->SetTitle("#bf{#it{p}_{T}^{assoc} > 1 GeV/#it{c}}"); 
  TLatex *obj2e=(TLatex*)(lc->At(6));
  obj2e->SetTitle("#bf{|#Delta#eta| < 1}");

  TH1D* h = (TH1D*)lc->At(1);
  h->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
  h->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetYaxis()->CenterTitle(kTRUE);

  c->SaveAs("PrelPlots_SQM2016/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0_99.0.root");
  c->SaveAs("PrelPlots_SQM2016/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0_99.0.eps");
  c->SaveAs("PrelPlots_SQM2016/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0_99.0.pdf");
  c->SaveAs("PrelPlots_SQM2016/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0_99.0.png");

} 

/*

  for(Int_t jl=0;jl<lc->GetEntries();jl++){
    TObject *obj=(TObject*)lc->At(jl);
    TString strName=obj->ClassName();
    printf("%d) %s\n",jl,strName.Data());
  }

*/
/*
    if(strName.Contains("TLatex")) {
      TLatex *tl=(TLatex*)obj;
      TString str=tl->GetTitle();
      str.ReplaceAll("#bf","#font[42]");
      if(str.Contains("#it{p}_{T}^{D^{0}}")) {tl->SetTitle(""); continue;}
   // if(jl==20) obj->Clear();
*/

