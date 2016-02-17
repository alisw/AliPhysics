{
  gROOT->LoadMacro("KMCDetector.cxx+");
  gROOT->LoadMacro("PrepSummaryKMC.C+");
  //
  TObjArray* NewSt100Perf       = PrepSummaryKMC("NewStEff100NoVtx.root",0,"NewStEff100Rerf");
  TObjArray* NewSt100Cor5       = PrepSummaryKMC("NewStEff100NoVtx.root",1,"NewStEff100Cor5");
  TObjArray* NewSt100Cor5hs2in  = PrepSummaryKMC("NewStEff100NoVtx.root",2,"NewStEff100Cor5hs2in");
  TObjArray* NewSt100Cor5hd2in  = PrepSummaryKMC("NewStEff100NoVtx.root",3,"NewStEff100Cor5hd2in");
  TObjArray* NewSt100Cor5hs2inG = PrepSummaryKMC("NewStEff100NoVtx.root",4,"NewStEff100Cor5hs2inG");
  TObjArray* NewSt100Cor5hd2inG = PrepSummaryKMC("NewStEff100NoVtx.root",5,"NewStEff100Cor5hd2inG");
  TObjArray* NewSt100Fake       = PrepSummaryKMC("NewStEff100NoVtx.root",6,"NewStEff100Fake");

  TObjArray* NewSt095Perf       = PrepSummaryKMC("NewStEff095NoVtx.root",0,"NewStEff095Rerf");
  TObjArray* NewSt095Cor5       = PrepSummaryKMC("NewStEff095NoVtx.root",1,"NewStEff095Cor5");
  TObjArray* NewSt095Cor5hs2in  = PrepSummaryKMC("NewStEff095NoVtx.root",2,"NewStEff095Cor5hs2in");
  TObjArray* NewSt095Cor5hd2in  = PrepSummaryKMC("NewStEff095NoVtx.root",3,"NewStEff095Cor5hd2in");
  TObjArray* NewSt095Cor5hs2inG = PrepSummaryKMC("NewStEff095NoVtx.root",4,"NewStEff095Cor5hs2inG");
  TObjArray* NewSt095Cor5hd2inG = PrepSummaryKMC("NewStEff095NoVtx.root",5,"NewStEff095Cor5hd2inG");
  TObjArray* NewSt095Fake       = PrepSummaryKMC("NewStEff095NoVtx.root",6,"NewStEff095Fake");

  //
  TObjArray* NewRS100Perf       = PrepSummaryKMC("NewRSEff100NoVtx.root",0,"NewRSEff100Rerf");
  TObjArray* NewRS100Cor5       = PrepSummaryKMC("NewRSEff100NoVtx.root",1,"NewRSEff100Cor5");
  TObjArray* NewRS100Cor5hs2in  = PrepSummaryKMC("NewRSEff100NoVtx.root",2,"NewRSEff100Cor5hs2in");
  TObjArray* NewRS100Cor5hd2in  = PrepSummaryKMC("NewRSEff100NoVtx.root",3,"NewRSEff100Cor5hd2in");
  TObjArray* NewRS100Cor5hs2inG = PrepSummaryKMC("NewRSEff100NoVtx.root",4,"NewRSEff100Cor5hs2inG");
  TObjArray* NewRS100Cor5hd2inG = PrepSummaryKMC("NewRSEff100NoVtx.root",5,"NewRSEff100Cor5hd2inG");
  TObjArray* NewRS100Fake       = PrepSummaryKMC("NewRSEff100NoVtx.root",6,"NewRSEff100Fake");

  TObjArray* NewRS095Perf       = PrepSummaryKMC("NewRSEff095NoVtx.root",0,"NewRSEff095Rerf");
  TObjArray* NewRS095Cor5       = PrepSummaryKMC("NewRSEff095NoVtx.root",1,"NewRSEff095Cor5");
  TObjArray* NewRS095Cor5hs2in  = PrepSummaryKMC("NewRSEff095NoVtx.root",2,"NewRSEff095Cor5hs2in");
  TObjArray* NewRS095Cor5hd2in  = PrepSummaryKMC("NewRSEff095NoVtx.root",3,"NewRSEff095Cor5hd2in");
  TObjArray* NewRS095Cor5hs2inG = PrepSummaryKMC("NewRSEff095NoVtx.root",4,"NewRSEff095Cor5hs2inG");
  TObjArray* NewRS095Cor5hd2inG = PrepSummaryKMC("NewRSEff095NoVtx.root",5,"NewRSEff095Cor5hd2inG");
  TObjArray* NewRS095Fake       = PrepSummaryKMC("NewRSEff095NoVtx.root",6,"NewRSEff095Fake");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas();
  c1->Draw();
  TH1F* heff = new TH1F("heff","",100,0.1,1.7);
  heff->GetYaxis()->SetNdivisions(521);
  heff->SetXTitle("p_{T}");
  heff->SetYTitle("Rate");  
  heff->SetMinimum(0);
  heff->SetMaximum(1.09);  
  heff->Draw();
  gPad->SetGrid(1,1);
  //
  TGraphErrors* gr;
  //
  Int_t col = 0;
  Int_t mtype = 20;
  Int_t msize = 1;
  TLegend *leg1 = new TLegend(0.20,0.25,0.95,0.65," ");
  leg1->SetNColumns(2);
  leg1->SetFillColor(0);
  TLegendEntry* le=0;
  //

  col = kRed;
  mtype = 20;
  gr = (TGraphErrors*)NewSt100Cor5->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 100%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kRed+2;
  mtype = 21;
  gr = (TGraphErrors*)NewSt100Cor5hs2in->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 100%, min 5h(2 inner)","lp");
  le->SetTextColor(col);
  // 
  //
  //-------------------------

  col = kRed;
  mtype = 24;
  gr = (TGraphErrors*)NewSt095Cor5->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kRed+2;
  mtype = 25;
  gr = (TGraphErrors*)NewSt095Cor5hs2in->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 95%, min 5h(2 inner)","lp");
  le->SetTextColor(col);
  // 
  //
  // 
  col = kRed;
  mtype = 20;
  gr = (TGraphErrors*)NewSt100Fake->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->SetLineStyle(2);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 100%, Fakes","lp");
  le->SetTextColor(col);
  //
  // 
  col = kRed+2;
  mtype = 21;
  gr = (TGraphErrors*)NewSt095Fake->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->SetLineStyle(2);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New St, eff 95%, Fakes","lp");
  le->SetTextColor(col);
  //
  //

  //-------------------------------------------------------------------------------------
  //
  
  col = kBlue;
  mtype = 20;
  gr = (TGraphErrors*)NewRS100Cor5->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 100%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kGreen+2;
  mtype = 21;
  gr = (TGraphErrors*)NewRS100Cor5hd2in->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 100%, min 5h(2 inner pair)","lp");
  le->SetTextColor(col);
  // 
  //
  //-------------------------

  col = kBlue;
  mtype = 24;
  gr = (TGraphErrors*)NewRS095Cor5->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kGreen+2;
  mtype = 25;
  gr = (TGraphErrors*)NewRS095Cor5hd2in->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 95%, min 5h(2 inner pair)","lp");
  le->SetTextColor(col);
  // 
  //
  // 
  col = kBlue;
  mtype = 20;
  gr = (TGraphErrors*)NewRS100Fake->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->SetLineStyle(2);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 100%, Fakes","lp");
  le->SetTextColor(col);
  //
  // 
  col = kGreen+2;
  mtype = 21;
  gr = (TGraphErrors*)NewRS095Fake->At(kEff);
  SetGStyle(gr,col,mtype,msize);
  gr->SetLineStyle(2);
  gr->Draw("pc");
  le = leg1->AddEntry(gr,"New+2in dbl, eff 95%, Fakes","lp");
  le->SetTextColor(col);
  //
  leg1->Draw();
  //
  SaveCanvas(c1,"fig/effFakes_new_st_2indbl","cg");

  //===============================================
  TCanvas* c2 = new TCanvas();
  c2->Draw();
  gPad->Modified();
  gPad->SetLeftMargin(0.15);
  gPad->Modified();
 
  TLegend *leg2 = new TLegend(0.45,0.5,0.9,0.85,"");
  leg2->SetFillColor(0);
  
  TH1F* hsigD = new TH1F("sigd","",100,0.1,1.7);
  hsigD->SetXTitle("p_{T}");
  hsigD->SetYTitle("resolution, #sigma_{r#phi}");  
  hsigD->SetMinimum(0);
  hsigD->SetMaximum(0.025);  
  hsigD->Draw();
  hsigD->GetYaxis()->SetTitleOffset(1.4);
  hsigD->GetYaxis()->SetTitleSize(0.05);
  // 
  gPad->SetGrid(1,1);
  //
  //------------------------------------------------------------

  col = kRed;
  mtype = 20;
  gr = (TGraphErrors*)NewSt095Cor5->At(kSigD);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg2->AddEntry(gr,"New St, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kRed+2;
  mtype = 25;
  gr = (TGraphErrors*)NewSt095Cor5hs2in->At(kSigD);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg2->AddEntry(gr,"New St, eff 95%, min 5h(2 inner)","lp");
  le->SetTextColor(col);
  // 
  //--------------------
  col = kBlue;
  mtype = 20;
  gr = (TGraphErrors*)NewRS095Cor5->At(kSigD);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg2->AddEntry(gr,"New+2in dbl, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kGreen+2;
  mtype = 25;
  gr = (TGraphErrors*)NewRS095Cor5hd2in->At(kSigD);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg2->AddEntry(gr,"New+2in dbl, eff 95%, min 5h(2 inner pair)","lp");
  le->SetTextColor(col);
  // 
  leg2->Draw();
  SaveCanvas(c2,"fig/resD_new_st_2indbl","cg");


  //----------------------------------------------
  TCanvas* c3 = new TCanvas();
  c3->Draw();
  gPad->Modified();
  gPad->SetLeftMargin(0.15);
  gPad->Modified();
 
  TLegend *leg3 = new TLegend(0.45,0.5,0.9,0.85,"");
  leg3->SetFillColor(0);
  
  TH1F* hsigZ = new TH1F("sigd","",100,0.1,1.7);
  hsigZ->SetXTitle("p_{T}");
  hsigZ->SetYTitle("resolution, #sigma_{Z}");  
  hsigZ->SetMinimum(0);
  hsigZ->SetMaximum(0.025);  
  hsigZ->Draw();
  hsigZ->GetYaxis()->SetTitleOffset(1.4);
  hsigZ->GetYaxis()->SetTitleSize(0.05);
  // 
  gPad->SetGrid(1,1);
  //


  col = kRed;
  mtype = 20;
  gr = (TGraphErrors*)NewSt095Cor5->At(kSigZ);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg3->AddEntry(gr,"New St, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kRed+2;
  mtype = 25;
  gr = (TGraphErrors*)NewSt095Cor5hs2in->At(kSigZ);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg3->AddEntry(gr,"New St, eff 95%, min 5h(2 inner)","lp");
  le->SetTextColor(col);
  // 
  //--------------------
  col = kBlue;
  mtype = 20;
  gr = (TGraphErrors*)NewRS095Cor5->At(kSigZ);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg3->AddEntry(gr,"New+2in dbl, eff 95%, min 5h","lp");
  le->SetTextColor(col);
  //
  col = kGreen+2;
  mtype = 25;
  gr = (TGraphErrors*)NewRS095Cor5hd2in->At(kSigZ);
  SetGStyle(gr,col,mtype,msize);
  gr->Draw("pc");
  le = leg3->AddEntry(gr,"New+2in dbl, eff 95%, min 5h(2 inner pair)","lp");
  le->SetTextColor(col);
  // 
  //
  leg3->Draw();
  SaveCanvas(c3,"fig/resZ_new_st_2indbl","cg");

}
