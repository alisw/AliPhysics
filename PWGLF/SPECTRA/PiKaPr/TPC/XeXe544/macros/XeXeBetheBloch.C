

void XeXeBetheBloch() {

  //-------contains BetheBloch Parameterization and Y-projection of dE/dx 
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

     TFile * inFile = TFile::Open("./Output_2018_03_17data/MergeData.root");        // open analysis results folder
    TList* list = (TList*) inFile->Get("chist");                   // make my output container to store results
 
  /*TH2D * HistPt = (TH2D*) list->FindObject("fHistPt");
    HistPt->GetXaxis()->SetTitle("p/z(GeV/c)");
    HistPt->GetYaxis()->SetTitle("counts");
    TCanvas * canvHistPt = new TCanvas("canvHistPt","canvHistPt");
    HistPt->Draw();

    TH2D * fPions = (TH2D*) list->FindObject("fPions");
    fPions->GetXaxis()->SetTitle("p/z(GeV/c)");
    fPions->GetYaxis()->SetTitle("dE/dx in TPC");
    TCanvas * canvfPions = new TCanvas("canvfPions","canvfPions");
    fPions->Draw("col");
    TLatex *   pion = new TLatex(0.1724, 116.0,"#pi");
    pion->SetTextSize(0.05);
    pion->SetTextFont(42);
    pion->SetLineWidth(2);
    pion->Draw();
    TLatex *   pionbar = new TLatex(-0.28, 116.0,"#bar#pi");
    pionbar->SetTextSize(0.05);
    pionbar->SetTextFont(42);
    pionbar->SetLineWidth(2);
    pionbar->Draw();


    TH2D * fKaons = (TH2D*) list->FindObject("fKaons");
    fKaons->GetXaxis()->SetTitle("p/z(GeV/c)");
    fKaons->GetYaxis()->SetTitle("dE/dx in TPC");
    TCanvas * canvfKaons = new TCanvas("canvfKaons","canvfKaons");
    fKaons->Draw("col");
    TLatex *   kaon = new TLatex(0.42, 140.0,"k");
    kaon->SetTextSize(0.05);
    kaon->SetTextFont(42);
    kaon->SetLineWidth(2);
    kaon->Draw();
    TLatex *   kaonbar = new TLatex(-0.50, 140.0,"#bark");
    kaonbar->SetTextSize(0.05);
    kaonbar->SetTextFont(42);
    kaonbar->SetLineWidth(2);
    kaonbar->Draw();

    TH2D * fProton = (TH2D*) list->FindObject("fProton");
    fProton->GetXaxis()->SetTitle("p/z(GeV/c)");
    fProton->GetYaxis()->SetTitle("dE/dx in TPC");
    TCanvas * canvfProton = new TCanvas("canvfProton","canvfProton");
    fProton->Draw("col");
    TLatex *   pro = new TLatex(0.67, 190.0,"p");
    pro->SetTextSize(0.05);
    pro->SetTextFont(42);
    pro->SetLineWidth(2);
    pro->Draw();
    TLatex *   probar = new TLatex(-0.8, 190.0,"#barp");
    probar->SetTextSize(0.05);
    probar->SetTextFont(42);
    probar->SetLineWidth(2);
    probar->Draw();

    TH2D * fDeuteron = (TH2D*) list->FindObject("fDeuteron");
    fDeuteron->GetXaxis()->SetTitle("p/z(GeV/c)");
    fDeuteron->GetYaxis()->SetTitle("dE/dx in TPC");
    TCanvas * canvfDeuteron = new TCanvas("canvfDeuteron","canvfDeuteron");
    fDeuteron->Draw("col");
    TLatex *   deut = new TLatex(0.90, 300.0,"d");
    deut->SetTextSize(0.05);
    deut->SetTextFont(42);
    deut->SetLineWidth(2);
    deut->Draw();
    TLatex *   deutbar = new TLatex(-1.05, 300.0,"#bard");
    deutbar->SetTextSize(0.05);
    deutbar->SetTextFont(42);
    deutbar->SetLineWidth(2);
    deutbar->Draw(); */
  //-----------------------------------------------------------------------------  
  TH2D * histdEdx = (TH2D*) list->FindObject("fHistdEdxData");
  histdEdx->SetNameTitle("histdEdx","histdEdx");
  
  histdEdx->GetXaxis()->SetRangeUser(-3,3);
  histdEdx->GetYaxis()->SetRangeUser(0.,1000.);
  histdEdx->GetXaxis()->SetTitle("p/z(GeV/c)");
  histdEdx->GetYaxis()->SetTitle("dE/dx in TPC");

  //--------------------------------------------------------------------------------------
  TF1 *foPion = new TF1("foPion","AliExternalTrackParam::BetheBlochAleph(x/0.138,[0],[1],[2],[3],[4])",0.01,100);
  foPion->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foAntiPion = new TF1("foAntiPion","AliExternalTrackParam::BetheBlochAleph(-x/0.138,[0],[1],[2],[3],[4])",-100,0);
  foAntiPion->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foKaon = new TF1("foKaon","AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.01,100);
  foKaon->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foAntiKaon = new TF1("foAntiKaon","AliExternalTrackParam::BetheBlochAleph(-x/0.493677,[0],[1],[2],[3],[4])",-100,0);
  foAntiKaon->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foProton = new TF1("foProton","AliExternalTrackParam::BetheBlochAleph(x/0.938,[0],[1],[2],[3],[4])",0.01,100);
  foProton->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foPbar = new TF1("foPbar","AliExternalTrackParam::BetheBlochAleph(-x/0.938,[0],[1],[2],[3],[4])",-100,0);
  foPbar->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foDeut = new TF1("foDeut","AliExternalTrackParam::BetheBlochAleph(x/1.8756,[0],[1],[2],[3],[4])",0.01,100);
  foDeut->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foAntiDeut = new TF1("foAntiDeut","AliExternalTrackParam::BetheBlochAleph(-x/1.8756,[0],[1],[2],[3],[4])",-100,0);
  foAntiDeut->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foTriton = new TF1("foTriton","AliExternalTrackParam::BetheBlochAleph(x/2.808921,[0],[1],[2],[3],[4])",0.01,100);
  foTriton->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foAntiTriton = new TF1("foAntiTriton","AliExternalTrackParam::BetheBlochAleph(-x/2.808921,[0],[1],[2],[3],[4])",-100,0);
  foAntiTriton->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foHe3 = new TF1("foHe3","4*AliExternalTrackParam::BetheBlochAleph(2*x/2.809413,[0],[1],[2],[3],[4])",0.01,100);
  foHe3->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);
  TF1 *foAntiHe3 = new TF1("foAntiHe3","4*AliExternalTrackParam::BetheBlochAleph(2*(-x)/2.809413,[0],[1],[2],[3],[4])",-100,0);
  foAntiHe3->SetParameters(1.45265,27.4992,4.003129e-15,2.46246,6.78938);

  //

  /*
    Double_t x[1000], y[1000], yErr[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    x[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (x[counter] < 0.3 || x[counter] > 1.0) continue;
    if (x[counter] > 0.3 && x[counter] < 1.0) {
    TH1D * proj = histdEdx->ProjectionY("proj",iBin,iBin);
    TF1 gausFunc("gausFunc","gaus",foProton->Eval(x[counter])*0.7, foProton->Eval(x[counter])*1.3);
    proj->Fit(&gausFunc,"QNR");
    y[counter]    = gausFunc.GetParameter(1);
    yErr[counter] = gausFunc.GetParError(1);
    delete proj;
    counter++;
    }
    } 
    TGraphErrors * grProton = new TGraphErrors(counter, x, y, 0,yErr); 
    //---------------------------------------------
 
    Double_t xPbar[1000], yPbar[1000], yErrPbar[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xPbar[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xPbar[counter] > -0.3 || xPbar[counter] < -1.0) continue;
    TH1D * projPbar = histdEdx->ProjectionY("projPbar",iBin,iBin);
    TF1 gausFuncPbar("gausFuncPbar","gaus",foPbar->Eval(xPbar[counter])*0.7, foPbar->Eval(xPbar[counter])*1.3);
    projPbar->Fit(&gausFuncPbar,"QNR");
    yPbar[counter]    = gausFuncPbar.GetParameter(1);
    yErrPbar[counter] = gausFuncPbar.GetParError(1);
    delete projPbar;
    counter++;
    }
  
    TGraphErrors * grPbar = new TGraphErrors(counter, xPbar, yPbar, 0,yErrPbar); 
    //-------------------
    Double_t xPion[1000], yPion[1000], yErrPion[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xPion[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xPion[counter] < 0.05 || xPion[counter] > 1.5) continue;
    TH1D * projPion = histdEdx->ProjectionY("projPion",iBin,iBin);
    TF1 gausFuncPion("gausFuncPion","gaus",foPion->Eval(xPion[counter])*0.7, foPion->Eval(xPion[counter])*1.3);
    projPion->Fit(&gausFuncPion,"QNR");
    yPion[counter]    = gausFuncPion.GetParameter(1);
    yErrPion[counter] = gausFuncPion.GetParError(1);
    delete projPion;
    counter++;
    }
    TGraphErrors * grPionFit = new TGraphErrors(counter, xPion, yPion, 0,yErrPion);
    //---------------------------
  
    Double_t xDeut[1000], yDeut[1000], yErrDeut[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xDeut[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xDeut[counter] < 0.6 || xDeut[counter] > 1.25) continue;
    TH1D * projDeut = histdEdx->ProjectionY("projDeut",iBin,iBin);
    TF1 gausFuncDeut("gausFuncDeut","gaus",foDeut->Eval(xDeut[counter])*0.7, foDeut->Eval(xDeut[counter])*1.3);
    projDeut->Fit(&gausFuncDeut,"QNR");
    yDeut[counter]    = gausFuncDeut.GetParameter(1);
    yErrDeut[counter] = gausFuncDeut.GetParError(1);
    delete projDeut;
    counter++;
    }
    TGraphErrors * grDeut = new TGraphErrors(counter, xDeut, yDeut, 0,yErrDeut); 
    //----------------------
    //
    Double_t xTriton[1000], yTriton[1000], yErrTriton[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xTriton[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xTriton[counter] < 0.8 || xTriton[counter] > 1.0) continue;
    TH1D * projTriton = histdEdx->ProjectionY("projTriton",iBin,iBin);
    TF1 gausFuncTriton("gausFuncTriton","gaus",foTriton->Eval(xTriton[counter])*0.7, foTriton->Eval(xTriton[counter])*1.3);
    projTriton->Fit(&gausFuncTriton,"QNR");
    yTriton[counter]    = gausFuncTriton.GetParameter(1);
    yErrTriton[counter] = gausFuncTriton.GetParError(1);
    delete projTriton;
    counter++;
    }
    TGraphErrors * grTriton = new TGraphErrors(counter, xTriton, yTriton, 0,yErrTriton);
    //----------------------
    //TF1 *fHe3 = new TF1("fHe3","4.0*AlephNew_He((2*x)/(2.8084))", 0.1,100);

 
  
    Double_t xHe3[1000], yHe3[1000], yErrHe3[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xHe3[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xHe3[counter] < 1.5 || xHe3[counter] > 2.0) continue;
    TH1D * projHe3 = histdEdx->ProjectionY("projHe3",iBin,iBin);
    TF1 gausFuncHe3("gausFuncHe3","gaus",foHe3->Eval(xHe3[counter])*0.7, foHe3->Eval(xHe3[counter])*1.3);
    projHe3->Fit(&gausFuncHe3,"QNR");
    yHe3[counter]    = gausFuncHe3.GetParameter(1);
    yErrHe3[counter] = gausFuncHe3.GetParError(1);
    delete projHe3;
    counter++;
    }
    TGraphErrors * grHe3 = new TGraphErrors(counter, xHe3, yHe3, 0,yErrHe3); 
    //----------------------
 
    Double_t xKaon[1000], yKaon[1000], yErrKaon[1000];
    Int_t counter = 0;
    for(Int_t iBin = 0; iBin < histdEdx->GetXaxis()->GetNbins(); iBin++) {
    xKaon[counter] = histdEdx->GetXaxis()->GetBinCenter(iBin);
    if (xKaon[counter] < 0.15 || xKaon[counter] > 0.55) continue;
    TH1D * projKaon = histdEdx->ProjectionY("projKaon",iBin,iBin);
    TF1 gausFuncKaon("gausFuncKaon","gaus",foKaon->Eval(xKaon[counter])*0.7, foKaon->Eval(xKaon[counter])*1.3);
    projKaon->Fit(&gausFuncKaon,"QNR");
    yKaon[counter]    = gausFuncKaon.GetParameter(1);
    yErrKaon[counter] = gausFuncKaon.GetParError(1);
    delete projKaon;
    counter++;
    }
    TGraphErrors * grKaon = new TGraphErrors(counter, xKaon, yKaon, 0,yErrKaon); */
  //----------------------
  TCanvas * canvdEdx = new TCanvas("canvdEdx","canvdEdx");
  histdEdx->DrawCopy("col");
  histdEdx->SetTitle("TPC PID for Xe-Xe");
  //gPad->SetLogx();
  //gPad->SetTickx(2);
  //gPad->SetTicky(2);
  foProton->Draw("SAME");
  //grProton->Draw("lP");
  //grPbar->Draw("lP");
  foPion->Draw("SAME");
  //grPionFit->Draw("lP");
  foDeut->Draw("SAME");
  //grDeut->Draw("lP");
  foTriton->Draw("SAME");
  //grTriton->Draw("P");
  foHe3->Draw("SAME");
  foAntiHe3->Draw("SAME");
  //grHe3->Draw("lP");
  foPbar->Draw("SAME");
  foAntiPion->Draw("SAME");
  foAntiDeut->Draw("SAME");
  foAntiTriton->Draw("SAME");
  foKaon->Draw("SAME");
  //grKaon->Draw("lp");
  foAntiKaon->Draw("SAME");
  

  //fHe3->Draw("SAME");
  
  TLatex *   pi = new TLatex(0.1724, 116.0,"#pi");
  pi->SetTextSize(0.05);
  pi->SetTextFont(42);
  pi->SetLineWidth(2);
  pi->Draw();

  TLatex *   pibar = new TLatex(-0.28, 116.0,"#bar#pi");
  pibar->SetTextSize(0.05);
  pibar->SetTextFont(42);
  pibar->SetLineWidth(2);
  pibar->Draw();
   
  TLatex *   k = new TLatex(0.42, 140.0,"k");
  k->SetTextSize(0.05);
  k->SetTextFont(42);
  k->SetLineWidth(2);
  k->Draw();

  TLatex *   kbar = new TLatex(-0.50, 140.0,"#bark");
  kbar->SetTextSize(0.05);
  kbar->SetTextFont(42);
  kbar->SetLineWidth(2);
  kbar->Draw();

  TLatex *   p = new TLatex(0.67, 190.0,"p");
  p->SetTextSize(0.05);
  p->SetTextFont(42);
  p->SetLineWidth(2);
  p->Draw();

  TLatex *   pbar = new TLatex(-0.8, 190.0,"#barp");
  pbar->SetTextSize(0.05);
  pbar->SetTextFont(42);
  pbar->SetLineWidth(2);
  pbar->Draw();

  TLatex *   d = new TLatex(0.90, 300.0,"d");
  d->SetTextSize(0.05);
  d->SetTextFont(42);
  d->SetLineWidth(2);
  d->Draw();

  TLatex *   dbar = new TLatex(-1.05, 300.0,"#bard");
  dbar->SetTextSize(0.05);
  dbar->SetTextFont(42);
  dbar->SetLineWidth(2);
  dbar->Draw();

  TLatex *   t = new TLatex(1.18, 350.0,"t");
  t->SetTextSize(0.05);
  t->SetTextFont(42);
  t->SetLineWidth(2);
  t->Draw();

  TLatex *   tbar = new TLatex(-1.20, 350.0,"#bart");
  tbar->SetTextSize(0.05);
  tbar->SetTextFont(42);
  tbar->SetLineWidth(2);
  tbar->Draw();

  TLatex *  He3 = new TLatex(1.40, 403.0,"He3");
  He3->SetTextSize(0.04);
  He3->SetTextFont(42);
  He3->SetLineWidth(2);
  He3->Draw();

  TLatex *  antiHe3 = new TLatex(-1.65, 403.0,"#barHe3");
  antiHe3->SetTextSize(0.04);
  antiHe3->SetTextFont(42);
  antiHe3->SetLineWidth(2);
  antiHe3->Draw();
  //--------------------------------------

  canvdEdx->SaveAs("FinalSpectra_forQM/canvdEdx.png");
  canvdEdx->SaveAs("FinalSpectra_forQM/canvdEdx.root");

  //---------------------------------------
 

}
