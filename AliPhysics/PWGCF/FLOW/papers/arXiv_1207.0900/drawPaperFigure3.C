Bool_t gPreliminary = kFALSE;

TString drawOptions="e0,p,same,X0";

//systematic error
Double_t gV2For30To40 = 0.0914567;
Double_t gV2For30To40Error = 0.0093375;
Double_t gSysErrorForPlusMinusOn3pCorrelator = 3.64e-06; 
Double_t gSysErrorForSameChargeOn3pCorrelator = 3.396e-06; 

static const Int_t numberOfBins = 25;
Double_t gDeltaPt[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gDeltaPtError[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInDeltaPt[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInDeltaPtError[numberOfBins] = {0.,0.,0.,0.,0.};

Double_t gAvgPt[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gAvgPtError[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInAvgPt[numberOfBins] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInAvgPtError[numberOfBins] = {0.,0.,0.,0.,0.};

Double_t gDeltaEta[numberOfBins*2] = {0.,0.,0.,0.,0.};
Double_t gDeltaEtaError[numberOfBins*2] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInDeltaEta[numberOfBins*2] = {0.,0.,0.,0.,0.};
Double_t gCorrelatorInDeltaEtaError[numberOfBins*2] = {0.,0.,0.,0.,0.};

void drawPaperFigure3() {
  //Draw the differential 3--particle correlator for the PRL
  gROOT->LoadMacro("SetFlowStyle.C");
  SetFlowStyle();

  //TGaxis::SetMaxDigits(2);
  int nRebin = 4;
  TFile *fIn = TFile::Open("differentialRP.root");
  //   fIn->ls();

  //Correlator vs pt diff
  TProfile *gHist3pvsPtDiffPlusMinus = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsPtDiffPlusMinus"));
  gHist3pvsPtDiffPlusMinus->Rebin(nRebin);
  //Printf("Bins: %d",gHist3pvsPtDiffPlusMinus->GetNbinsX());
  myTGraphSetUpPr(gHist3pvsPtDiffPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);

  Double_t sumYN = 0., sumN = 0.;
  Double_t dError = 0.;
  for(Int_t iBin = 1; iBin <= gHist3pvsPtDiffPlusMinus->GetNbinsX(); iBin++) {
    sumYN += gHist3pvsPtDiffPlusMinus->GetBinEntries(iBin)*gHist3pvsPtDiffPlusMinus->GetBinContent(iBin);
    sumN += gHist3pvsPtDiffPlusMinus->GetBinEntries(iBin);

    gDeltaPt[iBin-1] = gHist3pvsPtDiffPlusMinus->GetBinCenter(iBin);
    gCorrelatorInDeltaPt[iBin-1] = gHist3pvsPtDiffPlusMinus->GetBinContent(iBin);
    gCorrelatorInDeltaPtError[iBin-1] = (1.55/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForPlusMinusOn3pCorrelator,2) + TMath::Power((gCorrelatorInDeltaPt[iBin-1]*gV2For30To40Error),2));

    dError += gHist3pvsPtDiffPlusMinus->GetBinEntries(iBin)*gHist3pvsPtDiffPlusMinus->GetBinEntries(iBin)*gCorrelatorInDeltaPtError[iBin-1]*gCorrelatorInDeltaPtError[iBin-1];
    Printf("(Same charge) Bin: %d - Delta pt: %lf - Entries: %lf - Value: %lf - Error: %lf",iBin,gHist3pvsPtDiffPlusMinus->GetBinCenter(iBin),gHist3pvsPtDiffPlusMinus->GetBinEntries(iBin),gCorrelatorInDeltaPt[iBin-1],gCorrelatorInDeltaPtError[iBin-1]);
  }
  if(sumN != 0.)
    dError /= TMath::Power(sumN,2);
  dError = TMath::Sqrt(dError);
  cout<<"(PlusMinus) cos: "<<sumYN/sumN<<" - error: "<<dError<<endl;
  TGraphErrors *grHist3pvsPtDiffPlusMinusSystematic = new TGraphErrors(gHist3pvsPtDiffPlusMinus->GetNbinsX(),gDeltaPt,gCorrelatorInDeltaPt,gDeltaPtError,gCorrelatorInDeltaPtError);
  myTGraphSetUp(grHist3pvsPtDiffPlusMinusSystematic,24,mySysErrColorOpp,myMarkerSize,1,mySysErrColorOpp,10,1001,mySysErrColorOpp);

  TProfile *gHist3pvsPtDiffSameCharge = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsPtDiffSameCharge"));
  gHist3pvsPtDiffSameCharge->Rebin(nRebin);
  myTGraphSetUpPr(gHist3pvsPtDiffSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);

  sumYN = 0., sumN = 0.;
  dError = 0.;
  for(Int_t iBin = 1; iBin <= gHist3pvsPtDiffSameCharge->GetNbinsX(); iBin++) {
    sumYN += gHist3pvsPtDiffSameCharge->GetBinEntries(iBin)*gHist3pvsPtDiffSameCharge->GetBinContent(iBin);
    sumN += gHist3pvsPtDiffSameCharge->GetBinEntries(iBin);

    gDeltaPt[iBin-1] = gHist3pvsPtDiffSameCharge->GetBinCenter(iBin);
    gCorrelatorInDeltaPt[iBin-1] = gHist3pvsPtDiffSameCharge->GetBinContent(iBin);
    gCorrelatorInDeltaPtError[iBin-1] = (1.45/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForSameChargeOn3pCorrelator,2) + TMath::Power((gCorrelatorInDeltaPt[iBin-1]*gV2For30To40Error),2));
    //gCorrelatorInDeltaPtError[iBin-1] = 2.*TMath::Abs((gCorrelatorInDeltaPt[iBin-1]/gV2For30To40)*gV2For30To40Error);
    dError += gHist3pvsPtDiffSameCharge->GetBinEntries(iBin)*gHist3pvsPtDiffSameCharge->GetBinEntries(iBin)*gCorrelatorInDeltaPtError[iBin-1]*gCorrelatorInDeltaPtError[iBin-1];
    
    Printf("(Same charge) Bin: %d - Delta pt: %lf - Entries: %lf - Value: %lf - Error: %lf",iBin,gHist3pvsPtDiffSameCharge->GetBinCenter(iBin),gHist3pvsPtDiffSameCharge->GetBinEntries(iBin),gCorrelatorInDeltaPt[iBin-1],gCorrelatorInDeltaPtError[iBin-1]);
  }
  if(sumN != 0.)
    dError /= TMath::Power(sumN,2);
  dError = TMath::Sqrt(dError);
  cout<<"(Same charge) cos: "<<sumYN/sumN<<" - error: "<<dError<<endl;
  TGraphErrors *grHist3pvsPtDiffSameChargeSystematic = new TGraphErrors(gHist3pvsPtDiffSameCharge->GetNbinsX(),gDeltaPt,gCorrelatorInDeltaPt,gDeltaPtError,gCorrelatorInDeltaPtError);
  myTGraphSetUp(grHist3pvsPtDiffSameChargeSystematic,24,mySysErrColorSame,myMarkerSize,1,mySysErrColorSame,10,1001,mySysErrColorSame);
  
  //Correlator vs pt sum
  TProfile *gHist3pvsPtSumPlusMinus = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsPtSumPlusMinus"));
  gHist3pvsPtSumPlusMinus->Rebin(nRebin);
  myTGraphSetUpPr(gHist3pvsPtSumPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  for(Int_t iBin = 1; iBin <= gHist3pvsPtSumPlusMinus->GetNbinsX(); iBin++) {
    gAvgPt[iBin-1] = gHist3pvsPtSumPlusMinus->GetBinCenter(iBin);
    gCorrelatorInAvgPt[iBin-1] = gHist3pvsPtSumPlusMinus->GetBinContent(iBin);
    gCorrelatorInAvgPtError[iBin-1] = (1.55/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForPlusMinusOn3pCorrelator,2) + TMath::Power((gCorrelatorInAvgPt[iBin-1]*gV2For30To40Error),2));
    //gCorrelatorInAvgPtError[iBin-1] = 2.*TMath::Abs((gCorrelatorInAvgPt[iBin-1]/gV2For30To40)*gV2For30To40Error);
    //Printf("Bin: %d - Delta pt: %lf - Value: %lf (%lf) - Error: %lf",iBin,gHist3pvsPtSumPlusMinus->GetBinCenter(iBin),gCorrelatorInAvgPt[iBin-1],1.e+06*gHist3pvsPtSumPlusMinus->GetBinContent(iBin),gCorrelatorInAvgPtError[iBin-1]);
  }
  TGraphErrors *grHist3pvsPtSumPlusMinusSystematic = new TGraphErrors(gHist3pvsPtSumPlusMinus->GetNbinsX(),gAvgPt,gCorrelatorInAvgPt,gAvgPtError,gCorrelatorInAvgPtError);
  myTGraphSetUp(grHist3pvsPtSumPlusMinusSystematic,24,mySysErrColorOpp,myMarkerSize,1,mySysErrColorOpp,10,1001,mySysErrColorOpp);
   TProfile *gHist3pvsPtSumSameCharge = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsPtSumSameCharge"));
  gHist3pvsPtSumSameCharge->Rebin(nRebin);
  myTGraphSetUpPr(gHist3pvsPtSumSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  for(Int_t iBin = 1; iBin <= gHist3pvsPtSumSameCharge->GetNbinsX(); iBin++) {
    gAvgPt[iBin-1] = gHist3pvsPtSumSameCharge->GetBinCenter(iBin);
    gCorrelatorInAvgPt[iBin-1] = gHist3pvsPtSumSameCharge->GetBinContent(iBin);
    gCorrelatorInAvgPtError[iBin-1] = (1.45/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForSameChargeOn3pCorrelator,2) + TMath::Power((gCorrelatorInAvgPt[iBin-1]*gV2For30To40Error),2));
    //gCorrelatorInAvgPtError[iBin-1] = 2.*TMath::Abs((gCorrelatorInAvgPt[iBin-1]/gV2For30To40)*gV2For30To40Error);
    //Printf("Bin: %d - Delta pt: %lf - Value: %lf (%lf) - Error: %lf",iBin,gHist3pvsPtSumSameCharge->GetBinCenter(iBin),1.e+06*gCorrelatorInAvgPt[iBin-1],1.e+06*gHist3pvsPtSumSameCharge->GetBinContent(iBin),1.e+06*gCorrelatorInAvgPtError[iBin-1]);
  }
  TGraphErrors *grHist3pvsPtSumSameChargeSystematic = new TGraphErrors(gHist3pvsPtSumSameCharge->GetNbinsX(),gAvgPt,gCorrelatorInAvgPt,gAvgPtError,gCorrelatorInAvgPtError);
  myTGraphSetUp(grHist3pvsPtSumSameChargeSystematic,24,mySysErrColorSame,myMarkerSize,1,mySysErrColorSame,10,1001,mySysErrColorSame);

  //Correlator vs eta diff
  TProfile *gHist3pvsEtaDiffPlusMinus = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsEtaDiffPlusMinus"));
  gHist3pvsEtaDiffPlusMinus->Rebin(2);
  myTGraphSetUpPr(gHist3pvsEtaDiffPlusMinus,24,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  for(Int_t iBin = 1; iBin <= gHist3pvsEtaDiffPlusMinus->GetNbinsX(); iBin++) {
    gDeltaEta[iBin-1] = gHist3pvsEtaDiffPlusMinus->GetBinCenter(iBin);
    gCorrelatorInDeltaEta[iBin-1] = gHist3pvsEtaDiffPlusMinus->GetBinContent(iBin);
    gCorrelatorInDeltaEtaError[iBin-1] = (1.55/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForPlusMinusOn3pCorrelator,2) + TMath::Power((gCorrelatorInDeltaEta[iBin-1]*gV2For30To40Error),2));
    //gCorrelatorInDeltaEtaError[iBin-1] = 2.*TMath::Abs((gCorrelatorInDeltaEta[iBin-1]/gV2For30To40)*gV2For30To40Error);
    //Printf("Bin: %d - Delta pt: %lf - Value: %lf (%lf) - Error: %lf",iBin,gHist3pvsEtaDiffPlusMinus->GetBinCenter(iBin),gCorrelatorInDeltaEta[iBin-1],1.e+06*gHist3pvsEtaDiffPlusMinus->GetBinContent(iBin),gCorrelatorInDeltaEtaError[iBin-1]);
  }
  TGraphErrors *grHist3pvsEtaDiffPlusMinusSystematic = new TGraphErrors(gHist3pvsEtaDiffPlusMinus->GetNbinsX(),gDeltaEta,gCorrelatorInDeltaEta,gDeltaEtaError,gCorrelatorInDeltaEtaError);
  myTGraphSetUp(grHist3pvsEtaDiffPlusMinusSystematic,24,mySysErrColorOpp,myMarkerSize,1,mySysErrColorOpp,10,1001,mySysErrColorOpp);
   TProfile *gHist3pvsEtaDiffSameCharge = dynamic_cast<TProfile *>(fIn->Get("gHist3pvsEtaDiffSameCharge"));
  gHist3pvsEtaDiffSameCharge->Rebin(2);
  myTGraphSetUpPr(gHist3pvsEtaDiffSameCharge,20,myDarkRed,myMarkerSize,1,myDarkRed,2,1001,myDarkRed);
  for(Int_t iBin = 1; iBin <= gHist3pvsEtaDiffSameCharge->GetNbinsX(); iBin++) {
    gDeltaEta[iBin-1] = gHist3pvsEtaDiffSameCharge->GetBinCenter(iBin);
    gCorrelatorInDeltaEta[iBin-1] = gHist3pvsEtaDiffSameCharge->GetBinContent(iBin);
    gCorrelatorInDeltaEtaError[iBin-1] = (1.45/gV2For30To40)*TMath::Sqrt(TMath::Power(gSysErrorForSameChargeOn3pCorrelator,2) + TMath::Power((gCorrelatorInDeltaEta[iBin-1]*gV2For30To40Error),2));
    //gCorrelatorInDeltaEtaError[iBin-1] = 2.*TMath::Abs((gCorrelatorInDeltaEta[iBin-1]/gV2For30To40)*gV2For30To40Error);
    //Printf("Bin: %d - Delta pt: %lf - Value: %lf (%lf) - Error: %lf",iBin,gHist3pvsEtaDiffSameCharge->GetBinCenter(iBin),1.e+06*gCorrelatorInDeltaEta[iBin-1],1.e+06*gHist3pvsEtaDiffSameCharge->GetBinContent(iBin),1.e+06*gCorrelatorInDeltaEtaError[iBin-1]);
  }
  TGraphErrors *grHist3pvsEtaDiffSameChargeSystematic = new TGraphErrors(gHist3pvsEtaDiffSameCharge->GetNbinsX(),gDeltaEta,gCorrelatorInDeltaEta,gDeltaEtaError,gCorrelatorInDeltaEtaError);
  myTGraphSetUp(grHist3pvsEtaDiffSameChargeSystematic,24,mySysErrColorSame,myMarkerSize,1,mySysErrColorSame,10,1001,mySysErrColorSame);

  //======================================================//
  //gROOT->ForceStyle();
  TF1 *f1 = new TF1("f1","0",0,1000);
  f1->SetLineColor(1); f1->SetLineStyle(1); f1->SetLineWidth(1);

  TCanvas *c1 = new TCanvas("c1","Centrality dependence: nu_dyn",
			    0,0,1500,600);
//   c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);
//   c1->SetGridx(); c1->SetGridy();
  c1->Draw();
  c1->cd();
  //TPad *myPad1 = new TPad("myPad1","myPad1",0,0.0,0.33,1);
  TPad *myPad1 = new TPad("myPad1","myPad1",0,0.0,0.375,1);
  myPadSetUp(myPad1,0.25,0.08,0.00,0.17);
  TPad *myPad2 = new TPad("myPad2","myPad2",0.375,0,0.69,1);
  myPadSetUp(myPad2,0.0,0.08,0.00,0.17);
  TPad *myPad3 = new TPad("myPad3","myPad3",0.69,0,1,1);
  myPadSetUp(myPad3,0.0,0.08,0.04,0.17);
  myPad1->Draw();
  myPad2->Draw();
  myPad3->Draw();
  //======================================================//
  TLatex *captionLatex = new TLatex();
  captionLatex->SetNDC();
  captionLatex->SetTextFont(42);
//   captionLatex->SetTextColor(1);
  captionLatex->SetTextSize(0.08);

  //_____________________________________________________//
  //Draw the results
  //====================================//
  //Results vs centrality
  myPad1->cd();
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			   ";;#LT cos(#phi_{#alpha} + #phi_{#beta} - 2#Psi_{RP}) #GT",
			   1000,0,2.0,1000,-5e-04,3e-04);
//   gEmpty1->SetTitleSize(0.06,"xyz");
  gEmpty1->GetYaxis()->SetLabelSize(0.06);
  gEmpty1->GetYaxis()->SetTitleSize(0.07);
  gEmpty1->GetYaxis()->SetTitleOffset(1.5);
  gEmpty1->GetXaxis()->SetTitleSize(0.07);
  gEmpty1->GetXaxis()->SetLabelSize(0.06);
  gEmpty1->GetXaxis()->SetTitleOffset(1.1);
  gEmpty1->GetXaxis()->SetLabelOffset(0.0);
//   gEmpty1->GetYaxis()->SetFontSize(0.1);
  gEmpty1->SetAxisRange(0.,1.99,"X");
  gEmpty1->SetStats(kFALSE);
  //gEmpty1->GetYaxis()->SetTitleOffset(1.3);
  //gEmpty1->GetXaxis()->SetTitleOffset(1.5);
  gEmpty1->GetYaxis()->SetNdivisions(5);
  gEmpty1->GetXaxis()->SetNdivisions(5);
  gEmpty1->GetXaxis()->SetTitle("|p_{#font[42]{t},#alpha} - p_{#font[42]{t},#beta}| (GeV/#font[42]{c})");
  gEmpty1->Draw();
  captionLatex->DrawLatex(0.33,0.22,"(a)");

  //_____________________________________________________//

  //_____________________________________________________//
//   myPad1->SetTicks(0,0);
  f1->Draw("same");
  grHist3pvsPtDiffPlusMinusSystematic->Draw("Z");
  grHist3pvsPtDiffSameChargeSystematic->Draw("Z");
  gHist3pvsPtDiffPlusMinus->Draw(drawOptions);
  gHist3pvsPtDiffSameCharge->Draw(drawOptions);
  
  //_____________________________________________________//
  myPad2->cd();
  TH2F *gEmpty2 = (TH2F*)gEmpty1->Clone("gEmpty2");
  gEmpty2->GetYaxis()->SetLabelSize(0.1);
  gEmpty2->GetYaxis()->SetTitleSize(0.1);
  gEmpty2->GetYaxis()->SetTitleOffset(0.6);
  gEmpty2->GetXaxis()->SetTitleSize(0.08);
  gEmpty2->GetXaxis()->SetLabelSize(0.07);
  gEmpty2->GetXaxis()->SetTitleOffset(0.95);
  gEmpty2->GetXaxis()->SetLabelOffset(-0.01);
  gEmpty2->SetAxisRange(0.01,1.99,"X");
//   myPad2->SetTicks(1,0);
//   TH2F *gEmpty3 = (TH2F*)gEmpty1->Clone("gEmpty3");
  gEmpty2->GetYaxis()->SetNdivisions(0);
  gEmpty2->GetXaxis()->SetTitle("(p_{#font[42]{t},#alpha} + p_{#font[42]{t},#beta})/2 (GeV/#font[42]{c})");
  gEmpty2->Draw();
  f1->Draw("same");
  grHist3pvsPtSumPlusMinusSystematic->Draw("Z");
  grHist3pvsPtSumSameChargeSystematic->Draw("Z");
  gHist3pvsPtSumPlusMinus->Draw(drawOptions);
  gHist3pvsPtSumSameCharge->Draw(drawOptions);


  TLatex *captionLatex1 = new TLatex();
  captionLatex1->SetNDC();
  captionLatex1->SetTextFont(42);
//   captionLatex->SetTextColor(1);
  captionLatex1->SetTextSize(0.10);

  captionLatex1->DrawLatex(0.1,0.22,"(b)");

  //_____________________________________________________//
  myPad3->cd();
//   myPad3->SetTicks(1,1);
  TH2F *gEmpty3 = (TH2F*)gEmpty1->Clone("gEmpty3");
  gEmpty3->GetYaxis()->SetNdivisions(0);
  gEmpty3->GetXaxis()->SetTitle("#Delta #eta = |#eta_{#alpha} - #eta_{#beta}|");
  gEmpty3->SetAxisRange(0.01,1.6,"X");
  gEmpty3->GetYaxis()->SetLabelSize(0.1);
  gEmpty3->GetYaxis()->SetTitleSize(0.1);
  gEmpty3->GetYaxis()->SetTitleOffset(0.6);
  gEmpty3->GetXaxis()->SetTitleSize(0.08);
  gEmpty3->GetXaxis()->SetLabelSize(0.07);
  gEmpty3->GetXaxis()->SetTitleOffset(0.95);
  gEmpty3->GetXaxis()->SetLabelOffset(-0.01);
  gEmpty3->Draw();
  f1->Draw("same");
  grHist3pvsEtaDiffPlusMinusSystematic->Draw("Z");
  grHist3pvsEtaDiffSameChargeSystematic->Draw("Z");
  gHist3pvsEtaDiffPlusMinus->Draw(drawOptions);
  gHist3pvsEtaDiffSameCharge->Draw(drawOptions);
 
  TLatex *captionLatex2 = new TLatex();
  captionLatex2->SetNDC();
  captionLatex2->SetTextFont(42);
  captionLatex2->SetTextSize(0.11);
  captionLatex2->DrawLatex(0.1,0.22,"(c)");

  TLatex *myLatex = new TLatex();
  myLatex->SetNDC();
//   myLatex->SetTextColor(kRed+2);
  myLatex->SetTextSize(0.07);
  myLatex->SetLineWidth(2);

  TLegend *legend = new TLegend(0.29,0.73,0.7,0.9,"","brNDC");
  myLegendSetUp(legend,0.07);
  legend->AddEntry(gHist3pvsPtDiffSameCharge,"#color[1]{  same}","lp");
  legend->AddEntry(gHist3pvsPtDiffPlusMinus,"#color[1]{  opp.}","lp");
  myPad1->cd();
  legend->Draw();
  myPad2->cd();
  myLatex->DrawLatex(0.025,0.84,"ALICE Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV");
  myPad3->cd();
  myLatex->SetTextSize(0.08);
  myLatex->DrawLatex(0.2,0.84,"centrality 30-40%");
//   TLegend *legend1 = new TLegend(0.4,0.8,1.,0.95,
//                                 "","brNDC");
//   myLegendSetUp(legend1,0.08);
//   legend1->AddEntry(gHist3pvsPtDiffPlusMinus,"opp.","lp");
//   legend1->Draw();


  if(gPreliminary) {    
    TLatex *alice = new TLatex(0.69,0.9,"Preliminary");
    alice->SetNDC();
    alice->SetTextColor(kRed+2);
    alice->SetTextSize(0.055);
    alice->SetLineWidth(2);
    alice->Draw();
    
    TPad *myPadLogo = new TPad("myPadLogo",
			       "Pad for ALICE Logo",0.72,0.72,0.89,0.89);
    //myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !       
    myPadSetUp(myPadLogo,0,0,0,0);
    myPadLogo->Draw();
    myPadLogo->cd();
    TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
    myAliceLogo->Draw();

  }
 c1->SaveAs("figure3.eps");
 c1->SaveAs("figure3.pdf");
 c1->SaveAs("figure3.png");
}
