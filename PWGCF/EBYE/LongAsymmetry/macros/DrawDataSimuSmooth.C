TH1F *hMeanSlopeReg12WithStat, *hMeanSlopeReg12WithSys;
TH1F * h_c1Reg12_Fromfpol3par1_SimWithSys, *h_c1Reg12_Fromfpol3par1_SimWithStat;
TH1F *hy0Simu_Stat, *hy0Simu_Sys;
TH1F *hy0Data_Stat, *hy0Data_Sys;
TH1F *hy0SimuInput, *hy0SimuInputSys;
TH1F * hy0Actual;
TH1F * h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat, *h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys;
Float_t SystematicsSimuC1[8]={17.8,11.2,10.3,13.6,8.7,12.7,14.3,9.3};
// From y0 Distribution for alpha_ZDC_DATA
//Float_t y0Actual[8]={0.0140,0.0267, 0.0368, 0.0453, 0.0540, 0.0630, 0.0725, 0.0833};
//Float_t dy0Actual[8]={0.000026, 0.000047, 0.000063, 0.000078, 0.000093, 0.000109, 0.000125, 0.000145};
// From y0 Distribution for alpha_ZDC_Sim ( Glauber)
//Float_t y0Actual[8]={0.0152,0.0276, 0.0371, 0.0456, 0.0542, 0.0632, 0.0725, 0.0833};
//Float_t dy0Actual[8]={0.000028, 0.000048, 0.000064, 0.000079, 0.000093, 0.000109, 0.000125, 0.000145};
// From hlaf log A by B
Float_t y0Actual[8]={0.0162,0.0292, 0.0390, 0.0480, 0.0571, 0.0666, 0.0769, 0.0887};
Float_t dy0Actual[8]={0.000028, 0.000048, 0.000064, 0.000079, 0.000093, 0.000109, 0.000125, 0.000145};

TH1F *hws_c1Reg12_Fromfpol3par1_SimWithStat, *hws_c1Reg12_Fromfpol3par1_SimWithSys;
TH1F *h1g_c1Reg12_Fromfpol3par1_SimWithStat, *h1g_c1Reg12_Fromfpol3par1_SimWithSys;

Float_t SimuReg1fpol3_Par1[8], SimuReg1fpol3_dPar1[8], SimuReg2fpol3_Par1[8], SimuReg2fpol3_dPar1[8];
Float_t SimuReg12fpol3_Par1[8],SysErrReg12[8];
Float_t SimuReg1fpol3_Par2[8], SimuReg1fpol3_dPar2[8], SimuReg2fpol3_Par2[8], SimuReg2fpol3_dPar2[8];

TGraph * tg20SlopeData, *tg20y0Data, *tc1Simu, *ty0Simu, *ty0SimuInput;
//Float_t fpol3High[8], fpol3Low[8];
Float_t cent[8]={2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5};
void DrawDataSimuSmooth(){

  SetStyle();
  // Getting c1 Data points and their errors 
  TFile * fdata = new TFile("/Users/rashmiraniwala/ToALICE/ANALYSIS/dNdEtaRatios/DataSlopeWithStatSys.root");
  fdata->cd();
  hMeanSlopeReg12WithStat = (TH1F*)fdata->Get("hMeanSlopeReg12WithStat");
  hc1High = (TH1F*)fdata->Get("hc1High");
  hc1Low = (TH1F*)fdata->Get("hc1Low");
  // Tweeking errors for smooth fit
  Float_t mean = hMeanSlopeReg12WithStat->GetBinContent(6);
  Float_t err = mean*0.210;
  hc1High->SetBinContent(6,mean+err);

  mean = hMeanSlopeReg12WithStat->GetBinContent(6);
  err = mean*0.070;
  hc1Low->SetBinContent(6,mean-err);

  mean = hMeanSlopeReg12WithStat->GetBinContent(2);
  err = mean*0.065;
  hc1Low->SetBinContent(2,mean-err);

  

  tg20SlopeData =GetSysPlot(hc1High, hc1Low);

  // To simulated results now. 
  TFile * fsimu = new TFile("/Users/rashmiraniwala/ToALICE/NumericalSim/Sept2015/GaussWithb/WiderEta/SimuSlopeWithStatSys.root");
  fsimu->cd();
  h_c1Reg12_Fromfpol3par1_SimWithSys= (TH1F*)fsimu->Get("h_c1Reg12_Fromfpol3par1_SimWithSys");
  h_c1Reg12_Fromfpol3par1_SimWithStat= (TH1F*)fsimu->Get("h_c1Reg12_Fromfpol3par1_SimWithStat");
  /*
  TH1F * h_c1SimuHigh = (TH1F*)h_c1Reg12_Fromfpol3par1_SimWithSys->Clone("h_c1SimuHigh");
  TH1F * h_c1SimuLow = (TH1F*)h_c1Reg12_Fromfpol3par1_SimWithSys->Clone("h_c1SimuLow");
  h_c1SimuHigh->Reset();
  h_c1SimuLow->Reset();
  for(Int_t icent = 0;icent<8;icent++){
    Float_t c1high = h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinContent(icent+1)+h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinError(icent+1);
    h_c1SimuHigh->SetBinContent(icent+1, c1high);
    Float_t c1low = h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinContent(icent+1)-h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinError(icent+1);
    h_c1SimuLow->SetBinContent(icent+1, c1low);
    cout<<"Simu c1 = "<<h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinContent(icent+1)<<" +/- "<< h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinError(icent+1)<<" ("<< h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinError(icent+1)*100/h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinContent(icent+1) <<"%)"<<endl;
  }
  */
  // Tweeking c1 Simu errors for smooth fit
  Float_t mean = h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinContent(7);
  Float_t err = mean*0.175;
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetBinError(7,err);
  
  tc1Simu = GetSysPlot1(h_c1Reg12_Fromfpol3par1_SimWithSys);

  hy0Data_Stat = new TH1F("hy0Data_Stat","hy0Data_Stat",8,0.0,40.0);
  hy0Data_High = new TH1F("hy0Data_High","hy0Data_High",8,0.0,40.0);
  hy0Data_Low = new TH1F("hy0Data_Low","hy0Data_Low",8,0.0,40.0);
  
  Float_t y0, slope, sys, stat, y0high, syslow, dy0stat, dy0sys;  
  cout<<" Getting y0 for data"<<endl;
  Float_t EtaSigmaData[8]={3.85, 3.89, 3.86, 3.86, 3.89,3.89, 3.93, 3.93 };
  Float_t dEtaSigmaData[8]={0.22,0.22, 0.24, 0.24, 0.26, 0.26, 0.27, 0.27};
  Float_t syshigh, syslow, dy0syshigh, dy0syslow;
  for(Int_t icent = 0;icent<8;icent++){
    slope = hMeanSlopeReg12WithStat->GetBinContent(icent+1);
    stat = hMeanSlopeReg12WithStat->GetBinError(icent+1);
    //    sys = hMeanSlopeReg12WithSys->GetBinError(icent+1);
    syshigh = hc1High->GetBinContent(icent+1)-hMeanSlopeReg12WithStat->GetBinContent(icent+1);
    syslow = hc1Low->GetBinContent(icent+1)-hMeanSlopeReg12WithStat->GetBinContent(icent+1);
    y0 = slope*EtaSigmaData[icent]*EtaSigmaData[icent];
    dy0stat = (stat/slope)*(stat/slope)+2*(dEtaSigmaData[icent]/EtaSigmaData[icent])*(dEtaSigmaData[icent]/EtaSigmaData[icent]);
    dy0stat = y0*TMath::Sqrt(dy0stat);
    dy0syshigh = (syshigh/slope)*(syshigh/slope) + (0.05)*(0.05);
    dy0syshigh = y0*TMath::Sqrt(dy0syshigh);
    dy0syslow = (syslow/slope)*(syslow/slope) + (0.05)*(0.05);
    dy0syslow = y0*TMath::Sqrt(dy0syslow);
    cout<<"Data "<<icent<<" "<<y0<<" "<<dy0stat*100/y0<<" + "<<dy0syshigh*100/y0<<" - "<<dy0syslow*100/y0<<endl;
    hy0Data_Stat->SetBinContent(icent+1,y0);
    hy0Data_Stat->SetBinError(icent+1,dy0stat);
    hy0Data_High->SetBinContent(icent+1,y0+dy0syshigh);
    hy0Data_Low->SetBinContent(icent+1,y0-dy0syslow);
  }
  
  tg20y0Data = GetSysPlot(hy0Data_High, hy0Data_Low);
  
  hy0Simu_Stat = new TH1F("hy0Simu_Stat","hy0Simu_Stat",8,0.0,40.0);
  hy0Simu_Sys = new TH1F("hy0Simu_Sys","hy0Simu_Sys",8,0.0,40.0);
  
  hy0Actual = new TH1F("hy0Actual","hy0Actual",8,0.0,40.0);
  for(Int_t icent=0;icent<8;icent++){
    hy0Actual->SetBinContent(icent+1,y0Actual[icent]);
    hy0Actual->SetBinError(icent+1,dy0Actual[icent]);
  }
  
  GetMeanStatSimu1g();
  GetMeanStatSimuWS();
  
  
  //  Float_t SystematicSimuC1[8]={17.8,10.3,10.2,13.7,7.7,12.7,14.3,9.2};
  //  Float_t ySigma[8]={3.75,3.79,3.76,3.76,3.79,3.79,3.83, 3.83};
  Float_t EtaSigma[8]={3.836,3.879,3.844,3.844,3.876,3.877,3.92, 3.922};
  Float_t dEtaSigma[8]={0.000723,0.00074,0.00074,0.00074,0.00075,0.00075,0.00078, 0.00078};
  cout<<" Getting y0 for simulated data"<<endl;
  for(Int_t icent = 0;icent<8;icent++){
    slope = h_c1Reg12_Fromfpol3par1_SimWithStat->GetBinContent(icent+1);
    stat = h_c1Reg12_Fromfpol3par1_SimWithStat->GetBinError(icent+1);
    sys = h_c1Reg12_Fromfpol3par1_SimWithSys->GetBinError(icent+1);
    // sys = slope*SystematicsSimuC1[icent]/100.0;
    y0 = slope*EtaSigma[icent]*EtaSigma[icent];
    dy0stat = (stat/slope)*(stat/slope);
    dy0stat+=dEtaSigma[icent]/EtaSigma[icent]+dEtaSigma[icent]/EtaSigma[icent];
    dy0stat = y0*TMath::Sqrt(dy0stat);
    dy0sys = (sys/slope)*(sys/slope) + (0.05)*(0.05);
    //  dy0sys = dy0sys+
    dy0sys = y0*TMath::Sqrt(dy0sys);
    cout<<"Simu "<<icent<<" "<<y0<<" "<<dy0stat*100/y0<<" "<<dy0sys*100/y0<<endl;
    hy0Simu_Stat->SetBinContent(icent+1,y0);
    hy0Simu_Stat->SetBinError(icent+1,dy0stat);
    hy0Simu_Sys->SetBinContent(icent+1,y0);
    hy0Simu_Sys->SetBinError(icent+1,dy0sys);
  }
  ty0Simu = GetSysPlot1(hy0Simu_Sys);
  
  Gety0SimuInput();
  DrawSlope();
  Drawy0();
}

void Gety0SimuInput(){
  TH1F * hy0InReg1, *hy0InReg2;
  TFile * fin = new TFile("/Users/rashmiraniwala/ToALICE/tglauber/fSimuInputy0.root");
  cout<<" Reading file "<<fin->GetName()<<" for input y0"<<endl;
  hy0InReg1 = (TH1F*)fin->Get("hmeany0DataReg1");
  hy0InReg2 = (TH1F*)fin->Get("hmeany0DataReg2");
  hy0SimuInput = (TH1F*)fin->Get("hmeany0DataReg12WithStat");
  hy0SimuInputSys = (TH1F*)fin->Get("hmeany0DataReg12WithSys");
  ty0SimuInput = GetSysPlot1(hy0SimuInputSys);
  /*
  hy0SimuInput = new TH1F("hy0SimuInput","hy0SimuInput",8,0.0,40.0);
  for(Int_t icent =0;icent<8;icent++){
    Float_t y0Reg1 = hy0InReg1->GetBinContent(icent+1);
    Float_t y0Reg2 = hy0InReg2->GetBinContent(icent+1);
    Float_t dy0Reg1 = hy0InReg1->GetBinError(icent+1);
    Float_t dy0Reg2 = hy0InReg2->GetBinError(icent+1);
    Float_t mean = 0.5*(y0Reg1+y0Reg2);
    Float_t stat = mean*TMath::Sqrt(0.5*(dy0Reg1/y0Reg1)*(dy0Reg1/y0Reg1) + 0.5*(dy0Reg2/y0Reg2)*(dy0Reg2/y0Reg2));
    hy0SimuInput->SetBinContent(icent+1,mean);
    hy0SimuInput->SetBinError(icent+1,stat);
    cout<<"y0 = "<<mean<<" +/- "<<stat<<endl;
  }
  */
}

void Drawy0(){

  hy0Data_Stat->SetMarkerStyle(8);
  hy0Data_Stat->SetMarkerSize(1.2);
  hy0Data_Stat->SetLineWidth(2);
  hy0Data_Stat->SetMarkerColor(kRed+1);
  hy0Data_Stat->SetLineColor(kRed+1);
  /*  hy0Data_Sys->SetMarkerStyle(8);
  hy0Data_Sys->SetMarkerColor(kRed+1);
  hy0Data_Sys->SetFillColor(kRed+1-9);
  hy0Data_Sys->SetFillStyle(3001);
  hy0Data_Sys->SetLineColor(kRed+1);*/
  tg20y0Data->SetFillStyle(3001);
  tg20y0Data->SetFillColor(kRed-9);
  
  hy0SimuInput->SetMarkerStyle(29);
  hy0SimuInput->SetMarkerSize(1.6);
  hy0SimuInput->SetLineWidth(2);
  hy0SimuInput->SetMarkerColor(kCyan+2);
  hy0SimuInput->SetLineColor(kCyan+2);
  hy0SimuInput->SetFillColor(kCyan-8);
  /*hy0SimuInputSys->SetMarkerStyle(29);
  hy0SimuInputSys->SetMarkerColor(kCyan+2);
  hy0SimuInputSys->SetLineColor(kCyan+2);
  hy0SimuInputSys->SetFillColor(kCyan-8);*/
  ty0SimuInput->SetFillColor(kCyan-8);
  ty0SimuInput->SetFillStyle(3001);
  
  hy0Actual->SetMarkerStyle(34);
  hy0Actual->SetMarkerColor(kMagenta+1);
  hy0Actual->SetMarkerSize(1.5);
  hy0Actual->SetLineColor(kMagenta+1);
  
  hy0Simu_Stat->SetMarkerStyle(21);
  hy0Simu_Stat->SetMarkerColor(kBlue+1);
  hy0Simu_Stat->SetLineColor(kBlue+1);
  hy0Simu_Stat->SetLineWidth(2);
  hy0Simu_Sys->SetMarkerStyle(21);
  hy0Simu_Sys->SetMarkerColor(kBlue+1);
  hy0Simu_Sys->SetFillColor(kBlue-9);
  hy0Simu_Sys->SetFillStyle(3001);
  hy0Simu_Sys->SetLineColor(kBlue+1);
  ty0Simu->SetFillColor(kBlue+1-8);
  ty0Simu->SetFillStyle(3001);
  
  TCanvas * cy0 = new TCanvas("cy0","cy0",800,600);
  hy0Data_Stat->SetXTitle("Centrality %");
  hy0Data_Stat->SetYTitle("#LTy_{0}#GT Mean Rapidity Shift");
  hy0Data_Stat->SetMaximum(0.1);
  hy0Data_Stat->SetMinimum(0.0);
  hy0Data_Stat->Draw();
  //hy0Simu_Sys->Draw("same E3");
  //  hy0Data_Sys->Draw("same E3");
  ty0Simu->Draw("f");
  tg20y0Data->Draw("f");
  //hy0SimuInputSys->Draw("E3 same");
  ty0SimuInput->Draw("f");
  hy0Simu_Stat->Draw("same");
  hy0SimuInput->Draw("same");
  hy0Actual->Draw("same");
  hy0Data_Stat->Draw("same");
  
  TLatex * l = new TLatex(4,0.088,"Pb-Pb 2.76 GeV");
  l->SetTextColor(1);
  l->Draw();
  TLatex * l = new TLatex(20,0.088,"ALICE PRELIMINARY");
  l->SetTextColor(1);
  l->Draw();
  
  TMarker * m = new TMarker(18.0,0.017,8);
  m->SetMarkerColor(kRed+1);
  m->SetMarkerSize(1.2);
  m->Draw();
  l = new TLatex(19,0.016,"c_{1}#sigma_{#eta}^{2}: Data");
  l->SetTextColor(kRed+1);
  l->Draw();

  TMarker * m = new TMarker(4.0,0.072,29);
  m->SetMarkerSize(1.6);
  m->SetMarkerColor(kCyan+2);
  m->Draw();
  l = new TLatex(5,0.071,"#LTy_{0}#GT Using #alpha_{ZDC}^{Data}");
  l->SetTextColor(kCyan+2);
  l->Draw();
  
  TMarker * m = new TMarker(4.0,0.081,34);
  m->SetMarkerColor(kMagenta+1);
  m->SetMarkerSize(1.5);
  m->Draw();
  //  l = new TLatex(5,0.080,"#LTy_{0}#GT (#alpha_{Part}>0.0)");
  l = new TLatex(5,0.080,"#LT|y_{0}|#GT Glauber");
  l->SetTextColor(kMagenta+1);
  l->Draw();
  
  TMarker * m = new TMarker(18.0,0.009,21);
  m->SetMarkerColor(kBlue+1);
  m->Draw();
  //  l = new TLatex(19,0.016,"Simu (c_{1}#sigma_{#eta}^{2}),gaussian y");
  l = new TLatex(19,0.008,"c_{1}#sigma_{#eta}^{2}: Gaus(y_{0},#sigma_{0}+b(y-y_{0}))");
  l->SetTextColor(kBlue+1);
  l->Draw();

}
void DrawSlope(){
  TCanvas * cslope = new TCanvas("cslope","cslope",800,600);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetXTitle("Centrality %");
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetYTitle("c_{1} (From Fit To dN/d#eta Ratio)");
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetMinimum(0.0);
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetMaximum(0.007);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetMinimum(0.0);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetMaximum(0.007);
  hMeanSlopeReg12WithStat->SetMaximum(0.007);
  hMeanSlopeReg12WithStat->SetMinimum(0.00);

  h_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerSize(1.5);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetLineWidth(2);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerColor(kBlue+1);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerStyle(21);
  h_c1Reg12_Fromfpol3par1_SimWithStat->SetLineColor(kBlue+1);
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerSize(1.5);
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerColor(kBlue+1);
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerStyle(21);
  h_c1Reg12_Fromfpol3par1_SimWithSys->SetFillColor(kBlue-9);
  
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerSize(1.5);
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetLineWidth(2);
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerColor(kCyan+2);
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerStyle(34);
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetLineColor(kCyan+2);
  h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerSize(1.5);
  h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerColor(kCyan+2);
  h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetMarkerStyle(34);
  h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetFillColor(kCyan-8);

  tc1Simu->SetFillColor(kBlue-9);
  tc1Simu->SetFillStyle(3001);
  h_c1Reg12_Fromfpol3par1_SimWithStat->Draw();
  //  h_c1Reg12_Fromfpol3par1_SimWithSys->Draw("E4");
  
  h1g_c1Reg12_Fromfpol3par1_SimWithSys->Draw("E4 same");
    
  //hws_c1Reg12_Fromfpol3par1_SimWithSys->Draw("E4 same");
  hws_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerSize(1.6);
  hws_c1Reg12_Fromfpol3par1_SimWithStat->SetLineWidth(2);;
  hws_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerColor(kMagenta+1);
  hws_c1Reg12_Fromfpol3par1_SimWithStat->SetMarkerStyle(29);
  hws_c1Reg12_Fromfpol3par1_SimWithStat->SetLineColor(kMagenta+1);
  
  hws_c1Reg12_Fromfpol3par1_SimWithStat->Draw("same");
  
  

  hMeanSlopeReg12WithStat->SetMarkerSize(1.2);
  hMeanSlopeReg12WithStat->SetMarkerStyle(20);
  hMeanSlopeReg12WithStat->SetMarkerColor(kRed+1);
  hMeanSlopeReg12WithStat->SetLineColor(kRed+1);
  hMeanSlopeReg12WithStat->SetLineWidth(2);
  /*
  hMeanSlopeReg12WithSys->SetMarkerSize(1.5);
  hMeanSlopeReg12WithSys->SetMarkerStyle(20);
  hMeanSlopeReg12WithSys->SetMarkerColor(kRed+1);
  hMeanSlopeReg12WithSys->SetLineColor(kRed+1);
  hMeanSlopeReg12WithSys->SetFillColor(kRed-9);
  hMeanSlopeReg12WithSys->SetXTitle("Centrality %");
  hMeanSlopeReg12WithSys->SetYTitle("c_{1} (From Fit To dN/d#eta Ratio)");
  hMeanSlopeReg12WithSys->SetMaximum(0.007);
  hMeanSlopeReg12WithSys->SetMinimum(0.0);
  hMeanSlopeReg12WithSys->Draw("E3 same");
  */
  //hMeanSlopeReg12WithSys->Draw("E3");
  tg20SlopeData->SetFillColor(kRed-9);
 
  tg20SlopeData->SetFillStyle(3001);
  hMeanSlopeReg12WithStat->Draw();
  tc1Simu->Draw("f");
  tg20SlopeData->Draw("f same");
  h_c1Reg12_Fromfpol3par1_SimWithStat->Draw("same");
  h1g_c1Reg12_Fromfpol3par1_SimWithStat->Draw("same");
  hws_c1Reg12_Fromfpol3par1_SimWithStat->Draw("same");
  hMeanSlopeReg12WithStat->Draw("same");

  TLatex * l = new TLatex(4,0.0062,"Pb-Pb 2.76 GeV");
  l->SetTextColor(1);
  l->Draw();
  TLatex * l = new TLatex(20,0.0062,"ALICE PRELIMINARY");
  l->SetTextColor(1);
  l->Draw();
  
  l = new TLatex(5,0.0055,"Data");
  l->SetTextColor(kRed+1);
  l->Draw();
  TMarker * m = new TMarker(4, 0.0056, 20);
  m->SetMarkerColor(kRed+1);
  m->Draw();

  
  l = new TLatex(5,0.0049,"Gaus(y_{0},#sigma_{0}+b(y-y_{0}))");
  l->SetTextColor(kBlue+1);
  l->Draw();
  m = new TMarker(4, 0.0050, 21);
  m->SetMarkerColor(kBlue+1);
  m->Draw();
  
  l = new TLatex(5,0.0043,"Gaus(y_{0},#sigma_{0})");
  l->SetTextColor(kCyan+2);
  l->Draw();
  m = new TMarker(4, 0.0044, 34);
  m->SetMarkerColor(kCyan+3);
  m->SetMarkerSize(1.5);
  m->Draw();

  l = new TLatex(5,0.0037,"Wood-Saxon");
  l->SetTextColor(kMagenta+1);
  l->Draw();
  m = new TMarker(4, 0.0038, 29);
  m->SetMarkerColor(kMagenta+1);
  m->SetMarkerSize(1.6);
  m->Draw();
  
}

void SetStyle()
{
  gStyle->Reset();
  gROOT->SetStyle("Bold");

  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerSize(1.3);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetTextColor(1);
  gStyle->SetTextFont(42);
  gStyle->SetFillColor(0);
  gStyle->SetTitleTextColor(1);
  gStyle->SetLabelColor(1,"X");
  gStyle->SetLabelColor(1,"Y");
  gStyle->SetTitleOffset(1.1);  
  gStyle->SetTitleXOffset(1.1);  
  gStyle->SetTitleYOffset(1.5);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Z");
  gStyle->SetPadBottomMargin(0.19);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetLineWidth(2);
   //gROOT->ForceStyle();
}

void GetMeanStatSimu1gWithb(){

  Char_t fnameData1g[256];
  TH1F * h1gWithb_c1Reg1_Fromfpol3par1_Sim = new TH1F("h1gWithb_c1Reg1_Fromfpol3par1_Sim","h1gWithb_c1Reg1_Fromfpol3par1_Sim",8,0.0,40.0);
  TH1F * h1gWithb_c1Reg2_Fromfpol3par1_Sim = new TH1F("h1gWithb_c1Reg2_Fromfpol3par1_Sim","h1gWithb_c1Reg2_Fromfpol3par1_Sim",8,0.0,40.0);
  h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat = new TH1F("h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat","h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat",8,0.0,40.0);
  h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys = new TH1F("h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys","h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys",8,0.0,40.0);
 
  Float_t fdummy,StatErr[8];
  Float_t SimuReg1fpol3_Par1[8],SimuReg2fpol3_Par1[8],SimuReg1fpol3_dPar1[8],SimuReg2fpol3_dPar1[8];
  Int_t dummy;
  // Read Data2 files ( 8 centralities in a loop and fill histograms
  for(Int_t icent=0;icent<8;icent++){
    ifstream fin;
    sprintf(fnameData1g,"/Users/rashmiraniwala/ToALICE/NumericalSim/Sept2015/GaussWithb/GausWithb0.001yEta10K_41284_Cent%d_3.600.dat",icent);
    cout<<" Reading Simulated Data file "<<fnameData1g<<endl;
    fin.open(fnameData1g,ios::in);
    fin>>fdummy>>fdummy>>fdummy>>fdummy;
    fin>>fdummy>>fdummy>>fdummy;
    fin>>dummy>>dummy>>dummy;
    fin>>fdummy>>fdummy;
    fin>>fdummy>>fdummy;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;
    fin>>SimuReg1fpol3_Par1[icent]>>SimuReg1fpol3_dPar1[icent]>>SimuReg1fpol3_Par2[icent]>>SimuReg1fpol3_dPar2[icent];
    fin>>SimuReg2fpol3_Par1[icent]>>SimuReg2fpol3_dPar1[icent]>>SimuReg2fpol3_Par2[icent]>>SimuReg2fpol3_dPar2[icent];
    //cout<<SimuReg1fpol3_Par1[icent]<<" "<<SimuReg1fpol3_dPar1[icent]<<" "<<SimuReg1fpol3_Par2[icent]<<" "<<SimuReg1fpol3_dPar2[icent]<<endl;
    //cout<<SimuReg2fpol3_Par1[icent]<<" "<<SimuReg2fpol3_dPar1[icent]<<" "<<SimuReg2fpol3_Par2[icent]<<" "<<SimuReg2fpol3_dPar2[icent]<<endl;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;

    h1gWithb_c1Reg1_Fromfpol3par1_Sim->SetBinContent(icent+1,-1.0*SimuReg1fpol3_Par1[icent]);
    h1gWithb_c1Reg1_Fromfpol3par1_Sim->SetBinError(icent+1,-1.0*SimuReg1fpol3_dPar1[icent]);
    h1gWithb_c1Reg2_Fromfpol3par1_Sim->SetBinContent(icent+1,SimuReg2fpol3_Par1[icent]);
    h1gWithb_c1Reg2_Fromfpol3par1_Sim->SetBinError(icent+1,SimuReg2fpol3_dPar1[icent]);

    SimuReg12fpol3_Par1[icent] = 0.5*(TMath::Abs(SimuReg1fpol3_Par1[icent]) + TMath::Abs(SimuReg2fpol3_Par1[icent]));
    Float_t stat = (SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent])*(SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent]);
    stat += (SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent])*(SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent]);
    StatErr[icent] = SimuReg12fpol3_Par1[icent]*TMath::Sqrt(0.5*stat);
    
    h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    h1gWithb_c1Reg12_Fromfpol3par1_SimWithStat->SetBinError(icent+1,StatErr[icent]);
    h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    SysErrReg12[icent] = 0.5*(TMath::Abs(TMath::Abs(SimuReg1fpol3_Par1[icent]) - TMath::Abs(SimuReg2fpol3_Par1[icent])));
    cout<<"FromDiff = "<<-1.0*SysErrReg12[icent]*100/SimuReg1fpol3_Par1[icent]<<" "<<SysErrReg12[icent]<<endl;

    h1gWithb_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    cout<<"SimuReg12fpol3_Par1["<<icent<<"]="<<SimuReg12fpol3_Par1[icent]<<" stat%="<<StatErr[icent]*100.0/SimuReg12fpol3_Par1[icent]<<endl;
  }
  
}

void GetMeanStatSimu1g(){

  Char_t fnameData1g[256];
  TH1F * h1g_c1Reg1_Fromfpol3par1_Sim = new TH1F("h1g_c1Reg1_Fromfpol3par1_Sim","h1g_c1Reg1_Fromfpol3par1_Sim",8,0.0,40.0);
  TH1F *h1g_c1Reg2_Fromfpol3par1_Sim = new TH1F("h1g_c1Reg2_Fromfpol3par1_Sim","h1g_c1Reg2_Fromfpol3par1_Sim",8,0.0,40.0);
  h1g_c1Reg12_Fromfpol3par1_SimWithStat = new TH1F("h1g_c1Reg12_Fromfpol3par1_SimWithStat","h1g_c1Reg12_Fromfpol3par1_SimWithStat",8,0.0,40.0);
  h1g_c1Reg12_Fromfpol3par1_SimWithSys = new TH1F("h1g_c1Reg12_Fromfpol3par1_SimWithSys","h1g_c1Reg12_Fromfpol3par1_SimWithSys",8,0.0,40.0);
 
  Float_t fdummy,StatErr[8];
  Float_t SimuReg1fpol3_Par1[8],SimuReg2fpol3_Par1[8],SimuReg1fpol3_dPar1[8],SimuReg2fpol3_dPar1[8];
  Int_t dummy;
  // Read Data2 files ( 8 centralities in a loop and fill histograms
  Float_t ySigma[8]={3.75,3.79,3.76,3.76,3.79,3.79,3.83,3.83};
  for(Int_t icent=0;icent<8;icent++){
    ifstream fin;
    //    sprintf(fnameData1g,"/Users/rashmiraniwala/ToALICE/NumericalSim/Sept2015/Gaussian/GausyEta10K_41284_Cent%d_3.600.dat",icent);
    sprintf(fnameData1g,"/Users/rashmiraniwala/ToALICE/NumericalSim/Sept2015/Gaussian/WiderEta/GausyEta10K_41284_Cent%d_%4.3f.dat",icent,ySigma[icent]);
    cout<<" Reading Simulated Data file "<<fnameData1g<<endl;
    fin.open(fnameData1g,ios::in);
    fin>>fdummy>>fdummy>>fdummy>>fdummy;
    fin>>fdummy>>fdummy>>fdummy;
    fin>>dummy>>dummy>>dummy;
    fin>>fdummy>>fdummy;
    fin>>fdummy>>fdummy;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;
    fin>>SimuReg1fpol3_Par1[icent]>>SimuReg1fpol3_dPar1[icent]>>SimuReg1fpol3_Par2[icent]>>SimuReg1fpol3_dPar2[icent];
    fin>>SimuReg2fpol3_Par1[icent]>>SimuReg2fpol3_dPar1[icent]>>SimuReg2fpol3_Par2[icent]>>SimuReg2fpol3_dPar2[icent];
    //cout<<SimuReg1fpol3_Par1[icent]<<" "<<SimuReg1fpol3_dPar1[icent]<<" "<<SimuReg1fpol3_Par2[icent]<<" "<<SimuReg1fpol3_dPar2[icent]<<endl;
    //cout<<SimuReg2fpol3_Par1[icent]<<" "<<SimuReg2fpol3_dPar1[icent]<<" "<<SimuReg2fpol3_Par2[icent]<<" "<<SimuReg2fpol3_dPar2[icent]<<endl;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;

    h1g_c1Reg1_Fromfpol3par1_Sim->SetBinContent(icent+1,-1.0*SimuReg1fpol3_Par1[icent]);
    h1g_c1Reg1_Fromfpol3par1_Sim->SetBinError(icent+1,-1.0*SimuReg1fpol3_dPar1[icent]);
    h1g_c1Reg2_Fromfpol3par1_Sim->SetBinContent(icent+1,SimuReg2fpol3_Par1[icent]);
    h1g_c1Reg2_Fromfpol3par1_Sim->SetBinError(icent+1,SimuReg2fpol3_dPar1[icent]);

    SimuReg12fpol3_Par1[icent] = 0.5*(TMath::Abs(SimuReg1fpol3_Par1[icent]) + TMath::Abs(SimuReg2fpol3_Par1[icent]));
    Float_t stat = (SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent])*(SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent]);
    stat += (SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent])*(SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent]);
    StatErr[icent] = SimuReg12fpol3_Par1[icent]*TMath::Sqrt(0.5*stat);
    cout<<" mean = "<<SimuReg12fpol3_Par1[icent]<<" stat = "<<StatErr[icent]<<endl;
    h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    h1g_c1Reg12_Fromfpol3par1_SimWithStat->SetBinError(icent+1,StatErr[icent]);
    h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    SysErrReg12[icent] = 0.5*(TMath::Abs(TMath::Abs(SimuReg1fpol3_Par1[icent]) - TMath::Abs(SimuReg2fpol3_Par1[icent])));
    cout<<"FromDiff = "<<-1.0*SysErrReg12[icent]*100/SimuReg1fpol3_Par1[icent]<<" "<<SysErrReg12[icent]<<endl;

    h1g_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    cout<<"SimuReg12fpol3_Par1["<<icent<<"]="<<SimuReg12fpol3_Par1[icent]<<" stat%="<<StatErr[icent]*100.0/SimuReg12fpol3_Par1[icent]<<endl;
  }
  
  
}

void GetMeanStatSimuWS(){

  Char_t fnameDataWS[256];
  TH1F * hws_c1Reg1_Fromfpol3par1_Sim = new TH1F("hws_c1Reg1_Fromfpol3par1_Sim","hws_c1Reg1_Fromfpol3par1_Sim",8,0.0,40.0);
  TH1F * hws_c1Reg2_Fromfpol3par1_Sim = new TH1F("hws_c1Reg2_Fromfpol3par1_Sim","hws_c1Reg2_Fromfpol3par1_Sim",8,0.0,40.0);
  hws_c1Reg12_Fromfpol3par1_SimWithStat = new TH1F("hws_c1Reg12_Fromfpol3par1_SimWithStat","hws_c1Reg12_Fromfpol3par1_SimWithStat",8,0.0,40.0);
  hws_c1Reg12_Fromfpol3par1_SimWithSys = new TH1F("hws_c1Reg12_Fromfpol3par1_SimWithSys","hws_c1Reg12_Fromfpol3par1_SimWithSys",8,0.0,40.0);
 
  Float_t fdummy,StatErr[8];
  Float_t SimuReg1fpol3_Par1[8],SimuReg2fpol3_Par1[8],SimuReg1fpol3_dPar1[8],SimuReg2fpol3_dPar1[8];
  Int_t dummy;
  // Read Data2 files ( 8 centralities in a loop and fill histograms
  for(Int_t icent=0;icent<8;icent++){
    ifstream fin;
    sprintf(fnameDataWS,"/Users/rashmiraniwala/ToALICE/NumericalSim/Sept2015/WoodSaxon/WSyEta10K_41284_Cent%d_4.424_1.132.dat",icent);
    cout<<" Reading Simulated Data file "<<fnameDataWS<<endl;
    fin.open(fnameDataWS,ios::in);
    fin>>fdummy>>fdummy>>fdummy>>fdummy;
    fin>>fdummy>>fdummy>>fdummy;
    fin>>dummy>>dummy>>dummy;
    fin>>fdummy>>fdummy;
    fin>>fdummy>>fdummy;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;
    fin>>SimuReg1fpol3_Par1[icent]>>SimuReg1fpol3_dPar1[icent]>>SimuReg1fpol3_Par2[icent]>>SimuReg1fpol3_dPar2[icent];
    fin>>SimuReg2fpol3_Par1[icent]>>SimuReg2fpol3_dPar1[icent]>>SimuReg2fpol3_Par2[icent]>>SimuReg2fpol3_dPar2[icent];
    //cout<<SimuReg1fpol3_Par1[icent]<<" "<<SimuReg1fpol3_dPar1[icent]<<" "<<SimuReg1fpol3_Par2[icent]<<" "<<SimuReg1fpol3_dPar2[icent]<<endl;
    //cout<<SimuReg2fpol3_Par1[icent]<<" "<<SimuReg2fpol3_dPar1[icent]<<" "<<SimuReg2fpol3_Par2[icent]<<" "<<SimuReg2fpol3_dPar2[icent]<<endl;
    //Not reading the chisquare and ndf
    fin>>fdummy>>dummy;
    fin>>fdummy>>dummy;

    hws_c1Reg1_Fromfpol3par1_Sim->SetBinContent(icent+1,-1.0*SimuReg1fpol3_Par1[icent]);
    hws_c1Reg1_Fromfpol3par1_Sim->SetBinError(icent+1,-1.0*SimuReg1fpol3_dPar1[icent]);
    hws_c1Reg2_Fromfpol3par1_Sim->SetBinContent(icent+1,SimuReg2fpol3_Par1[icent]);
    hws_c1Reg2_Fromfpol3par1_Sim->SetBinError(icent+1,SimuReg2fpol3_dPar1[icent]);

    SimuReg12fpol3_Par1[icent] = 0.5*(TMath::Abs(SimuReg1fpol3_Par1[icent]) + TMath::Abs(SimuReg2fpol3_Par1[icent]));
    Float_t stat = (SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent])*(SimuReg1fpol3_dPar1[icent]/SimuReg1fpol3_Par1[icent]);
    stat += (SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent])*(SimuReg2fpol3_dPar1[icent]/SimuReg2fpol3_Par1[icent]);
    StatErr[icent] = SimuReg12fpol3_Par1[icent]*TMath::Sqrt(0.5*stat);
    
    hws_c1Reg12_Fromfpol3par1_SimWithStat->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    hws_c1Reg12_Fromfpol3par1_SimWithStat->SetBinError(icent+1,StatErr[icent]);
    hws_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    SysErrReg12[icent] = 0.5*(TMath::Abs(TMath::Abs(SimuReg1fpol3_Par1[icent]) - TMath::Abs(SimuReg2fpol3_Par1[icent])));
    cout<<"FromDiff = "<<-1.0*SysErrReg12[icent]*100/SimuReg1fpol3_Par1[icent]<<" "<<SysErrReg12[icent]<<endl;

    hws_c1Reg12_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,SimuReg12fpol3_Par1[icent]);
    cout<<"SimuReg12fpol3_Par1["<<icent<<"]="<<SimuReg12fpol3_Par1[icent]<<" stat%="<<StatErr[icent]*100.0/SimuReg12fpol3_Par1[icent]<<endl;
  }
  
  
}

TGraph* GetSysPlot( TH1F *hUp, TH1F * hDown) {

  Float_t ValHigh[8], ValLow[8];
  for(Int_t icent = 0;icent<8;icent++){
    ValHigh[icent] = hUp->GetBinContent(icent+1);
    ValLow[icent] = hDown->GetBinContent(icent+1);
  }
  //  TGraph * tg20 = new TGraph(20);
  TGraph * tg20 = new TGraph(16);
  tg20->SetNameTitle("tg20ForSys","tg20ForSys");
  TGraph * tg07 = new TGraph(8,cent,ValHigh);
  tg07->SetNameTitle("tg07","tg07");
  TGraph * tg815 = new TGraph(8,cent,ValLow);
  tg815->SetNameTitle("tg815","tg815");
  tg07->Fit("pol3");
  tg815->Fit("pol3");
  TF1 * f1 = tg07->GetFunction("pol3");
  TF1 * f2 = tg815->GetFunction("pol3");

  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f1(cent[icent]);
    // tg20->SetPoint(icent+1,cent[icent],y);
     tg20->SetPoint(icent,cent[icent],ValHigh[icent]);
  }
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f2(cent[icent]);
    //    tg20->SetPoint(18-icent,cent[icent],y);
    tg20->SetPoint(15-icent,cent[icent],ValLow[icent]);
  }
  /*
  tg20->SetPoint(0,0.0,f1(0.0));
  tg20->SetPoint(9,40.0,f1(40.0));
  tg20->SetPoint(10,40.0,f2(40.0));
  tg20->SetPoint(19,0.0,f2(0.0));
  */
  return tg20;
}
TGraph* GetSysPlot1(TH1F *hSys) {

  Float_t ValHigh[8], ValLow[8];
  for(Int_t icent = 0;icent<8;icent++){
    ValHigh[icent] = hSys->GetBinContent(icent+1)+hSys->GetBinError(icent+1);
    ValLow[icent] =  hSys->GetBinContent(icent+1)-hSys->GetBinError(icent+1);
  }
  TGraph * tg20 = new TGraph(16);
  tg20->SetNameTitle("tg20ForSys","tg20ForSys");
  TGraph * tg07 = new TGraph(8,cent,ValHigh);
  tg07->SetNameTitle("tg07","tg07");
  TGraph * tg815 = new TGraph(8,cent,ValLow);
  tg815->SetNameTitle("tg815","tg815");
  tg07->Fit("pol3");
  tg815->Fit("pol3");
  TF1 * f1 = tg07->GetFunction("pol3");
  TF1 * f2 = tg815->GetFunction("pol3");

  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f1(cent[icent]);
    //    tg20->SetPoint(icent+1,cent[icent],y);
    tg20->SetPoint(icent,cent[icent],ValHigh[icent]);
  }
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f2(cent[icent]);
    //    tg20->SetPoint(18-icent,cent[icent],y);
    tg20->SetPoint(15-icent,cent[icent],ValLow[icent]);
  }
  /*
  tg20->SetPoint(0,0.0,f1(0.0));
  tg20->SetPoint(9,40.0,f1(40.0));
  tg20->SetPoint(10,40.0,f2(40.0));
  tg20->SetPoint(19,0.0,f2(0.0));
  */
  return tg20;
}
