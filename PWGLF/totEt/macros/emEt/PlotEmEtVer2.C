Int_t colors[] = {TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};
Int_t CutSet = 0;//Defaults: 250 MeV for PHOS and 300 MeV for EMCal
//1:  350 MeV for both
Float_t fem = 0.246;
Float_t femerr = 0.028;
void SetStyles(TGraph *graph, Int_t marker, Int_t color){
  graph->SetMarkerStyle(marker);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerSize(1.5);
}
void WriteLatex();
Float_t finalemetCorrEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finalemetCorrPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finalemetErrorEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finalemetErrorPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finaltotaletCorrEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finaltotaletCorrPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finaltotaletErrorEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t finaltotaletErrorPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
//neutron
Float_t GetMostLikelyValue(TH1 *histo){
  Float_t max = 0;
  Float_t maxx=0;
  for(Int_t bin=1;bin<histo->GetNbinsX();bin++){
    if(histo->GetBinContent(bin)>max){
      max = histo->GetBinContent(bin);
      maxx = histo->GetBinCenter(bin);
    }
  }
  return maxx;
}

Float_t npartShort[10] =    {382.8,329.7,260.5,186.4,128.9, 85,52.8,30.0,15.8,7.48};
Float_t npartErrShort[10] = {    6,    6,  4.4,  3.9,  3.3,2.6,   2, 1.3, 0.6,0.29};
Float_t npart[20] = {382.7, 329.4, 281.2, 239, 202.1, 169.5, 141, 116, 94.11, 75.3, 59.24, 45.58, 34.33, 25.21, 17.96, 12.58, 8.812, 6.158, 4.376, 3.064};
Float_t npartErr[20] = {3, 4.3, 4.1, 3.5, 3.3, 3.3, 3.1, 2.8, 2.6, 2.3, 1.8, 1.4, 1.1, 0.87, 0.66, 0.45, 0.26, 0.19, 0.1, 0.059};

//========================Charged Pion Reference========================================
//Arrays for defining comparison plots
Float_t pionPlusEt[10] = {360.7,298.3,223.8,149.9,96.1, 58.1,32.4,16.4,7.3,2.7};
Float_t pionMinusEt[10] ={363.7,300.4,225.4,150.5,96.6, 58.4,32.5,16.5,7.4,2.8};
Float_t pionEtError[10] = {19.3, 15.3,11.3 ,7.5  , 4.8,  2.9, 1.6, 0.8,0.4,0.1};
Float_t pionEt[10] = {0,0,0,0,0, 0,0,0,0,0};
Float_t ypion[10]  = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t ypionerr[10]  = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPions[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsPerNpart[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsPerNpartErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t npartAlt1[20],npartAlt2[20],npartAlt3[20];

TGraphErrors *GetPionEtGraph(){
  for(int i=0;i<10;i++){
    pionEt[i] = (pionPlusEt[i]+pionMinusEt[i])/2.0;
    //emEtFromPions[i] = pionEt[i]*1.085;
    emEtFromPions[i] = pionEt[i]*1.171;
    emEtFromPionsErr[i] = emEtFromPions[i]*TMath::Sqrt(TMath::Power(0.11/1.171,2)+TMath::Power(pionEtError[i]/pionEt[i],2));
    ypion[i] = pionEt[i]/(npartShort[i]/2);
    ypionerr[i] = pionEtError[i]/(npartShort[i]/2);
    emEtFromPionsPerNpart[i] = emEtFromPions[i]/(npartShort[i]/2);
    emEtFromPionsPerNpartErr[i] = emEtFromPionsErr[i]/(npartShort[i]/2);
    npartAlt1[i] = npartShort[i]+2;
    npartAlt2[i] = npartShort[i]-2;
    npartAlt3[i] = npartShort[i]+4;
  }
  TGraphErrors *gr2 = new TGraphErrors(10,npartShort,ypion,npartErrShort,ypionerr);
  gr2->GetYaxis()->SetTitle("dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]");
  gr2->GetXaxis()->SetTitle("N_{part}");
  gr2->SetTitle("");
  gr2->GetXaxis()->SetRangeUser(0, 400);
  SetStyles(gr2,30,TColor::kBlue);

  return gr2;
}

TGraphErrors *GetPionEmEtGraph(){
    TGraphErrors *gr3 = new TGraphErrors(10,npartAlt3,emEtFromPionsPerNpart,npartErrShort,emEtFromPionsPerNpartErr);
    SetStyles(gr3,29,TColor::kBlue);
    return gr3;

}
//========================Reading in corrections========================================
Float_t energyscaleerror = 0.02;
Float_t nonLinError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t nonLinErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
//Float_t signalFraction[20] = {0.4,0.4,0.4,0.4,0.4, 0.4,0.4,0.4,0.4,0.4, 0.4,0.4,0.4,0.4,0.4, 0.4,0.4,0.4,0.4,0.4};
Float_t signalFraction[20] = {0.3,0.3,0.3,0.3,0.3, 0.3,0.3,0.3,0.3,0.3, 0.3,0.3,0.3,0.3,0.3, 0.3,0.3,0.3,0.3,0.3};
Float_t signalFractionError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
//Float_t averageEfficiency[20] = {0.5,0.5,0.5,0.5,0.5, 0.5,0.5,0.5,0.5,0.5, 0.5,0.5,0.5,0.5,0.5, 0.5,0.5,0.5,0.5,0.5};
Float_t averageEfficiency[20] = {1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0};
Float_t averageEfficiencyError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t efficiencyError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t efficiencyErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNChShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrPerNChShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNChShort[10] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNCh[20] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};


Float_t neutronCorrNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorShortNoEffCorr[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrShortNoEffCorr[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNChShortNoEffCorr[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrPerNChShortNoEffCorr[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorrPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNChShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNChShortNoEffCorr[10] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNChNoEffCorr[20] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNChShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNChShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadErrorShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorrShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadErrorNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorrNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNChShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNChShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNChNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};




Float_t kaonYield[10] = {109,90.5,68,46,30,18.2,10.2,5.1,2.3,0.855};
Float_t kaonYieldStatErr[10] = {0.3,0.2,0.1,0.1,0.1, 0.06,0.04,0.03,0.02,0.01};
Float_t kaonYieldSysErr[10] = {9,7,5,4,2, 1.5,0.8,0.4,0.2,0.09};
Float_t kaonYieldTotErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonYieldPerNCh[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonYieldPerNChErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonEtPerNCh[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonEtPerNChErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Double_t kaonPlusEt[2][10] = {{91.7712,75.8971,56.841,37.6962,23.9923,14.255,7.73469,3.74477,1.60505,0.578278},{6.61811,5.39337,3.9978,2.6337,1.67854,1.01849,0.557879,0.278199,0.125057,0.0592682}};
Double_t kaonMinusEt[2][10] = {{90.4723,74.9444,55.9463,37.286,23.6591,14.0413,7.63067,3.69337,1.59219,0.571019},{7.01588,5.76588,4.20933,2.80388,1.77983,1.06934,0.588003,0.292737,0.138191,0.0600075}};

Float_t averageHadEnergy = -1;
Float_t averageHadEnergyError = -1;

ofstream myfileHadCorrTable;


void ReadMinEtCorrections(){
  cout<<"Reading in min et corrections..."<<endl;
  string inline;
  float value = 0;
  float error = 0;
  float err = 0;
  int i=0;
  //file names like /tmp/MinEtEmcalShortCut7.dat
  //1 = 100 MeV/c
  //2 = 150 MeV/c
  //3 = 200 MeV/c
  //4 = 250 MeV/c
  //5 = 300 MeV/c
  //6 = 350 MeV/c
  //7 = 400 MeV/c
  //8 = 450 MeV/c
  //9 = 500 MeV/c
  //10 = 550 MeV/c
  TString cutstring = "";
  if(CutSet==1) cutstring = "Cut6";
  if(CutSet==2) cutstring = "Cut7";
  TString minetInfileNameShort = "MinEt"+detector+"Short"+cutstring+".dat";
  cout<<"Reading "<<minetInfileNameShort.Data()<<endl;
  ifstream myminetfileShort (minetInfileNameShort.Data());
  if (myminetfileShort.is_open()){
    while ( myminetfileShort.good() )
      {
	getline (myminetfileShort,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  minEtCorrShort[i] = value;
	  minEtErrorShort[i] = error;
	}
	//cout<<"min et corr cb "<<i<<" "<< minEtCorrShort[i]<<" +/- "<<minEtErrorShort[i]<<endl;
	i++;
      }
    myminetfileShort.close();
  }
  //doing the linear interpolation between the two data points we have
  TGraphErrors *graphMinEtCorrectionShort = GetMinEtCorrectionGraphShort();
  int shortbin = 0;
  for(int i=0;i<19;i++){
    //cout<<"Long bin "<<i<<" short bin "<<shortbin;
    if(i<2){//we have exact numbers so we don't need to interpolate
      minEtCorr[i] = minEtCorrShort[i];
      minEtError[i] = minEtErrorShort[i];
      shortbin++;
    }
    else{
      minEtCorr[i] = graphMinEtCorrectionShort->Eval(trackmultiplicity[i]);
      int altbin = shortbin-1;
      if(i%2==1){altbin = shortbin+1;}
      //cout<<" altbin "<<altbin;
      if(minEtErrorShort[shortbin]>minEtErrorShort[altbin]) minEtError[i] = minEtErrorShort[shortbin];
      else{minEtError[i] = minEtErrorShort[altbin];}
      if(i%2==1 && shortbin<10){shortbin++;}
    }
    // cout<<"min et corr cb "<<i<<" "<< minEtCorr[i]<<" +/- "<<minEtError[i]<<endl;
    //cout<<endl;
  }
  delete graphMinEtCorrectionShort;

}
void ReadInNeutronCorrections(){
  cout<<"Reading in neutron corrections..."<<endl;
  TString neutronInfileName = "Neutrons"+detector+".dat";
  ifstream myneutronfile (neutronInfileName.Data());
  string inline;
  float value = 0;
  float error = 0;
  int i=0;
  if (myneutronfile.is_open()){
    while ( myneutronfile.good() )
      {
	getline (myneutronfile,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<20){
	  neutronCorr[i] = value;
	  neutronError[i] = error;
	  if(trackmultiplicity[i]>0){
	    neutronCorrPerNCh[i] = value/(trackmultiplicity[i])/2.0;
	    neutronErrorPerNCh[i] = error/(trackmultiplicity[i])/2.0;
	  }
	}

	// cout<<"neutroncorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	i++;
      }
    myneutronfile.close();

  }
  TString neutronInfileNameShort = "Neutrons"+detector+"Short.dat";
  ifstream myneutronfileShort (neutronInfileNameShort.Data());
  i=0;
  if (myneutronfileShort.is_open()){
    while ( myneutronfileShort.good() )
      {
	getline (myneutronfileShort,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  neutronCorrShort[i] = value;
	  neutronErrorShort[i] = error;
	  neutronCorrPerNChShort[i] = value/(trackmultiplicityShort[i])/2.0;
	  neutronErrorPerNChShort[i] = error/(trackmultiplicityShort[i])/2.0;
	}
	//cout<<"neutroncorr cb "<<i<<" "<<neutronCorrShort[i]<<" +/- "<<neutronErrorShort[i]<<endl;
	i++;
      }
    myneutronfileShort.close();
  }
  //Begin reading in no eff corr corrections
  TString neutronInfileName3 = "Neutrons"+detector+"NoEffCorr.dat";
  ifstream myneutronfile3 (neutronInfileName3.Data());
  value = 0;
  error = 0;
  i=0;
  if (myneutronfile3.is_open()){
    while ( myneutronfile3.good() )
      {
	getline (myneutronfile3,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<20){
	  neutronCorrNoEffCorr[i] = value;
	  neutronErrorNoEffCorr[i] = error;
	  if(trackmultiplicity[i]>0){
	    neutronCorrPerNChNoEffCorr[i] = value/(trackmultiplicity[i])/2.0;
	    neutronErrorPerNChNoEffCorr[i] = error/(trackmultiplicity[i])/2.0;
	  }
	}

	//cout<<"neutroncorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	i++;
      }
    myneutronfile3.close();

  }
  TString neutronInfileNameShort3 = "Neutrons"+detector+"NoEffCorrShort.dat";
  ifstream myneutronfile3Short (neutronInfileNameShort3.Data());
  i=0;
  if (myneutronfile3Short.is_open()){
    while ( myneutronfile3Short.good() )
      {
	getline (myneutronfile3Short,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  neutronCorrShortNoEffCorr[i] = value;
	  neutronErrorShortNoEffCorr[i] = error;
	  neutronCorrPerNChShortNoEffCorr[i] = value/(trackmultiplicityShort[i])/2;
	  neutronErrorPerNChShortNoEffCorr[i] = error/(trackmultiplicityShort[i])/2;
	}
	//cout<<"neutroncorr cb "<<i<<" "<<neutronCorrShortNoEffCorr[i]<<" +/- "<<neutronErrorShortNoEffCorr[i]<<endl;
	i++;
      }
    myneutronfile3Short.close();
  }
  

}
void ReadInSecondaryCorrections(){
  cout<<"Reading in secondary corrections..."<<endl;

    TString secondaryInfileName = "Secondaries"+detector+".dat";
    ifstream mysecondaryfile (secondaryInfileName.Data());
    string inline;
    float value = 0;
    float error = 0;
    int i=0;
    if (mysecondaryfile.is_open()){
      while ( mysecondaryfile.good() )
	{
	  getline (mysecondaryfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    secondaryCorr[i] = value;
	    secondaryError[i] = error;
	    if(trackmultiplicity[i]>0){
	      secondaryCorrPerNCh[i] = value/(trackmultiplicity[i])/5.0;
	      secondaryErrorPerNCh[i] = error/(trackmultiplicity[i])/5.0;
	    }
	  }
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	  i++;
	}
        mysecondaryfile.close();
    }
    TString secondaryShortInfileName = "Secondaries"+detector+"Short.dat";
    ifstream mysecondaryShortfile (secondaryShortInfileName.Data());
    i=0;
    if (mysecondaryShortfile.is_open()){
      while ( mysecondaryShortfile.good() && i<10 )
	{
	  getline (mysecondaryShortfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<10){
	    secondaryCorrShort[i] = value;
	    secondaryErrorShort[i] = error;
	    secondaryCorrPerNChShort[i] = value/(trackmultiplicityShort[i])/5.0;
	    secondaryErrorPerNChShort[i] = error/(trackmultiplicityShort[i])/5.0;
	  }
 	  i++;
	}
       mysecondaryShortfile.close();
    }
    //Begin reading in no eff corr corrections
    TString secondaryInfileName = "Secondaries"+detector+"NoEffCorr.dat";
    ifstream mysecondaryfile3 (secondaryInfileName.Data());
    value = 0;
    error = 0;
    i=0;
    if (mysecondaryfile3.is_open()){
      while ( mysecondaryfile3.good() )
	{
	  getline (mysecondaryfile3,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    secondaryCorrNoEffCorr[i] = value;
	    secondaryErrorNoEffCorr[i] = error;
	    if(trackmultiplicity[i]>0){
	      secondaryCorrPerNChNoEffCorr[i] = value/(trackmultiplicity[i])/5.0;
	      secondaryErrorPerNChNoEffCorr[i] = error/(trackmultiplicity[i])/5.0;
	    }
	  }
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<" "<<secondaryCorrPerNChNoEffCorr[i]<<" +/- "<<secondaryErrorPerNChNoEffCorr[i]<<endl;
	  i++;
	}
        mysecondaryfile3.close();
    }
    TString secondaryShortInfileName = "Secondaries"+detector+"NoEffCorrShort.dat";
    ifstream mysecondaryShortfile3 (secondaryShortInfileName.Data());
    i=0;
    if (mysecondaryShortfile3.is_open()){
      while ( mysecondaryShortfile3.good() && i<10 )
	{
	  getline (mysecondaryShortfile3,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<10){
	    secondaryCorrShortNoEffCorr[i] = value;
	    secondaryErrorShortNoEffCorr[i] = error;
	    secondaryCorrPerNChShortNoEffCorr[i] = value/(trackmultiplicityShort[i])/5.0;
	    secondaryErrorPerNChShortNoEffCorr[i] = error/(trackmultiplicityShort[i])/5.0;
	  }
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<" "<<secondaryCorrPerNChShortNoEffCorr[i]<<" +/- "<<secondaryErrorPerNChShortNoEffCorr[i]<<endl;
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
 	  i++;
	}
       mysecondaryShortfile3.close();
    }

}
void ReadInKaonCorrections(){
  cout<<"Reading in kaon corrections..."<<endl;
  //junk.PHOS.CutNum6.txt
  TString cutstring = "7";
  if(detector.Contains("P")){ cutstring = "6";}
  if(CutSet==1){
    cutstring = "8";
  }
  if(CutSet==2){
    cutstring = "9";
  }
  TString kaonInfileName = "../spectrafits/KaonCut"+cutstring+"EMCalCorr."+year+".dat";
  if(detector.Contains("P")){
    kaonInfileName = "../spectrafits/KaonCut"+cutstring+"PHOSCorr."+year+".dat";
  }
  cout<<"Reading in "<<kaonInfileName<<endl;
    ifstream mykaonfile (kaonInfileName.Data());
    string inline;
    float value = 0;
    float error = 0;
    int i=0;
    if (mykaonfile.is_open()){
      while ( mykaonfile.good() )
	{
	  getline (mykaonfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<10){
	    kaonCorrShort[i] = value;
	    kaonErrorShort[i] = error;
	    kaonCorrPerNChShort[i] = value/(trackmultiplicityShort[i]);
	    kaonErrorPerNChShort[i] = error/(trackmultiplicityShort[i]);
	    // cout<<"kaoncorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	  }
	  i++;
	}
        mykaonfile.close();
    }


  TGraphErrors *graphKaonCorrectionShort = GetKaonCorrectionGraphShort();
  int shortbin = 0;
  for(int i=0;i<19;i++){
    //cout<<"Long bin "<<i<<" short bin "<<shortbin;
    if(i<2){//we have exact numbers so we don't need to interpolate
      kaonCorr[i] = kaonCorrShort[i];
      kaonError[i] = kaonErrorShort[i];
      shortbin++;
    }
    else{
      kaonCorr[i] = graphKaonCorrectionShort->Eval(trackmultiplicity[i]) * trackmultiplicity[i];
      int altbin = shortbin-1;
      if(i%2==1){altbin = shortbin+1;}
      //cout<<" altbin "<<altbin;
      if(kaonErrorPerNChShort[shortbin]>kaonErrorPerNChShort[altbin]) kaonError[i] = kaonErrorPerNChShort[shortbin] * trackmultiplicity[i];
      else{kaonError[i] =  kaonErrorPerNChShort[altbin] * trackmultiplicity[i];}
      if(i%2==1 && shortbin<10){shortbin++;}
    }
    //cout<<"kaoncorr cb "<<i<<" "<<kaonCorr[i]<<" +/- "<<kaonError[i]<<endl;
    kaonCorrPerNCh[i] = kaonCorr[i]/(trackmultiplicity[i]);
    kaonErrorPerNCh[i] = kaonError[i]/(trackmultiplicity[i]);
    //cout<<" track multiplicity "<<trackmultiplicity[i]<<" kaon corr/mult*1000 "<<kaonCorrPerNCh[i]*1000<<" +/- "<<kaonErrorPerNCh[i]*1000<<endl;
    //cout<<"min et corr cb "<<i<<" "<< kaonCorr[i]<<" +/- "<<kaonError[i]<<endl;
    //cout<<endl;
  }
  delete graphKaonCorrectionShort;


    for(int i=0;i<10;i++){
      kaonYieldTotErr[i] = TMath::Sqrt(kaonYieldStatErr[i]*kaonYieldStatErr[i]+kaonYieldSysErr[i]*kaonYieldSysErr[i]);
      kaonYieldPerNCh[i] = kaonYield[i]/(trackmultiplicityShort[i])/5;
      kaonYieldPerNChErr[i] = kaonYieldTotErr[i]/(trackmultiplicityShort[i])/5;
      float total = kaonPlusEt[0][i]+kaonMinusEt[0][i];
      float err = kaonPlusEt[1][i]+kaonMinusEt[1][i];
      kaonEtPerNCh[i] = total/(trackmultiplicityShort[i])/10;
      kaonEtPerNChErr[i] = err/(trackmultiplicityShort[i])/10;
    }
    //corrections for no eff corr
  TString kaonInfileName3 = "../spectrafits/KaonCut"+cutstring+"EMCalNoEffCorr."+year+".dat";
  if(detector.Contains("P")){
    kaonInfileName3 = "../spectrafits/KaonCut"+cutstring+"PHOSNoEffCorr.dat";
  }
  cout<<"Reading in "<<kaonInfileName3<<endl;
  ifstream mykaonfile3 (kaonInfileName3.Data());
  value = 0;
  error = 0;
  i=0;
  if (mykaonfile3.is_open()){
    while ( mykaonfile3.good() )
      {
	getline (mykaonfile3,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  kaonCorrShortNoEffCorr[i] = value;
	    kaonErrorShortNoEffCorr[i] = error;
	    kaonCorrPerNChShortNoEffCorr[i] = value/(trackmultiplicityShort[i]);
	    kaonErrorPerNChShortNoEffCorr[i] = error/(trackmultiplicityShort[i]);
	    //cout<<"kaoncorr cb "<<i<<" "<<value<<" +/- "<<error;
	    //cout<<endl;
	}
	i++;
      }
    mykaonfile3.close();
  }


  TGraphErrors *graphKaonCorrectionShort2 = GetKaonCorrectionGraphShortNoEffCorr();
  shortbin = 0;
  for(int i=0;i<19;i++){
    //cout<<"Long bin "<<i<<" short bin "<<shortbin;
    if(i<2){//we have exact numbers so we don't need to interpolate
      kaonCorrNoEffCorr[i] = kaonCorrShortNoEffCorr[i];
      kaonErrorNoEffCorr[i] = kaonErrorShortNoEffCorr[i];
      shortbin++;
    }
    else{
      kaonCorrNoEffCorr[i] = graphKaonCorrectionShort2->Eval(trackmultiplicity[i]) * trackmultiplicity[i];
      int altbin = shortbin-1;
      if(i%2==1){altbin = shortbin+1;}
      //cout<<" altbin "<<altbin;
      if(kaonErrorPerNChShortNoEffCorr[shortbin]>kaonErrorPerNChShortNoEffCorr[altbin]) kaonErrorNoEffCorr[i] = kaonErrorPerNChShortNoEffCorr[shortbin] * trackmultiplicity[i];
      else{kaonErrorNoEffCorr[i] =  kaonErrorPerNChShortNoEffCorr[altbin] * trackmultiplicity[i];}
      if(i%2==1 && shortbin<10){shortbin++;}
    }
    //cout<<"kaoncorr cb "<<i<<" "<<kaonCorrNoEffCorr[i]<<" +/- "<<kaonErrorNoEffCorr[i]<<endl;
    kaonCorrPerNChNoEffCorr[i] = kaonCorrNoEffCorr[i]/(trackmultiplicity[i]);
    kaonErrorPerNChNoEffCorr[i] = kaonErrorNoEffCorr[i]/(trackmultiplicity[i]);
    //cout<<"min et corr cb "<<i<<" "<< kaonCorrNoEffCorr[i]<<" +/- "<<kaonErrorNoEffCorr[i]<<endl;
    //cout<<endl;
  }
  delete graphKaonCorrectionShort2;




}
void CalculateHadronicCorrectionForOneBin(Int_t centbin1, Int_t centbin2, Bool_t isPhos, Bool_t isOver500MeV, Float_t &correction, Float_t &error, Bool_t effCorr, Bool_t writeTable){
  //cout<<"cb "<<centbin1<<" - "<<centbin2;
  //cout<<"entries "<<fHistMatchedTracksEvspTvsCentEffCorr->GetEntries()<<endl;
  fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(centbin1,centbin2);
  fHistPeripheralMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(centbin1,centbin2);
  fHistMatchedTracksEvspTvsCentEffCorr500MeV->GetZaxis()->SetRange(centbin1,centbin2);
  fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(centbin1,centbin2);
  TH1D *dataEffCorrTmp = NULL;
  TH1D *dataEffCorrTmp2 = NULL;
  TH1D *dataEffCorrPeripheralTmp = (TH1D*)fHistPeripheralMatchedTracksEvspTvsCentEffCorr->Project3D("y");
  float myDataEffCorrFromPeripheral = dataEffCorrPeripheralTmp->GetMean();
  TH1D *dataTmp = NULL;
  TH1D *foundTmp = NULL;
  TH1D *notfoundTmp = NULL;
  dataTmp = (TH1D*)fHistMatchedTracksEvspTvsCent->Project3D("y");
  dataTmp->SetName(Form("dataTmp%i",centbin1));
  if(isOver500MeV){
    dataEffCorrTmp =(TH1D*) fHistMatchedTracksEvspTvsCentEffCorr500MeV->Project3D("y");
    dataEffCorrTmp2 =(TH1D*) fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
    dataEffCorrTmp2->SetName("dataEffCorrNotOver500");
    dataEffCorrTmp->SetName("dataEffCorrOver500");
    foundTmp = fHistFoundHadronsvsCent500MeV->ProjectionX(Form("Found%iTmp",centbin1),centbin1,centbin2);
    notfoundTmp = fHistNotFoundHadronsvsCent500MeV->ProjectionX(Form("NotFound%iTmp",centbin1),centbin1,centbin2);
  }
  else{
    if(effCorr){            
      dataEffCorrTmp = (TH1D*)fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
      dataEffCorrTmp->SetName("dataEffCorr");
    }
    else{
      dataEffCorrTmp = (TH1D*)fHistMatchedTracksEvspTvsCent->Project3D("y");
      dataEffCorrTmp->SetName("dataNoEffCorr");
    }
    dataEffCorrTmp2 = dataEffCorrTmp;
    //cout<<" Using "<<dataEffCorrTmp2->GetName()<<" entries "<<dataEffCorrTmp2->GetEntries()<<endl;
    foundTmp = fHistFoundHadronsvsCent->ProjectionX(Form("Found%iTmp",centbin1),centbin1,centbin2);
    notfoundTmp = fHistNotFoundHadronsvsCent->ProjectionX(Form("NotFound%iTmp",centbin1),centbin1,centbin2);
  }
  float nfound = foundTmp->GetMean();//fHistFoundHadronsvsCent->GetBinContent(bin);
  float nnotfound = notfoundTmp->GetMean();//fHistNotFoundHadronsvsCent->GetBinContent(bin);
  //cout<<" nfound "<<nfound<<" nnotfound "<<nnotfound<<" total "<<nfound+nnotfound;
  float scaleLow = 0;
  float scaleHigh = 0;
  float scale1 = 1.0;
  float scale2 = 1.0;
  Float_t effectiveEfficiencyThis = 1;
  Float_t effectiveEfficiencyReference = 1;
  if(centbin1>=refBin){//for peripheral just rescale 
    scaleHigh = 1.01;
    scaleLow = 0.99;
  }
  else{
    //float refData = ((TH1D*)data[refBin])->GetMean();
    float refData =dataRefBin->GetMean();
    float myData = ((TH1D*)dataTmp)->GetMean();
    //float refDataEffCorr = ((TH1D*)dataEffCorr[refBin])->GetMean();
    //cout<<" ref bin n entries "<<((TH1D*)dataEffCorr[refBin])->GetEntries()<<", "<<dataEffCorrRefBin->GetEntries()<<" ";
    float refDataEffCorr = dataEffCorrRefBin->GetMean();
    float myDataEffCorr = ((TH1D*)dataEffCorrTmp2)->GetMean();
    //cout<<"ranges "<<refDataEffCorr<<", "<<myDataEffCorr<<", "<<myDataEffCorrFromPeripheral;//<<endl;
    effectiveEfficiencyThis = myData/myDataEffCorr;
    effectiveEfficiencyReference = refData/refDataEffCorr;
    Float_t scale3 = 0;
    //cout<<"efficiency this "<<effectiveEfficiencyThis<<" efficiency ref "<<effectiveEfficiencyReference;
    //cout<<"data Eff corr "<<myDataEffCorr<<" data eff corr ref "<<refDataEffCorr<<" data "<<myData<<" dataRef "<<refData<<endl;
    if(TMath::Abs(myData)>1e-5) scale1 = refData/myData;//scale without efficiency correction -> weights peripheral data by average efficiency of central bin
    if(TMath::Abs(myDataEffCorr)>1e-5){
      //dataEffCorrRefBin->GetMean()/dataEffCorrTmp->GetMean()* dataEffCorrTmp->GetMean();
      scale2 = refDataEffCorr/myDataEffCorr;//scale with efficiency correction ->actual data range
      scale3 = myDataEffCorrFromPeripheral/myDataEffCorr;
    }
    //cout<<" scale1 (uncorr ref) "<<scale1<<" scale2 (ref) "<<scale2<<" scale3 "<<scale3;
    //       if(scale1<scale2){
    // 	scaleLow = scale1;
    // 	scaleHigh = scale2;
    //       }
    //       else{
    // 	scaleLow = scale2;
    // 	scaleHigh = scale1;
    //       }
    if(scale3<scale2){
      scaleLow = scale3;
      scaleHigh = scale2;
    }
    else{
	scaleLow = scale2;
	scaleHigh = scale3;
    }
    //cout<<" scale low "<<scaleLow<<" scale High "<<scaleHigh<<" averge scale "<<(scaleLow+scaleHigh)/2.0;//<<" ref eff "<<effectiveEfficiencyReference<<" this ref "<<effectiveEfficiencyThis<<" refeff/thiseff*scale2 "<<effectiveEfficiencyReference/effectiveEfficiencyThis<<" this/ref "<<effectiveEfficiencyThis/effectiveEfficiencyReference<<endl;
    //cout<<"scale 1 "<<scale1<<" scale 2 "<<scale2<<endl;
    
  }
    //For EMCal the reference efficiency is lower than the central efficiency because of energy resolution
    //For PHOS the reference efficiency is higher than the central efficiency
    //The average background cluster has lower energy than the average signal cluster
    //For the EMCal the dominant effect is bad track matches.  For that reason the energy deposited in the current bin is an upper bound.
    //For the PHOS the dominant effect is the efficiency.  If the true cluster energy is lower than that observed, the estimate will overestimate the correction.
    if(isPhos){//for the PHOS - leave it alone, this is already a pretty good estimate.  One bound is the reference data (scale2) and the other is the current bin but taking into account that the efficiencies are different in the two bins (scale1).
      //scaleLow = scale2;
      //scaleHigh = scale1;
    }
    else{//For the EMCal the range is the reference range up to the current bin
      //scaleLow = scale1;
      //scaleHigh = 1;
    }
    //cout<<" scaleLow "<<scaleLow<<" scaleHigh "<<scaleHigh;
    //    if(averageHadEnergy<0){//if this hasn't been filled yet
      Float_t low = 100;
      Float_t high = -1;
      
//       for(int i=refBinHigh;i<=refBin;i++){
// 	Float_t val = ((TH1D*)dataEffCorr[i])->GetMean();
// 	//cout<<" i "<<val;
// 	if(val<low) low = val;
// 	if(val>high) high = val;
//       }
      //cout<<endl;
      //}
    float myavg = dataEffCorrTmp->GetMean();
    //if(myavg>high && !isPhos) high = myavg;
    averageHadEnergy = (low+high)/2.0;
    averageHadEnergyError = (high-low)/2.0;
    //cout<<" AVERAGE HAD ENERGY "<<averageHadEnergy<<"+/-"<<averageHadEnergyError;
    //this assumes that the reference bin is right

//     if(!isPhos){
//       scaleLow = 0.90;
//       scaleHigh = 1.0;
//     }
//    cout<<" ranges "<<scaleLow*myavg<<", "<<scaleHigh*myavg<<", "<<myDataEffCorrFromPeripheral<<endl;
    float avg = (scaleLow+scaleHigh)/2.0*myavg;
    float err = TMath::Abs((scaleLow-scaleHigh))/2.0*myavg;
//     float avg = (scaleHigh*myavg+myDataEffCorrFromPeripheral)/2.0;
//     float err = TMath::Abs((scaleHigh*myavg-myDataEffCorrFromPeripheral)/2.0);
    //cout<<" Nominal value "<<avg<<"+/-"<<err<<" this bin "<<myavg<<endl;
//     avg = averageHadEnergy;
//     err = averageHadEnergyError;
    // if(isPhos){//for the PHOS, false matches are actually rare, so we will use the value from this bin
    //averageHadEnergy = (0.97*low+high)/2.0;
      //averageHadEnergyError = (high-low*0.97)/2.0;
      // }
    //cout<<" avg "<<avg<<" +/- "<<err;
    if(TMath::Abs(avg)<1e-3){
      avg = 1e-3;
      cerr<<"WARNING:  ERROR NOT CALCULATED CORRECTLY!!"<<endl;//prevents a crash
    }
    //factor is the fraction to reduce the track-matched ET by to get the true background ET
    //corrfac is the factor to multiply by in order to get the fraction of hadrons leaving deposits which come from low pT
    float percentEfficiencyError = 0.05;
    //    float  factor = 1-0.04;
        float  factor = 1-0.03;
//         float corrfac = 1.275-1;
//     float corrfacerr = 0.059 ;
//         float corrfac = 0.208938;
// 	float corrfacerr = 0.0357719 ;
	//corrfac = 0.183584;
	//corrfacerr = 0.0393219;
// 	corrfac = 0.323685;
// 	corrfacerr = 0.0760131;
    float eLowAverage = avg;
    float eLowAverageErr = err;
    cout<<"eLowAverage "<<avg<<" +/- "<<err<<" % "<<err/avg*100.0<<endl;
    if(isPhos){
      factor = 1-0.02;
      //factor = 1-0.03;
//       corrfac = 0.183584;
//       corrfacerr = 0.0393219;
//       corrfac = 0.169776;
//       corrfacerr = 0.0468306;
//       corrfac = 1.300-1;
//       corrfacerr = 0.065;
    }
    if(isOver500MeV){
      eLowAverage = 1.0;
      eLowAverageErr = 0.05;
      if(isPhos){
	//fraction ranges from 5% - 26%
	corrfac = ((1.0/0.95 + 1.0/0.74)/2.0 - 1)*eLowAverage/avg;
	float corrfacerrtmp = (TMath::Abs((1.0/0.95 - 1.0/0.74)/2.0))*eLowAverage/avg;
	corrfacerr = corrfac * TMath::Sqrt(TMath::Power(corrfacerrtmp/corrfac,2)+TMath::Power(eLowAverageErr/eLowAverage,2));
      }
      else{
	//fraction ranges from 5% - 17%
	corrfac = ((1.0/0.95 + 1.0/0.83)/2.0 - 1)*eLowAverage/avg;
	float corrfacerrtmp = (TMath::Abs((1.0/0.95 - 1.0/0.83)/2.0))*eLowAverage/avg;
	corrfacerr = corrfac * TMath::Sqrt(TMath::Power(corrfacerrtmp/corrfac,2)+TMath::Power(eLowAverageErr/eLowAverage,2));
      }
    }
    // cout<<" nfound "<<nfound<<" nnotfound "<<nnotfound<<" total n "<<nfound+nnotfound<<" eff "<<nfound/(nfound+nnotfound);
    //the energy from low pT is the fraction of tracks that come from low pT tracks times the average energy of these tracks
    float eLow = corrfac * (nfound+nnotfound)*eLowAverage;
    float eLowErr = TMath::Sqrt(TMath::Power(corrfacerr*(nfound+nnotfound)*eLowAverage,2)+TMath::Power(eLowAverageErr*corrfac* (nfound+nnotfound),2)+TMath::Power(eLow*percentEfficiencyError,2));//error on the hadronic correction

    float eNotFound = nnotfound*avg;
    float eNotFoundErr = TMath::Sqrt(TMath::Power(err*nnotfound,2)+TMath::Power(percentEfficiencyError*eNotFound,2));//error on the hadronic correction

    float y = (corrfac * (nfound+nnotfound) +nnotfound)*avg;
    float finalerr = TMath::Sqrt(TMath::Power(corrfacerr*(nfound+nnotfound)*avg,2)+err*err*(TMath::Power(corrfac* (nfound+nnotfound),2)+nnotfound*nnotfound)+TMath::Power(percentEfficiencyError*y,2));//error on the hadronic correction
    correction = y;
    error = finalerr;
    //cout<<"error "<<finalerr/y<<endl;
    //error = 0.0;
    //cout<<" corr fac "<<corrfac<<" +/- "<<corrfacerr;
    //cout<<" AVERAGE eLow "<<eLow<<" +/- "<<eLowErr<<" % "<<eLowErr/y*100.0<<") eHigh "<<eNotFound<<" +/- "<<eNotFoundErr<<"("<<eNotFoundErr/y<<") total "<<y<<" +/- "<<finalerr<<"("<<finalerr/y <<")"<<endl;//<<" frac low "<<eLow/y<<endl;
    //cout<<endl;
    //correction = N(notfound)*<Enotfound> + fnotfound*Etotal
    //centrality, Nnotfound, <Enotfound>, fnotcound, Etotal, correction
    if(writeTable){
      myfileHadCorrTable<<Form("%i-%i",(centbin1-1)*5,centbin2*5)<<"\\% & ";
      myfileHadCorrTable<<Form("%2.1f $\\pm$ %2.1f",nnotfound,percentEfficiencyError*nnotfound)<<" & ";
      myfileHadCorrTable<<Form("%2.3f $\\pm$ %2.3f",avg,err)<<" & ";
      myfileHadCorrTable<<Form("%2.3f $\\pm$ %2.3f",corrfac,corrfacerr)<<" & ";
      myfileHadCorrTable<<Form("%2.1f $\\pm$ %2.1f",(nfound+nnotfound)*avg,(nfound+nnotfound)*avg*TMath::Sqrt(percentEfficiencyError*percentEfficiencyError+err*err/avg/avg))<<" & ";
      myfileHadCorrTable<<Form("%2.1f $\\pm$ %2.1f",correction,error)<<" \\\\ "<<endl;
    }


    //myfileHadCorrTable<<Form("%i-%i",(centbin1-1)*5,centbin2*5)<<"\\% & "<<Form("%2.3f",((TH1D *)data[centbin1])->GetMean())<<" & "<<Form("%2.3f",((TH1D *)dataEffCorr[centbin1])->GetMean())<<Form(" &  %2.3f $\\pm$ %2.3f",avg,err);
    //myfileHadCorrTable<<" & "<< Form("%2.1f $\\pm$ %2.1f",eLow,eLowErr);
    //myfileHadCorrTable<<" & "<< Form("%2.2f $\\pm$ %2.2f",eNotFound,eNotFoundErr);
    //myfileHadCorrTable<<"& "<< Form("%2.2f $\\pm$ %2.2f",y,finalerr);
    //myfileHadCorrTable<<" & "<< Form("%2.3f $\\pm$ %2.3f",y/npart[centbin1],finalerr/npart[centbin1]) <<"\\\\"<<endl;

    delete dataEffCorrTmp;
    delete foundTmp;
    delete notfoundTmp;
    delete dataEffCorrPeripheralTmp;
}
void CalculateHadronCorrections(Bool_t isPhos){
  
  TString texfilename = "/tmp/energydepositstable"+detector+".tex";
  myfileHadCorrTable.open (texfilename.Data());
  float plotscale = 5.0;
  for(int i=0;i<19;i++){
    Float_t correction = 0;
    Float_t error = 0;//not above 500 GeV, with eff corr
    CalculateHadronicCorrectionForOneBin(i+1,i+1,isPhos,kFALSE,correction,error,kTRUE,kTRUE);
    hadCorr[i] = correction;//hadCorrEmcal[i];
    hadError[i] = error;//hadErrorEmcal[i];
    hadronCorrPerNCh[i] = correction/(trackmultiplicity[i])/plotscale;//hadCorrEmcal[i];
    hadronErrorPerNCh[i] = error/(trackmultiplicity[i])/plotscale;//hadErrorEmcal[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicity[i])<< " +/- "<<  error/(trackmultiplicity[i]) <<endl;
  }

  int j=0;
  for(int i=0;i<10;i++){
    int centbinlow = i+1;
    int centbinhigh = i+1;
    if(i<2){//These bins are exactly what they should bin in the 20 bin binning
      j++;//i=0 j=0; i=1 j=1
      centbinlow = j;
      centbinhigh = j;
    }
    else{
      centbinlow = j+1;
      centbinhigh = j+2;
      j+=2;
    }
    Float_t correction = 0;
    Float_t error = 0;//not above 500 GeV, with eff corr, 10 bins
    CalculateHadronicCorrectionForOneBin(centbinlow,centbinhigh,isPhos,kFALSE,correction,error,kTRUE,kFALSE);
    hadCorrShort[i] = correction;//hadCorrEmcal[i];
    hadErrorShort[i] = error;//hadErrorEmcal[i];
    hadronCorrPerNChShort[i] = correction/(trackmultiplicityShort[i])/plotscale;//hadCorrEmcal[i];
    hadronErrorPerNChShort[i] = error/(trackmultiplicityShort[i])/plotscale;//hadErrorEmcal[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicityShort[i])<< " +/- "<<  error/(trackmultiplicityShort[i]) <<endl;
  }
  //No eff corr
  for(int i=0;i<19;i++){
    Float_t correction = 0;
    Float_t error = 0;//not above 500 GeV, without eff corr
    CalculateHadronicCorrectionForOneBin(i+1,i+1,isPhos,kFALSE,correction,error,kFALSE,kFALSE);
    hadCorrNoEffCorr[i] = correction;//hadCorrEmcalNoEffCorr[i];
    hadErrorNoEffCorr[i] = error;//hadErrorEmcalNoEffCorr[i];
    hadronCorrPerNChNoEffCorr[i] = correction/(trackmultiplicity[i])/plotscale;//hadCorrEmcalNoEffCorr[i];
    hadronErrorPerNChNoEffCorr[i] = error/(trackmultiplicity[i])/plotscale;//hadErrorEmcalNoEffCorr[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicityNoEffCorr[i])<< " +/- "<<  error/(trackmultiplicity[i]) <<endl;
  }

  int j=0;
  for(int i=0;i<10;i++){
    int centbinlow = i+1;
    int centbinhigh = i+1;
    if(i<2){//These bins are exactly what they should bin in the 20 bin binning
      j++;//i=0 j=0; i=1 j=1
      centbinlow = j;
      centbinhigh = j;
    }
    else{
      centbinlow = j+1;
      centbinhigh = j+2;
      j+=2;
    }
    Float_t correction = 0;
    Float_t error = 0;//not above 500 GeV, without eff corr, 10 bins
    CalculateHadronicCorrectionForOneBin(centbinlow,centbinhigh,isPhos,kFALSE,correction,error,kFALSE,kFALSE);
    hadCorrShortNoEffCorr[i] = correction;//hadCorrEmcalNoEffCorr[i];
    hadErrorShortNoEffCorr[i] = error;//hadErrorEmcalNoEffCorr[i];
    hadronCorrPerNChShortNoEffCorr[i] = correction/(trackmultiplicityShort[i])/plotscale;//hadCorrEmcalNoEffCorr[i];
    hadronErrorPerNChShortNoEffCorr[i] = error/(trackmultiplicityShort[i])/plotscale;//hadErrorEmcalNoEffCorr[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicityShortNoEffCorr[i])<< " +/- "<<  error/(trackmultiplicityShort[i]) <<endl;
  }

  myfileHadCorrTable.close();

}

TGraphErrors *GetSignalFractionGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,signalFraction,trackmultiplicityError,signalFractionError);
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,minEtCorrShort,npartErrShort,minEtErrorShort);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}

TGraphErrors *GetMinEtCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,minEtCorrShort,trackmultiplicityShortError,minEtErrorShort);
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,minEtCorrShort,npartErrShort,minEtErrorShort);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetMinEtCorrectionGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,minEtCorr,trackmultiplicityError,minEtError);
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,minEtCorrShort,npartErrShort,minEtErrorShort);
  SetStyles(gr3,30,TColor::kGreen+3);
    return gr3;

}


TGraphErrors *GetSecondaryCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,secondaryCorrPerNChShort,trackmultiplicityShortError,secondaryErrorPerNChShort);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetSecondaryCorrectionGraph(){
  TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,secondaryCorrPerNCh,trackmultiplicityError,secondaryErrorPerNCh);
  SetStyles(gr3,30,TColor::kGreen+3);
    return gr3;

}

TGraphErrors *GetSecondaryCorrectionGraphNoEffCorr(){
  TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,secondaryCorrPerNChNoEffCorr,trackmultiplicityError,secondaryErrorPerNChNoEffCorr);
  SetStyles(gr3,30,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetSecondaryCorrectionGraphShortNoEffCorr(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,secondaryCorrPerNChShortNoEffCorr,trackmultiplicityShortError,secondaryErrorPerNChShortNoEffCorr);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraph(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,neutronCorrPerNCh,trackmultiplicityError,neutronErrorPerNCh);
    SetStyles(gr3,24,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,npartErrShort,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,neutronCorrPerNChShort,trackmultiplicityShortError,neutronErrorPerNChShort);
    SetStyles(gr3,20,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraphNoEffCorr(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,neutronCorrPerNChNoEffCorr,trackmultiplicityError,neutronErrorPerNChNoEffCorr);
    SetStyles(gr3,24,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraphShortNoEffCorr(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,npartErrShort,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,neutronCorrPerNChShortNoEffCorr,trackmultiplicityShortError,neutronErrorPerNChShortNoEffCorr);
    SetStyles(gr3,20,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraph(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,hadronCorrPerNCh,trackmultiplicityError,hadronErrorPerNCh);
    SetStyles(gr3,25,1);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,xpionerr,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,hadronCorrPerNChShort,trackmultiplicityShortError,hadronErrorPerNChShort);
    SetStyles(gr3,21,1);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraphNoEffCorr(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,hadronCorrPerNChNoEffCorr,trackmultiplicityError,hadronErrorPerNChNoEffCorr);
    SetStyles(gr3,25,1);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraphShortNoEffCorr(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,xpionerr,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,hadronCorrPerNChShortNoEffCorr,trackmultiplicityShortError,hadronErrorPerNChShortNoEffCorr);
    SetStyles(gr3,21,1);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonCorrPerNChShort,trackmultiplicityShortError,kaonErrorPerNChShort);
    SetStyles(gr3,33,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,kaonCorrPerNCh,trackmultiplicityError,kaonErrorPerNCh);
    SetStyles(gr3,27,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraphShortNoEffCorr(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonCorrPerNChShortNoEffCorr,trackmultiplicityShortError,kaonErrorPerNChShortNoEffCorr);
    SetStyles(gr3,33,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraphNoEffCorr(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,kaonCorrPerNChNoEffCorr,trackmultiplicityError,kaonErrorPerNChNoEffCorr);
    SetStyles(gr3,27,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonYieldPerNCh,trackmultiplicityShortError,kaonYieldPerNChErr);
    SetStyles(gr3,33,TColor::kBlue);
    return gr3;
}
TGraphErrors *GetKaonEtGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);
    SetStyles(gr3,27,TColor::kBlue);
    return gr3;
}


//=====================READ IN DATA===================================
Float_t arrayofzeros[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityShortError[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityShortError[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totaltrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totaltrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtStatError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtStatErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairValuesNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairValuesShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairErrorNoEffCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totalCorrectionPerNPartPairErrorShortNoEffCorr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

Float_t rawEtNoEffCorrValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNPartValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNPartError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNPartValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNPartError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNPartValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNPartError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

Float_t corrEtValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtStatError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtStatErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totEtValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totEtValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtValuesFormulaC[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtErrorFormulaC[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValuesFormulaC[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairErrorFormulaC[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtValuesFormulaB[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtErrorFormulaB[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValuesFormulaB[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairErrorFormulaB[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

TH2F *fHistNominalRawEt;
TH2F *fHistNominalNonLinLowEt;
TH2F *fHistNominalNonLinHighEt;
TH2F *fHistNominalEffLowEt;
TH2F *fHistNominalEffHighEt;
TH2F *fHistTotRawEt;
TH2F *fHistTotAllRawEt;
TH2F *fHistTotRawEtNoEffCorr;
TH2F *fHistTotAllRawEtNoEffCorr;
TH2F *fHistTotRawEt500MeV;
TH3F  *fHistMatchedTracksEvspTvsCent;
TH3F  *fHistMatchedTracksEvspTvsCentEffCorr;
TH3F  *fHistPeripheralMatchedTracksEvspTvsCentEffCorr;
TH3F  *fHistMatchedTracksEvspTvsCentEffCorr500MeV;
TH2F *fHistFoundHadronsvsCent;
TH2F *fHistNotFoundHadronsvsCent;
TH2F *fHistFoundHadronsvsCent500MeV;
//old default: bin 18 = centrality bin 17
//ref bins in compiled code:  centrality bins 11-16
Int_t refBin = 12;//Reference bin for scaling
Int_t refBinHigh = 17;//Reference bin for scaling
Float_t corrfac = 0;//fraction of hadronic ET which is from low pT tracks
Float_t corrfacerr = 0;
TObjArray data(21);
TObjArray dataFits(21);
TObjArray dataEffCorr(21);
TObjArray dataEffCorr500MeV(21);
TH1D *dataRefBin = NULL;
TH1D *dataEffCorrRefBin = NULL;
TObjArray rawEt(21);
TObjArray rawEtShort(11);
TObjArray rawEtNoEffCorr(21);
TObjArray rawEtAll(21);
TObjArray rawEtAllNoEffCorr(21);

void ReadInData(char *filename,TString det){
  cout<<"Reading in data..."<<endl;
  Bool_t isPhos = kTRUE;
  if(det.Contains("Emc")){isPhos=kFALSE;}
  else{
    refBin = 12;
    cout<<"Decided this is PHOS and am using reference bin 12"<<endl;
  }

  TFile *f = TFile::Open(filename, "READ");
    if (!f)
    {
        std::cerr << "Could not open file: " << filename << std::endl;
    }

    TList *l = dynamic_cast<TList*>(f->Get("out1"));
    if (!l)
    {
        std::cerr << "Could not get object list from: " << filename << std::endl;
    }
    TString prefix = "fHistNominal";
    fHistNominalRawEt =(TH2F*) l->FindObject((prefix+"RawEt"+det+"Rec").Data());
    fHistNominalNonLinLowEt = (TH2F*)l->FindObject((prefix+"NonLinLowEt"+det+"Rec").Data());
    fHistNominalNonLinHighEt = (TH2F*)l->FindObject((prefix+"NonLinHighEt"+det+"Rec").Data());
    fHistNominalEffLowEt = (TH2F*)l->FindObject((prefix+"EffLowEt"+det+"Rec").Data());
    fHistNominalEffHighEt = (TH2F*)l->FindObject((prefix+"EffHighEt"+det+"Rec").Data());
    fHistTotRawEt =(TH2F*) l->FindObject("fHistTotRawEtEffCorr");
    fHistTotAllRawEt = (TH2F*) l->FindObject("fHistTotAllRawEtEffCorr");
    fHistTotRawEtNoEffCorr = (TH2F*)l->FindObject("fHistTotRawEt");
    fHistTotAllRawEtNoEffCorr = (TH2F*)l->FindObject("fHistTotAllRawEt");
    fHistTotRawEt500MeV =(TH2F*) l->FindObject("fHistTotRawEtEffCorr500MeV");

    fHistCentVsNchVsNcl =(TH3F*) l->FindObject("fHistCentVsNchVsNclReco");


    fHistMatchedTracksEvspTvsCent =(TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCent");
    fHistMatchedTracksEvspTvsCentEffCorr =(TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr");
    fHistPeripheralMatchedTracksEvspTvsCentEffCorr =(TH3F*) l->FindObject("fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr");
    fHistMatchedTracksEvspTvsCentEffCorr500MeV = (TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr500MeV");
    fHistFoundHadronsvsCent = (TH2F*)l->FindObject("fHistFoundHadronsvsCent");
    fHistNotFoundHadronsvsCent = (TH2F*)l->FindObject("fHistNotFoundHadronsvsCent");
    fHistFoundHadronsvsCent500MeV = (TH2F*)l->FindObject("fHistFoundHadronsvsCent500MeV");
    fHistNotFoundHadronsvsCent500MeV = (TH2F*)l->FindObject("fHistNotFoundHadronsvsCent500MeV");
    fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(refBin,refBinHigh);
    dataRefBin =  (TH1D*) fHistMatchedTracksEvspTvsCent->Project3D("y");
    dataRefBin->SetName("dataRefBin");
    fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(refBin,refBinHigh);
    dataEffCorrRefBin =(TH1D*) fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
    dataEffCorrRefBin->SetName("dataEffCorrRefBin");
    cout<<"Using reference centrality bins "<<refBin-1<<" - "<<refBinHigh-1<<endl;
    int nbins = 20;
    for(int bin = 1; bin<=nbins;bin++){
      fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(bin,bin);
      data[bin] = fHistMatchedTracksEvspTvsCent->Project3D("y");
      fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(bin,bin);
      dataEffCorr[bin] = fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
      fHistMatchedTracksEvspTvsCentEffCorr500MeV->GetZaxis()->SetRange(bin,bin);
      dataEffCorr500MeV[bin] = fHistMatchedTracksEvspTvsCentEffCorr500MeV->Project3D("y");
      ((TH1D*)data[bin])->SetName(Form("DataEff%i",bin));
      ((TH1D*)dataEffCorr[bin])->SetName(Form("DataEffCorr%i",bin));
      Int_t nentries = ((TH1D*)dataEffCorr[bin])->GetEntries();
      cout<<"centbin "<<bin<<" number of matched tracks "<<nentries<<endl;
      ((TH1D*)dataEffCorr500MeV[bin])->SetName(Form("DataEffCorr500MeV%i",bin));

      rawEt[bin]= fHistTotRawEt->ProjectionX(Form("RawEt%i",bin),bin,bin);
      ((TH1D*)rawEt[bin])->SetName(Form("rawEt%i",bin));
      partialCorrEtValues[bin-1] = (Float_t)((TH1D*)rawEt[bin])->GetMean();
      partialCorrEtError[bin-1] = (Float_t) ((TH1D*)rawEt[bin])->GetMeanError();
      partialCorrEtPerNPartPairValues[bin-1] = partialCorrEtValues[bin-1]/(npart[bin-1])/2.0*10;
      partialCorrEtPerNPartPairError[bin-1]  =  partialCorrEtError[bin-1]/(npart[bin-1])/2.0*10;

      rawEtNoEffCorr[bin]= fHistTotRawEtNoEffCorr->ProjectionX(Form("RawEtNoEffCorr%i",bin),bin,bin);
      ((TH1D*)rawEtNoEffCorr[bin])->SetName(Form("rawEtNoEffCorr%i",bin));
      rawEtNoEffCorrValues[bin-1] = (Float_t)((TH1D*)rawEtNoEffCorr[bin])->GetMean();
      rawEtNoEffCorrError[bin-1] = (Float_t) ((TH1D*)rawEtNoEffCorr[bin])->GetMeanError();

      rawEtAllNoEffCorr[bin]= fHistTotAllRawEtNoEffCorr->ProjectionX(Form("RawEtAllNoEffCorr%i",bin),bin,bin);
      ((TH1D*)rawEtAllNoEffCorr[bin])->SetName(Form("rawEtAllNoEffCorr%i",bin));
      rawEtAllNoEffCorrValues[bin-1] = (Float_t)((TH1D*)rawEtAllNoEffCorr[bin])->GetMean();
      rawEtAllNoEffCorrError[bin-1] = (Float_t) ((TH1D*)rawEtAllNoEffCorr[bin])->GetMeanError();

      rawEtAll[bin]= fHistTotAllRawEt->ProjectionX(Form("RawEtAll%i",bin),bin,bin);
      ((TH1D*)rawEtAll[bin])->SetName(Form("rawEtAll%i",bin));
      rawEtAllValues[bin-1] = (Float_t)((TH1D*)rawEtAll[bin])->GetMean();
      rawEtAllError[bin-1] = (Float_t) ((TH1D*)rawEtAll[bin])->GetMeanError();
      //cout<<"bin "<<bin<<" "<<partialCorrEtValues[bin-1]<<" "<<rawEtNoEffCorrValues[bin-1]<<" "<<rawEtAllNoEffCorrValues[bin-1]<<" "<<rawEtAllValues[bin-1]<<endl;
      //cout<<"bin "<<bin<<" eff corr "<<partialCorrEtValues[bin-1]<<" no eff corr "<<rawEtNoEffCorrValues[bin-1]<<" eff "<< rawEtNoEffCorrValues[bin-1]/partialCorrEtValues[bin-1] <<endl;
      if(partialCorrEtValues[bin-1]>0) averageEfficiency[bin-1] = rawEtNoEffCorrValues[bin-1]/partialCorrEtValues[bin-1];

      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",bin,bin);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",bin,bin);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",bin,bin);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",bin,bin);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",bin,bin);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = 0;
      if(nonlinhigh >0 || nonlinlow >0) nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
      float efffracerr = 0;
      if(effhigh >0 || efflow>0)efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      //cout<<"cb "<<bin-1<<" nonlinerr "<<nonlinfracerr<<" efficiencyError "<<efffracerr<<endl;
      nonLinError[bin-1] = nonlinfracerr;
      efficiencyError[bin-1] = efffracerr;

      //(not)matchedtrackmultiplicity, (not)matchedtrackmultiplicityPerNCh, (not)matchedtrackmultiplicityPerNCl
      temp = fHistFoundHadronsvsCent->ProjectionX(Form("Found%iTmp",bin),bin,bin);
      matchedtrackmultiplicity[bin-1] = temp->GetMean();
      delete temp;
      temp = fHistNotFoundHadronsvsCent->ProjectionX(Form("NotFound%iTmp",bin),bin,bin);
      notmatchedtrackmultiplicity[bin-1] = temp->GetMean();
      delete temp;



    }


    for(int cb=0;cb<20;cb++){
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb+1,cb+1);
      TH1D *trackmultiplicityHist = fHistCentVsNchVsNcl->Project3D("y");
      TH1D *clustermultiplicityHist = fHistCentVsNchVsNcl->Project3D("z");
      trackmultiplicity[cb] = (Float_t) trackmultiplicityHist->GetMean();
      clustermultiplicity[cb] = (Float_t)clustermultiplicityHist->GetMean();
      trackmultiplicityError[cb] = (Float_t) trackmultiplicityHist->GetMeanError();
      clustermultiplicityError[cb] = (Float_t)clustermultiplicityHist->GetMeanError();
      delete trackmultiplicityHist;
      delete clustermultiplicityHist;
    }
    int cb1 = 0;
    for(int cb=0;cb<10;cb++){
      int cb2 = cb1+1;
      if(cb1<2) cb2 = cb1;
      //cout<<"From "<<cb1<<" to "<<cb2<<endl;
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb1+1,cb2+1);
      TH1D *trackmultiplicityHistShort = fHistCentVsNchVsNcl->Project3D("y");
      TH1D *clustermultiplicityHistShort = fHistCentVsNchVsNcl->Project3D("z");
      trackmultiplicityShort[cb] = (Float_t) trackmultiplicityHistShort->GetMean();
      clustermultiplicityShort[cb] = (Float_t)clustermultiplicityHistShort->GetMean();
      trackmultiplicityShortError[cb] = (Float_t) trackmultiplicityHistShort->GetMeanError();
      clustermultiplicityShortError[cb] = (Float_t)clustermultiplicityHistShort->GetMeanError();
      delete trackmultiplicityHistShort;
      delete clustermultiplicityHistShort;
      if(cb1<2) cb1++;
      else{cb1+=2;}
      rawEtShort[cb]= fHistTotRawEt->ProjectionX(Form("RawEtShort%i",bin),cb1+1,cb2+1);
      ((TH1D*)rawEtShort[cb])->SetName(Form("rawEtShort%i",cb));
      partialCorrEtValuesShort[cb] = (Float_t)((TH1D*)rawEtShort[cb])->GetMean();
      partialCorrEtErrorShort[cb] = (Float_t) ((TH1D*)rawEtShort[cb])->GetMeanError();

      partialCorrEtPerNChValuesShort[cb] = partialCorrEtValuesShort[cb]/(trackmultiplicityShort[cb])/2.0;
      partialCorrEtPerNChErrorShort[cb]  =  partialCorrEtErrorShort[cb]/(trackmultiplicityShort[cb])/2.0;

      partialCorrEtPerNPartValuesShort[cb] = partialCorrEtValuesShort[cb]/(npart[cb])/2.0;
      partialCorrEtPerNPartErrorShort[cb]  =  partialCorrEtErrorShort[cb]/(npart[cb])/2.0;

      partialCorrEtPerNPartPairValuesShort[cb] = partialCorrEtValues[cb]/(npart[cb])/2.0*10;
      partialCorrEtPerNPartPairErrorShort[cb]  =  partialCorrEtError[cb]/(npart[cb])/2.0*10;

      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",cb1+1,cb2+1);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",cb1+1,cb2+1);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",cb1+1,cb2+1);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",cb1+1,cb2+1);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",cb1+1,cb2+1);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = 0;
      if(nonlinhigh >0 || nonlinlow >0) nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
//       if(isPhos){ 
// 	nonlinfracerr = 0.005;
//       }
      float efffracerr = 0;
      if(effhigh >0 || efflow>0)efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      nonLinErrorShort[cb] = nonlinfracerr;
      if(isPhos){
	efficiencyErrorShort[cb] = 0.005;
      }
      else{
	efficiencyErrorShort[cb] = 0.02;
      }



    }

}
void ApplyCorrections(Float_t scale){//scale takes into account the acceptance in eta and phi and the 1/etaacc
  //correlation of errors:
    //hadCorr - calculated using matched tracks, correcting with tracking efficiency, systematic error dominated by uncertainty in energy deposited in calorimeter
  //kaon correction - calculated using kaon spectra, systematic errors from spectra systematic errors
  //neutron correction - systematic errors a bit of a fudge.
  //secondary correction - systematic errors from Nch vs Ncl scaling
  //etmin - calculated from kinematics, simulation, and pion spectra
  //efficiency error - determined from different material budgets
  //nonlinearity error - difference between test beam data and simulation
  //arguably kaon correction and etmin correction are somewhat correlated

  for(int cb = 0;cb<19;cb++){
    if(trackmultiplicity[cb]>0){
      partialCorrEtPerNChValues[cb] = partialCorrEtValues[cb]/(trackmultiplicity[cb]);
      partialCorrEtPerNChError[cb]  =  partialCorrEtError[cb]/(trackmultiplicity[cb]);
      partialCorrEtPerNPartValues[cb] = partialCorrEtValues[cb]/(npart[cb]);
      partialCorrEtPerNPartError[cb]  =  partialCorrEtError[cb]/(npart[cb]);
      rawEtNoEffCorrPerNChValues[cb] = rawEtNoEffCorrValues[cb]/(trackmultiplicity[cb]);
      rawEtNoEffCorrPerNChError[cb]  =  rawEtNoEffCorrError[cb]/(trackmultiplicity[cb]);
      rawEtAllNoEffCorrPerNChValues[cb] = rawEtAllNoEffCorrValues[cb]/(trackmultiplicity[cb]);
      rawEtAllNoEffCorrPerNChError[cb]  =  rawEtAllNoEffCorrError[cb]/(trackmultiplicity[cb]);
      rawEtAllPerNChValues[cb] = rawEtAllValues[cb]/(trackmultiplicity[cb]);
      rawEtAllPerNChError[cb]  =  rawEtAllError[cb]/(trackmultiplicity[cb]);
      rawEtNoEffCorrPerNPartValues[cb] = rawEtNoEffCorrValues[cb]/(npart[cb]);
      rawEtNoEffCorrPerNPartError[cb]  =  rawEtNoEffCorrError[cb]/(npart[cb]);
      rawEtAllNoEffCorrPerNPartValues[cb] = rawEtAllNoEffCorrValues[cb]/(npart[cb]);
      rawEtAllNoEffCorrPerNPartError[cb]  =  rawEtAllNoEffCorrError[cb]/(npart[cb]);
      rawEtAllPerNPartValues[cb] = rawEtAllValues[cb]/(npart[cb]);
      rawEtAllPerNPartError[cb]  =  rawEtAllError[cb]/(npart[cb]);
      matchedtrackmultiplicityPerNCh[cb] = matchedtrackmultiplicity[cb]/(trackmultiplicity[cb]);
      notmatchedtrackmultiplicityPerNCh[cb] = notmatchedtrackmultiplicity[cb]/(trackmultiplicity[cb]);
    }
    if(clustermultiplicity[cb]>0){
      partialCorrEtPerNClValues[cb] = partialCorrEtValues[cb]/(clustermultiplicity[cb]);
      partialCorrEtPerNClError[cb]  =  partialCorrEtError[cb]/(clustermultiplicity[cb]);
      rawEtNoEffCorrPerNClValues[cb] = rawEtNoEffCorrValues[cb]/(clustermultiplicity[cb]);
      rawEtNoEffCorrPerNClError[cb]  =  rawEtNoEffCorrError[cb]/(clustermultiplicity[cb]);
      rawEtAllNoEffCorrPerNClValues[cb] = rawEtAllNoEffCorrValues[cb]/(clustermultiplicity[cb]);
      rawEtAllNoEffCorrPerNClError[cb]  =  rawEtAllNoEffCorrError[cb]/(clustermultiplicity[cb]);
      rawEtAllPerNClValues[cb] = rawEtAllValues[cb]/(clustermultiplicity[cb]);
      rawEtAllPerNClError[cb]  =  rawEtAllError[cb]/(clustermultiplicity[cb]);
      matchedtrackmultiplicityPerNCl[cb] = matchedtrackmultiplicity[cb]/(clustermultiplicity[cb]);
      notmatchedtrackmultiplicityPerNCl[cb] = notmatchedtrackmultiplicity[cb]/(clustermultiplicity[cb]);
    }

    totaltrackmultiplicityPerNCh[cb] = matchedtrackmultiplicityPerNCh[cb] + notmatchedtrackmultiplicityPerNCh[cb];
    totaltrackmultiplicityPerNCl[cb] = matchedtrackmultiplicityPerNCl[cb] + notmatchedtrackmultiplicityPerNCl[cb];
//     if(scale<4.0){
//       kaonCorr[cb] = 100.0/40.0*kaonCorr[cb];
//       kaonError[cb] = 100.0/40.0*kaonError[cb];
//     }
    //cout<<"cb "<<cb<<" "<<partialCorrEtValues[cb]<<" "<<rawEtNoEffCorrValues[cb]<<" "<<rawEtAllNoEffCorrValues[cb]<<" "<<rawEtAllValues[cb]<<endl;
    //cout<<"cb "<<cb<<" "<<partialCorrEtPerNChValues[cb]<<" "<<rawEtNoEffCorrPerNChValues[cb]<<" "<<rawEtAllNoEffCorrPerNChValues[cb]<<" "<<rawEtAllPerNChValues[cb]<<endl;
      //cout<<"cb "<<cb<<" "<<partialCorrEtPerNChValues[cb]<<" "<< partialCorrEtValues[cb]<<" "<< partialCorrEtError[cb]<<" "<<trackmultiplicity[cb] <<endl;
    corrEtValues[cb] = scale*(partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/minEtCorr[cb];
    if(rawEtAllNoEffCorr[cb]){
      float fracraw = (partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/rawEtAllNoEffCorrValues[cb];
      float fracraw2 = (partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/rawEtAllValues[cb];
      float fracraw3 = fracraw/fem*scale/minEtCorr[cb];//~0.28/0.23*6.5
      //if(cb>9) 
      fracraw3 = 15.0;
      //cout<<"cb "<<cb<<" frac raw "<<fracraw<<" frac raw 2 "<<fracraw2<<" frac raw 3 "<<fracraw3<<" raw et "<<rawEtAllNoEffCorrValues[cb]<<" test "<<fracraw3*rawEtAllValues[cb]<<endl;
    }
    //else{cerr<<"cannot calc fracraw"<<" raw et "<<rawEtAllNoEffCorrValues[cb]<<endl;}
    corrEtStatError[cb] = scale*(partialCorrEtError[cb])/minEtCorr[cb];
    totEtValues[cb] = scale*(partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/minEtCorr[cb]/fem;
    //cout<<"cb "<<cb<<"\t"<<corrEtValues[cb]<<" = \t"<<scale<<"*\t("<<partialCorrEtValues[cb]<<" -\t"<<hadCorr[cb]<<" -\t"<<kaonCorr[cb]<<" -\t"<<neutronCorr[cb]<<" -\t"<<secondaryCorr[cb]<<")\t/"<<minEtCorr[cb]<<endl;
    //cout<<"cb "<<cb<<"\t"<<corrEtValues[cb]<<" = \t"<<scale<<"*\t("<<partialCorrEtValues[cb]<<" -\t"<<hadCorr[cb]/partialCorrEtValues[cb]<<" -\t"<<kaonCorr[cb]/partialCorrEtValues[cb]<<" -\t"<<neutronCorr[cb]/partialCorrEtValues[cb]<<" -\t"<<secondaryCorr[cb]/partialCorrEtValues[cb]<<")\t/"<<minEtCorr[cb]<<endl;
    corrEtValuesFormulaB[cb] = scale*(rawEtNoEffCorrValues[cb] - hadCorrNoEffCorr[cb] - kaonCorrNoEffCorr[cb] - neutronCorrNoEffCorr[cb] - secondaryCorrNoEffCorr[cb])/minEtCorr[cb]/averageEfficiency[cb];
    //cout<<"cb "<<corrEtValuesFormulaB[cb]<<" = "<<scale<<"*("<<rawEtNoEffCorrValues[cb]<<" - "<<hadCorrNoEffCorr[cb]<<" - "<<kaonCorrNoEffCorr[cb]<<" - "<<neutronCorrNoEffCorr[cb]<<" - "<<secondaryCorrNoEffCorr[cb]<<")/"<<minEtCorr[cb]<<"/"<<averageEfficiency[cb]<<endl;
    //cout<<" fractions: had "<< hadCorrNoEffCorr[cb]/rawEtNoEffCorrValues[cb] <<" kaon "<<kaonCorrNoEffCorr[cb]/rawEtNoEffCorrValues[cb]<<" neutron "<<neutronCorrNoEffCorr[cb]/rawEtNoEffCorrValues[cb]<<" secondary "<<secondaryCorrNoEffCorr[cb]/rawEtNoEffCorrValues[cb]<<endl;
    //signalFraction[cb] = (1.0 - (hadCorrNoEffCorr[cb] + kaonCorrNoEffCorr[cb] + neutronCorrNoEffCorr[cb] + secondaryCorrNoEffCorr[cb])/rawEtNoEffCorrValues[cb]);
    //signalFraction[cb] = (1.0 - (hadCorr[cb] + kaonCorr[cb] + neutronCorr[cb] + secondaryCorr[cb])/partialCorrEtValues[cb]);
    //cout<<"cb "<<cb<<" fractions: \thad "<< hadCorr[cb]/partialCorrEtValues[cb]<<"+/-" << hadError[cb]/partialCorrEtValues[cb]<<"\tkaon "<<kaonCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<kaonError[cb]/partialCorrEtValues[cb]<<"\tneutron "<<neutronCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<neutronError[cb]/partialCorrEtValues[cb]<<"\tsecondary "<<secondaryCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<secondaryError[cb]/partialCorrEtValues[cb]<<endl;//rawEtValues
    //corrEtValuesFormulaC[cb] = scale*partialCorrEtValues[cb]*signalFraction[cb]/minEtCorr[cb];
//     signalFraction[cb] = 0.215;
//     corrEtValuesFormulaC[cb] = scale*rawEtAllValues[cb]*signalFraction[cb]/minEtCorr[cb];
//     //signalFraction[cb] = 0.32;
//     corrEtValuesFormulaC[cb] = scale*rawEtAllNoEffCorrValues[cb]*signalFraction[cb]/minEtCorr[cb];
    //signalFraction[cb] = (partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/rawEtAllNoEffCorrValues[cb];
    signalFraction[cb]= (partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/rawEtAllValues[cb];
    //signalFraction[cb]= (partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/rawEtAllValues[cb];
    //cout<<" f "<<
    //corrEtValuesFormulaC[cb] = scale*rawEtNoEffCorrValues[cb]*signalFraction[cb]/minEtCorr[cb]/0.22;
    //cout<<"cb "<<cb<<" "<<corrEtValuesFormulaC[cb]<<" = "<<scale<<"*"<<rawEtNoEffCorrValues[cb]<<"*"<<signalFraction[cb]<<"/"<<minEtCorr[cb]<<"/"<<averageEfficiency[cb]<<endl;

    totalCorrectionPerNPartPairValues[cb] = hadCorr[cb] + kaonCorr[cb] + neutronCorr[cb] + secondaryCorr[cb];
    totalCorrectionPerNPartPairValuesNoEffCorr[cb] = hadCorrNoEffCorr[cb] + kaonCorrNoEffCorr[cb] + neutronCorrNoEffCorr[cb] + secondaryCorrNoEffCorr[cb];
    //this comes up enough in the error calculations we'll just make a variable for it
    float  partialEt = scale*(partialCorrEtValues[cb])/minEtCorr[cb];
    float  partialEtNoEffCorr = scale*(rawEtNoEffCorrValues[cb])/minEtCorr[cb];
    //the lower the fraction EM ET is of the overall signal the higher the corrections are due to kaons, charged hadrons, neutrons, and secondaries
    float nominalTotEt = scale/minEtCorr[cb]*(partialCorrEtValues[cb]-hadCorr[cb]-kaonCorr[cb]-neutronCorr[cb]-secondaryCorr[cb])/(fem);
    float partialTotEtLow = (partialCorrEtValues[cb]-hadCorr[cb]-hadError[cb]-kaonCorr[cb]-kaonError[cb]-neutronCorr[cb]-neutronError[cb]-secondaryCorr[cb]-secondaryError[cb])/(fem-femerr);
    //and the higher it is the lower the corrections will be
    float partialTotEtHigh = (partialCorrEtValues[cb]-hadCorr[cb]+hadError[cb]-kaonCorr[cb]+kaonError[cb]-neutronCorr[cb]+neutronError[cb]-secondaryCorr[cb]+secondaryError[cb])/(fem+femerr);
    float partialTotEtAvg = (partialTotEtHigh+partialTotEtLow)/2.0;
    float partialTotEtErr = TMath::Abs(partialTotEtHigh-partialTotEtLow)/2.0;
    //cout<<"cb old "<<totEtValues[cb];
    totEtValues[cb] = scale/minEtCorr[cb]*partialTotEtAvg;
    totEtError[cb] = totEtValues[cb] * TMath::Sqrt(TMath::Power(minEtError[cb]/minEtCorr[cb],2)+TMath::Power(partialTotEtErr/partialTotEtAvg,2)+TMath::Power(nonLinError[cb],2)+TMath::Power(efficiencyError[cb],2));
    //cout<<" new "<<totEtValues[cb]<<endl;
    //add up the error squared
    float err = 0;
    float partialerr = 0;
    float totalcorrpartialerr = 0;
    float errNoEffCorr = 0;
    float partialerrNoEffCorr = 0;
    float totalcorrpartialerrNoEffCorr = 0;
    float totaletvariance = 0;
    float femkaoncorrelation = -1.0;//fractional correlation, says these are 100% anticorrelated
    float femhadroncorrelation = -0.0;//fractional correlation, says these are 100% anticorrelated
    float femneutroncorrelation = -1;//fractional correlation, says these are 100% anticorrelated
    float femsecondarycorrelation = 0;//fractional correlation, says these are 100% anticorrelated

    bool writeerror = false;
    bool writetoteterror = false;
    if(writeerror){
      cout<<"fraction errors: min et "<<minEtError[cb]/minEtCorr[cb];
    }
    //if(writeerror)cout<<"partialEt "<<partialEt<<" err^2 = ";
    //Et min correction
    partialerr = TMath::Power(minEtError[cb]/minEtCorr[cb]*corrEtValues[cb],2.0);
    partialerrNoEffCorr = TMath::Power(minEtError[cb]/minEtCorr[cb]*corrEtValuesFormulaB[cb],2.0);
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    partialerr = TMath::Power(minEtError[cb]/minEtCorr[cb]*totEtValues[cb],2.0);
    totaletvariance += partialerr;//the error from min et is not correlated with femet
    if(writetoteterror){cout<<"min et "<<partialerr;}
    //nonlinearity correction - this is saved as a fractional error
    partialerr = TMath::Power(nonLinError[cb]*corrEtValues[cb],2.0);
    partialerrNoEffCorr = TMath::Power(nonLinError[cb]*corrEtValuesFormulaB[cb],2.0);
    if(writeerror){cout<<" nonlinearity "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    partialerr = TMath::Power(nonLinError[cb]*totEtValues[cb],2.0);
    totaletvariance += partialerr;//the error from nonlinearity is not correlated with femet
    if(writetoteterror){cout<<" nonlinearity "<<partialerr;}
    //energy scale - this is saved as a fractional error
    partialerr = TMath::Power(energyscaleerror*corrEtValues[cb],2.0);
    partialerrNoEffCorr = TMath::Power(energyscaleerror*corrEtValuesFormulaB[cb],2.0);
    if(writeerror){cout<<" scale error "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    partialerr = TMath::Power(energyscaleerror*totEtValues[cb],2.0);
    totaletvariance += partialerr;//the error from nonlinearity is not correlated with femet
    if(writetoteterror){cout<<" scale error "<<partialerr;}
    //efficiency correction - this is also saved as a fractional error
    partialerr = TMath::Power(efficiencyError[cb]*corrEtValues[cb],2);
    partialerrNoEffCorr = TMath::Power(efficiencyError[cb]*corrEtValues[cb],2);
    if(writeerror){cout<<" efficiency "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    partialerr = TMath::Power(efficiencyError[cb]*totEtValues[cb],2);
    totaletvariance += partialerr;//the error from efficiency is not correlated with femet
    if(writetoteterror){cout<<" efficiency "<<partialerr;}
    //hadron correction
    partialerr = TMath::Power(hadError[cb]*scale/minEtCorr[cb],2);
    partialerrNoEffCorr = TMath::Power(hadErrorNoEffCorr[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerr += TMath::Power(hadError[cb],2);
    totalcorrpartialerrNoEffCorr += TMath::Power(hadErrorNoEffCorr[cb],2);
    if(writeerror){cout<<" hadronic corr "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    //correlation term is variance due to had ET + correlation fraction * df/dfem * error(fem) * df/dhadcorr * error(hadcorr)
    partialerr = TMath::Power(hadError[cb]*scale/minEtCorr[cb]/fem,2)+femhadroncorrelation * 2 * TMath::Abs(nominalTotEt*femerr/fem * scale/minEtCorr[cb]/fem*hadError[cb]);
    //cout<<" variance had "<<partialerr<<" = "<<TMath::Sqrt(TMath::Abs(partialerr))<<"^2 = "<<TMath::Power(hadError[cb]*scale/minEtCorr[cb]/fem,2)<<" + "<<femhadroncorrelation * TMath::Abs(nominalTotEt*femerr/fem * scale/minEtCorr[cb]/fem*hadError[cb])<<" first part sqrt "<<hadError[cb]*scale/minEtCorr[cb]/fem<<" second part "<<femhadroncorrelation <<" * "<<nominalTotEt*femerr/fem<<" * "<<scale/minEtCorr[cb]/fem*hadError[cb]<<" ";
    totaletvariance += partialerr;//the error from hadd corr IS correlated with femet
    if(writetoteterror){cout<<" hadronic corr "<<partialerr;}
    //neutron correction
    partialerr = TMath::Power(neutronError[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerr += TMath::Power(neutronError[cb],2);
    partialerrNoEffCorr = TMath::Power(neutronErrorNoEffCorr[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerrNoEffCorr += TMath::Power(neutronErrorNoEffCorr[cb],2);
    if(writeerror){cout<<" neutron corr "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    //if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    //correlation term is variance due to had ET + correlation fraction * df/dfem * error(fem) * df/dneutroncorr * error(neutroncorr)
    partialerr = TMath::Power(neutronError[cb]*scale/minEtCorr[cb]/fem,2)+femneutroncorrelation * 2 * TMath::Abs(nominalTotEt*femerr/fem * scale/minEtCorr[cb]/fem*neutronError[cb]);
    totaletvariance += partialerr;//the error from neutron corr IS correlated with femet
    if(writetoteterror){cout<<" neutron corr "<<partialerr;}
    //kaon correction
    partialerr = TMath::Power(kaonError[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerr += TMath::Power(kaonError[cb],2);
    partialerrNoEffCorr = TMath::Power(kaonErrorNoEffCorr[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerrNoEffCorr += TMath::Power(kaonErrorNoEffCorr[cb],2);
    if(writeerror){cout<<" kaon corr "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    //if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    //correlation term is variance due to had ET + correlation fraction * df/dfem * error(fem) * df/dkaoncorr * error(kaoncorr)
    partialerr =TMath::Power(kaonError[cb]*scale/minEtCorr[cb]/fem,2)+femkaoncorrelation * 2 * TMath::Abs(nominalTotEt*femerr/fem * scale/minEtCorr[cb]/fem*kaonError[cb]);
    totaletvariance += partialerr;//the error from kaon corr IS correlated with femet
    if(writetoteterror){cout<<" kaon corr "<<partialerr;}
    //secondary correction
    partialerr = TMath::Power(secondaryError[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerr += TMath::Power(secondaryError[cb],2);
    partialerrNoEffCorr = TMath::Power(secondaryErrorNoEffCorr[cb]*scale/minEtCorr[cb],2);
    totalcorrpartialerrNoEffCorr += TMath::Power(secondaryErrorNoEffCorr[cb],2);
    if(writeerror){cout<<" secondaries "<<TMath::Sqrt(partialerr)/corrEtValues[cb];}
    //if(writeerror)cout<<partialerr;
    err+=partialerr;
    errNoEffCorr+=partialerrNoEffCorr;
    //correlation term is variance due to had ET + correlation fraction * df/dfem * error(fem) * df/dkaoncorr * error(kaoncorr)
    partialerr = TMath::Power(secondaryError[cb]*scale/minEtCorr[cb]/fem,2)+femsecondarycorrelation * 2 * TMath::Abs(nominalTotEt*femerr/fem * scale/minEtCorr[cb]/fem*secondaryError[cb]);
    totaletvariance += partialerr;//the error from kaon corr IS correlated with femet
    if(writetoteterror){cout<<" secondary corr "<<partialerr;}

    partialerr = TMath::Power(nominalTotEt*femerr/fem,2);
    totaletvariance += partialerr;//the error from kaon corr IS correlated with femet
    if(writetoteterror){cout<<" fem "<<partialerr;}
    //And take the square root
    err = TMath::Sqrt(err);
    errNoEffCorr = TMath::Sqrt(errNoEffCorr);
    totalCorrectionPerNPartPairError[cb] = TMath::Sqrt(totalcorrpartialerr);
    totalCorrectionPerNPartPairErrorNoEffCorr[cb] = TMath::Sqrt(totalcorrpartialerrNoEffCorr);
    if(writeerror)cout<<" = "<<err/corrEtValues[cb]<<endl;
    signalFractionError[cb] = totalCorrectionPerNPartPairError[cb]/partialCorrEtValues[cb];
    cout<<"tot et "<<cb<<" "<<totEtValues[cb]<<" +/- "<<totEtError[cb]<<" vs ";
    totEtValues[cb] = nominalTotEt;
    if(totaletvariance>0) totEtError[cb] = TMath::Sqrt(totaletvariance);
    else{cout<<" var err "<<totaletvariance<<" ";}
    //totEtError[cb] = totaletvariance;
    cout<<totEtValues[cb]<<" +/- "<<totEtError[cb]<<endl;
    //cout<<"signal fraction "<<signalFraction[cb]<<" +/- ";
    //cout<<"cb "<<cb<<" fractions: \thad "<< hadCorr[cb]/partialCorrEtValues[cb]<<"+/-" << hadError[cb]/partialCorrEtValues[cb]<<"\tkaon "<<kaonCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<kaonError[cb]/partialCorrEtValues[cb]<<"\tneutron "<<neutronCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<neutronError[cb]/partialCorrEtValues[cb]<<"\tsecondary "<<secondaryCorr[cb]/partialCorrEtValues[cb]<<"+/-"<<secondaryError[cb]/partialCorrEtValues[cb]<<"\tsignalfrac "<<signalFraction[cb]<<"+/-"<<signalFractionError[cb]<<endl;//rawEtValues
    if(partialCorrEtValues[cb]>0)cout<<totalCorrectionPerNPartPairError[cb]/partialCorrEtValues[cb];
    //cout<<endl;
    corrEtError[cb] = err;
    corrEtErrorFormulaB[cb] = errNoEffCorr;
    corrEtPerNPartPairValues[cb] = corrEtValues[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairError[cb]  =  corrEtError[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairValuesFormulaC[cb] = corrEtValuesFormulaC[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairErrorFormulaC[cb]  =  corrEtErrorFormulaC[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairValuesFormulaB[cb] = corrEtValuesFormulaB[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairErrorFormulaB[cb]  =  corrEtErrorFormulaB[cb]/(npart[cb]/2.0);
    totalCorrectionPerNPartPairValues[cb] = totalCorrectionPerNPartPairValues[cb]/20.0/trackmultiplicity[cb];
    totalCorrectionPerNPartPairError[cb] = totalCorrectionPerNPartPairError[cb]/20.0/trackmultiplicity[cb];
    totalCorrectionPerNPartPairValuesNoEffCorr[cb] = totalCorrectionPerNPartPairValuesNoEffCorr[cb]/20.0/trackmultiplicity[cb];
    totalCorrectionPerNPartPairErrorNoEffCorr[cb] = totalCorrectionPerNPartPairErrorNoEffCorr[cb]/20.0/trackmultiplicity[cb];
    //cout<<"test cb "<<cb<<" total corr/npart pair "<<totalCorrectionPerNPartPairValues[cb]<<" "<<totalCorrectionPerNPartPairError[cb]<<endl;
    //cout<<"cb "<<cb <<" et "<< corrEtPerNPartPairValues[cb] <<" +/- "<<corrEtPerNPartPairError[cb];
    //cout<<"cb "<<cb <<" et "<< corrEtPerNPartPairValuesFormulaC[cb] <<" +/- "<<corrEtPerNPartPairErrorFormulaC[cb]<<endl;

    //cout<<" = "<<scale<<"*"<<"("<<partialCorrEtValues[cb]<<"+/-"<<partialCorrEtError[cb]<<" - "<<hadCorr[cb]<<"+/-"<<hadError[cb]<<" - "<<kaonCorr[cb]<<"+/-"<<kaonError[cb]<<" - "<<neutronCorr[cb]<<"+/-"<<neutronError[cb]<<" - "<<secondaryCorr[cb]<<"+/-"<<secondaryError[cb]<<")/"<<minEtCorr[cb]<<"+/-"<<minEtError[cb];
    //cout<<endl;
  }
//   for(int cb = 0;cb<10;cb++){
//   }

}

TGraphErrors *GetEtGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,corrEtPerNPartPairValues,npartErr,corrEtPerNPartPairError);
    SetStyles(gr3,25,1);
    return gr3;
}

TGraphErrors *GetEtGraphFormulaC(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,corrEtPerNPartPairValuesFormulaC,npartErr,corrEtPerNPartPairErrorFormulaC);
    for(int i=0;i<20;i++){
      //cout<<"i "<<i<<" "<<clustermultiplicity[i]<<": "<<rawEtNoEffCorrPerNClValues[i]<<"+/-"<<rawEtNoEffCorrPerNClError[i]<<endl;
    }
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetEtGraphFormulaB(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,corrEtPerNPartPairValuesFormulaB,npartErr,corrEtPerNPartPairErrorFormulaB);
    for(int i=0;i<20;i++){
      //cout<<"i "<<i<<" "<<clustermultiplicity[i]<<": "<<rawEtNoEffCorrPerNClValues[i]<<"+/-"<<rawEtNoEffCorrPerNClError[i]<<endl;
    }
    SetStyles(gr3,29,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,partialCorrEtPerNChValues,trackmultiplicityError,partialCorrEtPerNChError);
    //TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);

    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNPartGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,partialCorrEtPerNPartValues,npartErr,partialCorrEtPerNPartError);
    //TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);
    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,partialCorrEtPerNClValues,clustermultiplicityError,partialCorrEtPerNClError);
    //TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);

    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetRawEtNoEffCorrPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtNoEffCorrPerNChValues,trackmultiplicityError,rawEtNoEffCorrPerNChError);
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetRawEtNoEffCorrPerNPartGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,rawEtNoEffCorrPerNPartValues,npartErr,rawEtNoEffCorrPerNPartError);
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetRawEtNoEffCorrPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtNoEffCorrPerNClValues,clustermultiplicityError,rawEtNoEffCorrPerNClError);
    for(int i=0;i<20;i++){
      //cout<<"i "<<i<<" "<<clustermultiplicity[i]<<": "<<rawEtNoEffCorrPerNClValues[i]<<"+/-"<<rawEtNoEffCorrPerNClError[i]<<endl;
    }
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetRawEtAllNoEffCorrPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtAllNoEffCorrPerNChValues,trackmultiplicityError,rawEtAllNoEffCorrPerNChError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllNoEffCorrPerNPartGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,rawEtAllNoEffCorrPerNPartValues,npartErr,rawEtAllNoEffCorrPerNPartError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllNoEffCorrPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtAllNoEffCorrPerNClValues,clustermultiplicityError,rawEtAllNoEffCorrPerNClError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtAllPerNChValues,trackmultiplicityError,rawEtAllPerNChError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllPerNPartGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,rawEtAllPerNPartValues,npartErr,rawEtAllPerNPartError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtAllPerNClValues,clustermultiplicityError,rawEtAllPerNClError);
    SetStyles(gr3,22,1);
    return gr3;
}
TGraphErrors *GetMatchedTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,matchedtrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,20,1);
    return gr3;
}
TGraphErrors *GetMatchedTracksPerNPartGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,matchedtrackmultiplicityPerNPart,npartErr,arrayofzeros);
    SetStyles(gr3,20,1);
    return gr3;
}
TGraphErrors *GetMatchedTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,matchedtrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,20,1);
    return gr3;
}
TGraphErrors *GetTotalTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,totaltrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,34,1);
    return gr3;
}
TGraphErrors *GetTotalTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,totaltrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,34,1);
    return gr3;
}
TGraphErrors *GetNotMatchedTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,notmatchedtrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,24,1);
    return gr3;
}
TGraphErrors *GetNotMatchedTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,notmatchedtrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,24,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNPartPairGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,partialCorrEtPerNPartPairValues,npartErr,partialCorrEtPerNPartPairError);
    SetStyles(gr3,22,1);
    return gr3;
}
TGraphErrors *GetTotalCorrectionPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,totalCorrectionPerNPartPairValues,trackmultiplicityError,totalCorrectionPerNPartPairError);
    SetStyles(gr3,34,TColor::kOrange);
    return gr3;
}
TGraphErrors *GetTotalCorrectionPerNChGraphNoEffCorr(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,totalCorrectionPerNPartPairValuesNoEffCorr,trackmultiplicityError,totalCorrectionPerNPartPairErrorNoEffCorr);
    SetStyles(gr3,34,TColor::kOrange);
    return gr3;
}
//=================================Plotting code====================================
Bool_t sim = false;
//Bool_t isPhos = kFALSE;
TString detector = "";
TString year = "2010";
//void PlotEmEtVer2(TString filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root")
void PlotEmEtVer2(Bool_t isPhos = kFALSE, Bool_t isMC = kFALSE, Int_t cutset = 0)
{

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  CutSet = cutset;
  TString filename, simfilename;
  if(cutset==0){
    if(isPhos){
      if(isMC){
	//filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.root";
      }
      else{
	filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOS.LHC10hPass2.Run139465.root";
      }
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
    }
    else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
      if(isMC){
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
      }
      else{
	filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root";
      }
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
    }
  }
  if(cutset==1){
    if(isPhos){
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOS350MeVCut.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS350MeVCut.LHC11a10a_bis.Run139465.root";
    }
    else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal350MeVCut.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal350MeVCut.LHC11a10a_bis.Run139465.root";
    }
  }
  if(cutset==2){
    if(isPhos){
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOS400MeVCut.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS400MeVCut.LHC11a10a_bis.Run139465.root";
    }
    else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal400MeVCut.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal400MeVCut.LHC11a10a_bis.Run139465.root";
    }
  }
  if(cutset==3){
    if(isPhos){
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOSOldHadMethod.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOSOldHadMethod.LHC11a10a_bis.root";
      corrfac = 0.174236;
      corrfacerr = 0.0423712;
    }
    else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCalOldHadMethod.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCalOldHadMethod.LHC11a10a_bis.root";
      corrfac = 0.283921;
      corrfacerr = 0.0362492;
    }
  }
  if(cutset==4){
    year = "2011";
    if(isPhos){
      if(isMC){
	filename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.PHOS.LHC13e1abc.root";
      }
      else{
	filename = "rootFiles/LHC11hPass2/Et.ESD.realPbPb.PHOS.LHC11hPass2.Run168464.root";
      }
      corrfac = 0.16382;
      corrfacerr = 0.0830542;
      simfilename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.PHOS.LHC13e1abc.root";
    }
    else{
      if(isMC){
	filename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.EMCal.LHC13e1abc.root";
      }
      else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
	filename = "rootFiles/LHC11hPass2/Et.ESD.realPbPb.EMCal.LHC11hPass2.Run168464.root";
      }
      corrfac = 0.196084;
      corrfacerr = 0.0410244;
      simfilename = "rootFiles/LHC13e1abcCombined/Et.ESD.simPbPb.EMCal.LHC13e1abc.root";
    }
  }
  cout<<"data file name = "<<filename<<endl;
  cout<<" sim file name = "<<simfilename<<endl;
  TString detector = "Emcal";
  Float_t scale = 360.0/40.0/1.2*1.02;//Azimuthal acceptance over eta range, times energy scale correction
  if(filename.Contains("PHOS")){
    detector = "Phos";
    //isPhos = kTRUE;
    scale = 360.0/60/0.24;//Azimuthal acceptance over eta range
    energyscaleerror = 0.005;
  }
  else{
    if(cutset==4){//2011 data
      cout<<"scale "<<scale<<endl;
      scale = 360.0/100.0/1.2*1.02;//Azimuthal acceptance over eta range, times energy scale correction
      cout<<"scale "<<scale<<endl;
    }
  }
  TString det = detector;
  //gROOT->LoadMacro("macros/PlotSecondariesCorr.C");
  //PlotSecondariesCorr(simfilename,filename);
  ReadInData(filename,detector);
  ReadInNeutronCorrections();
  ReadInSecondaryCorrections();
  ReadInKaonCorrections();
  ReadMinEtCorrections();
  CalculateHadronCorrections(isPhos);
  ApplyCorrections(scale);

  TCanvas *c1 = new TCanvas("c1","Corrections",500,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.120968);
  //c1->SetRightMargin();
  //c1->SetTopMargin();
  c1->SetBottomMargin(0.15);
  c1->SetLogx();
  TGraphErrors *graphKaonCorrectionShort = GetKaonCorrectionGraphShort();
  graphKaonCorrectionShort->Draw("AP");
  //graphKaonCorrectionShort->SetMaximim(0.02);
  //graphKaonCorrectionShort->GetYaxis()->SetRange(1,graphKaonCorrectionShort->GetYaxis()->FindBin(0.02));
  float emcalmax = 0.006;
  if(cutset==4){//2011 data
    emcalmax = 7.65/scale*emcalmax;
  }
  graphKaonCorrectionShort->GetHistogram()->SetMaximum(emcalmax);
  //set scales the same within naive geometric scaling
  if(isPhos) graphKaonCorrectionShort->GetHistogram()->SetMaximum(emcalmax*0.24/1.2*60.0/40.0);
  graphKaonCorrectionShort->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  graphKaonCorrectionShort->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  graphKaonCorrectionShort->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  graphKaonCorrectionShort->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  graphKaonCorrectionShort->GetXaxis()->SetTitle("N_{Ch}");
  graphKaonCorrectionShort->GetYaxis()->SetTitle("Correction/N_{Ch}");
  TGraphErrors *graphKaonCorrection = GetKaonCorrectionGraph();
  graphKaonCorrection->Draw("P same");
  TGraphErrors *graphSecondaryCorrectionShort = GetSecondaryCorrectionGraphShort();
  graphSecondaryCorrectionShort->Draw("P same");
  TGraphErrors *graphNeutronCorrectionShort = GetNeutronCorrectionGraphShort();
  graphNeutronCorrectionShort->Draw("P same");
  TGraphErrors *graphSecondaryCorrection = GetSecondaryCorrectionGraph();
  graphSecondaryCorrection->Draw("P same");
  TGraphErrors *graphNeutronCorrection = GetNeutronCorrectionGraph();
  graphNeutronCorrection->Draw("P same");
  TGraphErrors *graphHadronCorrection = GetHadronCorrectionGraph();
  graphHadronCorrection->Draw("P same");
  TGraphErrors *graphHadronCorrectionShort = GetHadronCorrectionGraphShort();
  graphHadronCorrectionShort->Draw("P same");
  TGraphErrors *graphTotalCorrection = GetTotalCorrectionPerNChGraph();
  graphTotalCorrection->Draw("P same");
  TLegend *leg = new TLegend(0.538306,0.756684,0.739919,0.97861);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.038682);
  leg->AddEntry(graphKaonCorrectionShort,"Kaon correction/N_{ch}","P");
  leg->AddEntry(graphNeutronCorrectionShort,"Neutron correction/2.0 N_{ch}","P");
  leg->AddEntry(graphSecondaryCorrectionShort,"Secondary correction/5.0 N_{ch}","P");
  leg->AddEntry(graphHadronCorrectionShort,"Hadron correction/5.0 N_{ch}","P");
  leg->AddEntry(graphTotalCorrection,"(E_{T}^{kaon}+E_{T}^{n}+E_{T}^{secondary}+E_{T}^{h})/20.0 N_{ch}","P");
  leg->Draw();
  TString corrplotname1 = "/tmp/CorrectionsVsNch"+detector+".eps";
  c1->SaveAs(corrplotname1.Data());
//   TH1F *frame = new TH1F("frame","frame",1,0,2);
//   frame->GetYaxis()->SetTitle("dE_{T}/d#eta");
//   frame->GetXaxis()->SetTitle("N_{part}");
//   //fPion->SetRange(0,2);
//   frame->SetMinimum(0);
//   frame->SetMaximum(10);
//   TCanvas *c1a = new TCanvas("c1a","CorrectionsNoEffCorr",500,400);
//   c1a->SetTopMargin(0.02);
//   c1a->SetRightMargin(0.02);
//   c1a->SetBorderSize(0);
//   c1a->SetFillColor(0);
//   c1a->SetFillColor(0);
//   c1a->SetBorderMode(0);
//   c1a->SetFrameFillColor(0);
//   c1a->SetFrameBorderMode(0);
//   c1a->SetLeftMargin(0.120968);
//   //c1a->SetRightMargin();
//   //c1a->SetTopMargin();
//   c1a->SetBottomMargin(0.15);
//   c1a->SetLogx();
//   TGraphErrors *graphKaonCorrectionShortNoEffCorr = GetKaonCorrectionGraphShortNoEffCorr();
//   graphKaonCorrectionShortNoEffCorr->Draw("AP");
//   //graphKaonCorrectionShort->SetMaximim(0.02);
//   //graphKaonCorrectionShort->GetYaxis()->SetRange(1,graphKaonCorrectionShort->GetYaxis()->FindBin(0.02));
//   float emcalmax = 0.006;
//   graphKaonCorrectionShortNoEffCorr->GetHistogram()->SetMaximum(emcalmax);
//   //set scales the same within naive geometric scaling
//   if(isPhos) graphKaonCorrectionShortNoEffCorr->GetHistogram()->SetMaximum(emcalmax*0.24/1.2*60.0/40.0);
//   graphKaonCorrectionShortNoEffCorr->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
//   graphKaonCorrectionShortNoEffCorr->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
//   graphKaonCorrectionShortNoEffCorr->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
//   graphKaonCorrectionShortNoEffCorr->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
//   graphKaonCorrectionShortNoEffCorr->GetXaxis()->SetTitle("N_{Ch}");
//   graphKaonCorrectionShortNoEffCorr->GetYaxis()->SetTitle("Correction/N_{Ch}");
//   TGraphErrors *graphKaonCorrectionNoEffCorr = GetKaonCorrectionGraphNoEffCorr();
//   graphKaonCorrectionNoEffCorr->Draw("P same");
//   TGraphErrors *graphSecondaryCorrectionShortNoEffCorr = GetSecondaryCorrectionGraphShortNoEffCorr();
//   graphSecondaryCorrectionShortNoEffCorr->Draw("P same");
//   TGraphErrors *graphNeutronCorrectionShortNoEffCorr = GetNeutronCorrectionGraphShortNoEffCorr();
//   graphNeutronCorrectionShortNoEffCorr->Draw("P same");
//   TGraphErrors *graphSecondaryCorrectionNoEffCorr = GetSecondaryCorrectionGraphNoEffCorr();
//   graphSecondaryCorrectionNoEffCorr->Draw("P same");
//   TGraphErrors *graphNeutronCorrectionNoEffCorr = GetNeutronCorrectionGraphNoEffCorr();
//   graphNeutronCorrectionNoEffCorr->Draw("P same");
//   TGraphErrors *graphHadronCorrectionNoEffCorr = GetHadronCorrectionGraphNoEffCorr();
//   graphHadronCorrectionNoEffCorr->Draw("P same");
//   TGraphErrors *graphHadronCorrectionShortNoEffCorr = GetHadronCorrectionGraphShortNoEffCorr();
//   graphHadronCorrectionShortNoEffCorr->Draw("P same");
//   TGraphErrors *graphTotalCorrectionNoEffCorr = GetTotalCorrectionPerNChGraphNoEffCorr();
//   graphTotalCorrectionNoEffCorr->Draw("P same");
//   TLegend *leg2 = new TLegend(0.538306,0.756684,0.739919,0.97861);
//   leg2->SetFillStyle(0);
//   leg2->SetFillColor(0);
//   leg2->SetBorderSize(0);
//   leg2->SetTextSize(0.03);
//   leg2->SetTextSize(0.038682);
//   leg2->AddEntry(graphKaonCorrectionShortNoEffCorr,"Kaon correction/N_{ch}","P");
//   leg2->AddEntry(graphNeutronCorrectionShortNoEffCorr,"Neutron correction/N_{ch}","P");
//   leg2->AddEntry(graphSecondaryCorrectionShortNoEffCorr,"Secondary correction/5.0 N_{ch}","P");
//   leg2->AddEntry(graphHadronCorrectionShortNoEffCorr,"Hadron correction/5.0 N_{ch}","P");
//   leg2->AddEntry(graphTotalCorrectionNoEffCorr,"(E_{T}^{kaon}+E_{T}^{n}+E_{T}^{secondary}+E_{T}^{h})/20.0 N_{ch}","P");
//   leg2->Draw();
//   TString corrplotname1 = "/tmp/CorrectionsVsNch"+detector+"NoEffCorr.eps";
//   c1a->SaveAs(corrplotname1.Data());
//   TCanvas *c2 = new TCanvas("c2","Min Et Correction",600,400);
//   c2->SetTopMargin(0.02);
//   c2->SetRightMargin(0.02);
//   c2->SetBorderSize(0);
//   c2->SetFillColor(0);
//   c2->SetFillColor(0);
//   c2->SetBorderMode(0);
//   c2->SetFrameFillColor(0);
//   c2->SetFrameBorderMode(0);
//   c2->SetLogx();
//   TGraphErrors *graphMinEtCorrectionShort = GetMinEtCorrectionGraphShort();
//   graphMinEtCorrectionShort->GetXaxis()->SetTitle("f_{ETmin}");
//   graphMinEtCorrectionShort->GetYaxis()->SetTitle("Correction/N_{Ch}");
// //   graphMinEtCorrectionShort->GetHistogram()->SetMaximum(1.0);
// //   graphMinEtCorrectionShort->GetHistogram()->SetMinimum(0.0);
//   graphMinEtCorrectionShort->Draw("AP");
//   graphKaonCorrectionShort->GetXaxis()->SetTitle("N_{Ch}");
//   graphKaonCorrectionShort->GetYaxis()->SetTitle("f_{EtMin}");
//   TGraphErrors *graphMinEtCorrection = GetMinEtCorrectionGraph();
//   graphMinEtCorrection->Draw("P same");
//   TString corrplotname2 = "/tmp/EtMinVsNch"+detector+".eps";
//   c2->SaveAs(corrplotname2.Data());

//   TCanvas *c3 = new TCanvas("c3","Kaon Correction",600,400);
//   c3->SetTopMargin(0.02);
//   c3->SetRightMargin(0.02);
//   c3->SetBorderSize(0);
//   c3->SetFillColor(0);
//   c3->SetFillColor(0);
//   c3->SetBorderMode(0);
//   c3->SetFrameFillColor(0);
//   c3->SetFrameBorderMode(0);
//   c3->SetLogx();

//   TGraphErrors *graphKaonGraph = GetKaonGraph();
//   graphKaonGraph->Draw("AP");
//   //graphKaonGraph->GetHistogram()->SetMaximum();
//   graphKaonGraph->GetHistogram()->SetMinimum(0.0);
//   graphKaonGraph->GetXaxis()->SetTitle("N_{Ch}");
//   graphKaonGraph->GetYaxis()->SetTitle("val/N_{Ch}");
//   graphKaonCorrectionShort->Draw("P same");
//   TGraphErrors *graphKaonEtGraph = GetKaonEtGraph();
//   graphKaonEtGraph->Draw("P same");
//   TLegend *leg3 = new TLegend(0.602349,0.362903,0.805369,0.586022);
//   leg3->SetFillStyle(0);
//   leg3->SetFillColor(0);
//   leg3->SetBorderSize(0);
//   leg3->SetTextSize(0.03);
//   leg3->SetTextSize(0.038682);
//   leg3->AddEntry(graphKaonGraph,"N_{K^{#pm}}/N_{Ch}","P");
//   leg3->AddEntry(graphKaonCorrectionShort,"E_{T}^{kaon}","P");
//   leg3->AddEntry(graphKaonEtGraph,"E_{T}^{K^{#pm}}/N_{Ch}","P");
//   leg3->Draw();
//   TString corrplotname3 = "/tmp/KaonVsNch"+detector+".eps";
//   c3->SaveAs(corrplotname3.Data());


  TCanvas *c4 = new TCanvas("c4", "dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]",700, 600);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.02);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  TGraphErrors *graphEt = GetEtGraph();
  graphEt->GetHistogram()->GetXaxis()->SetTitle("N_{part}");
  graphEt->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/(N_{part}/2)");
  graphEt->GetHistogram()->SetMinimum(0.0);
  graphEt->GetHistogram()->SetMaximum(3.0);
  graphEt->Draw("AP");
  TGraphErrors *graphPartialCorr = GetPartialCorrEtPerNPartPairGraph();
  //graphPartialCorr->Draw("P same");
  TGraphErrors *graphPionEt = GetPionEtGraph();
  graphPionEt->Draw("P same");
  TGraphErrors *graphPionEmEt = GetPionEmEtGraph();
  graphPionEmEt->Draw("P same");
  TGraphErrors *graphEtFormulaC = GetEtGraphFormulaC();
  //graphEtFormulaC->Draw("P same");
  TGraphErrors *graphEtFormulaB = GetEtGraphFormulaB();
  //graphEtFormulaB->Draw("P same");
  TLegend *leg2 = new TLegend(0.607383,0.704301,0.810403,0.927419);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetTextSize(0.038682);
  leg2->AddEntry(graphEt,"E_{T}^{EM}","P");
  leg2->AddEntry(graphPionEt,"E_{T}^{#pi^{#pm}}","P");
  leg2->AddEntry(graphPionEmEt,"E_{T}^{EM calc. from #pi}","P");
  //leg2->AddEntry(graphEtFormulaC,"E_{T}^{EM} Formula C","P");
  //leg2->AddEntry(graphEtFormulaB,"E_{T}^{EM} Formula B","P");
  leg2->Draw();
  TString corrplotname4 = "/tmp/EtVsNpart"+detector+".eps";
  c4->SaveAs(corrplotname4.Data());
  TString corrplotname4a = "/tmp/EtVsNpart"+detector+".png";
  c4->SaveAs(corrplotname4a.Data());

//   TCanvas *c5 = new TCanvas("c5", "Raw Et Vs NCh",1200, 400);
//   c5->SetTopMargin(0.02);
//   c5->SetRightMargin(0.02);
//   c5->SetBorderSize(0);
//   c5->SetFillColor(0);
//   c5->SetFillColor(0);
//   c5->SetBorderMode(0);
//   c5->SetFrameFillColor(0);
//   c5->SetFrameBorderMode(0);
//   c5->Divide(2);
//   c5->cd(1);
//   TGraphErrors *partialCorrEtPerNCh = GetPartialCorrEtPerNChGraph();
//   partialCorrEtPerNCh->GetHistogram()->SetMaximum(0.1*7.5/scale);//so PHOS and EMCal scales are similar...  partialCorrEtPerNCh->GetHistogram()->SetMinimum(0.0);
//   partialCorrEtPerNCh->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
//   partialCorrEtPerNCh->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
//   partialCorrEtPerNCh->Draw("AP");
//   TGraphErrors *graphRawNoEffCorrCh = GetRawEtNoEffCorrPerNChGraph();
//   graphRawNoEffCorrCh->Draw("P same");
//   TGraphErrors *graphRawAllNoEffCorrCh = GetRawEtAllNoEffCorrPerNChGraph();
//   graphRawAllNoEffCorrCh->Draw("P same");
//   TGraphErrors *graphRawAllCh = GetRawEtAllPerNChGraph();
//   graphRawAllCh->Draw("P same");
//   TLegend *leg5 = new TLegend(0.416667,0.173077,0.62069,0.396853);
//   leg5->SetFillStyle(0);
//   leg5->SetFillColor(0);
//   leg5->SetBorderSize(0);
//   leg5->SetTextSize(0.03);
//   leg5->SetTextSize(0.038682);
//   leg5->AddEntry(partialCorrEtPerNCh,"E_{T}#delta_{TM}/#epsilon f_{nonlin}","P");
//   leg5->AddEntry(graphRawNoEffCorrCh,"E_{T}#delta_{TM}/ f_{nonlin}","P");
//   leg5->AddEntry(graphRawAllCh,"E_{T}/#epsilon f_{nonlin}","P");
//   leg5->AddEntry(graphRawAllNoEffCorrCh,"E_{T}/f_{nonlin}","P");
//   leg5->Draw();
//   TString corrplotname5 = "/tmp/RawEtVsNch"+detector+".eps";
//   c5->SaveAs(corrplotname5.Data());

// //   TCanvas *c6 = new TCanvas("c6", "Raw Et Vs NCl",700, 600);
// //   c6->SetTopMargin(0.02);
// //   c6->SetRightMargin(0.02);
// //   c6->SetBorderSize(0);
// //   c6->SetFillColor(0);
// //   c6->SetFillColor(0);
// //   c6->SetBorderMode(0);
// //   c6->SetFrameFillColor(0);
// //   c6->SetFrameBorderMode(0);

//   c5->cd(2);
//   TGraphErrors *graphRawAllCl = GetRawEtAllPerNClGraph();
//   graphRawAllCl->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
//   graphRawAllCl->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
//   graphRawAllCl->Draw("AP");
//   TGraphErrors *partialCorrEtPerNCl = GetPartialCorrEtPerNClGraph();
//   partialCorrEtPerNCl->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
//   partialCorrEtPerNCl->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
//   partialCorrEtPerNCl->Draw("P same");
//   TGraphErrors *graphRawNoEffCorrCl = GetRawEtNoEffCorrPerNClGraph();
//   graphRawNoEffCorrCl->Draw("P same");
//   TGraphErrors *graphRawAllNoEffCorrCl = GetRawEtAllNoEffCorrPerNClGraph();
//   graphRawAllNoEffCorrCl->Draw("P same");
//   TLegend *leg6 = new TLegend(0.416667,0.173077,0.62069,0.396853);
//   leg6->SetFillStyle(0);
//   leg6->SetFillColor(0);
//   leg6->SetBorderSize(0);
//   leg6->SetTextSize(0.03);
//   leg6->SetTextSize(0.038682);
//   leg6->AddEntry(partialCorrEtPerNCl,"E_{T}#delta_{TM}/#epsilon f_{nonlin}","P");
//   leg6->AddEntry(graphRawNoEffCorrCl,"E_{T}#delta_{TM}/ f_{nonlin}","P");
//   leg6->AddEntry(graphRawAllCl,"E_{T}/#epsilon f_{nonlin}","P");
//   leg6->AddEntry(graphRawAllNoEffCorrCl,"E_{T}/f_{nonlin}","P");
//   leg6->Draw();
//   TString corrplotname5 = "/tmp/RawEtVsNcl"+detector+".eps";
//   c5->SaveAs(corrplotname5.Data());

//   TCanvas *c7 = new TCanvas("c7", "Matched tracks",1200, 400);
//   c7->SetTopMargin(0.02);
//   c7->SetRightMargin(0.02);
//   c7->SetBorderSize(0);
//   c7->SetFillColor(0);
//   c7->SetFillColor(0);
//   c7->SetBorderMode(0);
//   c7->SetFrameFillColor(0);
//   c7->SetFrameBorderMode(0);
//   c7->Divide(2);
//   c7->cd(1);
//   TGraphErrors *graphTotalPerNCh = GetTotalTracksPerNChGraph();
//   graphTotalPerNCh->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
//   graphTotalPerNCh->GetHistogram()->GetYaxis()->SetTitle("Number of tracks/N_{Ch}");
//   graphTotalPerNCh->GetHistogram()->SetMaximum(0.06*7.5/scale);//so PHOS and EMCal scales are similar...
//   graphTotalPerNCh->Draw("AP");
//   TGraphErrors *graphMatchedPerNCh = GetMatchedTracksPerNChGraph();
//   graphMatchedPerNCh->Draw("P same");
//   TGraphErrors *graphNotMatchedPerNCh = GetNotMatchedTracksPerNChGraph();
//   graphNotMatchedPerNCh->Draw("P same");
//   TLegend *leg7 = new TLegend(0.176003,0.637152,0.379808,0.86208);
//   leg7->SetFillStyle(0);
//   leg7->SetFillColor(0);
//   leg7->SetBorderSize(0);
//   leg7->SetTextSize(0.03);
//   leg7->SetTextSize(0.038682);
//   leg7->AddEntry(graphMatchedPerNCh,"Number of tracks matched to clusters","P");
//   leg7->AddEntry(graphNotMatchedPerNCh,"Number of tracks not matched to clusters","P");
//   leg7->AddEntry(graphTotalPerNCh,"Total number of tracks","P");
//   leg7->Draw();
//   c7->cd(2);
//   TGraphErrors *graphTotalPerNCl = GetTotalTracksPerNClGraph();
//   graphTotalPerNCl->GetHistogram()->GetXaxis()->SetTitle("N_{Cl}");
//   graphTotalPerNCl->GetHistogram()->GetYaxis()->SetTitle("Number of tracks/N_{Cl}");
//   graphTotalPerNCl->GetHistogram()->SetMaximum(0.65);
//   graphTotalPerNCl->Draw("AP");
//   TGraphErrors *graphMatchedPerNCl = GetMatchedTracksPerNClGraph();
//   graphMatchedPerNCl->Draw("P same");
//   TGraphErrors *graphNotMatchedPerNCl = GetNotMatchedTracksPerNClGraph();
//   graphNotMatchedPerNCl->Draw("P same");
//   TLine *line = new TLine(graphTotalPerNCl->GetHistogram()->GetXaxis()->GetBinLowEdge(1),0.5,graphTotalPerNCl->GetHistogram()->GetXaxis()->GetBinLowEdge(graphTotalPerNCl->GetHistogram()->GetXaxis()->GetNbins()),0.5);
//   line->Draw();
//   TString corrplotname7 = "/tmp/TrackMatchCrossCheck"+detector+".eps";
//   c7->SaveAs(corrplotname7.Data());

  TCanvas *c8 = new TCanvas("c8","Signal Fraction",600,400);
  c8->SetTopMargin(0.02);
  c8->SetRightMargin(0.02);
  c8->SetBorderSize(0);
  c8->SetFillColor(0);
  c8->SetFillColor(0);
  c8->SetBorderMode(0);
  c8->SetFrameFillColor(0);
  c8->SetFrameBorderMode(0);
  //c8->SetLogx();
  TGraphErrors *graphSignalFraction = GetSignalFractionGraph();
  graphSignalFraction->GetXaxis()->SetTitle("f_{ETmin}");
  graphSignalFraction->GetYaxis()->SetTitle("Correction/N_{Ch}");
  graphSignalFraction->GetHistogram()->SetMaximum(0.4);//so PHOS and EMCal scales are similar...
  graphSignalFraction->GetHistogram()->SetMinimum(0.0);//so PHOS and EMCal scales are similar...
  //graphSignalFraction->GetHistogram()->SetMaximum(1.0);
  //graphSignalFraction->GetHistogram()->SetMinimum(0.0);
  graphSignalFraction->Draw("AP");
  TString corrplotname8 = "/tmp/SignalFraction"+detector+".eps";
  c8->SaveAs(corrplotname8.Data());

  TCanvas *c9 = new TCanvas("c9","Event by event ET distribution",600,400);
  c9->SetTopMargin(0.02);
  c9->SetRightMargin(0.02);
  c9->SetBorderSize(0);
  c9->SetFillColor(0);
  c9->SetFillColor(0);
  c9->SetBorderMode(0);
  c9->SetFrameFillColor(0);
  c9->SetFrameBorderMode(0);
  c9->SetLogy();
  c9->SetLogx();
  int nbins = 15;
  if(isPhos) nbins = 8;
  //nbins = 20;
  for(bin = 1; bin<=nbins;bin++){
    dataFits[bin] = new TF1(Form("Fit%i",bin),"[0]*exp(-(x-[1])*(x-[1])/[2])",0,250);
    ((TF1*)dataFits[bin])->SetParameter(0,1e3);
    ((TF1*)dataFits[bin])->SetParameter(1,partialCorrEtValues[bin-1]);
    ((TF1*)dataFits[bin])->SetParameter(2,20);
    //((TH1D*)rawEt[bin])->Fit(((TF1*)dataFits[bin]));
  }
  int bin = 1;
  ((TH1D*)rawEt[bin])->GetXaxis()->SetTitle("E_{T}^{raw}");
  //if(isPhos){((TH1D*)rawEt[bin])->GetXaxis()->SetRange(1,((TH1D*)rawEt[bin])->GetXaxis()->FindBin(100));}
  //else{((TH1D*)rawEt[bin])->GetXaxis()->SetRange(1,((TH1D*)rawEt[bin])->GetXaxis()->FindBin(200));}
  ((TH1D*)rawEt[bin])->Draw();
  ((TH1D*)rawEt[bin])->SetMaximum(5e4);
  for(bin = 1; bin<=nbins;bin+=1){
    ((TH1D*)rawEt[bin])->Draw("same");
    ((TH1D*)rawEt[bin])->SetLineColor(colors[bin-1]);
    ((TH1D*)rawEt[bin])->SetMarkerColor(colors[bin-1]);
    //((TF1*)dataFits[bin])->SetLineColor(colors[bin-1]);
    //Float_t mostlikely = GetMostLikelyValue((TH1D*)rawEt[bin]);
    //cout<<"cb "<<bin-1<<" average: "<<partialCorrEtValues[bin-1]<<" most likely "<<mostlikely<<" fit "<<((TF1*)dataFits[bin])->GetParameter(1)<<endl;
  }
  TString name9 = "ananotepics/emEt/ETDistribution"+detector+".eps";
  c9->SaveAs(name9.Data());


//   TCanvas *c10 = new TCanvas("c10", "Raw Et Vs NPart",1200, 400);
//   c10->SetTopMargin(0.02);
//   c10->SetRightMargin(0.02);
//   c10->SetBorderSize(0);
//   c10->SetFillColor(0);
//   c10->SetFillColor(0);
//   c10->SetBorderMode(0);
//   c10->SetFrameFillColor(0);
//   c10->SetFrameBorderMode(0);
//   c10->Divide(2);
//   c10->cd(1);
//   TGraphErrors *partialCorrEtPerNPart = GetPartialCorrEtPerNPartGraph();
//   partialCorrEtPerNPart->GetHistogram()->SetMaximum(0.1*7.5/scale);//so PHOS and EMCal scales are similar...  partialCorrEtPerNPart->GetHistogram()->SetMinimum(0.0);
//   partialCorrEtPerNPart->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
//   partialCorrEtPerNPart->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
//   partialCorrEtPerNPart->Draw("AP");
//   TGraphErrors *graphRawNoEffCorrNPart = GetRawEtNoEffCorrPerNPartGraph();
//   graphRawNoEffCorrNPart->Draw("P same");
//   TGraphErrors *graphRawAllNoEffCorrNPart = GetRawEtAllNoEffCorrPerNPartGraph();
//   graphRawAllNoEffCorrNPart->Draw("P same");
//   TGraphErrors *graphRawAllNPart = GetRawEtAllPerNPartGraph();
//   graphRawAllNPart->Draw("P same");
//   TLegend *leg5 = new TLegend(0.416667,0.173077,0.62069,0.396853);
//   leg5->SetFillStyle(0);
//   leg5->SetFillColor(0);
//   leg5->SetBorderSize(0);
//   leg5->SetTextSize(0.03);
//   leg5->SetTextSize(0.038682);
//   leg5->AddEntry(partialCorrEtPerNPart,"E_{T}#delta_{TM}/#epsilon f_{nonlin}","P");
//   leg5->AddEntry(graphRawNoEffCorrNPart,"E_{T}#delta_{TM}/ f_{nonlin}","P");
//   leg5->AddEntry(graphRawAllNPart,"E_{T}/#epsilon f_{nonlin}","P");
//   leg5->AddEntry(graphRawAllNoEffCorrNPart,"E_{T}/f_{nonlin}","P");
//   leg5->Draw();
//   TString corrplotname10 = "/tmp/RawEtVsNPart"+detector+".eps";
//   c10->SaveAs(corrplotname10.Data());


  TCanvas *c11 = new TCanvas("c11","Event by event ET distribution",600,400);
  c11->SetTopMargin(0.02);
  c11->SetRightMargin(0.02);
  c11->SetBorderSize(0);
  c11->SetFillColor(0);
  c11->SetFillColor(0);
  c11->SetBorderMode(0);
  c11->SetFrameFillColor(0);
  c11->SetFrameBorderMode(0);
  c11->SetLogy();
  c11->SetLogx();
  int nbins = 15;
  if(isPhos) nbins = 8;
  int bin = 1;
  ((TH1D*)rawEtAllNoEffCorr[bin])->GetXaxis()->SetTitle("E_{T}^{raw}");
  //if(isPhos){((TH1D*)rawEtAllNoEffCorr[bin])->GetXaxis()->SetRange(1,((TH1D*)rawEtAllNoEffCorr[bin])->GetXaxis()->FindBin(100));}
  //else{((TH1D*)rawEtAllNoEffCorr[bin])->GetXaxis()->SetRange(1,((TH1D*)rawEtAllNoEffCorr[bin])->GetXaxis()->FindBin(200));}
  ((TH1D*)rawEtAllNoEffCorr[bin])->Draw();
  ((TH1D*)rawEtAllNoEffCorr[bin])->SetMaximum(5e4);
  for(bin = 1; bin<=nbins;bin+=1){
    ((TH1D*)rawEtAllNoEffCorr[bin])->Draw("same");
    ((TH1D*)rawEtAllNoEffCorr[bin])->SetLineColor(colors[bin-1]);
    ((TH1D*)rawEtAllNoEffCorr[bin])->SetMarkerColor(colors[bin-1]);
    //((TF1*)dataFits[bin])->SetLineColor(colors[bin-1]);
    //Float_t mostlikely = GetMostLikelyValue((TH1D*)rawEtAllNoEffCorr[bin]);
    //cout<<"cb "<<bin-1<<" average: "<<partialCorrEtValues[bin-1]<<" most likely "<<mostlikely<<" fit "<<((TF1*)dataFits[bin])->GetParameter(1)<<endl;
  }
  TString name9 = "/tmp/RawETDistribution"+detector+".png";
  c11->SaveAs(name9.Data());



  //Create correction file with fractional corrections so that I can use them for cross checks
  ofstream myfile;
  TString texfilename = "SignalFrac"+detector+".dat";
  myfile.open (texfilename.Data());
  //neutron, hadron, kaon, secondaries
  for(int i=0;i<20;i++){
    myfile <<signalFraction[i]<<"\t";
    myfile <<signalFractionError[i]<<"\t";
    myfile <<endl;
  }
  myfile.close();

  //Create correction file with all corrections so that I can use them for cross checks
  ofstream myfile;
  TString texfilename = "allcorrections"+detector+".dat";
  myfile.open (texfilename.Data());
  //neutron, hadron, kaon, secondaries
  for(int i=0;i<20;i++){
    //ETem = 1/facc * scale & 1/fminet [sum(1/eff*1/nonlin) - neutron - hadron - kaon - secondary ]
    myfile <<partialCorrEtValues[i]<<"\t"<<partialCorrEtError[i]<<"\t";
    myfile <<scale<<"\t";
    myfile <<energyscaleerror<<"\t";
    myfile <<minEtCorr[i]<<"\t"<<minEtError[i]<<"\t";
    myfile <<nonLinError[i]<<"\t";
    myfile <<neutronCorr[i]<<"\t"<<neutronError[i]<<"\t";
    myfile <<hadCorr[i]<<"\t"<<hadError[i]<<"\t";
    myfile <<kaonCorr[i]<<"\t"<<kaonError[i]<<"\t";
    myfile <<secondaryCorr[i]<<"\t"<<secondaryError[i]<<endl;
  }
  myfile.close();
  //from spectra calc
  //0.232 $\pm$ 0.029
//   Float_t fem = 0.232;
//   Float_t femerr = 0.029;
  ofstream myfile2;
  ofstream myfile3;
  //EmEtFromSpectra.dat
  TString texfilename = "EmEtFrom"+detector+".dat";
  TString texfilename3 = "TotEtFrom"+detector+".dat";
  myfile2.open (texfilename.Data());
  myfile3.open (texfilename3.Data());
  Int_t maxcb = 10;
  if(isPhos) maxcb = 6;
  else{if(cutset==4){maxcb = 12;}}
  //neutron, hadron, kaon, secondaries
  for(int i=0;i<maxcb;i++){
    int cb1 = i*5;
    int cb2 = (i+1)*5;
    myfile2 <<cb1<<"\t"<<cb2<<"\t";
    myfile2 <<corrEtValues[i]<<"\t";
    myfile2 <<corrEtStatError[i]<<"\t";
    myfile2 <<corrEtError[i]<<"\t";
    myfile2 <<endl;
    
    //Float_t value = corrEtValues[i]/fem;
    //Float_t error = value * TMath::Sqrt(TMath::Power(corrEtError[i]/corrEtValues[i],2)+TMath::Power(femerr/fem,2));
    //cout<<"cb "<<i<<" old "<<value<<" +/- "<<error<<" ("<<error/value*100<<")"<<" new "<<totEtValues[i]<<" +/- "<<totEtError[i]<<" ("<<totEtError[i]/totEtValues[i]<<")"<<endl;
    myfile3 <<cb1<<"\t"<<cb2<<"\t";
    //myfile3 <<value<<"\t";
    //myfile3 <<error<<"\t";
    myfile3 <<totEtValues[i]<<"\t";
    myfile3 <<0<<"\t";//Will be statistical error
    myfile3 <<totEtError[i]<<"\t";
    myfile3 <<endl;
  }
  myfile2.close();
  myfile3.close();

  WriteLatex();
}


void WriteLatex(){
  TString detector = "Emcal";
  string inline;
  //FinalemetsEmcal.dat
  TString finalemetInfileName = "EmEtFrom"+detector+".dat";
  ifstream myfinalemetfile3 (finalemetInfileName.Data());
  Float_t value = 0;
  Float_t error = 0;
  Float_t junk1 = 0;
  Float_t junk2 = 0;
  Int_t i=0;
  if (myfinalemetfile3.is_open()){
    while ( myfinalemetfile3.good() )
      {
	getline (myfinalemetfile3,inline);
	istringstream tmp(inline);
	tmp >> junk1;
	tmp >> junk2;
	tmp >> value;
	tmp >> error;
	if(i<20){
	  //cout<<"value "<<value<<" +/- "<<error<<endl;
	  finalemetCorrEmcal[i] = value;
	  finalemetErrorEmcal[i] = error;
	}
	i++;
      }
    myfinalemetfile3.close();
  }

  detector = "Phos";
  finalemetInfileName = "EmEtFrom"+detector+".dat";
  ifstream myfinalemetfile4 (finalemetInfileName.Data());
  Float_t value = 0;
  Float_t error = 0;
  Int_t i=0;
  if (myfinalemetfile4.is_open()){
    while ( myfinalemetfile4.good() )
      {
	getline (myfinalemetfile4,inline);
	istringstream tmp(inline);
	tmp >> junk1;
	tmp >> junk2;
	tmp >> value;
	tmp >> error;
	if(i<20){
	  finalemetCorrPhos[i] = value;
	  finalemetErrorPhos[i] = error;
	}
	i++;
      }
    myfinalemetfile4.close();
  }

  ofstream myfile3;
  myfile3.open ("EmEt.tex");
  for(int i=0;i<20;i++){
    if(finalemetCorrPhos[i] > 1e-3 || finalemetCorrEmcal[i] > 1e-3){
      if(finalemetCorrPhos[i] > 1e-3){
	TString line = Form("%i-%i\\% & %2.3f $\\pm$ %2.3f & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,finalemetCorrPhos[i],finalemetErrorPhos[i],finalemetCorrEmcal[i],finalemetErrorEmcal[i]);
	myfile3<<line.Data()<<endl;
      }
      else{
	TString line = Form("%i-%i\\% & -- & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,finalemetCorrEmcal[i],finalemetErrorEmcal[i]);
	myfile3<<line.Data()<<endl;
      }
    }
  }
  myfile3.close();






  detector = "Emcal";
  string inline;
  //FinaltotaletsEmcal.dat
  TString finaltotaletInfileName = "TotEtFrom"+detector+".dat";
  ifstream myfinaltotaletfile5 (finaltotaletInfileName.Data());
  value = 0;
  error = 0;
  junk1 = 0;
  junk2 = 0;
  Int_t i=0;
  if (myfinaltotaletfile5.is_open()){
    while ( myfinaltotaletfile5.good() )
      {
	getline (myfinaltotaletfile5,inline);
	istringstream tmp(inline);
	tmp >> junk1;
	tmp >> junk2;
	tmp >> value;
	tmp >> error;
	if(i<20){
	  //cout<<"value "<<value<<" +/- "<<error<<endl;
	  finaltotaletCorrEmcal[i] = value;
	  finaltotaletErrorEmcal[i] = error;
	}
	i++;
      }
    myfinaltotaletfile5.close();
  }

  detector = "Phos";
  finaltotaletInfileName = "TotEtFrom"+detector+".dat";
  ifstream myfinaltotaletfile6 (finaltotaletInfileName.Data());
  value = 0;
  error = 0;
  Int_t i=0;
  if (myfinaltotaletfile6.is_open()){
    while ( myfinaltotaletfile6.good() )
      {
	getline (myfinaltotaletfile6,inline);
	istringstream tmp(inline);
	tmp >> junk1;
	tmp >> junk2;
	tmp >> value;
	tmp >> error;
	if(i<20){
	  finaltotaletCorrPhos[i] = value;
	  finaltotaletErrorPhos[i] = error;
	}
	i++;
      }
    myfinaltotaletfile6.close();
  }
  ofstream myfile4;
  myfile4.open ("TotEtFromCalorimeters.tex");
  for(int i=0;i<20;i++){
    if(finaltotaletCorrPhos[i] > 1e-3 || finaltotaletCorrEmcal[i] > 1e-3){
      if(finaltotaletCorrPhos[i] > 1e-3){
	TString line = Form("%i-%i\\% & %2.3f $\\pm$ %2.3f & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,finaltotaletCorrPhos[i],finaltotaletErrorPhos[i],finaltotaletCorrEmcal[i],finaltotaletErrorEmcal[i]);
	myfile4<<line.Data()<<endl;
      }
      else{
	TString line = Form("%i-%i\\% & -- & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,finaltotaletCorrEmcal[i],finaltotaletErrorEmcal[i]);
	myfile4<<line.Data()<<endl;
      }
    }
  }
  myfile4.close();





}
