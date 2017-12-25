TH1D* ReflectHisto(TH1D *h,Double_t scale);

void EvaluateMCClosModulations(Int_t binDmin, Int_t binDmax, Double_t thrMin, Double_t thrMax) {

  TFile fRatio(Form("Output_Root/MCClosure_Dzero_Canvas_PtIntBins%dto%d_PoolInt_thr%.1fto%.1f.root",binDmin,binDmax,thrMin,thrMax));
  TFile fAmplReco(Form("Output_Root/AzimCorrDistr_Dzero_Canvas_PtIntBins%dto%d_PoolInt_thr%.1fto%.1f_Superimposed_Reco.root",binDmin,binDmax,thrMin,thrMax));

  Double_t FPrompt = 0;
  if(binDmin==4) FPrompt = 0.935;
  if(binDmin==6) FPrompt = 0.934;
  if(binDmin==9) FPrompt = 0.927;
  if(binDmin==12) FPrompt = 0.916;

  TCanvas *cRatio = (TCanvas*)fRatio.Get(Form("cMCClosure_Dzero_%dto%d_%.1fto%.1f",binDmin,binDmax,thrMin,thrMax));
  TH1D *hRatio = (TH1D*)cRatio->FindObject("h1D_MCClosure_Orig5");

  TCanvas *cAmpl = (TCanvas*)fAmplReco.Get(Form("cOutSuperimp_%dto%d_%.1fto%.1f_Reco",binDmin,binDmax,thrMin,thrMax));
  TH1D *hAmplB = (TH1D*)cAmpl->FindObject("h1D_SubtrNorm_Orig5");
  TH1D *hAmplC = (TH1D*)cAmpl->FindObject("h1D_SubtrNorm_Orig4");

  TH1D *hAmplBRefl = ReflectHisto(hAmplB,0.5);
  TH1D *hAmplCRefl = ReflectHisto(hAmplC,0.5);
  Double_t relAmplC[6] = {0,0,0,0,0,0};
  Double_t relAmplB[6] = {0,0,0,0,0,0};
  Double_t recoKineVal[6] = {0,0,0,0,0,0};
  Double_t modul[6] = {0,0,0,0,0,0};

  printf("***BIN %d to %d, TRACK %1.1f TO %1.1f***\n",binDmin,binDmax,thrMin,thrMax);

  TF1 *funFit = new TF1("funFit","[0]",0,TMath::Pi()*6./18.);
  hRatio->Fit(funFit);
  Double_t fitVal = funFit->GetParameter(0);
  //Equation is:
  //CorrData = OrigData*modul = OrigData*[RelAmplC*fPrompt + RelAmplB*fPrompt/RecoKineVal]
  //Here we extract 'modul' in the first five dPhi bins, that we pass to the framework to be used for data scaling (and for uncertainty, i.e. its sqrt(12) relative amplitude bilateral)



  for(int i=0; i<6; i++) {
    recoKineVal[i] = hRatio->GetBinContent(i+1) - (fitVal-1);
    relAmplC[i] = hAmplCRefl->GetBinContent(i+1)/(hAmplCRefl->GetBinContent(i+1)*FPrompt + hAmplBRefl->GetBinContent(i+1)*(1-FPrompt));
    relAmplB[i] = hAmplBRefl->GetBinContent(i+1)/(hAmplCRefl->GetBinContent(i+1)*FPrompt + hAmplBRefl->GetBinContent(i+1)*(1-FPrompt));
    modul[i] = relAmplC[i]*FPrompt + relAmplB[i]*(1-FPrompt)/recoKineVal[i];
    
    printf("Bin%d) MODUL = %1.5f\t (Reco/Kine-fitVal = %1.4f, FPrompt = %1.3f, Ampl_ratio C,B = %1.4f, %1.4f)\n",i+1,modul[i],recoKineVal[i],FPrompt,relAmplC[i],relAmplB[i]);
  }

  return;
}


//________________________________________________________________________________________________
TH1D* ReflectHisto(TH1D *h,Double_t scale){
  
  TH1D *h2=new TH1D(Form("%sReflected",h->GetName()),Form("%sReflected",h->GetName()),h->GetNbinsX()/2.,0.,TMath::Pi());
  for(Int_t j=1;j<=h->GetNbinsX();j++){
    Double_t x=h->GetBinCenter(j);
    Double_t y0=h->GetBinContent(j);
    Double_t ey0=h->GetBinError(j);
    Int_t j2;
    if(x>0&&x<TMath::Pi())j2=h2->FindBin(x);
    else if(x<0)j2=h2->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h2->FindBin(2.*TMath::Pi()-x);
    else {
      printf("Point %d excluded \n",j);
      continue;
    }
    Double_t y=h2->GetBinContent(j2);
    Double_t ey=h2->GetBinError(j2);
    h2->SetBinContent(j2,(y+y0));
    h2->SetBinError(j2,TMath::Sqrt(ey0*ey0+ey*ey));  
  }
  h2->Scale(scale);


  return h2;
}
