
Double_t v2vsMass(Double_t *x, Double_t *par){
  // Fit function for signal+background
  // par[0] = S/B at the mass peak
  // par[1] = mass
  // par[2] = sigma
  // par[3] = v2sig
  // par[4] = v2back

  Double_t fracsig=par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
  Double_t fracbkg=1-fracsig;
  return par[3]*fracsig+par[4]*fracbkg;
}

void Extractv2from2Dhistos(){

  TFile* fil=new TFile("AnalysisResults.root");
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("PWG3_D2H_HFv2");
  TList* lst=(TList*)df->Get("coutputv2D0");
  Int_t minCent=30;
  Int_t maxCent=50;
  TH2F* hMassDphi=0x0;
  for(Int_t iHisC=minCent; iHisC<=maxCent-5; iHisC+=5){    
    TString hisname=Form("hMc2phicentr%d_%d",iHisC,iHisC+5);
    TH2F* htmp=(TH2F*)lst->FindObject(hisname.Data());
    if(iHisC==minCent) hMassDphi=(TH2F*)htmp->Clone("hMassCos2Dphi");
    else hMassDphi->Add(htmp);
    printf("Adding histogram %s\n",hisname.Data());
  }

  TH1F* hMass=(TH1F*)hMassDphi->ProjectionY();

  Double_t sigmaRangeForSig=3.;
  Double_t sigmaRangeForBkg=4.5;

  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  TCanvas* c1=new TCanvas("c1");
  c1->Divide(2,2);
  c1->cd(1);
  hMassDphi->Draw("colz");
  c1->cd(2);

  Int_t nMassBins=hMass->GetNbinsX();
  Int_t hMinBin=3;
  Int_t hMaxBin=nMassBins-2;
  Double_t hmin=hMass->GetBinLowEdge(hMinBin);
  Double_t hmax=hMass->GetBinLowEdge(hMaxBin)+hMass->GetBinWidth(hMaxBin);
  Int_t factor4refl=0;
  Float_t massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  AliHFMassFitter* fitter=new AliHFMassFitter(hMass,hmin,hmax,2,0,0);
  fitter->SetReflectionSigmaFactor(factor4refl);
  fitter->SetInitialGaussianMean(massD);
  Bool_t out=fitter->MassFitter(0);
  if(!out) return;
  fitter->DrawHere(gPad);
  Double_t sigfitter,esigfitter;
  fitter->Signal(sigmaRangeForSig, sigfitter,esigfitter);
  TH1F* hCos2PhiBkgLo=0x0;
  TH1F* hCos2PhiBkgHi=0x0;
  TH1F* hCos2PhiBkgLoScal=0x0;
  TH1F* hCos2PhiBkgHiScal=0x0;
  TH1F* hCos2PhiBkgAver=0x0;
  TH1F* hCos2PhiSigReg=0x0;
  TH1F* hCos2PhiSig=0x0;
  Double_t mass=fitter->GetMean();
  Double_t sigma=fitter->GetSigma();
  TF1* fB1=fitter->GetBackgroundFullRangeFunc();
  TF1* fB2=fitter->GetBackgroundRecalcFunc();
  TF1* fSB=fitter->GetMassFunc();
  Double_t minMassSig=mass-sigmaRangeForSig*sigma;
  Double_t maxMassSig=mass+sigmaRangeForSig*sigma;
  Int_t minBinSig=hMass->FindBin(minMassSig);
  Int_t maxBinSig=hMass->FindBin(maxMassSig);
  Double_t minMassSigBin=hMass->GetBinLowEdge(minBinSig);
  Double_t maxMassSigBin=hMass->GetBinLowEdge(maxBinSig)+hMass->GetBinWidth(maxBinSig);
  printf("Signal Fit Limits = %f %f\n",minMassSigBin,maxMassSigBin);
  Double_t maxMassBkgLow=mass-sigmaRangeForBkg*sigma;
  Int_t minBinBkgLow=hMinBin;
  Int_t maxBinBkgLow=hMass->FindBin(maxMassBkgLow);
  Double_t minMassBkgLowBin=hmin;
  Double_t maxMassBkgLowBin=hMass->GetBinLowEdge(maxBinBkgLow)+hMass->GetBinWidth(maxBinBkgLow);
  Double_t minMassBkgHi=mass+sigmaRangeForBkg*sigma;
  Int_t minBinBkgHi=hMass->FindBin(minMassBkgHi);
  Int_t maxBinBkgHi=hMaxBin;
  Double_t minMassBkgHiBin=hMass->GetBinLowEdge(minBinBkgHi);
  Double_t maxMassBkgHiBin=hmax;
  printf("BKG Fit Limits = %f %f  && %f %f\n",minMassBkgLowBin,maxMassBkgLowBin,minMassBkgHiBin,maxMassBkgHiBin);
  Double_t bkgSig=fB2->Integral(minMassSigBin,maxMassSigBin);
  Double_t bkgLow=fB2->Integral(minMassBkgLowBin,maxMassBkgLowBin);
  Double_t bkgHi=fB2->Integral(minMassBkgHiBin,maxMassBkgHiBin);
  printf("Background integrals = %f %f %f\n",bkgLow,bkgSig,bkgHi);
  TBox* bleft=new TBox(minMassBkgLowBin,0.,maxMassBkgLowBin,hMass->GetMaximum());
  bleft->SetFillColor(kRed+1);
  bleft->SetFillStyle(3002);
  bleft->Draw();
  TBox* bright=new TBox(minMassBkgHiBin,0.,maxMassBkgHiBin,hMass->GetMaximum());
  bright->SetFillColor(kBlue+1);
  bright->SetFillStyle(3002);
  bright->Draw();
  TBox* bsig=new TBox(minMassSigBin,0.,maxMassSigBin,hMass->GetMaximum()*2);
  bsig->SetFillColor(1);
  bsig->SetFillStyle(3002);
  bsig->Draw();

  TH1F* hCos2PhiBkgLo=(TH1F*)hMassDphi->ProjectionX("hCos2PhiBkgLo",minBinBkgLow,maxBinBkgLow);
  TH1F* hCos2PhiBkgHi=(TH1F*)hMassDphi->ProjectionX("hCos2PhiBkgHi",minBinBkgHi,maxBinBkgHi);
  TH1F* hCos2PhiSigReg=(TH1F*)hMassDphi->ProjectionX("hCos2PhiBkgSig",minBinSig,maxBinSig);
  hCos2PhiBkgLo->Rebin(4);
  hCos2PhiBkgHi->Rebin(4);
  hCos2PhiSigReg->Rebin(4);
  hCos2PhiSigReg->SetLineWidth(2);
  hCos2PhiBkgLo->SetLineWidth(2);
  hCos2PhiBkgHi->SetLineWidth(2);
  hCos2PhiBkgLo->SetLineColor(kRed+1);
  hCos2PhiBkgHi->SetLineColor(kBlue+1);
  TH1F* hCos2PhiBkgLoScal=(TH1F*)hCos2PhiBkgLo->Clone("hCos2PhiBkgLoScal");
  hCos2PhiBkgLoScal->Scale(bkgSig/bkgLow);
  TH1F* hCos2PhiBkgHiScal=(TH1F*)hCos2PhiBkgHi->Clone("hCos2PhiBkgHiScal");
  hCos2PhiBkgHiScal->Scale(bkgSig/bkgHi);
  hCos2PhiBkgLoScal->SetLineWidth(2);
  hCos2PhiBkgHiScal->SetLineWidth(2);
  hCos2PhiBkgLoScal->SetLineColor(kRed+1);
  hCos2PhiBkgHiScal->SetLineColor(kBlue+1);
  TH1F* hCos2PhiBkgAver=(TH1F*)hCos2PhiBkgLoScal->Clone("hCos2PhiBkgAver");
  hCos2PhiBkgAver->Add(hCos2PhiBkgHiScal);
  hCos2PhiBkgAver->Scale(0.5);
  hCos2PhiBkgAver->SetLineWidth(2);
  hCos2PhiBkgAver->SetLineColor(kGreen+1);
  TH1F* hCos2PhiSig=(TH1F*)hCos2PhiSigReg->Clone("hCos2PhiSig");
  hCos2PhiSig->Add(hCos2PhiBkgAver,-1.);   
  printf("v2=%f +- %f\n",hCos2PhiSig->GetMean(),hCos2PhiSig->GetMeanError());
  c1->cd(3);
  hCos2PhiSigReg->Draw();
  hCos2PhiBkgLoScal->Draw("same");
  hCos2PhiBkgHiScal->Draw("same");
  hCos2PhiBkgAver->Draw("same");
  TLegend* leg0=new TLegend(0.3,0.6,0.75,0.89);
  leg0->SetFillColor(0);
  TLegendEntry* ent=leg0->AddEntry(hCos2PhiBkgSig,"Signal region","L");
  ent->SetTextColor(hCos2PhiBkgSig->GetLineColor());
  ent=leg0->AddEntry(hCos2PhiBkgLoScal,"Left side band","L");
  ent->SetTextColor(hCos2PhiBkgLoScal->GetLineColor());
  ent=leg0->AddEntry(hCos2PhiBkgHiScal,"Right side band","L");
  ent->SetTextColor(hCos2PhiBkgHiScal->GetLineColor());
  ent=leg0->AddEntry(hCos2PhiBkgAver,"Average of side bands","L");
  ent->SetTextColor(hCos2PhiBkgAver->GetLineColor());
  leg0->Draw();
  c1->cd(4);
  hCos2PhiSig->Draw("EP");
  TPaveText* t0= new TPaveText(0.15,0.70,0.45,0.89,"NDC");
  t0->SetFillColor(0);
  t0->AddText(Form("v2=%.3f+-%.3f\n",hCos2PhiSig->GetMean(),hCos2PhiSig->GetMeanError()));
  t0->Draw();

  printf("Signal from mass fitter = %f  Signal from subracted histo= %f\n",
	 sigfitter,hCos2PhiSig->Integral());

  Int_t npars=fSB->GetNpar();
  Double_t sigma=fSB->GetParameter(npars-1);
  Double_t mass=fSB->GetParameter(npars-2);
  Double_t integr=fSB->GetParameter(npars-3);
  Double_t sOverAll=(fSB->Eval(mass)-fB2->Eval(mass))/fSB->Eval(mass);
  printf("mass=%f  S+B=%f   bkg=%f S/(S+B)=%f\n",mass,fSB->Eval(mass),fB2->Eval(mass),sOverAll);
  printf("Number of parameters: %d. Signal params: %f %f %f\n",npars,mass,sigma,integr);
  TF1* fSig=new TF1("fSig","[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])",hMassDphi->GetYaxis()->GetXmin(),hMassDphi->GetYaxis()->GetXmax());
  fSig->SetParameters(integr,mass,sigma);
 
  TH1F* hAverCos2Phi=new TH1F("hAverCos2Phi","",hMassDphi->GetNbinsY(),hMassDphi->GetYaxis()->GetXmin(),hMassDphi->GetYaxis()->GetXmax());
  TH1F* hFractionSig=new TH1F("hFractionSig","",hMassDphi->GetNbinsY(),hMassDphi->GetYaxis()->GetXmin(),hMassDphi->GetYaxis()->GetXmax());
  TH1F* hFractionBkg=new TH1F("hFractionBkg","",hMassDphi->GetNbinsY(),hMassDphi->GetYaxis()->GetXmin(),hMassDphi->GetYaxis()->GetXmax());

  for(Int_t iBin=1; iBin<=hMassDphi->GetNbinsY(); iBin++){
    TH1F* htemp=(TH1F*)hMassDphi->ProjectionX("htemp",iBin,iBin);
    hAverCos2Phi->SetBinContent(iBin,htemp->GetMean());
    hAverCos2Phi->SetBinError(iBin,htemp->GetMeanError());
    Double_t sig=fSig->Eval(hFractionSig->GetBinCenter(iBin));
    Double_t bkg=fB2->Eval(hFractionSig->GetBinCenter(iBin));
    if(bkg<1 && sig<1){
      hFractionSig->SetBinContent(iBin,0.);
      hFractionSig->SetBinError(iBin,0.);
      hFractionBkg->SetBinContent(iBin,1.);
      hFractionBkg->SetBinError(iBin,0.);
    }else{
      Double_t fracs=sig/(sig+bkg);
      Double_t fracb=bkg/(sig+bkg);
      Double_t efracs=0.;//TMath::Sqrt(fracs*(1.-fracs)/(sig+bkg));
      Double_t efracb=0.;//TMath::Sqrt(fracb*(1.-fracb)/(sig+bkg));
       
      hFractionSig->SetBinContent(iBin,fracs),
      hFractionSig->SetBinError(iBin,efracs);
      hFractionBkg->SetBinContent(iBin,fracb);      
      hFractionBkg->SetBinError(iBin,efracb);
    }
    delete htemp;
  }
  
  TF1* fv2=new TF1("fv2",v2vsMass,hMassDphi->GetYaxis()->GetXmin(),hMassDphi->GetYaxis()->GetXmax(),5);
  fv2->SetParameter(0,sOverAll);
  fv2->SetParameter(1,mass);
  fv2->SetParameter(2,sigma);
  fv2->SetParameter(3,0.2);
  fv2->SetParameter(4,0.2);
  fv2->FixParameter(0,sOverAll);
  fv2->FixParameter(1,mass);
  fv2->FixParameter(2,sigma);

  hAverCos2Phi->Rebin(2);
  hAverCos2Phi->Scale(0.5);

  TCanvas* c2=new TCanvas("c2");
  c2->Divide(2,2);
  c2->cd(1);
  hMassDphi->Draw("colz");
  c2->cd(2);
  hMass->Rebin(2);
  hMass->SetMinimum(0.);
  hMass->SetMarkerStyle(20);
  hMass->Draw("E");
  fSB->Draw("same");
  fSig->Draw("same");
  fB2->Draw("same");
  c2->cd(3);
  hFractionSig->SetMaximum(1.2);
  hFractionSig->Draw();
  hFractionSig->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  hFractionSig->GetYaxis()->SetTitle("Fraction");
  hFractionBkg->SetLineColor(2);
  hFractionBkg->Draw("same");
  TLegend* leg1=new TLegend(0.15,0.15,0.35,0.35);
  leg1->SetFillColor(0);
  ent=leg1->AddEntry(hFractionSig,"S/(S+B)","L");
  ent->SetTextColor(hFractionSig->GetLineColor());
  ent=leg1->AddEntry(hFractionBkg,"B/(S+B)","L");
  ent->SetTextColor(hFractionBkg->GetLineColor());
  leg1->Draw();
  c2->cd(4);
  hAverCos2Phi->Fit(fv2);
  hAverCos2Phi->GetXaxis()->SetTitle("Mass (GeV/c^2)");
  hAverCos2Phi->GetYaxis()->SetTitle("v_2^{obs}");
  TPaveText* t1= new TPaveText(0.55,0.70,0.89,0.89,"NDC");
  t1->SetFillColor(0);
  t1->AddText(Form("v2sig=%.3f+-%.3f\n",fv2->GetParameter(3),fv2->GetParError(3)));
  t1->AddText(Form("v2bkg=%.3f+-%.3f\n",fv2->GetParameter(4),fv2->GetParError(4)));
  t1->Draw();
  printf("v2(from fit)=%f+-%f\n",fv2->GetParameter(3),fv2->GetParError(3));

}
