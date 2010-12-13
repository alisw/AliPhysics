TH2F *hPtInvMass=0x0;
Int_t nptbins=0;
Double_t *ptbinlimits=0x0;
Double_t *rawsignal=0x0,*rawback=0x0,*sigma=0x0;
Int_t nbinsx;//=hPtInvMass->GetNbinsX();
Int_t nbinsy;//=hPtInvMass->GetNbinsY();
Double_t binwidthpt;//=hPtInvMass->GetYaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
Double_t binwidthInvMss;//=hPtInvMass->GetXaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
Double_t ptmin;//=hPtInvMass->GetYaxis()->GetBinLowEdge(1);
Double_t ptmax;//=hPtInvMass->GetYaxis()->GetBinUpEdge(nbinsy);
Bool_t useFitForSubtraction=kFALSE;
Bool_t subtractSB=kTRUE;
void SetSubtractSB(Bool_t subtract){subtractSB=subtract;}
void SetUseFitForSubtraction(Bool_t useFit){useFitForSubtraction=useFit;}
TCanvas **cPtDistrNoSubtr;
TF1 **fitfunc;
TH1D *hAvRawYieldSpectrum,*hSignal;
TGraphAsymmErrors *grAvRawYieldSpectrum;
TGraphAsymmErrors *grAvPtVSPtmean;

void SetPtBinLimits(const Int_t npt,Double_t ptbinlim[nptbin]){
  if(nptbins>0){
    delete ptbinlimits;ptbinlimits=0x0;
  }
  nptbins=npt;
  ptbinlimits=new Double_t[nptbins+1];
  for(Int_t bin=0;bin<=nptbins;bin++){
    ptbinlimits[bin]=ptbinlim[bin];
  }
  cPtDistrNoSubtr=new TCanvas*[nptbins];
  fitfunc=new TF1*[nptbins];
}
void SetPtBinLimits(TH1 *histo){
  
  if(nptbins>0){
    delete ptbinlimits;ptbinlimits=0x0;
  }
  TAxis *xax=histo->GetXaxis();
  nptbins=xax->GetNbins();
  ptbinlimits=new Double_t[nptbins+1];
  for(Int_t bin=0;bin<=nptbins;bin++){
    ptbinlimits[bin]=xax->GetBinLowEdge(bin+1);
  }
   cPtDistrNoSubtr=new TCanvas*[nptbins];
   fitfunc=new TF1*[nptbins];
}

void SetHistRawSignal(TH1D *hS){hSignal=hS;
 rawsignal=&(hSignal->GetArray())[1];
 return;}
void SetHistRawBack(TH1D *hB){hBack=hB;
 rawback=&(hBack->GetArray())[1];
 hAvRawYieldSpectrum=new TH1D(*hBack);
 hAvRawYieldSpectrum->SetName("hAvRawYieldDistrTotPeak");
 hAvRawYieldSpectrum->SetTitle("Average Raw Yield in bins after subtraction from bin counting"); 
 hAvRawYieldSpectrum->Reset(0);
 hAvRawYieldSpectrum->SetLineColor(kRed);
 return;
}
void SetHistSigma(TH1D *hSig){hSigma=hSig;
 sigma=&(hSigma->GetArray())[1];
 return;
}

void SetHistoMassPt(TH2F *h2){
  hPtInvMass=new TH2F(*h2);
  nbinsx=hPtInvMass->GetNbinsX();
  nbinsy=hPtInvMass->GetNbinsY();
  binwidthpt=hPtInvMass->GetYaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
  binwidthInvMss=hPtInvMass->GetXaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
  ptmin=hPtInvMass->GetYaxis()->GetBinLowEdge(1);
  ptmax=hPtInvMass->GetYaxis()->GetBinUpEdge(nbinsy);

}


TH1F *HistoPtShapeFromData(Int_t ptbin,Int_t rebin=1.){
  if(!hPtInvMass){
    printf("No histo set \n");
    return;
  }
  TString namehist="hPtDistrSB_Bin",funcname;
  namehist+=ptbin;

  TString cname="cPtDistrNoSubBin";
  cname+=ptbin;
  cPtDistrNoSubtr[ptbin]=new TCanvas(cname.Data(),"CPtDistrNoSub",700,700);
  
  printf("After Setting Canva \n");
  Int_t binleftpt=-1,binrightpt=-1;
  TH1F *hPtDistrSB=new TH1F(namehist.Data(),"Side band Pt distribution",nbinsy,ptmin,ptmax);
  namehist="hPtDistrTotPeak_Bin";
  namehist+=ptbin;
  printf("After Setting SB hist \n");
  TH1F *hPtDistrTotPeak=new TH1F(namehist.Data(),"Total Pt distribution in signal mass reagion",nbinsy,ptmin,ptmax);  
  
  printf("Before loop on xbins \n");
  for(Int_t j=1;j<nbinsx;j++){
    if((hPtInvMass->GetXaxis()->GetBinLowEdge(j)<=1.864-5.*sigma[ptbin])||(hPtInvMass->GetXaxis()->GetBinUpEdge(j)>=1.864+5.*sigma[ptbin])){
      for(Int_t k=1;k<nbinsy;k++){
	if(hPtInvMass->GetYaxis()->GetBinLowEdge(k)>=ptbinlimits[ptbin]&&hPtInvMass->GetYaxis()->GetBinUpEdge(k)<=ptbinlimits[ptbin+1]){
	  hPtDistrSB->SetBinContent(k,hPtDistrSB->GetBinContent(k)+hPtInvMass->GetBinContent(j,k));
	}
      }
    }
    else if((hPtInvMass->GetXaxis()->GetBinLowEdge(j)>=1.864-3.*sigma[ptbin])&&(hPtInvMass->GetXaxis()->GetBinUpEdge(j)<=1.864+3.*sigma[ptbin])){
      for(Int_t k=1;k<nbinsy;k++){
	if(hPtInvMass->GetYaxis()->GetBinLowEdge(k)>=ptbinlimits[ptbin]&&hPtInvMass->GetYaxis()->GetBinUpEdge(k)<=ptbinlimits[ptbin+1]){
	  if(binleftpt==-1)binleftpt=k;
	  hPtDistrTotPeak->SetBinContent(k,hPtDistrTotPeak->GetBinContent(k)+hPtInvMass->GetBinContent(j,k));
	  binrightpt=k;
	}
      }
    }    
  }
  hPtDistrSB->Sumw2();
  hPtDistrTotPeak->Sumw2();
  if(rebin!=1){//Rebin and recalculate needed info
    hPtDistrSB->Rebin(rebin);
    hPtDistrTotPeak->Rebin(rebin);
    binleftpt=-1;
    binleftpt=-1;
    for(Int_t k=1;k<=hPtDistrTotPeak->GetNbinsX();k++){      
      if(hPtDistrTotPeak->GetBinLowEdge(k)>=ptbinlimits[ptbin]&&hPtDistrTotPeak->GetBinLowEdge(k+1)<=ptbinlimits[ptbin+1]){
	if(binleftpt==-1)binleftpt=k;
	binrightpt=k;
      }
    }
    binwidthpt=hPtDistrTotPeak->GetBinWidth(1);    
  }

  TH1F *hPtSignalFromSubtraction=new TH1F(*hPtDistrTotPeak);
  namehist="hPtSignalFromSubtraction_Bin";
  namehist+=ptbin;
  hPtSignalFromSubtraction->SetName(namehist.Data());
  hPtSignalFromSubtraction->SetTitle("Signal Pt distr from subtraction");
  //hSignalFromSubtraction->Reset(0);
  if(subtractSB){
    if(!useFitForSubtraction)hPtSignalFromSubtraction->Add(hPtDistrSB,-rawback[ptbin]/hPtDistrSB->Integral());
    else{
      funcname="expofit";
      funcname+=ptbin;
      fitfunc[ptbin]=new TF1(funcname.Data(),"expo",ptbinlimits[ptbin],ptbinlimits[ptbin+1]);        
      TCanvas *cTempt=new TCanvas();
      cTempt->cd();
      hPtDistrSB->Fit(fitfunc[ptbin],"RLEV","",ptbinlimits[ptbin],ptbinlimits[ptbin+1]);
      delete cTempt;
      for(Int_t jb=binleftpt;jb<=binrightpt;jb++){
	hPtSignalFromSubtraction->SetBinContent(jb,hPtSignalFromSubtraction->GetBinContent(jb)-fitfunc[ptbin]->Eval(hPtSignalFromSubtraction->GetBinCenter(jb))/fitfunc[ptbin]->Integral(ptbinlimits[ptbin],ptbinlimits[ptbin+1])*binwidthpt*rawback[ptbin]);
      }    
    }
  }
  
  
  cPtDistrNoSubtr[ptbin]->Divide(1,2);
  cPtDistrNoSubtr[ptbin]->cd(1);
  hPtDistrSB->Draw();
  cPtDistrNoSubtr[ptbin]->cd(2);
  hPtDistrTotPeak->Draw();  
  hPtSignalFromSubtraction->Draw("sames");

  
  return hPtSignalFromSubtraction;

}

void CalculateAveragePt(TH2* hMassPt,TH1D *hB,TH1D *hSigm,TH1D *hS=0x0,Int_t rebin=1,Int_t firstbin=3,Int_t lastbin=8){
  
  SetPtBinLimits(hB);
  SetHistoMassPt((TH2F*)hMassPt);
  SetHistRawBack(hB);
  SetHistSigma(hSigm);
  if(hS)SetHistRawSignal(hS);  
  printf("Histos set \n");
  CalculateAveragePt(rebin,firstbin,lastbin);
}

void CalculateAveragePt(Int_t rebin=1,Int_t firstbin,Int_t lastbin){
  SetUseFitForSubtraction(kFALSE);
  TH1F *hPtSpectra=new TH1F("hPtDistrTotPeak","Total Pt distribution in signal mass reagion",nbinsy,ptmin,ptmax); 
  grAvRawYieldSpectrum=new TGraphAsymmErrors();
  grAvPtVSPtmean=new TGraphAsymmErrors();
  if(rebin!=1)hPtSpectra->Rebin(rebin);
  
  Double_t avptbin=0.;
  for(Int_t bin=firstbin;bin<=lastbin;bin++){    
    //    printf("Before Getting Subtracted Histo \n");
    avptbin=0.;
    TH1F *h=HistoPtShapeFromData(bin,rebin);
    hPtSpectra->Add(h);
    hAvRawYieldSpectrum->SetBinContent(bin+1,h->Integral("width")/hAvRawYieldSpectrum->GetBinWidth(bin+1));
    for(Int_t bbin=1;bbin<=h->GetNbinsX();bbin++){
      avptbin+=h->GetBinCenter(bbin)*h->GetBinContent(bbin)*h->GetBinWidth(bbin);
    }
    avptbin/=h->Integral("width");
    grAvRawYieldSpectrum->SetPoint(bin-firstbin,avptbin,h->Integral("width")/hAvRawYieldSpectrum->GetBinWidth(bin+1));
    grAvRawYieldSpectrum->SetPointError(bin-firstbin,avptbin-ptbinlimits[bin],ptbinlimits[bin+1]-avptbin,0.01,0.01);   

    grAvPtVSPtmean->SetPoint(bin-firstbin,0.5*(ptbinlimits[bin+1]+ptbinlimits[bin]),avptbin);
    grAvPtVSPtmean->SetPointError(bin-firstbin,0.5*(ptbinlimits[bin+1]-ptbinlimits[bin]),0.5*(ptbinlimits[bin+1]-ptbinlimits[bin]),0.01,0.01);
    

    hSignal->SetBinContent(bin+1,hSignal->GetBinContent(bin+1)*hPtSpectra->GetBinWidth(1)/hSignal->GetBinWidth(bin+1));
    hSignal->SetBinError(bin+1,hSignal->GetBinError(bin+1)*hPtSpectra->GetBinWidth(1)/hSignal->GetBinWidth(bin+1));
    //    hPtSpectra->Add(HistoPtShapeFromData(bin,rebin),1./efficiency[bin]/nev*factorTheoryCompare*primaryCharmCorrect[bin]/binwidthpt);//(ptbinlimits[bin+1]-ptbinlimits[bin])
    cPtDistrNoSubtr[bin]->Draw();
  }

   TString nameout="ptCorrection";
   TString namegrRawAvPt="grAvRawVsAvPt",namegrAvPtPtmean="grAvPtMeanPt";
   if(useFitForSubtraction){
     nameout.Append("FitSB");  
     namegrRawAvPt.Append("FitSB");  
     namegrAvPtPtmean.Append("FitSB"); 
   }
   if(!subtractSB){
     namegrRawAvPt.Append("NoSBsub.root");
     namegrAvPtPtmean.Append("NoSBsub.root");
     nameout.Append("NoSBsubtract.root");
   }
   else nameout.Append(".root");
  

   grAvRawYieldSpectrum->SetName(namegrRawAvPt.Data());
   grAvRawYieldSpectrum->SetMarkerStyle(21);
   grAvRawYieldSpectrum->SetMarkerSize(1.2);
   grAvRawYieldSpectrum->SetMarkerColor(kBlue);
   grAvRawYieldSpectrum->SetLineColor(kBlue);

   grAvPtVSPtmean->SetName(namegrAvPtPtmean.Data());
   grAvPtVSPtmean->SetMarkerStyle(21);
   grAvPtVSPtmean->SetMarkerSize(1.2);
   grAvPtVSPtmean->SetMarkerColor(kBlue);
   grAvPtVSPtmean->SetLineColor(kBlue);



  //  TFile *ftheory=TFile::Open("/home/andrea/ALICE/CHARM/ppData_2010/PtSelInBin/predictions_D0.root");
  // TH1D *histoThCompareCentral=(TH1D*)ftheory->Get("hD0Kpipred_central");
  // histoThCompareCentral->Rebin(5);  
 
  
 

  TCanvas *cPtSpectra=new TCanvas("cPtSpectra","cPtSpectra",700,700);
  cPtSpectra->cd();
  hPtSpectra->Draw();
  //hAvRawYieldSpectrum->Draw("same");
  hSignal->SetLineColor(kGreen);
  hSignal->Draw("same");
  grAvRawYieldSpectrum->Draw("p");

  TCanvas *cAvPtVSPtmean=new TCanvas("cAvPtVSPtmean","cAvPtVSPtmean",700,700);
  cAvPtVSPtmean->cd();
  grAvPtVSPtmean->Draw("ap");

  TString nameout="ptCorrection";
  if(useFitForSubtraction)nameout.Append("FitSB");  
  if(!subtractSB)nameout.Append("NoSBsubtract.root");
  else nameout.Append(".root");
  TFile *fout=new TFile(nameout.Data(),"RECREATE");
  fout->cd();
  grAvPtVSPtmean->Write();
  grAvRawYieldSpectrum->Write();
  fout->Close();

  //  histoThCompareCentral->Draw("sames");
}






