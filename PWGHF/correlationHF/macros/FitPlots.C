//--------------------------------------
//
//        Macro used to fit azimuthal correlations
//        some methods are obsolete (used before core of the fitting code was moved in a dedicated class)--> will be updates
//        A. Rossi (andrea.rossi@cern.ch)
//
//---------------------------------------
TF1 *fitFunction;
Double_t baseline=-999.;
Double_t errbaseline=-999.;
Bool_t isReflected=kFALSE;// obsolete parameter, will have to remove it
Double_t GetBaseline(Double_t &errbase){
  errbase=errbaseline;
  return baseline;
}

TF1 *GetFitFunction(){return fitFunction;}

void AddTextInfo(TPad *pd,Double_t ptDmin,Double_t ptDmax,Double_t ptAssmin,Double_t ptAssmax=-1,Double_t deltaEta=1.,Int_t system=0){
  AddTextInfo( ptDmin, ptDmax, ptAssmin, ptAssmax,deltaEta,system);
}

void AddTextInfo(Double_t ptDmin,Double_t ptDmax,Double_t ptAssmin,Double_t ptAssmax=-1,Double_t deltaEta=1.,Int_t system=0,TPad *pd = NULL){

  TLatex *tlTitle=new TLatex(0.15,0.85,"D meson-charged hadron azimuthal correlations");
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextSize(0.03);

  TLatex *tlMesons=new TLatex(0.17,0.8,"D^{0},D^{*+} average");
  tlMesons->SetNDC();
  tlMesons->Draw();
  tlMesons->SetTextSize(0.03);

  TLatex *tlcoll;
  if(system==0)tlcoll=new TLatex(0.17,0.76,"pp, #sqrt{s}=7 TeV");
  else if(system==1)tlcoll=new TLatex(0.17,0.76,"p-Pb, #sqrt{s}=5.02 TeV");
  else if(system==2)tlcoll=new TLatex(0.17,0.76,"Pb-Pb, #sqrt{s}=2.76 TeV");
  tlcoll->SetNDC();
  tlcoll->Draw();
  tlcoll->SetTextSize(0.03);
  

  TLatex *tlKine;
  if(ptAssmax>10)tlKine=new TLatex(0.17,0.71,Form("%.0f<p_{T}^{D}<%.0f GeV/c,p_{T}^{h}>%.1f GeV/c",ptDmin,ptDmax,ptAssmin));
  else if(ptAssmax>0.&&ptAssmax<10)tlKine=new TLatex(0.17,0.71,Form("%.0f<p_{T}^{D}<%.0f GeV/c,%.1f<p_{T}^{h}<%.1f GeV/c ",ptDmin,ptDmax,ptAssmin,ptAssmax));
  else tlKine=new TLatex(0.17,0.71,Form("%.0f<p_{T}^{D}<%.0f GeV/c, p_{T}^{h}>%.1f GeV/c ",ptDmin,ptDmax,ptAssmin));
  tlKine->Draw();  
  tlKine->SetTextSize(0.03);

//   TLatex *tlKineAss;
//   if(ptAssmax>0.)tlKineAss=new TLatex(0.17,0.67,Form("%.1f<p_{T}^{h}<%.1f",ptAssmin,ptAssmax));
//   else tlKineAss=new TLatex(0.17,0.67,Form("p_{T}^{h}>%.1f",ptAssmin));
//   tlKineAss->SetNDC();
//   tlKineAss->Draw();  
//   tlKineAss->SetTextSize(0.03);

  TLatex *tlKineDeltaEta=new TLatex(0.17,0.67,Form("|#Delta#eta|<%.1f",deltaEta));
  tlKineDeltaEta->SetNDC();
  tlKineDeltaEta->Draw();  
  tlKineDeltaEta->SetTextSize(0.03);
  
}

void DrawLegendWithParameters(TPad *pd,TF1 *f){
  DrawLegendWithParameters(f,pd);
}

void DrawLegendWithParameters(TF1 *f,TPad *pd = NULL){
  Int_t nsy=f->GetParNumber("NS Y");
  Int_t asy=f->GetParNumber("AS Y");
  Int_t nss=f->GetParNumber("NS #sigma");
  Int_t ass=f->GetParNumber("AS #sigma");
  Int_t bas=f->GetParNumber("ped");

  TPaveText *pvStatTests1=new TPaveText(0.61,0.6,0.91,0.82,"NDC");
  pvStatTests1->SetFillStyle(0);
  pvStatTests1->SetBorderSize(0);
  pvStatTests1->AddText(0.,0.87,Form("#chi^{2}/ndf = %.1f/%d ",f->GetChisquare(),f->GetNDF()));
  pvStatTests1->AddText(0.,0.69,Form("NS Y = %.2f#pm%.2f ",f->GetParameter(nsy),f->GetParError(nsy)));
  pvStatTests1->AddText(0.,0.51,Form("NS #sigma = %.2f#pm%.2f ",f->GetParameter(nss),f->GetParError(nss)));
  pvStatTests1->AddText(0.,0.33,Form("AS Y = %.2f#pm%.2f ",f->GetParameter(asy),f->GetParError(asy)));
  pvStatTests1->AddText(0.,0.15,Form("AS #sigma = %.2f#pm%.2f ",f->GetParameter(ass),f->GetParError(ass)));
  
  if(baseline<-998.){
    pvStatTests1->AddText(0.,0.,Form("baseline = %.2f#pm%.2f ",f->GetParameter(bas),f->GetParError(bas)));
  }
  else{     pvStatTests1->AddText(0.,0.,Form("baseline = %.2f#pm%.2f ",baseline,errbaseline));

}
  pvStatTests1->SetTextColor(kBlack);
  pvStatTests1->SetTextFont(42);
  //  pvStatTests1->SetNDC();
  pvStatTests1->SetTextSize(0.025);

  pvStatTests1->Draw("same");
}

TF1 *GetFitFunctionConst2Gaus(){
  TF1 *f=new TF1("myfunc","[0]+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))",-1.7,4.7);


  f->SetParameter(0,0.6);
  f->SetParameter(3,0.3);
  f->SetParLimits(3,0,3.5);

  f->SetParameter(6,0.3);
  f->SetParLimits(6,0,3.5);

  f->SetParameter(5,3.14);
  f->SetParLimits(2,2.14,4.14);
 
 f->SetParameter(2,0);
 f->SetParLimits(2,-1,1);
 f->SetParameter(2,0.);
 f->SetParameter(1,3);
 f->SetParLimits(1,0,999.);
 f->SetParameter(4,2);
 f->SetParLimits(4,0,999.);

 f->SetParName(0,"ped");
 f->SetParName(1,"NS Y");
 f->SetParName(2,"NS mean");
 f->SetParName(3,"NS #sigma");
 f->SetParName(4,"AS Y");
 f->SetParName(5,"AS mean");
 f->SetParName(6,"AS #sigma");
 return f;
}

TF1 *GetGaussianPeriodicity(){
  TF1 *f=new TF1("fGausPeriodicity","[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))",-1.7,5.2);

  f->SetParameter(0,3.);
  f->SetParameter(1,0.);  
  f->SetParameter(2,0.3);


  return f;
}

TF1 *GetGaussian(){
  TF1 *f=new TF1("fGaus","[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))",-1.7,5.2);

  f->SetParameter(0,3.);
  f->SetParameter(1,0.);  
  f->SetParameter(2,0.3);


  return f;
}

TF1 *GetFitFunctionConst2GausPeriodicity(){
  TF1 *f=new TF1("f2GausPeriodicity","[0]+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([3]*[3]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x+2.*TMath::Pi()-[5])*(x+2.*TMath::Pi()-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-2.*TMath::Pi()-[5])*(x-2.*TMath::Pi()-[5])/2./([6]*[6]))",-1.7,5.2);
  

  f->SetParLimits(0,0.,50);
  f->SetParameter(0,0.6);

  f->SetParameter(3,0.3);
  f->SetParLimits(3,0,3.5);

  f->SetParameter(6,0.3);
  f->SetParLimits(6,0,3.5);

  f->SetParLimits(5,2.85,3.55);
  f->SetParameter(5,TMath::Pi());

 f->SetParLimits(2,-0.55,0.55);
 f->SetParameter(2,0.);

 f->SetParameter(1,3);
 f->SetParameter(4,2);


 f->SetParName(0,"ped");
 f->SetParName(1,"NS Y");
 f->SetParName(2,"NS mean");
 f->SetParName(3,"NS #sigma");
 f->SetParName(4,"AS Y");
 f->SetParName(5,"AS mean");
 f->SetParName(6,"AS #sigma");
 f->SetDrawOption("9");
 return f;
}

TF1 *GetFitFunctionConst3GausPeriodicity(){
  TF1 *f=new TF1("f3GausPeriodicity","[0]+[1]*([7]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[7]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([3]*[3]))+[7]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([3]*[3]))+(1.-[7])/TMath::Sqrt(2.*TMath::Pi())/[8]*TMath::Exp(-(x-[2])*(x-[2])/2./([8]*[8]))+(1.-[7])/TMath::Sqrt(2.*TMath::Pi())/[8]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([8]*[8]))+(1.-[7])/TMath::Sqrt(2.*TMath::Pi())/[8]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([8]*[8])))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x+2.*TMath::Pi()-[5])*(x+2.*TMath::Pi()-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-2.*TMath::Pi()-[5])*(x-2.*TMath::Pi()-[5])/2./([6]*[6]))",-1.7,5.2);
  

  f->SetParLimits(0,0.,50);
  f->SetParameter(0,0.6);

  f->SetParameter(3,0.3);
  f->SetParLimits(3,0,3.5);

  f->SetParameter(6,0.3);
  f->SetParLimits(6,0,3.5);

  f->SetParameter(8,0.3);
  f->SetParLimits(8,0,3.5);

  f->SetParLimits(5,2.85,3.55);
  f->SetParameter(5,TMath::Pi());

  f->SetParLimits(2,-0.55,0.55);
  f->SetParameter(2,0.);
  
  f->SetParameter(1,3);
  f->SetParameter(4,2);
  
  f->SetParLimits(7,0.,1.);
  
  f->SetParName(0,"ped");
  f->SetParName(1,"NS Y");
  f->SetParName(2,"NS mean 1g");
  f->SetParName(3,"NS #sigma 1g");
  f->SetParName(4,"AS Y");
  f->SetParName(5,"AS mean");
  f->SetParName(6,"AS #sigma");
  f->SetParName(7,"fract 1g");
  f->SetParName(8,"NS #sigma 2g");

  f->SetDrawOption("9");
 return f;
}


TF1* FitPlotsShort(TH1D *h,Int_t fitFunc=1,Int_t fixBase=0,Int_t fixMean=0,Bool_t refl=kFALSE,Double_t rangeTransvMin=0.25*TMath::Pi(),Double_t rangeTransvMax=0.5*TMath::Pi(),Double_t minDpt=0,Double_t maxDpt=0,Double_t minAsspt=0,Double_t maxAsspt=0){
// - fitFunc=0: 2 gaussian + const baseline w/o periodicity (to be avoided, only for checks)
//               =1: 2 gauss + const baseline + periodicity
//               =2: 3 gaus (2 on the NS) + const baseline with periodicity. Useful for fitting MC templates
//
// - fixBase=0: baseline free
//                  =1 fix the baseline to the minimum of the histogram
//                  <0 fix the baseline to the weighted average of the abs(fixbaseline) lower points
//                  =2 :zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
//                  =4 :zyam at pi/2. Fix the baseline averaging the 4 points around +-pi/2 value
//
// - fixMean=0 : means free
//                   =1 : NS mean fixed to0, AS mean free
//                   =2: AS mean fixed to pi, NS mean free
//                 =3: NS mean fixed to 0, AS mean to pi

  Double_t nsybc, ensybc,asybc, easybc;
  return FitPlots(h,fitFunc, fixBase, fixMean, nsybc, ensybc, asybc, easybc,refl,rangeTransvMin,rangeTransvMax,minDpt,maxDpt,minAsspt,maxAsspt);
}

TF1 *FitPlots(TH1D *h,Int_t fitFunc=1,Int_t fixBase=0,Int_t fixMean=0,Double_t &nsybc,Double_t &ensybc,Double_t &asybc,Double_t &easybc,Bool_t refl,Double_t rangeTransvMin=0.25*TMath::Pi(),Double_t rangeTransvMax=0.5*TMath::Pi(),Double_t minDpt=0,Double_t maxDpt=0,Double_t minAsspt=0,Double_t maxAsspt=0){//
// - fitFunc=0: 2 gaussian + const baseline w/o periodicity (to be avoided, only for checks)
//               =1: 2 gauss + const baseline + periodicity
//               =2: 3 gaus (2 on the NS) + const baseline with periodicity. Useful for fitting MC templates
//
// - fixBase=0: baseline free
//                  =1 fix the baseline to the minimum of the histogram
//                  <0 fix the baseline to the weighted average of the abs(fixbaseline) lower points
//                  =2 :zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
//                  =4 :zyam at pi/2. Fix the baseline averaging the 4 points around +-pi/2 value
//                  =5 :pedestal fixed at avarage of the points in the transverse region |(pi/4 - pi2)|
//                  =6 :fit to external value
//
// - fixMean=0 : means free
//                   =1 : NS mean fixed to0, AS mean free
//                   =2: AS mean fixed to pi, NS mean free
//                 =3: NS mean fixed to 0, AS mean to pi
//
// nsybc, ensybc,asybc, easybc are the yield (and their errors) obtained via bin counting in the NS (-1.5<DeltaPhi<1.5) and AS (pi-1.5<deltaPhi<pi+1.5) regions

  Double_t min=0,max=TMath::Pi();
  if(!refl){
    min=-0.5*TMath::Pi();
    max=1.5*TMath::Pi();
  }
  AliHFCorrFitter *corrfitter=new AliHFCorrFitter((TH1F*)h,min,max,kTRUE);
  corrfitter->SetHistoIsReflected(refl);
  corrfitter->SetFuncType(fitFunc);
  corrfitter->SetFixBasetype(fixBase);
  if(fixBase==5){
    corrfitter->SetBaselineEstimationRange(rangeTransvMin,rangeTransvMax);
  }
  corrfitter->SetFixMeanType(fixMean);
  corrfitter->SetPtRanges(minDpt,maxDpt,minAsspt,maxAsspt);
  corrfitter->Fitting(kTRUE);
  nsybc=corrfitter->GetBinCountingYields(ensybc,asybc,easybc);
  fitFunction=corrfitter->GetFitFunction();// obsolete needs, will have to be removed
  baseline=corrfitter->GetPedestal();
  errbaseline=corrfitter->GetPedestalError();
  return fitFunction;
}


TF1* FitPlotsObsolete(TH1D *h,Int_t fitFunc=1,Int_t fixBase=0,Int_t fixMean=0,Double_t &nsybc,Double_t &ensybc,Double_t &asybc,Double_t &easybc,Bool_t refl){//
// - fitFunc=0: 2 gaussian + const baseline w/o periodicity (to be avoided, only for checks)
//               =1: 2 gauss + const baseline + periodicity
//               =2: 3 gaus (2 on the NS) + const baseline with periodicity. Useful for fitting MC templates
//
// - fixBase=0: baseline free
//                  =1 fix the baseline to the minimum of the histogram
//                  <0 fix the baseline to the weighted average of the abs(fixbaseline) lower points
//                  =2 :zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
//                  =4 :zyam at pi/2. Fix the baseline averaging the 4 points around +-pi/2 value
//                  =5 :pedestal fixed at avarage of the points in the transverse region |(pi/4 - pi2)|
//
// - fixMean=0 : means free
//                   =1 : NS mean fixed to0, AS mean free
//                   =2: AS mean fixed to pi, NS mean free
//                 =3: NS mean fixed to 0, AS mean to pi
//
// nsybc, ensybc,asybc, easybc are the yield (and their errors) obtained via bin counting in the NS (-1.5<DeltaPhi<1.5) and AS (pi-1.5<deltaPhi<pi+1.5) regions
  isReflected=refl;
  baseline=-999.;
  errbaseline=-999.;

  if(fitFunc==0){
    fitFunction=GetFitFunctionConst2Gaus();
  }
  else if(fitFunc==1)fitFunction=GetFitFunctionConst2GausPeriodicity();
  else if(fitFunc==2)fitFunction=GetFitFunctionConst3GausPeriodicity();
  
  if(fixBase==1){
    Double_t min=1.e10;
    Int_t k=-1;
    for(Int_t j=1;j<=h->GetNbinsX();j++){
      if(h->GetBinContent(j)<min){
	min=h->GetBinContent(j);
	k=j;
      }
    }
    fitFunction->FixParameter(0,min);
    baseline=min;
    errbaseline=h->GetBinError(k);      
  }
  if(fixBase<0){
    Int_t npointsAv=TMath::Abs(fixBase);
    Int_t *ind=new Int_t[h->GetNbinsX()];
    Double_t *hval=new Double_t[h->GetNbinsX()];// needed because problems were found with usage of fHist->GetArray();
    for(Int_t k=1;k<=h->GetNbinsX();k++){
      hval[k-1]=h->GetBinContent(k);
    }
    //    Double_t *hval=h->GetArray();
    Double_t errAv=0.,av=0.;
    TMath::Sort(h->GetNbinsX(),hval,ind,kFALSE);// need to exclude under and over flow bins
    // Average of abs(fixbase) lower points
    for(Int_t k=0;k<npointsAv;k++){
      //      Printf("Point %d, index %d,value: %f",k,ind[k],h->GetBinContent(ind[k]+1));
      av+=(h->GetBinContent(ind[k]+1)/(h->GetBinError(ind[k]+1)*h->GetBinError(ind[k]+1)));
      //printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
      errAv+=1./(h->GetBinError(ind[k]+1)*h->GetBinError(ind[k]+1));	  
    }
    av/=errAv;
    errAv=TMath::Sqrt(1./errAv);
    printf("Average baseline: %f +- %f \n",av,errAv);
    fitFunction->FixParameter(0,av);      
    baseline=av;
    errbaseline=errAv;
  }
  if(fixBase==2){// ZYAM, USE 2 POINTS AT +- pi/2
    Double_t errAv=0.,av=0.;
    // Average of abs(fixbase) lower points
    Int_t binPhi=h->FindBin(TMath::Pi()/2.);
    av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
    //	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
    errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  
    if(!isReflected){
      binPhi=h->FindBin(-TMath::Pi()/2.);
      if(binPhi<1)binPhi=1;
      av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
      //	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
      errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  
    }

    av/=errAv;
    errAv=TMath::Sqrt(1./errAv);
    printf("Average baseline: %f +- %f \n",av,errAv);
    fitFunction->FixParameter(0,av);      
    baseline=av;
    errbaseline=errAv;

  }

  if(fixBase==4){// ZYAM, USE 4 points around +- pi/2
    Double_t errAv=0.,av=0.;
    // Average of abs(fixbase) lower points
    Int_t binPhi=h->FindBin(TMath::Pi()/2.);
    av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
    //	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
    errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  

    Double_t binCentreCloseL=h->GetBinCenter(binPhi-1);
    Double_t binCentreCloseR=h->GetBinCenter(binPhi+1);

    if(TMath::Abs(binCentreCloseR-TMath::Pi()/2.)<TMath::Abs(binCentreCloseL-TMath::Pi()/2.))binPhi++;
    else  binPhi--;
    av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
    //	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
    errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  
    if!(isReflected){
	binPhi=h->FindBin(-TMath::Pi()/2.);
	if(binPhi<1)binPhi=h->GetNbinsX();
	av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
	//	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
	errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  

	Int_t binphiL,binphiR;
	if(binPhi==1){
	  binphiL=h->GetNbinsX();
	  binCentreCloseL=h->GetBinCenter(h->GetNbinsX())-TMath::Pi()*2.;
	}
	else {
	  binCentreCloseL=h->GetBinCenter(binPhi-1);
	  binphiL=binPhi-1;
	}
	
	if(binPhi==h->GetNbinsX()){
	  binphiR=1;
	  binCentreCloseR=h->GetBinCenter(1);
	}
	else {
	  binphiR=binPhi+1;
	  binCentreCloseR=h->GetBinCenter(binPhi+1);
	}
	
	if(TMath::Abs(binCentreCloseR+TMath::Pi()/2.)<TMath::Abs(binCentreCloseL+TMath::Pi()/2.)){
	  binPhi=binphiR;
	}
	else {
	  binPhi=binphiL;      
	}
	av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
	//	printf("havl: %f, hist :%f+-%f \n",hval[ind[k]+1],h->GetBinContent(ind[k]+1),h->GetBinError(ind[k]+1));
	errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));	  
      }
    
    av/=errAv;
    errAv=TMath::Sqrt(1./errAv);
    printf("Average baseline: %f +- %f \n",av,errAv);
    fitFunction->FixParameter(0,av);      
    baseline=av;
    errbaseline=errAv;
  }
    // begin change
    if(fixBase==5){// ZYAM, USE 4 points around +- pi/2
        Double_t errAv=0.,av=0.;
        
        
        Int_t nbins = 0;
        Double_t sum = 0;
        for(Int_t binPhi=1; binPhi<=h->GetNbinsX();binPhi++){
            
            if(h->GetBinLowEdge(binPhi)>=-0.5*TMath::Pi() && h->GetBinLowEdge(binPhi+1)<=-0.25*TMath::Pi()){
                cout << "iBin = " << binPhi << endl;
                av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
                errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));
            }
            
            if(h->GetBinLowEdge(binPhi)>=0.25*TMath::Pi() && h->GetBinLowEdge(binPhi+1)<=0.5*TMath::Pi()){
            cout << "iBin = " << binPhi << endl;
            av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
            errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));
            }
        }
        
        av/=errAv;
        errAv=TMath::Sqrt(1./errAv);
      //  av/=2;
      //  errAv/=2;
        printf("Average baseline: %f +- %f \n",av,errAv);
        fitFunction->FixParameter(0,av);
        baseline=av;
        errbaseline=errAv;
        
    } // end change

    
  if(fixMean==1||fixMean==3){
    fitFunction->FixParameter(2,0.);
  }
  if(fixMean==2||fixMean==3){
    if(fitFunc>0)fitFunction->FixParameter(5,TMath::Pi());
  }
  if(isReflected==0){
    h->Fit(fitFunction,"REMI","",-0.5*TMath::Pi(),1.5*TMath::Pi());
  }
  else {
    h->Fit(fitFunction,"REMI","",0*TMath::Pi(),TMath::Pi());
  }

  // Bin counting above baseline
  Double_t base=fitFunction->GetParameter(0);
  Double_t ebase=fitFunction->GetParError(0);
  nsybc=0.;
  ensybc=0.;
  asybc=0.;
  easybc=0.;

  Int_t binMinNS=h->FindBin(-1.5);// slightly more than -pi/2
  if(binMinNS<1)binMinNS=1;
  Int_t binMaxNS=h->FindBin(1.5);// slightly more than -pi/2

  Int_t binMinAS=h->FindBin(3.14-1.5);// slightly more than -pi/2
  Int_t binMaxAS=h->FindBin(3.14+1.5);// slightly more than -pi/2
  if(binMaxAS>h->GetNbinsX())binMaxNS=h->GetNbinsX();
  
  for(Int_t bm=binMinNS;bm<binMaxNS;bm++){
    nsybc+=(h->GetBinContent(bm)-base)*h->GetBinWidth(bm);
  }
  ensybc=TMath::Sqrt(nsybc+(binMaxNS-binMinNS+1)*ebase*ebase)*h->GetBinWidth(bm);

  for(Int_t bm=binMinAS;bm<binMaxAS;bm++){
    asybc+=(h->GetBinContent(bm)-base)*h->GetBinWidth(bm);
  }
  easybc=TMath::Sqrt(asybc+(binMaxAS-binMinAS+1)*ebase*ebase)*h->GetBinWidth(bm);

  printf("Bin counting results: NS y= %f =- %f \n AS y: %f =- %f \n",nsybc,ensybc,asybc,easybc);


  //  if(fitFunc==0)return fitFunction;

  Double_t *par;
  if(fitFunc==1||fitFunc==0){
    par=new Double_t*[7];
  }
  else if(fitFunc==2){
    par=new Double_t*[9];
  }

  fitFunction->GetParameters(par);
    
  TF1 *fGausNS;
  if(fitFunc==0)fGausNS=GetGaussian();
  else fGausNS=GetGaussianPeriodicity();
  fGausNS->SetName("fGausNS");
  if(fitFunc==1||fitFunc==0)fGausNS->SetParameter(0,par[1]);
  else if(fitFunc==2)fGausNS->SetParameter(0,par[1]*par[7]);
  fGausNS->SetParameter(1,par[2]);
  fGausNS->SetParameter(2,par[3]);
  fGausNS->SetLineStyle(2);
  fGausNS->SetLineColor(kBlue);

  TF1 *fGausNS2=0x0;
  if(fitFunc==2){
    TF1 *fGausNS2=GetGaussianPeriodicity();
    fGausNS2->SetName("fGausNS2");
    fGausNS2->SetParameter(0,par[1]*(1.-par[7]));
    fGausNS2->SetParameter(1,par[2]);
    fGausNS2->SetParameter(2,par[8]);
    fGausNS2->SetLineStyle(2);
    fGausNS2->SetLineColor(kCyan);
    
  }

  
  TF1 *fGausAS;
  if(fitFunc==0)fGausAS=GetGaussian();
  else fGausAS=GetGaussianPeriodicity();
  fGausAS->SetName("fGausAS");
  fGausAS->SetParameter(0,par[4]);
  fGausAS->SetParameter(1,par[5]);
  fGausAS->SetParameter(2,par[6]);
  fGausAS->SetLineStyle(2);
  fGausAS->SetLineColor(kGreen);
  
  TF1 *fPed=new TF1("fPed","[0]",-1.7,5.2);  
  fPed->SetParameter(0,par[0]);
  fPed->SetLineColor(6);//pink
  fPed->SetLineStyle(2);
  
  fGausAS->Draw("same");
  fGausNS->Draw("same");
  fPed->Draw("same");
  if(fGausNS2)fGausNS2->Draw("same");

  // TH1D *h2=new TH1D("h","h",500,-TMath::Pi()/2.,TMath::Pi()*3./2.);
  //   h2->Draw();
  //   fitFunction->Draw("same");


  return fitFunction;
}
