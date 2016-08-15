/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: $ */

/////////////////////////////////////////////////////////////
//
// Class for subtracting D-hadron correlations from secondary D from B meson decay from D-hadron correlations of prompt+secondary D mesons
//       n.b. requires the evaluation of the fraction of prompt D mesons using the same receipes used for D meson cross-section and RAA in D2H pag
//       (-> fprompt obtained from FONLL predictions of prompt and secondary D, prompt and secondary D meson efficiencies and a range of hypotheses for
//       RAA(DfromB)/RAA(promptD) (if needed)
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>

#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include "AliHFCorrelationFDsubtraction.h"


using std::cout;
using std::endl;

ClassImp(AliHFCorrelationFDsubtraction)

//************** constructors
AliHFCorrelationFDsubtraction::AliHFCorrelationFDsubtraction() : TNamed(),
  fptmin(),
  fptmax(),
  fhDataUncorrected(),
  fmethod(1),
  fgrConservativeFc(0x0),
  fgrConservativeNb(0x0),
  fnTemplates(0),
  fLastTemplAdded(0),
  fhTemplates(0x0),
  fhDataCorrected(0x0),
  fhRatioSameTemplate(0x0),
  fhRatioAllToCentralTempl(0x0),
  fcCompare(0x0),
  fcRatio(0x0),
  fCountTempl(0), 
  fcAllRatio(0),
  fhEnvelopeMax(0x0),
  fhEnvelopeMin(0x0),
  fhEnvelopeMaxRatio(0x0),
  fhEnvelopeMinRatio(0x0),
  fgrEnvelope(0x0),
  fgrEnvelopeRatio(0x0),
  fSystUseRMSfromFlatDistr(1)
{ // DEFAULT CONSTRUCTOR
  
}


AliHFCorrelationFDsubtraction::AliHFCorrelationFDsubtraction(const char* name, const char* title) : 
  TNamed(name,title),
  fptmin(),
  fptmax(),
  fhDataUncorrected(),
  fmethod(1),
  fgrConservativeFc(0x0),
  fgrConservativeNb(0x0),
  fnTemplates(0),
  fLastTemplAdded(0),
  fhTemplates(0x0),
  fhDataCorrected(0x0),
  fhRatioSameTemplate(0x0),
  fhRatioAllToCentralTempl(0x0),
  fcCompare(0x0),
  fcRatio(0x0),
  fCountTempl(0), 
  fcAllRatio(0),
  fhEnvelopeMax(0x0),
  fhEnvelopeMin(0x0),
  fhEnvelopeMaxRatio(0x0),
  fhEnvelopeMinRatio(0x0),
  fgrEnvelope(0x0),
  fgrEnvelopeRatio(0x0),
  fSystUseRMSfromFlatDistr(1)
{// default constructor
  
}

AliHFCorrelationFDsubtraction::~AliHFCorrelationFDsubtraction(){


  delete fhDataUncorrected;
  
  
  delete fgrConservativeFc;
  delete fgrConservativeNb;
  
  for(Int_t j=0;j<fnTemplates;j++){
    delete fhTemplates[j];
  }
  delete [] fhTemplates;
  
  for(Int_t j=0;j<3*fnTemplates;j++){
    delete fhDataCorrected[j];
  }
  delete [] fhDataCorrected;
  
  for(Int_t j=0;j<2*fnTemplates;j++){
    delete fhRatioSameTemplate[j];
  }
  delete [] fhRatioSameTemplate;
  

 for(Int_t j=0;j<3*fnTemplates-1;j++){
    delete fhRatioAllToCentralTempl[j];
  }
  delete [] fhRatioAllToCentralTempl;

  for(Int_t j=0;j<fnTemplates;j++){
    delete fcCompare[j];
  }
  delete [] fcCompare;

  for(Int_t j=0;j<fnTemplates;j++){
    delete fcRatio[j];
  }
  delete [] fcRatio;

  delete fcAllRatio;


  delete fhEnvelopeMax;
  delete fhEnvelopeMin;
  delete fhEnvelopeMaxRatio;
  delete fhEnvelopeMinRatio;
  delete fgrEnvelope;
  delete fgrEnvelopeRatio;
}


Bool_t AliHFCorrelationFDsubtraction::Init(){// init histo, graphs and canvases
  
  if(fnTemplates<=0){
    Printf("Number of templates not defined, aborting");
    return kFALSE;
  }
  fhTemplates=new TH1D*[fnTemplates];
  fhDataCorrected=new TH1D*[fnTemplates*3];
  fhRatioSameTemplate=new TH1D*[fnTemplates*2];
  fhRatioAllToCentralTempl=new TH1D*[fnTemplates*3-1];
  fcCompare=new TCanvas*[fnTemplates];
  fcRatio=new TCanvas*[fnTemplates];
  
  fCountTempl=0;  
  return kTRUE;
}
Bool_t AliHFCorrelationFDsubtraction::AddTemplateHisto(TH1D *h){
  if(fnTemplates<1){
    Printf("You first have to specify how many templates you want to use, aborting");
    return kFALSE;
  }
  if(fLastTemplAdded>=fnTemplates){
    Printf("You want to add more templates than those declared, aborting");
    return kFALSE;
  }
  fhTemplates[fLastTemplAdded]=(TH1D*)h->Clone();
  fLastTemplAdded++;
  return kTRUE;
}

void AliHFCorrelationFDsubtraction::SubtractFeedDown(TH1D *hFDtempl){
  TGraphAsymmErrors *gr;
  
  Double_t *xgr,*elxgr,*ehxgr,*fprompt,*elfprompt,*ehfprompt;
  if(fmethod==1||fmethod==3){
    gr=fgrConservativeFc;
  }
  if(fmethod==2){
    gr=fgrConservativeNb;
  }
  if(!gr){
    Printf("TGraph for feed-down subtraction not set. Aborting");
    return;
  }    
  
  xgr=gr->GetX();
  elxgr=gr->GetEXlow();
  ehxgr=gr->GetEXhigh();
  fprompt=gr->GetY();
  elfprompt=gr->GetEYlow();
  ehfprompt=gr->GetEYhigh();
  
  Int_t binmin=-1,binmax=-1,bincenter=-1;
  Double_t ptcentr=0.5*(fptmin+fptmax);
  Double_t fpromptmin=1.,fpromptmax=0.,fpromptcentr;
  // Look for right bins
  for(Int_t j=0;j<gr->GetN();j++){
    if(binmin<0&&fptmin*1.00001>(xgr[j]-elxgr[j])&&fptmin<(xgr[j]+ehxgr[j]))binmin=j;
    if(binmax<0&&fptmax>(xgr[j]-elxgr[j])&&fptmax*0.9999<(xgr[j]+ehxgr[j]))binmax=j;      
    if(ptcentr*1.00001>(xgr[j]-elxgr[j])&&ptcentr<(xgr[j]+ehxgr[j]))bincenter=j;
  }
  
  if(TMath::Abs(xgr[binmin]-elxgr[binmin]-fptmin)>0.0001){
    Printf("Bin mismatch: xmin=%f; aborting ",xgr[binmin]-elxgr[binmin]);
    return;
  }
  
  // central value for bin at centre of the bin (temporary)
  fpromptcentr=fprompt[bincenter];
  

  // Conservative approach: look for min and max values inside the pt range
  for(Int_t j=binmin;j<=binmax;j++){
    if(fprompt[j]+ehfprompt[j]>fpromptmax)fpromptmax=fprompt[j]+ehfprompt[j];
    if(fprompt[j]-elfprompt[j]<fpromptmin)fpromptmin=fprompt[j]-elfprompt[j];
  }
  
  
  TH1D *hFDCentr=(TH1D*)hFDtempl->Clone("hFDtemplCentr");
  hFDCentr->Scale(1.-fpromptcentr);
  printf("Centr fsec=%f \n",1.-fpromptcentr);

  fhDataCorrected[fCountTempl*3]=(TH1D*)fhDataUncorrected->Clone(Form("hDataCorrectedTempl%dCentrFprompt",fCountTempl));
  fhDataCorrected[fCountTempl*3]->Add(hFDCentr,-1.);// Treatment of Stat unc.
  fhDataCorrected[fCountTempl*3]->Scale(1./fpromptcentr);
  fhDataCorrected[fCountTempl*3]->SetTitle(Form("Templ %d central f_{prompt}",fCountTempl));

  TH1D *hFDmax=(TH1D*)hFDtempl->Clone("hFDtemplMax");
  hFDmax->Scale(1.-fpromptmin);
  Printf("Max fsec=%f ",1.-fpromptmin);
  fhDataCorrected[fCountTempl*3+1]=(TH1D*)fhDataUncorrected->Clone(Form("hDataCorrectedTempl%dMaxFprompt",fCountTempl)); 
  fhDataCorrected[fCountTempl*3+1]->Add(hFDmax,-1.);// Treatment of Stat unc.
  fhDataCorrected[fCountTempl*3+1]->Scale(1./fpromptmin);
  fhDataCorrected[fCountTempl*3+1]->SetTitle(Form("Templ %d min f_{prompt}",fCountTempl));


  TH1D *hFDmin=(TH1D*)hFDtempl->Clone("hFDtemplMin");
  hFDmin->Scale(1.-fpromptmax);
  Printf("Min fsec=%f ",1.-fpromptmax);
  fhDataCorrected[fCountTempl*3+2]=(TH1D*)fhDataUncorrected->Clone(Form("hDataCorrectedTempl%dMinFprompt",fCountTempl)); 
  fhDataCorrected[fCountTempl*3+2]->Add(hFDmin,-1.);// Treatment of Stat unc.
  fhDataCorrected[fCountTempl*3+2]->Scale(1./fpromptmax);
  fhDataCorrected[fCountTempl*3+2]->SetTitle(Form("Templ %d max f_{prompt}",fCountTempl));

 
  fcCompare[fCountTempl]=new TCanvas(Form("cCompareTempl%d",fCountTempl),Form("cCompare%d",fCountTempl),700,700);
  fcCompare[fCountTempl]->cd();

  TH1D *hData=(TH1D*)fhDataUncorrected->DrawCopy();
  hData->SetTitle("Uncorrected");
  hData->SetYTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#varphi} (rad^{-1})");
  hData->SetXTitle("#varphi_{ass}-#varphi_{trig} (rad)");
  hData->SetMarkerStyle(25);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerSize(1);
  hData->SetLineColor(kBlack);

  fhDataCorrected[fCountTempl*3]->SetLineColor(kRed);
  fhDataCorrected[fCountTempl*3]->SetYTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#varphi} (rad^{-1})");
  fhDataCorrected[fCountTempl*3]->SetXTitle("#varphi_{ass}-#varphi_{trig} (rad)");
  fhDataCorrected[fCountTempl*3]->SetLineWidth(2);
  fhDataCorrected[fCountTempl*3]->SetMarkerColor(kRed);
  fhDataCorrected[fCountTempl*3]->SetMarkerStyle(21);
  fhDataCorrected[fCountTempl*3]->SetMarkerSize(1);
  fhDataCorrected[fCountTempl*3]->Draw("sames");

  fhDataCorrected[fCountTempl*3+1]->SetLineColor(kBlue);
  fhDataCorrected[fCountTempl*3+1]->SetMarkerColor(kBlue);
  fhDataCorrected[fCountTempl*3+1]->SetMarkerStyle(23);
  fhDataCorrected[fCountTempl*3+1]->SetMarkerSize(1);
  fhDataCorrected[fCountTempl*3+1]->Draw("sames");

  fhDataCorrected[fCountTempl*3+2]->SetLineColor(kGreen-6);
  fhDataCorrected[fCountTempl*3+2]->SetMarkerColor(kGreen-6);
  fhDataCorrected[fCountTempl*3+2]->SetMarkerStyle(22);
  fhDataCorrected[fCountTempl*3+2]->SetMarkerSize(1);
  fhDataCorrected[fCountTempl*3+2]->Draw("sames");
  
  fcRatio[fCountTempl]=new TCanvas(Form("fcRatioTempl%d",fCountTempl),Form("fcRatioTempl%d",fCountTempl),700,700);
  fcRatio[fCountTempl]->cd();
  fhRatioSameTemplate[fCountTempl*2]=(TH1D*)fhDataCorrected[fCountTempl*3+2]->Clone("hRatioMinFSecToCentr");
  fhRatioSameTemplate[fCountTempl*2]->Divide(fhDataCorrected[fCountTempl*3]);
  for(Int_t j=1;j<=fhRatioSameTemplate[fCountTempl*2]->GetNbinsX();j++){
    fhRatioSameTemplate[fCountTempl*2]->SetBinError(j,0.001);
  }
  
  fhRatioSameTemplate[fCountTempl*2+1]=(TH1D*)fhDataCorrected[fCountTempl*3+1]->Clone("hRatioMaxFSecToCentr");
  fhRatioSameTemplate[fCountTempl*2+1]->Divide(fhDataCorrected[fCountTempl*3]);
  for(Int_t j=1;j<=fhRatioSameTemplate[fCountTempl*2+1]->GetNbinsX();j++){
    fhRatioSameTemplate[fCountTempl*2+1]->SetBinError(j,0.001);
  }

  fhRatioSameTemplate[fCountTempl*2]->Draw();
  fhRatioSameTemplate[fCountTempl*2+1]->Draw("same");

  fCountTempl++;  
  delete hFDCentr;
  delete hFDmax;
  delete hFDmin;

}

TGraphAsymmErrors* AliHFCorrelationFDsubtraction::GetUncGraphFromHistos(TH1D *hRef,TH1D *hMin,TH1D *hMax){
  
  //  Int_t npoints=hMin->GetNbinsX();
  Double_t ew=hMin->GetBinWidth(1)/2.;
  Double_t value,eyl,eym;
  
  TGraphAsymmErrors *gr=new TGraphAsymmErrors();
  for(Int_t j=1;j<=hMin->GetNbinsX();j++){
    if(hRef){
      value=hRef->GetBinContent(j);
      eyl=hMin->GetBinContent(j)-value;
      if(eyl<0.)eyl*=-1.;
      if(hMax){
	eym=hMax->GetBinContent(j)-value;
	if(eym<0.)eym*=-1.;
      }
      else eym=eyl;
    }
    else {
      value=0.;
      eyl=hMin->GetBinContent(j)-1;
      if(eyl<0.)eyl*=-1.;
      if(hMax){
	eym=hMax->GetBinContent(j)-1;
	if(eym<0.)eym*=-1.;
      }
      else eym=eyl;
    }
    
    gr->SetPoint(j-1,hMin->GetBinCenter(j),value);
    gr->SetPointError(j-1,ew,ew,eyl,eym);
  }
  
  return gr;
}


TH1D* AliHFCorrelationFDsubtraction::ReflectHisto(TH1D *h,Double_t scale){
  
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
    h2->SetBinContent(j2,scale*(y+y0));
    h2->SetBinError(j2,scale*TMath::Sqrt(ey0*ey0+ey*ey));
  }
  
  return h2;
}

TH1D* AliHFCorrelationFDsubtraction::DuplicateHistoTo2piRange(TH1D *h,Double_t scale){
  if(h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(h->GetNbinsX())-h->GetBinLowEdge(1)>1.01*TMath::Pi()){
    Printf("AliHFCorrelationFDsubtraction: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }

  if(h->GetBinLowEdge(h->GetNbinsX())>1.01*TMath::Pi()){
    Printf("AliHFCorrelationFDsubtraction: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }
  if(h->GetBinLowEdge(1)<-0.01){
    Printf("AliHFCorrelationFDsubtraction: DuplicateHistoTo2PiRange works for histos with x range between 0 and pi");
    return 0x0;
  }


  TH1D *h2=new TH1D(Form("%sTwoPiRange",h->GetName()),Form("%sTwoPiRange",h->GetName()),h->GetNbinsX()*2.,-TMath::Pi()/2.,1.5*TMath::Pi());
  for(Int_t j=1;j<=h2->GetNbinsX();j++){
    Double_t x=h2->GetBinCenter(j);
    Int_t j2;
    if(x>0&&x<TMath::Pi())j2=h->FindBin(x);
    else if(x<0)j2=h->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h->FindBin(2.*TMath::Pi()-x);
    else {
      printf("Point %d excluded \n",j);
      continue;
    }
    Double_t y=h->GetBinContent(j2);
    Double_t ey=h->GetBinError(j2);
    h2->SetBinContent(j,y*scale);
    h2->SetBinError(j,ey*scale);
  }
  
  return h2;
}


TH1D* AliHFCorrelationFDsubtraction::GetHistoRelSystUncMin(){// Method to extrac the rel syst unc to be used later on in the analysis from the envelopes. IMPORTANT: remember
  // that with our procedure the rel. unc. is variation/distribution, and thus its value is strongly dependent on the baseline: only the absolute unc (->variation) is really meaningful
  // Options: if fSystUseRMSfromFlatDistr
  // =0 the full spread (->DeltaMax, DeltaMin) of the envelopes is used to define the syst unc, 
  // =1 (default) DeltaMax/sqrt(3) and DeltaMin/sqrt(3) is used; 
  // =2: as 1 but then the result is reflected to the range (0,pi), to decrease stat unc; 
  // =3 as 1 then a smoothing is done (pol0 fit over 3 points), then the result is refelcted; other values not implemented yet (e.g. =2 could be taking the RMS)
  // =4 as 1 then a smoothing is done (linear fit over 3 points), then the result is refelcted; other values not implemented yet (e.g. =2 could be taking the RMS)
  // =13 as 3 + symmetrization
  // =14 as 4 + symmetrization

  if(!fhEnvelopeMinRatio){
    Printf("Cannot produce Rel syst unc histo: you first need to calculate the envelope");
    return 0x0;
  }
  TH1D *h=(TH1D*)fhEnvelopeMinRatio->Clone("hRelSystFeedDownMin");

  if(fSystUseRMSfromFlatDistr==0){
    for(Int_t j=1;j<=h->GetNbinsX();j++){
      h->SetBinContent(j,fhEnvelopeMinRatio->GetBinContent(j)-1.);
    }
    return h;
  }

  if(fSystUseRMSfromFlatDistr==1||fSystUseRMSfromFlatDistr==2||fSystUseRMSfromFlatDistr==3||fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==13||fSystUseRMSfromFlatDistr==14){
    Double_t sqrtThree=TMath::Sqrt(3.);
    for(Int_t j=1;j<=h->GetNbinsX();j++){
      h->SetBinContent(j,fhEnvelopeMinRatio->GetBinContent(j)/sqrtThree-1./sqrtThree);
    }
    if(fSystUseRMSfromFlatDistr==1)return h;
  
    if(fSystUseRMSfromFlatDistr==2){
      TH1D *hRefl=ReflectHisto(h);
      TH1D *hBackTo2Pi=DuplicateHistoTo2piRange(hRefl); // In this way we symmetrized to 0-pi    
      delete h;
      delete hRefl;      
      hBackTo2Pi->SetName("hRelSystFeedDownMin");
      return hBackTo2Pi;
    }

    // case fSystUseRMSfromFlatDistr==3 or 4
    // Do smoothing and then symmetrize
    // For the smoothing, assume linear trend every 3 points, starting from transverse region
    Int_t jstart=h->FindBin(TMath::Pi()*0.5)-1;
    TF1 *f;
    if(fSystUseRMSfromFlatDistr==3||fSystUseRMSfromFlatDistr==13){
      f=new TF1("mypol0","[0]",-TMath::Pi()*0.5,TMath::Pi()*1.5);
    }
    else if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14){
      f=new TF1("mypol1","[0]+x*[1]",-TMath::Pi()*0.5,TMath::Pi()*1.5);
    }

    TH1D *hOut=(TH1D*)h->Clone("hOut");

    Double_t median=0;
    Double_t rms=0,mean=0;
    Double_t *arr=hOut->GetArray();// note that this includes under and overflow
    Int_t *index=new Int_t[hOut->GetNbinsX()];
    TMath::Sort(hOut->GetNbinsX(),&arr[1],index,kTRUE);
    median=arr[index[(UInt_t)(hOut->GetNbinsX()/2)+2]];// the +2 comes from: +1 because the indices start from 0 but w.r.t. to the reduced array arr (w/o first bin); to this we add +1 because we prefer to take the point after the median (also to accounts for possible odd number of bins)
    // fast calculation of rms; mean computed excluding tails
    for(Int_t k=0;k<hOut->GetNbinsX()-4;k++){// n-2 points for cutting tails
      rms+=arr[index[k]+1]*arr[index[k]+1];// +1 due to underflow indexing
      mean+=arr[index[k]+1];
    }
    mean/=(hOut->GetNbinsX()-4.);
    rms=TMath::Sqrt(rms-mean*mean);


    for(Int_t j=jstart;j>=1;j--){// going towards 0

      if(j>=3){
	Double_t xmax=h->GetBinLowEdge(j+1)+h->GetBinWidth(j+1);
	Double_t xmin=h->GetBinLowEdge(j-1);

	Double_t y1=h->GetBinContent(j-1);
	Double_t y2=h->GetBinContent(j);
	Double_t y3=h->GetBinContent(j+1);
	if(y3<median-2.5*rms){
	  h->SetBinContent(j+1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",y3,median-1.5*rms);
	}
	else h->SetBinContent(j+1,y3);
	h->SetBinContent(j,y2);
	if(y1<median-2.5*rms){
	  h->SetBinContent(j-1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",y1,median-1.5*rms);
	}
	else h->SetBinContent(j-1,y1);

	f->SetParameter(0,h->GetBinContent(j));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(j+1)-h->GetBinContent(j-1))/(h->GetBinCenter(j+1)-h->GetBinCenter(j-1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(hOut->GetBinCenter(j)));
      }
      else if(j==2){
	Double_t xmax=h->GetBinLowEdge(3)+h->GetBinWidth(3);
	Double_t xmin=h->GetBinLowEdge(1);
	Double_t y1=h->GetBinContent(1);
	Double_t y2=h->GetBinContent(2);
	Double_t y3=h->GetBinContent(3);

	h->SetBinContent(3,y2);
	if(y1<median-2.5*rms){
	  h->SetBinContent(2,median-1.5*rms);
	  Printf("Moving to median: %f to %f",y1,median-1.5*rms);
	}
	else h->SetBinContent(2,y1);
	if(h->GetBinContent(h->GetNbinsX())<median-2.5*rms){
	  h->SetBinContent(1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(h->GetNbinsX()),median-1.5*rms);
	}
	else h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()));

	//	h->SetBinContent(3,y2);
	//	h->SetBinContent(2,y1);
	//	h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()));
	f->SetParameter(0,h->GetBinContent(2));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(3)-h->GetBinContent(1))/(h->GetBinCenter(3)-h->GetBinCenter(1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(2)));
	h->SetBinContent(3,y3);
	h->SetBinContent(2,y2);
	h->SetBinContent(1,y1);
      }
      else if(j==1){
	Double_t xmax=h->GetBinLowEdge(3)+h->GetBinWidth(3);
	Double_t xmin=h->GetBinLowEdge(1);
	Double_t y1=h->GetBinContent(1);
	Double_t y2=h->GetBinContent(2);
	Double_t y3=h->GetBinContent(3);

	h->SetBinContent(3,y1);
	if(h->GetBinContent(h->GetNbinsX())<median-2.5*rms){
	  h->SetBinContent(2,median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(h->GetNbinsX()),median-1.5*rms);
	}
	else h->SetBinContent(2,h->GetBinContent(h->GetNbinsX()));
	if(h->GetBinContent(h->GetNbinsX()-1)<median-2.5*rms){
	  h->SetBinContent(1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(h->GetNbinsX()-1),median-1.5*rms);
	}
	else h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()-1));

	//	h->SetBinContent(3,y1);
	//	h->SetBinContent(2,h->GetBinContent(h->GetNbinsX()));
	//	h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()-1));
	f->SetParameter(0,h->GetBinContent(2));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(3)-h->GetBinContent(1))/(h->GetBinCenter(3)-h->GetBinCenter(1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(2)));
	h->SetBinContent(3,y3);
	h->SetBinContent(2,y2);
	h->SetBinContent(1,y1);
      }

    }
    for(Int_t j=jstart+1;j<=hOut->GetNbinsX();j++){// now going towards 2pi
      
      if(j<=hOut->GetNbinsX()-2){
	Double_t xmax=h->GetBinLowEdge(j+1)+h->GetBinWidth(j+1);
	Double_t xmin=h->GetBinLowEdge(j-1);

	Double_t y1=h->GetBinContent(j-1);
	Double_t y2=h->GetBinContent(j);
	Double_t y3=h->GetBinContent(j+1);
	if(y3<median-2.5*rms){
	  h->SetBinContent(j+1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",y3,median-1.5*rms);
	}
	else h->SetBinContent(j+1,y3);
	h->SetBinContent(j,y2);
	if(y1<median-2.5*rms){
	  h->SetBinContent(j-1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",y1,median-1.5*rms);
	}
	else h->SetBinContent(j-1,y1);

	f->SetParameter(0,h->GetBinContent(j));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(j+1)-h->GetBinContent(j-1))/(h->GetBinCenter(j+1)-h->GetBinCenter(j-1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(hOut->GetBinCenter(j)));
      }
      else if(j==hOut->GetNbinsX()-1){
	Double_t xmax=h->GetBinLowEdge(hOut->GetNbinsX())+h->GetBinWidth(hOut->GetNbinsX());
	Double_t xmin=h->GetBinLowEdge(hOut->GetNbinsX()-2);
	Double_t y1=h->GetBinContent(hOut->GetNbinsX()-2);
	Double_t y2=h->GetBinContent(hOut->GetNbinsX()-1);
	Double_t y3=h->GetBinContent(hOut->GetNbinsX());


	h->SetBinContent(hOut->GetNbinsX()-2,y2);
	if(h->GetBinContent(hOut->GetNbinsX())<median-2.5*rms){
	  h->SetBinContent(hOut->GetNbinsX()-1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(hOut->GetNbinsX()),median-1.5*rms);
	}
	else h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(h->GetNbinsX()));
	if(h->GetBinContent(1)<median-2.5*rms){
	  Printf("Moving to median: %f to %f",h->GetBinContent(1),median-1.5*rms);
	  h->SetBinContent(hOut->GetNbinsX(),median-1.5*rms);
	}
	else h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(1));


	//	h->SetBinContent(hOut->GetNbinsX()-2,y2);
	//	h->SetBinContent(hOut->GetNbinsX()-1,y3);
	//	h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(1));
	f->SetParameter(0,h->GetBinContent(hOut->GetNbinsX()-1));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(hOut->GetNbinsX())-h->GetBinContent(hOut->GetNbinsX()-2))/(h->GetBinCenter(hOut->GetNbinsX())-h->GetBinCenter(hOut->GetNbinsX()-2)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(hOut->GetNbinsX()-1)));
	h->SetBinContent(hOut->GetNbinsX()-2,y1);
	h->SetBinContent(hOut->GetNbinsX()-1,y2);
	h->SetBinContent(hOut->GetNbinsX(),y3);
      }
      else if(j==hOut->GetNbinsX()){
	Double_t xmax=h->GetBinLowEdge(hOut->GetNbinsX())+h->GetBinWidth(hOut->GetNbinsX());
	Double_t xmin=h->GetBinLowEdge(hOut->GetNbinsX()-2);
	Double_t y1=h->GetBinContent(hOut->GetNbinsX()-2);
	Double_t y2=h->GetBinContent(hOut->GetNbinsX()-1);
	Double_t y3=h->GetBinContent(hOut->GetNbinsX());

	h->SetBinContent(hOut->GetNbinsX()-2,y3);
	if(h->GetBinContent(1)<median-2.5*rms){
	  h->SetBinContent(hOut->GetNbinsX()-1,median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(1),median-1.5*rms);
	}
	else h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(1));
	if(h->GetBinContent(2)<median-2.5*rms){
	  h->SetBinContent(hOut->GetNbinsX(),median-1.5*rms);
	  Printf("Moving to median: %f to %f",h->GetBinContent(2),median-1.5*rms);
	}
	else h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(2));

	//	h->SetBinContent(hOut->GetNbinsX()-2,y3);
	//	h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(1));
	//	h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(2));
	f->SetParameter(0,h->GetBinContent(hOut->GetNbinsX()-1));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(hOut->GetNbinsX())-h->GetBinContent(hOut->GetNbinsX()-2))/(h->GetBinCenter(hOut->GetNbinsX())-h->GetBinCenter(hOut->GetNbinsX()-2)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(hOut->GetNbinsX()-1)));
	h->SetBinContent(hOut->GetNbinsX()-2,y1);
	h->SetBinContent(hOut->GetNbinsX()-1,y2);
	h->SetBinContent(hOut->GetNbinsX(),y3);
      }
      
    }

    // Smoothing done, now simmetrize if requested
    if(fSystUseRMSfromFlatDistr==13||fSystUseRMSfromFlatDistr==14){
      TH1D *hRefl=ReflectHisto(hOut);
      TH1D *hBackTo2Pi=DuplicateHistoTo2piRange(hRefl); // In this way we symmetrized to 0-pi    
      delete[] index;
      delete h;
      delete hRefl;      
      delete hOut;
      delete f;
      hBackTo2Pi->SetName("hRelSystFeedDownMin");
      return hBackTo2Pi;
    }
    else{
      delete[] index;
      delete h;     
      delete f;
      hOut->SetName("hRelSystFeedDownMin");
      return hOut;
    }

    
  }
  
  Printf("AliHFCorrelationFDsubtraction: wrong option for getting syst unc");
  return 0x0;
  
  
}

TH1D* AliHFCorrelationFDsubtraction::GetHistoRelSystUncMax(){
  if(!fhEnvelopeMaxRatio){
    Printf("Cannot produce Rel syst unc histo: you first need to calculate the envelope");
    return 0x0;
  }
  TH1D *h=(TH1D*)fhEnvelopeMaxRatio->Clone("hRelSystFeedDownMax");
  if(fSystUseRMSfromFlatDistr==0){
    for(Int_t j=1;j<=h->GetNbinsX();j++){
      h->SetBinContent(j,fhEnvelopeMaxRatio->GetBinContent(j)-1);
    }
    return h;
  }
  
  if(fSystUseRMSfromFlatDistr==1||fSystUseRMSfromFlatDistr==2||fSystUseRMSfromFlatDistr==3||fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==13||fSystUseRMSfromFlatDistr==14){
    Double_t sqrtThree=TMath::Sqrt(3.);
    for(Int_t j=1;j<=h->GetNbinsX();j++){
      h->SetBinContent(j,fhEnvelopeMaxRatio->GetBinContent(j)/sqrtThree-1./sqrtThree);
    }
    if(fSystUseRMSfromFlatDistr==1)return h;
  
    if(fSystUseRMSfromFlatDistr==2){
      TH1D *hRefl=ReflectHisto(h);
      TH1D *hBackTo2Pi=DuplicateHistoTo2piRange(hRefl); // In this way we symmetrized to 0-pi    
      delete h;
      delete hRefl;      
      hBackTo2Pi->SetName("hRelSystFeedDownMin");
      return hBackTo2Pi;
    }
    
    // case fSystUseRMSfromFlatDistr==3 || ==4
    // Do smoothing and then symmetrize
    // For the smoothing, assume linear trend every 3 points, starting from transverse region
    Int_t jstart=h->FindBin(TMath::Pi()*0.5)-1;
    TF1 *f;
    if(fSystUseRMSfromFlatDistr==3||fSystUseRMSfromFlatDistr==13){
      f=new TF1("mypol0","[0]",-TMath::Pi()*0.5,TMath::Pi()*1.5);
    }
    else if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14){
      f=new TF1("mypol1","[0]+x*[1]",-TMath::Pi()*0.5,TMath::Pi()*1.5);
    }
    TH1D *hOut=(TH1D*)h->Clone("hOut");
    // Check median
    Double_t median=0;
    Double_t rms=0,mean=0;
    Double_t *arr=hOut->GetArray();// note that this includes under and overflow
    Int_t *index=new Int_t[hOut->GetNbinsX()];
    TMath::Sort(hOut->GetNbinsX(),&arr[1],index,kFALSE);
    median=arr[index[(UInt_t)(hOut->GetNbinsX()/2)+2]];// the +2 comes from: +1 because the indices start from 0 but w.r.t. to the reduced array arr (w/o first bin); to this we add +1 because we prefer to take the point after the median (also to accounts for possible odd number of bins)
    // fast calculation of rms; mean computed excluding tails
    for(Int_t k=0;k<hOut->GetNbinsX()-4;k++){// n-2 points for cutting tails
      rms+=arr[index[k]+1]*arr[index[k]+1];// +1 due to underflow indexing
      mean+=arr[index[k]+1];
    }
    mean/=(hOut->GetNbinsX()-4.);
    rms=TMath::Sqrt(rms-mean*mean);

    for(Int_t j=jstart;j>=1;j--){// going towards 0

      //      Printf("Value is for point %d : %f, around : %f, %f",j,h->GetBinContent(j),h->GetBinContent(j-1),h->GetBinContent(j+1));
      //      Printf("The median is: %f, the truncate mean is %f, the rms: %f",median,mean,rms);
    
      if(j>=3){

	Double_t xmax=h->GetBinLowEdge(j+1)+h->GetBinWidth(j+1);
	Double_t xmin=h->GetBinLowEdge(j-1);
	
	// now
	Double_t y1=h->GetBinContent(j-1);
	Double_t y2=h->GetBinContent(j);
	Double_t y3=h->GetBinContent(j+1);
	if(y3<median+2.5*rms)h->SetBinContent(j+1,y3);
	else h->SetBinContent(j+1,median+1.5*rms);
	h->SetBinContent(j,y2);
	if(y1<median+2.5*rms)h->SetBinContent(j-1,y1);
	else h->SetBinContent(j-1,median+1.5*rms);
	
	f->SetParameter(0,h->GetBinContent(j));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(j+1)-h->GetBinContent(j-1))/(h->GetBinCenter(j+1)-h->GetBinCenter(j-1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(2)));
	h->SetBinContent(j+1,y3);
	h->SetBinContent(j,y2);
	h->SetBinContent(j-1,y1);
	//	Printf("Value modif  : %f", f->Eval(hOut->GetBinCenter(j)));
	
      }
      else if(j==2){
	Double_t xmax=h->GetBinLowEdge(3)+h->GetBinWidth(3);
	Double_t xmin=h->GetBinLowEdge(1);
	Double_t y1=h->GetBinContent(1);
	Double_t y2=h->GetBinContent(2);
	Double_t y3=h->GetBinContent(3);


	h->SetBinContent(3,y2);
	if(y1<median+2.5*rms)h->SetBinContent(2,y1);
	else h->SetBinContent(2,median+1.5*rms);
	if(h->GetBinContent(h->GetNbinsX())<median+2.5*rms)h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()));
	else h->SetBinContent(1,median+1.5*rms);

	//	h->SetBinContent(3,y2);
	//	h->SetBinContent(2,y1);
	//	h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()));
	f->SetParameter(0,h->GetBinContent(2));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(3)-h->GetBinContent(1))/(h->GetBinCenter(3)-h->GetBinCenter(1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(2)));
	h->SetBinContent(3,y3);
	h->SetBinContent(2,y2);
	h->SetBinContent(1,y1);
      }
      else if(j==1){
	Double_t xmax=h->GetBinLowEdge(3)+h->GetBinWidth(3);
	Double_t xmin=h->GetBinLowEdge(1);
	Double_t y1=h->GetBinContent(1);
	Double_t y2=h->GetBinContent(2);
	Double_t y3=h->GetBinContent(3);
	h->SetBinContent(3,y1);
	if(h->GetBinContent(h->GetNbinsX())<median+2.5*rms)h->SetBinContent(2,h->GetBinContent(h->GetNbinsX()));
	else h->SetBinContent(2,median+1.5*rms);
	if(h->GetBinContent(h->GetNbinsX()-1)<median+2.5*rms)h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()-1));
	else h->SetBinContent(1,median+1.5*rms);
	//	h->SetBinContent(3,y1);
	//	h->SetBinContent(2,h->GetBinContent(h->GetNbinsX()));
	//	h->SetBinContent(1,h->GetBinContent(h->GetNbinsX()-1));
	f->SetParameter(0,h->GetBinContent(2));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)	f->SetParameter(1,(h->GetBinContent(3)-h->GetBinContent(1))/(h->GetBinCenter(3)-h->GetBinCenter(1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(2)));
	h->SetBinContent(3,y3);
	h->SetBinContent(2,y2);
	h->SetBinContent(1,y1);
      }
    }
    for(Int_t j=jstart+1;j<=hOut->GetNbinsX();j++){// now going towards 2pi
      
      if(j<=hOut->GetNbinsX()-2){

	// now
	Double_t xmax=h->GetBinLowEdge(j+1)+h->GetBinWidth(j+1);
	Double_t xmin=h->GetBinLowEdge(j-1);
	Double_t y1=h->GetBinContent(j-1);
	Double_t y2=h->GetBinContent(j);
	Double_t y3=h->GetBinContent(j+1);
	if(y3<median+2.5*rms)h->SetBinContent(j+1,y3);
	else h->SetBinContent(j+1,median+1.5*rms);
	h->SetBinContent(j,y2);
	if(y1<median+2.5*rms)h->SetBinContent(j-1,y1);
	else h->SetBinContent(j-1,median+1.5*rms);
	//	Double_t xmax=h->GetBinLowEdge(j+1)+h->GetBinWidth(j+1);
	//	Double_t xmin=h->GetBinLowEdge(j-1);
	f->SetParameter(0,h->GetBinContent(j));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(j+1)-h->GetBinContent(j-1))/(h->GetBinCenter(j+1)-h->GetBinCenter(j-1)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(hOut->GetBinCenter(j)));
      }
      else if(j==hOut->GetNbinsX()-1){
	Double_t xmax=h->GetBinLowEdge(hOut->GetNbinsX())+h->GetBinWidth(hOut->GetNbinsX());
	Double_t xmin=h->GetBinLowEdge(hOut->GetNbinsX()-2);
	Double_t y1=h->GetBinContent(hOut->GetNbinsX()-2);
	Double_t y2=h->GetBinContent(hOut->GetNbinsX()-1);
	Double_t y3=h->GetBinContent(hOut->GetNbinsX());



	h->SetBinContent(hOut->GetNbinsX()-2,y2);
	if(h->GetBinContent(hOut->GetNbinsX())<median+2.5*rms)h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(h->GetNbinsX()));
	else h->SetBinContent(hOut->GetNbinsX()-1,median+1.5*rms);
	if(h->GetBinContent(1)<median+2.5*rms)h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(1));
	else h->SetBinContent(hOut->GetNbinsX(),median+1.5*rms);

	//	h->SetBinContent(hOut->GetNbinsX()-2,y2);
	//	h->SetBinContent(hOut->GetNbinsX()-1,y3);
	//	h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(1));
	f->SetParameter(0,h->GetBinContent(hOut->GetNbinsX()-1));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)f->SetParameter(1,(h->GetBinContent(hOut->GetNbinsX())-h->GetBinContent(hOut->GetNbinsX()-2))/(h->GetBinCenter(hOut->GetNbinsX())-h->GetBinCenter(hOut->GetNbinsX()-2)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(hOut->GetNbinsX()-1)));
	h->SetBinContent(hOut->GetNbinsX()-2,y1);
	h->SetBinContent(hOut->GetNbinsX()-1,y2);
	h->SetBinContent(hOut->GetNbinsX(),y3);
      }
      else if(j==hOut->GetNbinsX()){
	Double_t xmax=h->GetBinLowEdge(hOut->GetNbinsX())+h->GetBinWidth(hOut->GetNbinsX());
	Double_t xmin=h->GetBinLowEdge(hOut->GetNbinsX()-2);
	Double_t y1=h->GetBinContent(hOut->GetNbinsX()-2);
	Double_t y2=h->GetBinContent(hOut->GetNbinsX()-1);
	Double_t y3=h->GetBinContent(hOut->GetNbinsX());

	h->SetBinContent(hOut->GetNbinsX()-2,y3);
	if(h->GetBinContent(1)<median+2.5*rms)h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(1));
	else h->SetBinContent(hOut->GetNbinsX()-1,median+1.5*rms);
	if(h->GetBinContent(2)<median+2.5*rms)h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(2));
	else h->SetBinContent(hOut->GetNbinsX(),median+1.5*rms);

// 	h->SetBinContent(hOut->GetNbinsX()-2,y3);
// 	h->SetBinContent(hOut->GetNbinsX()-1,h->GetBinContent(1));
// 	h->SetBinContent(hOut->GetNbinsX(),h->GetBinContent(2));
	f->SetParameter(0,h->GetBinContent(hOut->GetNbinsX()-1));
	if(fSystUseRMSfromFlatDistr==4||fSystUseRMSfromFlatDistr==14)	f->SetParameter(1,(h->GetBinContent(hOut->GetNbinsX())-h->GetBinContent(hOut->GetNbinsX()-2))/(h->GetBinCenter(hOut->GetNbinsX())-h->GetBinCenter(hOut->GetNbinsX()-2)));
	h->Fit(f,"RWW","N",xmin,xmax);
	hOut->SetBinContent(j,f->Eval(h->GetBinCenter(hOut->GetNbinsX()-1)));
	h->SetBinContent(hOut->GetNbinsX()-2,y1);
	h->SetBinContent(hOut->GetNbinsX()-1,y2);
	h->SetBinContent(hOut->GetNbinsX(),y3);
      }
      
    }    
    // Smoothing done, now simmetrize if requested
    if(fSystUseRMSfromFlatDistr==13||fSystUseRMSfromFlatDistr==14){
      TH1D *hRefl=ReflectHisto(hOut);
      TH1D *hBackTo2Pi=DuplicateHistoTo2piRange(hRefl); // In this way we symmetrized to 0-pi    
      delete[] index;
      delete h;
      delete hRefl;      
      delete hOut;
      delete f;
      hBackTo2Pi->SetName("hRelSystFeedDownMin");
      return hBackTo2Pi;
    }
    else{
      delete[] index;
      delete h;
      delete f;
      hOut->SetName("hRelSystFeedDownMin");
      return hOut;
    }
        
    
  }

  Printf("AliHFCorrelationFDsubtraction: wrong option for getting syst unc");
  return 0x0;
  
  
 
}


void AliHFCorrelationFDsubtraction::CalculateEnvelope(){

  if(fcAllRatio)delete fcAllRatio;
  if(fhEnvelopeMax)delete fhEnvelopeMax;
  if(fhEnvelopeMin)delete fhEnvelopeMin;
  if(fgrEnvelope)delete   fgrEnvelope;

  for(Int_t j=0;j<fLastTemplAdded;j++){
    SubtractFeedDown(fhTemplates[j]);
  }
  
  fhEnvelopeMaxRatio=(TH1D*)fhDataCorrected[0]->Clone("fhEnvelopeMaxRatio");
  fhEnvelopeMaxRatio->Reset();
  for(Int_t j=1;j<=  fhEnvelopeMaxRatio->GetNbinsX();j++)  fhEnvelopeMaxRatio->SetBinContent(j,1);
  fhEnvelopeMinRatio=(TH1D*)fhDataCorrected[0]->Clone("fhEnvelopeMinRatio");
  fhEnvelopeMinRatio->Reset();
  for(Int_t j=1;j<=  fhEnvelopeMinRatio->GetNbinsX();j++)  fhEnvelopeMinRatio->SetBinContent(j,1);
  fgrEnvelopeRatio=new TGraphAsymmErrors();

  fhEnvelopeMax=(TH1D*)fhDataCorrected[0]->Clone("fhEnvelopeMax");
  fhEnvelopeMin=(TH1D*)fhDataCorrected[0]->Clone("fhEnvelopeMin");
  fgrEnvelope=new TGraphAsymmErrors();
  

  Color_t col[10]={kGray,kRed,kBlue,kGreen,kOrange,kViolet,kCyan,kYellow,kSpring,kMagenta};
  fcAllRatio=new TCanvas("cRatioAll","cRatioAll",700,700);
  fcAllRatio->cd();
  // templ 0 is the reference 
  fhRatioAllToCentralTempl[0]=(TH1D*)fhDataCorrected[1]->Clone(Form("%sRatioToCentralPred",fhDataCorrected[1]->GetName()));
  fhRatioAllToCentralTempl[0]->Divide(fhDataCorrected[0]);
  
  fhRatioAllToCentralTempl[1]=(TH1D*)fhDataCorrected[2]->Clone(Form("%sRatioToCentralPred",fhDataCorrected[2]->GetName()));
  fhRatioAllToCentralTempl[1]->Divide(fhDataCorrected[0]);
  

  for(Int_t jb=1;jb<=fhRatioAllToCentralTempl[0]->GetNbinsX();jb++){
      
    fhRatioAllToCentralTempl[0]->SetBinError(jb,0.001);

    if(fhRatioAllToCentralTempl[0]->GetBinContent(jb)>1.){
      fhEnvelopeMaxRatio->SetBinContent(jb,fhRatioAllToCentralTempl[0]->GetBinContent(jb));    
      fhEnvelopeMax->SetBinContent(jb,fhDataCorrected[1]->GetBinContent(jb));
    }
    if(fhRatioAllToCentralTempl[0]->GetBinContent(jb)<1.){
      fhEnvelopeMinRatio->SetBinContent(jb,fhRatioAllToCentralTempl[0]->GetBinContent(jb));
      fhEnvelopeMin->SetBinContent(jb,fhDataCorrected[1]->GetBinContent(jb));
    }
    
    fhRatioAllToCentralTempl[1]->SetBinError(jb,0.001);
    if(fhRatioAllToCentralTempl[1]->GetBinContent(jb)>fhEnvelopeMaxRatio->GetBinContent(jb)){
	fhEnvelopeMaxRatio->SetBinContent(jb,fhRatioAllToCentralTempl[1]->GetBinContent(jb));
	fhEnvelopeMax->SetBinContent(jb,fhDataCorrected[2]->GetBinContent(jb));
      }
    if(fhRatioAllToCentralTempl[1]->GetBinContent(jb)<fhEnvelopeMinRatio->GetBinContent(jb)){
	fhEnvelopeMinRatio->SetBinContent(jb,fhRatioAllToCentralTempl[1]->GetBinContent(jb));
	fhEnvelopeMin->SetBinContent(jb,fhDataCorrected[2]->GetBinContent(jb));
    }
  }


  fhRatioAllToCentralTempl[0]->Draw();
  fhRatioAllToCentralTempl[1]->Draw("same");

  for(Int_t j=1;j<fCountTempl;j++){
    // templ j, central
    fhRatioAllToCentralTempl[j*3-1]=(TH1D*)fhDataCorrected[j*3]->Clone(Form("%sRatioToCentralPred",fhDataCorrected[j*3]->GetName()));
    fhRatioAllToCentralTempl[j*3-1]->Divide(fhDataCorrected[0]);
    // templ j, min and max
    fhRatioAllToCentralTempl[j*3]=(TH1D*)fhDataCorrected[j*3+1]->Clone(Form("%sRatioToCentralPred",fhDataCorrected[j*3+1]->GetName()));
    fhRatioAllToCentralTempl[j*3]->Divide(fhDataCorrected[0]);

    fhRatioAllToCentralTempl[j*3+1]=(TH1D*)fhDataCorrected[j*3+2]->Clone(Form("%sRatioToCentralPred",fhDataCorrected[j*3+2]->GetName()));
    fhRatioAllToCentralTempl[j*3+1]->Divide(fhDataCorrected[0]);

    for(Int_t jb=1;jb<=fhRatioAllToCentralTempl[j*3-1]->GetNbinsX();jb++){
      
      fhRatioAllToCentralTempl[j*3-1]->SetBinError(jb,0.001);
      if(fhRatioAllToCentralTempl[j*3-1]->GetBinContent(jb)>fhEnvelopeMaxRatio->GetBinContent(jb)){
	fhEnvelopeMaxRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3-1]->GetBinContent(jb));
	fhEnvelopeMax->SetBinContent(jb,fhDataCorrected[j*3-1]->GetBinContent(jb));
      }
      if(fhRatioAllToCentralTempl[j*3-1]->GetBinContent(jb)<fhEnvelopeMinRatio->GetBinContent(jb)){
	fhEnvelopeMinRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3-1]->GetBinContent(jb));
	fhEnvelopeMin->SetBinContent(jb,fhDataCorrected[j*3-1]->GetBinContent(jb));
      }

      fhRatioAllToCentralTempl[j*3]->SetBinError(jb,0.001);
      if(fhRatioAllToCentralTempl[j*3]->GetBinContent(jb)>fhEnvelopeMaxRatio->GetBinContent(jb)){
	fhEnvelopeMaxRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3]->GetBinContent(jb));
	fhEnvelopeMax->SetBinContent(jb,fhDataCorrected[j*3]->GetBinContent(jb));
      }
      if(fhRatioAllToCentralTempl[j*3]->GetBinContent(jb)<fhEnvelopeMinRatio->GetBinContent(jb)){
	fhEnvelopeMinRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3]->GetBinContent(jb));
	fhEnvelopeMin->SetBinContent(jb,fhDataCorrected[j*3]->GetBinContent(jb));
      }

      fhRatioAllToCentralTempl[j*3+1]->SetBinError(jb,0.001);
      if(fhRatioAllToCentralTempl[j*3+1]->GetBinContent(jb)>fhEnvelopeMaxRatio->GetBinContent(jb)){
	fhEnvelopeMaxRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3+1]->GetBinContent(jb));
	fhEnvelopeMax->SetBinContent(jb,fhDataCorrected[j*3+1]->GetBinContent(jb));
      }
      if(fhRatioAllToCentralTempl[j*3+1]->GetBinContent(jb)<fhEnvelopeMinRatio->GetBinContent(jb)){
	fhEnvelopeMinRatio->SetBinContent(jb,fhRatioAllToCentralTempl[j*3+1]->GetBinContent(jb));
	fhEnvelopeMin->SetBinContent(jb,fhDataCorrected[j*3+1]->GetBinContent(jb));
      }

     
    }
  
    if(j<=10){
      fhRatioAllToCentralTempl[j*3-1]->SetLineColor(col[j-1]);
      fhRatioAllToCentralTempl[j*3]->SetLineColor(col[j-1]+1);
      fhRatioAllToCentralTempl[j*3+1]->SetLineColor(col[j-1]+2);
    }
    else if(j<=20){
      fhRatioAllToCentralTempl[j*3-1]->SetLineColor(col[j-10]+3);
      fhRatioAllToCentralTempl[j*3]->SetLineColor(col[j-10]+4);
      fhRatioAllToCentralTempl[j*3+1]->SetLineColor(col[j-10]+5);
    }
    else{
      Printf("Too many templates, add more colors, aborting");
      return;
    }
    fhRatioAllToCentralTempl[j*3-1]->Draw("same");
    fhRatioAllToCentralTempl[j*3]->Draw("same");
    fhRatioAllToCentralTempl[j*3+1]->Draw("same");
  }
  
  for(Int_t j=1;j<=fhEnvelopeMin->GetNbinsX();j++){
    fhEnvelopeMin->SetBinError(j,0);
    fhEnvelopeMinRatio->SetBinError(j,0);
    fhEnvelopeMax->SetBinError(j,0);
    fhEnvelopeMaxRatio->SetBinError(j,0);
  }
  
  fgrEnvelope=GetUncGraphFromHistos(fhDataCorrected[0],fhEnvelopeMin,fhEnvelopeMax);
  fgrEnvelopeRatio=GetUncGraphFromHistos(0x0,fhEnvelopeMinRatio,fhEnvelopeMaxRatio);
  
  
}


// void AliHFCorrelationFDsubtraction::GetTemplateFromFit(TH1D *h,TH1D *hOut,TString strCanv="cFit",Int_t methodFD=0){

//   // Fit Template from B
//   TCanvas *cFit=new TCanvas(strCanv.Data(),strCanv.Data(),700,700);
//   cFit->cd();
//   h->Draw();
  
//   if(methodFD==0){
//     TF1 *fitFunction=FitPlotsShort(h,2,0,0);
//     for(Int_t j=1;j<=h->GetNbinsX();j++){
//       hOut->SetBinContent(j,fitFunction->Eval(hOut->GetBinCenter(j)));
//       hOut->SetBinError(j,0);
//     }
//     hOut->SetLineColor(kViolet);
//     hOut->Draw("same");
    
//   }
//   else if(methodFD==1){
    

//     Double_t nsybc, ensybc,asybc, easybc;
//     Double_t nsybcMC, ensybcMC,asybcMC, easybcMC;
//     cFit->cd();
//     TF1 *fitFunctionMC=FitPlots(h,2,2,3,nsybcMC, ensybcMC,asybcMC, easybcMC);
//     Double_t baseMC=fitFunctionMC->GetParameter(0);

//     TCanvas *cFitData=new TCanvas(Form("Data%s",strCanv.Data()),"cGetTemplatefromFit_FitData",700,700);
//     cFitData->cd();
    
//     TH1D *hDataForFit=(TH1D*)hOut->Clone("hDataForFit");
//     hDataForFit->Draw("same");

//     TF1 *fitFunctionData=FitPlots(hDataForFit,1,2,3,nsybc, ensybc,asybc, easybc);
//     Double_t baseData=fitFunctionData->GetParameter(0);

//     for(Int_t j=1;j<=h->GetNbinsX();j++){
//       hOut->SetBinContent(j,fitFunctionMC->Eval(hOut->GetBinCenter(j))+baseData-baseMC);
//       hOut->SetBinError(j,0);
//     }
    
//     cFit->cd();
//     hOut->SetLineColor(kViolet);
//     hOut->DrawClone("same");
    
//   }

//   return;
// }



// TH1D* AliHFCorrelationFDsubtraction::SubtracFDrough(Double_t ptmin,Double_t ptmax,Int_t methodSubtr=0,Int_t methodFD=0,Int_t rebin=1,Int_t testAssYieldExt=0,TString correlationDataFile="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2013Jul1PFapprovalNoPrelim/FabioD0/1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated_Bins10to11_Limits_2_4_TreshPt_0.3_Displ_0.00_Data.root", TString correlationFDtemplFile="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/MCtemplate/2013June14Fabio/Plots/1D_Signal_WithEMCorr_MCSample_Normal_Charg_OriginSuper_Integrated_Bins10to11_Limits_2_4_TreshPt_0.3_Displ_0.00_Kine.root",TString spectraMacroOutput="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/FeedDownSubtraction/ThirdTest/2013Jun17FabioD0properEff/Dati_CorrectEff/HFPtSpectrum_Nb.root"){

//   /*TString correlationDataFile="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/FeedDownSubtraction/SecondTest/PlotCorrelDeta1/1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated_Bins10to11_Limits_2_4_TreshPt_0.3_Displ_0.00_Data.root",TString correlationFDtemplFile="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/MCtemplate/2013June17Fabio/Plots/0.3/1D_Signal_WithEMCorr_MCSample_Normal_Charg_OriginSuper_Integrated_Bins10to11_Limits_2_4_TreshPt_0.3_Displ_0.00_Kine.root",TString spectraMacroOutput="../HFPtSpectrum_Nb.root"){
//    */
  
//   gROOT->LoadMacro("/Users/administrator/ALICE/CHARM/HFCJ/Macro/FitPlots.C");
  
//   TFile *fDataCorr=TFile::Open(correlationDataFile.Data(),"READ");
//   TCanvas *cData=fDataCorr->Get("c3");
//   TH1D *hData=(TH1D*)cData->FindObject("hsubtract_norm");
//   hData->Sumw2();
//   hData->Rebin(rebin);
//   hData->Scale(1./(Double_t)rebin);
  
//   TFile *fFDtemplCorr=TFile::Open(correlationFDtemplFile.Data(),"READ");
//   TCanvas *cFDtempl=fFDtemplCorr->Get("c3");
//   TH1D *hFDtemplFile=(TH1D*)cFDtempl->FindObject("hsubtract_norm_2");
//   TH1D *hPromptTempl=(TH1D*)cFDtempl->FindObject("hsubtract_norm_1");
//   TH1D *hFDtempl;
//   if(methodSubtr==0)hFDtempl=(TH1D*)hFDtemplFile->Clone("hTemplDfromB");
//   else if(methodSubtr==1){
//     hFDtempl=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
//     GetTemplateFromFit(hFDtemplFile,hFDtempl,"cFitFD");
//   }
//   else if(TMath::Abs(methodSubtr)==2){// Increase the baseline of MC template in order to match the baseline observed in data
//     hFDtempl=(TH1D*)hData->Clone("hTemplDfromBFitFunc");
//     GetTemplateFromFit(hFDtemplFile,hFDtempl,"cFitFD",1);
//   }

//   cData->Draw();
//   hFDtempl->SetLineColor(kBlue);
//   hFDtempl->SetMarkerColor(kBlue);
//   hFDtempl->Draw("same");


//   TFile *fSpectrum=TFile::Open(spectraMacroOutput.Data(),"READ");
//   TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fSpectrum->Get("gFcConservative");
//   //   TAxis *ax=gr->GetXaxis();
//   //   Int_t binmin=ax->FindBin(ptmin*1.0001);
//   //   Int_t binmax=ax->FindBin(ptmax*0.9999);
  
  

//   Double_t *xgr=gr->GetX();
//   Double_t *elxgr=gr->GetEXlow();
//   Double_t *ehxgr=gr->GetEXhigh();
//   Double_t *fprompt=gr->GetY();
//   Double_t *elfprompt=gr->GetEYlow();
//   Double_t *ehfprompt=gr->GetEYhigh();

//   Int_t binmin=-1,binmax=-1;
//   // Look for right bins
//   for(Int_t j=0;j<gr->GetN();j++){
//     if(binmin<0&&ptmin*1.00001>(xgr[j]-elxgr[j])&&ptmin<(xgr[j]+ehxgr[j]))binmin=j;
//     if(binmax<0&&ptmax>(xgr[j]-elxgr[j])&&ptmax*0.9999<(xgr[j]+ehxgr[j]))binmax=j;

//   }

//   if(TMath::Abs(xgr[binmin]-elxgr[binmin]-ptmin)>0.0001){
//     printf("Bin mismatch: xmin=%f \n",xgr[binmin]-elxgr[binmin]);
    
//     return;
//   }

  
//   // DEBUG BIAS:
//   //  fprompt[binmin]=1.;
//   //  elfprompt[binmin]=0.;
//   //  ehfprompt[binmin]=0.;

//   //  hFDtempl->Scale(1./5.);
  
//   TH1D *hFDCentr=(TH1D*)hFDtempl->Clone("hFDtemplCentr");
//   hFDCentr->Scale(1.-fprompt[binmin]);
//   printf("Centr fsec=%f \n",1.-fprompt[binmin]);

//   TH1D *hDataCentr=(TH1D*)hData->Clone("hDataCentr");
//   hDataCentr->Add(hFDCentr,-1.);// Treatment of Stat unc.
//   hDataCentr->Scale(1./fprompt[binmin]);
//   hDataCentr->SetTitle("central f_{prompt}");

//   TH1D *hFDmax=(TH1D*)hFDtempl->Clone("hFDtemplMax");
//   hFDmax->Scale(1.-(fprompt[binmin]-elfprompt[binmin]));
//   printf("Max fsec=%f \n",1.-(fprompt[binmin]-elfprompt[binmin]));
//   TH1D *hDataMaxFD=(TH1D*)hData->Clone("hDataMaxFD");
//   hDataMaxFD->Add(hFDmax,-1.);// Treatment of Stat unc.
//   hDataMaxFD->Scale(1./(fprompt[binmin]-elfprompt[binmin]));
//   hDataMaxFD->SetTitle("max f_{prompt}");

//   TH1D *hFDmin=(TH1D*)hFDtempl->Clone("hFDtemplMin");
//   hFDmin->Scale(1.-(fprompt[binmin]+ehfprompt[binmin]));
//   printf("Min fsec=%f \n",1.-(fprompt[binmin]+ehfprompt[binmin]));
//   TH1D *hDataMinFD=(TH1D*)hData->Clone("hDataMinFD");
//   hDataMinFD->Add(hFDmin,-1.);// Treatment of Stat unc.
//   hDataMinFD->Scale(1./(fprompt[binmin]+ehfprompt[binmin]));
//   hDataMinFD->SetTitle("min f_{prompt}");


    
//   TCanvas *cCompare=new TCanvas("cCompare","cCompare",700,700);
//   cCompare->cd();
//   hData->Draw();
//   hData->SetTitle("Uncorrected");
//   hData->SetYTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#varphi} (rad^{-1})");
//   hData->SetXTitle("#varphi_{ass}-#varphi_{trig} (rad)");
//   hData->SetMarkerStyle(25);
//   hData->SetMarkerColor(kBlack);
//   hData->SetMarkerSize(1);
//   hData->SetLineColor(kBlack);

//   hDataCentr->SetLineColor(kRed);
//   hDataCentr->SetYTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#varphi} (rad^{-1})");
//   hDataCentr->SetXTitle("#varphi_{ass}-#varphi_{trig} (rad)");
//   hDataCentr->SetLineWidth(2);
//   hDataCentr->SetMarkerColor(kRed);
//   hDataCentr->SetMarkerStyle(21);
//   hDataCentr->SetMarkerSize(1);
//   hDataCentr->Draw("sames");

//   hDataMinFD->SetLineColor(kBlue);
//   hDataMinFD->SetMarkerColor(kBlue);
//   hDataMinFD->SetMarkerStyle(23);
//   hDataMinFD->SetMarkerSize(1);
//   hDataMinFD->Draw("sames");

//   hDataMaxFD->SetLineColor(kGreen-6);
//   hDataMaxFD->SetMarkerColor(kGreen-6);
//   hDataMaxFD->SetMarkerStyle(22);
//   hDataMaxFD->SetMarkerSize(1);
//   hDataMaxFD->Draw("sames");
  
//   TCanvas *cRatio=new TCanvas("cRatio","cRatio",700,700);
//   cRatio->cd();
//   TH1D *hRatioMinToCentr=(TH1D*)hDataMinFD->Clone("hRatioMinToCentr");
//   hRatioMinToCentr->Divide(hDataCentr);
//   for(Int_t j=1;j<=hRatioMinToCentr->GetNbinsX();j++){
//     hRatioMinToCentr->SetBinError(j,0.001);
//   }
//   TH1D *hRatioMaxToCentr=(TH1D*)hDataMaxFD->Clone("hRatioMaxToCentr");
//   hRatioMaxToCentr->Divide(hDataCentr);
//   for(Int_t j=1;j<=hRatioMaxToCentr->GetNbinsX();j++){
//     hRatioMaxToCentr->SetBinError(j,0.001);
//   }

//   hRatioMinToCentr->Draw();
//   hRatioMaxToCentr->Draw("same");



//   TCanvas *cCompareToMC=new TCanvas("cCompareToMC","cCompareToMC",700,700);
//   cCompareToMC->cd();
//   TH1D *hDataCentrCompare=hDataCentr->DrawClone();
//   hDataCentrCompare->SetLineColor(kBlack);
//   hDataCentrCompare->SetMarkerColor(kBlack);
//   hPromptTempl->Draw("same");

//   TCanvas *cCompareToMCFit=new TCanvas("cCompareToMCFit","cCompareToMCFit",700,700);
//   cCompareToMCFit->cd();
//   TH1D *hDataCentrCompare2=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompare2->SetLineColor(kBlack);
//   hDataCentrCompare2->SetMarkerColor(kBlack);

//   TH1D *hPromptTemplFit=(TH1D*)hDataCentr->Clone("hPromptDtemplatFit");
//   GetTemplateFromFit(hPromptTempl,hPromptTemplFit,"cFitPrompt",1);
//   cCompareToMCFit->cd();
//   hPromptTemplFit->Draw("same");

//   TCanvas *cRatioToMCFit=new TCanvas("cRatioToMCFit","cRatioToMCFit",700,700);
//   cRatioToMCFit->cd();
//   TH1D *hDataCentrRatio=(TH1D*)hDataCentr->Clone("hDataCentrRatio");
//   hDataCentrRatio->SetLineColor(kBlack);
//   hDataCentrRatio->SetMarkerColor(kBlack);
//   hDataCentrRatio->Divide(hPromptTemplFit);
//   hDataCentrRatio->Draw();

//   TCanvas *cCompareFitToFit=new TCanvas("cCompareFitToFit","cCompareFitToFit",700,700);
//   cCompareFitToFit->cd();
//   TH1D *hDataCentrCompare3=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompare3->SetLineColor(kBlack);
//   hDataCentrCompare3->SetMarkerColor(kBlack);

//   Double_t nsybcSt=0.,ensybcSt=0.,asybcSt=0.,easybcSt=0.; 
//   Double_t nsySt=0., ensySt=0.,asySt=0., easySt=0.,sigmaNSSt=0.,esigmaNSSt=0.;
//   TF1 *fitSt=FitPlots(hDataCentrCompare2,1,0,0,   nsybcSt, ensybcSt,asybcSt, easybcSt);
//   nsySt=fitSt->GetParameter(1);
//   ensySt=fitSt->GetParError(1);
//   asySt=fitSt->GetParameter(4);
//   easySt=fitSt->GetParError(4);
//   sigmaNSSt=fitSt->GetParameter(3);
//   esigmaNSSt=fitSt->GetParameter(3);

//   hPromptTempl->Draw("sames");

//   if(!testAssYieldExt)  return hDataCentr;

//   // NOW SECTION ON ASSOCIATED YIELD EXTRACTION 
//   Double_t nsybc=0., ensybc=0.,asybc=0.,easybc=0.; 
//   Double_t nsy=0., ensy=0.,asy=0., easy=0.,sigmaNS=0.,esigmaNS=0.;

//   TH1D *hResidualNSyield=new TH1D("hResidualNSyield","hResidualNSyield",200,0.,10.);
//   TH1D *hResidualASyield=new TH1D("hResidualASyield","hResidualASyield",200,0.,10.);
//   TH1D *hResidualNSsigma=new TH1D("hResidualNSsigma","hResidualNSsigma",200,0.,4.);
//   TH1D *hResidualASsigma=new TH1D("hResidualASsigma","hResidualASsigma",200,0.,4.);

//   TH1D *hResidualNSyieldBC=new TH1D("hResidualNSyieldBC","hResidualNSyieldBC",200,0.,10.);
//   TH1D *hResidualASyieldBC=new TH1D("hResidualASyieldBC","hResidualASyieldBC",200,0.,10.);

//   TCanvas *cResiduals=new TCanvas("cResiduals","cResiduals",700,700);
//   cResiduals->Divide(2,2);
//   cResiduals->cd(1);

//   hResidualNSyield->Draw();
//   hResidualNSyieldBC->SetLineColor(kRed);
//   hResidualNSyieldBC->Draw("same");

//   cResiduals->cd(2);

//   hResidualASyield->Draw();
//   hResidualASyieldBC->SetLineColor(kRed);
//   hResidualASyieldBC->Draw("same");

//   cResiduals->cd(3);
//   hResidualNSsigma->Draw();

//   // START TESTS

//   // FIRST VALUES: STANDARD
//   hResidualNSyield->Fill(nsybcSt);  
//   hResidualASyield->Fill(asybcSt);  
//   hResidualNSyield->Fill(nsySt);  
//   hResidualASyield->Fill(asySt);  
//   hResidualNSsigma->Fill(sigmaNSSt);  


//   hResidualNSyieldBC->Fill(nsybcSt);  
//   hResidualASyieldBC->Fill(asybcSt);  


//   // FIX  Baseline at histo min 
//   TCanvas *cYieldExtrFixBaselineMin=new TCanvas("cYieldExtrFixBaselineMin","cYieldExtrFixBaselineMin",700,700);
//   cYieldExtrFixBaselineMin->cd();

//   TH1D *hDataCentrCompareBasMin=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareBasMin->SetLineColor(kBlack);
//   hDataCentrCompareBasMin->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareBasMin,1,1,0,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  

//   printf("Fix baseline at min ok \n");


//   // FIX  Baseline at average of lower 3 histo points 
//   TCanvas *cYieldExtrFixBaseAv3=new TCanvas("cYieldExtrFixBaseAv3","cYieldExtrFixBaseAv3",700,700);
//   cYieldExtrFixBaseAv3->cd();
//   TH1D *hDataCentrCompareBasAv3=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareBasAv3->SetLineColor(kBlack);
//   hDataCentrCompareBasAv3->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareBasAv3,1,-3,0,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  

//   printf("Fix baseline at av lower3 ok \n");

//   // FIX  Baseline at average of lower 2 histo points 
//   TCanvas *cYieldExtrFixBaseAv2=new TCanvas("cYieldExtrFixBaseAv2","cYieldExtrFixBaseAv2",700,700);
//   cYieldExtrFixBaseAv2->cd();
//   TH1D *hDataCentrCompareBasAv2=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareBasAv2->SetLineColor(kBlack);
//   hDataCentrCompareBasAv2->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareBasAv2,1,-2,0,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  


//   printf("Fix baseline at av lower2 ok \n");

//   // FIX  BOTH MEANS
//   TCanvas *cYieldExtrFixBothMeans=new TCanvas("cYieldExtrFixBothMeans","cYieldExtrFixBothMeans",700,700);
//   cYieldExtrFixBothMeans->cd();
//   TH1D *hDataCentrCompareFMB=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareFMB->SetLineColor(kBlack);
//   hDataCentrCompareFMB->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareFMB,1,0,3,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  

//   printf("Fix both means ok \n");
//   return;

//   // FIX  NS MEAN
//   TCanvas *cYieldExtrFixNSMean=new TCanvas("cYieldExtrFixNSMean","cYieldExtrFixNSMean",700,700);
//   cYieldExtrFixBothMeans->cd();
//   TH1D *hDataCentrCompareFMNS=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareFMNS->SetLineColor(kBlack);
//   hDataCentrCompareFMNS->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareFMNS,1,0,1,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  

//   printf("Fix NS mean ok \n");

//   // FIX  AS MEAN
//   TCanvas *cYieldExtrFixASMean=new TCanvas("cYieldExtrFixASMean","cYieldExtrFixASMean",700,700);
//   cYieldExtrFixBothMeans->cd();
//   TH1D *hDataCentrCompareFMAS=(TH1D*)hDataCentr->DrawClone();
//   hDataCentrCompareFMAS->SetLineColor(kBlack);
//   hDataCentrCompareFMAS->SetMarkerColor(kBlack);

//   TF1 *fitTest=FitPlots(hDataCentrCompareFMAS,1,0,2,nsybc,ensybc,asybc,easybc);
//   hPromptTempl->Draw("sames");

//   nsy=fitTest->GetParameter(1);
//   ensy=fitTest->GetParError(1);
//   asy=fitTest->GetParameter(4);
//   easy=fitTest->GetParError(4);
//   sigmaNS=fitTest->GetParameter(3);
//   esigmaNS=fitTest->GetParameter(3);

//   hResidualNSyield->Fill(nsybc);  
//   hResidualASyield->Fill(asybc);  
//   hResidualNSyield->Fill(nsy);  
//   hResidualASyield->Fill(asy);  
//   hResidualNSsigma->Fill(sigmaNS);  

//   hResidualNSyieldBC->Fill(nsybc);  
//   hResidualASyieldBC->Fill(asybc);  

//   printf("Fix AS mean ok \n");


// }
