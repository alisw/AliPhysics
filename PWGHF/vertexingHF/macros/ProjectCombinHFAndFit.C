#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>

#include "AliHFMassFitter.h"
#include "AliNormalizationCounter.h"
#endif

enum Method{kME,kRot,kLS,kSB};

// input files and pt binning
TString fileName="AnalysisResults_train543544.root";
TString suffix="c3SigPID_Pt400_EM1";
TString meson="Dzero";
const Int_t nPtBins=6;
Double_t binLims[nPtBins+1]={0.,1.,2.,3.,4.,6.,8.};
Double_t sigmas[nPtBins]={0.0075,0.008,0.009,0.009,0.010,0.012};//,0.012,0.015,0.017};// MC values up to 8

// fit configuration
Int_t rebin=5;
Bool_t fixSigma=kTRUE;
Bool_t fixMean=kFALSE;
Double_t minMass=1.74;
Double_t maxMass=2.04;
Int_t optForNorm=1;
Bool_t useEMwithLS=kTRUE;
Int_t typeb=2;
Double_t nsigmaBinCounting=4.;      // defines the mass interval over which the signal is bin counted
Double_t massD;

// objects and options related to side-band fit method
Bool_t tryDirectFit=kTRUE;
TCanvas *cDataSubtractedFit;
TCanvas *cDataSubtractedReFit;
TCanvas **cDataFitRes;
TCanvas *cPullsDirect;
TCanvas *cPullsDirectTrendVsMass;
TCanvas *cResidualsDirect;
TCanvas *cResidualDirectTrendVsMass;

TF1 **fFitSB;
TF1 *fback;// pointer to "current" background function
Int_t nparback=0;
Double_t fitrangelow[nPtBins]={1.74,1.74,1.74,1.72,1.72,1.72};
Double_t fitrangeup[nPtBins]={2.0,2.0,2.0,2.0,2.0,2.0};
Bool_t refitWithSignal=kTRUE;  // kFALSE: the background is obtained by fitting the side-bands only and then subtracted. The residual distribution is fittied with a pol0+gaus 
// using the min chi2 ROOT option
                               // kTRUE: the fit of the side-bands is used to initialize the background function parameters, then the whole distribution is fitted with sign+back, susing the max likelihood option accounting properly for poissonian uncertainties
Bool_t fixBackParFromFirstFit=kFALSE; // if kTRUE, in the case refitWithSignal=kTRUE the parameter related to the bacgkround part of the final fit function are fixed 
                                     // from the fit of the side-bands. (n.b. allows to use the max likelihood option but it is not that useful, better use refitWithSignal=kFALSE)
Double_t nsigmaSubBackFit=2.5;      // defines the mass window for excluding the signal region for the side-band fit
Int_t nDegreeBackPol=4;             // degree of polynomial function describing the background

// SMEARING PART

void WriteFitInfo(AliHFMassFitter *fitter, TH1D* histo);

AliHFMassFitter* ConfigureFitter(TH1D* histo, Int_t iPtBin){
  TH1F* histof=(TH1F*)histo->Clone(Form("%s_Fl",histo->GetName()));
  AliHFMassFitter* fitter=new AliHFMassFitter(histof,minMass,maxMass,1,typeb,0);
  fitter->SetReflectionSigmaFactor(0);
  fitter->SetUseChi2Fit();
  fitter->SetInitialGaussianMean(massD);
  fitter->SetInitialGaussianSigma(sigmas[iPtBin]);
  if(fixSigma) fitter->SetFixGaussianSigma(sigmas[iPtBin]);
  if(fixMean) fitter->SetFixGaussianMean(massD);
  return fitter;
}

Double_t GetSignalBinCounting(TH1 *h,TF1 *fbackground,Double_t &err,Double_t nsigmaBC=4.,Double_t sigmafit=0.010,Double_t minx=-999,Double_t maxx=-999){// default range for counting is pdg mass +- nsigmaBC x sigmafit, unless minx and maxx are specified
  
  Int_t binL,binR;
  if(minx>=0.){
    binL=h->FindBin(minx);
    if(binL<1){
      Printf("Left range for bin counting smaller than allowed by histogram axis, setting it to the minimum ");
      binL=1;
    }
  }
  else {
    Double_t lr=massD-sigmafit*nsigmaBC;
    if(minx<0&&minx>-999.){
      if(lr<-1.*minx){
	Printf("Left range for counting set to min of value defined");
	lr=-1.*minx;
      }
    }
    binL=h->FindBin(lr);

    if(binL<1){
      Printf("Left range for bin counting smaller than allowed by histogram axis, setting it to the minimum ");
      binL=1;
    }
  }

  if(maxx>=0.){
    binR=h->FindBin(maxx);
    if(binR>h->GetNbinsX()){
      Printf("Right range for bin counting larger than allowed by histogram axis, setting it to the maximum allowed ");
      binR=h->GetNbinsX();
    }
  }
  else {
    Double_t rr=massD+sigmafit*nsigmaBC;
    if(maxx<0. && maxx>-999){
      if(rr>-1.*maxx){
	Printf("Right range for counting set to max of value defined");
	rr=-1.*maxx;
      }
    }
    binR=h->FindBin(rr);
    if(binR>h->GetNbinsX()){
      Printf("Right range for bin counting larger than allowed by histogram axis, setting it to the maximum allowed ");
      binR=h->GetNbinsX();
    }
  }

  Double_t sign=0;
  err=0;
  for(Int_t j=binL;j<=binR;j++){
    sign+=h->GetBinContent(j)-fbackground->Integral(h->GetBinLowEdge(j),h->GetBinLowEdge(j)+h->GetBinWidth(j))/h->GetBinWidth(j);
    err+=h->GetBinError(j)*h->GetBinError(j);
  }
  err=TMath::Sqrt(err);
  return sign; 
}


TH1F* GetResidualsAndPulls(TH1 *h,TF1 *f,Double_t minrange=0,Double_t maxrange=-1,TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0){
  Int_t binmi=1,binma=h->GetNbinsX();

  if(maxrange>minrange){
    binmi=h->FindBin(minrange*1.001);
    binma=h->FindBin(maxrange*0.9999);
  }
  if(hResidualTrend){
    //h->Copy(hResidualTrend);
    hResidualTrend->SetBins(h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
    hResidualTrend->SetName(Form("%s_residualTrend",h->GetName()));
    hResidualTrend->SetTitle(Form("%s ; Invariant Mass (GeV/c^{2}) ; Residuals",h->GetTitle()));
    hResidualTrend->SetMarkerStyle(20);
    hResidualTrend->SetMarkerSize(1.0);
    hResidualTrend->Reset();
  }
  if(hPullsTrend){
    hPullsTrend->SetBins(h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
    hPullsTrend->Reset();
    hPullsTrend->SetName(Form("%s_pullTrend",h->GetName()));
    hPullsTrend->SetTitle(Form("%s ; Invariant Mass (GeV/c^{2}) ; Pulls",h->GetTitle()));
    hPullsTrend->SetMarkerStyle(20);
    hPullsTrend->SetMarkerSize(1.0);
  }
  if(hPulls){
    hPulls->SetName(Form("%s_pulls",h->GetName()));
    hPulls->SetTitle(Form("%s ; Pulls",h->GetTitle()));
    hPulls->SetBins(40,-10,10);
    hPulls->Reset();
  }

  Double_t res=-1.e-6,min=1.e+12,max=-1.e+12;
  TArrayD *arval=new TArrayD(binma-binmi+1);
  for(Int_t jst=binmi;jst<=binma;jst++){      
    
    res=h->GetBinContent(jst)-f->Integral(h->GetBinLowEdge(jst),h->GetBinLowEdge(jst)+h->GetBinWidth(jst))/h->GetBinWidth(jst);
    arval->AddAt(res,jst-binmi);
    if(res<min)min=res;
    if(res>max)max=res;
    //      Printf("Res = %f from %f - %f",res,h->GetBinContent(jst),f->Integral(h->GetBinLowEdge(jst),h->GetBinLowEdge(jst)+h->GetBinWidth(jst))/h->GetBinWidth(jst));
    if(hResidualTrend){
      hResidualTrend->SetBinContent(jst,res);
      hResidualTrend->SetBinError(jst,h->GetBinError(jst));
    }
    if(hPulls){
      hPulls->Fill(res/h->GetBinError(jst));
    }    
    if(hPullsTrend){
      hPullsTrend->SetBinContent(jst,res/h->GetBinError(jst));
      hPullsTrend->SetBinError(jst,0.0001);
    }
  }
  if(hPullsTrend){
    hPullsTrend->GetXaxis()->SetRange(binmi,binma);
    hPullsTrend->SetMinimum(-7);
    hPullsTrend->SetMaximum(+7);
  }
  if(TMath::Abs(min)>TMath::Abs(max))max=min;

  TH1F *hout=new TH1F(Form("%s_residuals",h->GetName()),Form("%s ; residuals",h->GetTitle()),25,-TMath::Abs(max)*1.5,TMath::Abs(max)*1.5);
  for(Int_t j=0;j<binma-binmi+1;j++){
    hout->Fill(arval->At(j));
  }
  hout->Sumw2();
  hout->Fit("gaus","LEMI","",-TMath::Abs(max)*1.2,TMath::Abs(max)*1.2);

  if(hPulls){
    hPulls->Sumw2();
    hPulls->Fit("gaus","RLEMI","",-3,3);
  }
  delete arval;
  return hout;
}
void SetStyleHisto(TH1 *h,Int_t method,Int_t isXpt=-1){
  if(isXpt==1)h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");// just for convenience, most of the time it's pt :)

  if(method==kRot){
    h->SetMarkerStyle(21);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(kBlack);
    h->SetMarkerSize(1.0);
    h->SetLineColor(kBlack);
    return;
  }

  if(method==kME){
    h->SetMarkerStyle(25);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(4);
    h->SetMarkerSize(1.0);
    h->SetLineColor(4);
    return;
  }

  if(method==kLS){
    h->SetMarkerStyle(22);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(kGreen+2);
    h->SetMarkerSize(1.0);
    h->SetLineColor(kGreen+2);
    return;
  }
  
  if(method==kSB){
    h->SetMarkerStyle(27);    
    h->GetYaxis()->SetTitleOffset(1.8);
    h->SetMarkerColor(6);
    h->SetMarkerSize(1.0);
    h->SetLineColor(6);
    return;
  }

  return;
}

void DivideCanvas(TCanvas *c,Int_t ndivisions){
  if(ndivisions>1){
    if(ndivisions<3){
      c->Divide(2,1);
    }
    else if(ndivisions<5){
      c->Divide(2,2);
    }
    else if(ndivisions<7){
      c->Divide(3,2);

    }
    else if(ndivisions<10){
      c->Divide(3,3);

    }
    else if(ndivisions<13){
      c->Divide(4,3);

    }
    else if(ndivisions<17){
      c->Divide(4,4);

    }
    else {
      c->Divide(5,5);

    }
  }
  return;
}



TF1 *GausPlusLine(Double_t minRange=1.72,Double_t maxRange=2.05){
  TF1 *f=new TF1("GausPlusLine","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]+[4]*x",minRange,maxRange);

  f->SetParameter(1,massD);
  f->SetParLimits(0,0.,1.e+7);
  f->SetParLimits(1,1.8,1.92);
  f->SetParameter(2,0.010);
  f->SetParLimits(2,0,0.1);
  f->SetParameter(3,0.);
  //  f->SetParLimits(3,-10000,10000);
  return f;

}


TH1D* FitBackgroundAndSubtract(TH1 *h,Int_t binNumber,Int_t binID=-1){
  // TO-DO: REWRITE IT IN A MORE CLEVER WAY (but inside Mass Fitter)
  
  if(binID<0)binID=binNumber;

  TF1 *fpol3 = new TF1(Form("mypol3Bin%d",binID),"pol3",0.,3.);
  fpol3->SetParameter(0,1.3e+11);
  fpol3->SetParameter(1,-0.35);

  TF1 *fpol4 = new TF1(Form("mypol4Bin%d",binID),"pol4",0.,3.);
  fpol4->SetParameter(0,1.3e+11);
  fpol4->SetParameter(1,-0.35);

  TF1 *fpol5 = new TF1(Form("mypol5Bin%d",binID),"pol5",0.,3.);
  fpol5->SetParameter(0,1.3e+11);
  fpol5->SetParameter(1,-0.35);

  TF1 *fpol6 = new TF1(Form("mypol6Bin%d",binID),"pol6",0.,3.);
  fpol6->SetParameter(0,1.3e+11);
  fpol6->SetParameter(1,-0.35);

  TF1 *fpol7 = new TF1(Form("mypol7Bin%d",binID),"pol7",0.,3.);
  fpol7->SetParameter(0,1.3e+11);
  fpol7->SetParameter(1,-0.35);

  TF1 *fpol8 = new TF1(Form("mypol8Bin%d",binID),"pol8",0.,3.);
  fpol8->SetParameter(0,1.3e+11);
  fpol8->SetParameter(1,-0.35);

  TF1 *fpol2mine = new TF1(Form("mypol2Bin%d",binID),"[0]+[1]*x+[2]*x*x",0.,3.);
  fpol2mine->SetParameter(0,1.3e+11);
  fpol2mine->SetParameter(1,-0.35);
  fpol2mine->SetParameter(2,0.001);

  TH1D *h2=(TH1D*)h->Clone(Form("%sSubtrBin%d",h->GetName(),binID));//new TH1D("h2","h2",50,1.7,2.2);
  TH1D *h3=(TH1D*)h->Clone(Form("%sCutBin%d",h->GetName(),binID));//,"h3",50,1.7,2.2);
  
  for(Int_t j=1;j<=h->GetNbinsX();j++){
    Double_t binCenter=h->GetBinCenter(j);
    Int_t j2=h2->FindBin(binCenter);
    if(binCenter>(massD-nsigmaSubBackFit*sigmas[binNumber]) && binCenter<(massD+nsigmaSubBackFit*sigmas[binNumber])){
      h3->SetBinContent(j2,0);
      h3->SetBinError(j2,0);
    }
  }
  
  if(nDegreeBackPol<2){
    Printf("%d order pol background fit not implemented: do you think it can work?",nDegreeBackPol);
    delete h2;
    delete h3;
    return 0x0;
  }
  
  if(nDegreeBackPol>=2){    
    h3->Fit(fpol2mine,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
    if(nDegreeBackPol==2){
	fback=fpol2mine;
	if(refitWithSignal)nparback=3;
    }
    if(nDegreeBackPol>=3){    
      fpol3->SetParameter(0,fpol2mine->GetParameter(0));
      fpol3->SetParameter(1,fpol2mine->GetParameter(1));
      fpol3->SetParameter(2,fpol2mine->GetParameter(2));

      h3->Fit(fpol3,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
      if(nDegreeBackPol==3){
	  fback=fpol3;
	  if(refitWithSignal)nparback=4;
      }
      
      if(nDegreeBackPol>=4){    
	fpol4->SetParameter(0,fpol3->GetParameter(0));
	fpol4->SetParameter(1,fpol3->GetParameter(1));
	fpol4->SetParameter(2,fpol3->GetParameter(2));
	fpol4->SetParameter(3,fpol3->GetParameter(3));
	h3->Fit(fpol4,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
	if(nDegreeBackPol==4){
	  fback=fpol4;
	  if(refitWithSignal)nparback=5;
	}
	
	if(nDegreeBackPol>=5){    
	  fpol5->SetParameter(0,fpol4->GetParameter(0));
	  fpol5->SetParameter(1,fpol4->GetParameter(1));
	  fpol5->SetParameter(2,fpol4->GetParameter(2));
	  fpol5->SetParameter(3,fpol4->GetParameter(3));
	  fpol5->SetParameter(4,fpol4->GetParameter(4));
	  h3->Fit(fpol5,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
	  if(nDegreeBackPol==5){
	    fback=fpol5;
	    if(refitWithSignal)nparback=6;
	  }
	  if(nDegreeBackPol>=6){    
	    fpol6->SetParameter(0,fpol5->GetParameter(0));
	    fpol6->SetParameter(1,fpol5->GetParameter(1));
	    fpol6->SetParameter(2,fpol5->GetParameter(2));
	    fpol6->SetParameter(3,fpol5->GetParameter(3));
	    fpol6->SetParameter(4,fpol5->GetParameter(4));
	    fpol6->SetParameter(5,fpol5->GetParameter(5));
	    h3->Fit(fpol6,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
	    if(nDegreeBackPol==6){
	      fback=fpol6;
		if(refitWithSignal)nparback=7;
	    }
	    if(nDegreeBackPol>=7){    
	      fpol7->SetParameter(0,fpol6->GetParameter(0));
	      fpol7->SetParameter(1,fpol6->GetParameter(1));
	      fpol7->SetParameter(2,fpol6->GetParameter(2));
	      fpol7->SetParameter(3,fpol6->GetParameter(3));
	      fpol7->SetParameter(4,fpol6->GetParameter(4));
	      fpol7->SetParameter(5,fpol6->GetParameter(5));
	      fpol7->SetParameter(6,fpol6->GetParameter(6));
	      h3->Fit(fpol7,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
	      if(nDegreeBackPol==7){
		fback=fpol7;
		if(refitWithSignal)nparback=8;
	      }
	      if(nDegreeBackPol>=8){    
		fpol8->SetParameter(0,fpol7->GetParameter(0));
		fpol8->SetParameter(1,fpol7->GetParameter(1));
		fpol8->SetParameter(2,fpol7->GetParameter(2));
		fpol8->SetParameter(3,fpol7->GetParameter(3));
		fpol8->SetParameter(4,fpol7->GetParameter(4));
		fpol8->SetParameter(5,fpol7->GetParameter(5));
		fpol8->SetParameter(6,fpol7->GetParameter(6));
		fpol8->SetParameter(7,fpol7->GetParameter(7));
		h3->Fit(fpol8,"REMIN0","",fitrangelow[binNumber],fitrangeup[binNumber]);
		
		if(nDegreeBackPol==8){
		  fback=fpol8;
		  if(refitWithSignal)nparback=9;
		}
		else{
		  Printf("Too high pol order for back fit (%d): not yet implemented",nDegreeBackPol);
		  delete h2;
		  delete h3;	       
		  return 0x0;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  for(Int_t k=1;k<=h2->GetNbinsX();k++){
    if(!refitWithSignal)h2->SetBinContent(k,h2->GetBinContent(k)-fback->Eval(h2->GetBinCenter(k)));
  }
  delete h3;
  return h2;
}




Bool_t DirectFit(TH1D *hUnSubtracted, Int_t iptbin, TH1F* hRawYieldSBfit){
  
  TH1D *hSubtractedSB=FitBackgroundAndSubtract(hUnSubtracted,iptbin);
  TH1F *hsubTemp;
  TH1F *hsubTempAllRange;
  if(nPtBins>1)cDataSubtractedFit->cd(iptbin+1);
  else cDataSubtractedFit->cd();
  hSubtractedSB->Draw();
  TF1 *fgl=0x0;
  if(refitWithSignal){
    fgl=GausPlusLine(fitrangelow[iptbin],fitrangeup[iptbin]);
    fgl->SetParameter(3,0);
    fgl->SetParameter(4,0);
    fgl->SetName(Form("fgl%d",iptbin));
    fFitSB[iptbin]=new TF1(Form("fFitSB%d",iptbin),Form("%s+%s",fback->GetName(),fgl->GetName()),fitrangelow[iptbin],fitrangeup[iptbin]);
    if(fixBackParFromFirstFit)for(Int_t j=0;j<9;j++)fFitSB[iptbin]->FixParameter(j,fback->GetParameter(j));
    fFitSB[iptbin]->SetParLimits(1+nparback,1.83,1.89);
    
    
    // fast bin counting to initialize fit param
    Int_t binlbc=hSubtractedSB->FindBin(massD-3.*sigmas[iptbin]);
    Int_t binrbc=hSubtractedSB->FindBin(massD+3.*sigmas[iptbin]);
    Double_t backest=fback->Integral(hSubtractedSB->GetBinLowEdge(binlbc),hSubtractedSB->GetBinLowEdge(binrbc+1))/hSubtractedSB->GetBinWidth(1);
    Double_t fastbcestimate=(hSubtractedSB->Integral(binlbc,binrbc)-backest)*hSubtractedSB->GetBinWidth(1);
    fFitSB[iptbin]->SetParameter(nparback,fastbcestimate);
    fFitSB[iptbin]->SetParLimits(nparback,0.1*fastbcestimate,fastbcestimate*10.);    
    if(fixSigma)fFitSB[iptbin]->FixParameter(2+nparback,sigmas[iptbin]);
    else fFitSB[iptbin]->SetParameter(2+nparback,sigmas[iptbin]);
    if(fixMean)fFitSB[iptbin]->FixParameter(1+nparback,massD);
    else fFitSB[iptbin]->SetParameter(1+nparback,massD);    
    fFitSB[iptbin]->FixParameter(3+nparback,0);// These are the linear terms in GausPlusLine
    fFitSB[iptbin]->FixParameter(4+nparback,0);
  }
  else{
    fFitSB[iptbin]=GausPlusLine(fitrangelow[iptbin],fitrangeup[iptbin]);
    fFitSB[iptbin]->SetName(Form("fFitSB%d",iptbin));
    if(fixSigma)fFitSB[iptbin]->FixParameter(2,sigmas[iptbin]);
    else fFitSB[iptbin]->SetParameter(2,sigmas[iptbin]);
    if(fixMean)fFitSB[iptbin]->FixParameter(1+nparback,massD);
    else fFitSB[iptbin]->SetParameter(1+nparback,massD);
  }
  fFitSB[iptbin]->SetRange(fitrangelow[iptbin],fitrangeup[iptbin]);
  if(refitWithSignal){
    hSubtractedSB->Fit(fFitSB[iptbin],"RLEM","");
    for(Int_t j=0;j<nparback;j++){
      fback->SetParameter(j,fFitSB[iptbin]->GetParameter(j));
      fback->SetParError(j,fFitSB[iptbin]->GetParError(j));
      fback->SetLineColor(kGreen);
      fback->SetLineStyle(2);
      //	  hSubtractedSB->GetListOfFunctions()->Add(fFitSB[iptbin]);
    }
    fback->Draw("same");
    fFitSB[iptbin]->Draw("same");
    if(nPtBins>1)cDataSubtractedReFit->cd(iptbin+1);
    else cDataSubtractedReFit->cd();
    hsubTemp=(TH1F*)hSubtractedSB->Clone(Form("%sSubBack",hSubtractedSB->GetName()));
    hsubTempAllRange=(TH1F*)hSubtractedSB->Clone(Form("%sSubBackAllRange",hSubtractedSB->GetName()));
    for(Int_t jst=1;jst<=hsubTemp->GetNbinsX();jst++){
      hsubTemp->SetBinContent(jst,hsubTemp->GetBinContent(jst)-fback->Integral(hsubTemp->GetBinLowEdge(jst),hsubTemp->GetBinLowEdge(jst)+hsubTemp->GetBinWidth(jst))/hsubTemp->GetBinWidth(jst));
      hsubTempAllRange->SetBinContent(jst,hsubTempAllRange->GetBinContent(jst)-fFitSB[iptbin]->Integral(hsubTemp->GetBinLowEdge(jst),hsubTempAllRange->GetBinLowEdge(jst)+hsubTempAllRange->GetBinWidth(jst))/hsubTempAllRange->GetBinWidth(jst));
    }
    hsubTemp->SetLineColor(kBlue);
    hsubTempAllRange->SetLineColor(kGray+2);
    // check better range for graphical purposes
    hsubTemp->GetXaxis()->SetRangeUser(fitrangelow[iptbin]-0.1*(fitrangeup[iptbin]-fitrangelow[iptbin]),fitrangeup[iptbin]+0.1*(fitrangeup[iptbin]-fitrangelow[iptbin]));
    Double_t minst,maxst;
    Int_t binst=hsubTemp->FindBin(fitrangelow[iptbin]-0.1*(fitrangeup[iptbin]-fitrangelow[iptbin]));
    Int_t binstm=hsubTemp->FindBin(fitrangeup[iptbin]+0.1*(fitrangeup[iptbin]-fitrangelow[iptbin]));
    minst=hsubTemp->GetBinContent(binst)-hsubTemp->GetBinError(binst);
    maxst=hsubTemp->GetBinContent(binstm)+hsubTemp->GetBinError(binstm);
    for(Int_t stbin=binst;stbin<=binstm;stbin++){
      Double_t valb=hsubTemp->GetBinContent(stbin)-hsubTemp->GetBinError(stbin);
      if(valb<minst)minst=valb;
      valb=hsubTemp->GetBinContent(stbin)+hsubTemp->GetBinError(stbin);
      if(valb>maxst)maxst=valb;
    }
    if(minst<0)minst*=1.2;
    else minst*=0.8;
    if(maxst<0)maxst*=0.8;
    else maxst*=1.2;
    hsubTemp->GetYaxis()->SetRangeUser(minst,maxst);
    hsubTemp->Draw();
    hsubTempAllRange->Draw("same");
    hsubTemp->Draw("same");
    fgl->SetParameter(0,fFitSB[iptbin]->GetParameter(nparback));
    fgl->SetParameter(1,fFitSB[iptbin]->GetParameter(nparback+1));
    fgl->SetParameter(2,fFitSB[iptbin]->GetParameter(nparback+2));
    fgl->Draw("same");
  
    if(nPtBins>1)cDataSubtractedFit->cd(iptbin+1);
    else cDataSubtractedFit->cd();
  }
  else hSubtractedSB->Fit(fFitSB[iptbin],"REMI","");
  
  hSubtractedSB->SetStats(kFALSE);
  hSubtractedSB->GetXaxis()->SetRangeUser(fitrangelow[iptbin],fitrangeup[iptbin]);
  
  // FITTING PART FINISHED
  // NOW DRAW PARAMETERS AND EXTRACT RAW YIELDS
  Double_t signif,signal,errsignal,total,sigma,mean;
  mean=fFitSB[iptbin]->GetParameter(1+nparback);
  sigma=fFitSB[iptbin]->GetParameter(2+nparback);
  TPaveText *pinfo=new TPaveText(0.5,0.65,0.9,0.9,"NDC");
  pinfo->SetBorderSize(0);
  pinfo->SetFillStyle(0);  
  pinfo->AddText(Form("#mu = %.3f#pm%.3f",mean,fFitSB[iptbin]->GetParError(1+nparback)));
  pinfo->AddText(Form("#sigma = %.3f#pm%.3f",sigma,fFitSB[iptbin]->GetParError(2+nparback)));
  pinfo->AddText(Form("#chi^{2}/NDOF = %.1f/%d",fFitSB[iptbin]->GetChisquare(),fFitSB[iptbin]->GetNDF()));
  pinfo->AddText(Form("prob =%.3f",fFitSB[iptbin]->GetProb()));
  pinfo->SetTextColor(4);
  pinfo->Draw();
  Int_t binMinTot=hUnSubtracted->FindBin((mean-3.*sigma)*1.0001);
  Int_t binMaxTot=hUnSubtracted->FindBin((mean+3.*sigma)*0.9999);
  total= hUnSubtracted->Integral(binMinTot,binMaxTot);
  signal=fFitSB[iptbin]->GetParameter(0+nparback)/hSubtractedSB->GetBinWidth(1);
  errsignal=fFitSB[iptbin]->GetParError(0+nparback)/hSubtractedSB->GetBinWidth(1);

  //  signalbc=GetSignalBinCounting(hSubtractedSB,fFitSB[iptbin],errsignbc,nsigmaBinCounting);
  
  //  hRawYieldBinCountBackSubFit->SetBinContent(ib+1,signalbc);
  //  hRawYieldBinCountBackSubFit->SetBinError(ib+1,errsignbc);

  //  hRawYieldBackSubFit->SetBinContent(ib+1,signal);
  //  hRawYieldBackSubFit->SetBinError(ib+1,errsignal);
  Double_t signifT=signal/TMath::Sqrt(total);
  signif=signifT;
  Double_t backfit=-1;
  if(refitWithSignal){
    backfit=fback->Integral(mean-3.*sigma,mean+3.*sigma)/hSubtractedSB->GetBinWidth(1);
    signif=signal/TMath::Sqrt(signal+backfit);
  }
  hRawYieldSBfit->SetBinContent(iptbin+1,signal);
  hRawYieldSBfit->SetBinError(iptbin+1,errsignal);
  
  // background values
  //  binMinTot=hUnSubtracted->FindBin((massD-3.*sigmas[iptbin])*1.00001);// the 1.00001 is useless...
  //  binMaxTot=hUnSubtracted->FindBin((massD+3.*sigmas[iptbin])*0.9999);// same for hte 0.9999
  //    Printf("Min bin and max bin are: ");
  //hTotYieldSignalRegion->SetBinContent(ib+1,hUnSubtracted->Integral(binMinTot,binMaxTot));// this is not backfit, I simply do not need more variables
  //  if(hTotYieldSignalRegion->GetBinContent(ib+1)>0.){
  //    hTotYieldSignalRegion->SetBinError(ib+1,TMath::Sqrt(hTotYieldSignalRegion->GetBinContent(ib+1)));
  //  }
  //  hBackBCSignalRegion->SetBinContent(ib+1,total-signal);
  //    if(total-signal>0.){
  //      hBackBCSignalRegion->SetBinError(ib+1,TMath::Sqrt(hBackBCSignalRegion->GetBinContent(ib+1)));
  //    }
  //    hBackfitSignalRegion->SetBinContent(ib+1,backfit);
  //    if(backfit>0.){
  //      hBackfitSignalRegion->SetBinError(ib+1,TMath::Sqrt(backfit));
  //    }

  TPaveText *pinfo2=new TPaveText(0.1,0.1,0.5,0.4,"NDC");
  pinfo2->SetBorderSize(0);
  pinfo2->SetFillStyle(0);
  TString str=Form("S(%d#sigma)= %.0f #pm %.0f ",3,signal,errsignal);
  pinfo2->AddText(str);
  str=Form("S/#sigma(S) = %.1f ",signal/errsignal);
  pinfo2->AddText(str);
  str=Form("S/#sqrt{Tot}= %.1f ",signifT);
  pinfo2->AddText(str);
  str=Form("signif(fit)= %.1f ",signif);
  pinfo2->AddText(str);
  if(backfit>1.e-9){
    str=Form("S/B(3#sigma)= %.1f#times 10^{-3} ",signal/backfit*1.e+3);
  }
  else {
    str=Form("S/B(3#sigma)= %.1f#times 10^{-3} ",signal/(total-signal)*1.e+3);
  }
  pinfo2->AddText(str);
  pinfo2->SetTextColor(2);
  pinfo2->Draw();
        
  cDataFitRes[iptbin]->cd(1);
  hSubtractedSB->Draw();
  hSubtractedSB->SetStats(kFALSE);
  TList *lf=hSubtractedSB->GetListOfFunctions();
  if(refitWithSignal){
    lf->Add(fgl);
    fback->Draw("same");
  }
  fFitSB[iptbin]->Draw("same");
  pinfo->Draw("same");
  pinfo2->Draw("same");
  cDataFitRes[iptbin]->cd(2);
  if(refitWithSignal){  
    hsubTemp->Draw("same");
    
    TLine *tlin=new TLine(hsubTemp->GetBinLowEdge(1),0.,fitrangeup[iptbin]+0.1*(fitrangeup[iptbin]-fitrangelow[iptbin]),0.);
    tlin->SetLineColor(kBlue);
    tlin->SetLineStyle(2);
    tlin->Draw("same");
    hsubTemp->SetMarkerColor(kBlack);
    hsubTemp->SetMarkerStyle(20);
    hsubTemp->SetMarkerSize(0.7);
    hsubTemp->SetLineColor(kBlack);
    hsubTemp->SetStats(kFALSE);
    hsubTempAllRange->SetStats(kFALSE);
    hsubTempAllRange->Draw("same");
    hsubTemp->Draw("same");
    fgl->Draw("same");
  }
  


  return kTRUE;
}



void ProjectCombinHFAndFit(){

  
  
  TString dirName=Form("PWG3_D2H_InvMass%sLowPt%s",meson.Data(),suffix.Data());
  TString lstName=Form("coutput%s%s",meson.Data(),suffix.Data());

  if(meson=="Dplus") massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  else if(meson=="Dzero") massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  TFile* fil=new TFile(fileName.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  
  AliNormalizationCounter *nC=(AliNormalizationCounter*)df->Get("NormalizationCounter");
  TH1F* hnEv=new TH1F("hEvForNorm","events for normalization",1,0,1);
  printf("Number of Ev. for norm=%d\n",(Int_t)nC->GetNEventsForNorm());
  hnEv->SetBinContent(1,nC->GetNEventsForNorm());

  TList* l=(TList*)df->Get(lstName.Data());
  l->ls();

  TH3F* h3d=(TH3F*)l->FindObject("hMassVsPtVsY");
  TH3F* h3dr=(TH3F*)l->FindObject("hMassVsPtVsYRot");
  TH3F* h3dme=(TH3F*)l->FindObject("hMassVsPtVsYME");
  TH3F* h3dmepp=(TH3F*)l->FindObject("hMassVsPtVsYMELSpp");
  TH3F* h3dmemm=(TH3F*)l->FindObject("hMassVsPtVsYMELSmm");
  TH2F* hpoolMix=(TH2F*)l->FindObject("hMixingsPerPool");
  TH2F* hpoolEv=(TH2F*)l->FindObject("hEventsPerPool");
  TH3F* h3dlsp=(TH3F*)l->FindObject("hMassVsPtVsYLSpp");
  TH3F* h3dlsm=(TH3F*)l->FindObject("hMassVsPtVsYLSmm");

  // Side-band fit method initialization
  if(tryDirectFit){
    fback=0x0;
    fFitSB=new TF1*[nPtBins];
    cDataSubtractedFit=new TCanvas("cDataSubtractedFit","cDataSubtractedFit",800,800);
    cDataSubtractedReFit=new TCanvas("cDataSubtractedReFit","cDataSubtractedReFit",800,800);
    cPullsDirect=new TCanvas("cPullsDirect","cPullsDirect",800,800);
    cPullsDirectTrendVsMass=new TCanvas("cPullsDirectTrendVsMass","cPullsDirectTrendVsMass",800,800);
    cResidualsDirect=new TCanvas("cResidualsDirect","cResidualsDirect",800,800);
    cResidualDirectTrendVsMass=new TCanvas("cResidualDirectTrendVsMass","cResidualDirectTrendVsMass",800,800);

    DivideCanvas(cDataSubtractedFit,nPtBins);
    DivideCanvas(cDataSubtractedReFit,nPtBins);
    DivideCanvas(cPullsDirect,nPtBins);
    DivideCanvas(cPullsDirectTrendVsMass,nPtBins);
    DivideCanvas(cResidualsDirect,nPtBins);
    DivideCanvas(cResidualDirectTrendVsMass,nPtBins);
    cDataFitRes=new TCanvas*[nPtBins];
    for(Int_t ib=0;ib<nPtBins;ib++){
      cDataFitRes[ib]=new TCanvas(Form("cDataSBFitResBin%d",ib),Form("cDataSBFitResBin%d",ib),800,800);
      cDataFitRes[ib]->Divide(1,2);
    }
    
  }

  TCanvas* cem=new TCanvas("cem","Pools",1200,600);
  cem->Divide(2,1);
  cem->cd(1);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  hpoolMix->GetXaxis()->SetTitle("zVertex");
  hpoolMix->GetYaxis()->SetTitle("Multiplicity");
  hpoolMix->GetYaxis()->SetTitleOffset(1.4);
  hpoolMix->Draw("colztext");
  cem->cd(2);
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.12);
  hpoolEv->Draw("colztext");
  hpoolEv->GetXaxis()->SetTitle("zVertex");
  hpoolEv->GetYaxis()->SetTitle("Multiplicity");
  hpoolEv->GetYaxis()->SetTitleOffset(1.4);


  TCanvas* c1=new TCanvas("c1","Mass",1200,800);
  c1->Divide(3,2);

  TCanvas* c2=new TCanvas("c2","Mass-Bkg Rot",1200,800);
  c2->Divide(3,2);
  TCanvas* c2pulls=new TCanvas("c2pulls","Mass-Bkg Rot pulls",1200,800);
  c2pulls->Divide(3,2);
  TCanvas* c2residuals=new TCanvas("c2residuals","Mass-Bkg Rot residuals",1200,800);
  c2residuals->Divide(3,2);
  TCanvas* c2residualTrend=new TCanvas("c2residualTrend","Mass-Bkg Rot residual trend vs. mass",1200,800);
  c2residualTrend->Divide(3,2);
  TCanvas* c2pullTrend=new TCanvas("c2pullTrend","Mass-Bkg Rot pull trend vs. mass",1200,800);
  c2pullTrend->Divide(3,2);

  TCanvas* c3=new TCanvas("c3","Mass-Bkg LS",1200,800);
  c3->Divide(3,2);
  TCanvas* c3pulls=new TCanvas("c3pulls","Mass-Bkg LS pulls",1200,800);
  c3pulls->Divide(3,2);
  TCanvas* c3residuals=new TCanvas("c3residuals","Mass-Bkg LS residuals",1200,800);
  c3residuals->Divide(3,2);
  TCanvas* c3residualTrend=new TCanvas("c3residualTrend","Mass-Bkg LS residual trend vs. mass",1200,800);
  c3residualTrend->Divide(3,2);
  TCanvas* c3pullTrend=new TCanvas("c3pullTrend","Mass-Bkg LS pull trend vs. mass",1200,800);
  c3pullTrend->Divide(3,2);

  TCanvas* c4=new TCanvas("c4","Mass-Bkg EM",1200,800);
  c4->Divide(3,2);
  TCanvas* c4pulls=new TCanvas("c4pulls","Mass-Bkg EM pulls",1200,800);
  c4pulls->Divide(3,2);
  TCanvas* c4residuals=new TCanvas("c4residuals","Mass-Bkg EM residuals",1200,800);
  c4residuals->Divide(3,2);
  TCanvas* c4residualTrend=new TCanvas("c4residualTrend","Mass-Bkg EM residual trend vs. mass",1200,800);
  c4residualTrend->Divide(3,2);
  TCanvas* c4pullTrend=new TCanvas("c4pullTrend","Mass-Bkg EM pull trend vs. mass",1200,800);
  c4pullTrend->Divide(3,2);

  TCanvas *cCompareResidualTrends=new TCanvas("cCompareResidualTrends","cCompareResidualTrends",1200,800);
  DivideCanvas(cCompareResidualTrends,nPtBins);  


  AliHFMassFitter* fitterRot[nPtBins];
  AliHFMassFitter* fitterLS[nPtBins];
  AliHFMassFitter* fitterME[nPtBins];

  TH1F* hRawYieldRot=new TH1F("hRawYieldRot","",nPtBins,binLims);
  TH1F* hRawYieldLS=new TH1F("hRawYieldLS","",nPtBins,binLims);
  TH1F* hRawYieldME=new TH1F("hRawYieldME","",nPtBins,binLims);
  TH1F* hRawYieldSBfit=new TH1F("hRawYieldSBfit","",nPtBins,binLims);

  TH1F* hRawYieldRotBC=new TH1F("hRawYieldRotBC","BC yield (rotational background)",nPtBins,binLims);
  TH1F* hRawYieldLSBC=new TH1F("hRawYieldLSBC","BC yield (like-sign background)",nPtBins,binLims);
  TH1F* hRawYieldMEBC=new TH1F("hRawYieldMEBC","BC yield (mixed-event background)",nPtBins,binLims);
  TH1F* hRawYieldSBfitBC=new TH1F("hRawYieldSBfitBC","BC yield (direct fit background)",nPtBins,binLims);


  TLatex* tME=new TLatex(0.65,0.82,"MixEv +- pairs");
  tME->SetNDC();
  TLatex* tMEpp=new TLatex(0.65,0.75,"MixEv ++ pairs");
  tMEpp->SetNDC();
  TLatex* tMEmm=new TLatex(0.65,0.68,"MixEv -- pairs");
  tMEmm->SetNDC();

  TDirectory *current = gDirectory;
  TFile* fout=new TFile(Form("outputMassFits_%s.root",suffix.Data()),"recreate");
  current->cd();

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){

    Int_t bin1=h3d->GetYaxis()->FindBin(binLims[iPtBin]);
    Int_t bin2=h3d->GetYaxis()->FindBin(binLims[iPtBin+1]-0.0001);
    printf("Bin %d   Pt range=%f %f\n",iPtBin,h3d->GetYaxis()->GetBinLowEdge(bin1),h3d->GetYaxis()->GetBinUpEdge(bin2));
    TH1D* hMassPtBin=h3d->ProjectionX(Form("hMassPtBin%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinr=h3dr->ProjectionX(Form("hMassPtBinr%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinme=h3dme->ProjectionX(Form("hMassPtBinme%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinmeLSpp=h3dmepp->ProjectionX(Form("hMassPtBinmeLSpp%d",iPtBin),bin1,bin2);
    TH1D* hMassPtBinmeLSmm=h3dmemm->ProjectionX(Form("hMassPtBinmeLSmm%d",iPtBin),bin1,bin2);

    TH1D* hMassPtBinlsp=0x0;
    TH1D* hMassPtBinlsm=0x0;
    TH1D* hMassPtBinls=0x0;
    if(h3dlsp){ 
      hMassPtBinlsp=h3dlsp->ProjectionX(Form("hMassPtBinlsp%d",iPtBin),bin1,bin2);
      hMassPtBinlsm=h3dlsm->ProjectionX(Form("hMassPtBinlsm%d",iPtBin),bin1,bin2);
      hMassPtBinls=(TH1D*)hMassPtBinlsp->Clone(Form("hMassPtBinls%d",iPtBin));
      hMassPtBinls->Reset();
      for(Int_t iBin=1; iBin<=hMassPtBinlsp->GetNbinsX(); iBin++){
	Double_t np=hMassPtBinlsp->GetBinContent(iBin);
	Double_t nm=hMassPtBinlsm->GetBinContent(iBin);
	Double_t tt=2*TMath::Sqrt(np*nm);
	Double_t enp=hMassPtBinlsp->GetBinError(iBin);
	Double_t enm=hMassPtBinlsm->GetBinError(iBin);
	Double_t ett=0;
	if(tt>0) ett=2./tt*TMath::Sqrt(np*np*enm*enm+nm*nm*enp*enp);
	hMassPtBinls->SetBinContent(iBin,tt);
	hMassPtBinls->SetBinError(iBin,ett);
      }
      hMassPtBinls->SetLineColor(kGreen+1);
    }

    // hMassPtBin->Sumw2();
    // hMassPtBinr->Sumw2();
    // hMassPtBinme->Sumw2();
    // hMassPtBinlsp->Sumw2();
    // hMassPtBinls->Sumw2();
    hMassPtBin->SetLineColor(1);
    hMassPtBinr->SetLineColor(4);
    hMassPtBinme->SetLineColor(kOrange+1);
    hMassPtBinmeLSpp->SetLineColor(kGreen+2);
    hMassPtBinmeLSmm->SetLineColor(kCyan);
    hMassPtBin->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    hMassPtBinr->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    hMassPtBinme->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hRatioMEpp=(TH1D*)hMassPtBinmeLSpp->Clone(Form("hRatioPtBinmeLSpp%d",iPtBin));
    hRatioMEpp->Divide(hMassPtBinme);
    TH1D* hRatioMEmm=(TH1D*)hMassPtBinmeLSmm->Clone(Form("hRatioPtBinmeLSmm%d",iPtBin));
    hRatioMEmm->Divide(hMassPtBinme);
    TH1D* hMassPtBinmeAll=(TH1D*)hMassPtBinme->Clone(Form("hRatioPtBinmeAll%d",iPtBin));
    hMassPtBinmeAll->Add(hMassPtBinmeLSpp);
    hMassPtBinmeAll->Add(hMassPtBinmeLSmm);
    hMassPtBinmeAll->SetLineColor(kRed);
    TH1D* hRatioME=(TH1D*)hMassPtBinme->Clone(Form("hRatioME%d",iPtBin));
    hRatioME->Divide(hMassPtBin);
    TH1D* hRatioMEAll=(TH1D*)hMassPtBinmeAll->Clone(Form("hRatioMEAll%d",iPtBin));
    hRatioMEAll->Divide(hMassPtBin);


    TCanvas* c0=new TCanvas(Form("CBin%d",iPtBin),Form("Bin%d norm",iPtBin),1000,700);
    c0->Divide(2,2);
    c0->cd(1);
    TH1D* hRatio=(TH1D*)hMassPtBinr->Clone("hRatio");
    hRatio->Divide(hMassPtBin);
    hRatio->Draw();
    hRatio->GetYaxis()->SetTitle("Rotational/All");
    hRatio->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    Double_t norm=hRatio->GetMaximum();
    if(optForNorm==1){
      norm=0.0001;
      for(Int_t iMassBin=1; iMassBin<hRatio->GetNbinsX(); iMassBin++){
	Double_t bce=hRatio->GetBinCenter(iMassBin);
	if(bce>minMass && bce<maxMass){
	  Double_t bco=hRatio->GetBinContent(iMassBin);
	  if(bco>norm) norm=bco;
	}
      }
    }else if(optForNorm==2){ 
      hRatio->Fit("pol0");
      TF1* func0=(TF1*)hRatio->GetListOfFunctions()->FindObject("pol0");
      norm=func0->GetParameter(0);
    }
    hMassPtBinr->Scale(1./norm);
    c0->cd(2);
    hMassPtBinme->GetYaxis()->SetTitle("Entries (EvMix)");
    hMassPtBinme->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassPtBinme->DrawCopy();
    hMassPtBinmeLSpp->Draw("same");
    hMassPtBinmeLSmm->Draw("same");
    tME->SetTextColor(hMassPtBinme->GetLineColor());
    tMEpp->SetTextColor(hMassPtBinmeLSpp->GetLineColor());
    tMEmm->SetTextColor(hMassPtBinmeLSmm->GetLineColor());
    tME->Draw();
    tMEpp->Draw();
    tMEmm->Draw();
    c0->cd(4);
    hRatioMEpp->Draw();
    hRatioMEpp->SetMinimum(0.4);
    hRatioMEpp->SetMaximum(0.6);
    hRatioMEpp->GetYaxis()->SetTitle("ME with LS / ME with OS");
    hRatioMEpp->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hRatioMEmm->Draw("same");
    c0->cd(3);
    hRatioME->Draw();
    hRatioME->SetMaximum(hRatioMEAll->GetMaximum()*1.05);
    hRatioMEAll->Draw("same");
    hRatioME->GetYaxis()->SetTitle("EvMix/All");
    hRatioME->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
 

    Double_t normME=hRatioME->GetMaximum();
    Double_t normMEAll=hRatioMEAll->GetMaximum();
    if(optForNorm==1){
      normME=0.0001;
      for(Int_t iMassBin=1; iMassBin<hRatioME->GetNbinsX(); iMassBin++){
	Double_t bce=hRatioME->GetBinCenter(iMassBin);
	if(bce>minMass && bce<maxMass){
	  Double_t bco=hRatioME->GetBinContent(iMassBin);
	  if(bco>normME) normME=bco;
	}
      }
      normMEAll=0.0001;
      for(Int_t iMassBin=1; iMassBin<hRatioME->GetNbinsX(); iMassBin++){
	Double_t bce=hRatioMEAll->GetBinCenter(iMassBin);
	if(bce>minMass && bce<maxMass){
	  Double_t bco=hRatioMEAll->GetBinContent(iMassBin);
	  if(bco>normMEAll) normMEAll=bco;
	}
      }
    }else if(optForNorm==2){ 
      hRatioME->Fit("pol0");
      TF1* func0ME=(TF1*)hRatioME->GetListOfFunctions()->FindObject("pol0");
      normME=func0ME->GetParameter(0);
      hRatioMEAll->Fit("pol0");
      TF1* func0MEAll=(TF1*)hRatioMEAll->GetListOfFunctions()->FindObject("pol0");
      normMEAll=func0MEAll->GetParameter(0);
    }
    hMassPtBinme->Scale(1./normME);
    hMassPtBinmeAll->Scale(1./normMEAll);



    c1->cd(iPtBin+1);
    hMassPtBin->Draw();
    hMassPtBin->GetYaxis()->SetTitle("Counts");
    hMassPtBin->GetYaxis()->SetTitleOffset(2.);
    hMassPtBin->GetXaxis()->SetTitle("Invariant mass (GeV/c^{2})");
    hMassPtBinr->Draw("same");
    if(!useEMwithLS) hMassPtBinme->Draw("same");
    else  hMassPtBinmeAll->Draw("same");
    if(hMassPtBinls) hMassPtBinls->Draw("same");
    if(iPtBin==0){
      TLegend* leg=new TLegend(0.5,0.6,0.89,0.89);
      leg->SetFillStyle(0);
      leg->AddEntry(hMassPtBin,"All candidates","L")->SetTextColor(hMassPtBin->GetLineColor());
      leg->AddEntry(hMassPtBinr,"Background (rotations)","L")->SetTextColor(hMassPtBinr->GetLineColor());
      if(!useEMwithLS) leg->AddEntry(hMassPtBinme,"Background (ME)","L")->SetTextColor(hMassPtBinme->GetLineColor());
      else leg->AddEntry(hMassPtBinmeAll,"Background (ME)","L")->SetTextColor(hMassPtBinmeAll->GetLineColor());
      if(hMassPtBinls) leg->AddEntry(hMassPtBinls,"Like-sign","L")->SetTextColor(hMassPtBinls->GetLineColor());
      leg->Draw();
    }
  
    TH1D* hMassSubRot=(TH1D*)hMassPtBin->Clone(Form("hMassSubRot_bin%d",iPtBin));
    hMassSubRot->Add(hMassPtBinr,-1);
    hMassSubRot->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c -- Rotational",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hMassSubME=(TH1D*)hMassPtBin->Clone(Form("hMassSubME_bin%d",iPtBin));
    if(!useEMwithLS) hMassSubME->Add(hMassPtBinme,-1);
    else hMassSubME->Add(hMassPtBinmeAll,-1);
    hMassSubME->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c -- Mixed Ev",binLims[iPtBin],binLims[iPtBin+1]));
    TH1D* hMassSubLS=0x0;
    if(hMassPtBinls){
      hMassSubLS=(TH1D*)hMassPtBin->Clone(Form("hMassSubLS_bin%d",iPtBin));
      hMassSubLS->Add(hMassPtBinls,-1);
      hMassSubLS->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c  -- Like Sign",binLims[iPtBin],binLims[iPtBin+1]));
    }
 
    fout->cd();
    hMassPtBin->Write();
    hMassSubRot->Write();
    hMassSubME->Write();
    if(hMassPtBinls) hMassSubLS->Write();
    current->cd();

    hMassSubRot->Rebin(rebin);
    hMassSubME->Rebin(rebin);
    if(hMassPtBinls) hMassSubLS->Rebin(rebin);
 
    fitterRot[iPtBin]=ConfigureFitter(hMassSubRot,iPtBin);
    if(hMassPtBinls) fitterLS[iPtBin]=ConfigureFitter(hMassSubLS,iPtBin);
    fitterME[iPtBin]=ConfigureFitter(hMassSubME,iPtBin);

    Bool_t out1=fitterRot[iPtBin]->MassFitter(0);
    Bool_t out2=kFALSE;
    if(hMassPtBinls) out2=fitterLS[iPtBin]->MassFitter(0);
    Bool_t out3=fitterME[iPtBin]->MassFitter(0);
    
    Bool_t out4=kFALSE;
    if(tryDirectFit){
      TH1D *hMassDirectFit=(TH1D*)hMassPtBin->Clone(Form("hMassDirectFit_bin%d",iPtBin));
      hMassDirectFit->Rebin(rebin);
      out4=DirectFit(hMassDirectFit,iPtBin,hRawYieldSBfit);
      if(out4){
	Double_t errbc;
	Double_t bc=GetSignalBinCounting(hMassDirectFit,fback,errbc,nsigmaBinCounting,fFitSB[iPtBin]->GetParameter(2+nparback),-1.*fitrangelow[iPtBin],-1.*fitrangeup[iPtBin]);
	hRawYieldSBfitBC->SetBinContent(iPtBin+1,bc);
	hRawYieldSBfitBC->SetBinError(iPtBin+1,errbc);
      

	cPullsDirect->cd(iPtBin+1);
	TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
	TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
	TH1F *hResiduals=(TH1F*)GetResidualsAndPulls(hMassDirectFit,fFitSB[iPtBin],fitrangelow[iPtBin],fitrangeup[iPtBin],hPulls,hResidualTrend,hPullsTrend);

	hPulls->Draw();
	
	cResidualsDirect->cd(iPtBin+1);
	hResiduals->Draw();
	
	SetStyleHisto(hResidualTrend,kSB);
	cCompareResidualTrends->cd(iPtBin+1);
	hResidualTrend->Draw();

	cResidualDirectTrendVsMass->cd(iPtBin+1);
	hResidualTrend->Draw();
	
	cPullsDirectTrendVsMass->cd(iPtBin+1);
	hPullsTrend->Draw();
	//TPaveStats *tp=(TPaveStats*)hPulls->FindObject("stats");
	//      tp->SetOptStat("remnpcev");
	delete hMassDirectFit;
      }
    }
    
    c2->cd(iPtBin+1);
    if(out1){
      fitterRot[iPtBin]->DrawHere(gPad,3,0);
      hRawYieldRot->SetBinContent(iPtBin+1,fitterRot[iPtBin]->GetRawYield());
      hRawYieldRot->SetBinError(iPtBin+1,fitterRot[iPtBin]->GetRawYieldError());
      WriteFitInfo(fitterRot[iPtBin],hMassPtBin);
      
      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubRot,fitterRot[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterRot[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldRotBC->SetBinContent(iPtBin+1,bc);
      hRawYieldRotBC->SetBinError(iPtBin+1,errbc);

      c2pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=(TH1F*)GetResidualsAndPulls(hMassSubRot,fitterRot[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);

      hPulls->Draw();

      c2residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kRot);
      c2residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();
      cCompareResidualTrends->cd(iPtBin+1);
      if(out4)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c2pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();

    }
    else hMassSubRot->Draw("");
    

    c3->cd(iPtBin+1);
    if(out2){
      fitterLS[iPtBin]->DrawHere(gPad,3,0); 
      hRawYieldLS->SetBinContent(iPtBin+1,fitterLS[iPtBin]->GetRawYield());
      hRawYieldLS->SetBinError(iPtBin+1,fitterLS[iPtBin]->GetRawYieldError());
      WriteFitInfo(fitterLS[iPtBin],hMassPtBin);

      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubLS,fitterLS[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterLS[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldLSBC->SetBinContent(iPtBin+1,bc);
      hRawYieldLSBC->SetBinError(iPtBin+1,errbc);

      c3pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=(TH1F*)GetResidualsAndPulls(hMassSubLS,fitterLS[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);
  
      hPulls->Draw();

      c3residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kLS);
      c3residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();
      
      cCompareResidualTrends->cd(iPtBin+1);
      if(out4||out1)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c3pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();
    }
    else if(hMassPtBinls) hMassSubLS->Draw("");

    c4->cd(iPtBin+1);
    if(out3){ 
      fitterME[iPtBin]->DrawHere(gPad,3,0); 
      hRawYieldME->SetBinContent(iPtBin+1,fitterME[iPtBin]->GetRawYield());
      hRawYieldME->SetBinError(iPtBin+1,fitterME[iPtBin]->GetRawYieldError());
      WriteFitInfo(fitterME[iPtBin],hMassPtBin);

      Double_t errbc;
      Double_t bc=GetSignalBinCounting(hMassSubME,fitterME[iPtBin]->GetBackgroundRecalcFunc(),errbc,nsigmaBinCounting,fitterME[iPtBin]->GetSigma(),-1.*minMass,-1.*maxMass);
      hRawYieldMEBC->SetBinContent(iPtBin+1,bc);
      hRawYieldMEBC->SetBinError(iPtBin+1,errbc);


      c4pulls->cd(iPtBin+1);
      TH1F *hPullsTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hPulls=new TH1F();//the name is changed afterwards, histo must not be deleted
      TH1F *hResidualTrend=new TH1F();// the name is changed afterwards, histo must not be deleted
      TH1F *hResiduals=(TH1F*)GetResidualsAndPulls(hMassSubME,fitterME[iPtBin]->GetMassFunc(),minMass,maxMass,hPulls,hResidualTrend,hPullsTrend);
      hPulls->Draw();

      c4residuals->cd(iPtBin+1);
      hResiduals->Draw();

      SetStyleHisto(hResidualTrend,kME);
      c4residualTrend->cd(iPtBin+1);
      hResidualTrend->Draw();

      cCompareResidualTrends->cd(iPtBin+1);
      if(out4||out1||out2)hResidualTrend->Draw("same");
      else hResidualTrend->Draw();

      c4pullTrend->cd(iPtBin+1);
      hPullsTrend->Draw();
    }
    else hMassSubME->Draw("");

  }
  TString path(gSystem->pwd());
  path.Append("/figures");
  if(gSystem->AccessPathName(path.Data())){ 
    gROOT->ProcessLine(Form(".!mkdir -p %s",path.Data()));  
  }

  c2->SaveAs(Form("figures/InvMassSpectra_%s_Rot.eps",suffix.Data()));
  c3->SaveAs(Form("figures/InvMassSpectra_%s_LS.eps",suffix.Data()));
  c4->SaveAs(Form("figures/InvMassSpectra_%s_EM.eps",suffix.Data()));

  c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_Rot.eps",suffix.Data()));
  c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_LS.eps",suffix.Data()));
  c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_EM.eps",suffix.Data()));

  c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_Rot.eps",suffix.Data()));
  c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_LS.eps",suffix.Data()));
  c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_EM.eps",suffix.Data()));

  c2pulls->SaveAs(Form("figures/PullDistribution_%s_Rot.eps",suffix.Data()));
  c3pulls->SaveAs(Form("figures/PullDistribution_%s_LS.eps",suffix.Data()));
  c4pulls->SaveAs(Form("figures/PullDistribution_%s_EM.eps",suffix.Data()));

  c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_Rot.eps",suffix.Data()));
  c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_LS.eps",suffix.Data()));
  c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_EM.eps",suffix.Data()));


  // save also .root
  c2->SaveAs(Form("figures/InvMassSpectra_%s_Rot.root",suffix.Data()));
  c3->SaveAs(Form("figures/InvMassSpectra_%s_LS.root",suffix.Data()));
  c4->SaveAs(Form("figures/InvMassSpectra_%s_EM.root",suffix.Data()));

  c2residuals->SaveAs(Form("figures/ResidualDistribution_%s_Rot.root",suffix.Data()));
  c3residuals->SaveAs(Form("figures/ResidualDistribution_%s_LS.root",suffix.Data()));
  c4residuals->SaveAs(Form("figures/ResidualDistribution_%s_EM.root",suffix.Data()));

  c2residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_Rot.root",suffix.Data()));
  c3residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_LS.root",suffix.Data()));
  c4residualTrend->SaveAs(Form("figures/residualTrendvsMass_%s_EM.root",suffix.Data()));

  c2pulls->SaveAs(Form("figures/PullDistribution_%s_Rot.root",suffix.Data()));
  c3pulls->SaveAs(Form("figures/PullDistribution_%s_LS.root",suffix.Data()));
  c4pulls->SaveAs(Form("figures/PullDistribution_%s_EM.root",suffix.Data()));

  c2pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_Rot.root",suffix.Data()));
  c3pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_LS.root",suffix.Data()));
  c4pullTrend->SaveAs(Form("figures/pullTrendvsMass_%s_EM.root",suffix.Data()));


  if(tryDirectFit){
    cDataSubtractedFit->SaveAs(Form("figures/InvMassSpectra_%s_SB.eps",suffix.Data()));
    cResidualsDirect->SaveAs(Form("figures/ResidualDistribution_%s_SB.eps",suffix.Data()));
    cResidualDirectTrendVsMass->SaveAs(Form("figures/residualTrendvsMass_%s_SB.eps",suffix.Data()));
    cPullsDirect->SaveAs(Form("figures/PullDistribution_%s_SB.eps",suffix.Data()));
    cPullsDirectTrendVsMass->SaveAs(Form("figures/pullTrendvsMass_%s_SB.eps",suffix.Data()));

    cDataSubtractedFit->SaveAs(Form("figures/InvMassSpectra_%s_SB.root",suffix.Data()));
    cResidualsDirect->SaveAs(Form("figures/ResidualDistribution_%s_SB.root",suffix.Data()));
    cResidualDirectTrendVsMass->SaveAs(Form("figures/residualTrendvsMass_%s_SB.root",suffix.Data()));
    cPullsDirect->SaveAs(Form("figures/PullDistribution_%s_SB.root",suffix.Data()));
    cPullsDirectTrendVsMass->SaveAs(Form("figures/pullTrendvsMass_%s_SB.root",suffix.Data()));
  }
  TCanvas* cry=new TCanvas("cry","RawYield",800,700);
  cry->SetLeftMargin(0.15);
  hRawYieldRot->SetMarkerStyle(21);
  hRawYieldRot->Draw("P");
  hRawYieldRot->SetMinimum(0);
  hRawYieldRot->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRawYieldRot->GetYaxis()->SetTitle("Raw Yield");
  hRawYieldRot->GetYaxis()->SetTitleOffset(1.8);
  Double_t max=hRawYieldRot->GetMaximum();
  if(hRawYieldLS->GetMaximum()>max)max=hRawYieldLS->GetMaximum();
  if(hRawYieldME->GetMaximum()>max)max=hRawYieldME->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSBfit->GetMaximum()>max)max=hRawYieldSBfit->GetMaximum();
  }
  hRawYieldRot->SetMaximum(max*1.2);

  hRawYieldLS->SetMarkerStyle(22);
  hRawYieldLS->SetMarkerColor(kGreen+2);
  hRawYieldLS->SetLineColor(kGreen+2);
  hRawYieldLS->Draw("PZSAME");
  hRawYieldME->SetMarkerStyle(25);
  hRawYieldME->SetMarkerColor(4);
  hRawYieldME->SetLineColor(4);
  hRawYieldME->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSBfit->SetMarkerStyle(27);
    hRawYieldSBfit->SetMarkerColor(6);
    hRawYieldSBfit->SetLineColor(6);
    hRawYieldSBfit->Draw("PSAME");
  }
  TLegend* legry=new TLegend(0.7,0.7,0.89,0.89);
  legry->SetFillStyle(0);
  legry->SetBorderSize(0);
  legry->AddEntry(hRawYieldRot,"Rotational","PL")->SetTextColor(1);
  legry->AddEntry(hRawYieldLS,"Like Sign","PL")->SetTextColor(kGreen+2);
  legry->AddEntry(hRawYieldME,"Ev Mix","PL")->SetTextColor(4);
  if(tryDirectFit){
    legry->AddEntry(hRawYieldSBfit,"Direct Fit","PL")->SetTextColor(6);
  }
  legry->Draw();
  cry->SaveAs(Form("figures/RawYield_%s.eps",suffix.Data()));
  cry->SaveAs(Form("figures/RawYield_%s.root",suffix.Data()));

  cCompareResidualTrends->cd(1);
  TLegend *legRT=new TLegend(*legry);
  legRT->SetX1NDC(0.12);
  legRT->SetY1NDC(0.7);
  legRT->SetX2NDC(0.3);
  legRT->SetY2NDC(0.9);
  legRT->Draw();

  // NOW COMPARE YIELDS OBTAINED WITH BC FROM DIFFERENT APPROACHES
  TCanvas* cryBC=new TCanvas("cryBC","RawYield with BC",800,700);
  cryBC->SetLeftMargin(0.15);
  hRawYieldRot->Draw("P");
  hRawYieldME->Draw("PSAME");
  hRawYieldLS->Draw("PZSAME");
  if(tryDirectFit){
    hRawYieldSBfit->Draw("PSAME");
  }

  if(hRawYieldRotBC->GetMaximum()>max)max=hRawYieldRotBC->GetMaximum();
  if(hRawYieldLSBC->GetMaximum()>max)max=hRawYieldLSBC->GetMaximum();
  if(hRawYieldMEBC->GetMaximum()>max)max=hRawYieldMEBC->GetMaximum();
  if(tryDirectFit){
    if(hRawYieldSBfitBC->GetMaximum()>max)max=hRawYieldSBfitBC->GetMaximum();
  }
  hRawYieldRot->SetMaximum(max*1.2);
  
  hRawYieldRotBC->SetMarkerStyle(21);
  hRawYieldRotBC->SetLineStyle(2);
  hRawYieldRotBC->Draw("PSAME");


  hRawYieldLSBC->SetMarkerStyle(22);
  hRawYieldLSBC->SetMarkerColor(kGreen+2);
  hRawYieldLSBC->SetLineColor(kGreen+2);
  hRawYieldLSBC->SetLineStyle(2);
  hRawYieldLSBC->Draw("PZSAME");

  hRawYieldMEBC->SetMarkerStyle(25);
  hRawYieldMEBC->SetMarkerColor(4);
  hRawYieldMEBC->SetLineColor(4);
  hRawYieldMEBC->SetLineStyle(2);
  hRawYieldMEBC->Draw("PSAME");
  if(tryDirectFit){
    hRawYieldSBfitBC->SetMarkerStyle(27);
    hRawYieldSBfitBC->SetMarkerColor(6);
    hRawYieldSBfitBC->SetLineColor(6);
    hRawYieldSBfitBC->SetLineStyle(2);
    hRawYieldSBfitBC->Draw("PSAME");
  }

  TLegend* legryBC=new TLegend(0.7,0.7,0.89,0.89);
  legryBC->SetFillStyle(0);
  legryBC->SetBorderSize(0);
  legryBC->AddEntry(hRawYieldRot,"Rotational","PL")->SetTextColor(1);
  legryBC->AddEntry(hRawYieldRotBC,"Rotational BC","PL")->SetTextColor(1);
  legryBC->AddEntry(hRawYieldLS,"Like Sign ","PL")->SetTextColor(kGreen+2);
  legryBC->AddEntry(hRawYieldLSBC,"Like Sign BC","PL")->SetTextColor(kGreen+2);
  legryBC->AddEntry(hRawYieldME,"Ev Mix","PL")->SetTextColor(4);
  legryBC->AddEntry(hRawYieldMEBC,"Ev Mix BC","PL")->SetTextColor(4);
  if(tryDirectFit){
    legryBC->AddEntry(hRawYieldSBfit,"Direct Fit","PL")->SetTextColor(6);
    legryBC->AddEntry(hRawYieldSBfitBC,"Direct Fit BC","PL")->SetTextColor(6);
  }
  legryBC->Draw();
  cryBC->SaveAs(Form("figures/RawYieldBC_%s.eps",suffix.Data()));
  cryBC->SaveAs(Form("figures/RawYieldBC_%s.root",suffix.Data()));




  fout->cd();
  hRawYieldRot->Write();
  hRawYieldLS->Write();
  hRawYieldME->Write();
  hRawYieldSBfit->Write();
  hnEv->Write();
  fout->Close();

  return;
}


void WriteFitInfo(AliHFMassFitter *fitter, TH1D* histo){
  Double_t sig,esig;
  fitter->Signal(3,sig,esig);
  Double_t mean=fitter->GetMean();
  Double_t emean=fitter->GetMeanUncertainty();
  Double_t sigma=fitter->GetSigma();
  Double_t esigma=fitter->GetSigmaUncertainty();
  Double_t minBin=histo->FindBin(mean-3.*sigma);
  Double_t maxBin=histo->FindBin(mean+3.*sigma);
  Double_t back=histo->Integral(minBin,maxBin);
  TPaveText* tpar=new TPaveText(0.5,0.7,0.89,.87,"NDC");
  tpar->SetBorderSize(0);
  tpar->SetFillStyle(0);
  tpar->AddText(Form("Mean = %.3f #pm %.3f",mean,emean));
  tpar->AddText(Form("Sigma = %.3f #pm %.3f",sigma,esigma));
  tpar->SetTextColor(4);
  tpar->Draw();

  TPaveText* tss=new TPaveText(0.15,0.15,0.5,0.4,"NDC");
  tss->SetBorderSize(0);
  tss->SetFillStyle(0);
  tss->AddText(Form("S(3#sigma) = %.0f #pm %.0f",sig,esig));
  tss->AddText(Form("B(3#sigma) = %.3g",back));
  tss->AddText(Form("S/B (3#sigma) = %.4f",sig/back));
  tss->AddText(Form("Significance(3#sigma) = %.2f",sig/TMath::Sqrt(back+sig)));
  tss->SetTextColor(1);
  tss->Draw();

}
