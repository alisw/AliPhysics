// MACRO TO COMPUTE THE <pt> from data for D mesons 
// works for D0, Ds and D+; D* case should be tested since the side bands definition is less clear
//
// Requires 
// a 2D (pt,mass) histogram, with fine pt bins in pt (finer than those used for the spectrum, e.g. 100 or 200 MeV/c)
// a histogram of background yield vs pt (with the pt bins in which the <pt> has to be computed)
// a histogram of the invariant mass sigma (with the pt bins in which the <pt> has to be computed)
// optionally a histogram of signal yield (not fully needed) (with the pt bins in which the <pt> has to be computed)
// optionally a histogram of the mean of the signal peak obtained from the fit (with the pt bins in which the <pt> has to be computed)
// The last 4 histograms are among the output of the standard macro used for fitting the mass distributions
// optionally (but strongly recommended) an histograms with prompt D reco efficiencies (in finer bins than those in which the <pt> is wanted, otherwise they will have no effect). 
//            The binning of the efficiency histo should match that of the pt axis of the 2D (pt,mass) plots x rebin (where rebin is used to rebin the data histogram with the method Rebin)
//            e.g. if you have data in 200 MeV/c and efficiency in 500 MeV/c bins than the minimum binning for which data and efficiencis match is 1 GeV/c
//                    -> you have to use rebin=5 (or multiples of 5) to have the macro working (it stops otherwise)
// Options: 
// rebin (see above)
// usefitforsubtraction: with this switched on the pt shape of the backgrond is fitted with an exponential and the fit function is used for the background subtraction in the signal region
//            (useful if there are stat fluctuations in the background pt shape distribution, but take care and check the fit quality in the control plots)
// correct for efficiency: if switched on, the pt distribution of the signal (i.e. that in the signal peak after back subtraction) is weighted by the efficiency
//                         it accounts for the pt dependence of the efficiency -> compulsory if a realistic (corrected) <pt> has to be computed 
//
// useParGenAccOverLimAcc: uses a parametrization for correcting the efficiencies by GenAcc/LimAcc. It is useful since the GenAcc/LimAcc factor might 
//                         not be available in very fine bins, especially at high pt  
//
//
// There are standard methods hard-coded for D0 (DoStandardForD0), D+ and Ds as an example (of course you have to change the paths of the files and the name of the histos!!).
//     To run them (e.g. for D0): 
//    -1 (my advice): copy the macro in your working directory and change the name of the paths and histo there (it helps to retrieve the results in the future). 
//                    Or copy the ingredient files in the working directory,
//     0- set the path of the files 
//     1- .L CalculateAveragePt.C++  (compile the macro, not really needed)
//     2- DoStandardForD0(5,kFALSE,kTRUE,kTRUE,0,3)
//    ... it should be easy :)

//
//
//  For problems write to: andrea.rossi@cern.ch
//
//////////////////////////////////////////////////////////////////////////////////// 

#include <TH2F.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TLegend.h>

TH2F *hPtInvMass=0x0;
Int_t nptbins=0;
Double_t *ptbinlimits=0x0;
Double_t *rawsignal=0x0,*rawback=0x0,*sigma=0x0,*meansignal=0x0;
Double_t nsigmaSignal=3.,nsigmaSBstart=5.;
Double_t mesonMass=1.8645;
Int_t nbinsx;//=hPtInvMass->GetNbinsX();
Int_t nbinsy;//=hPtInvMass->GetNbinsY();
Double_t binwidthpt;//=hPtInvMass->GetYaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
Double_t binwidthInvMss;//=hPtInvMass->GetXaxis()->GetBinWidth(1);// ASSUME ALL BINS ARE WITH THE SAME WIDTH
Double_t ptmin;//=hPtInvMass->GetYaxis()->GetBinLowEdge(1);
Double_t ptmax;//=hPtInvMass->GetYaxis()->GetBinUpEdge(nbinsy);
Bool_t useFitForSubtraction=kFALSE;
Bool_t useParGenAccOverLimAcc=kFALSE;
Bool_t subtractSB=kTRUE;
Bool_t corrForEff=kFALSE;
void SetSubtractSB(Bool_t subtract){subtractSB=subtract;}
void SetCorrForEff(Bool_t correff){corrForEff=correff;}
void SetUseFitForSubtraction(Bool_t useFit){useFitForSubtraction=useFit;}
void SetNsigmaForSignal(Double_t nsigm){nsigmaSignal=nsigm;}
void SetNsigmaStartSB(Double_t nsigm){nsigmaSBstart=nsigm;  
}
void CalculateAveragePt(Int_t rebin=1,Int_t firstbin=0,Int_t lastbin=0);
TF1* ParametricGenAccOverLimAccCorr(){// parametric form of GenAcc/LimAcc correction for D0->Kpi
  TF1 *fGenAccOverLimAcc=new TF1("fat7","[2]+([0]*TMath::ATan([1]*x+[3])+[4]*TMath::Log(1.+x/[5]))",0,24);
  fGenAccOverLimAcc->SetParameters(2.76153e-01,6.69658e-01,8.45865e-01,-1.87709e+00,2.89845e+03,1.39651e+05);  
  return fGenAccOverLimAcc;
}

Int_t standrebin[8]={2,2,2,2,2,2,4,5};
TCanvas **cPtDistrNoSubtr;
TF1 **fitfunc;
TH1D *hAvRawYieldSpectrum,*hSignal,*hMeanSignal,*hEfficNum,*hEfficDenum,*hBack,*hSigma;
void SetHistosEfficiency(TH1D *hNum,TH1D *hDenum){
  hEfficNum=(TH1D*)hNum;
  hEfficDenum=(TH1D*)hDenum;
  
}
TGraphAsymmErrors *grAvRawYieldSpectrum;
TGraphAsymmErrors *grAvPtVSPtmean;

TGraphAsymmErrors *grBackAvPtVSPtmean;

void SetPtBinLimits(const Int_t npt,Double_t *ptbinlim){
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

void SetHistMean(TH1D *hM){hMeanSignal=hM;
 meansignal=&(hMeanSignal->GetArray())[1]; 
 return;
}

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

TH1D* CheckBinningAndMerge(TH1D *hA,TH1D *hB,Double_t precision=0.001,Double_t minX=-9999.,Double_t maxX=-9999.){// Checks that histo hB binning is finer and fits hA binning
  // and return an Histo with hB bins integrated to match than hA bins  
  // TO BE ADDED A CHECK ON THE MIN AND MAX VALUES OF THE AXES IN THE 2 PLOTS
  // NO PARTICULAR TREATMENT OF THE ERROR IN THE INTEGRATION (SHOULD BE OK BUT FOR EFFIENCIES)

  TH1D *hWork=new TH1D(*hA);
  hWork->SetName(Form("%sNewBins",hB->GetName()));
  hWork->SetTitle(hB->GetTitle());
  
  Double_t lebinA,uebinA;
  Int_t binBleA=0,binBueA=0;
  Int_t binBleALast=0,binBueALast=0;
  Int_t maxbinA=hA->GetNbinsX();
  Int_t minbinA=1;
  if(maxX!=-9999.){
    maxbinA=hA->FindBin(maxX);
  }
  if(minX!=-9999.){
    minbinA=hA->FindBin(minX);
  }
  
  for(Int_t binA=minbinA;binA<=maxbinA;binA++){
    uebinA=hA->GetXaxis()->GetBinUpEdge(binA);
    lebinA=hA->GetXaxis()->GetBinLowEdge(binA);
    // Check that Low Edge of binA is close to the low or upper edge of a bin in hB
    // within precition
    binBleA=hB->FindBin(lebinA);    
    
    if(TMath::Abs(hB->GetXaxis()->GetBinLowEdge(binBleA)-lebinA)>precision&&TMath::Abs(hB->GetXaxis()->GetBinUpEdge(binBleA)-lebinA)>precision){
      printf("The bins of the efficiency histo do not match those of the pt spectrum histo (edges 1):\n Cannot correct for efficiensies \n");
      delete hWork;
      return 0x0;
    }
    if(TMath::Abs(hB->GetXaxis()->GetBinLowEdge(binBleA)-lebinA)>TMath::Abs(hB->GetXaxis()->GetBinUpEdge(binBleA)-lebinA))binBleA++;
    if(binBleA<=binBleALast){
      printf("Problems with rebinning Hist B, non compatibility with different bin of Hist A (1)\n");
      delete hWork;
      return 0x0;
    }
    if(hB->GetBinWidth(binBleA)-hA->GetBinWidth(binA)>precision){
      printf("The bins of the efficiency histo do not match those of the pt spectrum histo (bin width 1):\n Cannot correct for efficiensies \n");
      delete hWork;
      return 0x0;
    }
    // Check that Upper Edge of binA is close to the low or upper edge of a bin in hB
    // within precition
    binBueA=hB->FindBin(uebinA);
    if(TMath::Abs(hB->GetXaxis()->GetBinLowEdge(binBueA)-uebinA)>precision&&TMath::Abs(hB->GetXaxis()->GetBinUpEdge(binBueA)-uebinA)>precision){
      printf("The bins of the efficiency histo do not match those of the pt spectrum histo (edges 2) (bin %d Vs %d):\n Cannot correct for efficiensies \n",binA,binBueA);
      printf("%f Vs %f or %f \n",uebinA,hB->GetXaxis()->GetBinLowEdge(binBueA),hB->GetXaxis()->GetBinUpEdge(binBueA));
      delete hWork;
      return 0x0;
    }
    if(TMath::Abs(hB->GetXaxis()->GetBinLowEdge(binBueA)-uebinA)<TMath::Abs(hB->GetXaxis()->GetBinUpEdge(binBueA)-uebinA))binBueA--;
    if(binBueA<=binBueALast){
      printf("Problems with rebinning Hist B, non compatibility with different bin of Hist A (2) \n");
      delete hWork;
      return 0x0;
    }
    if(hB->GetBinWidth(binBueA)-hA->GetBinWidth(binA)>precision){
      printf("The bins of the efficiency histo do not match those of the pt spectrum histo (bin width 2):\n Cannot correct for efficiensies \n");
      delete hWork;
      return 0x0;
    }
    if(binBleA>binBueA){
      printf("Mathcing failure (probably a bug in the code \n");
      delete hWork;
      return 0x0;
    }
    
    hWork->SetBinContent(binA,hB->Integral(binBleA,binBueA));
  }
  
  return hWork;
 
}

Bool_t CorrectForEfficiency(TH1D *hPtHisto){
  TH1D *hWorkNum=0x0,*hWorkDenum=0x0;
  //  hWorkNum=CheckBinningAndMerge(hPtHisto,hEfficNum,0.01,1.,23.9);
  //  printf("ptmin and pt max: %f, %f \n",ptbinlimits[0],ptbinlimits[nptbins]);
  hWorkNum=CheckBinningAndMerge(hPtHisto,hEfficNum,0.01,ptbinlimits[0]*1.0001,ptbinlimits[nptbins]*0.9999);
  hWorkNum->Sumw2();
 
  if(!hWorkNum)return kFALSE;
  hWorkDenum=CheckBinningAndMerge(hPtHisto,hEfficDenum,0.01,ptbinlimits[0]*1.0001,ptbinlimits[nptbins]*0.9999);
  //  hWorkDenum=CheckBinningAndMerge(hPtHisto,hEfficDenum,0.01,1.,23.9);; 
  hWorkDenum->Sumw2();
  if(!hWorkDenum){
    delete hWorkNum;
    return kFALSE;
  }
  TH1D *hWorkNumCl=(TH1D*)hWorkDenum->Clone("hWorkNumCl");
  hWorkNumCl->Divide(hWorkNum,hWorkDenum,1,1,"B");
  if(useParGenAccOverLimAcc){
    TF1 *fGenAccoOverLimAcc=ParametricGenAccOverLimAccCorr();
    for(Int_t j=1;j<=hWorkNumCl->GetNbinsX();j++){
      printf(" EFFF for %f: %f * %f \n",hWorkNumCl->GetBinCenter(j), hWorkNumCl->GetBinContent(j),fGenAccoOverLimAcc->Eval(hWorkNumCl->GetBinCenter(j)));
      hWorkNumCl->SetBinContent(j,hWorkNumCl->GetBinContent(j)*fGenAccoOverLimAcc->Eval(hWorkNumCl->GetBinCenter(j)));
      hWorkNumCl->SetBinError(j,hWorkNumCl->GetBinError(j)*fGenAccoOverLimAcc->Eval(hWorkNumCl->GetBinCenter(j)));
   
    }
  }
  hPtHisto->Divide(hWorkNumCl);
  
  delete hWorkNum;
  delete hWorkNumCl;
  delete hWorkDenum;
  return kTRUE;
}

TH1D *HistoPtShapeFromData(Int_t ptbin,Int_t rebin=1.){
  if(!hPtInvMass){
    printf("No histo set \n");
    return 0x0;
  }
  Double_t mean=mesonMass;
  if(meansignal!=0x0){
    mean=meansignal[ptbin];    
    printf("USING MEAN FROM MASS FIT: %f\n",mean);
  }
  TString namehist="hPtDistrSB_Bin",funcname;
  namehist+=ptbin;

  TString cname="cPtDistrNoSubBin";
  cname+=ptbin;
  cPtDistrNoSubtr[ptbin]=new TCanvas(cname.Data(),"CPtDistrNoSub",700,700);
  
  printf("After Setting Canva \n");
  Int_t binleftpt=-1,binrightpt=-1;
  TH1D *hPtDistrSB=new TH1D(namehist.Data(),"Side band Pt distribution",nbinsy,ptmin,ptmax);
  namehist="hPtDistrTotPeak_Bin";
  namehist+=ptbin;
  printf("After Setting SB hist \n");
  TH1D *hPtDistrTotPeak=new TH1D(namehist.Data(),"Total Pt distribution in signal mass reagion",nbinsy,ptmin,ptmax);  
  for(Int_t jbi=1;jbi<hPtDistrTotPeak->GetNbinsX();jbi++){
    hPtDistrTotPeak->SetBinContent(jbi,0);
    hPtDistrTotPeak->SetBinError(jbi,0);
    hPtDistrSB->SetBinContent(jbi,0);
    hPtDistrSB->SetBinError(jbi,0);
  }
  


  printf("Before loop on xbins \n");
  for(Int_t j=1;j<nbinsx;j++){
    if((hPtInvMass->GetXaxis()->GetBinLowEdge(j)<=mean-nsigmaSBstart*sigma[ptbin])||(hPtInvMass->GetXaxis()->GetBinUpEdge(j)>=mean+nsigmaSBstart*sigma[ptbin])){
      for(Int_t k=1;k<nbinsy;k++){
	if(hPtInvMass->GetYaxis()->GetBinLowEdge(k)>=ptbinlimits[ptbin]&&hPtInvMass->GetYaxis()->GetBinUpEdge(k)<=ptbinlimits[ptbin+1]){
	  hPtDistrSB->SetBinContent(k,hPtDistrSB->GetBinContent(k)+hPtInvMass->GetBinContent(j,k));
	  hPtDistrSB->SetBinError(k,TMath::Sqrt(hPtDistrSB->GetBinError(k)*hPtDistrSB->GetBinError(k)+hPtInvMass->GetBinError(j,k)*hPtInvMass->GetBinError(j,k)));
	}
      }
    }
    else if((hPtInvMass->GetXaxis()->GetBinLowEdge(j)>=mean-nsigmaSignal*sigma[ptbin])&&(hPtInvMass->GetXaxis()->GetBinUpEdge(j)<=mean+nsigmaSignal*sigma[ptbin])){
      for(Int_t k=1;k<nbinsy;k++){
	if(hPtInvMass->GetYaxis()->GetBinLowEdge(k)>=ptbinlimits[ptbin]&&hPtInvMass->GetYaxis()->GetBinUpEdge(k)<=ptbinlimits[ptbin+1]){
	  if(binleftpt==-1)binleftpt=k;
	  hPtDistrTotPeak->SetBinContent(k,hPtDistrTotPeak->GetBinContent(k)+hPtInvMass->GetBinContent(j,k));
	  hPtDistrTotPeak->SetBinError(k,TMath::Sqrt(hPtDistrTotPeak->GetBinError(k)*hPtDistrTotPeak->GetBinError(k)+hPtInvMass->GetBinError(j,k)*hPtInvMass->GetBinError(j,k)));
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
  
  Double_t avptbin=0.,errA=0.,errB=0.,errC=0.,errisq=0.;
  // CALCULATE BACKGROUND AVERAGE PT
  for(Int_t bbin=1;bbin<=hPtDistrSB->GetNbinsX();bbin++){
    avptbin+=hPtDistrSB->GetBinCenter(bbin)*hPtDistrSB->GetBinContent(bbin)*hPtDistrSB->GetBinWidth(bbin);
    errisq=hPtDistrSB->GetBinError(bbin)*hPtDistrSB->GetBinError(bbin);
    errA+=hPtDistrSB->GetBinCenter(bbin)*hPtDistrSB->GetBinCenter(bbin)*errisq;
    errB+=hPtDistrSB->GetBinCenter(bbin)*errisq;
    errC+=errisq;
    // printf("ERROR VALUES: %f, %f \n",errA,errB);
  }
  avptbin/=hPtDistrSB->Integral("width");
  errisq=TMath::Sqrt(errA-2.*avptbin*errB+avptbin*avptbin*errC)/hPtDistrSB->Integral("width");
  if(!grBackAvPtVSPtmean){
    grBackAvPtVSPtmean=new TGraphAsymmErrors();
    grBackAvPtVSPtmean->SetName("grBackAvPtVSPtmean");
  }
  Int_t grbin=grBackAvPtVSPtmean->GetN();
  grBackAvPtVSPtmean->SetPoint(grbin,0.5*(ptbinlimits[ptbin+1]+ptbinlimits[ptbin]),avptbin);
  grBackAvPtVSPtmean->SetPointError(grbin,0.5*(ptbinlimits[ptbin+1]-ptbinlimits[ptbin]),0.5*(ptbinlimits[ptbin+1]-ptbinlimits[ptbin]),errisq,errisq);

  TH1D *hPtSignalFromSubtraction=new TH1D(*hPtDistrTotPeak);
  namehist="hPtSignalFromSubtraction_Bin";
  namehist+=ptbin;
  hPtSignalFromSubtraction->SetName(namehist.Data());
  hPtSignalFromSubtraction->SetTitle("Signal Pt distr from subtraction");


  if(subtractSB){
    if(!useFitForSubtraction)hPtSignalFromSubtraction->Add(hPtDistrSB,-rawback[ptbin]*(nsigmaSignal/3.)/hPtDistrSB->Integral());// ROUGH LINEAR SCALING OF BACKGOUND!!
    else{
      funcname="expofit";
      funcname+=ptbin;
      fitfunc[ptbin]=new TF1(funcname.Data(),"expo",ptbinlimits[ptbin],ptbinlimits[ptbin+1]);        
      TCanvas *cTempt=new TCanvas();
      cTempt->cd();
      hPtDistrSB->Fit(fitfunc[ptbin],"RLEV","",ptbinlimits[ptbin],ptbinlimits[ptbin+1]);
      delete cTempt;
      for(Int_t jb=binleftpt;jb<=binrightpt;jb++){
	printf("Err before sub: %f \n",hPtSignalFromSubtraction->GetBinError(jb));
	hPtSignalFromSubtraction->SetBinContent(jb,hPtSignalFromSubtraction->GetBinContent(jb)-fitfunc[ptbin]->Eval(hPtSignalFromSubtraction->GetBinCenter(jb))/fitfunc[ptbin]->Integral(ptbinlimits[ptbin],ptbinlimits[ptbin+1])*binwidthpt*rawback[ptbin]*(nsigmaSignal/3.));// ROUGH LINEAR SCALING OF BACKGOUND!!
	printf("Err after sub: %f \n\n",hPtSignalFromSubtraction->GetBinError(jb));
      }    
    }
  }
  
  
  cPtDistrNoSubtr[ptbin]->Divide(1,2);
  cPtDistrNoSubtr[ptbin]->cd(1);
  hPtDistrSB->Draw();
  cPtDistrNoSubtr[ptbin]->cd(2);
  if(corrForEff){
    if(!CorrectForEfficiency(hPtSignalFromSubtraction)){
      printf("Correction for efficieny failed \n");
    }    
  }
  hPtDistrTotPeak->Draw();  
  hPtSignalFromSubtraction->Draw("sames");

  
  return hPtSignalFromSubtraction;

}

void CalculateAveragePt(TH2* hMassPt,TH1D *hB,TH1D *hSigm,TH1D *hEffNum=0x0,TH1D *hEffDenum=0x0,TH1D *hS=0x0,TH1D *hMean=0x0,Int_t rebin=1,Int_t firstbin=3,Int_t lastbin=8){
  
  SetPtBinLimits(hB);
  SetHistoMassPt((TH2F*)hMassPt);
  SetHistRawBack(hB);
  SetHistSigma(hSigm);
  if(((!hEffNum)&&(hEffDenum))||((hEffNum)&&(!hEffDenum))){
    printf("Two histos are needed for the efficiency: numerator and denumerator (for rebinning) \n");
  }
  else if(hEffNum){
    SetHistosEfficiency(hEffNum,hEffDenum);
    //    SetHistEffNum(hEffNum);
    // SetHistEffDenum(hEffDenum);
  }
  if(hS)SetHistRawSignal(hS);  
  if(hMean)SetHistMean(hMean);  
  printf("Histos set \n");
  CalculateAveragePt(rebin,firstbin,lastbin);
}

void CalculateAveragePt(Int_t rebin,Int_t firstbin,Int_t lastbin){
  
  TH1D *hPtSpectra;
  if(rebin>0){
    hPtSpectra=new TH1D("hPtDistrTotPeak","Total Pt distribution in signal mass reagion",nbinsy,ptmin,ptmax); 
    if(rebin!=1)hPtSpectra->Rebin(rebin);
  }
  else if(rebin==-1){
    printf("Standard rebinning not implemented yet  \n");return;
  }
  else {
    printf("Wrong rebin numnber \n");return;
  }

  grAvRawYieldSpectrum=new TGraphAsymmErrors();
  grAvPtVSPtmean=new TGraphAsymmErrors();
  grBackAvPtVSPtmean=new TGraphAsymmErrors();
  grBackAvPtVSPtmean->SetName("grBackAvPtVSPtmean");
  

  Double_t avptbin=0.,errA,errB,errC,errisq;
  for(Int_t bin=firstbin;bin<=lastbin;bin++){    
    //    printf("Before Getting Subtracted Histo \n");
    avptbin=0.;
    errA=0.;
    errB=0.;
    errC=0.;
    errisq=0.;
    
    TH1D *h=HistoPtShapeFromData(bin,rebin);
    hPtSpectra->Add(h);
    hAvRawYieldSpectrum->SetBinContent(bin+1,h->Integral("width")/hAvRawYieldSpectrum->GetBinWidth(bin+1));
    for(Int_t bbin=1;bbin<=h->GetNbinsX();bbin++){
      avptbin+=h->GetBinCenter(bbin)*h->GetBinContent(bbin)*h->GetBinWidth(bbin);
      errisq=h->GetBinError(bbin)*h->GetBinError(bbin);
      errA+=h->GetBinCenter(bbin)*h->GetBinCenter(bbin)*errisq;
      errB+=h->GetBinCenter(bbin)*errisq;
      errC+=errisq;
      // printf("ERROR VALUES: %f, %f \n",errA,errB);
    }
    avptbin/=h->Integral("width");
    errisq=TMath::Sqrt(errA-2.*avptbin*errB+avptbin*avptbin*errC)/h->Integral("width");

    grAvRawYieldSpectrum->SetPoint(bin-firstbin,avptbin,h->Integral("width")/hAvRawYieldSpectrum->GetBinWidth(bin+1));
    grAvRawYieldSpectrum->SetPointError(bin-firstbin,avptbin-ptbinlimits[bin],ptbinlimits[bin+1]-avptbin,errisq/2.,errisq/2.);   

    grAvPtVSPtmean->SetPoint(bin-firstbin,0.5*(ptbinlimits[bin+1]+ptbinlimits[bin]),avptbin);
    grAvPtVSPtmean->SetPointError(bin-firstbin,0.5*(ptbinlimits[bin+1]-ptbinlimits[bin]),0.5*(ptbinlimits[bin+1]-ptbinlimits[bin]),errisq,errisq);
    

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



   //  TFile *ftheory=TFile::Open("/Users/administrator/ALICE/CHARM/ppData_2010/PtSelInBin/predictions_D0.root");
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
  grAvPtVSPtmean->GetYaxis()->SetTitle("<p_{t}> (GeV/c)");
  grAvPtVSPtmean->GetXaxis()->SetTitle("bin centre p_{t} (GeV/c)");

  grBackAvPtVSPtmean->SetMarkerColor(kRed);
  grBackAvPtVSPtmean->SetMarkerStyle(24);
  grBackAvPtVSPtmean->SetLineColor(kRed);
  grBackAvPtVSPtmean->GetYaxis()->SetTitle("<p_{t}> (GeV/c)");
  grBackAvPtVSPtmean->GetXaxis()->SetTitle("bin centre p_{t} (GeV/c)");
  grBackAvPtVSPtmean->Draw("p");
  

  TLegend *leg=new TLegend(0.7,0.5,0.95,0.2,"","NDC");
  leg->AddEntry(grAvPtVSPtmean,"signal","lp");
  leg->AddEntry(grBackAvPtVSPtmean,"background","lp");
  leg->Draw();

  // PRINTING RESULTS
  Double_t *ygrSignal,*ygrBack,*eygrSignal,*eygrBack;
  ygrSignal=grAvPtVSPtmean->GetY();
  eygrSignal=grAvPtVSPtmean->GetEYhigh();
  ygrBack=grBackAvPtVSPtmean->GetY();
  eygrBack=grBackAvPtVSPtmean->GetEYhigh();
  printf("Av Pt for signal: \n");
  for(Int_t in=0;in<grAvPtVSPtmean->GetN();in++){
    cout<<ygrSignal[in]<<endl;
  }  
  printf("Error on Av Pt for signal: \n");
  for(Int_t in=0;in<grAvPtVSPtmean->GetN();in++){
    cout<<eygrSignal[in]<<endl;
  }

 printf("Av Pt for back: \n");
  for(Int_t in=0;in<grAvPtVSPtmean->GetN();in++){
    cout<<ygrBack[in]<<endl;
  }  
  printf("Error on Av Pt for back: \n");
  for(Int_t in=0;in<grAvPtVSPtmean->GetN();in++){
    cout<<eygrBack[in]<<endl;
  }

  
  

  nameout="ptCorrection";
  if(useFitForSubtraction)nameout.Append("FitSB");  
  if(!subtractSB)nameout.Append("NoSBsubtract.root");
  else nameout.Append(".root");
  TFile *fout=new TFile(nameout.Data(),"RECREATE");
  fout->cd();
  grAvPtVSPtmean->Write();
  grAvRawYieldSpectrum->Write();
  grBackAvPtVSPtmean->Write();
  cPtSpectra->Write();
  for(Int_t bin=firstbin;bin<=lastbin;bin++){  
    cPtDistrNoSubtr[bin]->Write();
  }



  fout->Close();

  //  histoThCompareCentral->Draw("sames");
}


void DoStandardForD0(Int_t rebin,Bool_t usefit,Bool_t corrforeff=kTRUE,Bool_t useParGenAccLimacc=kTRUE,Int_t firstbin=0,Int_t lastbin=3){// firstbin and lastbin are w.r.t. the signal histo plot, starting from 0
  
  //N.B.  ADAPT THE PART BELOW WITH YOUR PATHS AND HISTO NAMES
  TFile *fData=TFile::Open("/Users/administrator/ALICE/CHARM/PbPBdata_2011/TestTrain/2013June4TrainData92_MC61/Data/Standard/RAAvsNPart/FinalMassPlots_v2/RawYieldBoth_tight.root");
  TFile *fCF=TFile::Open("/Users/administrator/ALICE/CHARM/PbPBdata_2011/2013Jun8MCptWeightDataMoreSplit/MergeWithMyCode/Standard/EffPrompt/fileEff_D0_from_c.root");//EffForHFspectrumLHC10f6NOTd4Corr.root");
  TFile *fDataList=TFile::Open("/Users/administrator/ALICE/CHARM/PbPBdata_2011/TestTrain/2013June4TrainData92_MC61/Data/AnalysisResults.root");
  TDirectory *fdir=(TDirectory*)fDataList->Get("PWG3_D2H_d0D0");
  TList *lslist=(TList*)fdir->Get("clistTGHCsign_d0D0");
  TH2F *hMassPt=(TH2F*)lslist->FindObject("hInvMassPtTGHCsign");
  TH1D *hSig=(TH1D*)fData->Get("hSigma");
  TH1D *hBackground=(TH1D*)fData->Get("hBackground");
  TH1D *hSign=(TH1D*)fData->Get("hSignal");
  TH1D *hMean=(TH1D*)fData->Get("hMass");

  mesonMass=1.86484;
  
  TH1D *hEffNum=(TH1D*)fCF->Get("hRecoPIDpt");
  TH1D *hEffDeNum=(TH1D*)fCF->Get("hMCAccpt");

  /////////////////////////////////////////////////////////////////  
  //   FOR DEBUGGING: CHECK EFFICIENCY PLOTS REBINNING
  //   TCanvas *cEff=new TCanvas();
  //   cEff->cd();
  //   TH1D *hRatio=new TH1D(*(TH1D*)hEffNum);
  //   hRatio->Divide(hEffDeNum);

  //   TH1D *hPt=(TH1D*)(hMassPt->ProjectionY());  
  
  //   hPt->Rebin(rebin);cout<< hPt->GetBinWidth(1)<<endl;
  //   cout<< hEffNum->GetBinWidth(1)<<endl;
  //   TH1D *hWorkNum=0x0;
  //   hWorkNum=CheckBinningAndMerge(hPt,(TH1D*)hEffNum,0.01,23.9);
  //   TH1D *hWorkDenum=0x0;
  //   hWorkDenum=CheckBinningAndMerge(hPt,(TH1D*)hEffDeNum,0.01.,23.9);
  //   TH1D *hWorkRatio=new TH1D(*hWorkNum);
  //   hWorkRatio->Divide(hWorkDenum);
  
  
  //   hWorkRatio->Draw();
  //   hPt->Draw("sames");
  
  //   hRatio->Draw("sames");
  /////////////////////////////////////////////////////////////////

  SetNsigmaForSignal(3.);
  SetNsigmaStartSB(5.);
  SetUseFitForSubtraction(usefit);
  corrForEff=corrforeff;
  useParGenAccOverLimAcc=useParGenAccLimacc;
  CalculateAveragePt(hMassPt,hBackground,hSig,hEffNum,hEffDeNum,hSign,hMean,rebin,firstbin,lastbin);
  



}


void DoStandardForDs(Int_t rebin,Bool_t usefit,Bool_t corrforeff=kTRUE,Int_t firstbin=0,Int_t lastbin=3){// firstbin and lastbin are w.r.t. the signal histo plot, starting from 0

  //N.B.  ADAPT THE PART BELOW WITH YOUR PATHS AND HISTO NAMES
TFile *fData=TFile::Open("RawYieldBoth.root"); 
TFile *fCF=TFile::Open("fileEff_Ds_CommonFramework_from_c_Enriched.root"); 

TH1D *hSig=(TH1D*)fData->Get("hSigma"); 
TH1D *hMean=(TH1D*)fData->Get("hMass"); 
TH1D *hBackground=(TH1D*)fData->Get("hBackground"); 
TH1D *hSign=(TH1D*)fData->Get("hSignal"); 

TH1D *hEffNum=(TH1D*)fCF->Get("hRecoPIDpt"); 
TH1D *hEffDeNum=(TH1D*)fCF->Get("hMCAccpt"); 
 if(!hEffNum)return;
 if(!hEffDeNum)return;
 TFile *file=new TFile("AnalysisResults.root"); 
 TDirectory *dirFile=(TDirectory*)file->Get("PWG3_D2H_InvMassDs"); 
 TList *cOutput = (TList*)dirFile->Get("coutputDs"); 
 TH2F *hMassPt=(TH2F*)cOutput->FindObject("hPtVsMass");
 //TH2F *hMassPt=(TH2F*)fData->Get("hPtVsMass"); 
 
 mesonMass=1.96847;

 SetUseFitForSubtraction(usefit);
 corrForEff=corrforeff;
 CalculateAveragePt(hMassPt,hBackground,hSig,hEffNum,hEffDeNum,hSign,hMean,rebin,firstbin,lastbin);
 

}
void DoStandardForDplus(Int_t rebin,Bool_t usefit,Bool_t corrforeff=kTRUE,Int_t firstbin=0,Int_t lastbin=2){
// firstbin and lastbin are w.r.t. the signal histo plot, starting from 0

  //N.B.  ADAPT THE PART BELOW WITH YOUR PATHS AND HISTO NAMES
  
  TFile *fData=TFile::Open("/Users/administrator/ALICE/CHARM/ppData_2010/PtSelInBin/testWithEff/Dplus/RawYieldBoth.root");
  TFile *fCF=TFile::Open("/Users/administrator/ALICE/CHARM/ppData_2010/PtSelInBin/testWithEff/Dplus/CFEfficiency.root");
  TH2F *hMassPt=(TH2F*)fData->Get("hPtVsMassTC");
  TH1D *hSig=(TH1D*)fData->Get("hSigma");
  TH1D *hMean=(TH1D*)fData->Get("hMean");
  TH1D *hBackground=(TH1D*)fData->Get("hBackground");
  TH1D *hSign=(TH1D*)fData->Get("hSignal");
  
  TH1D *hEffNum=(TH1D*)fCF->Get("RecoPID");
  TH1D *hEffDeNum=(TH1D*)fCF->Get("GenAcc");
  
  mesonMass=1.86962; 

  /////////////////////////////////////////////////////////////////
  //   FOR DEBUGGING: CHECK EFFICIENCY PLOTS REBINNING
  //   TCanvas *cEff=new TCanvas();
  //   cEff->cd();
  //   TH1D *hRatio=new TH1D(*(TH1D*)hEffNum);
  //   hRatio->Divide(hEffDeNum);

  //   TH1D *hPt=(TH1D*)(hMassPt->ProjectionY());  
  
  //   hPt->Rebin(rebin);cout<< hPt->GetBinWidth(1)<<endl;
  //   cout<< hEffNum->GetBinWidth(1)<<endl;
  //   TH1D *hWorkNum=0x0;
  //   hWorkNum=CheckBinningAndMerge(hPt,(TH1D*)hEffNum,0.01,23.9);
  //   TH1D *hWorkDenum=0x0;
  //   hWorkDenum=CheckBinningAndMerge(hPt,(TH1D*)hEffDeNum,0.01.,23.9);
  //   TH1D *hWorkRatio=new TH1D(*hWorkNum);
  //   hWorkRatio->Divide(hWorkDenum);
  //   hWorkRatio->Draw();
  //   hPt->Draw("sames");
  //   hRatio->Draw("sames");
  /////////////////////////////////////////////////////////////////  

  SetUseFitForSubtraction(usefit);
  corrForEff=corrforeff;
  CalculateAveragePt(hMassPt,hBackground,hSig,hEffNum,hEffDeNum,hSign,hMean,rebin,firstbin,lastbin);
  



}




TH1D *SmearEffHisto(TH1D *hInput,TString name="hEffNum",Double_t maxPt=40.,Double_t step=0.1){
  Int_t nbins=(Int_t)(maxPt/step);
  Int_t maxBinA=hInput->GetNbinsX();
  Int_t binBl,binBu;
  Double_t binlea,binuea,contInMin,contInMax;

  TH1D *h=new TH1D(name.Data(),name.Data(),nbins,0.,maxPt);
  for(Int_t bin=1;bin<=maxBinA;bin++){
    contInMin=hInput->GetBinContent(bin);
    if(bin==maxBinA) contInMax=contInMin*1.2;
    else contInMax=hInput->GetBinContent(bin+1);
    binlea=hInput->GetBinLowEdge(bin);
    binuea=hInput->GetXaxis()->GetBinUpEdge(bin);
    binBl=h->FindBin(binlea);
    binBu=h->FindBin(binuea);
    for(Int_t bb=binBl;bb<binBu;bb++){
      h->SetBinContent(bb,contInMin+(bb-binBl)*(contInMax-contInMin)/(binBu-binBl));
    } 
  }
  return h;
  
}


void AverageD0DplusResults(TString fileD0="/Users/administrator/ALICE/CHARM/ppData_2010/2011_Jul_05/data/LHC10bcdeAOD057/AvPt/MassRegSel3PkMore5SBEffCorrMeanFit/ptCorrectionFitSB.root",TString fileDplus="/Users/administrator/ALICE/CHARM/ppData_2010/2011_Jul_05/Dplus/AvPt/2011Jul26Renu/average_ptNew.root"){
  
  Double_t *yD0,*eyD0,*yDplus,*eyDplus,*x,*exLow,*exHigh;
  Double_t error=0.,average=0.;

  TFile *fD0=TFile::Open(fileD0.Data());
  // TCanvas *cD0=(TCanvas*)fD0->Get("cAvPtVSPtmean");
  //  TGraphAsymmErrors *grD0=(TGraphAsymmErrors*)cD0->FindObject("grAvPtMeanPtFitSB");
  TGraphAsymmErrors *grD0=(TGraphAsymmErrors*)fD0->Get("grAvPtMeanPtFitSB");

  grD0->SetName("grAvPtMeanPtFitSB_D0");
  yD0=grD0->GetY();
  eyD0=grD0->GetEYhigh();
  

  TFile *fDplus=TFile::Open(fileDplus.Data());
  TCanvas *cDplus=(TCanvas*)fDplus->Get("cAvPtVSPtmean");
  TGraphAsymmErrors *grDplus=(TGraphAsymmErrors*)cDplus->FindObject("grAvPtMeanPtFitSB");
  grDplus->SetName("grAvPtMeanPtFitSB_Dplus");
  yDplus=grDplus->GetY();
  eyDplus=grDplus->GetEYhigh();
  

  Int_t nlmax=grDplus->GetN();
  Int_t nlD0=grD0->GetN();
  Int_t nlDplus=grDplus->GetN();
  
  if(grD0->GetN()>nlmax){
    nlmax=grD0->GetN();
    x=grD0->GetX();
    exLow=grD0->GetEXlow();
    exHigh=grD0->GetEXhigh();
  }
  else {
    x=grDplus->GetX();
    exLow=grDplus->GetEXlow();
    exHigh=grDplus->GetEXhigh();
  }
  
  
  TGraphAsymmErrors *grAverage=new TGraphAsymmErrors();
  grAverage->SetName("grAverageD0Dplus");
  
  for(Int_t j=0;j<nlmax;j++){
    if(j<nlD0&&j<nlDplus){
      average=(yD0[j]/(eyD0[j]*eyD0[j])+yDplus[j]/(eyDplus[j]*eyDplus[j]))/(1./(eyD0[j]*eyD0[j])+1./(eyDplus[j]*eyDplus[j]));
      error=TMath::Sqrt(1./(1./(eyD0[j]*eyD0[j])+1./(eyDplus[j]*eyDplus[j])));
    }
    else if(j>=nlD0){
      average=yDplus[j];
      error=eyDplus[j];
    }
    else if(j>=nlDplus){
      average=yD0[j];
      error=eyD0[j];
    }
    grAverage->SetPoint(j,x[j],average);
    grAverage->SetPointError(j,exLow[j],exHigh[j],error,error);    
  }
  
  // Printing Info: D0
  printf("Average Value D0: \n");
  for(Int_t j=0;j<nlD0;j++){
    cout<<yD0[j]<<endl;
    
  }

 printf("Average Errors D0: \n");
  for(Int_t j=0;j<nlD0;j++){
    cout<<eyD0[j]<<endl;
    
  }
  
  // Printing Info: Dplus
  printf("Average Value Dplus: \n");
  for(Int_t j=0;j<nlDplus;j++){
    cout<<yDplus[j]<<endl;
    
  }

  printf("Average Errors Dplus: \n");
  for(Int_t j=0;j<nlDplus;j++){
    cout<<eyDplus[j]<<endl;
    
  }


  // Printing Info: Dplus D0 average
  yDplus=grAverage->GetY();
  eyDplus=grAverage->GetEYhigh();
  printf("Average Value D0 + Dplus: \n");
  for(Int_t j=0;j<nlDplus;j++){
    cout<<yDplus[j]<<endl;
    
  }

  printf("Average Errors D0+Dplus: \n");
  for(Int_t j=0;j<nlDplus;j++){
    cout<<eyDplus[j]<<endl;
    
  }
  
  


  TCanvas *cCompare=new TCanvas("cCompareD0Dplus","cCompareD0Dplus",700,700);
  cCompare->cd();
  
  grDplus->SetMarkerStyle(21);
  grDplus->SetMarkerSize(1.2);
  grDplus->SetMarkerColor(kBlue);
  grDplus->SetLineColor(kBlue);
  grDplus->Draw("ap");

  grD0->SetMarkerStyle(22);
  grD0->SetMarkerSize(1.2);
  grD0->SetMarkerColor(kRed);
  grD0->SetLineColor(kRed);
  grD0->Draw("p");

  grAverage->SetMarkerStyle(20);
  grAverage->SetMarkerSize(1.2);
  grAverage->SetMarkerColor(kBlack);
  grAverage->SetLineColor(kBlack);
  grAverage->Draw("p");
  
  

}
