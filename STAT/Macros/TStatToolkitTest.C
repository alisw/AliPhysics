/*
  .L $ALICE_ROOT/STAT/Macros/TStatToolkitTest.C+

*/

#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"
#include "TStatToolkit.h"
#include "TTreeStream.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TCanvas.h"

TObjArray arrayFit(3);
Int_t kMarkers[10]={25,24,20,21};
Int_t kColors[10]={1,2,4,3};
const char * names[10] = {"LTM","LThisto","LTMhisto1"};
const char * distNames[10] = {"Gauss","Gauss+flat","Gauss(m)+0.2*Gauss(m+5#sigma)","Gauss(m)+0.3*Gauss(m+3#sigma)"};


void TestLTM(Int_t nevents=10000){
  //
  // Goal test and benchamerk numerical stability and precission of the LTM method
  // Binned data and not binned data: 
  // Here we assume that binwidth<sigma
  // Distribution types examples used in test:
  //
  //   0 - gaussian
  //   1 - gauss+uniform background
  //   2 - gauss+second gaus 5 sigma away 
  //   
  Int_t npointsMax=1000;
  TTreeSRedirector * pcstream = new TTreeSRedirector("TStatToolkit_LTMtest.root","recreate");
  //
  TVectorD values(npointsMax);
  TVectorD vecLTM(6);
  TVectorD meanLTM(6);
  TVectorD sigmaLTM(6);
  TVectorD meanLTMHisto(6);
  TVectorD sigmaLTMHisto(6);
  TVectorD meanLTMHisto1(6);
  TVectorD sigmaLTMHisto1(6);
  TVectorD paramLTM(10);
  //
  for (Int_t iltm=0; iltm<6; iltm++) vecLTM[iltm]=0.99999*(0.5+iltm/10.);
  TF1 fg("fg","gaus");
  for (Int_t ievent = 0; ievent<nevents; ievent++){
    if (ievent%1000==0) printf("%d\n",ievent);
    Int_t distType=Int_t(gRandom->Rndm()*4);
    Int_t npoints= 50+(npointsMax-50)*gRandom->Rndm();
    Double_t mean  = 0.5+(gRandom->Rndm()-0.5)*0.2;
    Double_t sigma = mean*0.2*(1+2*gRandom->Rndm());
    TH1F histo("histo","histo",100,-0.5,1.5);
    //
    for (Int_t ipoint=0; ipoint<npoints; ipoint++){
      Double_t value=0;
      if (distType==0) value=gRandom->Gaus(mean,sigma); // gauss
      if (distType==1) {
	if (gRandom->Rndm()>0.2) {                     //  gauss + background
	  value=gRandom->Gaus(mean,sigma);
	}else{
	  value=mean+(gRandom->Rndm()-0.5)*2;
	}
      }
      if (distType==2) {
	if (gRandom->Rndm()>0.2) {                     //  gauss + second gaus 5 sigma away
	  value=gRandom->Gaus(mean,sigma);
	}else{
	  value=gRandom->Gaus(mean+5*sigma,sigma);
	}
      }
      if (distType==3) {
	if (gRandom->Rndm()>0.3) {                     //  gauss + second gaus 4 sigma away
	  value=gRandom->Gaus(mean,sigma);
	}else{
	  value=gRandom->Gaus(mean+5*sigma,sigma);
	}
      }
      values[ipoint]=value;
      histo.Fill(value);
    }
    //
    histo.Fit(&fg,"QN","QN");
    Double_t meanG = fg.GetParameter(1);
    Double_t rmsG = fg.GetParameter(2);
    Double_t meanA  = TMath::Mean(npoints,values.GetMatrixArray());
    Double_t rmsA   = TMath::Mean(npoints,values.GetMatrixArray());
    Double_t meanH  = histo.GetMean();
    Double_t rmsH  = histo.GetRMS();
    //
    for (Int_t iltm=0; iltm<6; iltm++){
      //    
      Double_t meanV,sigmaV=0;
      TStatToolkit::EvaluateUni(npoints,values.GetMatrixArray(), meanV,sigmaV, vecLTM[iltm]*npoints);
      meanLTM[iltm]=meanV;
      sigmaLTM[iltm]=sigmaV;
      //
      TStatToolkit::LTM(&histo, &paramLTM,  vecLTM[iltm]);
      meanLTMHisto1[iltm]=paramLTM[1];
      sigmaLTMHisto1[iltm]=paramLTM[2];
      //
      TStatToolkit::LTMHisto(&histo, paramLTM,  vecLTM[iltm]);
      meanLTMHisto[iltm]=paramLTM[1];
      sigmaLTMHisto[iltm]=paramLTM[2];
    }
    (*pcstream)<<"ltm"<<
      //Input
      "npoints="<<npoints<<
      "distType="<<distType<<
      "mean="<<mean<<
      "sigma="<<sigma<<
      // "Standart" statistic output
      "meanA="<<meanA<<
      "rmsA="<<rmsA<<
      "meanH="<<meanH<<
      "rmsH="<<rmsH<<
      "meanG="<<meanG<<
      "rmsG="<<rmsG<<
      // "LTM" output
      "vecLTM.="<<&vecLTM<<
      "meanLTM.="<<&meanLTM<<
      "meanLTMHisto.="<<&meanLTMHisto<<
      "meanLTMHisto1.="<<&meanLTMHisto1<<
      //
      "sigmaLTM.="<<&sigmaLTM<<
      "sigmaLTMHisto.="<<&sigmaLTMHisto<<
      "sigmaLTMHisto1.="<<&sigmaLTMHisto1<<      
      "\n";
  }
  delete pcstream;
  //
  //
  //
  TFile *fltm = TFile::Open("TStatToolkit_LTMtest.root","update");
  TTree * tree = (TTree*)fltm->Get("ltm");
  tree->SetMarkerSize(0.5);
  tree->SetMarkerStyle(25);
  //
  // 1. Get numerical error estimate of the LTM method for gaussian distribution  
  //
  TH2 * hisLTMTrunc[10]={0};  
  TH1 * hisResol[10]={0};
  TGraphErrors * grSigma[10]={0};
  TCanvas *canvasLTM= new TCanvas("canvasLTM","canvasLTM",800,700);
  canvasLTM->Divide(2,2);
  for (Int_t itype=0; itype<4; itype++){
    canvasLTM->cd(itype+1);
    TCut cutType=TString::Format("distType==0%d",itype).Data();
    tree->Draw("sqrt(npoints-2)*(meanLTM.fElements-mean)/sigma:vecLTM.fElements>>hisLTM(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgof");
    hisLTMTrunc[0]= (TH2*)(tree->GetHistogram()->Clone());
    tree->Draw("sqrt(npoints-2)*(meanLTMHisto.fElements-mean)/sigma:vecLTM.fElements>>hisLTMHisto(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[1]= (TH2*)tree->GetHistogram()->Clone();
    tree->Draw("sqrt(npoints-2)*(meanLTMHisto1.fElements-mean)/sigma:vecLTM.fElements>>hisLTMHist1(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[2]= (TH2*)tree->GetHistogram()->Clone();
    TLegend * legend = new TLegend(0.5,0.7,0.9,0.9,distNames[itype]);
    legend->SetBorderSize(0);
    for (Int_t ihis=0; ihis<3; ihis++){
      // MakeStat1D(TH2 * his, Int_t deltaBin, Double_t fraction, Int_t returnType, Int_t markerStyle, Int_t markerColor);
      grSigma[ihis]=TStatToolkit::MakeStat1D( hisLTMTrunc[ihis],0,0.99,5,kMarkers[ihis],kColors[ihis]);
      grSigma[ihis]->GetXaxis()->SetTitle("LTM fraction");
      grSigma[ihis]->GetYaxis()->SetTitle("#sigma_{meas}/sigma_{gauss}");
      if (ihis==0)grSigma[ihis]->Draw("alp");
      if (ihis>0)grSigma[ihis]->Draw("lp");
      legend->AddEntry(grSigma[ihis],names[ihis],"p");
    }
    legend->Draw();
  }
  canvasLTM->SaveAs("robustStatLTM_Performance.pdf");


}
