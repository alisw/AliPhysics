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
#include "TLinearFitter.h"
#include "TMultiGraph.h"
#include "TObjArray.h"

TObjArray arrayFit(3);
Int_t kMarkers[10]={25,24,20,21,22};
Int_t kColors[10]={1,2,4,3,6};
const char * names[10] = {"LTM","LThisto","LTMhistoPar","Gaus fit", "Gaus in range"};
const char * distNames[10] = {"Gauss","Gauss+flat","Gauss(m)+0.2*Gauss(m+5#sigma)","Gauss(m)+0.4*Gauss(m+5#sigma)"};
const char * normNames[TStatToolkit::kNNormType] = { "L1", "L2", "LP", "Max", "Hamming" };


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
  TVectorD meanLTMHistoPar(6);
  TVectorD sigmaLTMHistoPar(6);
  TVectorD meanLTMHistoGausLim(6);
  TVectorD sigmaLTMHistoGausLim(6);
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
	if (gRandom->Rndm()>0.4) {                     //  gauss + second gaus 5 sigma away
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
      //
      // graph fit
      //
      Int_t binC=histo.GetXaxis()->FindBin(paramLTM[1]);
      Int_t bin0=TMath::Max(histo.GetXaxis()->FindBin(paramLTM[1]-5*paramLTM[2]),1);
      Int_t bin1=TMath::Min(histo.GetXaxis()->FindBin(paramLTM[1]+5*paramLTM[2]),histo.GetXaxis()->GetNbins()) ;
      if ( (bin0-bin1) <3){
	bin0=TMath::Max(binC-2,1);
	bin1=TMath::Min(binC+2, histo.GetXaxis()->GetNbins());
      }
      //
      /* test input:
	 TH1F histo("aaa","aaa",100,-10,10);for (Int_t i=0; i<10000; i++) histo.Fill(gRandom->Gaus(2,1));
      */
      
      TLinearFitter fitter(3,"hyp2");
      for (Int_t ibin=bin0; ibin<bin1; ibin++){
	Double_t x=histo.GetXaxis()->GetBinCenter(ibin);
	Double_t y=histo.GetBinContent(ibin);
	Double_t sy=histo.GetBinError(ibin);
	Double_t xxx[3]={x,x*x};
	if (y>3.*sy){
	  fitter.AddPoint(xxx,TMath::Log(y),sy/y);
	}
      }
      if (fitter.GetNpoints()>3){
	fitter.Eval();
	// f(x)=[0]+[1]*x+[2]*x*x;
	// f(x)=s*(x0-x0)^2+s
	if (TMath::Abs(fitter.GetParameter(2))<1e-12) continue;
	Double_t parSigma=TMath::Sqrt(2*TMath::Abs(fitter.GetParameter(2)));
	Double_t parMean=-fitter.GetParameter(1)/(2.*fitter.GetParameter(2));  
	histo.Fit(&fg,"QN","QN", paramLTM[1]-5*paramLTM[2], paramLTM[1]+5*paramLTM[2]);
	meanLTMHistoPar[iltm]=parMean;
	sigmaLTMHistoPar[iltm]=parSigma;
	meanLTMHistoGausLim[iltm]=fg.GetParameter(1);
	sigmaLTMHistoGausLim[iltm]=fg.GetParameter(2);
      }else{
	meanLTMHistoPar[iltm]=paramLTM[1];
	sigmaLTMHistoPar[iltm]=paramLTM[2];
	meanLTMHistoGausLim[iltm]=meanG;
	sigmaLTMHistoGausLim[iltm]=rmsG;	
      }

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
      "meanLTMHistoPar.="<<&meanLTMHistoPar<<
      "meanLTMHistoGausLim.="<<&meanLTMHistoGausLim<<
      //
      "sigmaLTM.="<<&sigmaLTM<<
      "sigmaLTMHisto.="<<&sigmaLTMHisto<<
      "sigmaLTMHistoPar.="<<&sigmaLTMHistoPar<<      
      "sigmaLTMHistoGausLim.="<<&sigmaLTMHistoGausLim<<      
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
  //TH1 * hisResol[10]={0};
  TGraphErrors * grSigma[10]={0};
  TCanvas *canvasLTM= new TCanvas("canvasLTM","canvasLTM",800,700);
  canvasLTM->Divide(2,2,0,0);
  for (Int_t itype=0; itype<4; itype++){
    canvasLTM->cd(itype+1);
    TCut cutType=TString::Format("distType==0%d",itype).Data();
    tree->Draw("sqrt(npoints-2)*(meanLTM.fElements-mean)/sigma:vecLTM.fElements>>hisLTM(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgof");
    hisLTMTrunc[0]= (TH2*)(tree->GetHistogram()->Clone());
    tree->Draw("sqrt(npoints-2)*(meanLTMHisto.fElements-mean)/sigma:vecLTM.fElements>>hisLTMHisto(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[1]= (TH2*)tree->GetHistogram()->Clone();
    tree->Draw("sqrt(npoints-2)*(meanLTMHistoPar.fElements-mean)/sigma:vecLTM.fElements>>hisLTMHist1(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[2]= (TH2*)tree->GetHistogram()->Clone();
    tree->Draw("sqrt(npoints-2)*(meanG-mean)/sigma:vecLTM.fElements>>hisLTMHist1(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[3]= (TH2*)tree->GetHistogram()->Clone();
    tree->Draw("sqrt(npoints-2)*(meanLTMHistoGausLim.fElements-mean)/sigma:vecLTM.fElements>>hisLTMHist1(6,0.45,1.05,100,-10,10)",cutType+"npoints>50&&sigma>0.1","colzgoff");
    hisLTMTrunc[4]= (TH2*)tree->GetHistogram()->Clone();
    TLegend * legend = new TLegend(0.5,0.7,0.88,0.88,distNames[itype]);
    legend->SetBorderSize(0);
    for (Int_t ihis=0; ihis<5; ihis++){
      // MakeStat1D(TH2 * his, Int_t deltaBin, Double_t fraction, Int_t returnType, Int_t markerStyle, Int_t markerColor);
      grSigma[ihis]=TStatToolkit::MakeStat1D( hisLTMTrunc[ihis],0,0.99,1,kMarkers[ihis],kColors[ihis]);
      grSigma[ihis]->SetMaximum(5.0);
      grSigma[ihis]->SetMinimum(0.5);
      grSigma[ihis]->GetXaxis()->SetTitle("LTM fraction");
      grSigma[ihis]->GetYaxis()->SetTitle("#sigma_{meas}/#sigma_{Theor.}");
      if (ihis==0)grSigma[ihis]->Draw("alp");
      if (ihis>0)grSigma[ihis]->Draw("lp");
      legend->AddEntry(grSigma[ihis],names[ihis],"p");
    }
    legend->Draw();
  }
  canvasLTM->SaveAs("robustStatLTM_Performance.pdf");


}

void TestNormTypes() 
{
   //input arrays
   const Int_t nValues=5;

   Double_t v1=0.;
   Double_t v2=0.;
   Double_t v3=0.;

   TTree t1("t1","t1");
   TTree t2("t2","t2");

   t1.Branch("v1", &v1);
   t1.Branch("v2", &v2);
   t1.Branch("v3", &v3);
   t2.Branch("v1", &v1);
   t2.Branch("v2", &v2);
   t2.Branch("v3", &v3);

   for (Int_t ival=0; ival<nValues; ++ival) {
     v1=ival;
     v2=v1*v1;
     v3=2*ival; // same in both trees
     t1.Fill();

     v1=ival+1;
     v2=v1*v1;
     t2.Fill();
   }

   t1.AddFriend(&t2,"t2");

   // store norm for all type and different variable combinations
   const Int_t nVec = 6;
   const char* statGraphNames[nVec] = {"T1 v1", "T1 v2", "T1 v1&v2", "T2 v1", "T1-T2 v1", "T1-T2 v3"};
   TVectorD statTypesT1v1(TStatToolkit::kNNormType);
   TVectorD statTypesT1v2(TStatToolkit::kNNormType);
   TVectorD statTypesT1v1v2(TStatToolkit::kNNormType);
   TVectorD statTypesT2v1(TStatToolkit::kNNormType);
   TVectorD statTypesT1T2v1Diff(TStatToolkit::kNNormType);
   TVectorD statTypesT1T2v3Diff(TStatToolkit::kNNormType);
   TVectorD vx(TStatToolkit::kNNormType);

   TObjArray arrVect;
   arrVect.Add(&statTypesT1v1);
   arrVect.Add(&statTypesT1v2);
   arrVect.Add(&statTypesT1v1v2);
   arrVect.Add(&statTypesT2v1);
   arrVect.Add(&statTypesT1T2v1Diff);
   arrVect.Add(&statTypesT1T2v3Diff);

   // fill distance variables
   for (Int_t inormType=0; inormType<TStatToolkit::kNNormType; ++inormType) {
     vx(inormType) = inormType;
     statTypesT1v1(inormType)         = TStatToolkit::GetDistance(&t1, "v1",          "", (TStatToolkit::ENormType)inormType, 3);
     statTypesT1v2(inormType)         = TStatToolkit::GetDistance(&t1, "v2",          "", (TStatToolkit::ENormType)inormType, 3);
     statTypesT1v1v2(inormType)       = TStatToolkit::GetDistance(&t1, "v1:v2",       "", (TStatToolkit::ENormType)inormType, 3);
     statTypesT2v1(inormType)         = TStatToolkit::GetDistance(&t2, "v1",          "", (TStatToolkit::ENormType)inormType, 3);
     statTypesT1T2v1Diff(inormType)   = TStatToolkit::GetDistance(&t1, "v1-t2.v1",    "", (TStatToolkit::ENormType)inormType, 3);
     statTypesT1T2v3Diff(inormType)   = TStatToolkit::GetDistance(&t1, "v3-t2.v3",    "", (TStatToolkit::ENormType)inormType, 3);
   }

   // create graphs and draw
   TMultiGraph *mgr = new TMultiGraph;

   TLegend *leg = new TLegend(0.6, 0.6, 0.9, .9);
   leg->SetFillColor(10);
   leg->SetBorderSize(1);

   for (Int_t igraph=0; igraph<arrVect.GetEntriesFast(); ++igraph) {
     TVectorD vx2(vx);
     vx2+=1./(nVec+1)*(igraph+1);
     TGraph *gr = new TGraph(vx2, *((TVectorD*)arrVect.At(igraph)));
     gr->SetMarkerColor(igraph+1);
     gr->SetMarkerSize(1);
     gr->SetMarkerStyle(21);
     gr->SetTitle(statGraphNames[igraph]);
     leg->AddEntry(gr, statGraphNames[igraph], "p");
     mgr->Add(gr);
   }

   TCanvas *c=new TCanvas("statSummary","statSummary");
   TH1F *hDummy = new TH1F("hDummy", ";norm type;norm", TStatToolkit::kNNormType, 0., (Double_t)TStatToolkit::kNNormType);
   hDummy->SetMaximum(70);
   hDummy->Draw();
   mgr->Draw("p");
   TAxis *xAxis = hDummy->GetXaxis();
   for (Int_t inormType=0; inormType<TStatToolkit::kNNormType; ++inormType) {
     xAxis->SetBinLabel(inormType+1, normNames[inormType]);
   }
   leg->Draw();
   c->Modified();
   c->Update();
}

