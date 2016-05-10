/*
  Test of the TLinearFitter functionality.
  Written to test "Tikchonov" regularization patch for ROOT TLinearFitter
  //
  To run the test code:

  .L $ALICE_ROOT/../src/STAT/test/TLinearFitterTest.C+
  TestRobustLinearFitWithRegularizationGener(500,25);
  TestRobustLinearFitWithRegularizationGenerAnalyze();
  
*/

#include "TGraphErrors.h"
#include "TLinearFitter.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTreeStream.h"

void TestRobustLinearFitWithRegularizationGenerAnalyze();
void TestRobustLinearFitWithRegularizationGener(Int_t nSamples=500, Int_t nPoints=25){
  //
  // Test for TLinearFitter
  // Compared Standard linear fit with Robust fit 
  //   for robust fit  (and also for normal fit) Tichonov  regularization  (ridge regression) applied
  // 
  // Robust fitting tested using PDF with admixture f(5%)  of outliers (sigma_{out}=20*sigma_{base})
  //  
  //
  // 1.)  Generate data sample
  // 2.)  Test invariants
  // 
  TTreeSRedirector *pcstream = new TTreeSRedirector("testRobustWithRegularization.root","recreate");
  //
  //
  //
  Double_t snoise=0.1;
  Double_t outlierFraction=0.05;
  //
  //
  for (Int_t isample=0; isample<nSamples; isample++){
    Double_t fdeltaT0=gRandom->Rndm()-0.5;
    Double_t fdeltaZ=gRandom->Rndm()-0.5;
    Double_t fdeltaDrift=gRandom->Rndm()-0.5;
    //
    // 1.1) Add points to design matrix
    TLinearFitter fitterPlane(3,"hyp2");  // y=[0]+[1]*sqrt(i)
    Double_t xxx[2]={0};  
    Int_t mPoints=gRandom->Poisson(nPoints);
    for (Int_t ipoint=0; ipoint<mPoints; ipoint++){
      Bool_t   Aside=(gRandom->Rndm()>0.5);  
      Double_t drift=(gRandom->Rndm());
      xxx[0]=(Aside)? -1.:1.;
      xxx[1]=drift;
      //
      Double_t y=fdeltaT0+(xxx[0]*fdeltaZ+xxx[1]*fdeltaDrift);
      if (gRandom->Rndm()>outlierFraction) {
	y+= gRandom->Gaus()*snoise;
      }else{
	y+= gRandom->Gaus()*snoise/outlierFraction;
      }
      fitterPlane.AddPoint(xxx,y,snoise);
    }
    // 
    // 1.2) Export graphs   
    Double_t dummy[10];
    TGraphErrors * grFit[6];
    for (Int_t igr=0;igr<6; igr++) grFit[igr]= new TGraphErrors(10,dummy, dummy, dummy,dummy);
    //
    //
    for (Int_t type=0; type<=1; type++){
      for (Int_t ireg=0; ireg<10; ireg++){ 
	Float_t reg=(ireg*ireg+1);
	if (type==0) fitterPlane.Eval(reg);
	if (type==1) fitterPlane.EvalRobust(1-1.5*outlierFraction, reg);
	for (Int_t ipar=0;ipar<3; ipar++){      
	  TGraphErrors *gre = grFit[type*3+ipar];
	  gre->SetPoint(ireg, reg, fitterPlane.GetParameter(ipar));
	  gre->SetPointError(ireg,0,  fitterPlane.GetParError(ipar));
	}
      }
    }    
    // 1.3) Store data in trees
    //
    (*pcstream)<<"robustFit"<<
      "fdeltaT0="<<fdeltaT0<<          // deltaT0 shift
      "fdeltaZ="<<fdeltaZ<<            // deltaZ  shift
      "fdeltaDrift="<<fdeltaDrift<<    // delta drift shift
      "snoise="<<snoise<<              // noise
      "outlierFraction="<<outlierFraction<<
      "mPoints="<<mPoints<<
      // fit graphs 
      "grDef0.="<<grFit[0]<< // default param0
      "grRob0.="<<grFit[3]<< // robust param0
      "grDef1.="<<grFit[1]<< // default param1
      "grRob1.="<<grFit[4]<< // robust param1
      "grDef2.="<<grFit[2]<< // default param2
      "grRob2.="<<grFit[5]<< // robust param2
      "\n";
  }  
  delete pcstream; 
  TestRobustLinearFitWithRegularizationGenerAnalyze();
}
 

void TestRobustLinearFitWithRegularizationGenerAnalyze(){
  //
  //  Analyze results of robust and standard fit
  //  1.) Draw and save pull graph - compare robust and defaul mode
  //  2.) Make alarm in case of "bad performance" (to be added)
  //
  TFile *finput = TFile::Open("testRobustWithRegularization.root");
  TTree *robustFit= (TTree*)finput->Get("robustFit");
  robustFit->SetMarkerStyle(25);
  robustFit->SetMarkerSize(0.5);
  //
  //
  TH1 *hisPulls[6]={0};
  robustFit->Draw("(grRob0.fY[0]-fdeltaT0)/grDef0->GetErrorY(0)>>hisRob0(100,-20,20)","","goff");
  hisPulls[0]=robustFit->GetHistogram();
  robustFit->Draw("(grDef0.fY[0]-fdeltaT0)/grDef0->GetErrorY(0)>>hisDef0(100,-20,20)","","goff");
  hisPulls[1]=robustFit->GetHistogram();

  robustFit->Draw("(grRob1.fY[0]-fdeltaZ)/grDef1->GetErrorY(0)>>hisRob1(100,-20,20)","","goff");
  hisPulls[2]=robustFit->GetHistogram();
  robustFit->Draw("(grDef1.fY[0]-fdeltaZ)/grDef1->GetErrorY(0)>>hisDef1(100,-20,20)","","goff");
  hisPulls[3]=robustFit->GetHistogram();


  TCanvas *canvasRobustFit = new TCanvas("canvasRobustFit","canvasRobustFit",1000,800);
  canvasRobustFit->Divide(2,2);
  //
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TLatex latex;
  const char * desc[4]={"Robust fit #Delta_{T0}/#Sigma", "Default fit #Delta_{T0}/#Sigma", "Robust fit #Delta_{Z}/#Sigma", "Default fit #Delta_{Z}/#Sigma"};
  for (Int_t ipar=0; ipar<4; ipar++){
    canvasRobustFit->cd(1+ipar)->SetLogy();
    hisPulls[ipar]->GetXaxis()->SetTitle("pull");
    hisPulls[ipar]->Fit("gaus","","");
    Double_t sigma = hisPulls[ipar]->GetFunction("gaus")->GetParameter(2);
    Double_t rms = hisPulls[ipar]->GetRMS();
    printf("%d\t%f\t%f\t%f\n", ipar,sigma,rms, rms/sigma);
    latex.DrawLatexNDC(0.11,0.7,desc[ipar]);
  }
  canvasRobustFit->SaveAs("TestRobustLinearFitWithRegularizationGenerAnalyze.png");


}
