/*
  Fast toy MC for differnt purposes.
  testExtrapolation() - to show track extrapolation errror
  testdEdxGEM();  -to show the dEdx respolution detoraition as function of the electron trasnparency
 

  .L $ALICE_ROOT/TPC/Upgrade/macros/fastToyMCMI.C+
  testExtrapolation();

*/
#include "TMath.h"
#include "TVectorD.h"
#include "TLinearFitter.h"
#include "TGeoGlobalMagField.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TStatToolkit.h"

#include "TTreeStream.h"
#include "AliExternalTrackParam.h"
#include "AliCDBManager.h"
#include "AliMagF.h"
#include "AliTrackerBase.h"
#include "TAxis.h"

void testdEdxGEM(){
}

void testExtrapolation(const Int_t ntracks=10){
  //
  // check the extrapolation error
  // 2 scenarios 
  //     - ITS extrapolation
  //     - ITS + TRD interpolation
  //
  //
  const char *  ocdbpath = "local://$ALICE_ROOT/OCDB/";  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,       AliMagF::kBeamTypepp, 2.76/2.));
   Double_t   bz=AliTrackerBase::GetBz(); 

  Double_t rpoints[13]={2.2,   2.8,   3.6,   20,    22,    41,     43,       300,315,330,345,360,375};
  Double_t spoints[13]={0.0004,0.0004,0.0004,0.0004,0.0004,0.0004, 0.0004,   0.02,0.02,0.02,0.02,0.02,0.02}; // ITS layers R poition (http://arxiv.org/pdf/1304.1306v3.pdf - pixel scenario) 
  //
  // 
  //
  const Int_t nbins=40;  
  TFile * f = TFile::Open("testExtrapolationErr.root","update");
  if (f->Get("extrapol")==0){
    delete f;
    f=0;
    TTreeSRedirector *pcstream = new TTreeSRedirector("testExtrapolationErr.root","update");
    //
    TVectorD vecR(nbins);
    TVectorD vecITSErr0(nbins);
    TVectorD vecITSTRDErr0(nbins);
    
    for (Int_t itrack=0; itrack<ntracks; itrack++){
      //
      //Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
      Double_t phi = 0;// gRandom->Uniform(0.0, 2*TMath::Pi());
      Double_t eta = gRandom->Uniform(-1, 1);
      Double_t pt = 1/(gRandom->Rndm()*2+0.00001);   
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
      Double_t pxyz[3];
      pxyz[0]=pt*TMath::Cos(phi);
      pxyz[1]=pt*TMath::Sin(phi);
      pxyz[2]=pt*TMath::Tan(theta);
      Short_t psign=(gRandom->Rndm()>0.5)? -1.:1.;
      Double_t cv[21]={0}; 
      Double_t vertex[3]={0,0,0};
      AliExternalTrackParam *track= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
      track->Rotate(phi);
      //
      // 
      for (Int_t irbin=0; irbin<nbins; irbin++){
	Double_t x0 =85+Double_t(irbin)*(245.-85.)/Double_t(nbins);
	Double_t x;
	//
	// parabolic fit
	//
	TLinearFitter fitterITS(3,"pol2");
	TLinearFitter fitterITSTRD(3,"pol2");
	vecR[irbin]=x0;
	for (Int_t i=0; i<13; i++) { 
	  Double_t y = 0;
	  if (track->GetYAt(rpoints[irbin],bz,y)){
	    x=rpoints[i]-x0; 
	    if (i<7) fitterITS.AddPoint(&x,y,spoints[i]);
	    fitterITSTRD.AddPoint(&x,y,spoints[i]);      
	  }
	}
	fitterITS.Eval();
	fitterITSTRD.Eval();
	vecITSErr0[irbin]=fitterITS.GetParError(0);
	vecITSTRDErr0[irbin]=fitterITSTRD.GetParError(0);
	//
	//
      }
      //vecITSErr0.Print();
      (*pcstream)<<"extrapol"<<
	"itrack="<<itrack<<
	//"track.="<<track<<
	"vecR.="<<&vecR<<
	"vecITSErr0.="<<&vecITSErr0<<
	"vecITSTRDErr0.="<<&vecITSTRDErr0<<
	"\n";
    }
    delete pcstream;
  }
  delete f;
  f = TFile::Open("testExtrapolationErr.root","update");
  TTree * tree = (TTree*)f->Get("extrapol");
  //
  TGraphErrors * grITS0    = TStatToolkit::MakeGraphErrors(tree,"10*vecITSErr0.fElements:vecR.fElements","",25,4,0);
  TGraphErrors * grITSTRD0 = TStatToolkit::MakeGraphErrors(tree,"10*vecITSTRDErr0.fElements:vecR.fElements","",21,2,0);
  //
  grITS0->GetXaxis()->SetTitle("radius (cm)");
  grITS0->GetYaxis()->SetTitle("#sigma_{r#phi} (mm)");
  grITS0->Draw("ap");
  grITSTRD0->Draw("p");
  TLegend * legend  = new TLegend(0.11,0.5,0.5,0.89,"Track residuals at TPC (q/p_{t}=0)");
  legend->AddEntry(grITS0,"ITS extrapolation","p");
  legend->AddEntry(grITSTRD0,"ITS extrapolation","p");
  legend->Draw();


}
