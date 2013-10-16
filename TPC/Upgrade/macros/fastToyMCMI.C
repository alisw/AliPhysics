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
#include "TLegend.h"
#include "AliGeomManager.h"

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
  AliGeomManager::LoadGeometry();
  //
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
    Double_t cv[21]={0}; 
    Double_t covarSeed[15]={0};
    for (Int_t i=0; i<15; i++) covarSeed[i]=0;
    covarSeed[0]=0.1;
    covarSeed[2]=0.1;
    covarSeed[5]=0.001;
    covarSeed[9]=0.001;
    covarSeed[14]=0.02;
    //
    Double_t vertex[3]={0,0,0};
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
      Double_t pt = 2./(gRandom->Rndm()+0.00001);   
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
      Double_t pxyz[3];
      pxyz[0]=pt*TMath::Cos(phi);
      pxyz[1]=pt*TMath::Sin(phi);
      pxyz[2]=pt*TMath::Tan(theta);
      Short_t psign=(gRandom->Rndm()>0.5)? -1.:1.;
      AliExternalTrackParam *track= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
      //
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
      }
      //
      //  estimate q/pt resolution for the ITS+TPC, ITS+TPC+TRD and ITS+TRD scenario
      //
      AliExternalTrackParam * param = new AliExternalTrackParam(*track);
      Double_t *covar = (Double_t*)  param->GetCovariance();
      for (Int_t i=0; i<15; i++) covar[i]=covarSeed[i];
      AliTrackerBase::PropagateTrackToBxByBz(param, 370,0.13,1,kFALSE);
      
      AliExternalTrackParam *trackITSTPC= new AliExternalTrackParam(*param);   
      AliExternalTrackParam *trackITSTPCTRD= new AliExternalTrackParam(*param);   
      AliExternalTrackParam *trackITSTRD= new AliExternalTrackParam(*param);   
      AliExternalTrackParam *trackTPCTRD= new AliExternalTrackParam(*param);   
      //
      Bool_t tStatus=kTRUE;
      for (Int_t idet=2;idet>=0; idet--){
	Int_t nlayers=7;
	if (idet==1) nlayers=159;
	if (idet==2) nlayers=6;
	for (Int_t ilayer=nlayers; ilayer>=0; ilayer--){
	  Double_t rlayer=245.-ilayer;
	  Double_t slayer=0.1;
	  if (idet==0) {rlayer=rpoints[ilayer];  slayer=spoints[ilayer];}
	  if (idet==2) {rlayer=rpoints[ilayer+7]; slayer=spoints[ilayer+7];}
	  param->PropagateTo(rlayer,bz);
	  tStatus&=!AliTrackerBase::PropagateTrackToBxByBz(param, rlayer,0.13,1,kFALSE);
	  tStatus&=!AliTrackerBase::PropagateTrackToBxByBz(trackITSTPC, rlayer,0.13,1,kFALSE);
	  tStatus&=!AliTrackerBase::PropagateTrackToBxByBz(trackITSTPCTRD, rlayer,0.13,1,kFALSE);
	  tStatus&=!AliTrackerBase::PropagateTrackToBxByBz(trackITSTRD, rlayer,0.13,1,kFALSE);
	  tStatus&=!AliTrackerBase::PropagateTrackToBxByBz(trackTPCTRD, rlayer,0.13,1,kFALSE);
	  //if (tStatus==kFALSE) break;
	  Double_t pointPos[2]={param->GetY(),param->GetZ()};
	  Double_t pointCov[3]= { slayer*slayer,0,slayer*slayer};
	  tStatus&=!trackITSTPCTRD->Update(pointPos,pointCov);
	  if (idet!=2) {
	    trackITSTPC->Update(pointPos,pointCov);
	  }
	  if (idet!=1) {
	    trackITSTRD->Update(pointPos,pointCov);
	  }
	  if (idet!=0) {
	    trackTPCTRD->Update(pointPos,pointCov);
	  }
	}
      }
      

      //vecITSErr0.Print();
      (*pcstream)<<"extrapol"<<
	"itrack="<<itrack<<
	"track.="<<track<<
	"vecR.="<<&vecR<<
	"vecITSErr0.="<<&vecITSErr0<<
	"vecITSTRDErr0.="<<&vecITSTRDErr0<<
	"tStatus="<<tStatus<<
	"trackITSTPC.="<<trackITSTPC<<
	"trackITSTPCTRD.="<<trackITSTPCTRD<<
	"trackITSTRD.="<<trackITSTRD<<
	"trackTPCTRD.="<<trackTPCTRD<<
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
  TLegend * legend  = new TLegend(0.11,0.6,0.5,0.89,"Track residuals at TPC (q/p_{t}=0)");
  legend->AddEntry(grITS0,"ITS extrapolation","p");
  legend->AddEntry(grITSTRD0,"ITS-TRD interpolation","p");
  legend->Draw();
  //
  //
  //
  TCut cut="abs(track.fP[4])<0.25";
  TGraphErrors * grITSTPCqPt    = TStatToolkit::MakeGraphErrors(tree,"sqrt(trackITSTPC.fC[14]):abs(track.fP[4])",cut,20,1,0.8);
  TGraphErrors * grITSTRDqPt    = TStatToolkit::MakeGraphErrors(tree,"sqrt(trackITSTRD.fC[14]):abs(track.fP[4])",cut,21,4,0.8);
  TGraphErrors * grITSTPCTRDqPt = TStatToolkit::MakeGraphErrors(tree,"sqrt(trackITSTPCTRD.fC[14]):abs(track.fP[4])",cut,24,2,0.8);
  TGraphErrors * grTPCTRDqPt    = TStatToolkit::MakeGraphErrors(tree,"sqrt(trackTPCTRD.fC[14]):abs(track.fP[4])",cut,25,3,0.8);
  grITSTPCqPt->GetXaxis()->SetTitle("q/p_{t} (c/GeV)");
  grITSTPCqPt->GetYaxis()->SetTitle("#sigma_{q/p_{t}} (c/GeV)");
  grITSTPCqPt->SetMaximum(0.003);
  TCanvas * canvasResol = new TCanvas("canvasResol","canvasResol",800,600);
  canvasResol->cd();
  canvasResol->SetLeftMargin(0.15);
  grITSTPCqPt->GetYaxis()->SetTitleOffset(1.3);

  {
    grITSTPCqPt->Draw("ap");
    grITSTRDqPt->Draw("p");
    grITSTPCTRDqPt->Draw("p");
    grTPCTRDqPt->Draw("p");
    TLegend * legendQpt  = new TLegend(0.41,0.11,0.89,0.4,"q/p_{t} resolution (from covariance matrix) ");
    legendQpt->AddEntry(grITSTPCqPt,"ITS+TPC","p");
    legendQpt->AddEntry(grITSTRDqPt,"ITS+TRD","p");
    legendQpt->AddEntry(grITSTPCTRDqPt,"ITS+TPC+TRD","p");
    legendQpt->AddEntry(grTPCTRDqPt,"TPC+TRD","p");
    legendQpt->Draw();
  }

}
