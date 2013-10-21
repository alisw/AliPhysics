/*
  Fast toy MC for differnt purposes.
  testExtrapolation() - to show track extrapolation errror
  testdEdxGEM();  -to show the dEdx respolution detoraition as function of the electron trasnparency
 
  .x $HOME/NimStyle.C
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
#include "TCut.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"

void testdEdxGEM(const Int_t ntracks=10000, Double_t clNoise=2, Double_t corrNoiseAdd=0.15, Double_t sigmaLandau=0.26){
  //
  // test dEdx performance as function of the electron transparency.
  //
  // For simulation standard gRandom->Landau() generator with mean 15 is used
  // Charge between pad-rows are corelated via diffusion - qm = (0.15*q-1+0.7*q0+0.15*q1)
  // Electrons are randomly absorbed dependending on the transparency parameter
  // 
  // Parameters:
  //   clNoise       - noise integrated by cluster (snoise~1 x 4 bins ~ 2)
  //   corrNoiseAdd  - additional correlated noise to be added
  //   sigmaLandau   - relative sigma of Landau distirbution
  // Setup emulation
  //   pp   setup testdEdxGEM(20000,2,0.15,0.26)
  //   PbPb setup testdEdxGEM(20000,4,0.5,0.26)
  //   PbPb setup testdEdxGEM(20000,4,0.6,0.26)
  //   PbPb setup testdEdxGEM(20000,4,0.7,0.26)
  // 
 
  TTreeSRedirector *pcstream = new TTreeSRedirector("testdEdxResolution.root","update");
  TH2F * phisdEdx = (TH2F*)pcstream->GetFile()->Get("hisdEdx");
  if (!phisdEdx){

    TVectorD qVector(160);
    TVectorD qVectorCorr(160);
    TVectorD qVectorAcc(160);
    TVectorD qVectorSorted(160);
    Int_t indexes[160];
    Int_t ntrans=20;
    TVectorD vecdEdx(ntrans);
    TVectorD vecT(ntrans);
    //
    
    for (Int_t itrack=0; itrack<ntracks; itrack++){
      Double_t meanEl=15*(1.+1*gRandom->Rndm());
      Double_t sigmaEl=sigmaLandau*meanEl*TMath::Power(15./meanEl,0.25);
      for (Int_t irow=0; irow<160; irow++){
	qVector[irow]=gRandom->Landau(meanEl, sigmaEl);      
      }
      qVectorCorr[0]=qVector[0];
      qVectorCorr[158]=qVector[158];
      for (Int_t irow=1; irow<159; irow++){  //corralte measurement via diffusion
	qVectorCorr[irow]= 0.15*(qVector[irow-1]+ qVector[irow+1])+0.7*qVector[irow];
      }
      
      for (Int_t itrans=0; itrans<ntrans; itrans++){
	Double_t transparency=(itrans+1.)/ntrans;
	vecT[itrans]=transparency;
	for (Int_t irow=0; irow<160; irow++) {
	  qVectorAcc[irow]=0;
	  for (Int_t iel=0; iel<TMath::Nint(qVectorCorr[irow]); iel++) {
	    if (gRandom->Rndm()<transparency) 	qVectorAcc[irow]+=1./transparency;	
	  }
	}      
	for (Int_t irow=0; irow<160; irow++) {
	  qVectorAcc[irow]+=gRandom->Gaus(0, clNoise);
	  if (qVectorAcc[irow]<0) qVectorAcc[irow]=0;
	}
	TMath::Sort(160,qVectorAcc.GetMatrixArray(), indexes, kFALSE);
	for (Int_t irow=0; irow<160; irow++) {
	  qVectorSorted[irow]=qVectorAcc[indexes[irow]];
	}
	vecdEdx[itrans]=TMath::Mean(0.6*160.,	qVectorSorted.GetMatrixArray());    
	vecdEdx[itrans]+=gRandom->Gaus(0, corrNoiseAdd);
      }
      (*pcstream)<<"dEdx"<<
	"itrack="<<itrack<<              // itrack 
	"meanEl="<<meanEl<<              // 
	"sigmaEl="<<sigmaEl<<            //
	"vecT.="<<&vecT<<                //
	"vecdEdx.="<<&vecdEdx<<
	"qVector.="<<&qVector<<
	"\n";
    }
    TTree * tree = (TTree*)(pcstream->GetFile()->Get("dEdx"));
    tree->Draw("vecdEdx.fElements/(meanEl):vecT.fElements>>hisdEdx(16,0.199,1,100,0.70,1.50)","meanEl<100","colz");    
    phisdEdx= (TH2F*)tree->GetHistogram()->Clone();    
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    phisdEdx->Write("hisdEdx");
  }
  delete pcstream;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  pcstream = new TTreeSRedirector("testdEdxResolution.root","update"); 
  phisdEdx = (TH2F*)pcstream->GetFile()->Get("hisdEdx");
  TObjArray *arrFit = new TObjArray(3);
  phisdEdx->FitSlicesY(0,0,-1,0,"QNR",arrFit);
  TH1D * hisRes = (TH1D*)arrFit->At(2);
  hisRes->Divide((TH1D*)arrFit->At(1)); 
  hisRes->Scale(100);
  hisRes->SetTitle("dEdx resolution");
  hisRes->SetMarkerStyle(21);
  hisRes->GetXaxis()->SetTitle("Eff. electron  transparency");
  hisRes->GetYaxis()->SetTitle("#sigma_{dEdx}/dEdx (%)");
  hisRes->SetMinimum(4);
  hisRes->SetMaximum(8);
  //
  TF1 * ftrans = new TF1("ftrans","[0]*x**(-(abs([1])+0.000001))");
  ftrans->SetParName(0,"#sigma_{T1}");
  ftrans->SetParName(1,"#sigma slope");
  ftrans->SetParameters(0.05,1,0.05);
  hisRes->SetMarkerStyle(21);
  hisRes->SetMarkerSize(0.9);
  TCanvas * canvasTrans = new TCanvas("canvasTrans", "canvasTrans",600,500);
  canvasTrans->SetTicks(1,1);
  hisRes->Fit(ftrans,"","",0.2,1);
  canvasTrans->SaveAs("canvasElTrans.pdf");
  //TH
  
  
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
