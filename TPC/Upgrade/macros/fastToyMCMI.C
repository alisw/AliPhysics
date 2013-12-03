/*
  Fast toy MC for different purposes. primary goal prepare the motivation figures for the TPC TDR and
  expalantion internal note.
  
  testExtrapolationError() - to show track extrapolation errror
  testdEdxGEM();  -to show the dEdx respolution detoriation as function of the electron trasnparency
 

  .x $HOME/NimStyle.C
  .L $ALICE_ROOT/TPC/Upgrade/macros/fastToyMCMI.C+
 
*/
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
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



void testdEdxGEM(const Int_t ntracks=10000, Double_t clNoise=2, Double_t corrNoiseAdd=0.15, Double_t sigmaLandau=0.26);
void testExtrapolationError(Int_t ntracks, Int_t useGeo=kFALSE, Int_t seed=0);


void fastToyMCMI(Int_t action=1, Float_t arg0=200, Float_t arg1=0, Float_t arg2=0){
  //
  //
  //
  if (action==1)  testExtrapolationError(arg0,arg1,arg2); 
}


void testdEdxGEM(const Int_t ntracks, Double_t clNoise, Double_t corrNoiseAdd, Double_t sigmaLandau){
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

void testExtrapolationError(Int_t ntracks, Int_t useGeo, Int_t seed){
  //
  // check the extrapolation error
  // 2 scenarios 
  //     - ITS extrapolation
  //     - ITS + TRD interpolation
  //
  //
  gRandom->SetSeed(seed);
  const char *  ocdbpath = "local://$ALICE_ROOT/OCDB/";  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,       AliMagF::kBeamTypepp, 2.76/2.));
  Double_t   bz=AliTrackerBase::GetBz(); 
  if ( useGeo)   AliGeomManager::LoadGeometry("geometry.root");
  if (!useGeo)   AliGeomManager::LoadGeometry();
  //
  //  Double_t rpoints[13]={2.2,   2.8,   3.6,   20,    22,    41,     43,       300,315,330,345,360,375};
  Double_t spoints[13]={0.0004,0.0004,0.0004,0.0004,0.0004,0.0004, 0.0004,   0.02,0.02,0.02,0.02,0.02,0.02}; // ITS layers R poition (http://arxiv.org/pdf/1304.1306v3.pdf - pixel scenario) 
  //
  Double_t rpoints[13]={2.32,   3.13,   3.91,   19.41,    24.71,    35.33,     40.53,       300,315,330,345,360,375};

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
    covarSeed[14]=0.03;
    //
    Double_t vertex[3]={0,0,0};
    TTreeSRedirector *pcstream = new TTreeSRedirector("testExtrapolationErr.root","update");
    //
    TVectorD vecR(nbins);
    TVectorD vecITSErr0(nbins);
    TVectorD vecITSTRDErr0(nbins);
    //
    for (Int_t itrack=0; itrack<ntracks; itrack++){
      //
      //Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
      if (itrack%20==0) printf("processing\t%d\n", itrack);
      Double_t phi = (gRandom->Rndm()-1)*TMath::Pi()/18;
      Double_t eta = gRandom->Uniform(-1, 1);
      Double_t pt = 1./(gRandom->Rndm()+0.00001);   
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
      Double_t pxyz[3];
      pxyz[0]=pt*TMath::Cos(phi);
      pxyz[1]=pt*TMath::Sin(phi);
      pxyz[2]=pt*TMath::Tan(theta);
      Short_t psign=(gRandom->Rndm()>0.5)? -1.:1.;
      AliExternalTrackParam *track= new AliExternalTrackParam(vertex, pxyz, cv, psign);   
      Double_t alpha=TMath::ATan2(pxyz[1],pxyz[0]);
      track->Rotate(alpha);
      //
      // 0.) Estimate the error using the ITS extrapolation and ITS+TRD interpolation - neglecting the MS -inf. momanta tracks
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
      //  1.) estimate q/pt resolution for the ITS+TPC, ITS+TPC+TRD and ITS+TRD scenario
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
      //
      // 2.) Estimate propagation error at given radius
      //
      //  2.a) Fit ITS and TRD track (gauss smeared point error 0 not MS used)
      AliExternalTrackParam *trackITS= new AliExternalTrackParam(*track);
      AliExternalTrackParam *trackITS0= new AliExternalTrackParam(*track);
      covar = (Double_t*)  trackITS->GetCovariance();
      for (Int_t i=0; i<15; i++) covar[i]=covarSeed[i];
      AliExternalTrackParam *trackTRD= new AliExternalTrackParam(*param);      
      AliExternalTrackParam *trackTRD0= new AliExternalTrackParam(*param);      
      {for (Int_t ilayer=0; ilayer<7; ilayer++){ 
	  AliTrackerBase::PropagateTrackToBxByBz(trackITS, rpoints[ilayer],0.13,1,kFALSE);
	  AliTrackerBase::PropagateTrackToBxByBz(trackITS0, rpoints[ilayer],0.13,1,kFALSE);
	  Double_t pointPos[2]={trackITS0->GetY()+gRandom->Gaus(0,spoints[ilayer]),trackITS0->GetZ()+gRandom->Gaus(0,spoints[ilayer])};
	  Double_t pointCov[3]= {spoints[ilayer]*spoints[ilayer],0,spoints[ilayer]*spoints[ilayer]};
	  trackITS->Update(pointPos,pointCov);
	}}
      {for (Int_t ilayer=5; ilayer>=0; ilayer--){ 
	  AliTrackerBase::PropagateTrackToBxByBz(trackTRD, rpoints[ilayer+7],0.13,1,kFALSE);
	  AliTrackerBase::PropagateTrackToBxByBz(trackTRD0, rpoints[ilayer+7],0.13,1,kFALSE);
	  Double_t pointPos[2]={trackTRD0->GetY()+gRandom->Gaus(0,spoints[ilayer+7]),trackTRD0->GetZ()+gRandom->Gaus(0,spoints[ilayer+7])};
	  Double_t pointCov[3]= {spoints[ilayer+7]*spoints[ilayer+7],0,spoints[ilayer+7]*spoints[ilayer+7]};
	  trackTRD->Update(pointPos,pointCov);
	}}
      //
      // 2.b) get ITS extrapolation and TRD interpolation errors at random radisues positions
      //
      const Int_t npoints=20;
      TVectorD vecRKalman(npoints);
      TVectorD vecITSDeltaY(npoints),    vecITSTRDDeltaY(npoints);
      TVectorD vecITSErrY(npoints),      vecITSTRDErrY(npoints);
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
	Double_t rLayer= 85.+gRandom->Rndm()*(245.-85.);
	AliExternalTrackParam *trackITSP= new AliExternalTrackParam(*trackITS);
	AliExternalTrackParam *trackITSP0= new AliExternalTrackParam(*trackITS0);
	AliExternalTrackParam *trackTRDP= new AliExternalTrackParam(*trackTRD);
	AliTrackerBase::PropagateTrackToBxByBz(trackITSP, rLayer,0.13,1,kFALSE);
	AliTrackerBase::PropagateTrackToBxByBz(trackITSP0, rLayer,0.13,1,kFALSE);
	AliTrackerBase::PropagateTrackToBxByBz(trackTRDP, rLayer,0.13,1,kFALSE); 
	vecRKalman[ipoint]=rLayer;
	trackTRDP->Rotate(trackITSP->GetAlpha());
	vecITSDeltaY[ipoint]=trackITSP->GetY()-trackITSP0->GetY();
	vecITSErrY[ipoint]=TMath::Sqrt(trackITSP->GetSigmaY2());
	AliTrackerBase::UpdateTrack(*trackITSP, *trackTRDP);
	vecITSTRDDeltaY[ipoint]=trackITSP->GetY()-trackITSP0->GetY();
	vecITSTRDErrY[ipoint]=TMath::Sqrt(trackITSP->GetSigmaY2());
      }
      //

      //vecITSErr0.Print();
      (*pcstream)<<"extrapol"<<
	"itrack="<<itrack<<
	"track.="<<track<<
	"vecR.="<<&vecR<<
	//
	"vecITSErr0.="<<&vecITSErr0<<
	"vecITSTRDErr0.="<<&vecITSTRDErr0<<
	"tStatus="<<tStatus<<
	"trackITSTPC.="<<trackITSTPC<<
	"trackITSTPCTRD.="<<trackITSTPCTRD<<
	"trackITSTRD.="<<trackITSTRD<<
	"trackTPCTRD.="<<trackTPCTRD<<
	//
	// kalman extapoltation einterplotaltion studies 
	// extrapolation and ITSTRDitnerpolation at different radiuses 
	"verRKalman.="<<&vecRKalman<<            
	"vecITSDeltaY.="<<&vecITSDeltaY<<
	"vecITSErrY.="<<&vecITSErrY<<
	"vecITSTRDDeltaY.="<<&vecITSTRDDeltaY<<
	"vecITSTRDErrY.="<<&vecITSTRDErrY<<
	"\n";
    }
    delete pcstream;
  }
  delete f;

  gStyle->SetOptTitle(0);
  f = TFile::Open("testExtrapolationErr.root","update");
  TTree * tree = (TTree*)f->Get("extrapol");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.3);

  //
  TGraphErrors * grITS0    = TStatToolkit::MakeGraphErrors(tree,"10*vecITSErr0.fElements:vecR.fElements","",25,4,0);
  TGraphErrors * grITSTRD0 = TStatToolkit::MakeGraphErrors(tree,"10*vecITSTRDErr0.fElements:vecR.fElements","",21,2,0);
  //
  grITS0->GetXaxis()->SetTitle("radius (cm)");
  grITS0->GetYaxis()->SetTitle("#sigma_{r#varphi} (mm)");
  grITS0->Draw("ap");
  grITSTRD0->Draw("p");
  TLegend * legend  = new TLegend(0.11,0.65,0.55,0.89,"Track residuals at TPC (q/p_{t}=0)");
  legend->AddEntry(grITS0,"ITS extrapolation","p");
  legend->AddEntry(grITSTRD0,"ITS-TRD interpolation","p");
  legend->SetFillColor(10);
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
  //
  //
  TCanvas *canvasResolution = new TCanvas("canvasResolution","canvasResolution",600,600);
  canvasResolution->Divide(1,3);
  //
  canvasResolution->cd(1);
  tree->Draw("vecITSErrY.fElements:verRKalman.fElements:abs(track.fP[4])","","colz");
  canvasResolution->cd(2);
  tree->Draw("vecITSTRDErrY.fElements:verRKalman.fElements:abs(track.fP[4])","","colz");
  canvasResolution->cd(3);
  tree->Draw("vecITSTRDErrY.fElements/vecITSErrY.fElements:verRKalman.fElements:abs(track.fP[4])","","colz");
  canvasResolution->SaveAs("canvasResolution.pdf");
  canvasResolution->SaveAs("canvasResolution.png");
  //
  //
  //
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;  
  //
  TCut cut="";
  TString  fstringITS="";
  fstringITS+="abs(track.fP[4])*(verRKalman.fElements-40)^2++";   
  fstringITS+="(verRKalman.fElements-40)^2++";   
  TString * fitErrITS = TStatToolkit::FitPlane(tree,"vecITSErrY.fElements:(0.1+track.fP[4]**2)", fstringITS.Data(),"verRKalman.fElements>0", chi2,npoints,param,covar,-1,0,180000 , kFALSE);
  tree->SetAlias("fitErrITS",fitErrITS->Data());
  fitErrITS->Tokenize("++")->Print();
  tree->Draw("vecITSErrY.fElements/fitErrITS:verRKalman.fElements:abs(track.fP[4])","","colz");
  //
  //
  TString  fstringITSTRD="";
  fstringITSTRD+="abs(track.fP[4])++";   
  fstringITSTRD+="abs(track.fP[4]*verRKalman.fElements)++";   
  fstringITSTRD+="abs(track.fP[4]*verRKalman.fElements**2)++";   
  for (Int_t iter=0; iter<3; iter++){
    TCut cutOut="1";
    if (iter>0) cutOut="abs(vecITSTRDErrY.fElements/fitErrITSTRD-1)<0.8";
    TString * fitErrITSTRD = TStatToolkit::FitPlane(tree,"vecITSTRDErrY.fElements:(0.1+track.fP[4]**2)", fstringITSTRD.Data(),"verRKalman.fElements>0"+cutOut, chi2,npoints,param,covar,-1,0,180000 , kFALSE);
    tree->SetAlias("fitErrITSTRD",fitErrITSTRD->Data());
    fitErrITSTRD->Tokenize("++")->Print();
  }
  tree->Draw("vecITSTRDErrY.fElements/fitErrITSTRD:verRKalman.fElements:abs(track.fP[4])","","colz");

  
  TCanvas *canvasResolutionFit = new TCanvas("canvasResolutionFit","canvasResolutionFit",700,700);
  canvasResolutionFit->Divide(1,3);
  canvasResolutionFit->cd(1);
  tree->Draw("fitErrITS:verRKalman.fElements:abs(track.fP[4])>>hITS","","colz"); 
  htemp = (TH2F*)gPad->GetPrimitive("hITS")->Clone();
  htemp->GetXaxis()->SetTitle("#it{r} (cm)");
  htemp->GetYaxis()->SetTitle("#sigma_{r#phi} (cm) ITS");
  htemp->GetZaxis()->SetTitle("q/#it{p_{t}} (#it{c}/GeV)");
  htemp->Draw("colz");
  canvasResolutionFit->cd(2);
  tree->Draw("fitErrITSTRD:verRKalman.fElements:abs(track.fP[4])>>hITSTRD","","colz");
  htemp = (TH2F*)gPad->GetPrimitive("hITSTRD")->Clone();
  htemp->GetXaxis()->SetTitle("#it{r} (cm)");
  htemp->GetYaxis()->SetTitle("#sigma_{r#phi} (cm) ITS+TRD");
  htemp->GetZaxis()->SetTitle("q/#it{p_{t}} (#it{c}/GeV)");
  htemp->Draw("colz");
  canvasResolutionFit->cd(3);
  tree->Draw("fitErrITSTRD:verRKalman.fElements:abs(track.fP[4])>>hITSTRD","","colz");
  htemp = (TH2F*)gPad->GetPrimitive("hITSTRD")->Clone();
  htemp->GetXaxis()->SetTitle("#it{r} (cm)");
  htemp->GetYaxis()->SetTitle("#sigma_{r#phi} (cm) ITS+TRD");
  htemp->GetZaxis()->SetTitle("q/#it{p_{t}} (#it{c}/GeV)");
  htemp->Draw("colz");
  canvasResolutionFit->SaveAs("canvasResolutionFit.pdf");
  canvasResolutionFit->SaveAs("canvasResolutionFit.png");

}


void DrawInerpolationResolution(){

  // resRphi = 0.004390 + oneOverPt*(-0.136403) + oneOverPt*radius*(0.002266) + oneOverPt*radius*radius*(-0.000006);

  //TF1 f1("f1","[0]+[1]*(
}

void MakeRobustFitTest(){
  //
  //
  //
  // 
  TH2F *   his2D   = new TH2F("his2D", "his2D",50,0,1., 100,-2.,2.);
  TH2F *   his2DW   = new TH2F("his2DW", "his2DW",50,0,1., 100,-2.,2.);
  Double_t probRND = 0.1;
  Int_t    ntracks = 20*50000;
  Int_t    nclusters = 16;
  Double_t sigmaCluster=0.1;
  //
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    Double_t x     = gRandom->Rndm();
    Double_t widthTrack = 0.1/(0.5+gRandom->Exp(0.5));
    Double_t y     = 0;
    y= gRandom->Gaus(0.0,widthTrack);
    Bool_t isRandom=gRandom->Rndm()<probRND;
    //if (gRandom->Rndm()<probRND) y= -2+4*gRandom->Rndm();  
    //
    Double_t sigmaTrack= TMath::Sqrt(sigmaCluster*sigmaCluster/nclusters+widthTrack*widthTrack);
    for (Int_t icl=0; icl<nclusters; icl++){
      his2D->Fill(x,y+gRandom->Gaus(0,sigmaCluster));
      his2DW->Fill(x,y+gRandom->Gaus(0,sigmaCluster),1/sigmaTrack);
    }
  }
  {
    his2D->Draw("colz");
    his2DW->Draw("colz");
    his2D->FitSlicesY();
    his2DW->FitSlicesY();    
    his2D_1->Draw();
    his2DW_1->Draw("same");
  }
  TMath::RMS(50, &(his2D_1->GetArray()[1]));

  TH1D * phis1D = his2D->ProjectionY("aaa",5,5);
  Int_t nbinsY= phis1D->GetXaxis()->GetNbins();
  TVectorD vector(nbinsY, &(phis1D->GetArray()[1]));

  TStatToolkit::MakeStat1D((TH2*)his2D, 1,0.8,0,25,1)->Draw("alp");
  TStatToolkit::MakeStat1D((TH2*)his2D, 1,0.8,2,20,2)->Draw("lp");
  TStatToolkit::MakeStat1D((TH2*)his2D, 1,0.8,4,21,4)->Draw("lp");

}
