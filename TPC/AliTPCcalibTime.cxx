/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
    Comments to be written here:
    1. What do we calibrate.

Time dependence of gain and drift velocity in order to account for changes in: temperature, pressure, gas composition.

  AliTPCcalibTime *calibTime = new AliTPCcalibTime("cosmicTime","cosmicTime",0, 1213.9e+06, 1213.96e+06, 0.04e+04, 0.04e+04);


    2. How to interpret results
    3. Simple example

a.) determine the required time range:

AliXRDPROOFtoolkit tool;
TChain * chain = tool.MakeChain("pass2.txt","esdTree",0,6000);
chain->Draw("GetTimeStamp()")

b) analyse calibration object on Proof in calibration train 

AliTPCcalibTime *calibTime = new AliTPCcalibTime("cosmicTime","cosmicTime", StartTimeStamp, EndTimeStamp, IntegrationTimeVdrift, IntegrationTimeDeDx);
s
c) plot results

TFile f("CalibObjects.root");
AliTPCcalibTime *cal = (AliTPCcalibTime *)f->Get("TPCCalib")->FindObject("cosmicTime");
cal->GetHistVdrift()->Projection(1,0)->Draw()

    4. Analysis using debug streamers.    

    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
    AliXRDPROOFtoolkit tool;
    TChain * chainTime = tool.MakeChain("time.txt","timeInfo",0,10200);
    chainTime->Lookup();


*/


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibTime.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"

ClassImp(AliTPCcalibTime)


AliTPCcalibTime::AliTPCcalibTime() 
  :AliTPCcalibBase(), 
   fTriggerMask(0),
   fHistDeDxTgl(0),
   fHistDeDx(0),
   fHistVdrift(0),
   fIntegrationTimeDeDx(0),
   fIntegrationTimeVdrift(0),  
   fLaser(0),       // pointer to laser calibration
   fDz(0),          // current delta z
   fdEdx(0),        // current dEdx
   fdEdxRatio(0),   // current dEdx ratio
   fTl(0),          // current tan(lambda)
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutMaxDz(20),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products

{  
  AliInfo("Default Constructor");  
}


AliTPCcalibTime::AliTPCcalibTime(const Text_t *name, const Text_t *title, ULong64_t TriggerMask, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeDeDx, Int_t deltaIntegrationTimeVdrift)
  :AliTPCcalibBase(),
   fTriggerMask(0),
   fHistDeDxTgl(0),
   fHistDeDx(0),
   fHistVdrift(0),
   fIntegrationTimeDeDx(0),
   fIntegrationTimeVdrift(0),  
   fLaser(0),       // pointer to laser calibration
   fDz(0),          // current delta z
   fdEdx(0),        // current dEdx
   fdEdxRatio(0),   // current dEdx ratio
   fTl(0),          // current tan(lambda)
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutMaxDz(20),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{
  
  SetName(name);
  SetTitle(title);

  AliInfo("Non Default Constructor");

  fTriggerMask = TriggerMask;

  fIntegrationTimeDeDx = deltaIntegrationTimeDeDx;
  fIntegrationTimeVdrift = deltaIntegrationTimeVdrift;

  Double_t deltaTime = EndTime - StartTime;
  
  Int_t binsVdrift[2] = {TMath::Nint(deltaTime/deltaIntegrationTimeVdrift), 100};
  Double_t xminVdrift[2] = {StartTime, -20};
  Double_t xmaxVdrift[2] = {EndTime, 20};
  fHistVdrift = new THnSparseF("HistVdrift","vDrift; time;#Delta z",2,binsVdrift,xminVdrift,xmaxVdrift);

  Int_t binsDeDxTgl[3] = {TMath::Nint(deltaTime/deltaIntegrationTimeDeDx),30,100};
  Double_t xminDeDxTgl[3] = {StartTime,-1,0.7};
  Double_t xmaxDeDxTgl[3] = {EndTime,1,1.3};
  fHistDeDxTgl = new THnSparseF("HistDeDxTgl","dEdx vs tgl;time;tgl;dEdxUp/dEdxLow",3,binsDeDxTgl,xminDeDxTgl,xmaxDeDxTgl); 

  Int_t binsDeDx[2] = {TMath::Nint(deltaTime/deltaIntegrationTimeDeDx),100};
  Double_t xminDeDx[2] = {StartTime,1};
  Double_t xmaxDeDx[2] = {EndTime,100};
  fHistDeDx = new THnSparseF("HistDeDx","dEdx l;time;dEdx",2,binsDeDx,xminDeDx,xmaxDeDx);

}



AliTPCcalibTime::~AliTPCcalibTime(){
  //
  //
  //
}
void AliTPCcalibTime::ResetCurrent(){
  //
  // reset current values
  //
  fDz=0;          // current delta z
  fdEdx=0;        // current dEdx
  fdEdxRatio=0;   // current dEdx ratio
  fTl=0;          // current tan(lambda)

}


void AliTPCcalibTime::Process(AliESDEvent *event) {
  //
  //
  //
  Int_t ntracks=event->GetNumberOfTracks();
  if (ntracks<2) return;
  ResetCurrent();
  //
  ProcessCosmic(event);
  if (fTrigger==16){
    if (!fLaser) fLaser =  new AliTPCcalibLaser("laserTPC","laserTPC",kFALSE);
    fLaser->Process(event);
  }
  //  
  // fill debug streamer
  if (fStreamLevel>0 && fDz!=0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      TTimeStamp tstamp(fTime);
      Float_t valuePressure0  = AliTPCcalibDB::GetPressure(tstamp,fRun,0);
      Float_t valuePressure1 = AliTPCcalibDB::GetPressure(tstamp,fRun,1);
      Double_t ptrelative0   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
      Double_t ptrelative1   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
      Double_t temp0         = AliTPCcalibDB::GetTemperature(tstamp,fRun,0);
      Double_t temp1         = AliTPCcalibDB::GetTemperature(tstamp,fRun,1);
      TVectorD vecGoofie(20);
      AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(fRun);
      if (goofieArray) 
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor) vecGoofie[isensor]=gsensor->GetValue(tstamp);
	}

      TVectorD vdriftA, vdriftC,vdriftAC;
      if (fLaser && fTrigger==16) {
	if (fLaser->fFitAside)  vdriftA=*(fLaser->fFitAside);
	if (fLaser->fFitCside)  vdriftC=*(fLaser->fFitCside);
	if (fLaser->fFitACside) vdriftAC=*(fLaser->fFitACside);
      }
      (*cstream)<<"timeInfo"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	// Environment values
	"press0="<<valuePressure0<<
	"press1="<<valuePressure1<<
	"pt0="<<ptrelative0<<
	"pt1="<<ptrelative1<<
	"temp0="<<temp0<<
	"temp1="<<temp1<<
	"vecGoofie.=<<"<<&vecGoofie<<
	//
	// accumulated values
	//
	"fDz="<<fDz<<          //! current delta z
	"fdEdx="<<fdEdx<<        //! current dEdx
	"fdEdxRatio="<<fdEdxRatio<<   //! current dEdx ratio
	"fTl="<<fTl<<          //! current tan(lambda)
	//
	//laser
	//
	"laserA.="<<&vdriftA<<
	"laserC.="<<&vdriftC<<
	"laserAC.="<<&vdriftAC<<
	"\n";
    }
  }
}



void AliTPCcalibTime::ProcessCosmic(AliESDEvent *event) {

  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  if (event->GetTimeStamp() == 0 ) {
    Printf("no time stamp!");
    return;
  }

  if (fTriggerMask != 0 && event->GetTriggerMask() != fTriggerMask) return;

  UInt_t time = event->GetTimeStamp();

  //
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  Int_t ntracks=event->GetNumberOfTracks();
  if (ntracks==0) return;
  if (ntracks > 10) return; // temporary debug to remove LASER events


  if (GetDebugLevel()>1) printf("Hallo world: Im here\n");
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  
  TObjArray  tpcSeeds(ntracks);
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) {
   AliESDtrack *track = event->GetTrack(i);

   const AliExternalTrackParam * trackIn = track->GetInnerParam();
   const AliExternalTrackParam * trackOut = track->GetOuterParam();
   if (!trackIn) continue;
   if (!trackOut) continue;
   
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) { 
     tpcSeeds.AddAt(seed,i);
     if (track->GetTPCNcls() > 50) {
       Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
       Double_t TPCsignal = seed->CookdEdxNorm(0.0,0.6,1,0,159,0x0,kTRUE,kTRUE);
       Double_t vecDeDx[2] = {time, TPCsignal};
       if (meanP > 12) {
	 fdEdx = TPCsignal;
	 fHistDeDx->Fill(vecDeDx);
       }
     }
   }
  }
  
  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i);     
    // track0 - choosen upper part
    if (!track0) continue;
    if (!track0->GetOuterParam()) continue;
    if (track0->GetOuterParam()->GetAlpha()<0) continue;
    Double_t d1[3];
    track0->GetDirection(d1);    
    for (Int_t j=0;j<ntracks;++j) {
      if (i==j) continue;
       AliESDtrack *track1 = event->GetTrack(j);   
       //track 1 lower part
       if (!track1) continue;
       if (!track1->GetOuterParam()) continue;
       if (track1->GetOuterParam()->GetAlpha()>0) continue;
       //
       Double_t d2[3];
       track1->GetDirection(d2);
      
       AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
       AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
       if (! seed0) continue;
       if (! seed1) continue;
       Float_t dir = (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]);
       Float_t dist0  = track0->GetLinearD(0,0);
       Float_t dist1  = track1->GetLinearD(0,0);
       //
       // conservative cuts - convergence to be guarantied
       // applying before track propagation
       if (TMath::Abs(dist0+dist1)>fCutMaxD) continue;   // distance to the 0,0
       if (dir>fCutMinDir) continue;               // direction vector product
       Float_t bz = AliTracker::GetBz();
       Float_t dvertex0[2];   //distance to 0,0
       Float_t dvertex1[2];   //distance to 0,0 
       track0->GetDZ(0,0,0,bz,dvertex0);
       track1->GetDZ(0,0,0,bz,dvertex1);
       if (TMath::Abs(dvertex0[1])>250) continue;
       if (TMath::Abs(dvertex1[1])>250) continue;
       //
       //
       //
       Float_t dmax = TMath::Max(TMath::Abs(dist0),TMath::Abs(dist1));
       AliExternalTrackParam param0(*track0);
       AliExternalTrackParam param1(*track1);
       //
       // Propagate using Magnetic field and correct fo material budget
       //
       AliTracker::PropagateTrackTo(&param0,dmax+1,0.0005,3,kTRUE);
       AliTracker::PropagateTrackTo(&param1,dmax+1,0.0005,3,kTRUE);
       //
       // Propagate rest to the 0,0 DCA - z should be ignored
       //
       //Bool_t b0 = ;
       param0.PropagateToDCA(&vtx,bz,1000);
       //Bool_t b1 = 
       param1.PropagateToDCA(&vtx,bz,1000);
       //      
       param0.GetDZ(0,0,0,bz,dvertex0);
       param1.GetDZ(0,0,0,bz,dvertex1);
       //
       Double_t xyz0[3];//,pxyz0[3];
       Double_t xyz1[3];//,pxyz1[3];
       param0.GetXYZ(xyz0);
       param1.GetXYZ(xyz1);
       Bool_t isPair = IsPair(&param0,&param1);

       Double_t z0 = track0->GetOuterParam()->GetZ();
       Double_t z1 = track1->GetOuterParam()->GetZ();

       Double_t z0inner = track0->GetInnerParam()->GetZ();
       Double_t z1inner = track1->GetInnerParam()->GetZ();

       if (isPair && z0>0 && z0inner>0 && z1<0 && z1inner<0) {
	 Double_t vecVdrift[2] = {time, param0.GetZ() - param1.GetZ()};
	 
	 if (track0->GetTPCNcls() > 80) {
	   fHistVdrift->Fill(vecVdrift);
	   fDz = param0.GetZ() - param1.GetZ();
	 }
       }
       if (isPair && z0 > 0 && z1 > 0) {
	 if (track1->GetTPCNcls()> 110 && track0->GetTPCNcls()> 110 && seed1->CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE) > 0) {
	   Float_t dedxratio = seed0->CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE)/seed1->CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE);
	   Double_t vecDeDxTgl[3] = {time, track0->GetTgl(), seed0->CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE)/seed1->CookdEdxNorm(0,0.6,1,0,159,0,kFALSE,kTRUE)};
	   fHistDeDxTgl->Fill(vecDeDxTgl);
	   fdEdxRatio = dedxratio;
	   fTl = track0->GetTgl() ;
	 }
       }
       
    } // end 2nd order loop        
  } // end 1st order loop

}


void AliTPCcalibTime::Analyze() {
  //
  //
  //
  TH2D * hVdrift = GetHistVdrift()->Projection(1,0);
  hVdrift->Draw();
}


Long64_t AliTPCcalibTime::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibTime* cal = 0;

  while ((cal = (AliTPCcalibTime*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibTime::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }

    // add histograms here...
    fHistDeDxTgl->Add(cal->GetHistDeDxVsTgl());
    fHistVdrift->Add(cal->GetHistVdrift());
    fHistDeDx->Add(cal->GetHistDeDx());

  }
  
  return 0;
  
}



Bool_t  AliTPCcalibTime::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1){
  //
  //
  /*
  // 0. Same direction - OPOSITE  - cutDir +cutT    
  TCut cutDir("cutDir","dir<-0.99")
  // 1. 
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03")
  //
  // 2. The same rphi 
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5")
  //
  //
  //
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");  
  // 1/Pt diff cut
  */
  const Double_t *p0 = tr0->GetParameter();
  const Double_t *p1 = tr1->GetParameter();
  if (TMath::Abs(p0[3]+p1[3])>fCutTheta) return kFALSE;
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  if (TMath::Abs(p0[1]-p1[1])>fCutMaxDz)  return kFALSE;
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;
  //
  return kTRUE;  
}



/*
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
    AliXRDPROOFtoolkit tool;
    TChain * chainTime = tool.MakeChain("time.txt","timeInfo",0,10200);
    chainTime->Lookup();

TCut dzc("abs(fDz-500*(1-pt1))<4")
chainTime->SetMarkerSize(0.2);
chainTime->SetMarkerStyle(24);


chainTime->SetMarkerColor(2);
chainTime->Draw("fDz:time","trigger==4"+dzc);
chainTime->SetMarkerColor(1);
chainTime->Draw("fDz:time","trigger==8"+dzc,"same");
htemp->SetXTitle("time")
htemp->SetYTitle("#Delta_{z}(cm)")
gPad->SaveAs("~/Calibration/driftV/pic/cosmicdzraw_time.gif");


chainTime->SetMarkerColor(2);
chainTime->Draw("fDz:500*(1-pt1)","trigger==4"+dzc);
chainTime->SetMarkerColor(1);
chainTime->Draw("fDz:500*(1-pt1)","trigger==8"+dzc,"same");
htemp->SetXTitle("2L#frac{#Delta_{P/T}}{P/T}")
htemp->SetYTitle("#Delta_{z}(cm)")
gPad->SaveAs("~/Calibration/driftV/pic/cosmicdzraw_PT1.gif");




chainTime->SetLineColor(2);
chainTime->Draw("fDz-500*0.965*(1-pt1)","trigger==4"+dzc);
chainTime->SetLineColor(1);
chainTime->Draw("fDz-500*0.965*(1-pt1)","trigger==8"+dzc,"same");
htemp->SetXTitle("#Delta_{z}(cm)")
htemp->SetYTitle("")
gPad->SaveAs("~/Calibration/driftV/pic/cosmicdzcorr_1D.gif");


chainTime->SetMarkerColor(2);
chainTime->Draw("fDz-500*0.965*(1-pt1):time","trigger==4"+dzc);
chainTime->SetMarkerColor(1);
chainTime->Draw("fDz-500*0.965*(1-pt1):time","trigger==8"+dzc,"same");
htemp->SetXTitle("time")
htemp->SetYTitle("#Delta_{z}(cm)")
gPad->SaveAs("~/Calibration/driftV/pic/cosmicdzcorr_time.gif");



gSystem->Load("libSTAT.so");
TStatToolkit toolkit;
Double_t chi2=0;
Int_t    npoints=0;
TVectorD fitParam;
TMatrixD covMatrix;


chainTime->SetAlias("dpr","(1-press0/970)");
chainTime->SetAlias("dtr","(1-(temp0+273.15)/293.)");
chainTime->SetAlias("d1pt","(1-(temp0+273.15)/293.)");
chainTime->SetAlias("dptr","(press0/(273.15+temp0))/(970./293.15)");


chainTime->SetAlias("dvr","fDz");
TString *strvd = toolkit.FitPlane(chainTime,"fDz","dpr*500++dtr*500", "trigger==4"+dzc, chi2,npoints,fitParam,covMatrix,0.99);
strvd.Tokenize("++")->Print();
chainTime->SetAlias("vdcorr",strvd->Data());

TString *strvdpt1 = toolkit.FitPlane(chainTime,"fDz","(1-dptr)*500", "trigger==4"+dzc, chi2,npoints,fitParam,covMatrix,0.99);
strvdpt1.Tokenize("++")->Print();
chainTime->SetAlias("vdcorrpt1",strvdpt1->Data());


TString *strdedx = toolkit.FitPlane(chainTime,"1-fdEdx/27","(1-dptr)", "trigger==4&&run<62000"+dzc, chi2,npoints,fitParam,covMatrix,0.99);
strdedx.Tokenize("++")->Print();
chainTime->SetAlias("vdcorrpt1",strvdpt1->Data());

*/
