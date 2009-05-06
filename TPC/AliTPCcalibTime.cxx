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
 .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

TFile f("CalibObjects.root");
AliTPCcalibTime *cal = (AliTPCcalibTime *)f->Get("TPCCalib")->FindObject("calibTime");
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
   fCutMaxDz(20),      // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99),  // direction vector products
   fMapDz(0),          //NEW! Tmap of V drifts for different triggers
   fTimeBins(0),
   fTimeStart(0),
   fTimeEnd(0),
   fPtBins(0),
   fPtStart(0),
   fPtEnd(0),
   fVdriftBins(0),
   fVdriftStart(0),
   fVdriftEnd(0),
   fRunBins(0),
   fRunStart(0),
   fRunEnd(0)
//   fBinsVdrift(fTimeBins,fPtBins,fVdriftBins),
//   fXminVdrift(fTimeStart,fPtStart,fVdriftStart),
//   fXmaxVdrift(fTimeEnd,fPtEnd,fVdriftEnd)

{  
  AliInfo("Default Constructor");  
  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i]=0;
    fHistVdriftLaserC[i]=0;
  }
}


AliTPCcalibTime::AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeDeDx, Int_t deltaIntegrationTimeVdrift)
  :AliTPCcalibBase(),
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
   fCutMaxDz(20),      // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99),  // direction vector products
   fMapDz(0),          //NEW! Tmap of V drifts for different triggers
   fTimeBins(0),
   fTimeStart(0),
   fTimeEnd(0),
   fPtBins(0),
   fPtStart(0),
   fPtEnd(0),
   fVdriftBins(0),
   fVdriftStart(0),
   fVdriftEnd(0),
   fRunBins(0),
   fRunStart(0),
   fRunEnd(0)
{
  
  SetName(name);
  SetTitle(title);
  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i]=0;
    fHistVdriftLaserC[i]=0;
  }

  AliInfo("Non Default Constructor");


  fIntegrationTimeDeDx = deltaIntegrationTimeDeDx;
  fIntegrationTimeVdrift = deltaIntegrationTimeVdrift;

  Double_t deltaTime = EndTime - StartTime;
  
  Int_t binsVdrift[2] = {TMath::Nint(deltaTime/deltaIntegrationTimeVdrift), 100};
  Double_t xminVdrift[2] = {StartTime, -20};
  Double_t xmaxVdrift[2] = {EndTime, 20};
  fHistVdrift = new THnSparseF("HistVdrift","vDrift; time;#Delta z",2,binsVdrift,xminVdrift,xmaxVdrift);

  Int_t binsDeDxTgl[4] = {TMath::Nint(deltaTime/deltaIntegrationTimeDeDx),30,100,10000};
  Double_t xminDeDxTgl[4] = {StartTime,-1,0.7,-0.5};
  Double_t xmaxDeDxTgl[4] = {EndTime,1,1.3,999.5};
  fHistDeDxTgl = new THnSparseF("HistDeDxTgl","dEdx vs tgl;time;tgl;dEdxUp/dEdxLow",3,binsDeDxTgl,xminDeDxTgl,xmaxDeDxTgl); 

  Int_t binsDeDx[2] = {TMath::Nint(deltaTime/deltaIntegrationTimeDeDx),100};
  Double_t xminDeDx[2] = {StartTime,1};
  Double_t xmaxDeDx[2] = {EndTime,100};
  fHistDeDx = new THnSparseF("HistDeDx","dEdx l;time;dEdx",2,binsDeDx,xminDeDx,xmaxDeDx);

//  fLaser =  new AliTPCcalibLaser("laserTPC","laserTPC",kFALSE);

  fTimeBins   =TMath::Nint(deltaTime/deltaIntegrationTimeVdrift);
  fTimeStart  =StartTime; //(((TObjString*)(mapGRP->GetValue("fAliceStartTime")))->GetString()).Atoi();
  fTimeEnd    =EndTime;   //(((TObjString*)(mapGRP->GetValue("fAliceStopTime")))->GetString()).Atoi();
  fPtBins     = 200;
  fPtStart    = -0.04;
  fPtEnd      =  0.04;
  fVdriftBins = 200;
  fVdriftStart= -20.0/500.0;
  fVdriftEnd  =  20.0/500.0;
  fRunBins    = 10000;
  fRunStart   = -0.5;
  fRunEnd     = 0.5;

  Int_t    binsVdriftLaser[4] = {fTimeBins , fPtBins , fVdriftBins*20, fRunBins };
  Double_t xminVdriftLaser[4] = {fTimeStart, fPtStart, fVdriftStart  , fRunStart};
  Double_t xmaxVdriftLaser[4] = {fTimeEnd  , fPtEnd  , fVdriftEnd    , fRunEnd  };

  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i] = new THnSparseF("HistVdriftLaser","HistVdriftLaser;time;p/T ratio;Vdrift;run",4,binsVdriftLaser,xminVdriftLaser,xmaxVdriftLaser);
    fHistVdriftLaserC[i] = new THnSparseF("HistVdriftLaser","HistVdriftLaser;time;p/T ratio;Vdrift;run",4,binsVdriftLaser,xminVdriftLaser,xmaxVdriftLaser);
  }
  fBinsVdrift[0] = fTimeBins;
  fBinsVdrift[1] = fPtBins;
  fBinsVdrift[2] = fVdriftBins;
  fBinsVdrift[3] = fRunBins;
  fXminVdrift[0] = fTimeStart;
  fXminVdrift[1] = fPtStart;
  fXminVdrift[2] = fVdriftStart;
  fXminVdrift[3] = fRunStart;
  fXmaxVdrift[0] = fTimeEnd;
  fXmaxVdrift[1] = fPtEnd;
  fXmaxVdrift[2] = fVdriftEnd;
  fXmaxVdrift[3] = fRunEnd;

  fMapDz=new TMap();
}



AliTPCcalibTime::~AliTPCcalibTime(){
  //
  // Destructor
  //
  delete fHistDeDxTgl;
  delete fHistDeDx;
  delete fHistVdrift;
  for (Int_t i=0;i<3;i++){
    delete fHistVdriftLaserA[i];
    delete fHistVdriftLaserC[i];
  }
  fMapDz->SetOwner();
  fMapDz->Delete();
  delete fMapDz;  
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

Bool_t AliTPCcalibTime::IsLaser(AliESDEvent *event){
  return ((event->GetFiredTriggerClasses()).Contains("0LSR")==1);
}

void AliTPCcalibTime::Process(AliESDEvent *event){
  if(!event) return;
  if (event->GetNumberOfTracks()<2) return;
  ResetCurrent();
  
//  if(IsLaser(event))
  ProcessLaser (event);
//  else               
  ProcessCosmic(event);
}

void AliTPCcalibTime::ProcessLaser(AliESDEvent *event){
  //
  // Fit drift velocity using laser 
  // 
  // 0. cuts
  const Int_t    kMinTracks     = 20;    // minimal number of laser tracks
  const Int_t    kMinTracksSide = 10;    // minimal number of tracks per side
  const Float_t  kMaxRMS        = 0.1;   // maximal RMS of tracks
  const Float_t  kMaxDeltaZ     = 3.;    // maximal deltaZ A-C side
  const Float_t  kMaxDeltaV     = 0.01;  // maximal deltaV A-C side
  const Float_t  kMaxDeltaY     = 2.;    // maximal deltaY A-C side
  /*
    TCut cutRMS("sqrt(laserA.fElements[4])<0.1&&sqrt(laserC.fElements[4])<0.1");
    TCut cutZ("abs(laserA.fElements[0]-laserC.fElements[0])<3");
    TCut cutV("abs(laserA.fElements[1]-laserC.fElements[1])<0.01");
    TCut cutY("abs(laserA.fElements[2]-laserC.fElements[2])<2");
    TCut cutAll = cutRMS+cutZ+cutV+cutY;
  */
  if (event->GetNumberOfTracks()<kMinTracks) return;
  //
  if(!fLaser) fLaser = new AliTPCcalibLaser("laserTPC","laserTPC",kFALSE);
  fLaser->Process(event);
  if (fLaser->GetNtracks()<kMinTracks) return;   // small amount of tracks cut
  if (fLaser->fFitAside->GetNrows()==0) return;  // no fit A side
  if (fLaser->fFitCside->GetNrows()==0) return;  // no fit C side
  //
  // debug streamer  - activate stream level
  // Use it for tuning of the cuts
  //
  if (fStreamLevel>0){
    printf("Trigger: %s\n",event->GetFiredTriggerClasses().Data());

    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      TTimeStamp tstamp(fTime);
      Float_t valuePressure0 = AliTPCcalibDB::GetPressure(tstamp,fRun,0);
      Float_t valuePressure1 = AliTPCcalibDB::GetPressure(tstamp,fRun,1);
      Double_t ptrelative0   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
      Double_t ptrelative1   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
      Double_t temp0         = AliTPCcalibDB::GetTemperature(tstamp,fRun,0);
      Double_t temp1         = AliTPCcalibDB::GetTemperature(tstamp,fRun,1);
      TVectorD vecGoofie(20);
      AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(fRun);
      if (goofieArray){
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor) vecGoofie[isensor]=gsensor->GetValue(tstamp);
	}
      }
      (*cstream)<<"laserInfo"<<
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
	//laser
	"laserA.="<<fLaser->fFitAside<<
	"laserC.="<<fLaser->fFitCside<<
	"laserAC.="<<fLaser->fFitACside<<
	"trigger="<<event->GetFiredTriggerClasses()<<
	"\n";
    }
  }
  //
  // Apply custs
  //
  if ((*fLaser->fFitAside)[3] <kMinTracksSide) return; // enough tracks A side
  if ((*fLaser->fFitCside)[3]<kMinTracksSide) return;  // enough tracks C side
  //
  if (TMath::Abs((*fLaser->fFitAside)[0]-(*fLaser->fFitCside)[0])>kMaxDeltaZ) return;
  if (TMath::Abs((*fLaser->fFitAside)[2]-(*fLaser->fFitCside)[2])>kMaxDeltaY) return;
  if (TMath::Abs((*fLaser->fFitAside)[1]-(*fLaser->fFitCside)[1])>kMaxDeltaV) return;
  if (TMath::Sqrt(TMath::Abs((*fLaser->fFitAside)[4]))>kMaxRMS) return;
  if (TMath::Sqrt(TMath::Abs((*fLaser->fFitCside)[4]))>kMaxRMS) return;
  //
  // fill histos
  //
  TVectorD vdriftA(5), vdriftC(5),vdriftAC(5);
  vdriftA=*(fLaser->fFitAside);
  vdriftC=*(fLaser->fFitCside);
  vdriftAC=*(fLaser->fFitACside);
  Int_t npointsA=0, npointsC=0;
  Float_t chi2A=0, chi2C=0;
  npointsA= TMath::Nint(vdriftA[3]);
  chi2A= vdriftA[4];
  npointsC= TMath::Nint(vdriftC[3]);
  chi2C= vdriftC[4];
  //
  //

  //  
  TTimeStamp tstamp(fTime);
  Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
  Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
  //
  Double_t vecVdriftLaserA[4]={fTime,(ptrelative0+ptrelative1)/2.0,1./((*(fLaser->fFitAside))[1])-1,event->GetRunNumber()};
  Double_t vecVdriftLaserC[4]={fTime,(ptrelative0+ptrelative1)/2.0,1./((*(fLaser->fFitCside))[1])-1,event->GetRunNumber()};

  for (Int_t i=0;i<3;i++){
    if (i==0){ //z0 shift
      vecVdriftLaserA[3]=(*(fLaser->fFitAside))[0]/250.;
      vecVdriftLaserA[3]=(*(fLaser->fFitCside))[0]/250.;
    }
    if (i==1){ //vdrel shift
      vecVdriftLaserA[3]=1./(*(fLaser->fFitAside))[1]-1.;
      vecVdriftLaserA[3]=1./(*(fLaser->fFitCside))[1]-1.;
    }
    if (i==2){ //gy shift - full gy - full drift
      vecVdriftLaserA[3]=(*(fLaser->fFitAside))[2]/250.;
      vecVdriftLaserA[3]=(*(fLaser->fFitCside))[2]/250.;
    }
    fHistVdriftLaserA[i]->Fill(vecVdriftLaserA);
    fHistVdriftLaserC[i]->Fill(vecVdriftLaserC);
  }
  //
  //
  //


}

void AliTPCcalibTime::ProcessCosmic(AliESDEvent *event){

  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  if (event->GetTimeStamp() == 0 ) {
    Printf("no time stamp!");
    return;
  }


  UInt_t time = event->GetTimeStamp();

  //fd
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
	   fDz = param0.GetZ() - param1.GetZ();
           TTimeStamp tstamp(fTime);
           Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
           Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
           THnSparse* curHist   =(THnSparseF*)(fMapDz->GetValue(event->GetFiredTriggerClasses()));
           if(!curHist){
             curHist=new THnSparseF(event->GetFiredTriggerClasses(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
             fMapDz->Add(new TObjString(event->GetFiredTriggerClasses()),curHist);
           }
           Double_t vecVdrift2[4]={fTime,(ptrelative0+ptrelative1)/2.0,fDz/500.0,event->GetRunNumber()};
           curHist->Fill(vecVdrift2);

	   fHistVdrift->Fill(vecVdrift);
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

  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      TTimeStamp tstamp(fTime);
      Float_t valuePressure0 = AliTPCcalibDB::GetPressure(tstamp,fRun,0);
      Float_t valuePressure1 = AliTPCcalibDB::GetPressure(tstamp,fRun,1);
      Double_t ptrelative0   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
      Double_t ptrelative1   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
      Double_t temp0         = AliTPCcalibDB::GetTemperature(tstamp,fRun,0);
      Double_t temp1         = AliTPCcalibDB::GetTemperature(tstamp,fRun,1);
      TVectorD vecGoofie(20);
      AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(fRun);
      if (goofieArray){
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor) vecGoofie[isensor]=gsensor->GetValue(tstamp);
	}
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
	"trigger="<<event->GetFiredTriggerClasses()<<
	"\n";
    }
  }
  printf("Trigger: %s\n",event->GetFiredTriggerClasses().Data());
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

    for (Int_t imeas=0; imeas<3; imeas++){
      if (cal->GetHistVdriftLaserA(imeas) && cal->GetHistVdriftLaserA(imeas)){
	fHistVdriftLaserA[imeas]->Add(cal->GetHistVdriftLaserA(imeas));
	fHistVdriftLaserC[imeas]->Add(cal->GetHistVdriftLaserC(imeas));
      }
    }
    TMap * addMap=cal->GetMapDz();
    if(!addMap) return 0;
    TIterator* iterator = addMap->MakeIterator();
    iterator->Reset();
    TPair* addPair=0;
    while((addPair=(TPair *)(addMap->FindObject(iterator->Next())))){
      THnSparse* addHist=dynamic_cast<THnSparseF*>(addPair->Value());
      if (!addHist) continue;
      addHist->Print();
      THnSparse* localHist=dynamic_cast<THnSparseF*>(fMapDz->GetValue(addHist->GetName()));
      if(!localHist){
        localHist=new THnSparseF(addHist->GetName(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
        fMapDz->Add(new TObjString(addHist->GetName()),localHist);
      }
      localHist->Add(addHist);
    }
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


