
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

  a) determine the required time range:

  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("pass2.txt","esdTree",0,6000);
  chain->Draw("GetTimeStamp()")

  b) analyse calibration object on Proof in calibration train 

  AliTPCcalibTime *calibTime = new AliTPCcalibTime("cosmicTime","cosmicTime", StartTimeStamp, EndTimeStamp, IntegrationTimeVdrift);

  c) plot results
  .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

  TFile f("CalibObjectsTrain1.root");
  AliTPCcalibTime *calib = (AliTPCcalibTime *)f->Get("calibTime");
  calib->GetHistoDrift("all")->Projection(2,0)->Draw()
  calib->GetFitDrift("all")->Draw("lp")

4. Analysis using debug streamers.    

  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chainTime = tool.MakeChain("time.txt","timeInfo",0,10200);
  chainTime->Lookup();
  TChain * chainLaser = tool.MakeChain("time.txt","timeLaser",0,10200);
  chainLaser->Lookup();
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
   fLaser(0),       // pointer to laser calibration
   fDz(0),          // current delta z
   fCutMaxD(3),        // maximal distance in rfi ditection
   fCutMaxDz(25),      // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99),  // direction vector products
   fCutTracks(10),
   fArrayDz(0),          //NEW! Tmap of V drifts for different triggers
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
  for (Int_t i=0;i<10;i++) {
    fCosmiMatchingHisto[i]=0;
  }
}

AliTPCcalibTime::AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeVdrift)
  :AliTPCcalibBase(),
   fLaser(0),            // pointer to laser calibration
   fDz(0),               // current delta z
   fCutMaxD(5*0.5356),   // maximal distance in rfi ditection
   fCutMaxDz(40),   // maximal distance in rfi ditection
   fCutTheta(5*0.004644),// maximal distan theta
   fCutMinDir(-0.99),    // direction vector products
   fCutTracks(10),
   fArrayDz(0),            //Tmap of V drifts for different triggers
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
  fTimeBins   =(EndTime-StartTime)/deltaIntegrationTimeVdrift;
  fTimeStart  =StartTime; //(((TObjString*)(mapGRP->GetValue("fAliceStartTime")))->GetString()).Atoi();
  fTimeEnd    =EndTime;   //(((TObjString*)(mapGRP->GetValue("fAliceStopTime")))->GetString()).Atoi();
  fPtBins     = 400;
  fPtStart    = -0.04;
  fPtEnd      =  0.04;
  fVdriftBins = 500;
  fVdriftStart= -0.1;
  fVdriftEnd  =  0.1;
  fRunBins    = 100001;
  fRunStart   = -1.5;
  fRunEnd     = 99999.5;

  Int_t    binsVdriftLaser[4] = {fTimeBins , fPtBins , fVdriftBins*20, fRunBins };
  Double_t xminVdriftLaser[4] = {fTimeStart, fPtStart, fVdriftStart  , fRunStart};
  Double_t xmaxVdriftLaser[4] = {fTimeEnd  , fPtEnd  , fVdriftEnd    , fRunEnd  };
  TString axisTitle[4]={
    "T",
    "#delta_{P/T}",
    "value",
    "run"
  };
  TString histoName[3]={
    "Loffset",
    "Lcorr",
    "Lgy"
  };

  
  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i] = new THnSparseF("HistVdriftLaser","HistVdriftLaser;time;p/T ratio;Vdrift;run",4,binsVdriftLaser,xminVdriftLaser,xmaxVdriftLaser);
    fHistVdriftLaserC[i] = new THnSparseF("HistVdriftLaser","HistVdriftLaser;time;p/T ratio;Vdrift;run",4,binsVdriftLaser,xminVdriftLaser,xmaxVdriftLaser);
    fHistVdriftLaserA[i]->SetName(histoName[i]);
    fHistVdriftLaserC[i]->SetName(histoName[i]);
    for (Int_t iaxis=0; iaxis<4;iaxis++){
      fHistVdriftLaserA[i]->GetAxis(iaxis)->SetName(axisTitle[iaxis]);
      fHistVdriftLaserC[i]->GetAxis(iaxis)->SetName(axisTitle[iaxis]);
    }
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

  fArrayDz=new TObjArray();
  fArrayDz->AddLast(fHistVdriftLaserA[0]);
  fArrayDz->AddLast(fHistVdriftLaserA[1]);
  fArrayDz->AddLast(fHistVdriftLaserA[2]);
  fArrayDz->AddLast(fHistVdriftLaserC[0]);
  fArrayDz->AddLast(fHistVdriftLaserC[1]);
  fArrayDz->AddLast(fHistVdriftLaserC[2]);

  fCosmiMatchingHisto[0]=new TH1F("Cosmics matching","p0-all"   ,100,-10*0.5356  ,10*0.5356  );
  fCosmiMatchingHisto[1]=new TH1F("Cosmics matching","p1-all"   ,100,-10*4.541   ,10*4.541   );
  fCosmiMatchingHisto[2]=new TH1F("Cosmics matching","p2-all"   ,100,-10*0.01134 ,10*0.01134 );
  fCosmiMatchingHisto[3]=new TH1F("Cosmics matching","p3-all"   ,100,-10*0.004644,10*0.004644);
  fCosmiMatchingHisto[4]=new TH1F("Cosmics matching","p4-all"   ,100,-10*0.03773 ,10*0.03773 );
  fCosmiMatchingHisto[5]=new TH1F("Cosmics matching","p0-isPair",100,-10*0.5356  ,10*0.5356  );
  fCosmiMatchingHisto[6]=new TH1F("Cosmics matching","p1-isPair",100,-10*4.541   ,10*4.541   );
  fCosmiMatchingHisto[7]=new TH1F("Cosmics matching","p2-isPair",100,-10*0.01134 ,10*0.01134 );
  fCosmiMatchingHisto[8]=new TH1F("Cosmics matching","p3-isPair",100,-10*0.004644,10*0.004644);
  fCosmiMatchingHisto[9]=new TH1F("Cosmics matching","p4-isPair",100,-10*0.03773 ,10*0.03773 );
//  Char_t nameHisto[3]={'p','0','\n'};
//  for (Int_t i=0;i<10;i++){
//    fCosmiMatchingHisto[i]=new TH1F("Cosmics matching",nameHisto,8192,0,0);
//    nameHisto[1]++;
//    if(i==4) nameHisto[1]='0';
//  }
}

AliTPCcalibTime::~AliTPCcalibTime(){
  //
  // Destructor
  //
  for(Int_t i=0;i<3;i++){
    if(fHistVdriftLaserA[i]){
      delete fHistVdriftLaserA[i];
      fHistVdriftLaserA[i]=NULL;
    }
    if(fHistVdriftLaserC[i]){
      delete fHistVdriftLaserC[i];
      fHistVdriftLaserC[i]=NULL;
    }
  }
  if(fArrayDz){
    fArrayDz->SetOwner();
    fArrayDz->Delete();
    delete fArrayDz;
    fArrayDz=NULL;
  }
  for(Int_t i=0;i<5;i++){
    if(fCosmiMatchingHisto[i]){
      delete fCosmiMatchingHisto[i];
      fCosmiMatchingHisto[i]=NULL;
    }
  }
}

Bool_t AliTPCcalibTime::IsLaser(AliESDEvent */*event*/){
  return kTRUE; //More accurate creteria to be added
}
Bool_t AliTPCcalibTime::IsCosmics(AliESDEvent */*event*/){
  return kTRUE; //More accurate creteria to be added
}
Bool_t AliTPCcalibTime::IsBeam(AliESDEvent */*event*/){
  return kTRUE; //More accurate creteria to be added
}
void AliTPCcalibTime::ResetCurrent(){
  fDz=0; //Reset current dz
}
void AliTPCcalibTime::Process(AliESDEvent *event){
  if(!event) return;
  if (event->GetNumberOfTracks()<2) return;
  ResetCurrent();
  if(IsLaser  (event)) ProcessLaser (event);
  if(IsCosmics(event)) ProcessCosmic(event);
  if(IsBeam   (event)) ProcessBeam  (event);
}

void AliTPCcalibTime::ProcessLaser(AliESDEvent *event){
  //
  // Fit drift velocity using laser 
  // 
  // 0. cuts
  const Int_t    kMinTracks     = 40;    // minimal number of laser tracks
  const Int_t    kMinTracksSide = 20;    // minimal number of tracks per side
  const Float_t  kMaxDeltaZ     = 30.;   // maximal trigger delay
  const Float_t  kMaxDeltaV     = 0.05;  // maximal deltaV 
  const Float_t  kMaxRMS        = 0.1;   // maximal RMS of tracks
  //
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
  if (fLaser->fFitAside->GetNrows()==0  && fLaser->fFitCside->GetNrows()==0) return;  // no fit neither a or C side
  //
  // debug streamer  - activate stream level
  // Use it for tuning of the cuts
  //
  // cuts to be applied
  //
  Int_t isReject[2]={0,0};
  //
  // not enough tracks 
  if (TMath::Abs((*fLaser->fFitAside)[3]) < kMinTracksSide) isReject[0]|=1; 
  if (TMath::Abs((*fLaser->fFitCside)[3]) < kMinTracksSide) isReject[1]|=1; 
  // unreasonable z offset
  if (TMath::Abs((*fLaser->fFitAside)[0])>kMaxDeltaZ)  isReject[0]|=2;
  if (TMath::Abs((*fLaser->fFitCside)[0])>kMaxDeltaZ)  isReject[1]|=2;
  // unreasonable drift velocity
  if (TMath::Abs((*fLaser->fFitAside)[1]-1)>kMaxDeltaV)  isReject[0]|=4;
  if (TMath::Abs((*fLaser->fFitCside)[1]-1)>kMaxDeltaV)  isReject[1]|=4;
  // big chi2
  if (TMath::Sqrt(TMath::Abs((*fLaser->fFitAside)[4]))>kMaxRMS ) isReject[0]|=8;
  if (TMath::Sqrt(TMath::Abs((*fLaser->fFitCside)[4]))>kMaxRMS ) isReject[1]|=8;




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
	"rejectA="<<isReject[0]<<
	"rejectC="<<isReject[1]<<
	"laserA.="<<fLaser->fFitAside<<
	"laserC.="<<fLaser->fFitCside<<
	"laserAC.="<<fLaser->fFitACside<<
	"trigger="<<event->GetFiredTriggerClasses()<<
	"\n";
    }
  }
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

  TTimeStamp tstamp(fTime);
  Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
  Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
  Double_t driftA=0, driftC=0;
  if (vdriftA[1]>1.-kMaxDeltaV) driftA = 1./vdriftA[1]-1.;
  if (vdriftC[1]>1.-kMaxDeltaV) driftC = 1./vdriftC[1]-1.;
  //
  Double_t vecDriftLaserA[4]={fTime,(ptrelative0+ptrelative1)/2.0,driftA,event->GetRunNumber()};
  Double_t vecDriftLaserC[4]={fTime,(ptrelative0+ptrelative1)/2.0,driftC,event->GetRunNumber()};
  //  Double_t vecDrift[4]      ={fTime,(ptrelative0+ptrelative1)/2.0,1./((*(fLaser->fFitACside))[1])-1,event->GetRunNumber()};

  for (Int_t icalib=0;icalib<3;icalib++){
    if (icalib==0){ //z0 shift
      vecDriftLaserA[2]=vdriftA[0]/250.;
      vecDriftLaserC[2]=vdriftC[0]/250.;
    }
    if (icalib==1){ //vdrel shift
      vecDriftLaserA[2]=driftA;
      vecDriftLaserC[2]=driftC;
    }
    if (icalib==2){ //gy shift - full gy - full drift
      vecDriftLaserA[2]=vdriftA[2]/250.;
      vecDriftLaserC[2]=vdriftC[2]/250.;
    }
    if (npointsA>kMinTracks) fHistVdriftLaserA[icalib]->Fill(vecDriftLaserA);
    if (npointsC>kMinTracks) fHistVdriftLaserC[icalib]->Fill(vecDriftLaserC);
  }

//   THnSparse* curHist=new THnSparseF("","HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//   TString shortName=curHist->ClassName();
//   shortName+="_MEAN_DRIFT_LASER_";
//   delete curHist;
//   curHist=NULL;
//   TString name="";

//   name=shortName;
//   name+=event->GetFiredTriggerClasses();
//   name.ToUpper();
//   curHist=(THnSparseF*)fArrayDz->FindObject(name);
//   if(!curHist){
//     curHist=new THnSparseF(name,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//     fArrayDz->AddLast(curHist);
//   }
//   curHist->Fill(vecDrift);
	  
//   name=shortName;
//   name+="ALL";
//   name.ToUpper();
//   curHist=(THnSparseF*)fArrayDz->FindObject(name);
//   if(!curHist){
//     curHist=new THnSparseF(name,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//     fArrayDz->AddLast(curHist);
//   }
//   curHist->Fill(vecDrift);
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
  
  //fd
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  Int_t ntracks=event->GetNumberOfTracks();
  if (ntracks==0) return;
  if (ntracks > fCutTracks) return;
  
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
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    if (seed) tpcSeeds.AddAt(seed,i);
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
      param0.GetDZ(0,0,0,bz,dvertex0);
      param1.GetDZ(0,0,0,bz,dvertex1);
      Double_t xyz0[3];
      Double_t xyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair = IsPair(&param0,&param1);
      Bool_t isCross = IsCross(track0, track1);
      Bool_t isSame = IsSame(track0, track1);

      THnSparse* hist=new THnSparseF("","HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
      TString shortName=hist->ClassName();
      shortName+="_MEAN_VDRIFT_COSMICS_";
      delete hist;
      hist=NULL;

      if(isSame || (isCross && isPair)){
	if (track0->GetTPCNcls() > 80) {
	  fDz = param0.GetZ() - param1.GetZ();
	  if(track0->GetOuterParam()->GetZ()<0) fDz=-fDz;
	  TTimeStamp tstamp(fTime);
	  Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
	  Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
	  Double_t vecDrift[4]={fTime,(ptrelative0+ptrelative1)/2.0,fDz/500.0,event->GetRunNumber()};
          THnSparse* curHist=NULL;
          TString name="";

          name=shortName;
	  name+=event->GetFiredTriggerClasses();
	  name.ToUpper();
	  curHist=(THnSparseF*)fArrayDz->FindObject(name);
	  if(!curHist){
	    curHist=new THnSparseF(name,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
	    fArrayDz->AddLast(curHist);
	  }
//	  curHist=(THnSparseF*)(fMapDz->GetValue(event->GetFiredTriggerClasses()));
//	  if(!curHist){
//	    curHist=new THnSparseF(event->GetFiredTriggerClasses(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//	    fMapDz->Add(new TObjString(event->GetFiredTriggerClasses()),curHist);
//	  }
	  curHist->Fill(vecDrift);
	  
          name=shortName;
	  name+="ALL";
	  name.ToUpper();
	  curHist=(THnSparseF*)fArrayDz->FindObject(name);
	  if(!curHist){
	    curHist=new THnSparseF(name,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
	    fArrayDz->AddLast(curHist);
	  }
//	  curHist=(THnSparseF*)(fMapDz->GetValue("all"));
//	  if(!curHist){
//	    curHist=new THnSparseF("all","HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//	    fMapDz->Add(new TObjString("all"),curHist);
//	  }
	  curHist->Fill(vecDrift);
	}
      }
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (fStreamLevel>0){
        if (cstream){
        (*cstream)<<"trackInfo"<<
	  "tr0.="<<track0<<
	  "tr1.="<<track1<<
	  "p0.="<<&param0<<
	  "p1.="<<&param1<<
	  "isPair="<<isPair<<
	  "isCross="<<isCross<<
	  "isSame="<<isSame<<
	  "fDz="<<fDz<<
	  "fRun="<<fRun<<
	  "fTime="<<fTime<<
	  "\n";
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
	"trigger="<<event->GetFiredTriggerClasses()<<
	"\n";
    }
  }
  printf("Trigger: %s\n",event->GetFiredTriggerClasses().Data());
}

void AliTPCcalibTime::ProcessBeam(AliESDEvent */*event*/){
}

void AliTPCcalibTime::Analyze(){}

THnSparse* AliTPCcalibTime::GetHistoDrift(const char* name){
  TIterator* iterator = fArrayDz->MakeIterator();
  iterator->Reset();
  TString newName=name;
  newName.ToUpper();
  THnSparse* newHist=new THnSparseF(newName,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
  THnSparse* addHist=NULL;
  while((addHist=(THnSparseF*)iterator->Next())){
  if(!addHist) continue;
    TString histName=addHist->GetName();
    if(!histName.Contains(newName)) continue;
    addHist->Print();
    newHist->Add(addHist);
  }
  return newHist;
}

TObjArray* AliTPCcalibTime::GetHistoDrift(){
  return fArrayDz;
}

TGraphErrors* AliTPCcalibTime::GetGraphDrift(const char* name){
  THnSparse* histoDrift=GetHistoDrift(name);
  TGraphErrors* graphDrift=NULL;
  if(histoDrift){
    graphDrift=FitSlices(histoDrift,2,0,400,100,0.05,0.95, kTRUE);
    TString end=histoDrift->GetName();
    Int_t pos=end.Index("_");
    end=end(pos,end.Capacity()-pos);
    TString graphName=graphDrift->ClassName();
    graphName+=end;
    graphName.ToUpper();
    graphDrift->SetName(graphName);
  }
  return graphDrift;
}

TObjArray* AliTPCcalibTime::GetGraphDrift(){
  TObjArray* arrayGraphDrift=new TObjArray();
  TIterator* iterator=fArrayDz->MakeIterator();
  iterator->Reset();
  THnSparse* addHist=NULL;
  while((addHist=(THnSparseF*)iterator->Next())) arrayGraphDrift->AddLast(GetGraphDrift(addHist->GetName()));
  return arrayGraphDrift;
}

AliSplineFit* AliTPCcalibTime::GetFitDrift(const char* name){
  TGraph* graphDrift=GetGraphDrift(name);
  AliSplineFit* fitDrift=NULL;
  if(graphDrift && graphDrift->GetN()){
    fitDrift=new AliSplineFit();
    fitDrift->SetGraph(graphDrift);
    fitDrift->SetMinPoints(graphDrift->GetN()+1);
    fitDrift->InitKnots(graphDrift,2,0,0.001);
    fitDrift->SplineFit(0);
    TString end=graphDrift->GetName();
    Int_t pos=end.Index("_");
    end=end(pos,end.Capacity()-pos);
    TString fitName=fitDrift->ClassName();
    fitName+=end;
    fitName.ToUpper();
    //fitDrift->SetName(fitName);
    delete graphDrift;
    graphDrift=NULL;
  }
  return fitDrift;
}

//TObjArray* AliTPCcalibTime::GetFitDrift(){
//  TObjArray* arrayFitDrift=new TObjArray();
//  TIterator* iterator = fArrayDz->MakeIterator();
//  iterator->Reset();
//  THnSparse* addHist=NULL;
//  while((addHist=(THnSparseF*)iterator->Next())) arrayFitDrift->AddLast(GetFitDrift(addHist->GetName()));
//  return arrayFitDrift;
//}

Long64_t AliTPCcalibTime::Merge(TCollection *li) {
  TIterator* iter = li->MakeIterator();
  AliTPCcalibTime* cal = 0;

  while ((cal = (AliTPCcalibTime*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibTime::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    for (Int_t imeas=0; imeas<3; imeas++){
      if (cal->GetHistVdriftLaserA(imeas) && cal->GetHistVdriftLaserA(imeas)){
	fHistVdriftLaserA[imeas]->Add(cal->GetHistVdriftLaserA(imeas));
	fHistVdriftLaserC[imeas]->Add(cal->GetHistVdriftLaserC(imeas));
      }
    }
    TObjArray* addArray=cal->GetHistoDrift();
    if(!addArray) return 0;
    TIterator* iterator = addArray->MakeIterator();
    iterator->Reset();
    THnSparse* addHist=NULL;
    while((addHist=(THnSparseF*)iterator->Next())){
      if(!addHist) continue;
      addHist->Print();
      THnSparse* localHist=(THnSparseF*)fArrayDz->FindObject(addHist->GetName());
      if(!localHist){
        localHist=new THnSparseF(addHist->GetName(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
        fArrayDz->AddLast(localHist);
      }
      localHist->Add(addHist);
    }
//    TMap * addMap=cal->GetHistoDrift();
//    if(!addMap) return 0;
//    TIterator* iterator = addMap->MakeIterator();
//    iterator->Reset();
//    TPair* addPair=0;
//    while((addPair=(TPair *)(addMap->FindObject(iterator->Next())))){
//      THnSparse* addHist=dynamic_cast<THnSparseF*>(addPair->Value());
//      if (!addHist) continue;
//      addHist->Print();
//      THnSparse* localHist=dynamic_cast<THnSparseF*>(fMapDz->GetValue(addHist->GetName()));
//      if(!localHist){
//        localHist=new THnSparseF(addHist->GetName(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
//        fMapDz->Add(new TObjString(addHist->GetName()),localHist);
//      }
//      localHist->Add(addHist);
//    }
    for(Int_t i=0;i<10;i++) if (cal->GetCosmiMatchingHisto(i)) fCosmiMatchingHisto[i]->Add(cal->GetCosmiMatchingHisto(i));
  }
  return 0;
}

Bool_t  AliTPCcalibTime::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1){
  /*
  // 0. Same direction - OPOSITE  - cutDir +cutT    
  TCut cutDir("cutDir","dir<-0.99")
  // 1. 
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03")
  //
  // 2. The same rphi 
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5")
  //
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");  
  // 1/Pt diff cut
  */
  const Double_t *p0 = tr0->GetParameter();
  const Double_t *p1 = tr1->GetParameter();
  fCosmiMatchingHisto[0]->Fill(p0[0]+p1[0]);
  fCosmiMatchingHisto[1]->Fill(p0[1]-p1[1]);
  fCosmiMatchingHisto[2]->Fill(tr1->GetAlpha()-tr0->GetAlpha()+TMath::Pi());
  fCosmiMatchingHisto[3]->Fill(p0[3]+p1[3]);
  fCosmiMatchingHisto[4]->Fill(p0[4]+p1[4]);
  
  if (TMath::Abs(p0[3]+p1[3])>fCutTheta) return kFALSE;
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  if (TMath::Abs(p0[1]-p1[1])>fCutMaxDz)  return kFALSE;
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;

  fCosmiMatchingHisto[5]->Fill(p0[0]+p1[0]);
  fCosmiMatchingHisto[6]->Fill(p0[1]-p1[1]);
  fCosmiMatchingHisto[7]->Fill(tr1->GetAlpha()-tr0->GetAlpha()+TMath::Pi());
  fCosmiMatchingHisto[8]->Fill(p0[3]+p1[3]);
  fCosmiMatchingHisto[9]->Fill(p0[4]+p1[4]);

  return kTRUE;  
}
Bool_t AliTPCcalibTime::IsCross(AliESDtrack *tr0, AliESDtrack *tr1){
  return  tr0->GetOuterParam()->GetZ()*tr1->GetOuterParam()->GetZ()<0 && tr0->GetInnerParam()->GetZ()*tr1->GetInnerParam()->GetZ()<0 && tr0->GetOuterParam()->GetZ()*tr0->GetInnerParam()->GetZ()>0 && tr1->GetOuterParam()->GetZ()*tr1->GetInnerParam()->GetZ()>0;
}

Bool_t AliTPCcalibTime::IsSame(AliESDtrack */*tr0*/, AliESDtrack */*tr1*/){
  // To be implemented
  return kFALSE;
}

/*
chainDrift->Draw("p0.fP[0]+p1.fP[0]","isPair");
  mean ~-0.02  ~-0.03913
  RMS  ~ 0.5   ~ 0.5356    --> 3    (fCutMaxD)

chainDrift->Draw("p0.fP[1]-p1.fP[1]","isPair");
  mean         ~ 0.1855
  RMS          ~ 4.541     -->25    (fCutMaxDz)

chainDrift->Draw("p1.fAlpha-p0.fAlpha+pi","isPair");
//chainDrift->Draw("p1.fAlpha+p0.fAlpha","isPair");
//chainDrift->Draw("p1.fP[2]-p0.fP[2]+pi","isPair");
//chainDrift->Draw("p1.fP[2]+p0.fP[2]","isPair");
  mean ~ 0     ~ 0.001898
  RMS  ~ 0.009 ~ 0.01134   --> 0.06

chainDrift->Draw("p0.fP[3]+p1.fP[3]","isPair");
  mean ~ 0.0013 ~ 0.001539
  RMS  ~ 0.003  ~ 0.004644 --> 0.03  (fCutTheta)

chainDrift->Draw("p1.fP[4]+p0.fP[4]>>his(100,-0.2,0.2)","isPair")
  mean ~ 0.012  ~-0.0009729
  RMS  ~ 0.036  ~ 0.03773  --> 0.2
*/

