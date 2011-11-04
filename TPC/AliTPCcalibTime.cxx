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

*/

#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TVectorD.h"
//#include "TChain.h"
//#include "TFile.h"

#include "AliDCSSensor.h"
#include "AliDCSSensorArray.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDfriend.h"
#include "AliLog.h"
#include "AliRelAlignerKalman.h"
#include "AliTPCCalROC.h"
#include "AliTPCParam.h"
#include "AliTPCTracklet.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliTPCcalibTime.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTrackPointArray.h"
#include "AliTracker.h"
#include "AliKFVertex.h"
#include <AliLog.h>

ClassImp(AliTPCcalibTime)


AliTPCcalibTime::AliTPCcalibTime() 
  :AliTPCcalibBase(),  
   fMemoryMode(1), // 0 -do not fill THnSparse with residuals  1- fill only important QA THn 2 - Fill all THnsparse for calibration
   fLaser(0),       // pointer to laser calibration
   fDz(0),          // current delta z
   fCutMaxD(3),        // maximal distance in rfi ditection
   fCutMaxDz(25),      // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99),  // direction vector products
   fCutTracks(2500),
   fArrayLaserA(0),      //laser  fit parameters C
   fArrayLaserC(0),      //laser  fit parameters A
   fArrayDz(0),          //NEW! Tmap of V drifts for different triggers
   fAlignITSTPC(0),      //alignemnt array ITS TPC match
   fAlignTRDTPC(0),      //alignemnt array TRD TPC match 
   fAlignTOFTPC(0),      //alignemnt array TOF TPC match
   fTimeKalmanBin(60*15), //time bin width for kalman - 15 minutes default
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
  //
  // default constructor
  //
  AliInfo("Default Constructor");  
  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i]=0;
    fHistVdriftLaserC[i]=0;
  }
  for (Int_t i=0;i<10;i++) {
    fCosmiMatchingHisto[i]=0;
  }
  //
  for (Int_t i=0;i<5;i++) {
    fResHistoTPCCE[i]=0;
    fResHistoTPCITS[i]=0;
    fResHistoTPCTRD[i]=0;
    fResHistoTPCTOF[i]=0;
    fResHistoTPCvertex[i]=0;
    fTPCVertex[i]=0;
  }
  for (Int_t i=0;i<12;i++) {
    fTPCVertex[i]=0;
  }
  for (Int_t i=0;i<5;i++) {
    fTPCVertexCorrelation[i]=0;
  }
  static Int_t counter=0;
  if (1) {
    TTimeStamp s;
    Int_t time=s;
    AliInfo(Form("Counter Constructor\t%d\t%d",counter,time));
    counter++;
  }

}

AliTPCcalibTime::AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeVdrift, Int_t memoryMode)
  :AliTPCcalibBase(),
   fMemoryMode(memoryMode), // 0 -do not fill THnSparse with residuals  1- fill only important QA THn 2 - Fill all THnsparse for calibration
   fLaser(0),            // pointer to laser calibration
   fDz(0),               // current delta z
   fCutMaxD(5*0.5356),   // maximal distance in rfi ditection
   fCutMaxDz(40),   // maximal distance in rfi ditection
   fCutTheta(5*0.004644),// maximal distan theta
   fCutMinDir(-0.99),    // direction vector products
   fCutTracks(2500),
   fArrayLaserA(new TObjArray(1000)),      //laser  fit parameters C
   fArrayLaserC(new TObjArray(1000)),      //laser  fit parameters A
   fArrayDz(0),            //Tmap of V drifts for different triggers
   fAlignITSTPC(0),      //alignemnt array ITS TPC match
   fAlignTRDTPC(0),      //alignemnt array TRD TPC match 
   fAlignTOFTPC(0),      //alignemnt array TOF TPC match
   fTimeKalmanBin(60*15), //time bin width for kalman - 15 minutes default
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
  //
  // Non deafaul constructor - to be used in the Calibration setups 
  //

  SetName(name);
  SetTitle(title);
  for (Int_t i=0;i<3;i++) {
    fHistVdriftLaserA[i]=0;
    fHistVdriftLaserC[i]=0;
  }

  for (Int_t i=0;i<5;i++) {
    fResHistoTPCCE[i]=0;
    fResHistoTPCITS[i]=0;
    fResHistoTPCTRD[i]=0;
    fResHistoTPCTOF[i]=0;
    fResHistoTPCvertex[i]=0;
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
  fRunBins    = 1000001;
  fRunStart   = -1.5;
  fRunEnd     = 999999.5;

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
  fAlignITSTPC = new TObjArray;      //alignemnt array ITS TPC match
  fAlignTRDTPC = new TObjArray;      //alignemnt array ITS TPC match
  fAlignTOFTPC = new TObjArray;      //alignemnt array ITS TPC match
  fAlignITSTPC->SetOwner(kTRUE);
  fAlignTRDTPC->SetOwner(kTRUE);
  fAlignTOFTPC->SetOwner(kTRUE);
  

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
  for (Int_t i=0;i<12;i++) {
    fTPCVertex[i]=0;
  }
  for (Int_t i=0;i<5;i++) {
    fTPCVertexCorrelation[i]=0;
  }
  BookDistortionMaps();
  
}

AliTPCcalibTime::~AliTPCcalibTime(){
  //
  // Virtual Destructor
  //
  static Int_t counter=0;
  if (1) {
    TTimeStamp s;
    Int_t time=s;
    AliInfo(Form("Counter Destructor\t%s\t%d\t%d",GetName(),counter,time));
    counter++;
  }
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

  for (Int_t i=0;i<5;i++) {
    delete fResHistoTPCCE[i];
    delete fResHistoTPCITS[i];
    delete fResHistoTPCTRD[i];
    delete fResHistoTPCTOF[i];
    delete fResHistoTPCvertex[i];
    fResHistoTPCCE[i]=0;
    fResHistoTPCITS[i]=0;
    fResHistoTPCTRD[i]=0;
    fResHistoTPCTOF[i]=0;
    fResHistoTPCvertex[i]=0;
  }

  if (fTPCVertex[0]) {
    for (Int_t i=0;i<12;i++)  delete fTPCVertex[i];
  }
  if (fTPCVertexCorrelation[0]) {
    for (Int_t i=0;i<5;i++)  delete fTPCVertexCorrelation[i];
  }

  if (fAlignITSTPC){
    fAlignITSTPC->SetOwner(kTRUE);
    fAlignTRDTPC->SetOwner(kTRUE);
    fAlignTOFTPC->SetOwner(kTRUE);
    
    fAlignITSTPC->Delete();
    fAlignTRDTPC->Delete();
    fAlignTOFTPC->Delete();
    delete fAlignITSTPC;
    delete fAlignTRDTPC;
    delete fAlignTOFTPC;
  }
}

// Bool_t AliTPCcalibTime::IsLaser(const AliESDEvent *const /*event*/) const{
//   //
//   // Indicator is laser event not yet implemented  - to be done using trigger info or event specie
//   //
//   return kTRUE; //More accurate creteria to be added
// }
// Bool_t AliTPCcalibTime::IsCosmics(const AliESDEvent *const /*event*/){
//   //
//   // Indicator is cosmic event not yet implemented - to be done using trigger info or event specie
//   //

//   return kTRUE; //More accurate creteria to be added
// }
// Bool_t AliTPCcalibTime::IsBeam(const AliESDEvent *const /*event*/) const{
//   //
//   // Indicator is physic event not yet implemented - to be done using trigger info or event specie
//   //

//   return kTRUE; //More accurate creteria to be added
// }
void AliTPCcalibTime::ResetCurrent(){
  //
  //ResetCurrent
  //
  fDz=0; //Reset current dz
}



void AliTPCcalibTime::Process(AliESDEvent *event){
  //
  // main function to make calibration
  //
  if(!event) return;
  if (event->GetNumberOfTracks()<2) return; 
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
    return;
  }
  if (ESDfriend->TestSkipBit()) return;
  
  ResetCurrent();
  //if(IsLaser  (event)) 
  ProcessLaser (event);
  //if(IsCosmics(event)) 
  ProcessCosmic(event);
  //if(IsBeam   (event)) 
  ProcessBeam  (event);
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
      (*cstream)<<"laserInfo"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
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
  TVectorD vdriftA(5), vdriftC(5),vdriftAC(6);
  vdriftA=*(fLaser->fFitAside);
  vdriftC=*(fLaser->fFitCside);
  vdriftAC=*(fLaser->fFitACside);
  Int_t npointsA=0, npointsC=0;
  Float_t chi2A=0, chi2C=0;
  npointsA= TMath::Nint(vdriftA[3]);
  chi2A= vdriftA[4];
  npointsC= TMath::Nint(vdriftC[3]);
  chi2C= vdriftC[4];

  if (npointsA>kMinTracksSide || npointsC>kMinTracksSide){
    TVectorD *fitA = new TVectorD(6);
    TVectorD *fitC = new TVectorD(6);
    for (Int_t ipar=0; ipar<5; ipar++){
      (*fitA)[ipar]=vdriftA[ipar];
      (*fitC)[ipar]=vdriftC[ipar];
    }
    (*fitA)[5]=fTime;
    (*fitC)[5]=fTime;
    fArrayLaserA->AddLast(fitA);
    fArrayLaserC->AddLast(fitC);
  }
  //
  
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
    //if (isReject[0]==0) fHistVdriftLaserA[icalib]->Fill(vecDriftLaserA);
    //if (isReject[1]==0) fHistVdriftLaserC[icalib]->Fill(vecDriftLaserC);
    fHistVdriftLaserA[icalib]->Fill(vecDriftLaserA);
    fHistVdriftLaserC[icalib]->Fill(vecDriftLaserC);
  }
}

void AliTPCcalibTime::ProcessCosmic(const AliESDEvent *const event){
  //
  // process Cosmic event - track matching A side C side
  //
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
  const Int_t kMinClustersCross =30;
  const Int_t kMinClusters      =80;
  Int_t ntracks=event->GetNumberOfTracks();
  if (ntracks==0) return;
  if (ntracks > fCutTracks) return;
  
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  AliESDfriend *esdFriend=(AliESDfriend*)(((AliESDEvent*)event)->FindListObject("AliESDfriend"));
  
  TObjArray  tpcSeeds(ntracks);
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  // track loop
  //
  TArrayI clusterSideA(ntracks);
  TArrayI clusterSideC(ntracks);
  for (Int_t i=0;i<ntracks;++i) {
    clusterSideA[i]=0;
    clusterSideC[i]=0;
    AliESDtrack *track = event->GetTrack(i);
    
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = track->GetOuterParam();
    if (!trackIn) continue;
    if (!trackOut) continue;
    
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
    if (friendTrack) ProcessSame(track,friendTrack,event);
    if (friendTrack) ProcessAlignITS(track,friendTrack,event,esdFriend);
    if (friendTrack) ProcessAlignTRD(track,friendTrack);
    if (friendTrack) ProcessAlignTOF(track,friendTrack);
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    if (seed) {
      tpcSeeds.AddAt(seed,i);
      Int_t nA=0, nC=0;
      for (Int_t irow=159;irow>0;irow--) {
	AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
	if (!cl) continue;
	if ((cl->GetDetector()%36)<18) nA++;
	if ((cl->GetDetector()%36)>=18) nC++;
      }
      clusterSideA[i]=nA;
      clusterSideC[i]=nC;
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
      if (track0->GetTPCNcls()+ track1->GetTPCNcls()< kMinClusters) continue;
      Int_t nAC = TMath::Max( TMath::Min(clusterSideA[i], clusterSideC[j]), 
			      TMath::Min(clusterSideC[i], clusterSideA[j]));
      if (nAC<kMinClustersCross) continue; 
      Int_t nA0=clusterSideA[i];
      Int_t nC0=clusterSideC[i];
      Int_t nA1=clusterSideA[j];
      Int_t nC1=clusterSideC[j];
      //      if (track1->GetOuterParam()->GetAlpha()>0) continue;
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
      if (TMath::Abs(TMath::Abs(dist0)-TMath::Abs(dist1))>fCutMaxD) continue;   // distance to the 0,0
      if (TMath::Abs(dir)<TMath::Abs(fCutMinDir)) continue;               // direction vector product
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
      AliTracker::PropagateTrackTo(&param0,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE);
      AliTracker::PropagateTrackTo(&param1,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE);
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

      if((isSame) || (isCross && isPair)){
	if (track0->GetTPCNcls()+ track1->GetTPCNcls()> 80) {
	  fDz = param0.GetZ() - param1.GetZ();
	  Double_t sign=(nA0>nA1)? 1:-1; 
	  fDz*=sign;
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
	  "nAC="<<nAC<<
	  "nA0="<<nA0<<
	  "nA1="<<nA1<<
	  "nC0="<<nC0<<
	  "nC1="<<nC1<<
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
      (*cstream)<<"timeInfo"<<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field
	// Environment values
	//
	// accumulated values
	//
	"fDz="<<fDz<<          //! current delta z
	"trigger="<<event->GetFiredTriggerClasses()<<
	"\n";
    }
  }
  if (GetDebugLevel()>20) printf("Trigger: %s\n",event->GetFiredTriggerClasses().Data());
}

void AliTPCcalibTime::ProcessBeam(const AliESDEvent *const event){
  //
  // Process beam data - calculates vartex
  //                     from A side and C side
  // Histogram the differences                 
  // 
  const Int_t kMinClusters  =80;
  const Int_t kMinTracks    =2;      // minimal number of tracks to define the vertex
  const Int_t kMinTracksVertex=30;   // minimal number of tracks to define the cumulative vertex
  const Double_t kMaxTgl    =1.2;    // maximal Tgl (z angle)
  const Double_t kMinPt     =0.2;    // minimal pt
  const Double_t kMaxD0     =5.;     // cut on distance to the primary vertex first guess
  const Double_t kMaxZ0     =20; 
  const Double_t kMaxD      =2.5;    // cut on distance to the primary vertex 
  const Double_t kMaxZ      =4;      // maximal z distance between tracks form the same side
  const Double_t kMaxChi2   =15;     // maximal chi2 of the TPCvertex 
  const Double_t kCumulCovarXY=0.003; //increase the error of cumul vertex 30 microns profile
  const Double_t kCumulCovarZ=250.;  //increase the error of cumul vertex
  const Double_t kMaxDvertex = 1.0;  // cut to accept the vertex; 
  //
  Int_t  flags=0;
  const  Int_t  kBuffSize=100;
  static Double_t deltaZ[kBuffSize]={0};
  static Int_t counterZ=0;
  static AliKFVertex cumulVertexA, cumulVertexC, cumulVertexAC; // cumulative vertex 
  AliKFVertex vertexA, vertexC;

  Float_t  dca0[2]={0,0};
  Double_t dcaVertex[2]={0,0};
  Int_t ntracks=event->GetNumberOfTracks();
  if (ntracks==0) return;
  if (ntracks > fCutTracks) return;
  //
  AliESDfriend *esdFriend=(AliESDfriend*)(((AliESDEvent*)event)->FindListObject("AliESDfriend"));
  //
  // Divide tracks to A and C side tracks - using the cluster indexes
  TObjArray tracksA(ntracks);  
  TObjArray tracksC(ntracks);  
  //
  AliESDVertex *vertexSPD =  (AliESDVertex *)event->GetPrimaryVertexSPD();
  AliESDVertex *vertex    =  (AliESDVertex *)event->GetPrimaryVertex();
  AliESDVertex *vertexTracks =  (AliESDVertex *)event->GetPrimaryVertexTracks();
  Double_t vertexZA[10000], vertexZC[10000];
  //
  Int_t ntracksA= 0;
  Int_t ntracksC= 0;
  //
  for (Int_t itrack=0;itrack<ntracks;itrack++) {
    AliESDtrack *track = event->GetTrack(itrack);
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(itrack);
    if (!friendTrack) continue;
    if (TMath::Abs(track->GetTgl())>kMaxTgl) continue;
    if (TMath::Abs(track->Pt())<kMinPt) continue;
    const AliExternalTrackParam * trackIn  = track->GetInnerParam();
    TObject *calibObject=0;
    AliTPCseed *seed = 0;
    Int_t nA=0, nC=0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    if (seed) {
      for (Int_t irow=159;irow>0;irow--) {
	AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
	if (!cl) continue;
	if ((cl->GetDetector()%36)<18) nA++;
	if ((cl->GetDetector()%36)>=18) nC++;
      }
      if ((nA>kMinClusters || nC>kMinClusters) && (nA*nC==0) ){
	track->GetImpactParameters(dca0[0],dca0[1]);
	if (TMath::Abs(dca0[0])>kMaxD0) continue;
	if (TMath::Abs(dca0[1])>kMaxZ0) continue;
	AliExternalTrackParam pTPCvertex(*trackIn);
	if (!AliTracker::PropagateTrackToBxByBz(&pTPCvertex,4.+4.*TMath::Abs(dca0[0]),0.1,2,kTRUE)) continue;
	pTPCvertex.PropagateToDCA(vertex,AliTracker::GetBz(), kMaxD, dcaVertex,0);
	if (TMath::Abs(dcaVertex[0])>kMaxD) continue;
	if (nA>kMinClusters &&nC==0) { tracksA.AddLast(pTPCvertex.Clone()); vertexZA[ntracksA++] = pTPCvertex.GetZ();}
	if (nC>kMinClusters &&nA==0) {tracksC.AddLast(pTPCvertex.Clone());  vertexZC[ntracksC++] = pTPCvertex.GetZ();}
      }
    }
  }
  Double_t medianZA=TMath::Median(ntracksA, vertexZA);  // tracks median
  Double_t medianZC=TMath::Median(ntracksC, vertexZC);  // tracks median
  //
  ntracksA= tracksA.GetEntriesFast();
  ntracksC= tracksC.GetEntriesFast();
  if (ntracksA>kMinTracks && ntracksC>kMinTracks){
    deltaZ[counterZ%kBuffSize]=medianZA-medianZC;
    counterZ+=1;
    Double_t medianDelta=(counterZ>=kBuffSize)? TMath::Median(kBuffSize, deltaZ): TMath::Median(counterZ, deltaZ);
    if (TMath::Abs(medianDelta-(medianZA-medianZC))>kMaxZ) flags+=16;
    // increse the error of cumulative vertex at the beginning of event
    cumulVertexA.Covariance(0,0)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexA.Covariance(1,1)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexA.Covariance(2,2)+=kCumulCovarZ*kCumulCovarZ;
    cumulVertexC.Covariance(0,0)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexC.Covariance(1,1)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexC.Covariance(2,2)+=kCumulCovarZ*kCumulCovarZ;
    cumulVertexAC.Covariance(0,0)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexAC.Covariance(1,1)+=kCumulCovarXY*kCumulCovarXY;
    cumulVertexAC.Covariance(2,2)+=kCumulCovarZ*kCumulCovarZ;
    //
    for (Int_t iA=0; iA<ntracksA; iA++){
      if (flags!=0) continue;
      AliExternalTrackParam *aliTrack =  (AliExternalTrackParam *)tracksA.At(iA);
      if (TMath::Abs(aliTrack->GetZ()-medianZA)>kMaxZ) continue;
      AliKFParticle part(*aliTrack,211);
      vertexA+=part;
    }	
    for (Int_t iC=0; iC<ntracksC; iC++){
      if (flags!=0) continue;
      AliExternalTrackParam *aliTrack =  (AliExternalTrackParam *)tracksC.At(iC);
      if (TMath::Abs(aliTrack->GetZ()-medianZC)>kMaxZ) continue;
      AliKFParticle part(*aliTrack,211);
      vertexC+=part;
    }	 
    //
    if (vertexA.GetNDF()<kMinTracks) flags+=32;
    if (vertexC.GetNDF()<kMinTracks) flags+=32;
    if (TMath::Abs(vertexA.Z()-medianZA)>kMaxZ) flags+=1;   //apply cuts
    if (TMath::Abs(vertexC.Z()-medianZC)>kMaxZ) flags+=2;
    if (TMath::Abs(vertexA.GetChi2()/vertexA.GetNDF()+vertexC.GetChi2()/vertexC.GetNDF())> kMaxChi2) flags+=4;
    //
    if (flags==0){
      for (Int_t iA=0; iA<ntracksA; iA++){
	if (flags!=0) continue;
	AliExternalTrackParam *aliTrack =  (AliExternalTrackParam *)tracksA.At(iA);
	if (TMath::Abs(aliTrack->GetZ()-medianZA)>kMaxZ) continue;
	AliKFParticle part(*aliTrack,211);
	cumulVertexA+=part;
	cumulVertexAC+=part;
      }	
      for (Int_t iC=0; iC<ntracksC; iC++){
	if (flags!=0) continue;
	AliExternalTrackParam *aliTrack =  (AliExternalTrackParam *)tracksC.At(iC);
	if (TMath::Abs(aliTrack->GetZ()-medianZC)>kMaxZ) continue;
	AliKFParticle part(*aliTrack,211);
	cumulVertexC+=part;
	cumulVertexAC+=part;
      }	     
      //
      if (TMath::Abs(cumulVertexA.X()-vertexA.X())>kMaxDvertex) flags+=64;
      if (TMath::Abs(cumulVertexA.Y()-vertexA.Y())>kMaxDvertex) flags+=64;
      if (TMath::Abs(cumulVertexA.Z()-vertexA.Z())>kMaxDvertex) flags+=64;
      //
      if (TMath::Abs(cumulVertexC.X()-vertexC.X())>kMaxDvertex) flags+=64;
      if (TMath::Abs(cumulVertexC.Y()-vertexC.Y())>kMaxDvertex) flags+=64;
      if (TMath::Abs(cumulVertexC.Z()-vertexC.Z())>kMaxDvertex) flags+=64;
      
      
      if ( flags==0 && cumulVertexC.GetNDF()>kMinTracksVertex&&cumulVertexA.GetNDF()>kMinTracksVertex){
	Double_t cont[2]={0,fTime};
	//
	cont[0]= cumulVertexA.X();
	fTPCVertex[0]->Fill(cont);
	cont[0]= cumulVertexC.X();
	fTPCVertex[1]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.X()-cumulVertexC.X());
	fTPCVertex[2]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.X()+cumulVertexC.X())-vertexSPD->GetX();
	fTPCVertex[3]->Fill(cont);
	//
	cont[0]= cumulVertexA.Y();
	fTPCVertex[4]->Fill(cont);
	cont[0]= cumulVertexC.Y();
	fTPCVertex[5]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.Y()-cumulVertexC.Y());
	fTPCVertex[6]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.Y()+cumulVertexC.Y())-vertexSPD->GetY();
	fTPCVertex[7]->Fill(cont);
	//
	//
	cont[0]= 0.5*(cumulVertexA.Z()+cumulVertexC.Z());
	fTPCVertex[8]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.Z()-cumulVertexC.Z());
	fTPCVertex[9]->Fill(cont);
	cont[0]= 0.5*(cumulVertexA.Z()-cumulVertexC.Z());
	fTPCVertex[10]->Fill(cont);      
	cont[0]= 0.5*(cumulVertexA.Z()+cumulVertexC.Z())-vertexSPD->GetZ();
	fTPCVertex[11]->Fill(cont);
	//
	Double_t correl[2]={0,0};
	//
	correl[0]=cumulVertexC.Z();
	correl[1]=cumulVertexA.Z();
	fTPCVertexCorrelation[0]->Fill(correl);   // fill A side :TPC
	correl[0]=cumulVertexA.Z();
	correl[1]=cumulVertexC.Z(); 
	fTPCVertexCorrelation[1]->Fill(correl);   // fill C side :TPC
	//
	correl[0]=vertexSPD->GetZ();
	correl[1]=cumulVertexA.Z()-correl[0];
	fTPCVertexCorrelation[2]->Fill(correl);   // fill A side :ITS
	correl[1]=cumulVertexC.Z()-correl[0]; 
	fTPCVertexCorrelation[3]->Fill(correl);   // fill C side :ITS
	correl[1]=0.5*(cumulVertexA.Z()+cumulVertexC.Z())-correl[0]; 
	fTPCVertexCorrelation[4]->Fill(correl);   // fill C side :ITS
      }
    }        
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      /*
	TCut cutChi2= "sqrt(vA.fChi2/vA.fNDF+vC.fChi2/vC.fNDF)<10";  // chi2 Cut e.g 10 	
	TCut cutXY= "sqrt((vA.fP[0]-vC.fP[0])^2+(vA.fP[0]-vC.fP[1])^2)<5";   // vertex Cut	
	TCut cutZ= "abs(vA.fP[2]-mZA)<3&&abs(vC.fP[2]-mZC)<5";           // vertex Cut	
	tree->Draw("sqrt(vA.fChi2/vA.fNDF)","sqrt(vA.fChi2/vA.fNDF)<100","")
	
      */
      //vertexA.Print();
      //vertexC.Print();      
      (*cstream)<<"vertexTPC"<<
	"flags="<<flags<<        // rejection flags
	"vSPD.="<<vertexSPD<<    // SPD vertex
	"vT.="<<vertexTracks<<   // track vertex
	"v.="<<vertex<<          // esd vertex
	"mZA="<<medianZA<<       // median Z position at vertex A side
	"mZC="<<medianZC<<       // median Z position at vertex C side
	"mDelta="<<medianDelta<< // median delta A side -C side
	"counter="<<counterZ<<    // counter Z
	//
	"vA.="<<&vertexA<<       // vertex A side
	"vC.="<<&vertexC<<       // vertex C side
	"cvA.="<<&cumulVertexA<<       // cumulative vertex A side
	"cvC.="<<&cumulVertexC<<       // cumulative vertex C side
	"cvAC.="<<&cumulVertexAC<<       // cumulative vertex A+C side
	"nA="<<ntracksA<<        // contributors
	"nC="<<ntracksC<<        // contributors
	"\n";
    }      
  }
  tracksA.Delete();
  tracksC.Delete();
}

void AliTPCcalibTime::Analyze(){
  //
  // Special macro to analyze result of calibration and extract calibration entries
  // Not yet ported to the Analyze function yet
  //
}

THnSparse* AliTPCcalibTime::GetHistoDrift(const char* name) const
{
  //
  // Get histogram for given trigger mask
  //
  TIterator* iterator = fArrayDz->MakeIterator();
  iterator->Reset();
  TString newName=name;
  newName.ToUpper();
  THnSparse* newHist=new THnSparseF(newName,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
  THnSparse* addHist=NULL;
  while((addHist=(THnSparseF*)iterator->Next())){
    //  if(!addHist) continue;
    TString histName=addHist->GetName();
    if(!histName.Contains(newName)) continue;
    addHist->Print();
    newHist->Add(addHist);
  }
  return newHist;
}

TObjArray* AliTPCcalibTime::GetHistoDrift() const
{
  //
  // return array of histograms
  //
  return fArrayDz;
}

TGraphErrors* AliTPCcalibTime::GetGraphDrift(const char* name){
  //
  // Make a drift velocity (delta Z) graph
  //
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
  //
  // make a array of drift graphs
  //
  TObjArray* arrayGraphDrift=new TObjArray();
  TIterator* iterator=fArrayDz->MakeIterator();
  iterator->Reset();
  THnSparse* addHist=NULL;
  while((addHist=(THnSparseF*)iterator->Next())) arrayGraphDrift->AddLast(GetGraphDrift(addHist->GetName()));
  return arrayGraphDrift;
}

AliSplineFit* AliTPCcalibTime::GetFitDrift(const char* name){
  //
  // Make a fit AliSplinefit  of drift velocity
  //
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


Long64_t AliTPCcalibTime::Merge(TCollection *const li) {
  //
  // Object specific merging procedure
  //
  TIterator* iter = li->MakeIterator();
  AliTPCcalibTime* cal = 0;
  //
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
    //
    if (fTPCVertexCorrelation[0] && cal->fTPCVertexCorrelation[0]){
      for (Int_t imeas=0; imeas<5; imeas++){
	if (fTPCVertexCorrelation[imeas] && cal->fTPCVertexCorrelation[imeas]) fTPCVertexCorrelation[imeas]->Add(cal->fTPCVertexCorrelation[imeas]);
      }
    }
      
    if (fTPCVertex[0] && cal->fTPCVertex[0]) 
      for (Int_t imeas=0; imeas<12; imeas++){
	if (fTPCVertex[imeas] && cal->fTPCVertex[imeas]) fTPCVertex[imeas]->Add(cal->fTPCVertex[imeas]);
      }
    
    if (fMemoryMode>0) for (Int_t imeas=0; imeas<5; imeas++){
      if (fMemoryMode>1){
	if ( cal->GetResHistoTPCCE(imeas) && cal->GetResHistoTPCCE(imeas)){
	  fResHistoTPCCE[imeas]->Add(cal->fResHistoTPCCE[imeas]);
	}else{
	  fResHistoTPCCE[imeas]=(THnSparse*)cal->fResHistoTPCCE[imeas]->Clone();
	}
      }
      //
      if ((fMemoryMode>0) &&cal->GetResHistoTPCITS(imeas) && cal->GetResHistoTPCITS(imeas)){
	if (fMemoryMode>1 || (imeas%2)==1) fResHistoTPCITS[imeas]->Add(cal->fResHistoTPCITS[imeas]);
	if (fMemoryMode>1) fResHistoTPCvertex[imeas]->Add(cal->fResHistoTPCvertex[imeas]);
      }
      //
      if ((fMemoryMode>1) && cal->fResHistoTPCTRD[imeas]){
	if (fResHistoTPCTRD[imeas])
	  fResHistoTPCTRD[imeas]->Add(cal->fResHistoTPCTRD[imeas]);
	else
	  fResHistoTPCTRD[imeas]=(THnSparse*)cal->fResHistoTPCTRD[imeas]->Clone();
      }
      //
      if  ((fMemoryMode>1) && cal->fResHistoTPCTOF[imeas]){
	if (fResHistoTPCTOF[imeas])
	  fResHistoTPCTOF[imeas]->Add(cal->fResHistoTPCTOF[imeas]);
	else
	  fResHistoTPCTOF[imeas]=(THnSparse*)cal->fResHistoTPCTOF[imeas]->Clone();      
      }
      //
      if (cal->fArrayLaserA){
	fArrayLaserA->Expand(fArrayLaserA->GetEntriesFast()+cal->fArrayLaserA->GetEntriesFast());
	fArrayLaserC->Expand(fArrayLaserC->GetEntriesFast()+cal->fArrayLaserC->GetEntriesFast());
	for (Int_t ical=0; ical<cal->fArrayLaserA->GetEntriesFast(); ical++){
	  if (cal->fArrayLaserA->UncheckedAt(ical)) fArrayLaserA->AddLast(cal->fArrayLaserA->UncheckedAt(ical)->Clone());
	  if (cal->fArrayLaserC->UncheckedAt(ical)) fArrayLaserC->AddLast(cal->fArrayLaserC->UncheckedAt(ical)->Clone());
	}
      }

    }
    TObjArray* addArray=cal->GetHistoDrift();
    if(!addArray) return 0;
    TIterator* iterator = addArray->MakeIterator();
    iterator->Reset();
    THnSparse* addHist=NULL;
    if ((fMemoryMode>1)) while((addHist=(THnSparseF*)iterator->Next())){
      //      if(!addHist) continue;
      addHist->Print();
      THnSparse* localHist=(THnSparseF*)fArrayDz->FindObject(addHist->GetName());
      if(!localHist){
        localHist=new THnSparseF(addHist->GetName(),"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
        fArrayDz->AddLast(localHist);
      }
      localHist->Add(addHist);
    }

    for(Int_t i=0;i<10;i++) if (cal->GetCosmiMatchingHisto(i)) fCosmiMatchingHisto[i]->Add(cal->GetCosmiMatchingHisto(i));
    //
    // Merge alignment
    //
    for (Int_t itype=0; itype<3; itype++){
      //
      //
      TObjArray *arr0= 0;
      TObjArray *arr1= 0;
      if (itype==0) {arr0=fAlignITSTPC; arr1=cal->fAlignITSTPC;}
      if (itype==1) {arr0=fAlignTRDTPC; arr1=cal->fAlignTRDTPC;}
      if (itype==2) {arr0=fAlignTOFTPC; arr1=cal->fAlignTOFTPC;}
      if (!arr1) continue;
      if (!arr0) arr0=new TObjArray(arr1->GetEntriesFast());
      if (arr1->GetEntriesFast()>arr0->GetEntriesFast()){
	arr0->Expand(arr1->GetEntriesFast());
      }
      for (Int_t i=0;i<arr1->GetEntriesFast(); i++){
	AliRelAlignerKalman *kalman1 = (AliRelAlignerKalman *)arr1->UncheckedAt(i);
	AliRelAlignerKalman *kalman0 = (AliRelAlignerKalman *)arr0->UncheckedAt(i);
	if (!kalman1)  continue;
	if (!kalman0) {arr0->AddAt(new AliRelAlignerKalman(*kalman1),i); continue;}
	kalman0->SetRejectOutliers(kFALSE);
	kalman0->Merge(kalman1);
      }
    }

  }
  return 0;
}

Bool_t  AliTPCcalibTime::IsPair(const AliExternalTrackParam *tr0, const AliExternalTrackParam *tr1){
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
Bool_t AliTPCcalibTime::IsCross(const AliESDtrack *const tr0, const AliESDtrack *const tr1){
  //
  // check if the cosmic pair of tracks crossed A/C side
  // 
  Bool_t result= tr0->GetOuterParam()->GetZ()*tr1->GetOuterParam()->GetZ()<0;
  if (result==kFALSE) return result;
  result=kTRUE;
  return result;
}

Bool_t AliTPCcalibTime::IsSame(const AliESDtrack *const tr0, const AliESDtrack *const tr1){
  // 
  // track crossing the CE
  // 0. minimal number of clusters 
  // 1. Same sector +-1
  // 2. Inner and outer track param on opposite side
  // 3. Outer and inner track parameter close each to other
  // 3. 
  Bool_t result=kTRUE;
  //
  // inner and outer on opposite sides in z
  //
  const Int_t knclCut0  = 30;
  const Double_t kalphaCut = 0.4;
  //
  // 0. minimal number of clusters
  //
  if (tr0->GetTPCNcls()<knclCut0) return kFALSE;
  if (tr1->GetTPCNcls()<knclCut0) return kFALSE;
  //
  // 1. alpha cut - sector+-1
  //
  if (TMath::Abs(tr0->GetOuterParam()->GetAlpha()-tr1->GetOuterParam()->GetAlpha())>kalphaCut) return kFALSE;
  //
  // 2. Z crossing
  //
  if (tr0->GetOuterParam()->GetZ()*tr0->GetInnerParam()->GetZ()>0) result&=kFALSE;
  if (tr1->GetOuterParam()->GetZ()*tr1->GetInnerParam()->GetZ()>0) result&=kFALSE;
  if (result==kFALSE){
    return result;
  }
  //
  //
  const Double_t *p0I = tr0->GetInnerParam()->GetParameter();
  const Double_t *p1I = tr1->GetInnerParam()->GetParameter();
  const Double_t *p0O = tr0->GetOuterParam()->GetParameter();
  const Double_t *p1O = tr1->GetOuterParam()->GetParameter();
  //
  if (TMath::Abs(p0I[0]-p1I[0])>fCutMaxD)  result&=kFALSE;
  if (TMath::Abs(p0I[1]-p1I[1])>fCutMaxDz) result&=kFALSE;
  if (TMath::Abs(p0I[2]-p1I[2])>fCutTheta) result&=kFALSE;
  if (TMath::Abs(p0I[3]-p1I[3])>fCutTheta) result&=kFALSE;
  if (TMath::Abs(p0O[0]-p1O[0])>fCutMaxD)  result&=kFALSE;
  if (TMath::Abs(p0O[1]-p1O[1])>fCutMaxDz) result&=kFALSE;
  if (TMath::Abs(p0O[2]-p1O[2])>fCutTheta) result&=kFALSE;
  if (TMath::Abs(p0O[3]-p1O[3])>fCutTheta) result&=kFALSE;
  if (result==kTRUE){
    result=kTRUE; // just to put break point here
  }
  return result;
}


void  AliTPCcalibTime::ProcessSame(const AliESDtrack *const track, AliESDfriendTrack *const friendTrack, const AliESDEvent *const event){
  //
  // Process  TPC tracks crossing CE
  //
  // 0. Select only track crossing the CE
  // 1. Cut on the track length
  // 2. Refit the the track on A and C side separatelly
  // 3. Fill time histograms
  const Int_t kMinNcl=100;
  const Int_t kMinNclS=25;  // minimul number of clusters on the sides
  const Double_t pimass=TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  const Double_t kMaxDy=1;  // maximal distance in y
  const Double_t kMaxDsnp=0.05;  // maximal distance in snp
  const Double_t kMaxDtheta=0.05;  // maximal distance in theta
  
  if (!friendTrack->GetTPCOut()) return;
  //
  // 0. Select only track crossing the CE
  //
  if (track->GetInnerParam()->GetZ()*friendTrack->GetTPCOut()->GetZ()>0) return;
  //
  // 1. cut on track length
  //
  if (track->GetTPCNcls()<kMinNcl) return;
  //
  // 2. Refit track sepparatel on A and C side
  //
  TObject *calibObject;
  AliTPCseed *seed = 0;
  for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
    if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
  }
  if (!seed) return;
  //
  AliExternalTrackParam trackIn(*track->GetInnerParam());
  AliExternalTrackParam trackOut(*track->GetOuterParam());
  Double_t cov[3]={0.01,0.,0.01}; //use the same errors
  Double_t xyz[3]={0,0.,0.0};  
  Double_t bz   =0;
  Int_t nclIn=0,nclOut=0;
  trackIn.ResetCovariance(1000.);
  trackOut.ResetCovariance(1000.);
  //
  //2.a Refit inner
  // 
  Int_t sideIn=0;
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    if (sideIn==0){
      if (cl->GetDetector()%36<18) sideIn=1;
      if (cl->GetDetector()%36>=18) sideIn=-1;
    }
    if (sideIn== -1 && (cl->GetDetector()%36)<18) break;
    if (sideIn==  1 &&(cl->GetDetector()%36)>=18) break;
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackIn.GetAlpha();
    if (TMath::Abs(dalpha)>0.01){
      if (!trackIn.Rotate(TMath::DegToRad()*(sector%18*20.+10.))) break;
    }
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    trackIn.GetXYZ(xyz);
    bz = AliTracker::GetBz(xyz);
    AliTracker::PropagateTrackToBxByBz(&trackIn,r[0],pimass,1.,kFALSE);
    if (!trackIn.PropagateTo(r[0],bz)) break;
    nclIn++;
    trackIn.Update(&r[1],cov);    
  }
  //
  //2.b Refit outer
  //
  Int_t sideOut=0;
  for (Int_t irow=159;irow>0;irow--) {
    AliTPCclusterMI *cl=seed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()<80) continue;
    if (sideOut==0){
      if (cl->GetDetector()%36<18) sideOut=1;
      if (cl->GetDetector()%36>=18) sideOut=-1;
      if (sideIn==sideOut) break;
    }
    if (sideOut== -1 && (cl->GetDetector()%36)<18) break;
    if (sideOut==  1 &&(cl->GetDetector()%36)>=18) break;
    //
    Int_t sector = cl->GetDetector();
    Float_t dalpha = TMath::DegToRad()*(sector%18*20.+10.)-trackOut.GetAlpha();
    if (TMath::Abs(dalpha)>0.01){
      if (!trackOut.Rotate(TMath::DegToRad()*(sector%18*20.+10.))) break;
    }
    Double_t r[3]={cl->GetX(),cl->GetY(),cl->GetZ()};
    trackOut.GetXYZ(xyz);
    bz = AliTracker::GetBz(xyz);
    AliTracker::PropagateTrackToBxByBz(&trackOut,r[0],pimass,1.,kFALSE);
    if (!trackOut.PropagateTo(r[0],bz)) break;
    nclOut++;
    trackOut.Update(&r[1],cov);    
  }
  trackOut.Rotate(trackIn.GetAlpha());
  Double_t meanX = (trackIn.GetX()+trackOut.GetX())*0.5;
  trackIn.PropagateTo(meanX,bz); 
  trackOut.PropagateTo(meanX,bz); 
  if (TMath::Abs(trackIn.GetY()-trackOut.GetY())>kMaxDy) return;
  if (TMath::Abs(trackIn.GetSnp()-trackOut.GetSnp())>kMaxDsnp) return;
  if (TMath::Abs(trackIn.GetTgl()-trackOut.GetTgl())>kMaxDtheta) return;
  if (TMath::Min(nclIn,nclOut)>kMinNclS){
    FillResHistoTPCCE(&trackIn,&trackOut);
  }
  TTreeSRedirector *cstream = GetDebugStreamer();
  if (cstream){
    TVectorD gxyz(3);
    trackIn.GetXYZ(gxyz.GetMatrixArray());
    TTimeStamp tstamp(fTime);
    (*cstream)<<"tpctpc"<<
      "run="<<fRun<<              //  run number
      "event="<<fEvent<<          //  event number
      "time="<<fTime<<            //  time stamp of event
      "trigger="<<fTrigger<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      //
      "sideIn="<<sideIn<<         // side at inner part
      "sideOut="<<sideOut<<         // side at puter part
      "xyz.="<<&gxyz<<             // global position
      "tIn.="<<&trackIn<<         // refitterd track in 
      "tOut.="<<&trackOut<<       // refitter track out
      "nclIn="<<nclIn<<           // 
      "nclOut="<<nclOut<<         //
      "\n";  
  }
  //
  // 3. Fill time histograms
  // Debug stremaer expression
  // chainTPCTPC->Draw("(tIn.fP[1]-tOut.fP[1])*sign(-tIn.fP[3]):tIn.fP[3]","min(nclIn,nclOut)>30","")
  if (TMath::Min(nclIn,nclOut)>kMinNclS){
    fDz = trackOut.GetZ()-trackIn.GetZ();
    if (trackOut.GetTgl()<0) fDz*=-1.;
    TTimeStamp tstamp(fTime);
    Double_t ptrelative0 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
    Double_t ptrelative1 = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
    Double_t vecDrift[4]={fTime,(ptrelative0+ptrelative1)/2.0,fDz/500.0,event->GetRunNumber()};
    //
    // fill histograms per trigger class and itegrated
    //
    THnSparse* curHist=NULL;
    for (Int_t itype=0; itype<2; itype++){
      TString name="MEAN_VDRIFT_CROSS_";  
      if (itype==0){
	name+=event->GetFiredTriggerClasses();
	name.ToUpper();
      }else{
	name+="ALL";
      }
      curHist=(THnSparseF*)fArrayDz->FindObject(name);
      if(!curHist){
	curHist=new THnSparseF(name,"HistVdrift;time;p/T ratio;Vdrift;run",4,fBinsVdrift,fXminVdrift,fXmaxVdrift);
	fArrayDz->AddLast(curHist);
      }
      curHist->Fill(vecDrift);
    }
  }

}

void  AliTPCcalibTime::ProcessAlignITS(AliESDtrack *const track, const AliESDfriendTrack *const friendTrack, const AliESDEvent *const event, AliESDfriend *const esdFriend){
  //
  // Process track - Update TPC-ITS alignment
  // Updates: 
  // 0. Apply standartd cuts 
  // 1. Recalucluate the current statistic median/RMS
  // 2. Apply median+-rms cut
  // 3. Update kalman filter
  //
  const Int_t    kMinTPC  = 80;    // minimal number of TPC cluster
  const Int_t    kMinITS  = 3;     // minimal number of ITS cluster
  const Double_t kMinZ    = 10;    // maximal dz distance
  const Double_t kMaxDy   = 2.;    // maximal dy distance
  const Double_t kMaxAngle= 0.07;  // maximal angular distance
  const Double_t kSigmaCut= 5;     // maximal sigma distance to median
  const Double_t kVdErr   = 0.1;  // initial uncertainty of the vd correction 
  const Double_t kT0Err   = 3.;  // initial uncertainty of the T0 time
  const Double_t kVdYErr  = 0.05;  // initial uncertainty of the vd correction 
  const Double_t kOutCut  = 1.0;   // outlyer cut in AliRelAlgnmentKalman
  const Double_t kMinPt   = 0.3;   // minimal pt
  const Double_t kMax1Pt=0.5;        //maximal 1/pt distance
  const  Int_t     kN=50;         // deepnes of history
  static Int_t     kglast=0;
  static Double_t* kgdP[4]={new Double_t[kN], new Double_t[kN], new Double_t[kN], new Double_t[kN]};
  //
  // 0. Apply standard cuts
  // 
  Int_t dummycl[1000];
  if (track->GetTPCNcls()<kMinTPC) return;  // minimal amount of clusters cut
  if (!track->IsOn(AliESDtrack::kTPCrefit)) return;
  if (!track->GetInnerParam())   return;
  if (!track->GetOuterParam())   return;
  if (track->GetInnerParam()->Pt()<kMinPt)  return;
  // exclude crossing track
  if (track->GetOuterParam()->GetZ()*track->GetInnerParam()->GetZ()<0)   return;
  if (TMath::Abs(track->GetInnerParam()->GetZ())<kMinZ/3.)   return;
  if (track->GetInnerParam()->GetX()>90)   return;
  //
  AliExternalTrackParam &pTPC=(AliExternalTrackParam &)(*(track->GetInnerParam()));
  //  
  AliExternalTrackParam pITS;   // ITS standalone if possible
  AliExternalTrackParam pITS2;  //TPC-ITS track
  if (friendTrack->GetITSOut()){
    pITS2=(*(friendTrack->GetITSOut()));  //TPC-ITS track - snapshot ITS out
    pITS2.Rotate(pTPC.GetAlpha());
    AliTracker::PropagateTrackToBxByBz(&pITS2,pTPC.GetX(),0.1,0.1,kFALSE);
  }

  AliESDfriendTrack *itsfriendTrack=0;
  //
  // try to find standalone ITS track corresponing to the TPC if possible
  //
  Bool_t hasAlone=kFALSE;
  Int_t ntracks=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntracks; i++){
    AliESDtrack * trackITS = event->GetTrack(i); 
    if (!trackITS) continue;
    if (trackITS->GetITSclusters(dummycl)<kMinITS) continue;  // minimal amount of clusters
    itsfriendTrack = esdFriend->GetTrack(i);
    if (!itsfriendTrack) continue;
    if (!itsfriendTrack->GetITSOut()) continue;
     
    if (TMath::Abs(pTPC.GetTgl()-itsfriendTrack->GetITSOut()->GetTgl())> kMaxAngle) continue;
    if (TMath::Abs(pTPC.GetSigned1Pt()-itsfriendTrack->GetITSOut()->GetSigned1Pt())> kMax1Pt) continue;
    pITS=(*(itsfriendTrack->GetITSOut()));
    //
    pITS.Rotate(pTPC.GetAlpha());
    AliTracker::PropagateTrackToBxByBz(&pITS,pTPC.GetX(),0.1,0.1,kFALSE);
    if (TMath::Abs(pTPC.GetY()-pITS.GetY())> kMaxDy) continue;
    if (TMath::Abs(pTPC.GetSnp()-pITS.GetSnp())> kMaxAngle) continue;
    hasAlone=kTRUE;
  }
  if (!hasAlone) {
    if (track->GetITSclusters(dummycl)<kMinITS) return;
    pITS=pITS2;  // use combined track if it has ITS
  }
  //
  if (TMath::Abs(pITS.GetY()-pTPC.GetY())    >kMaxDy)    return;
  if (TMath::Abs(pITS.GetSnp()-pTPC.GetSnp())>kMaxAngle) return;
  if (TMath::Abs(pITS.GetTgl()-pTPC.GetTgl())>kMaxAngle) return;
  //
  // 1. Update median and RMS info
  //
  TVectorD vecDelta(5),vecMedian(5), vecRMS(5);
  TVectorD vecDeltaN(5);
  Double_t sign=(pITS.GetParameter()[1]>0)? 1.:-1.;
  vecDelta[4]=0;
  for (Int_t i=0;i<4;i++){
    vecDelta[i]=(pITS.GetParameter()[i]-pTPC.GetParameter()[i])*sign;
    kgdP[i][kglast%kN]=vecDelta[i];
  }
  kglast=(kglast+1);
  Int_t entries=(kglast<kN)?kglast:kN;
  for (Int_t i=0;i<4;i++){
    vecMedian[i] = TMath::Median(entries,kgdP[i]);
    vecRMS[i]    = TMath::RMS(entries,kgdP[i]);
    vecDeltaN[i] = 0;
    if (vecRMS[i]>0.){
      vecDeltaN[i] = (vecDelta[i]-vecMedian[i])/vecRMS[i];
      vecDeltaN[4]+= TMath::Abs(vecDeltaN[i]);  //sum of abs residuals
    }
  }
  //
  // 2. Apply median+-rms cut
  //
  if (kglast<3)  return;   //median and RMS to be defined
  if ( vecDeltaN[4]/4.>kSigmaCut) return;
  //
  // 3. Update alignment
  //
  Int_t htime = (fTime-fTimeKalmanBin/2)/fTimeKalmanBin; //time bins number
  if (fAlignITSTPC->GetEntriesFast()<htime){
    fAlignITSTPC->Expand(htime*2+20);
  }
  AliRelAlignerKalman* align =  (AliRelAlignerKalman*)fAlignITSTPC->At(htime);
  if (!align){
    // make Alignment object if doesn't exist
    align=new AliRelAlignerKalman(); 
    align->SetRunNumber(fRun);
    (*align->GetStateCov())(6,6)=kVdErr*kVdErr;
    (*align->GetStateCov())(7,7)=kT0Err*kT0Err;
    (*align->GetStateCov())(8,8)=kVdYErr*kVdYErr;
    align->SetOutRejSigma(kOutCut+kOutCut*kN);
    align->SetRejectOutliers(kFALSE);

    align->SetTPCvd(AliTPCcalibDB::Instance()->GetParameters()->GetDriftV()/1000000.);
    align->SetMagField(fMagF); 
    fAlignITSTPC->AddAt(align,htime);
  }
  align->AddTrackParams(&pITS,&pTPC);
  Double_t averageTime =  fTime;
  if (align->GetTimeStamp()>0&&align->GetNUpdates()>0){
    averageTime=((Double_t(align->GetTimeStamp())*Double_t(align->GetNUpdates())+Double_t(fTime)))/(Double_t(align->GetNUpdates())+1.);
  }
  align->SetTimeStamp(Int_t(averageTime));

  align->SetRunNumber(fRun );
  Float_t dca[2],cov[3];
  track->GetImpactParameters(dca,cov);
  if (TMath::Abs(dca[0])<kMaxDy){
    FillResHistoTPCITS(&pTPC,&pITS);
    FillResHistoTPC(track);
  }
  //
  Int_t nupdates=align->GetNUpdates();
  align->SetOutRejSigma(kOutCut+kOutCut*kN/Double_t(nupdates));
  align->SetRejectOutliers(kFALSE);
  TTreeSRedirector *cstream = GetDebugStreamer();  
  if (cstream && align->GetState() && align->GetState()->GetNrows()>2 ){
    TVectorD gpTPC(3), gdTPC(3);
    TVectorD gpITS(3), gdITS(3);
    pTPC.GetXYZ(gpTPC.GetMatrixArray());
    pTPC.GetDirection(gdTPC.GetMatrixArray());
    pITS.GetXYZ(gpITS.GetMatrixArray());
    pITS.GetDirection(gdITS.GetMatrixArray());
    (*cstream)<<"itstpc"<<
      "run="<<fRun<<              //  run number
      "event="<<fEvent<<          //  event number
      "time="<<fTime<<            //  time stamp of event
      "trigger="<<fTrigger<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      //
      "hasAlone="<<hasAlone<<    // has ITS standalone ?
      "track.="<<track<<  // track info
      "nmed="<<kglast<<        // number of entries to define median and RMS
      "vMed.="<<&vecMedian<<    // median of deltas
      "vRMS.="<<&vecRMS<<       // rms of deltas
      "vDelta.="<<&vecDelta<<   // delta in respect to median
      "vDeltaN.="<<&vecDeltaN<< // normalized delta in respect to median
      "a.="<<align<<            // current alignment
      "pITS.="<<&pITS<<         // track param ITS 
      "pITS2.="<<&pITS2<<       // track param ITS+TPC
      "pTPC.="<<&pTPC<<         // track param TPC
      "gpTPC.="<<&gpTPC<<       // global position  TPC
      "gdTPC.="<<&gdTPC<<       // global direction TPC
      "gpITS.="<<&gpITS<<       // global position  ITS
      "gdITS.="<<&gdITS<<       // global position  ITS
      "\n";
  }
}




void  AliTPCcalibTime::ProcessAlignTRD(AliESDtrack *const track, const AliESDfriendTrack *const friendTrack){
  //
  // Process track - Update TPC-TRD alignment
  // Updates: 
  // 0. Apply standartd cuts 
  // 1. Recalucluate the current statistic median/RMS
  // 2. Apply median+-rms cut
  // 3. Update kalman filter
  //
  const Int_t    kMinTPC  = 80;    // minimal number of TPC cluster
  const Int_t    kMinTRD  = 50;    // minimal number of TRD cluster
  const Double_t kMinZ    = 20;    // maximal dz distance
  const Double_t kMaxDy   = 5.;    // maximal dy distance
  const Double_t kMaxAngle= 0.1;  // maximal angular distance
  const Double_t kSigmaCut= 10;     // maximal sigma distance to median
  const Double_t kVdErr   = 0.1;  // initial uncertainty of the vd correction 
  const Double_t kT0Err   = 3.;  // initial uncertainty of the T0 time
  const Double_t kVdYErr  = 0.05;  // initial uncertainty of the vd correction 
  const Double_t kOutCut  = 1.0;   // outlyer cut in AliRelAlgnmentKalman
  const Double_t kRefX    = 275;   // reference X
  const  Int_t     kN=50;         // deepnes of history
  static Int_t     kglast=0;
  static Double_t* kgdP[4]={new Double_t[kN], new Double_t[kN], new Double_t[kN], new Double_t[kN]};
  //
  // 0. Apply standard cuts
  //
  Int_t dummycl[1000];
  if (track->GetTRDclusters(dummycl)<kMinTRD) return;  // minimal amount of clusters
  if (track->GetTPCNcls()<kMinTPC) return;  // minimal amount of clusters cut
  if (!friendTrack->GetTRDIn()) return;  
  if (!track->IsOn(AliESDtrack::kTRDrefit)) return;  
  if (!track->IsOn(AliESDtrack::kTRDout)) return;  
  if (!track->GetInnerParam())   return;
  if (!friendTrack->GetTPCOut())   return;
  // exclude crossing track
  if (friendTrack->GetTPCOut()->GetZ()*track->GetInnerParam()->GetZ()<0)   return;
  if (TMath::Abs(track->GetInnerParam()->GetZ())<kMinZ)   return;
  //
  AliExternalTrackParam &pTPC=(AliExternalTrackParam &)(*(friendTrack->GetTPCOut()));
  AliTracker::PropagateTrackToBxByBz(&pTPC,kRefX,0.1,0.1,kFALSE);
  AliExternalTrackParam pTRD(*(friendTrack->GetTRDIn()));
  pTRD.Rotate(pTPC.GetAlpha());
  //  pTRD.PropagateTo(pTPC.GetX(),fMagF);
  AliTracker::PropagateTrackToBxByBz(&pTRD,pTPC.GetX(),0.1,0.1,kFALSE);

  ((Double_t*)pTRD.GetCovariance())[2]+=3.*3.;   // increas sys errors
  ((Double_t*)pTRD.GetCovariance())[9]+=0.1*0.1; // increse sys errors

  if (TMath::Abs(pTRD.GetY()-pTPC.GetY())    >kMaxDy)    return;
  if (TMath::Abs(pTRD.GetSnp()-pTPC.GetSnp())>kMaxAngle) return;
  //  if (TMath::Abs(pTRD.GetTgl()-pTPC.GetTgl())>kMaxAngle) return;
  //
  // 1. Update median and RMS info
  //
  TVectorD vecDelta(5),vecMedian(5), vecRMS(5);
  TVectorD vecDeltaN(5);
  Double_t sign=(pTRD.GetParameter()[1]>0)? 1.:-1.;
  vecDelta[4]=0;
  for (Int_t i=0;i<4;i++){
    vecDelta[i]=(pTRD.GetParameter()[i]-pTPC.GetParameter()[i])*sign;
    kgdP[i][kglast%kN]=vecDelta[i];
  }
  kglast=(kglast+1);
  Int_t entries=(kglast<kN)?kglast:kN;
  for (Int_t i=0;i<4;i++){
    vecMedian[i] = TMath::Median(entries,kgdP[i]);

    vecRMS[i]    = TMath::RMS(entries,kgdP[i]);
    vecDeltaN[i] = 0;
    if (vecRMS[i]>0.){
      vecDeltaN[i] = (vecDelta[i]-vecMedian[i])/vecRMS[i];
      vecDeltaN[4]+= TMath::Abs(vecDeltaN[i]);  //sum of abs residuals
    }
  }
  //
  // 2. Apply median+-rms cut
  //
  if (kglast<3)  return;   //median and RMS to be defined
  if ( vecDeltaN[4]/4.>kSigmaCut) return;
  //
  // 3. Update alignment
  //
  //Int_t htime = fTime/3600; //time in hours
  Int_t htime = (Int_t)(fTime-fTimeKalmanBin/2)/fTimeKalmanBin; //time in half hour
  if (fAlignTRDTPC->GetEntriesFast()<htime){
    fAlignTRDTPC->Expand(htime*2+20);
  }
  AliRelAlignerKalman* align =  (AliRelAlignerKalman*)fAlignTRDTPC->At(htime);
  if (!align){
    // make Alignment object if doesn't exist
    align=new AliRelAlignerKalman(); 
    align->SetRunNumber(fRun);
    (*align->GetStateCov())(6,6)=kVdErr*kVdErr;
    (*align->GetStateCov())(7,7)=kT0Err*kT0Err;
    (*align->GetStateCov())(8,8)=kVdYErr*kVdYErr;
    align->SetOutRejSigma(kOutCut+kOutCut*kN);
    align->SetRejectOutliers(kFALSE);
    align->SetTPCvd(AliTPCcalibDB::Instance()->GetParameters()->GetDriftV()/1000000.);
    align->SetMagField(fMagF); 
    fAlignTRDTPC->AddAt(align,htime);
  }
  align->AddTrackParams(&pTRD,&pTPC);
  //align->SetTimeStamp(fTime);
  Double_t averageTime =  fTime;
  if (align->GetTimeStamp()>0 && align->GetNUpdates()>0) {
    averageTime = (((Double_t)fTime) + ((Double_t)align->GetTimeStamp())*align->GetNUpdates()) / (align->GetNUpdates() + 1.);
    //printf("align->GetTimeStamp() %d, align->GetNUpdates() %d \n", align->GetTimeStamp(), align->GetNUpdates());
  }
  align->SetTimeStamp((Int_t)averageTime);

  //printf("fTime %d, averageTime %d \n", fTime, (Int_t)averageTime);

  align->SetRunNumber(fRun );
  Float_t dca[2],cov[3];
  track->GetImpactParameters(dca,cov);
  if (TMath::Abs(dca[0])<kMaxDy){
    FillResHistoTPCTRD(&pTPC,&pTRD);  //only primaries
  }
  //
  Int_t nupdates=align->GetNUpdates();
  align->SetOutRejSigma(kOutCut+kOutCut*kN/Double_t(nupdates));
  align->SetRejectOutliers(kFALSE);
  TTreeSRedirector *cstream = GetDebugStreamer();  
  if (cstream && align->GetState() && align->GetState()->GetNrows()>2 ){
    TVectorD gpTPC(3), gdTPC(3);
    TVectorD gpTRD(3), gdTRD(3);
    pTPC.GetXYZ(gpTPC.GetMatrixArray());
    pTPC.GetDirection(gdTPC.GetMatrixArray());
    pTRD.GetXYZ(gpTRD.GetMatrixArray());
    pTRD.GetDirection(gdTRD.GetMatrixArray());
    (*cstream)<<"trdtpc"<<
      "run="<<fRun<<              //  run number
      "event="<<fEvent<<          //  event number
      "time="<<fTime<<            //  time stamp of event
      "trigger="<<fTrigger<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      //
      "nmed="<<kglast<<        // number of entries to define median and RMS
      "vMed.="<<&vecMedian<<    // median of deltas
      "vRMS.="<<&vecRMS<<       // rms of deltas
      "vDelta.="<<&vecDelta<<   // delta in respect to median
      "vDeltaN.="<<&vecDeltaN<< // normalized delta in respect to median
      "t.="<<track<<            // ful track - find proper cuts
      "a.="<<align<<            // current alignment
      "pTRD.="<<&pTRD<<         // track param TRD
      "pTPC.="<<&pTPC<<         // track param TPC
      "gpTPC.="<<&gpTPC<<       // global position  TPC
      "gdTPC.="<<&gdTPC<<       // global direction TPC
      "gpTRD.="<<&gpTRD<<       // global position  TRD
      "gdTRD.="<<&gdTRD<<       // global position  TRD
      "\n";
  }
}


void  AliTPCcalibTime::ProcessAlignTOF(AliESDtrack *const track, const AliESDfriendTrack *const friendTrack){
  //
  //
  // Process track - Update TPC-TOF alignment
  // Updates: 
  // -1. Make a TOF "track"
  // 0. Apply standartd cuts 
  // 1. Recalucluate the current statistic median/RMS
  // 2. Apply median+-rms cut
  // 3. Update kalman filter
  //
  const Int_t      kMinTPC  = 80;    // minimal number of TPC cluster
  //  const Double_t   kMinZ    = 10;    // maximal dz distance
  const Double_t   kMaxDy   = 5.;    // maximal dy distance
  const Double_t   kMaxAngle= 0.05;  // maximal angular distance
  const Double_t   kSigmaCut= 5;     // maximal sigma distance to median
  const Double_t   kVdErr   = 0.1;  // initial uncertainty of the vd correction 
  const Double_t   kT0Err   = 3.;  // initial uncertainty of the T0 time
  const Double_t   kVdYErr  = 0.05;  // initial uncertainty of the vd correction 

  const Double_t   kOutCut  = 1.0;   // outlyer cut in AliRelAlgnmentKalman
  const  Int_t     kN=50;         // deepnes of history
  static Int_t     kglast=0;
  static Double_t* kgdP[4]={new Double_t[kN], new Double_t[kN], new Double_t[kN], new Double_t[kN]};
  //
  // -1. Make a TOF track-
  //     Clusters are not in friends - use alingment points
  //
  if (track->GetTOFsignal()<=0)  return;
  if (!friendTrack->GetTPCOut()) return;
  if (!track->GetInnerParam())   return;
  if (!friendTrack->GetTPCOut())   return;
  const AliTrackPointArray *points=friendTrack->GetTrackPointArray();
  if (!points) return;
  AliExternalTrackParam pTPC(*(friendTrack->GetTPCOut()));
  AliExternalTrackParam pTOF(pTPC);
  Double_t mass = TDatabasePDG::Instance()->GetParticle("mu+")->Mass();
  Int_t npoints = points->GetNPoints();
  AliTrackPoint point;
  Int_t naccept=0;
  //
  for (Int_t ipoint=0;ipoint<npoints;ipoint++){
    points->GetPoint(point,ipoint);
    Float_t xyz[3];
    point.GetXYZ(xyz);
    Double_t r=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    if (r<350)  continue;
    if (r>400)  continue;
    AliTracker::PropagateTrackToBxByBz(&pTPC,r,mass,2.,kTRUE);
    AliTracker::PropagateTrackToBxByBz(&pTPC,r,mass,0.1,kTRUE);    
    AliTrackPoint lpoint = point.Rotate(pTPC.GetAlpha());
    pTPC.PropagateTo(lpoint.GetX(),fMagF);
    pTOF=pTPC;
    ((Double_t*)pTOF.GetParameter())[0] =lpoint.GetY();
    ((Double_t*)pTOF.GetParameter())[1] =lpoint.GetZ();
    ((Double_t*)pTOF.GetCovariance())[0]+=3.*3./12.;
    ((Double_t*)pTOF.GetCovariance())[2]+=3.*3./12.;
    ((Double_t*)pTOF.GetCovariance())[5]+=0.1*0.1;
    ((Double_t*)pTOF.GetCovariance())[9]+=0.1*0.1;
    naccept++;
  }
  if (naccept==0) return;  // no tof match clusters
  //
  // 0. Apply standard cuts
  //
  if (track->GetTPCNcls()<kMinTPC) return;  // minimal amount of clusters cut
  // exclude crossing track
  if (friendTrack->GetTPCOut()->GetZ()*track->GetInnerParam()->GetZ()<0)   return;
  //
  if (TMath::Abs(pTOF.GetY()-pTPC.GetY())    >kMaxDy)    return;
  if (TMath::Abs(pTOF.GetSnp()-pTPC.GetSnp())>kMaxAngle) return;
  if (TMath::Abs(pTOF.GetTgl()-pTPC.GetTgl())>kMaxAngle) return;
  //
  // 1. Update median and RMS info
  //
  TVectorD vecDelta(5),vecMedian(5), vecRMS(5);
  TVectorD vecDeltaN(5);
  Double_t sign=(pTOF.GetParameter()[1]>0)? 1.:-1.;
  vecDelta[4]=0;
  for (Int_t i=0;i<4;i++){
    vecDelta[i]=(pTOF.GetParameter()[i]-pTPC.GetParameter()[i])*sign;
    kgdP[i][kglast%kN]=vecDelta[i];
  }
  kglast=(kglast+1);
  Int_t entries=(kglast<kN)?kglast:kN;
  Bool_t isOK=kTRUE;
  for (Int_t i=0;i<4;i++){
    vecMedian[i] = TMath::Median(entries,kgdP[i]);
    vecRMS[i]    = TMath::RMS(entries,kgdP[i]);
    vecDeltaN[i] = 0;
    if (vecRMS[i]>0.){
      vecDeltaN[i] = (vecDelta[i]-vecMedian[i])/(vecRMS[i]+1.);
      vecDeltaN[4]+= TMath::Abs(vecDeltaN[i]);  //sum of abs residuals
      if (TMath::Abs(vecDeltaN[i])>kSigmaCut) isOK=kFALSE;
    }
  }
  //
  // 2. Apply median+-rms cut
  //
  if (kglast<10)  return;   //median and RMS to be defined
  if (!isOK) return;
  //
  // 3. Update alignment
  //
  //Int_t htime = fTime/3600; //time in hours
  Int_t htime = (Int_t)(fTime-fTimeKalmanBin)/fTimeKalmanBin; //time bin
  if (fAlignTOFTPC->GetEntriesFast()<htime){
    fAlignTOFTPC->Expand(htime*2+20);
  }
  AliRelAlignerKalman* align =  (AliRelAlignerKalman*)fAlignTOFTPC->At(htime);
  if (!align){
    // make Alignment object if doesn't exist
    align=new AliRelAlignerKalman(); 
    align->SetRunNumber(fRun);
    (*align->GetStateCov())(6,6)=kVdErr*kVdErr;
    (*align->GetStateCov())(7,7)=kT0Err*kT0Err;
    (*align->GetStateCov())(8,8)=kVdYErr*kVdYErr;
    align->SetOutRejSigma(kOutCut+kOutCut*kN);
    align->SetRejectOutliers(kFALSE);
    align->SetTPCvd(AliTPCcalibDB::Instance()->GetParameters()->GetDriftV()/1000000.);
    align->SetMagField(fMagF); 
    fAlignTOFTPC->AddAt(align,htime);
  }
  align->AddTrackParams(&pTOF,&pTPC);
  Float_t dca[2],cov[3];
  track->GetImpactParameters(dca,cov);
  if (TMath::Abs(dca[0])<kMaxDy){
    FillResHistoTPCTOF(&pTPC,&pTOF);
  }
  //align->SetTimeStamp(fTime);
  Double_t averageTime =  fTime;
  if (align->GetTimeStamp()>0 && align->GetNUpdates()>0) {
    averageTime = (((Double_t)fTime) + ((Double_t)align->GetTimeStamp())*align->GetNUpdates()) / (align->GetNUpdates() + 1.);
    //printf("align->GetTimeStamp() %d, align->GetNUpdates() %d \n", align->GetTimeStamp(), align->GetNUpdates());
  }
  align->SetTimeStamp((Int_t)averageTime);

  //printf("fTime %d, averageTime %d \n", fTime, (Int_t)averageTime);

  align->SetRunNumber(fRun );
  //
  Int_t nupdates=align->GetNUpdates();
  align->SetOutRejSigma(kOutCut+kOutCut*kN/Double_t(nupdates));
  align->SetRejectOutliers(kFALSE);
  TTreeSRedirector *cstream = GetDebugStreamer();  
  if (cstream && align->GetState() && align->GetState()->GetNrows()>2 ){
    TVectorD gpTPC(3), gdTPC(3);
    TVectorD gpTOF(3), gdTOF(3);
    pTPC.GetXYZ(gpTPC.GetMatrixArray());
    pTPC.GetDirection(gdTPC.GetMatrixArray());
    pTOF.GetXYZ(gpTOF.GetMatrixArray());
    pTOF.GetDirection(gdTOF.GetMatrixArray());
    (*cstream)<<"toftpc"<<
      "run="<<fRun<<              //  run number
      "event="<<fEvent<<          //  event number
      "time="<<fTime<<            //  time stamp of event
      "trigger="<<fTrigger<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      //
      "nmed="<<kglast<<        // number of entries to define median and RMS
      "vMed.="<<&vecMedian<<    // median of deltas
      "vRMS.="<<&vecRMS<<       // rms of deltas
      "vDelta.="<<&vecDelta<<   // delta in respect to median
      "vDeltaN.="<<&vecDeltaN<< // normalized delta in respect to median
      "t.="<<track<<            // ful track - find proper cuts
      "a.="<<align<<            // current alignment
      "pTOF.="<<&pTOF<<         // track param TOF
      "pTPC.="<<&pTPC<<         // track param TPC
      "gpTPC.="<<&gpTPC<<       // global position  TPC
      "gdTPC.="<<&gdTPC<<       // global direction TPC
      "gpTOF.="<<&gpTOF<<       // global position  TOF
      "gdTOF.="<<&gdTOF<<       // global position  TOF
      "\n";
  }
}


void  AliTPCcalibTime::BookDistortionMaps(){
  //
  //   Book ndimensional histograms of distortions/residuals
  //   Only primary tracks are selected for analysis
  //
 
  Double_t xminTrack[5], xmaxTrack[5];
  Int_t binsTrack[5];
  TString axisName[5];
  TString axisTitle[5];
  //
  binsTrack[0]  =50;
  axisName[0]   ="#Delta";
  axisTitle[0]  ="#Delta";
  //
  binsTrack[1] =44;
  xminTrack[1] =-1.1; xmaxTrack[1]=1.1;
  axisName[1]  ="tanTheta";
  axisTitle[1]  ="tan(#Theta)";
  //
  binsTrack[2] =180;
  xminTrack[2] =-TMath::Pi(); xmaxTrack[2]=TMath::Pi(); 
  axisName[2]  ="phi";
  axisTitle[2]  ="#phi";
  //
  binsTrack[3] =20;
  xminTrack[3] =-1.; xmaxTrack[3]=1.;   // 0.33 GeV cut 
  axisName[3]  ="snp";
  axisTitle[3]  ="snp";
  //
  binsTrack[4] =10;
  xminTrack[4] =120.; xmaxTrack[4]=215.;   // crossing radius for CE only 
  axisName[4]  ="r";
  axisTitle[4] ="r(cm)";
  //
  // delta y
  xminTrack[0] =-1.5; xmaxTrack[0]=1.5;  // 
  fResHistoTPCCE[0] = new THnSparseS("TPCCE#Delta_{Y} (cm)","#Delta_{Y} (cm)",    5, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCITS[0] = new THnSparseS("TPCITS#Delta_{Y} (cm)","#Delta_{Y} (cm)",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCvertex[0]    = new THnSparseS("TPCVertex#Delta_{Y} (cm)","#Delta_{Y} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-1.5; xmaxTrack[0]=1.5;  // 
  fResHistoTPCTRD[0] = new THnSparseS("TPCTRD#Delta_{Y} (cm)","#Delta_{Y} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5; xmaxTrack[0]=5;  // 
  fResHistoTPCTOF[0] = new THnSparseS("TPCTOF#Delta_{Y} (cm)","#Delta_{Y} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  //
  // delta z
  xminTrack[0] =-3.; xmaxTrack[0]=3.;  // 
  fResHistoTPCCE[1] = new THnSparseS("TPCCE#Delta_{Z} (cm)","#Delta_{Z} (cm)",    5, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCITS[1] = new THnSparseS("TPCITS#Delta_{Z} (cm)","#Delta_{Z} (cm)",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCvertex[1]    = new THnSparseS("TPCVertex#Delta_{Z} (cm)","#Delta_{Z} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTRD[1] = new THnSparseS("TPCTRD#Delta_{Z} (cm)","#Delta_{Z} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  xminTrack[0] =-5.; xmaxTrack[0]=5.;  // 
  fResHistoTPCTOF[1] = new THnSparseS("TPCTOF#Delta_{Z} (cm)","#Delta_{Z} (cm)", 4, binsTrack,xminTrack, xmaxTrack);
  //
  // delta snp-P2
  xminTrack[0] =-0.015; xmaxTrack[0]=0.015;  // 
  fResHistoTPCCE[2] = new THnSparseS("TPCCE#Delta_{#phi}","#Delta_{#phi}",    5, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCITS[2] = new THnSparseS("TPCITS#Delta_{#phi}","#Delta_{#phi}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCvertex[2] = new THnSparseS("TPCITSv#Delta_{#phi}","#Delta_{#phi}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTRD[2] = new THnSparseS("TPCTRD#Delta_{#phi}","#Delta_{#phi}", 4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTOF[2] = new THnSparseS("TPCTOF#Delta_{#phi}","#Delta_{#phi}", 4, binsTrack,xminTrack, xmaxTrack);
  //
  // delta theta-P3
  xminTrack[0] =-0.025; xmaxTrack[0]=0.025;  // 
  fResHistoTPCCE[3] = new THnSparseS("TPCCE#Delta_{#theta}","#Delta_{#theta}",    5, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCITS[3] = new THnSparseS("TPCITS#Delta_{#theta}","#Delta_{#theta}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCvertex[3] = new THnSparseS("TPCITSv#Delta_{#theta}","#Delta_{#theta}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTRD[3] = new THnSparseS("TPCTRD#Delta_{#theta}","#Delta_{#theta}", 4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTOF[3] = new THnSparseS("TPCTOF#Delta_{#theta}","#Delta_{#theta}", 4, binsTrack,xminTrack, xmaxTrack);
  //
  // delta theta-P4
  xminTrack[0] =-0.2; xmaxTrack[0]=0.2;  // 
  fResHistoTPCCE[4] = new THnSparseS("TPCCE#Delta_{1/pt}","#Delta_{1/pt}",    5, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCITS[4] = new THnSparseS("TPCITS#Delta_{1/pt}","#Delta_{1/pt}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCvertex[4] = new THnSparseS("TPCITSv#Delta_{1/pt}","#Delta_{1/pt}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTRD[4] = new THnSparseS("TPCTRD#Delta_{1/pt}","#Delta_{1/pt}",    4, binsTrack,xminTrack, xmaxTrack);
  fResHistoTPCTOF[4] = new THnSparseS("TPCTOF#Delta_{1/pt}","#Delta_{1/pt}",    4, binsTrack,xminTrack, xmaxTrack);
  //
  for (Int_t ivar=0;ivar<4;ivar++){
    for (Int_t ivar2=0;ivar2<5;ivar2++){      
      fResHistoTPCCE[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      fResHistoTPCCE[ivar]->GetAxis(ivar2)->SetTitle(axisTitle[ivar2].Data());
      if (ivar2<4){
	fResHistoTPCITS[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
	fResHistoTPCITS[ivar]->GetAxis(ivar2)->SetTitle(axisTitle[ivar2].Data());
	fResHistoTPCTRD[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
	fResHistoTPCTRD[ivar]->GetAxis(ivar2)->SetTitle(axisTitle[ivar2].Data());
	fResHistoTPCvertex[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
	fResHistoTPCvertex[ivar]->GetAxis(ivar2)->SetTitle(axisTitle[ivar2].Data());
      }
    }
  }
  //
  // Book vertex: time histograms
  //
  Int_t    binsVertex[2]={500, fTimeBins};
  Double_t    aminVertex[2]={-5,fTimeStart};
  Double_t    amaxVertex[2]={5, fTimeEnd};
  const char* hnames[12]={"TPCXAside", "TPCXCside","TPCXACdiff","TPCXAPCdiff",
			  "TPCYAside", "TPCYCside","TPCYACdiff","TPCYAPCdiff",
			  "TPCZAPCside", "TPCZAMCside","TPCZACdiff","TPCZAPCdiff"}; 
  const char* anames[12]={"x (cm) - A side ", "x (cm) - C side","#Delta_{x} (cm) - TPC-A-C","#Delta_{x} (cm) - TPC-Common",
			  "y (cm) - A side ", "y (cm) - C side","#Delta_{x} (cm) - TPC-A-C","#Delta_{y} (cm) - TPC-Common",
			  "z (cm)", "#Delta_{Z} (cm) A-C side","#Delta_{x} (cm) - TPC-A-C","#Delta_{Z} (cm) TPC-common"}; 
  for (Int_t ihis=0; ihis<12; ihis++) {
    if (ihis>=8) aminVertex[0]=-20.;
    if (ihis>=8) amaxVertex[0]=20.;
    fTPCVertex[ihis]=new THnSparseF(hnames[ihis],hnames[ihis],2,binsVertex,aminVertex,amaxVertex);
    fTPCVertex[ihis]->GetAxis(1)->SetTitle("Time");
    fTPCVertex[ihis]->GetAxis(0)->SetTitle(anames[ihis]);
  }
  
  Int_t    binsVertexC[2]={40, 300};
  Double_t aminVertexC[2]={-20,-30};
  Double_t amaxVertexC[2]={20,30};
  const char* hnamesC[5]={"TPCA_TPC","TPCC_TPC","TPCA_ITS","TPCC_ITS","TPC_ITS"};
  for (Int_t ihis=0; ihis<5; ihis++) {
    fTPCVertexCorrelation[ihis]=new THnSparseF(hnamesC[ihis],hnamesC[ihis],2,binsVertexC,aminVertexC,amaxVertexC);
    fTPCVertexCorrelation[ihis]->GetAxis(1)->SetTitle("z (cm)");
    fTPCVertexCorrelation[ihis]->GetAxis(0)->SetTitle("z (cm)");
  }
}


void        AliTPCcalibTime::FillResHistoTPCCE(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pTPCOut ){
  //
  // fill residual histograms pTPCOut-pTPCin - trac crossing CE
  // Histogram 
  //
  if (fMemoryMode<2) return;
  Double_t histoX[5];
  Double_t xyz[3];
  pTPCIn->GetXYZ(xyz);
  Double_t phi= TMath::ATan2(xyz[1],xyz[0]);
  histoX[1]= pTPCIn->GetTgl();
  histoX[2]= phi;
  histoX[3]= pTPCIn->GetSnp();
  histoX[4]= pTPCIn->GetX();
  AliExternalTrackParam lout(*pTPCOut);
  lout.Rotate(pTPCIn->GetAlpha());
  lout.PropagateTo(pTPCIn->GetX(),fMagF);
  //
  for (Int_t ihisto=0; ihisto<5; ihisto++){
    histoX[0]=lout.GetParameter()[ihisto]-pTPCIn->GetParameter()[ihisto];
    fResHistoTPCCE[ihisto]->Fill(histoX);
  }
}  
void        AliTPCcalibTime::FillResHistoTPCITS(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pITSOut ){
  //
  // fill residual histograms pTPCIn-pITSOut
  // Histogram is filled only for primary tracks
  //
  Double_t histoX[4];
  Double_t xyz[3];
  pTPCIn->GetXYZ(xyz);
  Double_t phi= TMath::ATan2(xyz[1],xyz[0]);
  histoX[1]= pTPCIn->GetTgl();
  histoX[2]= phi;
  histoX[3]= pTPCIn->GetSnp();
  AliExternalTrackParam lits(*pITSOut);
  lits.Rotate(pTPCIn->GetAlpha());
  lits.PropagateTo(pTPCIn->GetX(),fMagF);
  //
  for (Int_t ihisto=0; ihisto<5; ihisto++){
    histoX[0]=pTPCIn->GetParameter()[ihisto]-lits.GetParameter()[ihisto];
    fResHistoTPCITS[ihisto]->Fill(histoX);
  }
}  

     
void        AliTPCcalibTime::FillResHistoTPC(const AliESDtrack * pTrack){
  //
  // fill residual histograms pTPC - vertex
  // Histogram is filled only for primary tracks
  //
  if (fMemoryMode<2) return;
  Double_t histoX[4];
  const AliExternalTrackParam * pTPCIn = pTrack->GetInnerParam();
  AliExternalTrackParam pTPCvertex(*(pTrack->GetInnerParam()));
  //
  AliExternalTrackParam lits(*pTrack);
  if (TMath::Abs(pTrack->GetY())>3) return;  // beam pipe
  pTPCvertex.Rotate(lits.GetAlpha());
  //pTPCvertex.PropagateTo(pTPCvertex->GetX(),fMagF);
  AliTracker::PropagateTrackToBxByBz(&pTPCvertex,lits.GetX(),0.1,2,kFALSE);
  AliTracker::PropagateTrackToBxByBz(&pTPCvertex,lits.GetX(),0.1,0.1,kFALSE);
  Double_t xyz[3];
  pTPCIn->GetXYZ(xyz);
  Double_t phi= TMath::ATan2(xyz[1],xyz[0]);
  histoX[1]= pTPCIn->GetTgl();
  histoX[2]= phi;
  histoX[3]= pTPCIn->GetSnp();
  //
  Float_t dca[2], cov[3];
  pTrack->GetImpactParametersTPC(dca,cov);
  for (Int_t ihisto=0; ihisto<5; ihisto++){
    histoX[0]=pTPCvertex.GetParameter()[ihisto]-lits.GetParameter()[ihisto];
    //    if (ihisto<2) histoX[0]=dca[ihisto];
    fResHistoTPCvertex[ihisto]->Fill(histoX);
  }
}


void        AliTPCcalibTime::FillResHistoTPCTRD(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTRDIn ){
  //
  // fill resuidual histogram TPCout-TRDin
  //
  if (fMemoryMode<2) return;
  Double_t histoX[4];
  Double_t xyz[3];
  pTPCOut->GetXYZ(xyz);
  Double_t phi= TMath::ATan2(xyz[1],xyz[0]);
  histoX[1]= pTPCOut->GetTgl();
  histoX[2]= phi;
  histoX[3]= pTPCOut->GetSnp();
  //
  AliExternalTrackParam ltrd(*pTRDIn);
  ltrd.Rotate(pTPCOut->GetAlpha());
  //  ltrd.PropagateTo(pTPCOut->GetX(),fMagF);
  AliTracker::PropagateTrackToBxByBz(&ltrd,pTPCOut->GetX(),0.1,0.1,kFALSE);

  for (Int_t ihisto=0; ihisto<5; ihisto++){
    histoX[0]=pTPCOut->GetParameter()[ihisto]-ltrd.GetParameter()[ihisto];
    fResHistoTPCTRD[ihisto]->Fill(histoX);
  }

}

void        AliTPCcalibTime::FillResHistoTPCTOF(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTOFIn ){
  //
  // fill resuidual histogram TPCout-TOFin
  // track propagated to the TOF position
  if (fMemoryMode<2) return;
  Double_t histoX[4];
  Double_t xyz[3];

  AliExternalTrackParam ltpc(*pTPCOut);
  ltpc.Rotate(pTOFIn->GetAlpha());
  AliTracker::PropagateTrackToBxByBz(&ltpc,pTOFIn->GetX(),0.1,0.1,kFALSE);
  //
  ltpc.GetXYZ(xyz);
  Double_t phi= TMath::ATan2(xyz[1],xyz[0]);
  histoX[1]= ltpc.GetTgl();
  histoX[2]= phi;
  histoX[3]= ltpc.GetSnp();
  //
  for (Int_t ihisto=0; ihisto<2; ihisto++){
    histoX[0]=ltpc.GetParameter()[ihisto]-pTOFIn->GetParameter()[ihisto];
    fResHistoTPCTOF[ihisto]->Fill(histoX);
  }

}
