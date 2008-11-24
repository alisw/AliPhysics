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

//--------------------------------------------------------------------//
//                                                                    //
// AliTOFtrackerV1 Class                                              //
// Task: Perform association of the ESD tracks to TOF Clusters        //
// and Update ESD track with associated TOF Cluster parameters        //
//                                                                    //
// -- Authors : S. Arcelli, C. Zampolli (Bologna University and INFN) //
// -- Contacts: Annalisa.De.Caro@cern.ch                              //
// --         : Chiara.Zampolli@bo.infn.it                            //
// --         : Silvia.Arcelli@bo.infn.it                             //
//                                                                    //
//--------------------------------------------------------------------//

#include <Rtypes.h>
#include <TROOT.h>

#include <TClonesArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSeqCollection.h>

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"

#include "AliTOFRecoParam.h"
#include "AliTOFReconstructor.h"
#include "AliTOFcluster.h"
#include "AliTOFGeometry.h"
#include "AliTOFtrackerV1.h"
#include "AliTOFtrack.h"
#include "AliTOFpidESD.h"

extern TROOT *gROOT;

ClassImp(AliTOFtrackerV1)

//_____________________________________________________________________________
AliTOFtrackerV1::AliTOFtrackerV1():
  fRecoParam(0x0),
  fPid(0x0),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fngoodmatch(0),
  fnbadmatch(0),
  fnunmatch(0),
  fnmatch(0),
  fTracks(new TClonesArray("AliTOFtrack")),
  fSeeds(new TClonesArray("AliESDtrack")),
  fHDigClusMap(0x0),
  fHDigNClus(0x0),
  fHDigClusTime(0x0),
  fHDigClusToT(0x0),
  fHRecNClus(0x0),
  fHRecChi2(0x0),
  fHRecDistZ(0x0),
  fHRecSigYVsP(0x0),
  fHRecSigZVsP(0x0),
  fHRecSigYVsPWin(0x0),
  fHRecSigZVsPWin(0x0)
 { 
  //AliTOFtrackerV1 main Ctor

   InitCheckHists();

}
//_____________________________________________________________________________
AliTOFtrackerV1::AliTOFtrackerV1(const AliTOFtrackerV1 &t):
  AliTracker(),
  fRecoParam(0x0),
  fPid(0x0),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fngoodmatch(0),
  fnbadmatch(0),
  fnunmatch(0),
  fnmatch(0),
  fTracks(new TClonesArray("AliTOFtrack")),
  fSeeds(new TClonesArray("AliESDtrack")),
  fHDigClusMap(0x0),
  fHDigNClus(0x0),
  fHDigClusTime(0x0),
  fHDigClusToT(0x0),
  fHRecNClus(0x0),
  fHRecChi2(0x0),
  fHRecDistZ(0x0),
  fHRecSigYVsP(0x0),
  fHRecSigZVsP(0x0),
  fHRecSigYVsPWin(0x0),
  fHRecSigZVsPWin(0x0)
 { 
  //AliTOFtrackerV1 copy Ctor

  fNseeds=t.fNseeds;
  fNseeds=t.fNseeds;
  fNseedsTOF=t.fNseedsTOF;
  fngoodmatch=t.fngoodmatch;
  fnbadmatch=t.fnbadmatch;
  fnunmatch=t.fnunmatch;
  fnmatch=t.fnmatch;
  fRecoParam=t.fRecoParam;
  fPid=t.fPid;
  fSeeds=t.fSeeds;
  fTracks=t.fTracks;
  fN=t.fN;
}

//_____________________________________________________________________________
AliTOFtrackerV1& AliTOFtrackerV1::operator=(const AliTOFtrackerV1 &t)
{ 
  //AliTOFtrackerV1 assignment operator

  this->fNseeds=t.fNseeds;
  this->fNseedsTOF=t.fNseedsTOF;
  this->fngoodmatch=t.fngoodmatch;
  this->fnbadmatch=t.fnbadmatch;
  this->fnunmatch=t.fnunmatch;
  this->fnmatch=t.fnmatch;
  this->fRecoParam = t.fRecoParam;
  this->fPid = t.fPid;
  this->fSeeds=t.fSeeds;
  this->fTracks=t.fTracks;
  this->fN=t.fN;
  return *this;

}
//_____________________________________________________________________________
AliTOFtrackerV1::~AliTOFtrackerV1() {
  //
  // Dtor
  //

  SaveCheckHists();

  if(!(AliCDBManager::Instance()->GetCacheFlag())){
    delete fRecoParam;
  }
  delete fPid; 
  delete fHDigClusMap;
  delete fHDigNClus;
  delete fHDigClusTime;
  delete fHDigClusToT;
  delete fHRecNClus;
  delete fHRecChi2;
  delete fHRecDistZ;
  delete fHRecSigYVsP;
  delete fHRecSigZVsP;
  delete fHRecSigYVsPWin;
  delete fHRecSigZVsPWin;
  if (fTracks){
    fTracks->Delete();
    delete fTracks;
    fTracks=0x0;
  }
  if (fSeeds){
    fSeeds->Delete();
    delete fSeeds;
    fSeeds=0x0;
  }
}
//_____________________________________________________________________________
Int_t AliTOFtrackerV1::PropagateBack(AliESDEvent* event) {
  //
  // Gets seeds from ESD event and Match with TOF Clusters
  //

  // initialize RecoParam for current event

  AliInfo("Initializing params for TOF... ");

  fRecoParam = AliTOFReconstructor::GetRecoParam();  // instantiate reco param from STEER...

  if (fRecoParam == 0x0) { 
    AliFatal("No Reco Param found for TOF!!!");
  }
  //fRecoParam->Dump();
  //if(fRecoParam->GetApplyPbPbCuts())fRecoParam=fRecoParam->GetPbPbparam();
  //fRecoParam->PrintParameters();

  Double_t parPID[2];   
  parPID[0]=fRecoParam->GetTimeResolution();
  parPID[1]=fRecoParam->GetTimeNSigma();
  fPid=new AliTOFpidESD(parPID);

  //Initialise some counters

  fNseeds=0;
  fNseedsTOF=0;
  fngoodmatch=0;
  fnbadmatch=0;
  fnunmatch=0;
  fnmatch=0;

  Int_t ntrk=event->GetNumberOfTracks();
  fNseeds = ntrk;
  TClonesArray &aESDTrack = *fSeeds;


  //Load ESD tracks into a local Array of ESD Seeds

  for (Int_t i=0; i<fNseeds; i++) {
    AliESDtrack *t=event->GetTrack(i);
    new(aESDTrack[i]) AliESDtrack(*t);
  }

  //Prepare ESD tracks candidates for TOF Matching
  CollectESD();

  //Matching Step
  MatchTracks();

  AliInfo(Form("Number of matched tracks: %d",fnmatch));
  AliInfo(Form("Number of good matched tracks: %d",fngoodmatch));
  AliInfo(Form("Number of bad  matched tracks: %d",fnbadmatch));

  //Update the matched ESD tracks

  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    AliESDtrack *seed =(AliESDtrack*)fSeeds->UncheckedAt(i);
    if(seed->GetTOFsignal()>0){
      t->SetTOFsignal(seed->GetTOFsignal());
      t->SetTOFcluster(seed->GetTOFcluster());
      t->SetTOFsignalToT(seed->GetTOFsignalToT());
      t->SetTOFsignalRaw(seed->GetTOFsignalRaw());
      t->SetTOFsignalDz(seed->GetTOFsignalDz());
      t->SetTOFCalChannel(seed->GetTOFCalChannel());
      Int_t tlab[3]; seed->GetTOFLabel(tlab);    
      t->SetTOFLabel(tlab);
      AliTOFtrack *track = new AliTOFtrack(*seed); 
      t->UpdateTrackParams(track,AliESDtrack::kTOFout);   
      delete track;
    }
  }

  //Handle Time Zero information

  Double_t timeZero=0.;
  Double_t timeZeroMax=99999.;
  Bool_t usetimeZero     = fRecoParam->UseTimeZero();
  Bool_t timeZeroFromT0  = fRecoParam->GetTimeZerofromT0();
  Bool_t timeZeroFromTOF = fRecoParam->GetTimeZerofromTOF();

  AliDebug(1,Form("Use Time Zero?: %d",usetimeZero));
  AliDebug(1,Form("Time Zero from T0? : %d",timeZeroFromT0));
  AliDebug(1,Form("Time Zero From TOF? : %d",timeZeroFromTOF));

  if(usetimeZero){
    if(timeZeroFromT0){
      timeZero=GetTimeZerofromT0(event); 
    }
    if(timeZeroFromTOF && (timeZero>timeZeroMax || !timeZeroFromT0)){
      timeZero=GetTimeZerofromTOF(event); 
    }
  }
  AliDebug(2,Form("time Zero used in PID: %f",timeZero));
  //Make TOF PID
  fPid->MakePID(event,timeZero);

  fSeeds->Clear();
  fTracks->Clear();
  return 0;
  
}
//_________________________________________________________________________
void AliTOFtrackerV1::CollectESD() {
   //prepare the set of ESD tracks to be matched to clusters in TOF

  Int_t seedsTOF1=0;
  Int_t seedsTOF2=0;
 
  TClonesArray &aTOFTrack = *fTracks;
  for (Int_t i=0; i<fNseeds; i++) {

    AliESDtrack *t =(AliESDtrack*)fSeeds->UncheckedAt(i);
    if ((t->GetStatus()&AliESDtrack::kTPCout)==0)continue;

    // TRD 'good' tracks, already propagated at 372 cm

    AliTOFtrack *track = new AliTOFtrack(*t); // New
    Double_t x = track->GetX(); //New

    if (((t->GetStatus()&AliESDtrack::kTRDout)!=0 ) && 
	 ( x >= AliTOFGeometry::RinTOF()) ){
      track->SetSeedIndex(i);
      t->UpdateTrackParams(track,AliESDtrack::kTOFout);    
      new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
      fNseedsTOF++;
      seedsTOF1++;
      delete track;
    }

    // Propagate the rest of TPCbp  

    else {
      if(track->PropagateToInnerTOF()){ 
      	track->SetSeedIndex(i);
	t->UpdateTrackParams(track,AliESDtrack::kTOFout);    
 	new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	fNseedsTOF++;
	seedsTOF2++;
      }
      delete track;
    }
  }

  AliInfo(Form("Number of TOF seeds %i",fNseedsTOF));
  AliInfo(Form("Number of TOF seeds Type 1 %i",seedsTOF1));
  AliInfo(Form("Number of TOF seeds Type 2 %i",seedsTOF2));

  // Sort according uncertainties on track position 
  fTracks->Sort();

}
//_________________________________________________________________________
void AliTOFtrackerV1::MatchTracks( ){
  //
  //Match ESD tracks to clusters in TOF
  //


  // Parameters regulating the reconstruction
  Float_t dY=AliTOFGeometry::XPad(); 
  Float_t dZ=AliTOFGeometry::ZPad(); 

  const Int_t kncmax = 100;
  Float_t sensRadius = fRecoParam->GetSensRadius();
  Float_t scaleFact   = fRecoParam->GetWindowScaleFact();
  Float_t dyMax=fRecoParam->GetWindowSizeMaxY(); 
  Float_t dzMax=fRecoParam->GetWindowSizeMaxZ();
  Double_t maxChi2=fRecoParam->GetMaxChi2();
  Bool_t timeWalkCorr    = fRecoParam->GetTimeWalkCorr();
  AliDebug(1,"++++++++++++++TOF Reconstruction Parameters:++++++++++++ \n");
  AliDebug(1,Form("TOF sens radius: %f",sensRadius));
  AliDebug(1,Form("TOF Window scale factor: %f",scaleFact));
  AliDebug(1,Form("TOF Window max dy: %f",dyMax));
  AliDebug(1,Form("TOF Window max dz: %f",dzMax));
  AliDebug(1,Form("TOF Max Chi2: %f",maxChi2));
  AliDebug(1,Form("Time Walk Correction? : %d",timeWalkCorr));   


  //The matching loop

  for (Int_t iseed=0; iseed<fNseedsTOF; iseed++) {

    AliTOFtrack *track =(AliTOFtrack*)fTracks->UncheckedAt(iseed);
    AliESDtrack *t =(AliESDtrack*)fSeeds->UncheckedAt(track->GetSeedIndex());
    if(t->GetTOFsignal()>0. ) continue;
    AliTOFtrack *trackTOFin =new AliTOFtrack(*track);
     
    // Determine a window around the track
    Double_t x,par[5]; trackTOFin->GetExternalParameters(x,par);
    Double_t cov[15]; trackTOFin->GetExternalCovariance(cov);

    Double_t z    = par[1];   
    Double_t dz   =  scaleFact*3.*TMath::Sqrt(cov[2]+dZ*dZ/12.);
    Double_t dphi =  scaleFact*3.*TMath::Sqrt(cov[0]+dY*dY/12.)/sensRadius; 

    Double_t phi=TMath::ATan2(par[0],x) + trackTOFin->GetAlpha();
    if (phi<-TMath::Pi())phi+=2*TMath::Pi();
    if (phi>=TMath::Pi())phi-=2*TMath::Pi();

    //upper limit on window's size.
    if(dz> dzMax) dz=dzMax;
    if(dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;


    // find the clusters inside the selected window 
    Int_t nc=0;
    AliTOFcluster *clusters[kncmax]; // pointers to the clusters in the window
    Int_t index[kncmax];//to keep track of the cluster index
    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {  
      AliTOFcluster *c=fClusters[k];
      //      if(nc>kncmax)break; /* R+ fix (buffer overflow) */
      if(nc>=kncmax)break; /* R+ fix (buffer overflow protection) */
      if(c->GetZ() > z+dz) break;
      if(c->IsUsed()) continue;      
      if(!c->GetStatus()) {
	      AliDebug(1,"Cluster in channel declared bad!");
	      continue; // skip bad channels as declared in OCDB  
      }
      Float_t xyz[3]; c->GetGlobalXYZ(xyz);
      Double_t clPhi=TMath::ATan2(xyz[1],xyz[0]);
      Double_t dph=TMath::Abs(clPhi-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;
      clusters[nc]=c;
      index[nc] = k;      
      nc++;  
    }

    //start propagation: go to the average TOF pad middle plane at ~379.5 cm

    Float_t  xTOF = sensRadius;
    Double_t ymax = xTOF*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
    Bool_t skip = kFALSE;
    Double_t ysect = trackTOFin->GetYat(xTOF,skip);
    if (skip) break;
    if (ysect > ymax) {
      if (!trackTOFin->Rotate(AliTOFGeometry::GetAlpha())) {
	break;
      }
    } else if (ysect <-ymax) {
      if (!trackTOFin->Rotate(AliTOFGeometry::GetAlpha())) {
	break;
      }
    }
    if(!trackTOFin->PropagateTo(xTOF)) {
      break;
    }


    AliTOFcluster *bestCluster=0;
    Double_t bestChi2=maxChi2; 
    Int_t idclus=-1;
    //    for (Int_t i=0; i<nc; i++){ /* R+ fix (unsafe) */
    for (Int_t i=0; i<nc && i<kncmax; i++){ /* R+ fix (buffer overflow protection) */
      AliTOFcluster *c=clusters[i];  // one of the preselected clusters     
      Double_t chi2=trackTOFin->GetPredictedChi2((AliCluster3D*)c); 
      if (chi2 >= bestChi2) continue;
      bestChi2=chi2;
      bestCluster=c;
      idclus=index[i];
    }
    
    if (!bestCluster) {  // no matching , go to the next track 
      fnunmatch++;
      delete trackTOFin;
      continue;
    }

    fnmatch++;

    AliDebug(2, Form("%7i     %7i     %10i     %10i  %10i  %10i      %7i",
		     iseed,
		     fnmatch-1,
		     TMath::Abs(trackTOFin->GetLabel()),
		     bestCluster->GetLabel(0), 
		     bestCluster->GetLabel(1), 
		     bestCluster->GetLabel(2),
		     idclus)); // AdC

    bestCluster->Use(); 
    if (
	(bestCluster->GetLabel(0)==TMath::Abs(trackTOFin->GetLabel()))
	||
	(bestCluster->GetLabel(1)==TMath::Abs(trackTOFin->GetLabel()))
	||
	(bestCluster->GetLabel(2)==TMath::Abs(trackTOFin->GetLabel()))
	) {
      fngoodmatch++;
       AliDebug(2,Form(" track label good %5i",trackTOFin->GetLabel()));

    }
    else{
      fnbadmatch++;
      AliDebug(2,Form(" track label bad %5i",trackTOFin->GetLabel()));
    }

    //Propagate the track to the best matched cluster
    trackTOFin->PropagateTo(bestCluster);


    //now take the local distance in Z from the pad center for time walk correction
    Float_t tiltangle = AliTOFGeometry::GetAngles(bestCluster->GetDetInd(1),bestCluster->GetDetInd(2))*TMath::DegToRad();
    Double_t dzTW=trackTOFin->GetZ()-bestCluster->GetZ(); // in cm
    dzTW/=TMath::Cos(tiltangle);

    //update the ESD track and delete the TOFtrack
    t->UpdateTrackParams(trackTOFin,AliESDtrack::kTOFout);    
    delete trackTOFin;

    //  Store quantities to be used in the TOF Calibration
    Float_t tToT=AliTOFGeometry::ToTBinWidth()*bestCluster->GetToT()*1E-3; // in ns
    t->SetTOFsignalToT(tToT);
    Float_t rawTime=AliTOFGeometry::TdcBinWidth()*bestCluster->GetTDCRAW()+32; // RAW time,in ps
    t->SetTOFsignalRaw(rawTime);
    t->SetTOFsignalDz(dzTW);
    AliDebug(2,Form(" Setting TOF raw time: %f  z distance: %f time: %f = ",rawTime,dzTW));    
    Int_t ind[5];
    ind[0]=bestCluster->GetDetInd(0);
    ind[1]=bestCluster->GetDetInd(1);
    ind[2]=bestCluster->GetDetInd(2);
    ind[3]=bestCluster->GetDetInd(3);
    ind[4]=bestCluster->GetDetInd(4);
    Int_t calindex = AliTOFGeometry::GetIndex(ind);
    t->SetTOFCalChannel(calindex);

    // keep track of the track labels in the matched cluster
    Int_t tlab[3];
    tlab[0]=bestCluster->GetLabel(0);
    tlab[1]=bestCluster->GetLabel(1);
    tlab[2]=bestCluster->GetLabel(2);
    AliDebug(2,Form(" tdc time of the matched track %i = ",bestCluster->GetTDC()));    
    Double_t tof=AliTOFGeometry::TdcBinWidth()*bestCluster->GetTDC()+32; // in ps
    AliDebug(2,Form(" tof time of the matched track: %f = ",tof));    
    Double_t tofcorr=tof;
    if(timeWalkCorr)tofcorr=CorrectTimeWalk(dzTW,tof);
    AliDebug(2,Form(" tof time of the matched track, after TW corr: %f = ",tofcorr));    
    //Set TOF time signal and pointer to the matched cluster
    t->SetTOFsignal(tofcorr);
    t->SetTOFcluster(idclus); // pointing to the recPoints tree
    t->SetTOFLabel(tlab);

    Double_t mom=t->GetP();
    // Fill Reco-QA histos for Reconstruction
    fHRecNClus->Fill(nc);
    fHRecChi2->Fill(bestChi2);
    fHRecDistZ->Fill(dzTW);
    fHRecSigYVsP->Fill(mom,TMath::Sqrt(cov[0]));
    fHRecSigZVsP->Fill(mom,TMath::Sqrt(cov[2]));
    fHRecSigYVsPWin->Fill(mom,dphi*sensRadius);
    fHRecSigZVsPWin->Fill(mom,dz);

    // Fill Tree for on-the-fly offline Calibration
    // no longer there - all info is in the ESDs now

  }

}
//_________________________________________________________________________
Int_t AliTOFtrackerV1::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------

  Int_t npadX = AliTOFGeometry::NpadX();
  Int_t npadZ = AliTOFGeometry::NpadZ();
  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  TBranch *branch=cTree->GetBranch("TOF");
  if (!branch) { 
    AliError("can't get the branch with the TOF clusters !");
    return 1;
  }

  static TClonesArray dummy("AliTOFcluster",10000);
  dummy.Clear();
  TClonesArray *clusters=&dummy;
  branch->SetAddress(&clusters);

  cTree->GetEvent(0);
  Int_t nc=clusters->GetEntriesFast();
  fHDigNClus->Fill(nc);

  AliInfo(Form("Number of clusters: %d",nc));

  for (Int_t i=0; i<nc; i++) {
    AliTOFcluster *c=(AliTOFcluster*)clusters->UncheckedAt(i);
    fClusters[i]=new AliTOFcluster(*c); fN++;

  // Fill Digits QA histos
 
    Int_t isector = c->GetDetInd(0);
    Int_t iplate = c->GetDetInd(1);
    Int_t istrip = c->GetDetInd(2);
    Int_t ipadX = c->GetDetInd(4);
    Int_t ipadZ = c->GetDetInd(3);

    Float_t time =(AliTOFGeometry::TdcBinWidth()*c->GetTDC())*1E-3; // in ns
    Float_t tot = (AliTOFGeometry::TdcBinWidth()*c->GetToT())*1E-3;//in ns
 
    Int_t stripOffset = 0;
    switch (iplate) {
    case 0:
      stripOffset = 0;
      break;
    case 1:
      stripOffset = nStripC;
      break;
    case 2:
      stripOffset = nStripC+nStripB;
      break;
    case 3:
      stripOffset = nStripC+nStripB+nStripA;
      break;
    case 4:
      stripOffset = nStripC+nStripB+nStripA+nStripB;
      break;
    default:
      AliError(Form("Wrong plate number in TOF (%d) !",iplate));
      break;
    };
    Int_t zindex=npadZ*(istrip+stripOffset)+(ipadZ+1);
    Int_t phiindex=npadX*isector+ipadX+1;
    fHDigClusMap->Fill(zindex,phiindex);
    fHDigClusTime->Fill(time);
    fHDigClusToT->Fill(tot);
  }


  return 0;
}
//_________________________________________________________________________
void AliTOFtrackerV1::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads TOF clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) {
    delete fClusters[i];
    fClusters[i] = 0x0;
  }
  fN=0;
}

//_________________________________________________________________________
Int_t AliTOFtrackerV1::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  //MOD
  //Here we need to get the Z in the tracking system

  if (fN==0) return 0;
  if (z <= fClusters[0]->GetZ()) return 0;
  if (z > fClusters[fN-1]->GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}

//_________________________________________________________________________
Bool_t AliTOFtrackerV1::GetTrackPoint(Int_t index, AliTrackPoint& p) const
{
  // Get track space point with index i
  // Coordinates are in the global system
  AliTOFcluster *cl = fClusters[index];
  Float_t xyz[3];
  cl->GetGlobalXYZ(xyz);
  Float_t phi=TMath::ATan2(xyz[1],xyz[0]);
  Float_t phiangle = (Int_t(phi*TMath::RadToDeg()/20.)+0.5)*20.*TMath::DegToRad();
  Float_t sinphi = TMath::Sin(phiangle), cosphi = TMath::Cos(phiangle);
  Float_t tiltangle = AliTOFGeometry::GetAngles(cl->GetDetInd(1),cl->GetDetInd(2))*TMath::DegToRad();
  Float_t sinth = TMath::Sin(tiltangle), costh = TMath::Cos(tiltangle);
  Float_t sigmay2 = AliTOFGeometry::XPad()*AliTOFGeometry::XPad()/12.;
  Float_t sigmaz2 = AliTOFGeometry::ZPad()*AliTOFGeometry::ZPad()/12.;
  Float_t cov[6];
  cov[0] = sinphi*sinphi*sigmay2 + cosphi*cosphi*sinth*sinth*sigmaz2;
  cov[1] = -sinphi*cosphi*sigmay2 + sinphi*cosphi*sinth*sinth*sigmaz2;
  cov[2] = -cosphi*sinth*costh*sigmaz2;
  cov[3] = cosphi*cosphi*sigmay2 + sinphi*sinphi*sinth*sinth*sigmaz2;
  cov[4] = -sinphi*sinth*costh*sigmaz2;
  cov[5] = costh*costh*sigmaz2;
  p.SetXYZ(xyz[0],xyz[1],xyz[2],cov);

  // Detector numbering scheme
  Int_t nSector = AliTOFGeometry::NSectors();
  Int_t nPlate  = AliTOFGeometry::NPlates();
  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  Int_t isector = cl->GetDetInd(0);
  if (isector >= nSector)
    AliError(Form("Wrong sector number in TOF (%d) !",isector));
  Int_t iplate = cl->GetDetInd(1);
  if (iplate >= nPlate)
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
  Int_t istrip = cl->GetDetInd(2);

  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
    break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC+nStripB;
    break;
  case 3:
    stripOffset = nStripC+nStripB+nStripA;
    break;
  case 4:
    stripOffset = nStripC+nStripB+nStripA+nStripB;
    break;
  default:
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
    break;
  };

  Int_t idet = (2*(nStripC+nStripB)+nStripA)*isector +
               stripOffset +
               istrip;
  UShort_t volid = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,idet);
  p.SetVolumeID((UShort_t)volid);
  return kTRUE;
}
//_________________________________________________________________________
void AliTOFtrackerV1::InitCheckHists() {

  //Init histos for Digits/Reco QA and Calibration

  TDirectory *dir = gDirectory;
  TFile *logFileTOF = 0;

  TSeqCollection *list = gROOT->GetListOfFiles();
  int n = list->GetEntries();
  Bool_t isThere=kFALSE;
  for(int i=0; i<n; i++) {
    logFileTOF = (TFile*)list->At(i);
    if (strstr(logFileTOF->GetName(), "TOFQA.root")){
      isThere=kTRUE;
      break;
    } 
  }

  if(!isThere)logFileTOF = new TFile( "TOFQA.root","RECREATE");
  logFileTOF->cd(); 

  //Digits "QA" 
  fHDigClusMap = new TH2F("TOFDig_ClusMap", "",182,0.5,182.5,864, 0.5,864.5);  
  fHDigNClus = new TH1F("TOFDig_NClus", "",200,0.5,200.5);  
  fHDigClusTime = new TH1F("TOFDig_ClusTime", "",2000,0.,200.);  
  fHDigClusToT = new TH1F("TOFDig_ClusToT", "",500,0.,100);  

  //Reco "QA"
  fHRecNClus =new TH1F("TOFRec_NClusW", "",50,0.5,50.5);
  fHRecDistZ=new TH1F("TOFRec_DistZ", "",50,0.5,10.5);
  fHRecChi2=new TH1F("TOFRec_Chi2", "",100,0.,10.);
  fHRecSigYVsP=new TH2F("TOFDig_SigYVsP", "",40,0.,4.,100, 0.,5.);
  fHRecSigZVsP=new TH2F("TOFDig_SigZVsP", "",40,0.,4.,100, 0.,5.);
  fHRecSigYVsPWin=new TH2F("TOFDig_SigYVsPWin", "",40,0.,4.,100, 0.,50.);
  fHRecSigZVsPWin=new TH2F("TOFDig_SigZVsPWin", "",40,0.,4.,100, 0.,50.);

  dir->cd();

}

//_________________________________________________________________________
void AliTOFtrackerV1::SaveCheckHists() {

  //write histos for Digits/Reco QA and Calibration

  TDirectory *dir = gDirectory;
  TFile *logFile = 0;
  TFile *logFileTOF = 0;

  TSeqCollection *list = gROOT->GetListOfFiles();
  int n = list->GetEntries();
  for(int i=0; i<n; i++) {
    logFile = (TFile*)list->At(i);
    if (strstr(logFile->GetName(), "AliESDs.root")) break;
  }

  Bool_t isThere=kFALSE;
  for(int i=0; i<n; i++) {
    logFileTOF = (TFile*)list->At(i);
    if (strstr(logFileTOF->GetName(), "TOFQA.root")){
      isThere=kTRUE;
      break;
    } 
  }
   
  if(!isThere) {
	  AliError(Form("File TOFQA.root not found!! not wring histograms...."));
	  return;
  }
  logFile->cd();
  fHDigClusMap->Write(fHDigClusMap->GetName(), TObject::kOverwrite);
  fHDigNClus->Write(fHDigNClus->GetName(), TObject::kOverwrite);
  fHDigClusTime->Write(fHDigClusTime->GetName(), TObject::kOverwrite);
  fHDigClusToT->Write(fHDigClusToT->GetName(), TObject::kOverwrite);
  fHRecNClus->Write(fHRecNClus->GetName(), TObject::kOverwrite);
  fHRecChi2->Write(fHRecChi2->GetName(), TObject::kOverwrite);
  fHRecDistZ->Write(fHRecDistZ->GetName(), TObject::kOverwrite);
  fHRecSigYVsP->Write(fHRecSigYVsP->GetName(), TObject::kOverwrite);
  fHRecSigZVsP->Write(fHRecSigZVsP->GetName(), TObject::kOverwrite);
  fHRecSigYVsPWin->Write(fHRecSigYVsPWin->GetName(), TObject::kOverwrite);
  fHRecSigZVsPWin->Write(fHRecSigZVsPWin->GetName(), TObject::kOverwrite);
  logFile->Flush();  

  dir->cd();

  }
//_________________________________________________________________________
Float_t AliTOFtrackerV1::CorrectTimeWalk( Float_t dist, Float_t tof) {

  //dummy, for the moment
  Float_t tofcorr=0.;
  if(dist<AliTOFGeometry::ZPad()*0.5){
    tofcorr=tof;
    //place here the actual correction
  }else{
    tofcorr=tof; 
  } 
  return tofcorr;
}
//_________________________________________________________________________
Float_t AliTOFtrackerV1::GetTimeZerofromT0(AliESDEvent *event) const {

  //Returns TimeZero as measured by T0 detector

  return event->GetT0();
}
//_________________________________________________________________________
Float_t AliTOFtrackerV1::GetTimeZerofromTOF(AliESDEvent * /*event*/) const {

  //dummy, for the moment. T0 algorithm using tracks on TOF
  {
    //place T0 algo here...
  }
  return 0.;
}
//_________________________________________________________________________

void AliTOFtrackerV1::FillClusterArray(TObjArray* arr) const
{
  //
  // Returns the TOF cluster array
  //

  if (fN==0)
    arr = 0x0;
  else
    for (Int_t i=0; i<fN; ++i) arr->Add(fClusters[i]);

}
