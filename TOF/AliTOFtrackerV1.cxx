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
#include <TObjArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSeqCollection.h>

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
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

extern TROOT *gROOT;

ClassImp(AliTOFtrackerV1)

//_____________________________________________________________________________
AliTOFtrackerV1::AliTOFtrackerV1():
  fkRecoParam(0x0),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fngoodmatch(0),
  fnbadmatch(0),
  fnunmatch(0),
  fnmatch(0),
  fTracks(new TClonesArray("AliTOFtrack")),
  fSeeds(new TObjArray(100)),
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

   for (Int_t ii=0; ii<kMaxCluster; ii++) fClusters[ii]=0x0;

   InitCheckHists();

}
//_____________________________________________________________________________
AliTOFtrackerV1::~AliTOFtrackerV1() {
  //
  // Dtor
  //

  SaveCheckHists();

  if(!(AliCDBManager::Instance()->GetCacheFlag())){
    delete fkRecoParam;
  }
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


  for (Int_t ii=0; ii<kMaxCluster; ii++)
    if (fClusters[ii]) fClusters[ii]->Delete();

}
//_____________________________________________________________________________
void AliTOFtrackerV1::GetPidSettings(AliESDpid *esdPID) {
  // 
  // Sets TOF resolution from RecoParams
  //
  if (fkRecoParam)
    esdPID->GetTOFResponse().SetTimeResolution(fkRecoParam->GetTimeResolution());
  else
    AliWarning("fkRecoParam not yet set; cannot set PID settings");
}
//_____________________________________________________________________________
Int_t AliTOFtrackerV1::PropagateBack(AliESDEvent * const event) {
  //
  // Gets seeds from ESD event and Match with TOF Clusters
  //

  // initialize RecoParam for current event
  AliDebug(1,"Initializing params for TOF");

  fkRecoParam = AliTOFReconstructor::GetRecoParam();  // instantiate reco param from STEER...

  if (fkRecoParam == 0x0) { 
    AliFatal("No Reco Param found for TOF!!!");
  }
  //fkRecoParam->Dump();
  //if(fkRecoParam->GetApplyPbPbCuts())fkRecoParam=fkRecoParam->GetPbPbparam();
  //fkRecoParam->PrintParameters();

  //Initialise some counters

  fNseeds=0;
  fNseedsTOF=0;
  fngoodmatch=0;
  fnbadmatch=0;
  fnunmatch=0;
  fnmatch=0;

  Int_t ntrk=event->GetNumberOfTracks();
  fNseeds = ntrk;

  //Load ESD tracks into a local Array of ESD Seeds
  for (Int_t i=0; i<fNseeds; i++)
    fSeeds->AddLast(event->GetTrack(i));

  //Prepare ESD tracks candidates for TOF Matching
  CollectESD();

  //Matching Step
  MatchTracks();

  AliInfo(Form("Number of matched tracks = %d (good = %d, bad = %d)",fnmatch,fngoodmatch,fnbadmatch));

  //Update the matched ESD tracks

  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    AliESDtrack *seed =(AliESDtrack*)fSeeds->At(i);

    if ( (seed->GetStatus()&AliESDtrack::kTOFin)!=0 ) {
      t->SetStatus(AliESDtrack::kTOFin);
      //if(seed->GetTOFsignal()>0){
      if ( (seed->GetStatus()&AliESDtrack::kTOFout)!=0 ) {
	t->SetStatus(AliESDtrack::kTOFout);
	t->SetTOFsignal(seed->GetTOFsignal());
	t->SetTOFcluster(seed->GetTOFcluster());
	t->SetTOFsignalToT(seed->GetTOFsignalToT());
	t->SetTOFsignalRaw(seed->GetTOFsignalRaw());
	t->SetTOFsignalDz(seed->GetTOFsignalDz());
	t->SetTOFsignalDx(seed->GetTOFsignalDx());
	t->SetTOFDeltaBC(seed->GetTOFDeltaBC());
	t->SetTOFL0L1(seed->GetTOFL0L1());
	t->SetTOFCalChannel(seed->GetTOFCalChannel());
	Int_t tlab[3]; seed->GetTOFLabel(tlab);
	t->SetTOFLabel(tlab);

	Double_t alphaA = (Double_t)t->GetAlpha();
	Double_t xA = (Double_t)t->GetX();
	Double_t yA = (Double_t)t->GetY();
	Double_t zA = (Double_t)t->GetZ();
	Double_t p1A = (Double_t)t->GetSnp();
	Double_t p2A = (Double_t)t->GetTgl();
	Double_t p3A = (Double_t)t->GetSigned1Pt();
	const Double_t *covA = (Double_t*)t->GetCovariance();

	// Make attention, please:
	//      AliESDtrack::fTOFInfo array does not be stored in the AliESDs.root file
	//      it is there only for a check during the reconstruction step.
	Float_t info[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	seed->GetTOFInfo(info);
	t->SetTOFInfo(info);
	AliDebug(2,Form(" distance=%f; residual in the pad reference frame: dX=%f, dZ=%f", info[0],info[1],info[2]));

	// Check done:
	//       by calling the AliESDtrack::UpdateTrackParams,
	//       the current track parameters are changed
	//       and it could cause refit problems.
	//       We need to update only the following track parameters:
        //            the track length and expected times.
	//       Removed AliESDtrack::UpdateTrackParams call
	//       Called AliESDtrack::SetIntegratedTimes(...) and
	//       AliESDtrack::SetIntegratedLength() routines.
	/*
	  AliTOFtrack *track = new AliTOFtrack(*seed);
	  t->UpdateTrackParams(track,AliESDtrack::kTOFout); // to be checked - AdC
	  delete track;
	  Double_t time[10]; t->GetIntegratedTimes(time);
	*/

	Double_t time[10]; seed->GetIntegratedTimes(time);
	t->SetIntegratedTimes(time);

	Double_t length =  seed->GetIntegratedLength();
	t->SetIntegratedLength(length);

	Double_t alphaB = (Double_t)t->GetAlpha();
	Double_t xB = (Double_t)t->GetX();
	Double_t yB = (Double_t)t->GetY();
	Double_t zB = (Double_t)t->GetZ();
	Double_t p1B = (Double_t)t->GetSnp();
	Double_t p2B = (Double_t)t->GetTgl();
	Double_t p3B = (Double_t)t->GetSigned1Pt();
	const Double_t *covB = (Double_t*)t->GetCovariance();
	AliDebug(2,"Track params -now(before)-:");
	AliDebug(2,Form("    X: %f(%f), Y: %f(%f), Z: %f(%f) --- alpha: %f(%f)",
			xB,xA,
			yB,yA,
			zB,zA,
			alphaB,alphaA));
	AliDebug(2,Form("    p1: %f(%f), p2: %f(%f), p3: %f(%f)",
			p1B,p1A,
			p2B,p2A,
			p3B,p3A));
	AliDebug(2,Form("    cov1: %f(%f), cov2: %f(%f), cov3: %f(%f)"
			" cov4: %f(%f), cov5: %f(%f), cov6: %f(%f)"
			" cov7: %f(%f), cov8: %f(%f), cov9: %f(%f)"
			" cov10: %f(%f), cov11: %f(%f), cov12: %f(%f)"
			" cov13: %f(%f), cov14: %f(%f), cov15: %f(%f)",
			covB[0],covA[0],
			covB[1],covA[1],
			covB[2],covA[2],
			covB[3],covA[3],
			covB[4],covA[4],
			covB[5],covA[5],
			covB[6],covA[6],
			covB[7],covA[7],
			covB[8],covA[8],
			covB[9],covA[9],
			covB[10],covA[10],
			covB[11],covA[11],
			covB[12],covA[12],
			covB[13],covA[13],
			covB[14],covA[14]
			));
	AliDebug(3,Form(" %6d  %f %f %f %f %f %6d %3d %f  %f %f %f %f %f",
			i,
			t->GetTOFsignalRaw(),
			t->GetTOFsignal(),
			t->GetTOFsignalToT(),
			t->GetTOFsignalDz(),
			t->GetTOFsignalDx(),
			t->GetTOFCalChannel(),
			t->GetTOFcluster(),
			t->GetIntegratedLength(),
			time[0], time[1], time[2], time[3], time[4]
			)
		 );
      }
    }
  }

  fSeeds->Clear();
  fTracks->Delete();
  return 0;
  
}
//_________________________________________________________________________
void AliTOFtrackerV1::CollectESD() {
   //prepare the set of ESD tracks to be matched to clusters in TOF

  Int_t seedsTOF1=0;
  Int_t seedsTOF3=0;
  Int_t seedsTOF2=0;
 
  TClonesArray &aTOFTrack = *fTracks;
  for (Int_t i=0; i<fNseeds; i++) {

    AliESDtrack *t =(AliESDtrack*)fSeeds->At(i);
    if ((t->GetStatus()&AliESDtrack::kTPCout)==0)continue;

    AliTOFtrack *track = new AliTOFtrack(*t); // New
    Float_t x = (Float_t)track->GetX(); //New

    // TRD 'good' tracks
    if ( ( (t->GetStatus()&AliESDtrack::kTRDout)!=0 ) )  {

      AliDebug(1,Form(" Before propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));

      // TRD 'good' tracks, already propagated at 371 cm
      if ( x >= AliTOFGeometry::Rmin() ) {

	if ( track->PropagateToInnerTOF() ) {

	  AliDebug(1,Form(" TRD propagated track till rho = %fcm."
			  " And then the track has been propagated till rho = %fcm.",
			  x, (Float_t)track->GetX()));

	  track->SetSeedIndex(i);
	  t->UpdateTrackParams(track,AliESDtrack::kTOFin);
	  new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	  fNseedsTOF++;
	  seedsTOF1++;

	  AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));
	}
	delete track;

      }
      else { // TRD 'good' tracks, propagated rho<371cm

	if  ( track->PropagateToInnerTOF() ) {

	  AliDebug(1,Form(" TRD propagated track till rho = %fcm."
			  " And then the track has been propagated till rho = %fcm.",
			  x, (Float_t)track->GetX()));

	  track->SetSeedIndex(i);
	  t->UpdateTrackParams(track,AliESDtrack::kTOFin);
	  new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	  fNseedsTOF++;
	  seedsTOF3++;

	  AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));
	}
	delete track;

      }

    }

    else { // Propagate the rest of TPCbp

      AliDebug(1,Form(" Before propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));

      if ( track->PropagateToInnerTOF() ) {

	AliDebug(1,Form(" TRD propagated track till rho = %fcm."
			" And then the track has been propagated till rho = %fcm.",
			x, (Float_t)track->GetX()));

      	track->SetSeedIndex(i);
	t->UpdateTrackParams(track,AliESDtrack::kTOFin);
 	new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	fNseedsTOF++;
	seedsTOF2++;
      }
      delete track;
    }
  }

  AliInfo(Form("Number of TOF seeds = %d (kTRDout371 = %d, kTRDoutLess371 = %d, !kTRDout = %d)",fNseedsTOF,seedsTOF1,seedsTOF3,seedsTOF2));

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

  const Float_t kTimeOffset = 0.; // time offset for tracking algorithm [ps]

  const Int_t kncmax = 100;
  Float_t sensRadius = fkRecoParam->GetSensRadius();
  Float_t scaleFact   = fkRecoParam->GetWindowScaleFact();
  Float_t dyMax=fkRecoParam->GetWindowSizeMaxY(); 
  Float_t dzMax=fkRecoParam->GetWindowSizeMaxZ();
  Double_t maxChi2=fkRecoParam->GetMaxChi2();
  Bool_t timeWalkCorr    = fkRecoParam->GetTimeWalkCorr();
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
    AliESDtrack *t =(AliESDtrack*)fSeeds->At(track->GetSeedIndex());
    //if ( t->GetTOFsignal()>0. ) continue;
    if ( (t->GetStatus()&AliESDtrack::kTOFout)!=0 ) continue;
    AliTOFtrack *trackTOFin = new AliTOFtrack(*track);
     
    // Determine a window around the track
    Double_t x,par[5]; trackTOFin->GetExternalParameters(x,par);
    Double_t cov[15]; trackTOFin->GetExternalCovariance(cov);

    if (cov[0]<0. || cov[2]<0.) {
      AliWarning(Form("Very strange track (%d)! At least one of its covariance matrix diagonal elements is negative!",iseed));
      //delete trackTOFin;
      //continue;
    }

    Double_t z    = par[1];   
    Double_t dz   = scaleFact*3.*TMath::Sqrt(TMath::Abs(cov[2])+dZ*dZ/12.);
    Double_t dphi = scaleFact*3.*TMath::Sqrt(TMath::Abs(cov[0])+dY*dY/12.)/sensRadius; 

    Double_t phi=TMath::ATan2(par[0],x) + trackTOFin->GetAlpha();
    if (phi<-TMath::Pi())phi+=2*TMath::Pi();
    if (phi>=TMath::Pi())phi-=2*TMath::Pi();

    //upper limit on window's size.
    if (dz> dzMax) dz=dzMax;
    if (dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;

    // find the clusters inside the selected window 
    Int_t nc=0;
    AliTOFcluster *clusters[kncmax]; // pointers to the clusters in the window
    Int_t index[kncmax];//to keep track of the cluster index
    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {  
      AliTOFcluster *c=fClusters[k];
      //      if(nc>kncmax)break; /* R+ fix (buffer overflow) */
      if (nc>=kncmax) {
 	AliWarning("No more matchable clusters can be stored! Please, increase the corresponding vectors size.");
	break; /* R+ fix (buffer overflow protection) */
      }
      if (c->GetZ() > z+dz) break;
      if (c->IsUsed()) continue;      
      if (!c->GetStatus()) {
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

    AliDebug(1,Form(" Number of matchable TOF clusters for the track number %d: %d",iseed,nc));

    //start propagation: go to the average TOF pad middle plane at ~379.5 cm

    // First of all, propagate the track...
    Float_t xTOF = sensRadius;
    if (!(trackTOFin->PropagateTo(xTOF))) {
      delete trackTOFin;
      continue;
    }

    // ...and then, if necessary, rotate the track
    Double_t ymax = xTOF*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
    Double_t ysect = trackTOFin->GetY();
    if (ysect > ymax) {
      if (!(trackTOFin->Rotate(AliTOFGeometry::GetAlpha()))) {
	delete trackTOFin;
	continue;
      }
    } else if (ysect <-ymax) {
      if (!(trackTOFin->Rotate(-AliTOFGeometry::GetAlpha()))) {
	delete trackTOFin;
	continue;
      }
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
       AliDebug(2,Form(" track label good %5d",trackTOFin->GetLabel()));

    } else {
      fnbadmatch++;
      AliDebug(2,Form(" track label bad %5d",trackTOFin->GetLabel()));
    }

    //Propagate the track to the best matched cluster
    if (!(trackTOFin->PropagateTo(bestCluster))) {
      delete trackTOFin;
      continue;
    }

    // If necessary, rotate the track
    Double_t yATxMax = trackTOFin->GetX()*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
    Double_t yATx = trackTOFin->GetY();
    if (yATx > yATxMax) {
      if (!(trackTOFin->Rotate(AliTOFGeometry::GetAlpha()))) {
	delete trackTOFin;
	continue;
      }
    } else if (yATx <-yATxMax) {
      if (!(trackTOFin->Rotate(-AliTOFGeometry::GetAlpha()))) {
	delete trackTOFin;
	continue;
      }
    }

    // Fill the track residual histograms.
    FillResiduals(trackTOFin,bestCluster,kFALSE);

    //now take the local distance in Z from the pad center for time walk correction
    Float_t tiltangle = AliTOFGeometry::GetAngles(bestCluster->GetDetInd(1),bestCluster->GetDetInd(2))*TMath::DegToRad();
    Double_t dzTW=trackTOFin->GetZ()-bestCluster->GetZ(); // in cm - in the ALICE RF -
    dzTW/=TMath::Cos(tiltangle); // from ALICE/tracking RF to pad RF (1)
    dzTW=-dzTW; // from ALICE/tracking RF to pad RF (2)
    if (tiltangle!=0.) AliDebug(2,Form(" rho_track = %f --- rho_cluster = %f ",trackTOFin->GetX(),bestCluster->GetX()));

    //update the ESD track and delete the TOFtrack
    t->UpdateTrackParams(trackTOFin,AliESDtrack::kTOFout);

    //  Store quantities to be used in the TOF Calibration
    Float_t tToT=AliTOFGeometry::ToTBinWidth()*bestCluster->GetToT()*1E-3; // in ns
    t->SetTOFsignalToT(tToT);
    Float_t rawTime=AliTOFGeometry::TdcBinWidth()*bestCluster->GetTDCRAW()+kTimeOffset; // RAW time,in ps
    t->SetTOFsignalRaw(rawTime);
    t->SetTOFsignalDz(dzTW);

    Float_t deltaY = trackTOFin->GetY()-bestCluster->GetY();
    t->SetTOFsignalDx(deltaY);

    t->SetTOFDeltaBC(bestCluster->GetDeltaBC());
    t->SetTOFL0L1(bestCluster->GetL0L1Latency());

    Float_t distR = (trackTOFin->GetX()-bestCluster->GetX())*
      (trackTOFin->GetX()-bestCluster->GetX());
    distR+=deltaY*deltaY;
    distR+=dzTW*dzTW;
    distR = TMath::Sqrt(distR);
    Float_t info[10] = {distR, deltaY, dzTW,
			0.,0.,0.,0.,0.,0.,0.};
    t->SetTOFInfo(info);

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
    AliDebug(2,Form(" tdc time of the matched track %6d = ",bestCluster->GetTDC()));    
    Double_t tof=AliTOFGeometry::TdcBinWidth()*bestCluster->GetTDC()+kTimeOffset; // in ps
    AliDebug(2,Form(" tof time of the matched track: %f = ",tof));
    Double_t tofcorr=tof;
    if(timeWalkCorr)tofcorr=CorrectTimeWalk(dzTW,tof);
    AliDebug(2,Form(" tof time of the matched track, after TW corr: %f = ",tofcorr));    
    //Set TOF time signal and pointer to the matched cluster
    t->SetTOFsignal(tofcorr);
    t->SetTOFcluster(idclus); // pointing to the recPoints tree
    t->SetTOFLabel(tlab);

    AliDebug(2,Form(" Setting TOF raw time: %f  z distance: %f  corrected time: %f",rawTime,dzTW,tofcorr));

    Double_t mom=t->GetP();
    AliDebug(2,Form(" Momentum for track %d -> %f", iseed,mom));
    // Fill Reco-QA histos for Reconstruction
    fHRecNClus->Fill(nc);
    fHRecChi2->Fill(bestChi2);
    fHRecDistZ->Fill(dzTW);
    if (cov[0]>=0.)
      fHRecSigYVsP->Fill(mom,TMath::Sqrt(cov[0]));
    else
      fHRecSigYVsP->Fill(mom,-TMath::Sqrt(-cov[0]));
    if (cov[2]>=0.)
      fHRecSigZVsP->Fill(mom,TMath::Sqrt(cov[2]));
    else
      fHRecSigZVsP->Fill(mom,-TMath::Sqrt(-cov[2]));
    fHRecSigYVsPWin->Fill(mom,dphi*sensRadius);
    fHRecSigZVsPWin->Fill(mom,dz);

    // Fill Tree for on-the-fly offline Calibration
    // no longer there - all info is in the ESDs now

    delete trackTOFin;

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
//PH    fClusters[i]=new AliTOFcluster(*c); fN++;
    fClusters[i]=c; fN++;

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
//PH    delete fClusters[i];
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
  //TFile *logFile = 0;
  TFile *logFileTOF = 0;

  TSeqCollection *list = gROOT->GetListOfFiles();
  int n = list->GetEntries();
  /*
  for(int i=0; i<n; i++) {
    logFile = (TFile*)list->At(i);
    if (strstr(logFile->GetName(), "AliESDs.root")) break;
  }
  */
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
  //logFile->cd();
  logFileTOF->cd();
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
  //logFile->Flush();  
  logFileTOF->Flush();  

  dir->cd();

  }
//_________________________________________________________________________
Float_t AliTOFtrackerV1::CorrectTimeWalk( Float_t dist, Float_t tof) const {

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
Float_t AliTOFtrackerV1::GetTimeZerofromT0(const AliESDEvent * const event) const {

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
