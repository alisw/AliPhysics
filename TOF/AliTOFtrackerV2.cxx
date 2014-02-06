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
// AliTOFtrackerV2 Class                                              //
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
#include <TGeoManager.h>
#include <TTree.h>

#include "AliGeomManager.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDTOFcluster.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"

#include "AliTOFRecoParam.h"
#include "AliTOFReconstructor.h"
#include "AliTOFcluster.h"
#include "AliTOFGeometry.h"
#include "AliTOFtrackerV2.h"
#include "AliTOFtrack.h"

extern TGeoManager *gGeoManager;

ClassImp(AliTOFtrackerV2)

//_____________________________________________________________________________
AliTOFtrackerV2::AliTOFtrackerV2():
  fkRecoParam(0x0),
  fGeom(0x0),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fnunmatch(0),
  fnmatch(0),
  fTracks(new TClonesArray("AliTOFtrack")),
  fSeeds(new TObjArray(100)),
  fClusters(0x0)
{
  //AliTOFtrackerV2 main Ctor

  // Getting the geometry
  fGeom = new AliTOFGeometry();

}
//_____________________________________________________________________________
AliTOFtrackerV2::~AliTOFtrackerV2() {
  //
  // Dtor
  //

  if(!(AliCDBManager::Instance()->GetCacheFlag())){
    delete fkRecoParam;
  }
  delete fGeom; 
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

  if(fClusters){
    delete[] fClusters;
    fClusters = NULL;
  }

}
//_____________________________________________________________________________
void AliTOFtrackerV2::GetPidSettings(AliESDpid *esdPID) {
  // 
  // Sets TOF resolution from RecoParams
  //
  if (fkRecoParam)
    esdPID->GetTOFResponse().SetTimeResolution(fkRecoParam->GetTimeResolution());
  else
    AliWarning("fkRecoParam not yet set; cannot set PID settings");
} 
//_____________________________________________________________________________
Int_t AliTOFtrackerV2::PropagateBack(AliESDEvent * const event) {
  //
  // Gets seeds from ESD event and Match with TOF Clusters
  //

  //Update the matched ESD tracks
  // needed in case of call of TOF info before of the selection of matching and in case of no clusters available at all
  if(fN==0)
    event->SetTOFcluster(1,fClusters); 
  else
    event->SetTOFcluster(fN,fClusters);

  if (fN==0) {
    AliInfo("No TOF recPoints to be matched with reconstructed tracks");
    return 0;
  }

  // initialize RecoParam for current event
  AliDebug(1,"Initializing params for TOF");

  fkRecoParam = AliTOFReconstructor::GetRecoParam();  // instantiate reco param from STEER...

  if (fkRecoParam == 0x0) { 
    AliFatal("No Reco Param found for TOF!!!");
  }

  //Initialise some counters

  fNseeds=0;
  fNseedsTOF=0;
  fnunmatch=0;
  fnmatch=0;

  Int_t ntrk=event->GetNumberOfTracks();
  fNseeds = ntrk;

  //Load ESD tracks into a local Array of ESD Seeds
  for (Int_t i=0; i<fNseeds; i++){
    fSeeds->AddLast(event->GetTrack(i));
    event->GetTrack(i)->SetESDEvent(event);
  }
  //Prepare ESD tracks candidates for TOF Matching
  CollectESD();

  if (fNseeds==0 || fNseedsTOF==0) {
    AliInfo("No seeds to try TOF match");
    return 0;
  }

  // clusterize before of matching
  Clusterize();

  //Second Step with Looser Matching Criterion
  MatchTracks();

  AliInfo(Form("Number of matched tracks = %d",fnmatch));

  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    //t->SetESDEvent(event);
    AliESDtrack *seed =(AliESDtrack*)fSeeds->At(i);

    if ( (seed->GetStatus()&AliESDtrack::kTOFin)!=0 ) {
      t->SetStatus(AliESDtrack::kTOFin);
      if ( (seed->GetStatus()&AliESDtrack::kTOFout)!=0 ) {
	t->SetStatus(AliESDtrack::kTOFout);
	//t->SetTOFclusterArray(seed->GetNTOFclusters(),seed->GetTOFclusterArray());
	t->SortTOFcluster();

	// Make attention, please:
	//      AliESDtrack::fTOFInfo array does not be stored in the AliESDs.root file
	//      it is there only for a check during the reconstruction step.
	Float_t info[10]; seed->GetTOFInfo(info);
	t->SetTOFInfo(info);
	AliDebug(3,Form(" distance=%f; residual in the pad reference frame: dX=%f, dZ=%f", info[0],info[1],info[2]));

	/*
	Double_t alphaA = (Double_t)t->GetAlpha();
	Double_t xA = (Double_t)t->GetX();
	Double_t yA = (Double_t)t->GetY();
	Double_t zA = (Double_t)t->GetZ();
	Double_t p1A = (Double_t)t->GetSnp();
	Double_t p2A = (Double_t)t->GetTgl();
	Double_t p3A = (Double_t)t->GetSigned1Pt();
	const Double_t *covA = (Double_t*)t->GetCovariance();

	// Check done:
	//       by calling the AliESDtrack::UpdateTrackParams,
	//       the current track parameters are changed
	//       and it could cause refit problems.
	//       We need to update only the following track parameters:
        //            the track length and expected times.
	//       Removed AliESDtrack::UpdateTrackParams call
	//       Called AliESDtrack::SetIntegratedTimes(...) and
	//       AliESDtrack::SetIntegratedLength() routines.

	AliTOFtrack *track = new AliTOFtrack(*seed);
	t->UpdateTrackParams(track,AliESDtrack::kTOFout); // to be checked - AdC
	delete track;
	Double_t time[AliPID::kSPECIESC]; t->GetIntegratedTimes(time);
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
	*/
	Double_t time[AliPID::kSPECIESC]; t->GetIntegratedTimes(time,AliPID::kSPECIESC);
	AliDebug(2,Form(" TOF params: %6d  %f %f %f %f %f %6d %3d %f",
			i,
			t->GetTOFsignalRaw(),t->GetTOFsignal(),t->GetTOFsignalToT(),
			t->GetTOFsignalDz(),t->GetTOFsignalDx(),t->GetTOFCalChannel(),
			t->GetTOFcluster(),t->GetIntegratedLength()));
	AliDebug(2,Form("  %f %f %f %f %f %f %f %f %f",
			time[0], time[1], time[2], time[3], time[4], time[5], time[6], time[7], time[8]));
      }
    }
  }
  
  fSeeds->Clear();
  fTracks->Clear();
  
  AliInfo(Form("Number of cluster to be checked = %d",fN));
  if(fN){
    Int_t *matchmap = new Int_t[fN];
    event->SetTOFcluster(fN,fClusters,matchmap);
    for (Int_t i=0; i<ntrk; i++) { // remapping after TOF matching selection
      AliESDtrack *t=event->GetTrack(i);
      t->ReMapTOFcluster(fN,matchmap);
    }

    delete[] matchmap;
  }
  
  

  return 0;
  
}
//_________________________________________________________________________
void AliTOFtrackerV2::CollectESD() {
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
    if ( ( (t->GetStatus()&AliESDtrack::kTRDout)!=0 ) ) {

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
      //delete track;
    }

    else { // Propagate the rest of TPCbp

      AliDebug(1,Form(" Before propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));

      if ( track->PropagateToInnerTOF() ) { 

	AliDebug(1,Form(" TPC propagated track till rho = %fcm."
			" And then the track has been propagated till rho = %fcm.",
			x, (Float_t)track->GetX()));

      	track->SetSeedIndex(i);
	t->UpdateTrackParams(track,AliESDtrack::kTOFin);
 	new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	fNseedsTOF++;
	seedsTOF2++;

	AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track->GetIntegratedLength()));
      }
      delete track;
    }
  }

  AliInfo(Form("Number of TOF seeds = %d (kTRDout371 = %d, kTRDoutLess371 = %d, !kTRDout = %d)",fNseedsTOF,seedsTOF1,seedsTOF3,seedsTOF2));

  // Sort according uncertainties on track position 
  fTracks->Sort();

}

//_________________________________________________________________________
void AliTOFtrackerV2::MatchTracks() {
  //
  //Match ESD tracks to clusters in TOF
  //

  // Parameters used/regulating the reconstruction
  static Float_t detDepth=18.;
  static Float_t padDepth=0.5;

  const Float_t kSpeedOfLight= 2.99792458e-2; // speed of light [cm/ps]

  Float_t dY=AliTOFGeometry::XPad(); 
  Float_t dZ=AliTOFGeometry::ZPad(); 

  Float_t sensRadius = fkRecoParam->GetSensRadius();
  Float_t stepSize   = fkRecoParam->GetStepSize();
  Float_t scaleFact  = fkRecoParam->GetWindowScaleFact();
  Float_t dyMax=fkRecoParam->GetWindowSizeMaxY(); 
  Float_t dzMax=fkRecoParam->GetWindowSizeMaxZ();
  Float_t dCut=10.;//fkRecoParam->GetDistanceCut(); // This is to be loaded by OCDB. It should be 10cm always.
  Double_t maxChi2=fkRecoParam->GetMaxChi2TRD();
  Bool_t timeWalkCorr = fkRecoParam->GetTimeWalkCorr();
  AliDebug(1,"++++++++++++++TOF Reconstruction Parameters:++++++++++++++");
  AliDebug(1,Form("TOF sens radius: %f",sensRadius));
  AliDebug(1,Form("TOF step size: %f",stepSize));
  AliDebug(1,Form("TOF Window scale factor: %f",scaleFact));
  AliDebug(1,Form("TOF Window max dy: %f",dyMax));
  AliDebug(1,Form("TOF Window max dz: %f",dzMax));
  AliDebug(1,Form("TOF distance Cut: %f",dCut));
  AliDebug(1,Form("TOF Max Chi2: %f",maxChi2));
  AliDebug(1,Form("Time Walk Correction? : %d",timeWalkCorr));   

  //Match ESD tracks to clusters in TOF

  // Get the number of propagation steps
  Int_t nSteps=(Int_t)(detDepth/stepSize);
  AliDebug(1,Form(" Number of steps to be done %d",nSteps));

  AliDebug(1,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

  //PH Arrays (moved outside of the loop)
  Float_t * trackPos[4];
  for (Int_t ii=0; ii<4; ii++) trackPos[ii] = new Float_t[nSteps];
  Int_t * clind = new Int_t[fN];
  
  // Some init
  const Int_t kNclusterMax = 1000; // related to fN value
  TGeoHMatrix global[kNclusterMax];

  //The matching loop
  for (Int_t iseed=0; iseed<fNseedsTOF; iseed++) {

    for (Int_t ii=0; ii<fN; ii++) clind[ii]=-1;
    for (Int_t ii=0; ii<kNclusterMax; ii++) global[ii] = 0x0;
    for (Int_t ii=0; ii<4; ii++)
      for (Int_t jj=0; jj<nSteps; jj++) trackPos[ii][jj]=0.;

    AliTOFtrack *track =(AliTOFtrack*)fTracks->UncheckedAt(iseed);
    AliESDtrack *t =(AliESDtrack*)fSeeds->At(track->GetSeedIndex());
    AliTOFtrack *trackTOFin = new AliTOFtrack(*track);

    Double_t timesOr[AliPID::kSPECIESC]; t->GetIntegratedTimes(timesOr,AliPID::kSPECIESC); // in ps

    // Determine a window around the track
    Double_t x,par[5]; 
    trackTOFin->GetExternalParameters(x,par);
    Double_t cov[15]; 
    trackTOFin->GetExternalCovariance(cov);

    if (cov[0]<0. || cov[2]<0.) {
      AliWarning(Form("Very strange track (%d)! At least one of its covariance matrix diagonal elements is negative!",iseed));
      delete trackTOFin;
      continue;
    }

    Double_t dphi=
      scaleFact*
      ((5*TMath::Sqrt(TMath::Abs(cov[0])) + 0.5*dY + 2.5*TMath::Abs(par[2]))/sensRadius); 
    Double_t dz=
       scaleFact*
       (5*TMath::Sqrt(TMath::Abs(cov[2])) + 0.5*dZ + 2.5*TMath::Abs(par[3]));

    Double_t phi=TMath::ATan2(par[0],x) + trackTOFin->GetAlpha();
    if (phi<-TMath::Pi())phi+=2*TMath::Pi();
    if (phi>=TMath::Pi())phi-=2*TMath::Pi();
    Double_t z=par[1];   

    //upper limit on window's size.
    if (dz> dzMax) dz=dzMax;
    if (dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;


    // find the clusters in the window of the track
    Int_t nc=0;
    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {

      if (nc>=kNclusterMax) {
 	AliWarning("No more matchable clusters can be stored! Please, increase the corresponding vectors size.");
 	break;
      }

      AliESDTOFcluster *c=&(fClusters[k]);
      if (c->GetZ() > z+dz) break;
      if (!c->GetStatus()) {
	AliDebug(1,"Cluster in channel declared bad!");
	continue; // skip bad channels as declared in OCDB
      }

      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;

      Double_t yc=(c->GetPhi() - trackTOFin->GetAlpha())*c->GetR();
      Double_t p[2]={yc, c->GetZ()};
      Double_t cov2[3]= {dY*dY/12., 0., dZ*dZ/12.};
      if (trackTOFin->AliExternalTrackParam::GetPredictedChi2(p,cov2) > maxChi2)continue;

      clind[nc] = k;      
      Char_t path[200];
      Int_t ind[5]; fGeom->GetVolumeIndices(c->GetTOFchannel(),ind);
      fGeom->GetVolumePath(ind,path);
      gGeoManager->cd(path);
      global[nc] = *gGeoManager->GetCurrentMatrix();
      nc++;
    }

    if (nc == 0 ) {
      AliDebug(1,Form("No available clusters for the track number %d",iseed));
      fnunmatch++;
      delete trackTOFin;
      continue;
    }

    AliDebug(1,Form(" Number of available TOF clusters for the track number %d: %d",iseed,nc));

    //start fine propagation 

    Double_t *times[AliPID::kSPECIESC];
    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
      times[isp] = new Double_t[nSteps];
    }

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<nSteps; istep++){ 
      
      // First of all, propagate the track...
      Float_t xs = AliTOFGeometry::RinTOF()+istep*stepSize;
      if (!(trackTOFin->PropagateTo(xs))) break;

      //  ...and then, if necessary, rotate the track
      Double_t ymax = xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
      Double_t ysect = trackTOFin->GetY();
      if (ysect > ymax) {
	if (!(trackTOFin->Rotate(AliTOFGeometry::GetAlpha()))) break;
      } else if (ysect <-ymax) {
	if (!(trackTOFin->Rotate(-AliTOFGeometry::GetAlpha()))) break;
      }

      Double_t mom = trackTOFin->P();

      if(istep == 0){
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  times[isp][nStepsDone] = stepSize/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }
      else{
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  times[isp][nStepsDone] = times[isp][nStepsDone-1] + (trackTOFin->GetIntegratedLength()-trackPos[3][nStepsDone-1])/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }

      // store the running point (Globalrf) - fine propagation     

      Double_t r[3]; trackTOFin->GetXYZ(r);
      trackPos[0][nStepsDone]= (Float_t) r[0];
      trackPos[1][nStepsDone]= (Float_t) r[1];
      trackPos[2][nStepsDone]= (Float_t) r[2];   
      trackPos[3][nStepsDone]= trackTOFin->GetIntegratedLength();

      nStepsDone++;
      AliDebug(3,Form(" current step %d (%d) - nStepsDone=%d",istep,nSteps,nStepsDone));
    }

    if ( nStepsDone == 0 ) {
      AliDebug(1,Form(" No track points for track number %d",iseed));
      fnunmatch++;
      delete trackTOFin;
      continue;
    }

    AliDebug(3,Form(" Number of steps done for the track number %d: %d",iseed,nStepsDone));

    Int_t *isClusterMatchable = NULL;
    if(nc){
      isClusterMatchable = new Int_t[nc];
      for (Int_t i=0; i<nc; i++) isClusterMatchable[i] = kFALSE;	  	
    }

    Int_t nfound = 0;
    Bool_t accept = kFALSE;
    Bool_t isInside = kFALSE;
    for (Int_t istep=0; istep<nStepsDone; istep++) {

      Bool_t gotInsideCluster = kFALSE;
      Int_t trackInsideCluster = -1;

      Float_t ctrackPos[3];     
      ctrackPos[0] = trackPos[0][istep];
      ctrackPos[1] = trackPos[1][istep];
      ctrackPos[2] = trackPos[2][istep];

      //now see whether the track matches any of the TOF clusters            

      Float_t dist3d[3]={0.,0.,0.};
      accept = kFALSE;
     
      for (Int_t i=0; i<nc; i++) {

        // ***** NEW *****
        /* check whether track was inside another cluster
         * and in case inhibit this cluster.
         * this will allow to only go on and add track points for
         * that cluster where the track got inside first */
        if (gotInsideCluster && trackInsideCluster != i) {
	  AliDebug(3,Form(" A - istep=%d ~ %d %d ~ nfound=%d",istep,trackInsideCluster,i,nfound));
          continue;
	}
	AliDebug(3,Form(" B - istep=%d ~ %d %d ~ nfound=%d",istep,trackInsideCluster,i,nfound));

        /* check whether track is inside this cluster */
	for (Int_t hh=0; hh<3; hh++) dist3d[hh]=0.;
	isInside = fGeom->IsInsideThePad((TGeoHMatrix*)(&global[i]),ctrackPos,dist3d);

        // ***** NEW *****
        /* if track is inside this cluster set flags which will then
         * inhibit to add track points for the other clusters */
        if (isInside) {
          gotInsideCluster = kTRUE;
          trackInsideCluster = i;
        }

	Float_t yLoc = dist3d[1];
	Float_t rLoc = TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[2]*dist3d[2]);
	accept = (TMath::Abs(yLoc)<padDepth*0.5 && rLoc<dCut);

	//***** NEW *****
	/* add point everytime that:
	 * - the track is inside the cluster
	 * - the track got inside the cluster, even when it eventually exited the cluster
	 * - the tracks is within dCut from the cluster
	 */
        if (accept || isInside || gotInsideCluster) {

	  Double_t timesCurrent[AliPID::kSPECIESC];
	  AliDebug(3,Form(" Momentum for track %d -> %f", iseed,t->P()));
	  for (Int_t j=0;j<AliPID::kSPECIESC;j++) {
	    timesCurrent[j] = timesOr[j] + times[j][istep];
	  }

	  if (TMath::Abs(dist3d[1])<stepSize && !isClusterMatchable[i]) {
	    isClusterMatchable[i] = kTRUE;
	    fClusters[clind[i]].Update(t->GetID(),dist3d[1],dist3d[0],dist3d[2],trackPos[3][istep],timesCurrent);//x,y,z -> tracking RF
	    t->AddTOFcluster(clind[i]);
	    t->SetStatus(AliESDtrack::kTOFout);
	  }
          AliDebug(2,Form(" dist3dLoc[0] = %f, dist3dLoc[1] = %f, dist3dLoc[2] = %f ",dist3d[0],dist3d[1],dist3d[2]));

          nfound++;

	  AliDebug(3,Form(" C - istep=%d ~ %d %d ~ nfound=%d",istep,trackInsideCluster,i,nfound));
        
          // ***** NEW *****
        }//end if accept
        
      } //end for on the clusters
    } //end for on the steps     
    if(nc) delete[] isClusterMatchable;

    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
      delete[] times[isp];
    }


    if (nfound == 0 ) {
      AliDebug(1,Form(" No matchable track points for the track number %d",iseed));
      fnunmatch++;
      delete trackTOFin;
      continue;
    }

    AliDebug(1,Form(" Number of track points for the track number %d: %d",iseed,nfound));

    Int_t nMatchedClusters = t->GetNTOFclusters();
 
    if (nMatchedClusters==0) {
      AliDebug(1,Form("Reconstructed track %d doesn't match any TOF cluster", iseed));
      fnunmatch++;
      delete trackTOFin;
      continue;
    }

    AliDebug(1,Form(" %d - matched (%d)",track->GetSeedIndex()/*iseed*/,nMatchedClusters));

    fnmatch++;

    /*
    AliTOFcluster cTOF = AliTOFcluster(volIdClus,
    (Float_t)posClus[0],(Float_t)posClus[1],(Float_t)posClus[2],
    (Float_t)covClus[0],(Float_t)covClus[1],(Float_t)covClus[2],
    (Float_t)covClus[3],(Float_t)covClus[4],(Float_t)covClus[5],
    tofLabels,volIndices,parClu,kTRUE,index[i]);

    // Fill the track residual histograms.
    FillResiduals(trackTOFin,c,kFALSE);
    */

    delete trackTOFin;

  } // loop on fSeeds

  for (Int_t ii=0; ii<4; ii++) delete [] trackPos[ii];
  delete [] clind;
 
}
//_________________________________________________________________________
Int_t AliTOFtrackerV2::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------

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
  fN=clusters->GetEntriesFast();
  AliInfo(Form("Number of clusters: %d",fN));

  if(fClusters){
    delete[] fClusters;
    fClusters = NULL;
  }

  if(fN)
    fClusters = new AliESDTOFcluster[fN];
  else{
    fClusters = new AliESDTOFcluster[1];
    fN = 1;
    return 0;
  }

  for (Int_t i=0; i<fN; i++) {
    AliTOFcluster *c=(AliTOFcluster*)clusters->UncheckedAt(i);
    Int_t ind[5];
    ind[0]=c->GetDetInd(0);
    ind[1]=c->GetDetInd(1);
    ind[2]=c->GetDetInd(2);
    ind[3]=c->GetDetInd(3);
    ind[4]=c->GetDetInd(4);
    Int_t calindex = AliTOFGeometry::GetIndex(ind);
    Int_t tofLabels[3]={c->GetLabel(0),c->GetLabel(1),c->GetLabel(2)};
    AliESDTOFcluster esdTOFclus(i,calindex,
				AliTOFGeometry::TdcBinWidth()*c->GetTDC()/*ps*/,
				AliTOFGeometry::TdcBinWidth()*c->GetTDCRAW()/*ps*/,
				AliTOFGeometry::ToTBinWidth()*c->GetToT()*1E-3/*ns*/,
				tofLabels,
				c->GetDeltaBC(),c->GetL0L1Latency(),
				c->GetStatus(),c->GetZ(),c->GetPhi(),c->GetR());

    fClusters[i] = esdTOFclus;

  }

  return 0;
}
//_________________________________________________________________________
void AliTOFtrackerV2::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads TOF clusters
  //--------------------------------------------------------------------

  // don't delete TOF clusters here because they should be written
}

//_________________________________________________________________________
Int_t AliTOFtrackerV2::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  if (fN==0) return 0;
  if (z <= fClusters[0].GetZ()) return 0;
  if (z > fClusters[fN-1].GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m].GetZ()) b=m+1;

    else e=m; 
  }
  return m;
}

//_________________________________________________________________________
Bool_t AliTOFtrackerV2::GetTrackPoint(Int_t index, AliTrackPoint& p) const
{
  // Get track space point with index i
  // Coordinates are in the global system
  AliESDTOFcluster *cl = &(fClusters[index]);
  Float_t xyz[3];
  xyz[0] = cl->GetR()*TMath::Cos(cl->GetPhi());
  xyz[1] = cl->GetR()*TMath::Sin(cl->GetPhi());
  xyz[2] = cl->GetZ();
  Float_t phiangle = (Int_t(cl->GetPhi()*TMath::RadToDeg()/20.)+0.5)*20.*TMath::DegToRad();
  Float_t sinphi = TMath::Sin(phiangle), cosphi = TMath::Cos(phiangle);
  Int_t tofChannel=cl->GetTOFchannel();
  Int_t ind[5]; fGeom->GetVolumeIndices(tofChannel,ind);
  Float_t tiltangle = AliTOFGeometry::GetAngles(ind[1],ind[2])*TMath::DegToRad();
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

  Int_t isector = ind[0];//cl->GetDetInd(0);
  if (isector >= nSector)
    AliError(Form("Wrong sector number in TOF (%d) !",isector));
  Int_t iplate = ind[1];//cl->GetDetInd(1);
  if (iplate >= nPlate)
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
  Int_t istrip = ind[2];//cl->GetDetInd(2);

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

void AliTOFtrackerV2::FillClusterArray(TObjArray* arr) const
{
  //
  // Returns the TOF cluster array
  //

  if (fN==0)
    arr = 0x0;
  else
    for (Int_t i=0; i<fN; ++i) arr->Add(&(fClusters[i]));

}
//_________________________________________________________________________
Float_t AliTOFtrackerV2::CorrectTimeWalk( Float_t dist, Float_t tof) const {

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
void AliTOFtrackerV2::Clusterize(){
  Int_t detId[5];
  for(Int_t i=0; i < fN-1;i++){
    AliESDTOFcluster *c1=&(fClusters[i]);
    if(!c1->GetStatus()) continue;

    Int_t chan1 = c1->GetTOFchannel();
    AliTOFGeometry::GetVolumeIndices(chan1, detId); // Get volume index from channel index

    Int_t ieta = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
    if(detId[1]/*module*/ == 0) ieta += 0;
    else if(detId[1] == 1) ieta += 38;
    else if(detId[1] == 2) ieta += 76;
    else if(detId[1] == 3) ieta += 106;
    else if(detId[1] == 4) ieta += 144;
    Int_t iphi = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;
    

    for(Int_t j=i+1; j < fN;j++){
      AliESDTOFcluster *c2=&(fClusters[j]);
      if(!c2->GetStatus()) continue;

      Int_t chan2 = c2->GetTOFchannel();

      // check if the two TOF hits are in the same strip
      if(chan1/96 != chan2/96) continue;

      AliTOFGeometry::GetVolumeIndices(chan2, detId); // Get volume index from channel index
      Int_t ieta2 = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
      if(detId[1]/*module*/ == 0) ieta2 += 0;
      else if(detId[1] == 1) ieta2 += 38;
      else if(detId[1] == 2) ieta2 += 76;
      else if(detId[1] == 3) ieta2 += 106;
      else if(detId[1] == 4) ieta2 += 144;
      Int_t iphi2 = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;
      
      // check if the fired pad are close in space
      if(TMath::Abs(iphi-iphi2)>1 || TMath::Abs(ieta-ieta2)>1) continue;

      // check if the TOF time are close enough to be merged
      if(TMath::Abs(c1->GetTime() - c2->GetTime()) > 500/*in ps*/) continue;

      // merge them
      Int_t label[3] = {c2->GetLabel(0),c2->GetLabel(1),c2->GetLabel(2)};
      fClusters[i].AddTOFhit(c2->GetClusterIndex(),chan2,c2->GetTime(),c2->GetTimeRaw(),c2->GetTOT(),label,
                                       c2->GetDeltaBC(),c2->GetL0L1Latency(),1,c2->GetZ(),c2->GetPhi(),c2->GetR());
      
      c2->SetStatus(0); // only the merged one should be used
      j = fN; // cluster "i" merged go to the next one ("i+1")
    }
  }

  // second step of clusterization
  for(Int_t i=0; i < fN-1;i++){
    AliESDTOFcluster *c1=&(fClusters[i]);
    if(!c1->GetStatus()) continue;

    Int_t chan1 = c1->GetTOFchannel(0);
    AliTOFGeometry::GetVolumeIndices(chan1, detId); // Get volume index from channel index

    Int_t ieta = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
    if(detId[1]/*module*/ == 0) ieta += 0;
    else if(detId[1] == 1) ieta += 38;
    else if(detId[1] == 2) ieta += 76;
    else if(detId[1] == 3) ieta += 106;
    else if(detId[1] == 4) ieta += 144;
    Int_t iphi = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;

    Int_t ieta2,iphi2,chan2=chan1;
    if(c1->GetNTOFhits() > 1){
      chan2 = c1->GetTOFchannel(1);
      AliTOFGeometry::GetVolumeIndices(chan2, detId); // Get volume index from channel index

      ieta2 = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
      if(detId[1]/*module*/ == 0) ieta2 += 0;
      else if(detId[1] == 1) ieta2 += 38;
      else if(detId[1] == 2) ieta2 += 76;
      else if(detId[1] == 3) ieta2 += 106;
      else if(detId[1] == 4) ieta2 += 144;
      iphi2 = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;
    }
    else{
      iphi2=iphi;
      ieta2=ieta;
    }

    for(Int_t j=i+1; j < i;j++){
      AliESDTOFcluster *c2=&(fClusters[j]);
      if(!c2->GetStatus()) continue;

      Int_t chan3 = c2->GetTOFchannel();

      // check if the two TOF hits are in the same strip
      if(chan1/96 != chan3/96) continue;

      AliTOFGeometry::GetVolumeIndices(chan3, detId); // Get volume index from channel index
      Int_t ieta3 = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
      if(detId[1]/*module*/ == 0) ieta3 += 0;
      else if(detId[1] == 1) ieta3 += 38;
      else if(detId[1] == 2) ieta3 += 76;
      else if(detId[1] == 3) ieta3 += 106;
      else if(detId[1] == 4) ieta3 += 144;
      Int_t iphi3 = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;
      
      // check if the fired pad are close in space
      if((TMath::Abs(iphi-iphi3)>1 && TMath::Abs(iphi2-iphi3)>1) || (TMath::Abs(ieta-ieta3)>1 && TMath::Abs(ieta2-ieta3)>1)) 
continue;
      
      // check if the TOF time are close enough to be merged
      if(TMath::Abs(c1->GetTime() - c2->GetTime()) > 500/*in ps*/) continue;
      
      // merge them
      Int_t label[3] = {c2->GetLabel(0),c2->GetLabel(1),c2->GetLabel(2)};
      fClusters[i].AddTOFhit(c2->GetClusterIndex(),chan2,c2->GetTime(),c2->GetTimeRaw(),c2->GetTOT(),label,
                                       c2->GetDeltaBC(),c2->GetL0L1Latency(),1,c2->GetZ(),c2->GetPhi(),c2->GetR());

      if(c2->GetNTOFhits() > 1){ // in case also the second cluster has two hits
        Int_t label2[3] = {c2->GetLabel(0,1),c2->GetLabel(1,1),c2->GetLabel(2,1)};
        fClusters[i].AddTOFhit(c2->GetClusterIndex(1),c2->GetTOFchannel(1),c2->GetTime(1),c2->GetTimeRaw(1),
                                         c2->GetTOT(1),label2,c2->GetDeltaBC(2),c2->GetL0L1Latency(2),1,c2->GetZ(),
                                         c2->GetPhi(),c2->GetR());
      }

      c1->SetStatus(0); // only the merged one should be used
      j = fN; // cluster "i" merged go to the next one ("i+1")
    }
  }  
}

