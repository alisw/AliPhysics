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
#include "AliESDTOFCluster.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"

#include "AliTOFRecoParam.h"
#include "AliTOFReconstructor.h"
#include "AliTOFGeometry.h"
#include "AliTOFtrackerV2.h"
#include "AliTOFtrack.h"

#include "AliESDTOFHit.h"

extern TGeoManager *gGeoManager;

ClassImp(AliTOFtrackerV2)

//_____________________________________________________________________________
AliTOFtrackerV2::AliTOFtrackerV2():
  fkRecoParam(0x0),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fnunmatch(0),
  fnmatch(0),
  fSeeds(new TObjArray(100)),
  fClustersESD(new TClonesArray("AliESDTOFCluster")),
  fHitsESD(new TClonesArray("AliESDTOFHit")),
  fEvent(0),
  fNsteps(0)
{
  //AliTOFtrackerV2 main Ctor
  for (Int_t ii=0; ii<kMaxCluster; ii++){
    fClusters[ii]=0x0;
    fWrittenInPos[ii] = -1;
  }

  for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++)
    fTimesAr[isp] = NULL;

  for (Int_t ii=0; ii<4; ii++)
    fTrackPos[ii] = NULL;
}
//_____________________________________________________________________________
AliTOFtrackerV2::~AliTOFtrackerV2() {
  //
  // Dtor
  //

  if(!(AliCDBManager::Instance()->GetCacheFlag())){
    delete fkRecoParam;
  }
  if (fSeeds){
    fSeeds->Delete();
    delete fSeeds;
    fSeeds=0x0;
  }

  if (fClustersESD){
    fClustersESD->Delete();
    delete fClustersESD;
    fClustersESD=0x0;
  }
  if (fHitsESD){
    fHitsESD->Delete();
    delete fHitsESD;
    fHitsESD=0x0;
  }

  for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
    if(fTimesAr[isp]) delete[] fTimesAr[isp];
    fTimesAr[isp] = NULL;
  }


  for (Int_t ii=0; ii<4; ii++)
    if(fTrackPos[ii])
      delete [] fTrackPos[ii];
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

  fEvent = event;

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
    fSeeds->Clear();
    fClustersESD->Clear();
    fHitsESD->Clear();
    return 0;
  }

  // clusterize before of matching
  Clusterize(); // fN might change

  //Second Step with Looser Matching Criterion
  MatchTracks(); 
  
  // switch array from ALL to filter for ESD (moving from all fClusterESD to filtered AliESDEvent->GetESDTOFClusters())
  TClonesArray* esdTOFHitArr = event->GetESDTOFHits();
  TClonesArray* esdTOFClArr = event->GetESDTOFClusters();

  AliInfo(Form("TOF before the matching: hits = %i and clusters = %i",esdTOFHitArr->GetEntriesFast(),fClustersESD->GetEntriesFast()));

  while(esdTOFHitArr->GetEntriesFast()){ // remove all hits
    delete esdTOFHitArr->RemoveAt(esdTOFHitArr->GetEntriesFast()-1);
  }
  for(Int_t i=0;i < fHitsESD->GetEntriesFast();i++){
    AliESDTOFHit *hitToStored = (AliESDTOFHit *) fHitsESD->At(i);
    AliESDTOFHit *hitNew = new ( (*esdTOFHitArr)[esdTOFHitArr->GetEntriesFast()] ) AliESDTOFHit(*hitToStored);
    hitNew->SetESDTOFClusterIndex(hitToStored->GetESDTOFClusterIndex());
  }

  AliInfo(Form("TOF after the matching: hits = %i and clusters = %i",esdTOFHitArr->GetEntriesFast(),esdTOFClArr->GetEntriesFast()));

  // Sort tof cluster
  for (Int_t i=0; i<fNseeds; i++) {
    AliESDtrack *t =(AliESDtrack*)fSeeds->At(i);
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0)continue;
    t->SortTOFcluster();
  }


  fSeeds->Clear();
  fClustersESD->Clear();
  fHitsESD->Clear();
  return 0;
  
}
//_________________________________________________________________________
void AliTOFtrackerV2::CollectESD() {
   //prepare the set of ESD tracks to be matched to clusters in TOF

  Int_t seedsTOF1=0;
  Int_t seedsTOF3=0;
  Int_t seedsTOF2=0;
 
  AliTOFtrack track;

  for (Int_t i=0; i<fNseeds; i++) {

    AliESDtrack *t =(AliESDtrack*)fSeeds->At(i);
    if ((t->GetStatus()&AliESDtrack::kTPCout)==0)continue;

    track = *t; // New
    Float_t x = (Float_t)track.GetX(); //New

    // TRD 'good' tracks
    if ( ( (t->GetStatus()&AliESDtrack::kTRDout)!=0 ) ) {

      AliDebug(1,Form(" Before propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track.GetIntegratedLength()));

      // TRD 'good' tracks, already propagated at 371 cm
      if ( x >= AliTOFGeometry::Rmin() ) {

	if ( track.PropagateToInnerTOF() ) {

	  AliDebug(1,Form(" TRD propagated track till rho = %fcm."
			  " And then the track has been propagated till rho = %fcm.",
			  x, (Float_t)track.GetX()));

	  track.SetSeedIndex(i);
	  t->UpdateTrackParams(&track,AliESDtrack::kTOFin);
	  fNseedsTOF++;
	  seedsTOF1++;

	  AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track.GetIntegratedLength()));
	}
      }
      else { // TRD 'good' tracks, propagated rho<371cm

	if  ( track.PropagateToInnerTOF() ) {

	  AliDebug(1,Form(" TRD propagated track till rho = %fcm."
			  " And then the track has been propagated till rho = %fcm.",
			  x, (Float_t)track.GetX()));

	  track.SetSeedIndex(i);
	  t->UpdateTrackParams(&track,AliESDtrack::kTOFin);
	  fNseedsTOF++;
	  seedsTOF3++;

	  AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track.GetIntegratedLength()));
	}
      }
    }

    else { // Propagate the rest of TPCbp

      AliDebug(1,Form(" Before propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track.GetIntegratedLength()));

      if ( track.PropagateToInnerTOF() ) { 

	AliDebug(1,Form(" TPC propagated track till rho = %fcm."
			" And then the track has been propagated till rho = %fcm.",
			x, (Float_t)track.GetX()));

      	track.SetSeedIndex(i);
	t->UpdateTrackParams(&track,AliESDtrack::kTOFin);
	fNseedsTOF++;
	seedsTOF2++;

	AliDebug(1,Form(" After propagation till inner TOF radius, ESDtrackLength=%f, TOFtrackLength=%f",t->GetIntegratedLength(),track.GetIntegratedLength()));
      }
    }
  }

  AliInfo(Form("Number of TOF seeds = %d (kTRDout371 = %d, kTRDoutLess371 = %d, !kTRDout = %d)",fNseedsTOF,seedsTOF1,seedsTOF3,seedsTOF2));

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
  TClonesArray* TOFClArr = fClustersESD; // use temporary array
  TClonesArray* esdTOFClArr = fEvent->GetESDTOFClusters();
  TClonesArray* esdTOFHitArr = fEvent->GetESDTOFHits();

  if(Int_t(detDepth/stepSize) > fNsteps){ // create array for each step
    // Get the number of propagation steps
    fNsteps =(Int_t)(detDepth/stepSize);

    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
      if(fTimesAr[isp]) delete[] fTimesAr[isp];
    }

    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
      fTimesAr[isp] = new Double_t[fNsteps];
    }

    for (Int_t ii=0; ii<4; ii++)
      if(fTrackPos[ii])
	delete [] fTrackPos[ii];
    
    for (Int_t ii=0; ii<4; ii++) fTrackPos[ii] = new Float_t[fNsteps];
  }



  AliDebug(1,Form(" Number of steps to be done %d",fNsteps));

  AliDebug(1,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

  // Some init
  const Int_t kNclusterMax = 1000; // related to fN value
  TGeoHMatrix global[kNclusterMax];
  Int_t clind[kNclusterMax];
  Bool_t isClusterMatchable[kNclusterMax]; // true if track and cluster were already matched (set to false below upto nc < kNclusterMax)

  AliTOFtrack trackTOFin;

  //The matching loop
  for (Int_t iseed=0; iseed<fSeeds->GetEntriesFast(); iseed++) {
    AliESDtrack *t =(AliESDtrack*)fSeeds->At(iseed); // ciao replace with loop on ESD + kTOFin
    if( (t->GetStatus()&AliESDtrack::kTOFin) == 0 ) continue;

    trackTOFin = *t;

    for (Int_t ii=0; ii<4; ii++)
      for (Int_t jj=0; jj<fNsteps; jj++) fTrackPos[ii][jj]=0.;

    for (Int_t ii=0; ii<kNclusterMax; ii++) clind[ii]=-1;
    for (Int_t ii=0; ii<kNclusterMax; ii++) global[ii] = 0x0;
    for (Int_t ii=0; ii<kNclusterMax; ii++) isClusterMatchable[ii] = kFALSE;	  	

    Double_t timesOr[AliPID::kSPECIESC]; t->GetIntegratedTimes(timesOr,AliPID::kSPECIESC); // in ps

    // Determine a window around the track
    Double_t x,par[5]; 
    trackTOFin.GetExternalParameters(x,par);
    Double_t cov[15]; 
    trackTOFin.GetExternalCovariance(cov);

    if (cov[0]<0. || cov[2]<0.) {
      AliWarning(Form("Very strange track (%d)! At least one of its covariance matrix diagonal elements is negative!",iseed));
      continue;
    }

    Double_t dphi=
      scaleFact*
      ((5*TMath::Sqrt(TMath::Abs(cov[0])) + 0.5*dY + 2.5*TMath::Abs(par[2]))/sensRadius); 
    Double_t dz=
       scaleFact*
       (5*TMath::Sqrt(TMath::Abs(cov[2])) + 0.5*dZ + 2.5*TMath::Abs(par[3]));

    Double_t phi=TMath::ATan2(par[0],x) + trackTOFin.GetAlpha();
    if (phi<-TMath::Pi())phi+=2*TMath::Pi();
    if (phi>=TMath::Pi())phi-=2*TMath::Pi();
    Double_t z=par[1];   

    //upper limit on window's size.
    if (dz> dzMax) dz=dzMax;
    if (dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;


    // find the clusters in the window of the track
    Int_t nc=0;
    for (Int_t k=FindClusterIndex(z-dz); k<TOFClArr->GetEntriesFast(); k++) {

      if (nc>=kNclusterMax) {
 	AliWarning("No more matchable clusters can be stored! Please, increase the corresponding vectors size.");
 	break;
      }

      AliESDTOFCluster *c=(AliESDTOFCluster *) TOFClArr->At(k);
      if (c->GetZ() > z+dz) break;
      if (!c->GetStatus()) {
	AliDebug(1,"Cluster in channel declared bad!");
	continue; // skip bad channels as declared in OCDB
      }


      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;

      Double_t yc=(c->GetPhi() - trackTOFin.GetAlpha())*c->GetR();
      Double_t p[2]={yc, c->GetZ()};
      Double_t cov2[3]= {dY*dY/12., 0., dZ*dZ/12.};
      if (trackTOFin.AliExternalTrackParam::GetPredictedChi2(p,cov2) > maxChi2)continue;

      clind[nc] = k;      
      Char_t path[200];
      Int_t ind[5]; AliTOFGeometry::GetVolumeIndices(c->GetTOFchannel(),ind);
      AliTOFGeometry::GetVolumePath(ind,path);
      gGeoManager->cd(path);
      global[nc] = *gGeoManager->GetCurrentMatrix();
      nc++;
    }


    if (nc == 0 ) {
      AliDebug(1,Form("No available clusters for the track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" Number of available TOF clusters for the track number %d: %d",iseed,nc));

    //start fine propagation 

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<fNsteps; istep++){ 
      
      // First of all, propagate the track...
      Float_t xs = AliTOFGeometry::RinTOF()+istep*stepSize;
      if (!(trackTOFin.PropagateTo(xs))) break;

      //  ...and then, if necessary, rotate the track
      Double_t ymax = xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
      Double_t ysect = trackTOFin.GetY();
      if (ysect > ymax) {
	if (!(trackTOFin.Rotate(AliTOFGeometry::GetAlpha()))) break;
      } else if (ysect <-ymax) {
	if (!(trackTOFin.Rotate(-AliTOFGeometry::GetAlpha()))) break;
      }

      Double_t mom = trackTOFin.P();

      if(istep == 0){
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = stepSize/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }
      else{
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = fTimesAr[isp][nStepsDone-1] + (trackTOFin.GetIntegratedLength()-fTrackPos[3][nStepsDone-1])/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }

      // store the running point (Globalrf) - fine propagation     

      Double_t r[3]; trackTOFin.GetXYZ(r);
      fTrackPos[0][nStepsDone]= (Float_t) r[0];
      fTrackPos[1][nStepsDone]= (Float_t) r[1];
      fTrackPos[2][nStepsDone]= (Float_t) r[2];   
      fTrackPos[3][nStepsDone]= trackTOFin.GetIntegratedLength();

      nStepsDone++;
      AliDebug(3,Form(" current step %d (%d) - nStepsDone=%d",istep,fNsteps,nStepsDone));
    }

    if ( nStepsDone == 0 ) {
      AliDebug(1,Form(" No track points for track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(3,Form(" Number of steps done for the track number %d: %d",iseed,nStepsDone));

    if(nc){
      for (Int_t i=0; i<nc; i++) isClusterMatchable[i] = kFALSE;	  	
    }

    Int_t nfound = 0;
    Bool_t accept = kFALSE;
    for (Int_t istep=0; istep<nStepsDone; istep++) {
      Float_t ctrackPos[3];     
      ctrackPos[0] = fTrackPos[0][istep];
      ctrackPos[1] = fTrackPos[1][istep];
      ctrackPos[2] = fTrackPos[2][istep];

      //now see whether the track matches any of the TOF clusters            

      Float_t dist3d[3]={0.,0.,0.};
      accept = kFALSE;

      for (Int_t i=0; i<nc; i++) {

	AliTOFGeometry::IsInsideThePad((TGeoHMatrix*)(&global[i]),ctrackPos,dist3d);

	// check multiple hit cases
	AliESDTOFCluster *cmatched=(AliESDTOFCluster *) TOFClArr->At(clind[i]);

	if(cmatched->GetNTOFhits() > 1){ // correct residual for mean position of the clusters (w.r.t. the first pad/hit)
	  Float_t zmain = cmatched->GetTOFchannel(0)/48;
	  Float_t xmain = cmatched->GetTOFchannel(0)%48;
	  for(Int_t ihit=1;ihit < cmatched->GetNTOFhits();ihit++){
	    Float_t deltaz = (cmatched->GetTOFchannel(ihit)/48 - zmain) * 3.5;
	    Float_t deltax = (cmatched->GetTOFchannel(ihit)%48 - xmain) * 2.5;
	    dist3d[0] -= deltax / cmatched->GetNTOFhits();
	    dist3d[2] -= deltaz / cmatched->GetNTOFhits();
	  }
	}

        // ***** NEW *****
        /* if track is inside this cluster set flags which will then
         * inhibit to add track points for the other clusters */
	Float_t yLoc = dist3d[1];
	Float_t rLoc = TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[2]*dist3d[2]);
	accept = (TMath::Abs(yLoc)<padDepth*0.5 && rLoc<dCut);

	//***** NEW *****
	/* add point everytime that:
	 * - the tracks is within dCut from the cluster
	 */
        if (accept) {

	  Double_t timesCurrent[AliPID::kSPECIESC];
	  AliDebug(3,Form(" Momentum for track %d -> %f", iseed,t->P()));
	  for (Int_t j=0;j<AliPID::kSPECIESC;j++) {
	    timesCurrent[j] = timesOr[j] + fTimesAr[j][istep];
	  }


	  if (TMath::Abs(dist3d[1])<stepSize && !isClusterMatchable[i]) {
	    isClusterMatchable[i] = kTRUE;
	    
	    Int_t currentpos = esdTOFClArr->GetEntriesFast(); // position of cluster in ESD
	    if(fWrittenInPos[clind[i]] != -1){
	      currentpos = fWrittenInPos[clind[i]];
	      cmatched = (AliESDTOFCluster *) esdTOFClArr->At(currentpos); // update the new one in the ESDEvent
	    }
	    else{ // add as a new cluster in the ESD TClonesArray
	      AliESDTOFCluster *clnew =  new( (*esdTOFClArr)[currentpos] ) AliESDTOFCluster(*cmatched);
	      clnew->SetEvent(fEvent);
	      clnew->SetESDID(currentpos);

	      // remap also TOF hit in the filtered array
	      for(Int_t ii=0;ii < cmatched->GetNTOFhits();ii++){
		Int_t index = cmatched->GetHitIndex(ii);
		AliESDTOFHit *hitOld = (AliESDTOFHit *) esdTOFHitArr->At(index);
		Int_t index_new = fHitsESD->GetEntriesFast();
		AliESDTOFHit *hitNew = new( (*fHitsESD)[index_new] ) AliESDTOFHit(*hitOld);
		hitNew->SetESDTOFClusterIndex(currentpos);
		clnew->SetHitIndex(ii,index_new);
	      }

	      fWrittenInPos[clind[i]] = currentpos;
	      cmatched = clnew; // update the new one added to the ESDEvent
	    }

	    if(cmatched->GetNMatchableTracks() < AliESDTOFCluster::kMaxMatches){
	      cmatched->Update(t->GetID(),dist3d[0],dist3d[1],dist3d[2],fTrackPos[3][istep],timesCurrent);//x,y,z -> tracking RF
	      t->AddTOFcluster(currentpos);
	      t->SetStatus(AliESDtrack::kTOFout);
	    }
	  }
          AliDebug(2,Form(" dist3dLoc[0] = %f, dist3dLoc[1] = %f, dist3dLoc[2] = %f ",dist3d[0],dist3d[1],dist3d[2]));

          nfound++;
        
          // ***** NEW *****
        }//end if accept
        
      } //end for on the clusters
    } //end for on the steps     


    if (nfound == 0 ) {
      AliDebug(1,Form(" No matchable track points for the track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" Number of track points for the track number %d: %d",iseed,nfound));

    Int_t nMatchedClusters = t->GetNTOFclusters();
 
    if (nMatchedClusters==0) {
      AliDebug(1,Form("Reconstructed track %d doesn't match any TOF cluster", iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" %d - matched (%d)",iseed,nMatchedClusters));

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
  } // loop on fSeeds
 
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
  Int_t ncl =clusters->GetEntriesFast();
  AliInfo(Form("Number of clusters: %d",ncl));

  fN = ncl; // set cluster counter

  for(Int_t i=0;i < ncl;i++) // reset position of clusters in ESD
    fWrittenInPos[i] = -1;

  if(ncl==0){
    return 0;
  }

  for (Int_t i=0; i<ncl; i++) {
    AliTOFcluster *c=(AliTOFcluster*)clusters->UncheckedAt(i);

    /*
    Int_t ind[5];
    ind[0]=c->GetDetInd(0);
    ind[1]=c->GetDetInd(1);
    ind[2]=c->GetDetInd(2);
    ind[3]=c->GetDetInd(3);
    ind[4]=c->GetDetInd(4);
    Int_t calindex = AliTOFGeometry::GetIndex(ind);
    Int_t tofLabels[3]={c->GetLabel(0),c->GetLabel(1),c->GetLabel(2)};
    */
      
    fClusters[i] = c;
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
  TClonesArray* TOFClArr = fClustersESD;; // use temporary array
  Int_t n = TOFClArr->GetEntriesFast();

  if (n==0) return 0;
  if (z <= ((AliESDTOFCluster *) TOFClArr->At(0))->GetZ()) return 0;
  if (z > ((AliESDTOFCluster *) TOFClArr->At(n-1))->GetZ()) return n;
  Int_t b=0, e=n-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > ((AliESDTOFCluster *) TOFClArr->At(m))->GetZ()) b=m+1;

    else e=m; 
  }
  return m;
} 
//_________________________________________________________________________
Bool_t AliTOFtrackerV2::GetTrackPoint(Int_t index, AliTrackPoint& p) const
{
  // Get track space point with index i
  // Coordinates are in the global system
  TClonesArray* esdTOFClArr = fEvent->GetESDTOFClusters();
  AliESDTOFCluster *cl = (AliESDTOFCluster *) esdTOFClArr->At(index);

  Float_t xyz[3];
  xyz[0] = cl->GetR()*TMath::Cos(cl->GetPhi());
  xyz[1] = cl->GetR()*TMath::Sin(cl->GetPhi());
  xyz[2] = cl->GetZ();
  Float_t phiangle = (Int_t(cl->GetPhi()*TMath::RadToDeg()/20.)+0.5)*20.*TMath::DegToRad();
  Float_t sinphi = TMath::Sin(phiangle), cosphi = TMath::Cos(phiangle);
  Int_t tofChannel=cl->GetTOFchannel();
  Int_t ind[5]; AliTOFGeometry::GetVolumeIndices(tofChannel,ind);
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

  // Load 1 hit in 1 cluster (ESD)
  TClonesArray* TOFClArr = fClustersESD;// use a temporary copy //fEvent->GetESDTOFClusters();
  TClonesArray* esdTOFHitArr = fEvent->GetESDTOFHits();

  if(TOFClArr->GetEntriesFast()) TOFClArr->Clear();
  if(esdTOFHitArr->GetEntriesFast()) esdTOFHitArr->Clear();

  AliESDTOFCluster *c1 = NULL;
  AliESDTOFCluster *c2 = NULL;

  for(Int_t i=0; i < fN;i++){
    AliTOFcluster *c = fClusters[i];
    Int_t ind[5];
    ind[0]=c->GetDetInd(0);
    ind[1]=c->GetDetInd(1);
    ind[2]=c->GetDetInd(2);
    ind[3]=c->GetDetInd(3);
    ind[4]=c->GetDetInd(4);
    Int_t calindex = AliTOFGeometry::GetIndex(ind);
    Int_t tofLabels[3]={c->GetLabel(0),c->GetLabel(1),c->GetLabel(2)};
    
    new ( (*esdTOFHitArr)[i] ) AliESDTOFHit( AliTOFGeometry::TdcBinWidth()*c->GetTDC(),
                              AliTOFGeometry::TdcBinWidth()*c->GetTDCRAW(),
                              AliTOFGeometry::ToTBinWidth()*c->GetToT()*1E-3,
                              calindex,tofLabels,c->GetL0L1Latency(),
                              c->GetDeltaBC(),i,c->GetZ(),c->GetR(),c->GetPhi() );
    
    c1 =  new( (*TOFClArr)[i] ) AliESDTOFCluster(i);
    c1->SetEvent(fEvent);
    c1->SetStatus( c->GetStatus() );
    c1->SetESDID(i);
    // 
    // register hits in the cluster
    c1->AddESDTOFHitIndex(i);

  }
  // start to merge clusters
  Int_t chan1,chan2,chan3;
  Int_t strip1,strip2;
  Int_t iphi,iphi2,iphi3;
  Int_t ieta,ieta2,ieta3;

  for(Int_t i=0; i <  TOFClArr->GetEntriesFast()-1;i++){
    c1=(AliESDTOFCluster *) TOFClArr->At(i);
    if(!c1->GetStatus()) continue;

    chan1 = c1->GetTOFchannel(0);
    AliTOFGeometry::GetVolumeIndices(chan1, detId); // Get volume index from channel index

    ieta = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
    if(detId[1]/*module*/ == 0) ieta += 0;
    else if(detId[1] == 1) ieta += 38;
    else if(detId[1] == 2) ieta += 76;
    else if(detId[1] == 3) ieta += 106;
    else if(detId[1] == 4) ieta += 144;
    iphi = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;

    chan2=chan1;
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

    // 1 and 2 belong now to the first cluster, 3 to the second one
    
    strip1 = chan1/96;
    for(Int_t j=i+1; j < TOFClArr->GetEntriesFast();j++){
      c2=(AliESDTOFCluster *) TOFClArr->At(j);
      if(!c2->GetStatus()) continue;

      chan3 = c2->GetTOFchannel();

      // check if the two TOF hits are in the same strip
      strip2 = chan3/96;
      if(strip1 != strip2) continue;

      AliTOFGeometry::GetVolumeIndices(chan3, detId); // Get volume index from channel index
      ieta3 = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
      if(detId[1]/*module*/ == 0) ieta3 += 0;
      else if(detId[1] == 1) ieta3 += 38;
      else if(detId[1] == 2) ieta3 += 76;
      else if(detId[1] == 3) ieta3 += 106;
      else if(detId[1] == 4) ieta3 += 144;
      iphi3 = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;
      

      if(ieta3-ieta > 2) j = fN; // because cluster are order along Z so you can skip all the rest, go to the next one ("i+1")

      // check if the fired pad are close in space
      if((TMath::Abs(iphi-iphi3)>1 && TMath::Abs(iphi2-iphi3)>1) || (TMath::Abs(ieta-ieta3)>1 && TMath::Abs(ieta2-ieta3)>1)) 
	continue; // double checks
      
      // check if the TOF time are close enough to be merged
      if(TMath::Abs(c1->GetTime() - c2->GetTime()) > 500/*in ps*/) continue;
      
      // merge them
      MergeClusters(i,j);

      // new hit is added as a second hit for the first cluster 
      iphi2 = iphi3;
      ieta2 = ieta3;
    }
  }  
}

void AliTOFtrackerV2::MergeClusters(Int_t i,Int_t j){
  TClonesArray* TOFClArr = fClustersESD;// use a temporary copy //fEvent->GetESDTOFClusters();

  if(i == j){
    AliInfo("No TOF cluster mergine possible (cannot merge a cluster with itself)");
    return;
  }

  if(i > j){ // check right order
    Int_t k=i;
    i=j;
    j=k;
  }

  Int_t last = TOFClArr->GetEntriesFast()-1;

  if(j > last){
    AliInfo("No TOF cluster mergine possible (cluster not available)");
    return;
  }
  
  AliESDTOFCluster *c1 = (AliESDTOFCluster *) TOFClArr->At(i);
  AliESDTOFCluster *c2 = (AliESDTOFCluster *) TOFClArr->At(j);

  if(c2->GetNMatchableTracks()){
    AliInfo("No TOF cluster mergine possible (cluster already matched)");
    return; // cannot merge a cluster already matched
  }

  Int_t nhit1 = c1->GetNTOFhits();
  Int_t nhit2 = c2->GetNTOFhits();

  if(nhit1+nhit2 >= AliESDTOFCluster::kMaxHits) 
    {
      AliInfo("No TOF cluster mergine possible (too many hits)");
      return;
    }

  for(Int_t k=0;k < nhit2;k++){// add hits in c2 to c1
    c1->AddESDTOFHitIndex(c2->GetHitIndex(k));

    // ID re-setting for hits not needed (done when matching is found)
  }

  // remove c2 from array
  if(j == last) delete TOFClArr->RemoveAt(j);
  else{
    for(Int_t ii=j;ii < last;ii++){
      AliESDTOFCluster *old= (AliESDTOFCluster *) TOFClArr->At(ii);
      if (!old) {AliFatal(Form("NULL pointer for TOF cluster %d",ii));}
      AliESDTOFCluster *replace= (AliESDTOFCluster *) TOFClArr->At(ii+1);
      if (!replace) {AliFatal(Form("NULL pointer for TOF cluster %d",ii+1));}
      *old = *replace;
      old->SetESDID(j);
    }
    delete TOFClArr->RemoveAt(last);
  }

}
