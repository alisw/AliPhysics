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
// AliTOFtracker Class
// Task: Perform association of the ESD tracks to TOF Clusters
// and Update ESD track with associated TOF Cluster parameters 
//
// -- Authors : S. Arcelli, C. Zampolli (Bologna University and INFN) 
// -- Contacts: Annalisa.De.Caro@cern.ch
// --         : Chiara.Zampolli@bo.infn.it
// --         : Silvia.Arcelli@bo.infn.it
//--------------------------------------------------------------------

#include <Rtypes.h>
#include "AliTOFtracker.h"
#include "AliTOFtrack.h"
#include "TClonesArray.h"
#include "TError.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"
#include "AliTOF.h"
#include "AliRun.h"
#include "AliModule.h"

ClassImp(AliTOFtracker)

//_____________________________________________________________________________
AliTOFtracker::AliTOFtracker(AliTOFGeometry * geom, Double_t parPID[2]) { 
  //AliTOFtracker main Ctor

  fHoles=true;
  fNseeds=0;
  fNseedsTOF=0;
  fngoodmatch=0;
  fnbadmatch=0;
  fnunmatch=0;
  fnmatch=0;
  fGeom = geom;
  fTOFpid = new AliTOFpidESD(parPID);
  fR=378.; 
  fTOFHeigth=15.3;  
  fdCut=3.; 
  fDy=AliTOFGeometry::XPad(); 
  fDz=AliTOFGeometry::ZPad(); 
  fDx=1.5; 
  fSeeds=0x0;
  fTracks=0x0;
  fN=0;
  Init(); // temporary solution to know about Holes/no Holes
}
//_____________________________________________________________________________
AliTOFtracker::AliTOFtracker(const AliTOFtracker &t):AliTracker() { 
  //AliTOFtracker copy Ctor

  fHoles=t.fHoles;
  fNseeds=t.fNseeds;
  fNseedsTOF=t.fNseedsTOF;
  fngoodmatch=t.fngoodmatch;
  fnbadmatch=t.fnbadmatch;
  fnunmatch=t.fnunmatch;
  fnmatch=t.fnmatch;
  fGeom = t.fGeom;
  fTOFpid = t.fTOFpid;
  fR=t.fR; 
  fTOFHeigth=t.fTOFHeigth;  
  fdCut=t.fdCut; 
  fDy=t.fDy; 
  fDz=t.fDz; 
  fDx=1.5; 
  fSeeds=t.fSeeds;
  fTracks=t.fTracks;
  fN=t.fN;
}
//_____________________________________________________________________________
void AliTOFtracker::Init() { 

// temporary solution to know about Holes/no Holes, will be implemented as 
// an AliTOFGeometry getter

  AliModule* frame=gAlice->GetModule("FRAME"); 

  if(!frame) {
    Error("Init","Could Not load FRAME! Assume Frame with Holes \n");
    fHoles=true;
  } else{
    if(frame->IsVersion()==1) {fHoles=false;}    
    else {fHoles=true;}      
  }
}
//_____________________________________________________________________________
Int_t AliTOFtracker::PropagateBack(AliESD* event) {
  //
  // Gets seeds from ESD event and Match with TOF Clusters
  //


  //Initialise some counters

  fNseeds=0;
  fNseedsTOF=0;
  fngoodmatch=0;
  fnbadmatch=0;
  fnunmatch=0;
  fnmatch=0;

  Int_t ntrk=event->GetNumberOfTracks();
  fNseeds = ntrk;
  fSeeds= new TClonesArray("AliESDtrack");
  TClonesArray &aESDTrack = *fSeeds;


  //Load ESD tracks into a local Array of ESD Seeds

  for (Int_t i=0; i<fNseeds; i++) {
    AliESDtrack *t=event->GetTrack(i);
    new(aESDTrack[i]) AliESDtrack(*t);
  }

  //Prepare ESD tracks candidates for TOF Matching
  CollectESD();

  //First Step with Strict Matching Criterion
  MatchTracks(kFALSE);

  //Second Step with Looser Matching Criterion
  MatchTracks(kTRUE);

  Info("PropagateBack","Number of matched tracks: %d",fnmatch);
  Info("PropagateBack","Number of good matched tracks: %d",fngoodmatch);
  Info("PropagateBack","Number of bad  matched tracks: %d",fnbadmatch);

  //Update the matched ESD tracks

  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    AliESDtrack *seed =(AliESDtrack*)fSeeds->UncheckedAt(i);
    if(seed->GetTOFsignal()>0){
      t->SetTOFsignal(seed->GetTOFsignal());
      t->SetTOFcluster(seed->GetTOFcluster());
      AliTOFtrack *track = new AliTOFtrack(*seed); 
      t->UpdateTrackParams(track,AliESDtrack::kTOFout);    
      delete track;
    }
  }


  //Make TOF PID
  fTOFpid->MakePID(event);

  if (fSeeds) {
    fSeeds->Delete();
    delete fSeeds;
    fSeeds = 0x0;
  }
  if (fTracks) {
    fTracks->Delete();
    delete fTracks;
    fTracks = 0x0;
  }
  return 0;
  
}
//_________________________________________________________________________
void AliTOFtracker::CollectESD() {
   //prepare the set of ESD tracks to be matched to clusters in TOF
 
  fTracks= new TClonesArray("AliTOFtrack");
  TClonesArray &aTOFTrack = *fTracks;
  for (Int_t i=0; i<fNseeds; i++) {

    AliESDtrack *t =(AliESDtrack*)fSeeds->UncheckedAt(i);
    if ((t->GetStatus()&AliESDtrack::kTPCout)==0)continue;

    // TRD good tracks, already propagated at 371 cm

    AliTOFtrack *track = new AliTOFtrack(*t); // New
    Double_t x = track->GetX(); //New

    if (((t->GetStatus()&AliESDtrack::kTRDout)!=0 ) && 
	 ( x >= AliTOFGeometry::RinTOF()) ){
      track->SetSeedIndex(i);
      t->UpdateTrackParams(track,AliESDtrack::kTOFout);    
      new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
      fNseedsTOF++;
      delete track;
    }

    // Propagate the rest of TPCbp  

    else {
      if(track->PropagateToInnerTOF(fHoles)){ // temporary solution
	//      if(track->PropagateToInnerTOF(fGeom->GetHoles())){
      	track->SetSeedIndex(i);
	t->UpdateTrackParams(track,AliESDtrack::kTOFout);    
 	new(aTOFTrack[fNseedsTOF]) AliTOFtrack(*track);
	fNseedsTOF++;
      }
      delete track;
    }
  }

  // Sort according uncertainties on track position 
  fTracks->Sort();

}
//_________________________________________________________________________
void AliTOFtracker::MatchTracks( Bool_t mLastStep){

  //Match ESD tracks to clusters in TOF

  static const Double_t kMasses[]={
    0.000511, 0.105658, 0.139570, 0.493677, 0.938272, 1.875613
  };
  
  Int_t nSteps=(Int_t)(fTOFHeigth/0.1);

  //PH Arrays (moved outside of the loop)
  Float_t * trackPos[4];
  for (Int_t ii=0; ii<4; ii++) trackPos[ii] = new Float_t[nSteps];
  Int_t * clind[6];
  for (Int_t ii=0;ii<6;ii++) clind[ii] = new Int_t[fN];
  
  for (Int_t i=0; i<fNseedsTOF; i++) {

    AliTOFtrack *track =(AliTOFtrack*)fTracks->UncheckedAt(i);
    AliESDtrack *t =(AliESDtrack*)fSeeds->UncheckedAt(track->GetSeedIndex());
    if(t->GetTOFsignal()>0. ) continue;
    AliTOFtrack *trackTOFin =new AliTOFtrack(*track);

    // Some init 

    Int_t         index[10000];
    Float_t        dist[10000];
    Float_t       cxpos[10000];
    Float_t       crecL[10000];
     
    // Determine a window around the track

    Double_t x,par[5]; 
    trackTOFin->GetExternalParameters(x,par);
    Double_t cov[15]; 
    trackTOFin->GetExternalCovariance(cov);
    Float_t scalefact=3.;    
    Double_t dphi=
      scalefact*
      ((5*TMath::Sqrt(cov[0]) + 0.5*fDy + 2.5*TMath::Abs(par[2]))/fR); 
    Double_t dz=
      scalefact*
      (5*TMath::Sqrt(cov[2]) + 0.5*fDz + 2.5*TMath::Abs(par[3]));
    
    Double_t phi=TMath::ATan2(par[0],x) + trackTOFin->GetAlpha();
    if (phi<-TMath::Pi())phi+=2*TMath::Pi();
    if (phi>=TMath::Pi())phi-=2*TMath::Pi();
    Double_t z=par[1];   

    Int_t nc=0;
    
    // find the clusters in the window of the track

    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {
      AliTOFcluster *c=fClusters[k];
      if (c->GetZ() > z+dz) break;
      if (c->IsUsed()) continue;
      
      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;
    
      clind[0][nc] = c->GetDetInd(0);
      clind[1][nc] = c->GetDetInd(1);
      clind[2][nc] = c->GetDetInd(2);
      clind[3][nc] = c->GetDetInd(3);
      clind[4][nc] = c->GetDetInd(4);
      clind[5][nc] = k;      
      nc++;
    }

    //start fine propagation 

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<nSteps; istep++){ 

      Float_t xs=AliTOFGeometry::RinTOF()+istep*0.1;
      Double_t ymax=xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());

      Bool_t skip=kFALSE;
      Double_t ysect=trackTOFin->GetYat(xs,skip);
      if(skip)break;
      if (ysect > ymax) {
	if (!trackTOFin->Rotate(AliTOFGeometry::GetAlpha())) {
	  break;
	}
      } else if (ysect <-ymax) {
	if (!trackTOFin->Rotate(-AliTOFGeometry::GetAlpha())) {
	  break;
	}
      }

      if(!trackTOFin->PropagateTo(xs)) {
	break;
      }

      nStepsDone++;

      // store the running point (Globalrf) - fine propagation     

      Double_t x,y,z;
      trackTOFin->GetGlobalXYZ(x,y,z);
      trackPos[0][istep]= (Float_t) x;
      trackPos[1][istep]= (Float_t) y;
      trackPos[2][istep]= (Float_t) z;   
      trackPos[3][istep]= trackTOFin->GetIntegratedLength();
    }


    Int_t nfound = 0;
    for (Int_t istep=0; istep<nStepsDone; istep++) {

      Bool_t isInside =kFALSE;
      Float_t ctrackPos[3];	

      ctrackPos[0]= trackPos[0][istep];
      ctrackPos[1]= trackPos[1][istep];
      ctrackPos[2]= trackPos[2][istep];

      //now see whether the track matches any of the TOF clusters            

      for (Int_t i=0; i<nc; i++){
	Int_t cind[5];
	cind[0]= clind[0][i];
	cind[1]= clind[1][i];
	cind[2]= clind[2][i];
	cind[3]= clind[3][i];
	cind[4]= clind[4][i];
        Bool_t accept = kFALSE;
        if( mLastStep)accept = (fGeom->DistanceToPad(cind,ctrackPos)<fdCut);
        if(!mLastStep)accept = (fGeom->IsInsideThePad(cind,ctrackPos));
	if(accept){
	  if(!mLastStep)isInside=kTRUE;
	  dist[nfound]=fGeom->DistanceToPad(cind,ctrackPos);
	  crecL[nfound]=trackPos[3][istep];
	  index[nfound]=clind[5][i]; // store cluster id 	    
	  cxpos[nfound]=AliTOFGeometry::RinTOF()+istep*0.1; //store prop.radius
	  nfound++;
	  if(isInside)break;
	}//end if accept
      } //end for on the clusters


      if(isInside)break;
    } //end for on the steps     



    if (nfound == 0 ) {
      fnunmatch++;
      delete trackTOFin;
      continue;
    }
    
    fnmatch++;

    // now choose the cluster to be matched with the track.

    Int_t idclus=0;
    Float_t  recL = 0.;
    Float_t  xpos=0.;
    Float_t  mindist=1000.;
    for (Int_t iclus= 0; iclus<nfound;iclus++){
      if (dist[iclus]< mindist){
      	mindist = dist[iclus];
      	xpos = cxpos[iclus];
        idclus =index[iclus]; 
        recL=crecL[iclus]+fDx*0.5;
      }
    }

    AliTOFcluster *c=fClusters[idclus];
    c->Use();

    // Track length correction for matching Step 2 

    if(mLastStep){
      Float_t rc=TMath::Sqrt(c->GetR()*c->GetR() + c->GetZ()*c->GetZ());
      Float_t rt=TMath::Sqrt(trackPos[0][70]*trackPos[0][70]
			     +trackPos[1][70]*trackPos[1][70]
			     +trackPos[2][70]*trackPos[2][70]);
      Float_t dlt=rc-rt;      
      recL=trackPos[3][70]+dlt;
    }    

    if (
	(c->GetLabel(0)==TMath::Abs(trackTOFin->GetLabel()))
	||
	(c->GetLabel(1)==TMath::Abs(trackTOFin->GetLabel()))
	||
	(c->GetLabel(2)==TMath::Abs(trackTOFin->GetLabel()))
	) {
      fngoodmatch++;
    }
    else{
      fnbadmatch++;
    }

    delete trackTOFin;

    Double_t tof=50*c->GetTDC()+32; // in ps
    t->SetTOFsignal(tof);
    t->SetTOFcluster(c->GetIndex());
    Double_t time[10]; t->GetIntegratedTimes(time);
    Double_t mom=t->GetP();
    for(Int_t j=0;j<=5;j++){
      Double_t mass=kMasses[j];
      time[j]+=(recL-trackPos[3][0])/3e-2*TMath::Sqrt(mom*mom+mass*mass)/mom;
    }

    AliTOFtrack *trackTOFout = new AliTOFtrack(*t); 
    trackTOFout->PropagateTo(xpos);
    t->UpdateTrackParams(trackTOFout,AliESDtrack::kTOFout);    
    t->SetIntegratedLength(recL);
    t->SetIntegratedTimes(time);

    delete trackTOFout;
  }
  for (Int_t ii=0; ii<4; ii++) delete [] trackPos[ii];
  for (Int_t ii=0;ii<6;ii++) delete [] clind[ii];
}
//_________________________________________________________________________
Int_t AliTOFtracker::LoadClusters(TTree *dTree) {
  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------

  TBranch *branch=dTree->GetBranch("TOF");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch with the TOF digits !\n");
    return 1;
  }

  TClonesArray dummy("AliTOFdigit",10000), *digits=&dummy;
  branch->SetAddress(&digits);

  dTree->GetEvent(0);
  Int_t nd=digits->GetEntriesFast();
  Info("LoadClusters","number of digits: %d",nd);

  for (Int_t i=0; i<nd; i++) {
    AliTOFdigit *d=(AliTOFdigit*)digits->UncheckedAt(i);
    Int_t dig[5]; Float_t g[3];
    dig[0]=d->GetSector();
    dig[1]=d->GetPlate();
    dig[2]=d->GetStrip();
    dig[3]=d->GetPadz();
    dig[4]=d->GetPadx();

    fGeom->GetPos(dig,g);

    Double_t h[5];
    h[0]=TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
    h[1]=TMath::ATan2(g[1],g[0]); h[2]=g[2]; 
    h[3]=d->GetTdc(); h[4]=d->GetAdc();

    AliTOFcluster *cl=new AliTOFcluster(h,d->GetTracks(),dig,i);
    InsertCluster(cl);
  }  

  return 0;
}
//_________________________________________________________________________
void AliTOFtracker::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads TOF clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  fN=0;
}

//_________________________________________________________________________
Int_t AliTOFtracker::InsertCluster(AliTOFcluster *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to the array of clusters sorted in Z
  //--------------------------------------------------------------------
  if (fN==kMaxCluster) {
    Error("InsertCluster","Too many clusters !\n");
    return 1;
  }

  if (fN==0) {fClusters[fN++]=c; return 0;}
  Int_t i=FindClusterIndex(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTOFcluster*));
  fClusters[i]=c; fN++;

  return 0;
}

//_________________________________________________________________________
Int_t AliTOFtracker::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
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

