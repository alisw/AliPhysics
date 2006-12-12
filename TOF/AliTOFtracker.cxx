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
// AliTOFtracker Class                                                //
// Task: Perform association of the ESD tracks to TOF Clusters        //
// and Update ESD track with associated TOF Cluster parameters        //
//                                                                    //
// -- Authors : S. Arcelli, C. Zampolli (Bologna University and INFN) //
// -- Contacts: Annalisa.De.Caro@cern.ch                              //
// --         : Chiara.Zampolli@bo.infn.it                            //
// --         : Silvia.Arcelli@bo.infn.it                             //
//                                                                    //
//--------------------------------------------------------------------//

#include "Rtypes.h"

#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TTree.h"

#include "AliAlignObj.h"
#include "AliESDtrack.h"
#include "AliESD.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"

#include "AliTOFcalib.h"
#include "AliTOFcluster.h"
#include "AliTOFGeometry.h"
#include "AliTOFtracker.h"
#include "AliTOFtrack.h"

extern TGeoManager *gGeoManager;

ClassImp(AliTOFtracker)

//_____________________________________________________________________________
AliTOFtracker::AliTOFtracker(AliTOFGeometry * geom, Double_t parPID[2]):
  fGeom(geom),
  fTOFpid(new AliTOFpidESD(parPID)),
  fHoles(kFALSE),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fngoodmatch(0),
  fnbadmatch(0),
  fnunmatch(0),
  fnmatch(0),
  fR(378.), 
  fTOFHeigth(15.3),  
  fdCut(3.), 
  fDx(1.5), 
  fDy(0), 
  fDz(0), 
  fDzMax(35.), 
  fDyMax(50.), 
  fTracks(0x0),
  fSeeds(0x0)
 { 
  //AliTOFtracker main Ctor

  fDy=AliTOFGeometry::XPad(); 
  fDz=AliTOFGeometry::ZPad(); 
}
//_____________________________________________________________________________
AliTOFtracker::AliTOFtracker(const AliTOFtracker &t):
  AliTracker(),
  fGeom(0x0),
  fTOFpid(0x0),
  fHoles(kFALSE),
  fN(0),
  fNseeds(0),
  fNseedsTOF(0),
  fngoodmatch(0),
  fnbadmatch(0),
  fnunmatch(0),
  fnmatch(0),
  fR(378.), 
  fTOFHeigth(15.3),  
  fdCut(3.), 
  fDx(1.5), 
  fDy(0), 
  fDz(0), 
  fDzMax(35.), 
  fDyMax(50.), 
  fTracks(0x0),
  fSeeds(0x0)
 { 
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
  fDx=t.fDx; 
  fDzMax=t.fDzMax; 
  fDyMax=t.fDyMax; 
  fSeeds=t.fSeeds;
  fTracks=t.fTracks;
  fN=t.fN;
}

//_____________________________________________________________________________
AliTOFtracker& AliTOFtracker::operator=(const AliTOFtracker &t)
{ 
  //AliTOFtracker assignment operator

  this->fHoles=t.fHoles;
  this->fNseeds=t.fNseeds;
  this->fNseedsTOF=t.fNseedsTOF;
  this->fngoodmatch=t.fngoodmatch;
  this->fnbadmatch=t.fnbadmatch;
  this->fnunmatch=t.fnunmatch;
  this->fnmatch=t.fnmatch;
  this->fGeom = t.fGeom;
  this->fTOFpid = t.fTOFpid;
  this->fR=t.fR; 
  this->fTOFHeigth=t.fTOFHeigth;  
  this->fdCut=t.fdCut; 
  this->fDy=t.fDy; 
  this->fDz=t.fDz; 
  this->fDx=t.fDx; 
  this->fDzMax=t.fDzMax; 
  this->fDyMax=t.fDyMax; 
  this->fSeeds=t.fSeeds;
  this->fTracks=t.fTracks;
  this->fN=t.fN;
  return *this;

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
  fSeeds= new TClonesArray("AliESDtrack",ntrk);
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
  /*
  for (Int_t ijk=0; ijk<fN; ijk++) {
    AliInfo(Form("%4i %4i  %f %f %f  %f %f   %2i %1i %2i %1i %2i",ijk, fClusters[ijk]->GetIndex(),fClusters[ijk]->GetZ(),fClusters[ijk]->GetR(),fClusters[ijk]->GetPhi(), fClusters[ijk]->GetTDC(),fClusters[ijk]->GetADC(),fClusters[ijk]->GetDetInd(0),fClusters[ijk]->GetDetInd(1),fClusters[ijk]->GetDetInd(2),fClusters[ijk]->GetDetInd(3),fClusters[ijk]->GetDetInd(4)));
  }
  */

  //Second Step with Looser Matching Criterion
  MatchTracks(kTRUE);

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
      t->SetTOFCalChannel(seed->GetTOFCalChannel());
      Int_t tlab[3]; seed->GetTOFLabel(tlab);    
      t->SetTOFLabel(tlab);
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
	 ( x >= fGeom->RinTOF()) ){
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

  AliInfo(Form("Number of TOF seedds %i",fNseedsTOF));

  // Sort according uncertainties on track position 
  fTracks->Sort();

}
//_________________________________________________________________________
void AliTOFtracker::MatchTracks( Bool_t mLastStep){

  //Match ESD tracks to clusters in TOF


  Int_t nSteps=(Int_t)(fTOFHeigth/0.1);

  AliTOFcalib *calib = new AliTOFcalib(fGeom);
  //PH Arrays (moved outside of the loop)
  Float_t * trackPos[4];
  for (Int_t ii=0; ii<4; ii++) trackPos[ii] = new Float_t[nSteps];
  Int_t * clind[6];
  for (Int_t ii=0;ii<6;ii++) clind[ii] = new Int_t[fN];
  
  for (Int_t iseed=0; iseed<fNseedsTOF; iseed++) {

    AliTOFtrack *track =(AliTOFtrack*)fTracks->UncheckedAt(iseed);
    AliESDtrack *t =(AliESDtrack*)fSeeds->UncheckedAt(track->GetSeedIndex());
    if(t->GetTOFsignal()>0. ) continue;
    AliTOFtrack *trackTOFin =new AliTOFtrack(*track);

    // Some init

    Int_t         index[10000];
    Float_t        dist[10000];
    Float_t       cxpos[10000];
    Float_t       crecL[10000];
    TGeoHMatrix   global[1000];
     
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

    //upper limit on window's size.

    if(dz> fDzMax) dz=fDzMax;
    if(dphi*fR>fDyMax) dphi=fDyMax/fR;


    Int_t nc=0;

    // find the clusters in the window of the track

    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {

      AliTOFcluster *c=fClusters[k];
      if (c->GetZ() > z+dz) break;
      if (c->IsUsed()) continue;

      if (!c->GetStatus()) continue; // skip bad channels as declared in OCDB
      
      //AliInfo(Form(" fClusters[k]->GetZ() (%f) z-dz (%f)   %4i ", fClusters[k]->GetZ(), z-dz, k));

      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;

      {
      Double_t maxChi2=150.; // "calibration constant". Needs to be tuned.
      Double_t yc=(c->GetPhi() - trackTOFin->GetAlpha())*c->GetR();
      Double_t p[2]={yc, c->GetZ()};
      Double_t cov[3]={fDy*fDy/12., 0., fDz*fDz/12.};
      if (trackTOFin->AliExternalTrackParam::GetPredictedChi2(p,cov) > 150.) 
         continue;
      }

      clind[0][nc] = c->GetDetInd(0);
      clind[1][nc] = c->GetDetInd(1);
      clind[2][nc] = c->GetDetInd(2);
      clind[3][nc] = c->GetDetInd(3);
      clind[4][nc] = c->GetDetInd(4);
      clind[5][nc] = k;      
      Char_t path[100];
      Int_t ind[5];
      ind[0]=clind[0][nc];
      ind[1]=clind[1][nc];
      ind[2]=clind[2][nc];
      ind[3]=clind[3][nc];
      ind[4]=clind[4][nc];
      fGeom->GetVolumePath(ind,path);
      gGeoManager->cd(path);
      global[nc] = *gGeoManager->GetCurrentMatrix();
      nc++;
    }

    //if (nc) AliInfo(Form("seed for TOF %4i and number of clusters in the track window %4i (cluster index %4i)     %4i",i,nc, clind[5][0], fN));

    //start fine propagation 

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<nSteps; istep++){ 

      Float_t xs=fGeom->RinTOF()+istep*0.1;
      Double_t ymax=xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());

      Bool_t skip=kFALSE;
      Double_t ysect=trackTOFin->GetYat(xs,skip);
      if (skip) break;
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

      Double_t r[3];
      trackTOFin->GetXYZ(r);
      trackPos[0][istep]= (Float_t) r[0];
      trackPos[1][istep]= (Float_t) r[1];
      trackPos[2][istep]= (Float_t) r[2];   
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
	if( mLastStep)accept = (fGeom->DistanceToPad(cind,global[i],ctrackPos)<fdCut);
	if(!mLastStep)accept = (fGeom->IsInsideThePad(cind,global[i],ctrackPos));
	if(accept){
	  if(!mLastStep)isInside=kTRUE;
	  dist[nfound]=fGeom->DistanceToPad(cind,global[i],ctrackPos);
	  crecL[nfound]=trackPos[3][istep];
	  index[nfound]=clind[5][i]; // store cluster id 	    
	  cxpos[nfound]=fGeom->RinTOF()+istep*0.1; //store prop.radius
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
    c->Use(); //AliInfo(Form("I am using the cluster"));

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

      //AliInfo(Form(" track label good %5i",trackTOFin->GetLabel()));

    }
    else{
      fnbadmatch++;

      //AliInfo(Form(" track label  bad %5i",trackTOFin->GetLabel()));

    }

    delete trackTOFin;

    //  Store quantities to be used in the TOF Calibration
    Float_t tToT=c->GetToT(); // in ps
    t->SetTOFsignalToT(tToT);
    Int_t ind[5];
    ind[0]=c->GetDetInd(0);
    ind[1]=c->GetDetInd(1);
    ind[2]=c->GetDetInd(2);
    ind[3]=c->GetDetInd(3);
    ind[4]=c->GetDetInd(4);
    Int_t calindex = calib->GetIndex(ind);
    t->SetTOFCalChannel(calindex);

    // keep track of the track labels in the matched cluster
    Int_t tlab[3];
    tlab[0]=c->GetLabel(0);
    tlab[1]=c->GetLabel(1);
    tlab[2]=c->GetLabel(2);
    
    Double_t tof=AliTOFGeometry::TdcBinWidth()*c->GetTDC()+32; // in ps
    t->SetTOFsignal(tof);
    //t->SetTOFcluster(c->GetIndex()); // pointing to the digits tree
    t->SetTOFcluster(idclus); // pointing to the recPoints tree
    Double_t time[10]; t->GetIntegratedTimes(time);
    Double_t mom=t->GetP();
    for(Int_t j=0;j<=AliPID::kSPECIES;j++){
      Double_t mass=AliPID::ParticleMass(j);
      time[j]+=(recL-trackPos[3][0])/3e-2*TMath::Sqrt(mom*mom+mass*mass)/mom;
    }

    AliTOFtrack *trackTOFout = new AliTOFtrack(*t); 
    trackTOFout->PropagateTo(xpos);
    t->UpdateTrackParams(trackTOFout,AliESDtrack::kTOFout);    
    t->SetIntegratedLength(recL);
    t->SetIntegratedTimes(time);
    t->SetTOFLabel(tlab);

    delete trackTOFout;
  }
  for (Int_t ii=0; ii<4; ii++) delete [] trackPos[ii];
  for (Int_t ii=0;ii<6;ii++) delete [] clind[ii];
  delete calib;
}
//_________________________________________________________________________
Int_t AliTOFtracker::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------

  TBranch *branch=cTree->GetBranch("TOF");
  if (!branch) { 
    AliError("can't get the branch with the TOF clusters !");
    return 1;
  }

  TClonesArray dummy("AliTOFcluster",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  cTree->GetEvent(0);
  Int_t nc=clusters->GetEntriesFast();
  AliInfo(Form("Number of clusters: %d",nc));

  for (Int_t i=0; i<nc; i++) {
    AliTOFcluster *c=(AliTOFcluster*)clusters->UncheckedAt(i);
    fClusters[i]=new AliTOFcluster(*c); fN++;
    //AliInfo(Form("%4i %4i  %f %f %f  %f %f   %2i %1i %2i %1i %2i",i, fClusters[i]->GetIndex(),fClusters[i]->GetZ(),fClusters[i]->GetR(),fClusters[i]->GetPhi(), fClusters[i]->GetTDC(),fClusters[i]->GetADC(),fClusters[i]->GetDetInd(0),fClusters[i]->GetDetInd(1),fClusters[i]->GetDetInd(2),fClusters[i]->GetDetInd(3),fClusters[i]->GetDetInd(4)));
    //AliInfo(Form("%i %f",i, fClusters[i]->GetZ()));
  }

  //AliInfo(Form("Number of clusters: %d",fN));

  return 0;
}
//_________________________________________________________________________
void AliTOFtracker::UnloadClusters() {
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

//_________________________________________________________________________
Bool_t AliTOFtracker::GetTrackPoint(Int_t index, AliTrackPoint& p) const
{
  // Get track space point with index i
  // Coordinates are in the global system
  AliTOFcluster *cl = fClusters[index];
  Float_t xyz[3];
  xyz[0] = cl->GetR()*TMath::Cos(cl->GetPhi());
  xyz[1] = cl->GetR()*TMath::Sin(cl->GetPhi());
  xyz[2] = cl->GetZ();
  Float_t phiangle = (Int_t(cl->GetPhi()*TMath::RadToDeg()/20.)+0.5)*20.*TMath::DegToRad();
  Float_t sinphi = TMath::Sin(phiangle), cosphi = TMath::Cos(phiangle);
  Float_t tiltangle = fGeom->GetAngles(cl->GetDetInd(1),cl->GetDetInd(2))*TMath::DegToRad();
  Float_t sinth = TMath::Sin(tiltangle), costh = TMath::Cos(tiltangle);
  Float_t sigmay2 = fGeom->XPad()*fGeom->XPad()/12.;
  Float_t sigmaz2 = fGeom->ZPad()*fGeom->ZPad()/12.;
  Float_t cov[6];
  cov[0] = sinphi*sinphi*sigmay2 + cosphi*cosphi*sinth*sinth*sigmaz2;
  cov[1] = -sinphi*cosphi*sigmay2 + sinphi*cosphi*sinth*sinth*sigmaz2;
  cov[2] = -cosphi*sinth*costh*sigmaz2;
  cov[3] = cosphi*cosphi*sigmay2 + sinphi*sinphi*sinth*sinth*sigmaz2;
  cov[4] = -sinphi*sinth*costh*sigmaz2;
  cov[5] = costh*costh*sigmaz2;
  p.SetXYZ(xyz[0],xyz[1],xyz[2],cov);

  // Detector numbering scheme
  Int_t nSector = fGeom->NSectors();
  Int_t nPlate  = fGeom->NPlates();
  Int_t nStripA = fGeom->NStripA();
  Int_t nStripB = fGeom->NStripB();
  Int_t nStripC = fGeom->NStripC();

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
  UShort_t volid = AliAlignObj::LayerToVolUID(AliAlignObj::kTOF,idet);
  p.SetVolumeID((UShort_t)volid);
  return kTRUE;
}
