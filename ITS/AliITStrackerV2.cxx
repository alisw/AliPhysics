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

//-------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    It reads AliITSRecPoint clusters and creates AliITStrackV2 tracks
//                   and fills with them the ESD
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------

#include <new>

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TGeoMatrix.h>

#include "AliITSgeom.h"
#include "AliAlignObj.h"
#include "AliITSRecPoint.h"
#include "AliESD.h"
#include "AliITSRecPoint.h"
#include "AliITStrackerV2.h"

ClassImp(AliITStrackerV2)

AliITStrackerV2::AliITSlayer AliITStrackerV2::fgLayers[kMaxLayer]; //ITS layers

AliITStrackerV2::AliITStrackerV2(): 
  AliTracker(), 
  fI(kMaxLayer),
  fBestTrack(),
  fTrackToFollow(),
  fPass(0),
  fLastLayerToTrackTo(kLastLayerToTrackTo)
{
  //--------------------------------------------------------------------
  //This is the AliITStrackerV2 default constructor
  //--------------------------------------------------------------------

  for (Int_t i=1; i<kMaxLayer+1; i++) new(fgLayers+i-1) AliITSlayer();

  fConstraint[0]=1; fConstraint[1]=0;

  Double_t xyz[]={kXV,kYV,kZV}, ers[]={kSigmaXV,kSigmaYV,kSigmaZV}; 
  SetVertex(xyz,ers);

  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=kLayersNotToSkip[i];

}

AliITStrackerV2::AliITStrackerV2(const AliITStrackerV2 &t): 
  AliTracker(t), 
  fI(t.fI),
  fBestTrack(t.fBestTrack),
  fTrackToFollow(t.fTrackToFollow),
  fPass(t.fPass),
  fLastLayerToTrackTo(t.fLastLayerToTrackTo)
{
  //--------------------------------------------------------------------
  //This is the AliITStrackerV2 copy constructor
  //--------------------------------------------------------------------

  //for (Int_t i=1; i<kMaxLayer+1; i++) new(fgLayers+i-1) AliITSlayer();

  fConstraint[0]=t.fConstraint[0]; fConstraint[1]=t.fConstraint[1];

  Double_t xyz[]={kXV,kYV,kZV}, ers[]={kSigmaXV,kSigmaYV,kSigmaZV};
  xyz[0]=t.GetX(); xyz[1]=t.GetY(); xyz[2]=t.GetZ(); 
  ers[0]=t.GetSigmaX(); ers[1]=t.GetSigmaY(); ers[2]=t.GetSigmaZ(); 
  SetVertex(xyz,ers);

  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=t.fLayersNotToSkip[i];

}

AliITStrackerV2::AliITStrackerV2(const AliITSgeom *geom) : 
  AliTracker(), 
  fI(kMaxLayer),
  fBestTrack(),
  fTrackToFollow(),
  fPass(0),
  fLastLayerToTrackTo(kLastLayerToTrackTo)
{
  //--------------------------------------------------------------------
  //This is the AliITStrackerV2 constructor
  //--------------------------------------------------------------------
  AliITSgeom *g=(AliITSgeom*)geom;

  Float_t x,y,z;  Int_t i;
  for (i=1; i<kMaxLayer+1; i++) {
    Int_t nlad=g->GetNladders(i);
    Int_t ndet=g->GetNdetectors(i);

    g->GetTrans(i,1,1,x,y,z); 
    Double_t poff=TMath::ATan2(y,x);
    Double_t zoff=z;
    Double_t r=TMath::Sqrt(x*x + y*y);

    g->GetTrans(i,1,2,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    g->GetTrans(i,2,1,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    g->GetTrans(i,2,2,x,y,z);
    r += TMath::Sqrt(x*x + y*y);
    r*=0.25;

    new (fgLayers+i-1) AliITSlayer(r,poff,zoff,nlad,ndet);

    for (Int_t j=1; j<nlad+1; j++) {
      for (Int_t k=1; k<ndet+1; k++) { //Fill this layer with detectors
        Float_t x,y,zshift; g->GetTrans(i,j,k,x,y,zshift); 
        Double_t rot[9]; g->GetRotMatrix(i,j,k,rot);

        Double_t phi=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
        phi+=TMath::Pi()/2;
        if (i==1) phi+=TMath::Pi();

        if (phi<0) phi+=TMath::TwoPi();
        else if (phi>=TMath::TwoPi()) phi-=TMath::TwoPi();

        Double_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);
        Double_t r=x*cp+y*sp;

        AliITSdetector &det=fgLayers[i-1].GetDetector((j-1)*ndet + k-1); 
        new(&det) AliITSdetector(r,phi); 
      } 
    }  

  }

  fConstraint[0]=1; fConstraint[1]=0;

  Double_t xyz[]={kXV,kYV,kZV}, ers[]={kSigmaXV,kSigmaYV,kSigmaZV}; 
  SetVertex(xyz,ers);

  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=kLayersNotToSkip[i];

}

void AliITStrackerV2::SetLayersNotToSkip(Int_t *l) {
  //--------------------------------------------------------------------
  //This function set masks of the layers which must be not skipped
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=l[i];
}

Int_t AliITStrackerV2::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads ITS clusters
  //--------------------------------------------------------------------
  TBranch *branch=cTree->GetBranch("ITSRecPoints");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch !\n");
    return 1;
  }

  TClonesArray dummy("AliITSRecPoint",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  Int_t j=0;
  for (Int_t i=0; i<kMaxLayer; i++) {
    Int_t ndet=fgLayers[i].GetNdetectors();
    Int_t jmax = j + fgLayers[i].GetNladders()*ndet;

    Double_t r=fgLayers[i].GetR();
    Double_t circ=TMath::TwoPi()*r;

    for (; j<jmax; j++) {           
      if (!cTree->GetEvent(j)) continue;
      Int_t ncl=clusters->GetEntriesFast();
      while (ncl--) {
        AliITSRecPoint *c=(AliITSRecPoint*)clusters->UncheckedAt(ncl);

        Int_t idx=c->GetDetectorIndex();
        AliITSdetector &det=fgLayers[i].GetDetector(idx);
   
        //Shift the cluster to the misaligned position
        Double_t x=det.GetR(); //y=...;      z=...;
        c->SetX(x);            //c->SetY(y); c->SetZ(z); 

        Double_t y=r*det.GetPhi()+c->GetY();
        if (y>circ) y-=circ; else if (y<0) y+=circ;
        c->SetPhiR(y);

        fgLayers[i].InsertCluster(new AliITSRecPoint(*c));
      }
      clusters->Delete();
    }
    fgLayers[i].ResetRoad(); //road defined by the cluster density
  }

  return 0;
}

void AliITStrackerV2::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer; i++) fgLayers[i].ResetClusters();
}

static Int_t CorrectForDeadZoneMaterial(AliITStrackV2 *t) {
  //--------------------------------------------------------------------
  // Correction for the material between the TPC and the ITS
  // (should it belong to the TPC code ?)
  //--------------------------------------------------------------------
  Double_t riw=80., diw=0.0053, x0iw=30; // TPC inner wall ? 
  Double_t rcd=61., dcd=0.0053, x0cd=30; // TPC "central drum" ?
  Double_t yr=12.8, dr=0.03; // rods ?
  Double_t zm=0.2, dm=0.40;  // membrane
  //Double_t rr=52., dr=0.19, x0r=24., yyr=7.77; //rails
  Double_t rs=50., ds=0.001; // something belonging to the ITS (screen ?)

  if (t->GetX() > riw) {
     if (!t->PropagateTo(riw,diw,x0iw)) return 1;
     if (TMath::Abs(t->GetY())>yr) t->CorrectForMaterial(dr);
     if (TMath::Abs(t->GetZ())<zm) t->CorrectForMaterial(dm);
     if (!t->PropagateTo(rcd,dcd,x0cd)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,dr,x0r); 
     if (!t->PropagateTo(rs,ds)) return 1;
  } else if (t->GetX() < rs) {
     if (!t->PropagateTo(rs,-ds)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,-dr,x0r); 
     if (!t->PropagateTo(rcd,-dcd,x0cd)) return 1;
     if (!t->PropagateTo(riw+0.001,-diw,x0iw)) return 1;
  } else {
  ::Error("CorrectForDeadZoneMaterial","track is already in the dead zone !");
    return 1;
  }
  
  return 0;
}

Int_t AliITStrackerV2::Clusters2Tracks(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions reconstructs ITS tracks
  // The clusters must be already loaded !
  //--------------------------------------------------------------------
  TObjArray itsTracks(15000);

  {/* Read ESD tracks */
    Int_t nentr=event->GetNumberOfTracks();
    Info("Clusters2Tracks", "Number of ESD tracks: %d\n", nentr);
    while (nentr--) {
      AliESDtrack *esd=event->GetTrack(nentr);

      if ((esd->GetStatus()&AliESDtrack::kTPCin)==0) continue;
      if (esd->GetStatus()&AliESDtrack::kTPCout) continue;
      if (esd->GetStatus()&AliESDtrack::kITSin) continue;

      AliITStrackV2 *t=0;
      try {
        t=new AliITStrackV2(*esd);
      } catch (const Char_t *msg) {
        Warning("Clusters2Tracks",msg);
        delete t;
        continue;
      }
      if (TMath::Abs(t->GetD(GetX(),GetY()))>4) {
	delete t;
	continue;
      }

      if (CorrectForDeadZoneMaterial(t)!=0) {
         Warning("Clusters2Tracks",
                 "failed to correct for the material in the dead zone !\n");
         delete t;
         continue;
      }
      itsTracks.AddLast(t);
    }
  } /* End Read ESD tracks */

  itsTracks.Sort();
  Int_t nentr=itsTracks.GetEntriesFast();

  Int_t ntrk=0;
  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
       AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
       if (t==0) continue;           //this track has been already tracked
       Int_t tpcLabel=t->GetLabel(); //save the TPC track label

       ResetTrackToFollow(*t);
       ResetBestTrack();

       for (FollowProlongation(); fI<kMaxLayer; fI++) {
          while (TakeNextProlongation()) FollowProlongation();
       }

       if (fBestTrack.GetNumberOfClusters() == 0) continue;

       if (fConstraint[fPass]) {
          ResetTrackToFollow(*t);
          if (!RefitAt(3.7, &fTrackToFollow, &fBestTrack)) continue;
          ResetBestTrack();
       }

       fBestTrack.SetLabel(tpcLabel);
       fBestTrack.CookdEdx();
       CookLabel(&fBestTrack,0.); //For comparison only
       fBestTrack.UpdateESDtrack(AliESDtrack::kITSin);
       UseClusters(&fBestTrack);
       delete itsTracks.RemoveAt(i);
       ntrk++;
     }
  }

  itsTracks.Delete();

  Info("Clusters2Tracks","Number of prolonged tracks: %d\n",ntrk);

  return 0;
}

Int_t AliITStrackerV2::PropagateBack(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions propagates reconstructed ITS tracks back
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  Int_t nentr=event->GetNumberOfTracks();
  Info("PropagateBack", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
     AliESDtrack *esd=event->GetTrack(i);

     if ((esd->GetStatus()&AliESDtrack::kITSin)==0) continue;
     if (esd->GetStatus()&AliESDtrack::kITSout) continue;

     AliITStrackV2 *t=0;
     try {
        t=new AliITStrackV2(*esd);
     } catch (const Char_t *msg) {
        Warning("PropagateBack",msg);
        delete t;
        continue;
     }

     ResetTrackToFollow(*t);

     // propagete to vertex [SR, GSI 17.02.2003]
     // Start Time measurement [SR, GSI 17.02.2003], corrected by I.Belikov
     if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {
       if (fTrackToFollow.PropagateToVertex(event->GetVertex())) {
          fTrackToFollow.StartTimeIntegral();
       }
       fTrackToFollow.PropagateTo(3.,-0.0028,65.19);
     }

     fTrackToFollow.ResetCovariance(10.); fTrackToFollow.ResetClusters();
     if (RefitAt(49.,&fTrackToFollow,t)) {
        if (CorrectForDeadZoneMaterial(&fTrackToFollow)!=0) {
          Warning("PropagateBack",
                  "failed to correct for the material in the dead zone !\n");
          delete t;
          continue;
        }
        fTrackToFollow.SetLabel(t->GetLabel());
        //fTrackToFollow.CookdEdx();
        CookLabel(&fTrackToFollow,0.); //For comparison only
        fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSout);
        UseClusters(&fTrackToFollow);
        ntrk++;
     }
     delete t;
  }

  Info("PropagateBack","Number of back propagated ITS tracks: %d\n",ntrk);

  return 0;
}

Int_t AliITStrackerV2::RefitInward(AliESD *event) {
  //--------------------------------------------------------------------
  // This functions refits ITS tracks using the 
  // "inward propagated" TPC tracks
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  Int_t nentr=event->GetNumberOfTracks();
  Info("RefitInward", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);

    if ((esd->GetStatus()&AliESDtrack::kITSout) == 0) continue;
    if (esd->GetStatus()&AliESDtrack::kITSrefit) continue;
    if (esd->GetStatus()&AliESDtrack::kTPCout)
    if ((esd->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;

    AliITStrackV2 *t=0;
    try {
        t=new AliITStrackV2(*esd);
    } catch (const Char_t *msg) {
        Warning("RefitInward",msg);
        delete t;
        continue;
    }

    if (CorrectForDeadZoneMaterial(t)!=0) {
       Warning("RefitInward",
               "failed to correct for the material in the dead zone !\n");
       delete t;
       continue;
    }

    ResetTrackToFollow(*t);
    fTrackToFollow.ResetClusters();

    //Refitting...
    if (RefitAt(3.7, &fTrackToFollow, t, kTRUE)) {
       fTrackToFollow.SetLabel(t->GetLabel());
       fTrackToFollow.CookdEdx();
       CookLabel(&fTrackToFollow,0.); //For comparison only

       if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {//The beam pipe 
         AliESDtrack  *esdTrack =fTrackToFollow.GetESDtrack();
         esdTrack->UpdateTrackParams(&fTrackToFollow,AliESDtrack::kITSrefit);
         Float_t r[3]={0.,0.,0.};
         Double_t maxD=3.;
	 esdTrack->RelateToVertex(event->GetVertex(),GetBz(r),maxD);
         ntrk++;
       }
    }
    delete t;
  }

  Info("RefitInward","Number of refitted tracks: %d\n",ntrk);

  return 0;
}

AliCluster *AliITStrackerV2::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fgLayers[l].GetCluster(c);
}


void AliITStrackerV2::FollowProlongation() {
  //--------------------------------------------------------------------
  //This function finds a track prolongation 
  //--------------------------------------------------------------------
  while (fI>fLastLayerToTrackTo) {
    Int_t i=fI-1;

    AliITSlayer &layer=fgLayers[i];
    AliITStrackV2 &track=fTracks[i];

    Double_t r=layer.GetR();

    if (i==3 || i==1) {
       Double_t rs=0.5*(fgLayers[i+1].GetR() + r);
       Double_t d=0.0034, x0=38.6;
       if (i==1) {rs=9.; d=0.0097; x0=42;}
       if (!fTrackToFollow.PropagateTo(rs,d,x0)) {
	 //Warning("FollowProlongation","propagation failed !\n");
         return;
       }
    }

    //find intersection
    Double_t x,y,z;  
    if (!fTrackToFollow.GetGlobalXYZat(r,x,y,z)) {
      //Warning("FollowProlongation","failed to estimate track !\n");
      return;
    }
    Double_t phi=TMath::ATan2(y,x);

    Int_t idet=layer.FindDetectorIndex(phi,z);
    if (idet<0) {
      //Warning("FollowProlongation","failed to find a detector !\n");
      return;
    }

    //propagate to the intersection
    const AliITSdetector &det=layer.GetDetector(idet);
    phi=det.GetPhi();
    if (!fTrackToFollow.Propagate(phi,det.GetR())) {
      //Warning("FollowProlongation","propagation failed !\n");
      return;
    }
    fTrackToFollow.SetDetectorIndex(idet);

    //Select possible prolongations and store the current track estimation
    track.~AliITStrackV2(); new(&track) AliITStrackV2(fTrackToFollow);
    Double_t dz=7*TMath::Sqrt(track.GetSigmaZ2() + kSigmaZ2[i]);
    Double_t dy=7*TMath::Sqrt(track.GetSigmaY2() + kSigmaY2[i]);
    Double_t road=layer.GetRoad();
    if (dz*dy>road*road) {
       Double_t dd=TMath::Sqrt(dz*dy), scz=dz/dd, scy=dy/dd;
       dz=road*scz; dy=road*scy;
    } 

    //Double_t dz=4*TMath::Sqrt(track.GetSigmaZ2() + kSigmaZ2[i]);
    if (dz < 0.5*TMath::Abs(track.GetTgl())) dz=0.5*TMath::Abs(track.GetTgl());
    if (dz > kMaxRoad) {
      //Warning("FollowProlongation","too broad road in Z !\n");
      return;
    }

    if (TMath::Abs(fTrackToFollow.GetZ()-GetZ()) > r+dz) return;

    //Double_t dy=4*TMath::Sqrt(track.GetSigmaY2() + kSigmaY2[i]);
    if (dy < 0.5*TMath::Abs(track.GetSnp())) dy=0.5*TMath::Abs(track.GetSnp());
    if (dy > kMaxRoad) {
      //Warning("FollowProlongation","too broad road in Y !\n");
      return;
    }

    fI--;

    Double_t zmin=track.GetZ() - dz; 
    Double_t zmax=track.GetZ() + dz;
    Double_t ymin=track.GetY() + r*phi - dy;
    Double_t ymax=track.GetY() + r*phi + dy;
    if (layer.SelectClusters(zmin,zmax,ymin,ymax)==0) 
       if (fLayersNotToSkip[fI]) return;  

    if (!TakeNextProlongation()) 
       if (fLayersNotToSkip[fI]) return;

  } 

  //deal with the best track
  Int_t ncl=fTrackToFollow.GetNumberOfClusters();
  Int_t nclb=fBestTrack.GetNumberOfClusters();
  if (ncl)
  if (ncl >= nclb) {
     Double_t chi2=fTrackToFollow.GetChi2();
     if (chi2/ncl < kChi2PerCluster) {        
        if (ncl > nclb || chi2 < fBestTrack.GetChi2()) {
           ResetBestTrack();
        }
     }
  }

}

Int_t AliITStrackerV2::TakeNextProlongation() {
  //--------------------------------------------------------------------
  // This function takes another track prolongation 
  //
  //  dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch 
  //--------------------------------------------------------------------
  AliITSlayer &layer=fgLayers[fI];
  ResetTrackToFollow(fTracks[fI]);

  Double_t dz=7*TMath::Sqrt(fTrackToFollow.GetSigmaZ2() + kSigmaZ2[fI]);
  Double_t dy=7*TMath::Sqrt(fTrackToFollow.GetSigmaY2() + kSigmaY2[fI]);
  Double_t road=layer.GetRoad();
  if (dz*dy>road*road) {
     Double_t dd=TMath::Sqrt(dz*dy), scz=dz/dd, scy=dy/dd;
     dz=road*scz; dy=road*scy;
  } 

  const AliITSRecPoint *c=0; Int_t ci=-1;
  const AliITSRecPoint *cc=0; Int_t cci=-1;
  Double_t chi2=kMaxChi2;
  while ((c=layer.GetNextCluster(ci))!=0) {
    Int_t idet=c->GetDetectorIndex();

    if (fTrackToFollow.GetDetectorIndex()!=idet) {
       const AliITSdetector &det=layer.GetDetector(idet);
       ResetTrackToFollow(fTracks[fI]);
       if (!fTrackToFollow.Propagate(det.GetPhi(),det.GetR())) {
         //Warning("TakeNextProlongation","propagation failed !\n");
         continue;
       }
       fTrackToFollow.SetDetectorIndex(idet);
       if (TMath::Abs(fTrackToFollow.GetZ()-GetZ())>layer.GetR()+dz) continue;
    }

    if (TMath::Abs(fTrackToFollow.GetZ() - c->GetZ()) > dz) continue;
    if (TMath::Abs(fTrackToFollow.GetY() - c->GetY()) > dy) continue;

    Double_t ch2=fTrackToFollow.GetPredictedChi2(c); 
    if (ch2 > chi2) continue;
    chi2=ch2;
    cc=c; cci=ci;
    break;
  }

  if (!cc) return 0;

  if (!fTrackToFollow.PropagateTo(cc->GetX(),0.,0.)) return 0; // Alignment
  if (!fTrackToFollow.Update(cc,chi2,(fI<<28)+cci)) {
     //Warning("TakeNextProlongation","filtering failed !\n");
     return 0;
  }

  if (fTrackToFollow.GetNumberOfClusters()>1)
    if (TMath::Abs(fTrackToFollow.GetD(GetX(),GetY()))>4) return 0;

  fTrackToFollow.
    SetSampledEdx(cc->GetQ(),fTrackToFollow.GetNumberOfClusters()-1); //b.b.

  {
  Double_t x0;
 Double_t d=layer.GetThickness(fTrackToFollow.GetY(),fTrackToFollow.GetZ(),x0);
  fTrackToFollow.CorrectForMaterial(d,x0);
  }

  if (fConstraint[fPass]) {
    Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
    Double_t xyz[]={GetX(),GetY(),GetZ()};
    Double_t ers[]={GetSigmaX(),GetSigmaY(),GetSigmaZ()};
    fTrackToFollow.Improve(d,xyz,ers);
  }

  return 1;
}


AliITStrackerV2::AliITSlayer::AliITSlayer():
  fR(0.),
  fPhiOffset(0.),
  fNladders(0),
  fZOffset(0.),
  fNdetectors(0),
  fDetectors(0),
  fNsel(0),
  fRoad(2*fR*TMath::Sqrt(3.14/1.)) //assuming that there's only one cluster
{
  //--------------------------------------------------------------------
  //default AliITSlayer constructor
  //--------------------------------------------------------------------
  
  for (Int_t i=0; i<kNsector; i++) fN[i]=0;

}

AliITStrackerV2::AliITSlayer::
AliITSlayer(Double_t r,Double_t p,Double_t z,Int_t nl,Int_t nd): 
  fR(r), 
  fPhiOffset(p), 
  fNladders(nl),
  fZOffset(z),
  fNdetectors(nd),
  fDetectors(new AliITSdetector[nl*nd]),
  fNsel(0),
  fRoad(2*r*TMath::Sqrt(3.14/1.)) //assuming that there's only one cluster
{
  //--------------------------------------------------------------------
  //main AliITSlayer constructor
  //--------------------------------------------------------------------

  for (Int_t i=0; i<kNsector; i++) fN[i]=0;

  for (Int_t i=0; i<kMaxClusterPerLayer; i++) fClusters[i]=0;

}

AliITStrackerV2::AliITSlayer::~AliITSlayer() {
  //--------------------------------------------------------------------
  // AliITSlayer destructor
  //--------------------------------------------------------------------
  delete[] fDetectors;
  ResetClusters();
}

void AliITStrackerV2::AliITSlayer::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------
   for (Int_t s=0; s<kNsector; s++) {
       Int_t &n=fN[s];
       while (n) {
          n--;
          delete fClusters[s*kMaxClusterPerSector+n];
       }
   }
}

void AliITStrackerV2::AliITSlayer::ResetRoad() {
  //--------------------------------------------------------------------
  // This function calculates the road defined by the cluster density
  //--------------------------------------------------------------------
  Int_t n=0;
  for (Int_t s=0; s<kNsector; s++) {
    Int_t i=fN[s];
    while (i--) 
       if (TMath::Abs(fClusters[s*kMaxClusterPerSector+i]->GetZ())<fR) n++;
  }
  if (n>1) fRoad=2*fR*TMath::Sqrt(3.14/n);
}

Int_t AliITStrackerV2::AliITSlayer::InsertCluster(AliITSRecPoint *c) {
  //--------------------------------------------------------------------
  // This function inserts a cluster to this layer in increasing
  // order of the cluster's fZ
  //--------------------------------------------------------------------
  Float_t circ=TMath::TwoPi()*fR;
  Int_t sec=Int_t(kNsector*c->GetPhiR()/circ);
  if (sec>=kNsector) {
     ::Error("InsertCluster","Wrong sector !\n");
     return 1;
  }
  Int_t &n=fN[sec];
  if (n>=kMaxClusterPerSector) {
     ::Error("InsertCluster","Too many clusters !\n");
     return 1;
  }
  if (n==0) fClusters[sec*kMaxClusterPerSector]=c;
  else {
     Int_t i=FindClusterIndex(c->GetZ(),sec);
     Int_t k=n-i+sec*kMaxClusterPerSector;
     memmove(fClusters+i+1 ,fClusters+i,k*sizeof(AliITSRecPoint*));
     fClusters[i]=c;
  }
  n++;
  return 0;
}

Int_t 
AliITStrackerV2::AliITSlayer::FindClusterIndex(Float_t z,Int_t s) const {
  //--------------------------------------------------------------------
  // For the sector "s", this function returns the index of the first 
  // with its fZ >= "z". 
  //--------------------------------------------------------------------
  Int_t nc=fN[s];
  if (nc==0) return kMaxClusterPerSector*s;

  Int_t b=kMaxClusterPerSector*s;
  if (z <= fClusters[b]->GetZ()) return b;

  Int_t e=b+nc-1;
  if (z > fClusters[e]->GetZ()) return e+1;

  Int_t m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}

Int_t AliITStrackerV2::AliITSlayer::
SelectClusters(Float_t zmin,Float_t zmax,Float_t ymin, Float_t ymax) {
  //--------------------------------------------------------------------
  // This function selects clusters within the "window"
  //--------------------------------------------------------------------
    Float_t circ=fR*TMath::TwoPi();

    if (ymin>circ) ymin-=circ; else if (ymin<0) ymin+=circ;
    if (ymax>circ) ymax-=circ; else if (ymax<0) ymax+=circ;

    Int_t i1=Int_t(kNsector*ymin/circ); if (i1==kNsector) i1--;
    if (fN[i1]!=0) {
       Float_t ym = (ymax<ymin) ? ymax+circ : ymax;
       Int_t i=FindClusterIndex(zmin,i1), imax=i1*kMaxClusterPerSector+fN[i1];
       for (; i<imax; i++) {
           AliITSRecPoint *c=fClusters[i];
           if (c->IsUsed()) continue;
           if (c->GetZ()>zmax) break;
           if (c->GetPhiR()<=ymin) continue;
           if (c->GetPhiR()>ym) continue;
           fIndex[fNsel++]=i;
       }
    }

    Int_t i2=Int_t(kNsector*ymax/circ); if (i2==kNsector) i2--;
    if (i2==i1) return fNsel;

    if (fN[i2]!=0) {
       Float_t ym = (ymin>ymax) ? ymin-circ : ymin;
       Int_t i=FindClusterIndex(zmin,i2), imax=i2*kMaxClusterPerSector+fN[i2];
       for (; i<imax; i++) {
           AliITSRecPoint *c=fClusters[i];
           if (c->IsUsed()) continue;
           if (c->GetZ()>zmax) break;
           if (c->GetPhiR()<=ym) continue;
           if (c->GetPhiR()>ymax) continue;
           fIndex[fNsel++]=i;
       }
    }

    return fNsel;
}

const AliITSRecPoint *AliITStrackerV2::AliITSlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------
  AliITSRecPoint *c=0;
  ci=-1;
  if (fNsel) {
     fNsel--;
     ci=fIndex[fNsel]; 
     c=fClusters[ci];
  }
  return c; 
}

Int_t AliITStrackerV2::AliITSlayer::GetNumberOfClusters() const {
  Int_t n=0;
  for (Int_t s=0; s<kNsector; s++) n+=fN[s];
  return n; 
}

Int_t 
AliITStrackerV2::AliITSlayer::FindDetectorIndex(Double_t phi,Double_t z)const {
  //--------------------------------------------------------------------
  //This function finds the detector crossed by the track
  //--------------------------------------------------------------------
  Double_t dphi=-(phi-fPhiOffset);
  if      (dphi <  0) dphi += 2*TMath::Pi();
  else if (dphi >= 2*TMath::Pi()) dphi -= 2*TMath::Pi();
  Int_t np=Int_t(dphi*fNladders*0.5/TMath::Pi()+0.5);
  if (np>=fNladders) np-=fNladders;
  if (np<0)          np+=fNladders;

  Double_t dz=fZOffset-z;
  Int_t nz=Int_t(dz*(fNdetectors-1)*0.5/fZOffset+0.5);
  if (nz>=fNdetectors) return -1;
  if (nz<0)            return -1;

  return np*fNdetectors + nz;
}

Double_t 
AliITStrackerV2::AliITSlayer::GetThickness(Double_t y,Double_t z,Double_t &x0)
const {
  //--------------------------------------------------------------------
  //This function returns the layer thickness at this point (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0085;
  x0=21.82;

  if (43<fR&&fR<45) { //SSD2
     Double_t dd=0.0034;
     d=dd;
     if (TMath::Abs(y-0.00)>3.40) d+=dd;
     if (TMath::Abs(y-1.90)<0.45) {d+=(0.013-0.0034);}
     if (TMath::Abs(y+1.90)<0.45) {d+=(0.013-0.0034);}
     for (Int_t i=0; i<12; i++) {
       if (TMath::Abs(z-3.9*(i+0.5))<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=0.0034; 
          break;
       }
       if (TMath::Abs(z+3.9*(i+0.5))<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=0.0034; 
          break;
       }         
       if (TMath::Abs(z-3.4-3.9*i)<0.50) {d+=(0.016-0.0034); break;}
       if (TMath::Abs(z+0.5+3.9*i)<0.50) {d+=(0.016-0.0034); break;}
     }
  } else 
  if (37<fR&&fR<41) { //SSD1
     Double_t dd=0.0034;
     d=dd;
     if (TMath::Abs(y-0.00)>3.40) d+=dd;
     if (TMath::Abs(y-1.90)<0.45) {d+=(0.013-0.0034);}
     if (TMath::Abs(y+1.90)<0.45) {d+=(0.013-0.0034);}
     for (Int_t i=0; i<11; i++) {
       if (TMath::Abs(z-3.9*i)<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=dd; 
          break;
       }
       if (TMath::Abs(z+3.9*i)<0.15) {
          if (TMath::Abs(y-0.00)>3.40) d+=dd;
          d+=dd; 
          break;
       }         
       if (TMath::Abs(z-1.85-3.9*i)<0.50) {d+=(0.016-0.0034); break;}
       if (TMath::Abs(z+2.05+3.9*i)<0.50) {d+=(0.016-0.0034); break;}         
     }
  } else
  if (13<fR&&fR<26) { //SDD
     Double_t dd=0.0033;
     d=dd;
     if (TMath::Abs(y-0.00)>3.30) d+=dd;

     if (TMath::Abs(y-1.80)<0.55) {
        d+=0.016;
        for (Int_t j=0; j<20; j++) {
          if (TMath::Abs(z+0.7+1.47*j)<0.12) {d+=0.08; x0=9.; break;}
          if (TMath::Abs(z-0.7-1.47*j)<0.12) {d+=0.08; x0=9.; break;}
        } 
     }
     if (TMath::Abs(y+1.80)<0.55) {
        d+=0.016;
        for (Int_t j=0; j<20; j++) {
          if (TMath::Abs(z-0.7-1.47*j)<0.12) {d+=0.08; x0=9.; break;}
          if (TMath::Abs(z+0.7+1.47*j)<0.12) {d+=0.08; x0=9.; break;}
        } 
     }

     for (Int_t i=0; i<4; i++) {
       if (TMath::Abs(z-7.3*i)<0.60) {
          d+=dd;
          if (TMath::Abs(y-0.00)>3.30) d+=dd; 
          break;
       }
       if (TMath::Abs(z+7.3*i)<0.60) {
          d+=dd; 
          if (TMath::Abs(y-0.00)>3.30) d+=dd; 
          break;
       }
     }
  } else
  if (6<fR&&fR<8) {   //SPD2
     Double_t dd=0.0063; x0=21.5;
     d=dd;
     if (TMath::Abs(y-3.08)>0.5) d+=dd;
     //if (TMath::Abs(y-3.08)>0.45) d+=dd;
     if (TMath::Abs(y-3.03)<0.10) {d+=0.014;}
  } else
  if (3<fR&&fR<5) {   //SPD1
     Double_t dd=0.0063; x0=21.5;
     d=dd;
     if (TMath::Abs(y+0.21)>0.6) d+=dd;
     //if (TMath::Abs(y+0.21)>0.45) d+=dd;
     if (TMath::Abs(y+0.10)<0.10) {d+=0.014;}
  }

  return d;
}

Double_t AliITStrackerV2::GetEffectiveThickness(Double_t y,Double_t z) const
{
  //--------------------------------------------------------------------
  //Returns the thickness between the current layer and the vertex (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0028*3*3; //beam pipe
  Double_t x0=0;

  Double_t xn=fgLayers[fI].GetR();
  for (Int_t i=0; i<fI; i++) {
    Double_t xi=fgLayers[i].GetR();
    d+=fgLayers[i].GetThickness(y,z,x0)*xi*xi;
  }

  if (fI>1) {
    Double_t xi=9.;
    d+=0.0097*xi*xi;
  }

  if (fI>3) {
    Double_t xi=0.5*(fgLayers[3].GetR()+fgLayers[4].GetR());
    d+=0.0034*xi*xi;
  }

  return d/(xn*xn);
}

Bool_t AliITStrackerV2::RefitAt(Double_t xx,AliITStrackV2 *t,
                                const AliITStrackV2 *c, Bool_t extra) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"
  // If "extra"==kTRUE, 
  //    the clusters from overlapped modules get attached to "t" 
  //--------------------------------------------------------------------
  Int_t index[kMaxLayer];
  Int_t k;
  for (k=0; k<kMaxLayer; k++) index[k]=-1;
  Int_t nc=c->GetNumberOfClusters();
  for (k=0; k<nc; k++) { 
    Int_t idx=c->GetClusterIndex(k),nl=(idx&0xf0000000)>>28;
    index[nl]=idx; 
  }

  Int_t from, to, step;
  if (xx > t->GetX()) {
      from=0; to=kMaxLayer;
      step=+1;
  } else {
      from=kMaxLayer-1; to=-1;
      step=-1;
  }

  for (Int_t i=from; i != to; i += step) {
     AliITSlayer &layer=fgLayers[i];
     Double_t r=layer.GetR();
 
     {
     Double_t hI=i-0.5*step; 
     if (TMath::Abs(hI-1.5)<0.01 || TMath::Abs(hI-3.5)<0.01) {             
        Double_t rs=0.5*(fgLayers[i-step].GetR() + r);
        Double_t d=0.0034, x0=38.6; 
        if (TMath::Abs(hI-1.5)<0.01) {rs=9.; d=0.0097; x0=42;}
        if (!t->PropagateTo(rs,-step*d,x0)) {
          return kFALSE;
        }
     }
     }

     // remember old position [SR, GSI 18.02.2003]
     Double_t oldX=0., oldY=0., oldZ=0.;
     if (t->IsStartedTimeIntegral() && step==1) {
        t->GetGlobalXYZat(t->GetX(),oldX,oldY,oldZ);
     }
     //

     Double_t x,y,z;
     if (!t->GetGlobalXYZat(r,x,y,z)) { 
       return kFALSE;
     }
     Double_t phi=TMath::ATan2(y,x);
     Int_t idet=layer.FindDetectorIndex(phi,z);
     if (idet<0) { 
       return kFALSE;
     }
     const AliITSdetector &det=layer.GetDetector(idet);
     phi=det.GetPhi();
     if (!t->Propagate(phi,det.GetR())) {
       return kFALSE;
     }
     t->SetDetectorIndex(idet);

     const AliITSRecPoint *cl=0;
     Double_t maxchi2=kMaxChi2;

     Int_t idx=index[i];
     if (idx>0) {
        const AliITSRecPoint *c=(AliITSRecPoint *)GetCluster(idx); 
        if (idet != c->GetDetectorIndex()) {
           idet=c->GetDetectorIndex();
           const AliITSdetector &det=layer.GetDetector(idet);
           if (!t->Propagate(det.GetPhi(),det.GetR())) {
             return kFALSE;
           }
           t->SetDetectorIndex(idet);
        }
        Double_t chi2=t->GetPredictedChi2(c);
        if (chi2<maxchi2) { 
	  cl=c; 
	  maxchi2=chi2; 
	} else {
	  return kFALSE;
	}
     }
 
     if (cl) {
       if (!t->PropagateTo(cl->GetX(),0.,0.)) return kFALSE; //Alignment
       if (!t->Update(cl,maxchi2,idx)) {
          return kFALSE;
       }
       t->SetSampledEdx(cl->GetQ(),t->GetNumberOfClusters()-1);
     }

     {
     Double_t x0;
     Double_t d=layer.GetThickness(t->GetY(),t->GetZ(),x0);
     t->CorrectForMaterial(-step*d,x0);
     }
                 
     if (extra) { //search for extra clusters
        AliITStrackV2 tmp(*t);
        Double_t dz=4*TMath::Sqrt(tmp.GetSigmaZ2()+kSigmaZ2[i]);
        if (dz < 0.5*TMath::Abs(tmp.GetTgl())) dz=0.5*TMath::Abs(tmp.GetTgl());
        Double_t dy=4*TMath::Sqrt(t->GetSigmaY2()+kSigmaY2[i]);
        if (dy < 0.5*TMath::Abs(tmp.GetSnp())) dy=0.5*TMath::Abs(tmp.GetSnp());
        Double_t zmin=t->GetZ() - dz;
        Double_t zmax=t->GetZ() + dz;
        Double_t ymin=t->GetY() + phi*r - dy;
        Double_t ymax=t->GetY() + phi*r + dy;
        layer.SelectClusters(zmin,zmax,ymin,ymax);

        const AliITSRecPoint *c=0; Int_t ci=-1,cci=-1;
        Double_t maxchi2=1000.*kMaxChi2, tolerance=0.1;
        while ((c=layer.GetNextCluster(ci))!=0) {
           if (idet == c->GetDetectorIndex()) continue;

	   const AliITSdetector &det=layer.GetDetector(c->GetDetectorIndex());

	   if (!tmp.Propagate(det.GetPhi(),det.GetR())) continue;
           
	   if (TMath::Abs(tmp.GetZ() - c->GetZ()) > tolerance) continue;
           if (TMath::Abs(tmp.GetY() - c->GetY()) > tolerance) continue;

           Double_t chi2=tmp.GetPredictedChi2(c);
           if (chi2<maxchi2) { maxchi2=chi2; cci=ci; }
        }
        if (cci>=0) t->SetExtraCluster(i,(i<<28)+cci);
     }

     // track time update [SR, GSI 17.02.2003]
     if (t->IsStartedTimeIntegral() && step==1) {
        Double_t newX, newY, newZ;
        t->GetGlobalXYZat(t->GetX(),newX,newY,newZ);
        Double_t dL2 = (oldX-newX)*(oldX-newX) + (oldY-newY)*(oldY-newY) + 
                       (oldZ-newZ)*(oldZ-newZ);
        t->AddTimeStep(TMath::Sqrt(dL2));
     }
     //

  }

  if (!t->PropagateTo(xx,0.,0.)) return kFALSE;
  return kTRUE;
}

void AliITStrackerV2::UseClusters(const AliKalmanTrack *t, Int_t from) const {
  //--------------------------------------------------------------------
  // This function marks clusters assigned to the track
  //--------------------------------------------------------------------
  AliTracker::UseClusters(t,from);

  Int_t clusterIndex = t->GetClusterIndex(0);
  AliITSRecPoint *c= 0x0;

  if (clusterIndex>-1)
    c = (AliITSRecPoint *)GetCluster(clusterIndex);
  if (c && c->GetSigmaZ2()>0.1) c->UnUse();

  c = 0x0;
  clusterIndex = t->GetClusterIndex(1);
  if (clusterIndex>-1)
    c=(AliITSRecPoint *)GetCluster(clusterIndex);
  if (c && c->GetSigmaZ2()>0.1) c->UnUse();

}
