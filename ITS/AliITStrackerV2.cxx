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
//    It reads AliITSclusterV2 clusters and creates AliITStrackV2 tracks
//                   and fills with them the ESD
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------

#include <new>

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>

#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliESD.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerV2.h"

ClassImp(AliITStrackerV2)

AliITStrackerV2::AliITSlayer AliITStrackerV2::fgLayers[kMaxLayer]; // ITS layers

AliITStrackerV2::AliITStrackerV2(const AliITSgeom *geom) : AliTracker() {
  //--------------------------------------------------------------------
  //This is the AliITStrackerV2 constructor
  //--------------------------------------------------------------------
  AliITSgeom *g=(AliITSgeom*)geom;

  Float_t x,y,z;
  Int_t i;
  for (i=1; i<kMaxLayer+1; i++) {
    Int_t nlad=g->GetNladders(i);
    Int_t ndet=g->GetNdetectors(i);

    g->GetTrans(i,1,1,x,y,z); 
    Double_t r=TMath::Sqrt(x*x + y*y);
    Double_t poff=TMath::ATan2(y,x);
    Double_t zoff=z;

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
        Double_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);
        Double_t r=x*cp+y*sp;

        AliITSdetector &det=fgLayers[i-1].GetDetector((j-1)*ndet + k-1); 
        new(&det) AliITSdetector(r,phi); 
      } 
    }  

  }

  fI=kMaxLayer;

  fPass=0;
  fConstraint[0]=1; fConstraint[1]=0;

  Double_t xyz[]={kXV,kYV,kZV}, ers[]={kSigmaXV,kSigmaYV,kSigmaZV}; 
  SetVertex(xyz,ers);

  for (Int_t i=0; i<kMaxLayer; i++) fLayersNotToSkip[i]=kLayersNotToSkip[i];
  fLastLayerToTrackTo=kLastLayerToTrackTo;

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
  TBranch *branch=cTree->GetBranch("Clusters");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch !\n");
    return 1;
  }

  TClonesArray dummy("AliITSclusterV2",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  Int_t j=0;
  Int_t detector=0;
  for (Int_t i=0; i<kMaxLayer; i++) {
    Int_t ndet=fgLayers[i].GetNdetectors();
    Int_t jmax = j + fgLayers[i].GetNladders()*ndet;
    for (; j<jmax; j++) {           
      if (!cTree->GetEvent(j)) continue;
      Int_t ncl=clusters->GetEntriesFast();
      while (ncl--) {
        AliITSclusterV2 *c=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
	detector = c->GetDetectorIndex();
        fgLayers[i].InsertCluster(new AliITSclusterV2(*c));
      }
      clusters->Delete();
      //add dead zone virtual "cluster"
      
      if (i<1){
	for (Float_t ydead = -2.; ydead < 2. ; ydead+=0.05){     
	  Int_t lab[4] = {0,0,0,detector};
	  Int_t info[3] = {0,0,0};
	  Float_t hit[5]={ydead,0,1,0.01,0};
	  fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	  //
	  hit[1]=-7.;
	  fgLayers[i].InsertCluster(new AliITSclusterV2(lab, hit, info));
	}	
      }
      
    }
    //
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
      if (TMath::Abs(t->GetD())>5) {
	delete t;
	continue;
      }

      if (CorrectForDeadZoneMaterial(t)!=0) {
         Warning("Clusters2Tracks",
                 "failed to correct for the material in the dead zone !\n");
         delete t;
         continue;
      }
      t->SetReconstructed(kFALSE);
      itsTracks.AddLast(t);
    }
  } /* End Read ESD tracks */

  itsTracks.Sort();
  Int_t nentr=itsTracks.GetEntriesFast();
  fTrackHypothesys.Expand(nentr);
  Int_t ntrk=0;
  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
       fCurrentEsdTrack = i;
       AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
       if (t==0) continue;              //this track has been already tracked
       if (t->GetReconstructed()) continue;  //this track was  already  "succesfully" reconstructed
       if ( (TMath::Abs(t->GetD(GetX(),GetY()))  >2.) && fConstraint[fPass]) continue;
       if ( (TMath::Abs(t->GetZat(GetX())-GetZ())>2.) && fConstraint[fPass]) continue;

       Int_t tpcLabel=t->GetLabel(); //save the TPC track label

       ResetTrackToFollow(*t);
       ResetBestTrack();

       for (FollowProlongation(); fI<kMaxLayer; fI++) {
	 while (TakeNextProlongation()) {
	    FollowProlongation();
	 }
       }

       SortTrackHypothesys(fCurrentEsdTrack,0.9);  //MI change
       CompressTrackHypothesys(fCurrentEsdTrack,0.90,2);  //MI change

       /*
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
       //UseClusters(&fBestTrack);
       delete itsTracks.RemoveAt(i);
       ntrk++;
       */
       AliITStrackV2 * besttrack = GetBestHypothesys(fCurrentEsdTrack,t,5);
       if (!besttrack) continue;
       besttrack->SetLabel(tpcLabel);
       besttrack->CookdEdx();
       besttrack->SetFakeRatio(1.);
       CookLabel(besttrack,0.); //For comparison only
       //       besttrack->UpdateESDtrack(AliESDtrack::kITSin);
       
       if (besttrack->GetChi2()/besttrack->GetNumberOfClusters()>3.5){
	 if ( (TMath::Abs(besttrack->GetD(GetX(),GetY()))>0.4) && fConstraint[fPass]) {
	   CompressTrackHypothesys(fCurrentEsdTrack,0.0,0);
	   continue;
	 }
	 if ( (TMath::Abs(besttrack->GetZat(GetX()) -GetZ() )>0.4) && fConstraint[fPass]){
	   CompressTrackHypothesys(fCurrentEsdTrack,0.0,0);
	   continue;
	 }
       }
       
       //delete itsTracks.RemoveAt(i);
       t->SetReconstructed();
       ntrk++;              
       
     }
  }
  
  for (Int_t i=0; i<nentr; i++) {
    fCurrentEsdTrack = i;
    AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
    if (t==0) continue;         
    Int_t tpcLabel=t->GetLabel(); //save the TPC track label
    AliITStrackV2 * besttrack = GetBestHypothesysMIP(fCurrentEsdTrack,t);
    if (!besttrack) continue;

    besttrack->SetLabel(tpcLabel);
    besttrack->CookdEdx();
    besttrack->SetFakeRatio(1.);
    CookLabel(besttrack,0.); //For comparison only
    besttrack->UpdateESDtrack(AliESDtrack::kITSin);
  }
  //

  itsTracks.Delete();
  //
  Int_t entries = fTrackHypothesys.GetEntriesFast();
  for (Int_t ientry=0;ientry<entries;ientry++){
    TObjArray * array =(TObjArray*)fTrackHypothesys.UncheckedAt(ientry);
    if (array) array->Delete();
    delete fTrackHypothesys.RemoveAt(ientry); 
  }

  fTrackHypothesys.Delete();
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
       if (fTrackToFollow.PropagateToVertex()) {
          fTrackToFollow.StartTimeIntegral();
       }
       fTrackToFollow.PropagateTo(3.,-0.0028,65.19);
     }

     fTrackToFollow.ResetCovariance(); fTrackToFollow.ResetClusters();
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
        //UseClusters(&fTrackToFollow);
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

    if ((esd->GetStatus()&AliESDtrack::kTPCin)==0)
      fTrackToFollow.ResetCovariance();

    //Refitting...
    if (RefitAt(3.7, &fTrackToFollow, t)) {
       fTrackToFollow.SetLabel(t->GetLabel());
       fTrackToFollow.CookdEdx();
       CookLabel(&fTrackToFollow,0.0); //For comparison only

       if (fTrackToFollow.PropagateTo(3.,0.0028,65.19)) {//The beam pipe    
         Double_t a=fTrackToFollow.GetAlpha();
         Double_t cs=TMath::Cos(a),sn=TMath::Sin(a);
         Double_t xv= GetX()*cs + GetY()*sn;
         Double_t yv=-GetX()*sn + GetY()*cs;
         
         Double_t c=fTrackToFollow.GetC(), snp=fTrackToFollow.GetSnp();
         Double_t x=fTrackToFollow.GetX(), y=fTrackToFollow.GetY();
         Double_t tgfv=-(c*(x-xv)-snp)/(c*(y-yv) + TMath::Sqrt(1.-snp*snp));
         Double_t fv=TMath::ATan(tgfv);

         cs=TMath::Cos(fv); sn=TMath::Sin(fv);
         x = xv*cs + yv*sn;
         yv=-xv*sn + yv*cs; xv=x;

	 if (fTrackToFollow.Propagate(fv+a,xv)) {
            fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSrefit);
            //UseClusters(&fTrackToFollow);
            {
            AliITSclusterV2 c; c.SetY(yv); c.SetZ(GetZ());
            c.SetSigmaY2(GetSigmaY()*GetSigmaY());
            c.SetSigmaZ2(GetSigmaZ()*GetSigmaZ());
            Double_t chi2=fTrackToFollow.GetPredictedChi2(&c);
            if (chi2<kMaxChi2)
	       if (fTrackToFollow.Update(&c,-chi2,0))
                   fTrackToFollow.SetConstrainedESDtrack(chi2);            
            }
            ntrk++;
         }
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

    //    if (TMath::Abs(fTrackToFollow.GetZ()-GetZ()) > r+dz) return;

    //Double_t dy=4*TMath::Sqrt(track.GetSigmaY2() + kSigmaY2[i]);
    if (TMath::Abs(track.GetSnp()>kMaxSnp)) {
      fI--;
      continue;   // MI
    }
    if (dy < 0.5*TMath::Abs(track.GetSnp())) dy=0.5*TMath::Abs(track.GetSnp());
    if (dy > kMaxRoad) {
      //Warning("FollowProlongation","too broad road in Y !\n");
      return;
    }

    Double_t zmin=track.GetZ() - dz; 
    Double_t zmax=track.GetZ() + dz;
    Double_t ymin=track.GetY() + r*phi - dy;
    Double_t ymax=track.GetY() + r*phi + dy;
    layer.SelectClusters(zmin,zmax,ymin,ymax); 
    fI--;

    //take another prolongation
    if (!TakeNextProlongation()){ 
      //skipped++;
      fTrackToFollow.IncrementNSkipped();
      if (fLayersNotToSkip[fI]||fTrackToFollow.GetNSkipped()>1) return;
    }
    if (fTrackToFollow.GetNUsed()>1) return;
    if (fTrackToFollow.GetNUsed()+fTrackToFollow.GetNSkipped()>1) return;    
    if ( (fI<3) && ( fTrackToFollow.GetChi2()/fTrackToFollow.GetNumberOfClusters()>kChi2PerCluster)) return;
  } 

  //deal with the best track
  Int_t ncl=fTrackToFollow.GetNumberOfClusters();
  Int_t nclb=fBestTrack.GetNumberOfClusters();
  if (ncl){
    if (ncl<4) return;
    if ( (ncl<6) && (fTrackToFollow.GetChi2()/float(ncl))>3) return;
    if (fTrackToFollow.GetChi2()/ncl>5.5) return;
    fTrackToFollow.CookdEdx();
    if (fTrackToFollow.GetESDtrack()->GetTPCsignal()>80.)
      if ((fTrackToFollow.GetdEdx()/fTrackToFollow.GetESDtrack()->GetTPCsignal())<0.35){
	// mismatch in dedx
	return;
      }
    //
    fTrackToFollow.SetLabel(fTrackToFollow.GetESDtrack()->GetLabel());
    CookLabel(&fTrackToFollow,0.); //
    //
    if ( (nclb>3) && ((fTrackToFollow.GetChi2()/ncl)<(3*fBestTrack.GetChi2()/(nclb))))
      AddTrackHypothesys(new AliITStrackV2(fTrackToFollow), fCurrentEsdTrack);
    else 
      if (ncl>3) AddTrackHypothesys(new AliITStrackV2(fTrackToFollow), fCurrentEsdTrack);
    if (ncl >= nclb) {
      Double_t chi2=fTrackToFollow.GetChi2();
      if (chi2/ncl < kChi2PerCluster) {        
        if (ncl > nclb ) {
	  ResetBestTrack();
        }
	if ( (ncl == nclb) && chi2 < fBestTrack.GetChi2()) {
	  ResetBestTrack();
        }	
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

  const AliITSclusterV2 *c=0; Int_t ci=-1;
  Double_t chi2=12345.;
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

    chi2=fTrackToFollow.GetPredictedChi2(c); if (chi2<kMaxChi2) break;
  }

  if (chi2>=kMaxChi2) return 0;
  if (!c) return 0;
  if (c->IsUsed()&&c->GetNy()<5) {  //shared factor    
    chi2+=1;
    chi2*=2*(1./(TMath::Max(c->GetNy(),1)));
  }
  if (c->GetQ()==0){  //dead zone virtual cluster
    chi2*=4.;
    chi2+=20;
    fTrackToFollow.IncrementNUsed();
    return 1;
  }
  //if ((fI<2)&&chi2>kMaxChi2In) return 0;

  if (!fTrackToFollow.Update(c,chi2,(fI<<28)+ci)) {
     //Warning("TakeNextProlongation","filtering failed !\n");
     return 0;
  }
  if (c->IsUsed()&&c->GetNy()<5) fTrackToFollow.IncrementNUsed();
  if (fTrackToFollow.GetNumberOfClusters()>1)
  if (TMath::Abs(fTrackToFollow.GetD())>4) return 0;

  fTrackToFollow.
    SetSampledEdx(c->GetQ(),fTrackToFollow.GetNumberOfClusters()-1); //b.b.

  {
  Double_t x0;
 Double_t d=layer.GetThickness(fTrackToFollow.GetY(),fTrackToFollow.GetZ(),x0);
  fTrackToFollow.CorrectForMaterial(d,x0);
  }

  if (fConstraint[fPass]) {
    Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!
    Double_t xyz[]={GetX(),GetY(),GetZ()};
    Double_t ers[]={GetSigmaX(),GetSigmaY(),GetSigmaZ()};
    Double_t deltad = TMath::Abs(fTrackToFollow.GetD(GetX(),GetY()));
    Double_t deltaz = TMath::Abs(fTrackToFollow.GetZat(GetX())-GetZ());

    if ( (fI==4) &&  (deltad>2.0 || deltaz>1.5))  return 0; // don't "improve" secondaries     
    if ( (fI==3) &&  (deltad>1.5 || deltaz>0.9))  return 0; // don't "improve" secondaries 
    if ( (fI==2) &&  (deltad>0.9 || deltaz>0.6))  return 1; // don't "improve" secondaries 
    if ( (fI==1) &&  (deltad>0.3 || deltaz>0.3))  return 1; // don't "improve" secondaries 
    if ( (fI==0) &&  (deltad>0.1 || deltaz>0.1))  return 1; // don't "improve" secondaries 

    fTrackToFollow.Improve(d,xyz,ers);
  }

  return 1;
}


AliITStrackerV2::AliITSlayer::AliITSlayer() {
  //--------------------------------------------------------------------
  //default AliITSlayer constructor
  //--------------------------------------------------------------------
  fN=0;
  fDetectors=0;
}

AliITStrackerV2::AliITSlayer::
AliITSlayer(Double_t r,Double_t p,Double_t z,Int_t nl,Int_t nd) {
  //--------------------------------------------------------------------
  //main AliITSlayer constructor
  //--------------------------------------------------------------------
  fR=r; fPhiOffset=p; fZOffset=z;
  fNladders=nl; fNdetectors=nd;
  fDetectors=new AliITSdetector[fNladders*fNdetectors];

  fN=0;
  fI=0;

  fRoad=2*fR*TMath::Sqrt(3.14/1.);//assuming that there's only one cluster
}

AliITStrackerV2::AliITSlayer::~AliITSlayer() {
  //--------------------------------------------------------------------
  // AliITSlayer destructor
  //--------------------------------------------------------------------
  delete[] fDetectors;
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) fClusterWeight[i]=0;
}

void AliITStrackerV2::AliITSlayer::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) fClusterWeight[i]=0;
  fN=0;
  fI=0;
}

void AliITStrackerV2::AliITSlayer::ResetWeights() {
  //--------------------------------------------------------------------
  // This function reset weights of the clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxClusterPerLayer;i++) fClusterWeight[i]=0;
}

void AliITStrackerV2::AliITSlayer::ResetRoad() {
  //--------------------------------------------------------------------
  // This function calculates the road defined by the cluster density
  //--------------------------------------------------------------------
  Int_t n=0;
  for (Int_t i=0; i<fN; i++) {
     if (TMath::Abs(fClusters[i]->GetZ())<fR) n++;
  }
  if (n>1) fRoad=2*fR*TMath::Sqrt(3.14/n);
}

Int_t AliITStrackerV2::AliITSlayer::InsertCluster(AliITSclusterV2 *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to this layer
  //--------------------------------------------------------------------
  if (fN==kMaxClusterPerLayer) {
    ::Error("InsertCluster","Too many clusters !\n");
    return 1;
  }

  if (fN==0) {fClusters[fN++]=c; return 0;}
  Int_t i=FindClusterIndex(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliITSclusterV2*));
  fClusters[i]=c; fN++;

  return 0;
}

Int_t AliITStrackerV2::AliITSlayer::FindClusterIndex(Double_t z) const {
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

void AliITStrackerV2::AliITSlayer::
SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin, Double_t ymax) {
  //--------------------------------------------------------------------
  // This function sets the "window"
  //--------------------------------------------------------------------
  fI=FindClusterIndex(zmin); fZmax=zmax;
  Double_t circle=2*TMath::Pi()*fR;
  if (ymax>circle) { ymax-=circle; ymin-=circle; }
  fYmin=ymin; fYmax=ymax;
}

const AliITSclusterV2 *AliITStrackerV2::AliITSlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------
  const AliITSclusterV2 *cluster=0;
  for (Int_t i=fI; i<fN; i++) {
    const AliITSclusterV2 *c=fClusters[i];
    if (c->GetZ() > fZmax) break;
    //    if (c->IsUsed()) continue;
    const AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + c->GetY();

    if (y>2.*fR*TMath::Pi()) y -= 2*fR*TMath::Pi();
    if (y>1.*fR*TMath::Pi() && fYmax<y) y -= 2*fR*TMath::Pi();

    if (y<fYmin) continue;
    if (y>fYmax) continue;
    cluster=c; ci=i;
    fI=i+1;
    break; 
  }

  return cluster;
}

Int_t AliITStrackerV2::AliITSlayer::
FindDetectorIndex(Double_t phi, Double_t z) const {
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

Int_t AliITStrackerV2::AliITSlayer::InRoad() const {
  //--------------------------------------------------------------------
  // This function returns number of clusters within the "window" 
  //--------------------------------------------------------------------
  Int_t ncl=0;
  for (Int_t i=fI; i<fN; i++) {
    const AliITSclusterV2 *c=fClusters[i];
    if (c->GetZ() > fZmax) break;
    if (c->IsUsed()) continue;
    const AliITSdetector &det=GetDetector(c->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + c->GetY();

    if (y>2.*fR*TMath::Pi()) y -= 2*fR*TMath::Pi();
    if (y>1.*fR*TMath::Pi() && fYmax<y) y -= 2*fR*TMath::Pi();

    if (y<fYmin) continue;
    if (y>fYmax) continue;
    ncl++;
  }
  return ncl;
}

Bool_t 
AliITStrackerV2::RefitAt(Double_t xx,AliITStrackV2 *t,const AliITStrackV2 *c) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"
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

     const AliITSclusterV2 *cl=0;
     Double_t maxchi2=kMaxChi2;

     Int_t idx=index[i];
     if (idx>0) {
        const AliITSclusterV2 *c=(AliITSclusterV2 *)GetCluster(idx); 
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
     /*
     if (cl==0)
     if (t->GetNumberOfClusters()>2) {
        Double_t dz=4*TMath::Sqrt(t->GetSigmaZ2()+kSigmaZ2[i]);
        Double_t dy=4*TMath::Sqrt(t->GetSigmaY2()+kSigmaY2[i]);
        Double_t zmin=t->GetZ() - dz;
        Double_t zmax=t->GetZ() + dz;
        Double_t ymin=t->GetY() + phi*r - dy;
        Double_t ymax=t->GetY() + phi*r + dy;
        layer.SelectClusters(zmin,zmax,ymin,ymax);

        const AliITSclusterV2 *c=0; Int_t ci=-1;
        while ((c=layer.GetNextCluster(ci))!=0) {
           if (idet != c->GetDetectorIndex()) continue;
           Double_t chi2=t->GetPredictedChi2(c);
           if (chi2<maxchi2) { cl=c; maxchi2=chi2; idx=ci; }
        }
     }
     */
     if (cl) {
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


Float_t  *AliITStrackerV2::GetWeight(Int_t index) {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fgLayers[l].GetWeight(c);
}



void AliITStrackerV2::UseClusters(const AliKalmanTrack *t, Int_t from) const {
  //--------------------------------------------------------------------
  // This function marks clusters assigned to the track
  //--------------------------------------------------------------------
  AliTracker::UseClusters(t,from);

  AliITSclusterV2 *c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(0));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();
  c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(1));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();

}


void AliITStrackerV2::AddTrackHypothesys(AliITStrackV2 * track, Int_t esdindex)
{
  //------------------------------------------------------------------
  // add track to the list of hypothesys
  //------------------------------------------------------------------

  if (esdindex>=fTrackHypothesys.GetEntriesFast()) fTrackHypothesys.Expand(esdindex*2+10);
  //
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) {
    array = new TObjArray(10);
    fTrackHypothesys.AddAt(array,esdindex);
  }
  array->AddLast(track);
}

void AliITStrackerV2::SortTrackHypothesys(Int_t esdindex, Float_t likelihoodlevel)
{
  //-------------------------------------------------------------------
  // compress array of track hypothesys
  // keep only maxsize best hypothesys
  //-------------------------------------------------------------------
  if (esdindex>fTrackHypothesys.GetEntriesFast()) return;
  if (! (fTrackHypothesys.At(esdindex)) ) return;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  Int_t entries = array->GetEntriesFast();

  Float_t * chi2        = new Float_t[entries];
  Float_t * probability = new Float_t[entries];
  Float_t sumprobabilityall=0;
  Int_t * index         = new Int_t[entries];
  //
  //
  for (Int_t itrack=0;itrack<array->GetEntriesFast();itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    //
    if (track && track->GetNumberOfClusters()>(track->GetNUsed()+track->GetNSkipped())){
      //
      chi2[itrack] = track->GetChi2()/(track->GetNumberOfClusters()-track->GetNUsed()-track->GetNSkipped());      
      if (track->GetESDtrack())
	if (track->GetESDtrack()->GetTPCsignal()>80){
	  track->CookdEdx();
	  if ((track->GetdEdx()/track->GetESDtrack()->GetTPCsignal())<0.4){
	    Float_t factor = 2.+10.*(0.6-track->GetdEdx()/track->GetESDtrack()->GetTPCsignal());
	    chi2[itrack]*= factor;  //mismatch in dEdx
	  }
	}      
    }
    else
      chi2[itrack] = 10000000.;
    probability[itrack] = 1./(0.3+chi2[itrack]);
    sumprobabilityall+=probability[itrack];
  }

  TMath::Sort(entries,chi2,index,kFALSE);
  TObjArray * newarray = new TObjArray();
  Float_t     sumprobability = 0.;
  for (Int_t i=0;i<entries;i++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(index[i]);
    if (!track) break;
    if (chi2[index[i]]<30){
      newarray->AddLast(array->RemoveAt(index[i])); 
      track->SetChi2MIP(0,chi2[index[i]]);     // base chi 2    
      sumprobability+= probability[index[i]]/sumprobabilityall;
      if (sumprobability> likelihoodlevel) break;
    }
    else{
      delete array->RemoveAt(index[i]);
    }
  }

  array->Delete();
  delete fTrackHypothesys.RemoveAt(esdindex);
  fTrackHypothesys.AddAt(newarray,esdindex);

  delete [] chi2;
  delete [] index;

}


void AliITStrackerV2::CompressTrackHypothesys(Int_t esdindex, Float_t likelihoodlevel, Int_t maxsize)
{
  //
  //
  if (esdindex>fTrackHypothesys.GetEntriesFast()) return;
  if (! (fTrackHypothesys.At(esdindex)) ) return;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  Int_t entries = array->GetEntriesFast();
  //
  if (likelihoodlevel>0.000001){
    Float_t *probability = new Float_t[entries];
    for (Int_t i=0;i<entries;i++) probability[i]=0.;
    Float_t  sumprobabilityall=0;
    Int_t   *index         = new Int_t[entries];
    
    for (Int_t itrack=0;itrack<entries;itrack++){
      AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
      probability[itrack]=0;
      if (!track) continue;
      probability[itrack] = 1./(0.3+track->GetChi2MIP(0));
      sumprobabilityall  += probability[itrack];
      //    
    }
    if (sumprobabilityall<=0.000000000000001){
      return;
    }
    for (Int_t itrack=0;itrack<entries;itrack++){
      probability[itrack]/=sumprobabilityall;
    }
    
    TMath::Sort(entries, probability, index, kTRUE);
    //
    TObjArray * newarray = new TObjArray();
    Float_t     sumprobability = 0.;
    for (Int_t i=0;i<entries;i++){
      AliITStrackV2 * track = (AliITStrackV2*)array->At(index[i]);
      if (!track) continue;    
      newarray->AddLast(array->RemoveAt(index[i]));      
      sumprobability+= probability[index[i]];
      if (sumprobability> likelihoodlevel) break;
      if (i>maxsize) break;
    }
    
    array->Delete();
    delete fTrackHypothesys.RemoveAt(esdindex);
    fTrackHypothesys.AddAt(newarray,esdindex);
    //    
    delete []index;
    delete []probability;
  }
  else{
    array->Delete();
    delete fTrackHypothesys.RemoveAt(esdindex);
  }
}


AliITStrackV2 * AliITStrackerV2::GetBestHypothesys(Int_t esdindex, AliITStrackV2 * original, Int_t checkmax)
{
  //-------------------------------------------------------------
  // try to find best hypothesy
  // currently - minimal chi2 of track+backpropagated track+matching to the tpc track
  //-------------------------------------------------------------
  if (fTrackHypothesys.GetEntriesFast()<=esdindex) return 0;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) return 0;
  Int_t entries = array->GetEntriesFast();
  if (!entries) return 0;  
  Float_t minchi2 = 100000;
  Int_t   maxn    = 3;
  AliITStrackV2 * besttrack=0;
  Int_t accepted =0;
  Int_t maxindex=0;
  //
  Float_t sumz2=0;
  Float_t sumy2=0;
  Float_t sumall=0;
  for (Int_t itrack=0;itrack<entries; itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track) continue;
    sumall++;
    sumz2+=track->GetSigmaZ2();
    sumy2+=track->GetSigmaY2();
  }
  sumz2/=sumall;
  sumy2/=sumall;

  Float_t dedxmismatch=1;
  for (Int_t i=0;i<entries;i++){    
    maxindex = i;
    AliITStrackV2 * track = (AliITStrackV2*)array->At(i);    
    if (!track) continue;
    track->SetChi2MIP(1,1000000);
    track->SetChi2MIP(2,1000000);

    if ( (track->GetNumberOfClusters()-track->GetNSkipped()-track->GetNUsed())<2) continue;
    //
    if (track->GetESDtrack())
      if (track->GetESDtrack()->GetTPCsignal()>80){
	track->CookdEdx(); 	
	if ((track->GetdEdx()/track->GetESDtrack()->GetTPCsignal())<0.4){
	  //mismatch in dEdx 
	  dedxmismatch= 2.+ 10.*(0.6-track->GetdEdx()/track->GetESDtrack()->GetTPCsignal());
	}
      }

    //    track->SetLabel(original->GetLabel());
    //CookLabel(track,0.0);
    //if (track->GetFakeRatio()>0.01) continue;
    //
    //
    // backtrack
    AliITStrackV2 * backtrack = new AliITStrackV2(*track);
    backtrack->ResetCovariance();
    backtrack->ResetClusters();
    Double_t x = original->GetX();
    if (!RefitAt(x,backtrack,track)){
      delete backtrack;
      delete array->RemoveAt(i);
      continue;
    }
    if (  (backtrack->GetChi2() / float(backtrack->GetNumberOfClusters()-track->GetNSkipped()-track->GetNUsed()-0.5))>6) 
      {
	delete backtrack; 
	delete array->RemoveAt(i);
	continue;
      }
    Double_t deltac   = backtrack->GetC()-original->GetC();
    Double_t deltatgl = backtrack->GetTgl()-original->GetTgl();
    //
    Double_t poolc2      = (deltac*deltac)/(original->GetCov44()+backtrack->GetCov44());
    Double_t pooltgl2    = (deltatgl*deltatgl)/(original->GetCov33()+backtrack->GetCov33());
    if ((poolc2+pooltgl2)>32){  //4 sigma
      delete backtrack; 
      delete array->RemoveAt(i);
      continue;      
    }
    //Double_t bpoolc      = (deltac*deltac)/(original->GetCov44());
    //Double_t bpooltgl    = (deltatgl*deltatgl)/(original->GetCov33());

    //
    //forward track - without constraint
    AliITStrackV2 * forwardtrack = new AliITStrackV2(*original);
    //    forwardtrack->ResetCovariance();
    forwardtrack->ResetClusters();
    x = track->GetX();
    if (!RefitAt(x,forwardtrack,track)){
      delete forwardtrack;
      delete backtrack;
      delete array->RemoveAt(i);
      continue;
    }
    if ( (forwardtrack->GetChi2()/float(forwardtrack->GetNumberOfClusters()-track->GetNSkipped()-track->GetNUsed()))>6) 
      {
	delete forwardtrack; 
	delete array->RemoveAt(i);
	continue;
      }
    //
    accepted++;
    if (accepted>checkmax){
      delete backtrack;
      delete forwardtrack;
      break;
    }
    Double_t chi2 = (backtrack->GetChi2()/(backtrack->GetNumberOfClusters()-1-track->GetNSkipped()-track->GetNUsed())+
		     forwardtrack->GetChi2()/(forwardtrack->GetNumberOfClusters()-track->GetNSkipped()-track->GetNUsed()));
      //                     bpoolc+bpooltgl;
    //    chi2 *= (forwardtrack->GetSigmaZ2()/sumz2+forwardtrack->GetSigmaY2()/sumy2);
    chi2 *= dedxmismatch;
    //
    //
    track->SetChi2MIP(1,backtrack->GetChi2()/(backtrack->GetNumberOfClusters()-1-track->GetNSkipped()-track->GetNUsed()));
    track->SetChi2MIP(2,forwardtrack->GetChi2()/(forwardtrack->GetNumberOfClusters()-track->GetNSkipped()-track->GetNUsed()));
    track->SetChi2MIP(3,poolc2+pooltgl2);
    //
    
    if (track->GetNumberOfClusters()>maxn){
      besttrack =  new AliITStrackV2(*forwardtrack);
      maxn      =  track->GetNumberOfClusters();
      minchi2   =  chi2;
      delete backtrack;
      delete forwardtrack;
      continue;
    }   
    //
    if (chi2 < minchi2){
      besttrack = new AliITStrackV2(*forwardtrack);
      minchi2   = chi2;      
    }    
    delete backtrack;
    delete forwardtrack;
  }
  //
  //
  if (!besttrack || besttrack->GetNumberOfClusters()<4) {
    return 0;
  }

  //
  besttrack->SetLabel(original->GetLabel());
  CookLabel(besttrack,0.0);
  //
  // calculate "weight of the cluster"
  //

  {
    //sign usage information for clusters
    Int_t clusterindex[6][100];
    Double_t clusterweight[6][100];
    for (Int_t ilayer=0;ilayer<6;ilayer++)
      for (Int_t icluster=0;icluster<100;icluster++){
	clusterindex[ilayer][icluster]   = -1;
	clusterweight[ilayer][icluster]  = 0;
      }
    //printf("%d\t%d\n",esdindex, entries);
    //
    Float_t sumchi2=0;
    for (Int_t itrack=0;itrack<entries; itrack++){
      AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
      if (!track) continue;
      if (track->GetChi2MIP(1)>1000) continue;  //not accepted
      sumchi2 +=1./(0.3+track->GetChi2MIP(1)+track->GetChi2MIP(2));
    }
    for (Int_t itrack=0;itrack<entries;itrack++){
      AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
      if (!track) continue;
      if (track->GetChi2MIP(1)>1000) continue;  //not accepted  
      for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){     
	Int_t tindex = track->GetClusterIndex(icluster);
	Int_t ilayer = (tindex & 0xf0000000) >> 28;
	if (tindex<0) continue;
	Int_t cindex =0;
	//
	for (cindex=0;cindex<100;cindex++){
	  if (clusterindex[ilayer][cindex]<0) break;
	  if (clusterindex[ilayer][cindex]==tindex) break;
	}
	if (cindex>100) break;
	if (clusterindex[ilayer][cindex]!=tindex) clusterindex[ilayer][cindex] = tindex;
	clusterweight[ilayer][cindex]+= (1./(0.3+track->GetChi2MIP(1)+track->GetChi2MIP(2)))* (1./sumchi2); 
	Float_t *weight = GetWeight(tindex);
	
	if (weight){
	  *weight+= (1./(0.3+track->GetChi2MIP(1)+track->GetChi2MIP(2)))* (1./sumchi2);
	}
      }
    }
    
    if (besttrack->GetChi2()/besttrack->GetNumberOfClusters()>3.5) return besttrack; //don't sign clusters
    Int_t current=0;
    Double_t deltad = besttrack->GetD(GetX(),GetY());
    Double_t deltaz = besttrack->GetZat(GetX()) - GetZ();
    Double_t deltaprim = TMath::Sqrt(deltad*deltad+deltaz*deltaz);

    for (Int_t icluster=0;icluster<besttrack->GetNumberOfClusters();icluster++){
      Int_t index = besttrack->GetClusterIndex(icluster);
      Int_t ilayer =  (index & 0xf0000000) >> 28;
      AliITSclusterV2 *c = (AliITSclusterV2*)GetCluster(index);
      if (!c) continue;
      // 
      for (Int_t icluster=0;icluster<100;icluster++){
	// sign non "doubt" clusters as used
	if (clusterindex[ilayer][icluster]!=index) continue;
	//	Float_t * weight = GetWeight(index);
	//if (weight) if (*weight>1){
	//  if (c->IsUsed()) continue;
	//  c->Use();
	//}      
	if ( (ilayer*0.2+0.2)<deltaprim) continue; // secondaries
	if (c->GetNy()>4) continue; // don sign cluster
	if ( (ilayer>1&&clusterweight[ilayer][icluster]>0.7) || (ilayer<2&&clusterweight[ilayer][icluster]>0.8) ){
	  current++;
	  if (c->IsUsed()) continue;
	  c->Use();
	}
      }
    }
  }
 
  //
  return besttrack;
} 


AliITStrackV2 * AliITStrackerV2::GetBestHypothesysMIP(Int_t esdindex, AliITStrackV2 * original)
{
  //-------------------------------------------------------------
  // try to find best hypothesy
  // currently - minimal chi2 of track+backpropagated track+matching to the tpc track
  //-------------------------------------------------------------
  if (fTrackHypothesys.GetEntriesFast()<=esdindex) return 0;
  TObjArray * array = (TObjArray*) fTrackHypothesys.At(esdindex);
  if (!array) return 0;
  Int_t entries = array->GetEntriesFast();
  if (!entries) return 0;  
  AliITStrackV2 * besttrack=0;
  //
  //sign usage information for clusters
  Int_t clusterindex[6][100];
  Double_t clusterweight[6][100];
  for (Int_t ilayer=0;ilayer<6;ilayer++)
    for (Int_t icluster=0;icluster<100;icluster++){
      clusterindex[ilayer][icluster]   = -1;
      clusterweight[ilayer][icluster]  = 0;
    }
  //printf("%d\t%d\n",esdindex, entries);
  //
  Float_t sumchi2=0;
  for (Int_t itrack=0;itrack<entries; itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track) continue;
    if (track->GetChi2MIP(1)>1000) continue;  //not accepted - before
    sumchi2 +=1./(0.3+track->GetChi2MIP(1)+track->GetChi2MIP(2));
  }
  //
  // get cluster weight 
  for (Int_t itrack=0;itrack<entries;itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track) continue;
    if (track->GetChi2MIP(1)>1000) continue;  // track not accepted in previous iterration  
    for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){     
      Int_t tindex = track->GetClusterIndex(icluster);
      Int_t ilayer = (tindex & 0xf0000000) >> 28;
      if (tindex<0) continue;
      Int_t cindex =0;
      //
      for (cindex=0;cindex<100;cindex++){
	if (clusterindex[ilayer][cindex]<0) break;
	if (clusterindex[ilayer][cindex]==tindex) break;
      }
      if (cindex>100) break;
      if (clusterindex[ilayer][cindex]!=tindex) clusterindex[ilayer][cindex] = tindex;
      clusterweight[ilayer][cindex]+= (1./(0.3+track->GetChi2MIP(1)+track->GetChi2MIP(2)))* (1./sumchi2); 
    }
  }
  //  
  // get cluster relative sharing - factor
  //
  // 
  for (Int_t ilayer=0;ilayer<6;ilayer++)
    for (Int_t cindex=0;cindex<100;cindex++){
      if (clusterindex[ilayer][cindex]<0) continue;
      Int_t tindex = clusterindex[ilayer][cindex];
      Float_t *weight = GetWeight(tindex);
      if (!weight){
	printf("Problem 1\n");  // not existing track
	continue;
      }
      if (*weight<(clusterweight[ilayer][cindex]-0.00001)){
	printf("Problem 2\n");  // not normalized probability
	continue;
      }   
      AliITSclusterV2 *c = (AliITSclusterV2*)GetCluster(tindex);
      if (c->GetNy()<5){
      	clusterweight[ilayer][cindex]/= *weight; 
      }
      else{
	Float_t weight2 = TMath::Max(*weight-0.5*clusterweight[ilayer][cindex],0.0000001);
	clusterweight[ilayer][cindex]/= weight2;
	if (	clusterweight[ilayer][cindex]>1)   clusterweight[ilayer][cindex] =1.;
      }
    }
  //
  //take to the account sharing factor
  //
  Float_t chi2 =10000000;
  Float_t sharefactor=0;
  Float_t minchi2 = 100000000;
  Float_t secchi2 = 100000000;
  Float_t norm=0;
  for (Int_t itrack=0;itrack<entries; itrack++){
    AliITStrackV2 * track = (AliITStrackV2*)array->At(itrack);
    if (!track) continue;
    if (track->GetChi2MIP(1)>1000) continue;  //not accepted - before
    chi2 = track->GetChi2MIP(1);
    Float_t newchi2=0;
    sharefactor=0;
    norm =0;
    //
    for (Int_t icluster=0;icluster<track->GetNumberOfClusters();icluster++){     
      Int_t tindex = track->GetClusterIndex(icluster);
      Int_t ilayer = (tindex & 0xf0000000) >> 28;
      if (tindex<0) continue;
      Int_t cindex =0;
      Float_t cchi2 = (track->GetDy(ilayer)*track->GetDy(ilayer))/(track->GetSigmaY(ilayer)*track->GetSigmaY(ilayer)) +
	(track->GetDz(ilayer)*track->GetDz(ilayer))/(track->GetSigmaZ(ilayer)*track->GetSigmaZ(ilayer)) ;
      //
      for (cindex=0;cindex<100;cindex++){
	if (clusterindex[ilayer][cindex]<0) break;
	if (clusterindex[ilayer][cindex]==tindex) break;
      }
      if (cindex>100) continue;
      if (clusterweight[ilayer][cindex]>0.00001){
	sharefactor+= clusterweight[ilayer][cindex];
	cchi2/=clusterweight[ilayer][cindex];
	norm++;
      }
      newchi2+=cchi2;
    }
    newchi2/=(norm-track->GetNSkipped()-track->GetNUsed());
    track->SetChi2MIP(4,newchi2);
    //if (norm>0) sharefactor/=norm;
    //if (sharefactor<0.5) return 0;
    //    chi2*=1./(0.5+sharefactor);
    if (newchi2<minchi2){
      besttrack = track;
      minchi2   = newchi2;
    }
    else{
      if (newchi2<secchi2){
	secchi2 = newchi2;
      }
    }
    //
  }

  //
  //
  if (!besttrack || besttrack->GetNumberOfClusters()<4) {
    return 0;
  }
  //
  if ((minchi2/secchi2)<0.7){
    //
    //increase the weight for clusters if the probability of other hypothesys is small
    Double_t deltad = besttrack->GetD(GetX(),GetY());
    Double_t deltaz = besttrack->GetZat(GetX()) - GetZ();
    Double_t deltaprim = TMath::Sqrt(deltad*deltad+deltaz*deltaz);
    //
    for (Int_t icluster=0;icluster<besttrack->GetNumberOfClusters();icluster++){
      Int_t index = besttrack->GetClusterIndex(icluster);
      Int_t ilayer =  (index & 0xf0000000) >> 28;
      AliITSclusterV2 *c = (AliITSclusterV2*)GetCluster(index);
      if (!c) continue; 
      if (c->GetNy()>3) continue;
      if ( (ilayer*0.2+0.2)<deltaprim) continue; // secondaries     
      Float_t * weight = GetWeight(index);
      *weight*=2.;
      *weight+=1.;
    }   
  }
  

  besttrack->SetLabel(original->GetLabel());
  CookLabel(besttrack,0.0);
  
  //
  return besttrack;
} 

