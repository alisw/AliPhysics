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
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------

#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <Riostream.h>

#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliTPCtrack.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerV2.h"

//#define DEBUG

#ifdef DEBUG
Int_t LAB=38;
#endif

AliITStrackerV2::AliITSlayer AliITStrackerV2::fLayers[kMaxLayer]; // ITS layers

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

    new (fLayers+i-1) AliITSlayer(r,poff,zoff,nlad,ndet);

    for (Int_t j=1; j<nlad+1; j++) {
      for (Int_t k=1; k<ndet+1; k++) { //Fill this layer with detectors
        Float_t x,y,zshift; g->GetTrans(i,j,k,x,y,zshift); 
        Double_t rot[9]; g->GetRotMatrix(i,j,k,rot);

        Double_t r =-x*rot[1] + y*rot[0];         if (i==1) r=-r;
        Double_t phi=TMath::ATan2(rot[1],rot[0]); if (i==1) phi-=3.1415927;
        phi+=0.5*TMath::Pi(); if (phi<0) phi += 2*TMath::Pi();
        AliITSdetector &det=fLayers[i-1].GetDetector((j-1)*ndet + k-1); 

        new(&det) AliITSdetector(r,phi); 
      } 
    }  

  }

  fI=kMaxLayer;

  fPass=0;
  fConstraint[0]=1; fConstraint[1]=0;
}

void AliITStrackerV2::LoadClusters() {
  //--------------------------------------------------------------------
  //This function loads ITS clusters
  //--------------------------------------------------------------------
  char   cname[100]; 
  sprintf(cname,"TreeC_ITS_%d",GetEventNumber());
  TTree *cTree=(TTree*)gDirectory->Get(cname);
  if (!cTree) { 
    cerr<<"AliITStrackerV2::LoadClusters can't get cTree !\n";
    exit(1);
  }
  TBranch *branch=cTree->GetBranch("Clusters");
  if (!branch) { 
    cerr<<"AliITStrackerV2::LoadClusters can't get Clusters branch !\n";
    exit(1);
  }

  TClonesArray dummy("AliITSclusterV2",10000), *clusters=&dummy;
  branch->SetAddress(&clusters);

  Int_t j=0;
  for (Int_t i=0; i<kMaxLayer; i++) {
    Int_t ndet=fLayers[i].GetNdetectors();
    Int_t jmax = j + fLayers[i].GetNladders()*ndet;
    for (; j<jmax; j++) {           
      if (!cTree->GetEvent(j)) continue;
      Int_t ncl=clusters->GetEntriesFast();
      while (ncl--) {
        AliITSclusterV2 *c=(AliITSclusterV2*)clusters->UncheckedAt(ncl);

#ifdef DEBUG
if (c->GetLabel(2)!=LAB)
if (c->GetLabel(1)!=LAB)
if (c->GetLabel(0)!=LAB) continue;
Int_t idx=c->GetDetectorIndex();
Int_t lay=i, lad=idx/ndet, det=idx-lad*ndet;
cout<<lay<<' '<<lad<<' '<<det<<' '<<c->GetY()<<' '<<c->GetZ()<<endl;
#endif

        fLayers[i].InsertCluster(new AliITSclusterV2(*c));
      }
      clusters->Delete();
    }
    fLayers[i].ResetRoad(); //road defined by the cluster density
  }
  delete cTree; //Thanks to Mariana Bondila
}

void AliITStrackerV2::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer; i++) fLayers[i].ResetClusters();
}

//#ifdef DEBUG
static Int_t lbl;
//#endif

Int_t AliITStrackerV2::Clusters2Tracks(const TFile *inp, TFile *out) {
  //--------------------------------------------------------------------
  //This functions reconstructs ITS tracks
  //--------------------------------------------------------------------
  TFile *in=(TFile*)inp;
  TDirectory *savedir=gDirectory; 

  LoadClusters();

  if (!in->IsOpen()) {
     cerr<<"AliITStrackerV2::Clusters2Tracks(): ";
     cerr<<"file with TPC tracks is not open !\n";
     return 1;
  }

  if (!out->IsOpen()) {
     cerr<<"AliITStrackerV2::Clusters2Tracks(): ";
     cerr<<"file for ITS tracks is not open !\n";
     return 2;
  }

  in->cd();
 
  Char_t tname[100];
  Int_t nentr=0; TObjArray itsTracks(15000);

  {/* Read TPC tracks */ 
    sprintf(tname,"TreeT_TPC_%d",GetEventNumber());
    TTree *tpcTree=(TTree*)in->Get(tname);
    if (!tpcTree) {
       cerr<<"AliITStrackerV2::Clusters2Tracks(): "
             "can't get a tree with TPC tracks !\n";
       return 3;
    }
    AliTPCtrack *itrack=new AliTPCtrack; 
    tpcTree->SetBranchAddress("tracks",&itrack);
    nentr=(Int_t)tpcTree->GetEntries();
    for (Int_t i=0; i<nentr; i++) {
       tpcTree->GetEvent(i);
       AliITStrackV2 *t=0;
       try {
           t=new AliITStrackV2(*itrack);
       } catch (const Char_t *msg) {
           cerr<<msg<<endl;
           delete t;
           continue;
       }
       if (TMath::Abs(t->GetD())>4) continue;

       t->PropagateTo(80.,0.0053,30);
       if (TMath::Abs(t->GetY())>12.8) t->CorrectForMaterial(0.03);
       if (TMath::Abs(t->GetZ())<0.2) t->CorrectForMaterial(0.40);
       t->PropagateTo(61.,0.0053,30);
       //Double_t xk=52.,x,y,z; t->GetGlobalXYZat(xk,x,y,z);
       //if (TMath::Abs(y)<7.77) t->PropagateTo(xk,0.19,24.); 
       t->PropagateTo(50.,0.001);

       itsTracks.AddLast(t);
    }
    delete tpcTree; //Thanks to Mariana Bondila
    delete itrack;
  }
  itsTracks.Sort();

  out->cd();

  sprintf(tname,"TreeT_ITS_%d",GetEventNumber());
  TTree itsTree(tname,"Tree with ITS tracks");
  AliITStrackV2 *otrack=&fBestTrack;
  itsTree.Branch("tracks","AliITStrackV2",&otrack,32000,0);

  for (fPass=0; fPass<2; fPass++) {
     Int_t &constraint=fConstraint[fPass]; if (constraint<0) continue;
     for (Int_t i=0; i<nentr; i++) {
       if (i%10==0) cerr<<nentr-i<<" \r";
       AliITStrackV2 *t=(AliITStrackV2*)itsTracks.UncheckedAt(i);
       if (t==0) continue;           //this track has been already tracked
       Int_t tpcLabel=t->GetLabel(); //save the TPC track label

lbl=tpcLabel;
#ifdef DEBUG
lbl=tpcLabel;
if (TMath::Abs(tpcLabel)!=LAB) continue;
cout<<tpcLabel<<" *****************\n";
#endif

       ResetTrackToFollow(*t);
       ResetBestTrack();

       for (FollowProlongation(); fI<kMaxLayer; fI++) {
          while (TakeNextProlongation()) FollowProlongation();
       }

#ifdef DEBUG
cout<<fBestTrack.GetNumberOfClusters()<<" number of clusters\n\n";
#endif

      if (fBestTrack.GetNumberOfClusters() < kMaxLayer-kLayersToSkip) continue;

        if (fConstraint[fPass]) {
	   Int_t index[kMaxLayer];
           Int_t k;
           for (k=0; k<kMaxLayer; k++) index[k]=-1;
           Int_t nc=fBestTrack.GetNumberOfClusters();
           for (k=0; k<nc; k++) { 
              Int_t idx=fBestTrack.GetClusterIndex(k),nl=(idx&0xf0000000)>>28;
              index[nl]=idx; 
           }
           fBestTrack.~AliITStrackV2(); new(&fBestTrack) AliITStrackV2(*t);
	   if (!RefitAt(3.7, &fBestTrack, index)) continue;
	}

        fBestTrack.SetLabel(tpcLabel);
	fBestTrack.CookdEdx();
        CookLabel(&fBestTrack,0.); //For comparison only
        itsTree.Fill();
        UseClusters(&fBestTrack);
        delete itsTracks.RemoveAt(i);

     }
  }

  itsTree.Write();

  itsTracks.Delete();

  UnloadClusters();

  savedir->cd();
  cerr<<"Number of TPC tracks: "<<nentr<<endl;
  cerr<<"Number of prolonged tracks: "<<itsTree.GetEntries()<<endl;

  return 0;
}

Int_t AliITStrackerV2::PropagateBack(const TFile *inp, TFile *out) {
  //--------------------------------------------------------------------
  //This functions propagates reconstructed ITS tracks back
  //--------------------------------------------------------------------
  TFile *in=(TFile*)inp;
  TDirectory *savedir=gDirectory; 

  LoadClusters();

  if (!in->IsOpen()) {
     cerr<<"AliITStrackerV2::PropagateBack(): ";
     cerr<<"file with ITS tracks is not open !\n";
     return 1;
  }

  if (!out->IsOpen()) {
     cerr<<"AliITStrackerV2::PropagateBack(): ";
     cerr<<"file for back propagated ITS tracks is not open !\n";
     return 2;
  }

  in->cd();

  Char_t tname[100];
  sprintf(tname,"TreeT_ITS_%d",GetEventNumber());
  TTree *itsTree=(TTree*)in->Get(tname);
  if (!itsTree) {
     cerr<<"AliITStrackerV2::PropagateBack() ";
     cerr<<"can't get a tree with ITS tracks !\n";
     return 3;
  }
  AliITStrackV2 *itrack=new AliITStrackV2; 
  itsTree->SetBranchAddress("tracks",&itrack);

  out->cd();

  sprintf(tname,"TreeT_ITSb_%d",GetEventNumber());
  TTree backTree(tname,"Tree with back propagated ITS tracks");
  AliTPCtrack *otrack=0;
  backTree.Branch("tracks","AliTPCtrack",&otrack,32000,2);

  Int_t ntrk=0;

  Int_t nentr=(Int_t)itsTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
    itsTree->GetEvent(i);
    ResetTrackToFollow(*itrack);

    // propagete to vertex [SR, GSI 17.02.2003]
    fTrackToFollow.PropagateTo(3.,0.0028,65.19);
    fTrackToFollow.PropagateToVertex();

    // Start Time measurement [SR, GSI 17.02.2003]
    fTrackToFollow.StartTimeIntegral();
    fTrackToFollow.PropagateTo(3.,-0.0028,65.19);
    //

    fTrackToFollow.ResetCovariance(); fTrackToFollow.ResetClusters();

    Int_t itsLabel=fTrackToFollow.GetLabel(); //save the ITS track label

#ifdef DEBUG
if (TMath::Abs(itsLabel)!=LAB) continue;
cout<<itsLabel<<" *****************\n";
#endif

    try {
       Int_t nc=itrack->GetNumberOfClusters();
#ifdef DEBUG
for (Int_t k=0; k<nc; k++) {
    Int_t index=itrack->GetClusterIndex(k);
    AliITSclusterV2 *c=(AliITSclusterV2*)GetCluster(index);
    cout<<c->GetLabel(0)<<' '<<c->GetLabel(1)<<' '<<c->GetLabel(2)<<endl;
}
#endif       
       Int_t idx=0, l=0; 
       const  AliITSclusterV2 *c=0; 
       if (nc--) {
          idx=itrack->GetClusterIndex(nc); l=(idx&0xf0000000)>>28;
          c=(AliITSclusterV2*)GetCluster(idx);
       }
       for (fI=0; fI<kMaxLayer; fI++) {
         AliITSlayer &layer=fLayers[fI];
         Double_t r=layer.GetR();
         if (fI==2 || fI==4) {             
            Double_t rs=0.5*(fLayers[fI-1].GetR() + r);
            Double_t d=0.0034, x0=38.6;
            if (fI==2) {rs=9.; d=0.0097; x0=42.;}
            if (!fTrackToFollow.PropagateTo(rs,-d,x0)) throw "";
         }

         // remember old position [SR, GSI 18.02.2003]
	 Double_t oldX, oldY, oldZ;
	 fTrackToFollow.GetGlobalXYZat(fTrackToFollow.GetX(),oldX,oldY,oldZ);
	 //

         Double_t x,y,z;
         if (!fTrackToFollow.GetGlobalXYZat(r,x,y,z)) 
            throw "AliITStrackerV2::PropagateBack: failed to estimate track !";
         Double_t phi=TMath::ATan2(y,x);
         Int_t idet=layer.FindDetectorIndex(phi,z);
         if (idet<0) 
         throw "AliITStrackerV2::PropagateBack: failed to find a detector !\n";
         const AliITSdetector &det=layer.GetDetector(idet);
         r=det.GetR(); phi=det.GetPhi();
         if (!fTrackToFollow.Propagate(phi,r)) throw "";
         fTrackToFollow.SetDetectorIndex(idet);

         const AliITSclusterV2 *cl=0;
         Int_t index=0;
         Double_t maxchi2=kMaxChi2;

         if (l==fI) {
           idet=c->GetDetectorIndex();
           if (idet != fTrackToFollow.GetDetectorIndex()) {
             const AliITSdetector &det=layer.GetDetector(idet);
             r=det.GetR(); phi=det.GetPhi();
             if (!fTrackToFollow.Propagate(phi,r)) throw "";
             fTrackToFollow.SetDetectorIndex(idet);
           }
           Double_t chi2=fTrackToFollow.GetPredictedChi2(c);
           if (chi2<kMaxChi2) {
              cl=c; maxchi2=chi2; index=idx;
           }
           if (nc--) {
              idx=itrack->GetClusterIndex(nc); l=(idx&0xf0000000)>>28;
              c=(AliITSclusterV2*)GetCluster(idx);
           } 
         }

         if (fTrackToFollow.GetNumberOfClusters()>2) {
           Double_t dz=3*TMath::Sqrt(fTrackToFollow.GetSigmaZ2()+kSigmaZ2[fI]);
           Double_t dy=3*TMath::Sqrt(fTrackToFollow.GetSigmaY2()+kSigmaY2[fI]);
           Double_t zmin=fTrackToFollow.GetZ() - dz;
           Double_t zmax=fTrackToFollow.GetZ() + dz;
           Double_t ymin=fTrackToFollow.GetY() + phi*r - dy;
           Double_t ymax=fTrackToFollow.GetY() + phi*r + dy;
           layer.SelectClusters(zmin,zmax,ymin,ymax);

           const AliITSclusterV2 *cc=0; Int_t ci;
           while ((cc=layer.GetNextCluster(ci))!=0) {
              idet=cc->GetDetectorIndex();
              if (idet != fTrackToFollow.GetDetectorIndex()) continue;
              Double_t chi2=fTrackToFollow.GetPredictedChi2(cc);
              if (chi2<maxchi2) {
                 cl=cc; index=(fI<<28)+ci; maxchi2=chi2;
              }
           }
         }

         if (cl) {
            if (!fTrackToFollow.Update(cl,maxchi2,index)) 
              cerr<<"AliITStrackerV2::PropagateBack: filtering failed !\n";
         }
         {
          Double_t x0;
          x=layer.GetThickness(fTrackToFollow.GetY(),fTrackToFollow.GetZ(),x0);
          fTrackToFollow.CorrectForMaterial(-x,x0); 
         }
         	 
         // track time update [SR, GSI 17.02.2003]
         Double_t newX, newY, newZ;
	 fTrackToFollow.GetGlobalXYZat(fTrackToFollow.GetX(),newX,newY,newZ);
         Double_t dL2 = (oldX-newX)*(oldX-newX)+(oldY-newY)*(oldY-newY)+(oldZ-newZ)*(oldZ-newZ);
	 fTrackToFollow.AddTimeStep(TMath::Sqrt(dL2));
         //

       }

       fTrackToFollow.PropagateTo(50.,-0.001);
       //Double_t xk=52.,x,y,z; fTrackToFollow.GetGlobalXYZat(xk,x,y,z);
       //if (TMath::Abs(y)<7.77) fTrackToFollow.PropagateTo(xk,-0.19,24.); 
       fTrackToFollow.PropagateTo(61,-0.0053,30);
       fTrackToFollow.PropagateTo(80.,-0.0053,30);

       fTrackToFollow.SetLabel(itsLabel);
       otrack=new AliTPCtrack(fTrackToFollow,fTrackToFollow.GetAlpha()); 
       backTree.Fill(); delete otrack;
       UseClusters(&fTrackToFollow);
       cerr << ++ntrk << "                \r";
    }
    catch (const Char_t *msg) {
       cerr<<msg<<endl;
    }
  }

  backTree.Write();

  cerr<<"Number of ITS tracks: "<<nentr<<endl;
  cerr<<"Number of back propagated ITS tracks: "<<ntrk<<endl;

  delete itrack;

  delete itsTree; //Thanks to Mariana Bondila

  UnloadClusters();

  savedir->cd();

  return 0;
}

AliCluster *AliITStrackerV2::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 00;
  return fLayers[l].GetCluster(c);
}


void AliITStrackerV2::FollowProlongation() {
  //--------------------------------------------------------------------
  //This function finds a track prolongation 
  //--------------------------------------------------------------------
  Int_t tryAgain=kLayersToSkip;

  while (fI) {
    Int_t i=fI-1;
#ifdef DEBUG
cout<<i<<' ';
#endif
    AliITSlayer &layer=fLayers[i];
    AliITStrackV2 &track=fTracks[i];

    Double_t r=layer.GetR();

    if (i==3 || i==1) {
       Double_t rs=0.5*(fLayers[i+1].GetR() + r);
       Double_t d=0.0034, x0=38.6;
       if (i==1) {rs=9.; d=0.0097; x0=42;}
       if (!fTrackToFollow.PropagateTo(rs,d,x0)) {
	 //cerr<<"AliITStrackerV2::FollowProlongation: "
         //"propagation failed !\n";
         break;
       }
    }

    //find intersection
    Double_t x,y,z;  
    if (!fTrackToFollow.GetGlobalXYZat(r,x,y,z)) {
      //cerr<<"AliITStrackerV2::FollowProlongation: "
      //"failed to estimate track !\n";
      break;
    }
    Double_t phi=TMath::ATan2(y,x);
    Int_t idet=layer.FindDetectorIndex(phi,z);
    if (idet<0) {
      //cerr<<"AliITStrackerV2::FollowProlongation: "
      //"failed to find a detector !\n";
      break;
    }

    //propagate to the intersection
    const AliITSdetector &det=layer.GetDetector(idet);
    phi=det.GetPhi();
    if (!fTrackToFollow.Propagate(phi,det.GetR())) {
      //cerr<<"AliITStrackerV2::FollowProlongation: propagation failed !\n";
      break;
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
      //cerr<<"AliITStrackerV2::FollowProlongation: too broad road in Z !\n";
      break;
    }

    if (TMath::Abs(fTrackToFollow.GetZ()-GetZ()) > r+dz) break;

    //Double_t dy=4*TMath::Sqrt(track.GetSigmaY2() + kSigmaY2[i]);
    if (dy < 0.5*TMath::Abs(track.GetSnp())) dy=0.5*TMath::Abs(track.GetSnp());
    if (dy > kMaxRoad) {
      //cerr<<"AliITStrackerV2::FollowProlongation: too broad road in Y !\n";
      break;
    }

    Double_t zmin=track.GetZ() - dz; 
    Double_t zmax=track.GetZ() + dz;
    Double_t ymin=track.GetY() + r*phi - dy;
    Double_t ymax=track.GetY() + r*phi + dy;
    layer.SelectClusters(zmin,zmax,ymin,ymax); 
    fI--;

    //take another prolongation
    if (!TakeNextProlongation()) if (!tryAgain--) break;
    tryAgain=kLayersToSkip;

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
  AliITSlayer &layer=fLayers[fI];
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
    //if (c->GetLabel(0)!=TMath::Abs(lbl)) continue; 
    Int_t idet=c->GetDetectorIndex();

    if (fTrackToFollow.GetDetectorIndex()!=idet) {
       const AliITSdetector &det=layer.GetDetector(idet);
       ResetTrackToFollow(fTracks[fI]);
       if (!fTrackToFollow.Propagate(det.GetPhi(),det.GetR())) {
         //cerr<<"AliITStrackerV2::TakeNextProlongation: "
         //"propagation failed !\n";
         continue;
       }
       fTrackToFollow.SetDetectorIndex(idet);
       if (TMath::Abs(fTrackToFollow.GetZ()-GetZ())>layer.GetR()+dz) continue;

#ifdef DEBUG
cout<<fI<<" change detector !\n";
#endif

    }

    if (TMath::Abs(fTrackToFollow.GetZ() - c->GetZ()) > dz) continue;
    if (TMath::Abs(fTrackToFollow.GetY() - c->GetY()) > dy) continue;

    chi2=fTrackToFollow.GetPredictedChi2(c); if (chi2<kMaxChi2) break;
  }

#ifdef DEBUG
cout<<fI<<" chi2="<<chi2<<' '
    <<fTrackToFollow.GetY()<<' '<<fTrackToFollow.GetZ()<<' '
    <<dy<<' '<<dz<<endl;
{
Double_t phi=TMath::ATan(fTrackToFollow.GetY()/fTrackToFollow.GetX())
            +fTrackToFollow.GetAlpha(); phi=phi*180./3.141593;
cout<<phi<<endl;
}
#endif

  if (chi2>=kMaxChi2) return 0;
  if (!c) return 0;

  if (!fTrackToFollow.Update(c,chi2,(fI<<28)+ci)) {
     //cerr<<"AliITStrackerV2::TakeNextProlongation: filtering failed !\n";
     return 0;
  }

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
    fTrackToFollow.Improve(d,GetY(),GetZ());
  }

#ifdef DEBUG
cout<<"accepted lab="<<c->GetLabel(0)<<' '
    <<fTrackToFollow.GetNumberOfClusters()<<' '
    <<fTrackToFollow.GetY()<<' '<<fTrackToFollow.GetZ()<<' '
    <<fTrackToFollow.Get1Pt()<<endl<<endl;
#endif

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
}

void AliITStrackerV2::AliITSlayer::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  fN=0;
  fI=0;
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
     cerr<<"AliITStrackerV2::AliITSlayer::InsertCluster(): "
           "Too many clusters !\n";
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
    if (c->IsUsed()) continue;
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
  Double_t dphi=phi-fPhiOffset;
  if      (dphi <  0) dphi += 2*TMath::Pi();
  else if (dphi >= 2*TMath::Pi()) dphi -= 2*TMath::Pi();
  Int_t np=Int_t(dphi*fNladders*0.5/TMath::Pi()+0.5);
  if (np>=fNladders) np-=fNladders;
  if (np<0)          np+=fNladders;

  Double_t dz=fZOffset-z;
  Int_t nz=Int_t(dz*(fNdetectors-1)*0.5/fZOffset+0.5);
  if (nz>=fNdetectors) return -1;
  if (nz<0)            return -1;

#ifdef DEBUG
cout<<np<<' '<<nz<<endl;
#endif

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

#ifdef DEBUG
  cout<<"d="<<d<<endl;
#endif

  return d;
}

Double_t AliITStrackerV2::GetEffectiveThickness(Double_t y,Double_t z) const
{
  //--------------------------------------------------------------------
  //Returns the thickness between the current layer and the vertex (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0028*3*3; //beam pipe
  Double_t x0=0;

  Double_t xn=fLayers[fI].GetR();
  for (Int_t i=0; i<fI; i++) {
    Double_t xi=fLayers[i].GetR();
    d+=fLayers[i].GetThickness(y,z,x0)*xi*xi;
  }

  if (fI>1) {
    Double_t xi=9.;
    d+=0.0097*xi*xi;
  }

  if (fI>3) {
    Double_t xi=0.5*(fLayers[3].GetR()+fLayers[4].GetR());
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

Bool_t AliITStrackerV2::RefitAt(Double_t x, AliITStrackV2 *t, Int_t *index) {
  //--------------------------------------------------------------------
  // This function refits a track at a given position
  //--------------------------------------------------------------------
  Int_t from, to, step;
  if (x > t->GetX()) {
      from=0; to=kMaxLayer;
      step=+1;
  } else {
      from=kMaxLayer-1; to=-1;
      step=-1;
  }

  for (Int_t i=from; i != to; i += step) {
     AliITSlayer &layer=fLayers[i];
     Double_t r=layer.GetR();
 
     {
     Double_t hI=i-0.5*step; 
     if (hI==1.5 || hI==3.5) {             
        Double_t rs=0.5*(fLayers[i-step].GetR() + r);
        Double_t d=0.0034, x0=38.6; 
        if (hI==1.5) {rs=9.; d=0.0097; x0=42;}
        if (!t->PropagateTo(rs,d,x0)) {
          return kFALSE;
        }
     }
     }

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
        if (chi2<maxchi2) { cl=c; maxchi2=chi2; }
        else return kFALSE;
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
     }

     {
     Double_t x0;
     Double_t d=layer.GetThickness(t->GetY(),t->GetZ(),x0);
     t->CorrectForMaterial(-step*d,x0);
     }

  }

  if (!t->PropagateTo(x,0.,0.)) return kFALSE;
  return kTRUE;
}

void AliITStrackerV2::UseClusters(const AliKalmanTrack *t, Int_t from) const {
  //--------------------------------------------------------------------
  // This function marks clusters assigned to the track
  //--------------------------------------------------------------------
  AliTracker::UseClusters(t,from);
  /*
  AliITSclusterV2 *c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(0));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();
  c=(AliITSclusterV2 *)GetCluster(t->GetClusterIndex(1));
  //if (c->GetQ()>2) c->Use();
  if (c->GetSigmaZ2()>0.1) c->Use();
  */
}
