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
//-------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <iostream.h>

#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "../TPC/AliTPCtrack.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerV2.h"

//#define DEBUG

#ifdef DEBUG
Int_t LAB=236;
#endif

extern TRandom *gRandom;

AliITStrackerV2::AliITSlayer AliITStrackerV2::fLayers[kMaxLayer]; // ITS layers

AliITStrackerV2::
AliITStrackerV2(const AliITSgeom *geom) throw (const Char_t *) {
  //--------------------------------------------------------------------
  //This is the AliITStracker constructor
  //It also reads clusters from a root file
  //--------------------------------------------------------------------
  fYV=fZV=0.;

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
        phi+=0.5*TMath::Pi(); 
        AliITSdetector &det=fLayers[i-1].GetDetector((j-1)*ndet + k-1); 

if (phi<0) phi += 2*TMath::Pi();

        new(&det) AliITSdetector(r,phi); 
      } 
    }  

  }

  try {
     //Read clusters
     TTree *cTree=(TTree*)gDirectory->Get("cTree");
     if (!cTree) throw 
        ("AliITStrackerV2::AliITStrackerV2 can't get cTree !\n");

     TBranch *branch=cTree->GetBranch("Clusters");
     if (!branch) throw 
        ("AliITStrackerV2::AliITStrackerV2 can't get Clusters branch !\n");

     TClonesArray dummy("AliITSclusterV2",10000), *clusters=&dummy;
     branch->SetAddress(&clusters);

     Int_t nentr=(Int_t)cTree->GetEntries();
     for (i=0; i<nentr; i++) {
       if (!cTree->GetEvent(i)) continue;
       Int_t lay,lad,det; g->GetModuleId(i-1,lay,lad,det);
       Int_t ncl=clusters->GetEntriesFast();
       while (ncl--) {
         AliITSclusterV2 *c=(AliITSclusterV2*)clusters->UncheckedAt(ncl);

#ifdef DEBUG
if (c->GetLabel(2)!=LAB)
if (c->GetLabel(1)!=LAB)
if (c->GetLabel(0)!=LAB) continue;
cout<<lay-1<<' '<<lad-1<<' '<<det-1<<' '<<c->GetY()<<' '<<c->GetZ()<<endl;
#endif

         fLayers[lay-1].InsertCluster(new AliITSclusterV2(*c));
       }
       clusters->Delete();
     }
  }
  catch (const Char_t *msg) {
    cerr<<msg<<endl;
    throw;
  }

  fI=kMaxLayer;
}

static Int_t lbl;

Int_t AliITStrackerV2::Clusters2Tracks(const TFile *inp, TFile *out) {
  //--------------------------------------------------------------------
  //This functions reconstructs ITS tracks
  //--------------------------------------------------------------------
  TFile *in=(TFile*)inp;
  TDirectory *savedir=gDirectory; 

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
  TTree *tpcTree=(TTree*)in->Get("TreeT");
  if (!tpcTree) {
     cerr<<"AliITStrackerV2::Clusters2Tracks() ";
     cerr<<"can't get a tree with TPC tracks !\n";
     return 3;
  }

  AliTPCtrack *itrack=new AliTPCtrack; 
  tpcTree->SetBranchAddress("tracks",&itrack);

  out->cd();
  TTree itsTree("TreeT","Tree with ITS tracks");
  AliITStrackV2 *otrack=&fBestTrack;
  itsTree.Branch("tracks","AliITStrackV2",&otrack,32000,0);

  Int_t ntrk=0;

  Int_t nentr=(Int_t)tpcTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {

    if (!tpcTree->GetEvent(i)) continue;
    Int_t tpcLabel=itrack->GetLabel(); //save the TPC track label
lbl=tpcLabel;

#ifdef DEBUG
if (TMath::Abs(tpcLabel)!=LAB) continue;
cout<<tpcLabel<<" *****************\n";
#endif

    try {
      ResetTrackToFollow(AliITStrackV2(*itrack));
    } catch (const Char_t *msg) {
      cerr<<msg<<endl;
      continue;
    }
    ResetBestTrack();

    fYV=gRandom->Gaus(0.,kSigmaYV); fZV=gRandom->Gaus(0.,kSigmaZV);

    Double_t r2=fTrackToFollow.GetX()*fTrackToFollow.GetX();
    Double_t x0=0.082/21.82*2.33*(45*45+40*40)/r2+2.000/41.28*0.03*75*75/r2;
    fTrackToFollow.Improve(x0,fYV,fZV);

    //Double_t xk=77.2;
    Double_t xk=88.;
    fTrackToFollow.PropagateTo(xk,0.,0.); //Ne if it's still there

    xk-=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/44.77*1.71, 0.005*1.71); //Tedlar
    xk-=0.02;
    fTrackToFollow.PropagateTo(xk, 0.02/44.86*1.45, 0.02*1.45);   //Kevlar
       xk-=2.0;
       fTrackToFollow.PropagateTo(xk, 2.0/41.28*0.029, 2.0*0.029);//Nomex
    xk-=0.02;
    fTrackToFollow.PropagateTo(xk, 0.02/44.86*1.45, 0.02*1.45);   //Kevlar
    xk-=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/44.77*1.71, 0.005*1.71); //Tedlar

    xk-=14.5;
    fTrackToFollow.PropagateTo(xk,0.,0.); //C02

    xk -=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/24.01*2.7, 0.005*2.7);    //Al    
    xk -=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/44.77*1.71, 0.005*1.71);  //Tedlar
    xk -=0.02;
    fTrackToFollow.PropagateTo(xk, 0.02/44.86*1.45, 0.02*1.45);    //Kevlar
       xk -=0.5;
       fTrackToFollow.PropagateTo(xk, 0.5/41.28*.029, 0.5*0.029);  //Nomex    
    xk -=0.02;
    fTrackToFollow.PropagateTo(xk, 0.02/44.86*1.45, 0.02*1.45);    //Kevlar
    xk -=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/44.77*1.71, 0.005*1.71);  //Tedlar
    xk -=0.005;
    fTrackToFollow.PropagateTo(xk, 0.005/24.01*2.7, 0.005*2.7);    //Al    


    for (FollowProlongation(); fI<kMaxLayer; fI++) {
       while (TakeNextProlongation()) FollowProlongation();
    }

#ifdef DEBUG
cout<<fBestTrack.GetNumberOfClusters()<<" number of clusters\n\n";
#endif

    if (fBestTrack.GetNumberOfClusters() >= kMaxLayer-kLayersToSkip) {
       cerr << ++ntrk << "                \r";
       fBestTrack.SetLabel(tpcLabel);
       UseClusters(&fBestTrack);
       itsTree.Fill();
    }

  }

  itsTree.Write();
  savedir->cd();
  cerr<<"Number of TPC tracks: "<<nentr<<endl;
  cerr<<"Number of prolonged tracks: "<<ntrk<<endl;

  delete itrack;

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
    fI--;

#ifdef DEBUG
cout<<fI<<' ';
#endif

    AliITSlayer &layer=fLayers[fI];
    AliITStrackV2 &track=fTracks[fI];

    if (fI==3 || fI==1) {
      Double_t rs=0.5*(fLayers[fI+1].GetR() + layer.GetR());
      Double_t ds=0.034; if (fI==3) ds=0.039;
      Double_t dx0r=ds/21.82*2.33, dr=ds*2.33;
      fTrackToFollow.Propagate(fTrackToFollow.GetAlpha(),rs,1*dx0r,dr);
    }

    //find intersection
    Double_t x,y,z;  
    if (!fTrackToFollow.GetGlobalXYZat(layer.GetR(),x,y,z)) {
     cerr<<"AliITStrackerV2::FollowProlongation: failed to estimate track !\n";
     break;
    }
    Double_t phi=TMath::ATan2(y,x);
    Double_t d=layer.GetThickness(phi,z);
    Int_t idet=layer.FindDetectorIndex(phi,z);
    if (idet<0) {
    cerr<<"AliITStrackerV2::FollowProlongation: failed to find a detector !\n";
      break;
    }

    //propagate to the intersection
    const AliITSdetector &det=layer.GetDetector(idet);
    if (!fTrackToFollow.Propagate(det.GetPhi(),det.GetR(),1*d/21.82*2.33,d*2.33))
    {
      cerr<<"AliITStrackerV2::FollowProlongation: propagation failed !\n";
      break;
    }
    fTrackToFollow.SetDetectorIndex(idet);

    //Select possible prolongations and store the current track estimation
    track.~AliITStrackV2(); new(&track) AliITStrackV2(fTrackToFollow);
    Double_t dz=3*TMath::Sqrt(track.GetSigmaZ2() + kSigmaZ2[fI]);
    if (dz<0.5*TMath::Abs(track.GetTgl())) dz=0.5*TMath::Abs(track.GetTgl());
    if (dz > kMaxRoad/4) {
      //cerr<<"AliITStrackerV2::FollowProlongation: too broad road in Z !\n";
      break;
    }
    Double_t dy=4*TMath::Sqrt(track.GetSigmaY2() + kSigmaY2[fI]);
    if (dy<0.5*TMath::Abs(track.GetSnp())) dy=0.5*TMath::Abs(track.GetSnp());
    if (dy > kMaxRoad/4) {
      //cerr<<"AliITStrackerV2::FollowProlongation: too broad road in Y !\n";
      break;
    }

    Double_t zmin=track.GetZ() - dz;
    Double_t zmax=track.GetZ() + dz;
    Double_t ymin=track.GetY() + layer.GetR()*det.GetPhi() - dy;
    Double_t ymax=track.GetY() + layer.GetR()*det.GetPhi() + dy;
    if (ymax>layer.GetR()*2*TMath::Pi()) {
       ymax-=layer.GetR()*2*TMath::Pi();
       ymin-=layer.GetR()*2*TMath::Pi();
    }
    layer.SelectClusters(zmin,zmax,ymin,ymax);

    //if (1/TMath::Abs(track.Get1Pt())<0.2)
    //cout<<fI<<' '<<1/TMath::Abs(track.Get1Pt())<<' '
    //    <<dy<<' '<<dz<<' '<<layer.InRoad()<<endl;

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

  if (fI) fI++;
}


Int_t AliITStrackerV2::TakeNextProlongation() {
  //--------------------------------------------------------------------
  //This function takes another track prolongation 
  //--------------------------------------------------------------------
  //Double_t m[20];
  Double_t d=GetEffectiveThickness(0,0); //Think of this !!!!

  AliITSlayer &layer=fLayers[fI];
  AliITStrackV2 &t=fTracks[fI];

  Double_t dz=4*TMath::Sqrt(t.GetSigmaZ2() + kSigmaZ2[fI]);
  Double_t dy=4*TMath::Sqrt(t.GetSigmaY2() + kSigmaY2[fI]);

  const AliITSclusterV2 *c=0; Int_t ci=-1;
  Double_t chi2=12345.;
  while ((c=layer.GetNextCluster(ci))!=0) {
    //if (fI==0)
    //if (c->GetLabel(0)!=TMath::Abs(lbl)) continue; //88888888888888888888
    Int_t idet=c->GetDetectorIndex();

    if (t.GetDetectorIndex()!=idet) {
       const AliITSdetector &det=layer.GetDetector(idet);
       if (!t.Propagate(det.GetPhi(),det.GetR(),0.,0.)) {
         cerr<<"AliITStrackerV2::TakeNextProlongation: propagation failed !\n";
         continue;
       }
       t.SetDetectorIndex(idet);

#ifdef DEBUG
cout<<fI<<" change detector !\n";
#endif

    }

    if (TMath::Abs(t.GetZ() - c->GetZ()) > dz) continue;
    if (TMath::Abs(t.GetY() - c->GetY()) > dy) continue;

    //m[0]=fYV; m[1]=fZV;
    //chi2=t.GetPredictedChi2(c,m,d/21.82*2.33);
    chi2=t.GetPredictedChi2(c);

    if (chi2<kMaxChi2) break;
  }

#ifdef DEBUG
cout<<fI<<" chi2="<<chi2<<' '<<t.GetY()<<' '<<t.GetZ()<<' '<<dy<<' '<<dz<<endl;
#endif

  if (chi2>=kMaxChi2) return 0;
  if (!c) return 0;

  ResetTrackToFollow(t);

  //if (!fTrackToFollow.Update(m,chi2,(fI<<28)+ci)) {
  if (!fTrackToFollow.Update(c,chi2,(fI<<28)+ci)) {
     cerr<<"AliITStrackerV2::TakeNextProlongation: filtering failed !\n";
     return 0;
  }
  fTrackToFollow.Improve(d/21.82*2.33,fYV,fZV);

#ifdef DEBUG
cout<<"accepted lab="<<c->GetLabel(0)<<' '
    <<fTrackToFollow.GetNumberOfClusters()<<' '
    <<fTrackToFollow.GetY()<<' '<<fTrackToFollow.GetZ()<<endl<<endl;
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
}

AliITStrackerV2::AliITSlayer::~AliITSlayer() {
  //--------------------------------------------------------------------
  // AliITSlayer destructor
  //--------------------------------------------------------------------
  delete[] fDetectors;
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
}

Int_t AliITStrackerV2::AliITSlayer::InsertCluster(AliITSclusterV2 *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to this layer
  //--------------------------------------------------------------------
  if (fN==kMaxClusterPerLayer) {
 cerr<<"AliITStrackerV2::AliITSlayer::InsertCluster(): Too many clusters !\n"; 
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

  Double_t dz=fZOffset-z;
  Int_t nz=Int_t(dz*(fNdetectors-1)*0.5/fZOffset+0.5);

  if (np>=fNladders) np-=fNladders;
  if (np<0)          np+=fNladders;

#ifdef DEBUG
cout<<np<<' '<<nz<<endl;
#endif

  return np*fNdetectors + nz;
}

Double_t 
AliITStrackerV2::AliITSlayer::GetThickness(Double_t phi, Double_t z) const {
  //--------------------------------------------------------------------
  //This function returns the thickness of this layer
  //--------------------------------------------------------------------
  //-pi<phi<+pi
  if (3 <fR&&fR<8 ) return 1.1*0.096;
  if (13<fR&&fR<26) return 1.1*0.088;
  if (37<fR&&fR<41) return 1.1*0.085;
  return 1.1*0.081;
}


Double_t AliITStrackerV2::GetEffectiveThickness(Double_t phi,Double_t z) const
{
  //--------------------------------------------------------------------
  //Returns the thickness between the current layer and the vertex
  //--------------------------------------------------------------------
  //-pi<phi<+pi
  Double_t d=0.1;

  Double_t xn=fLayers[fI].GetR();
  for (Int_t i=0; i<fI; i++) {
    Double_t xi=fLayers[i].GetR();
    d+=fLayers[i].GetThickness(phi,z)*xi*xi;
  }

  if (fI>1) {
    Double_t xi=0.5*(fLayers[1].GetR()+fLayers[2].GetR());
    d+=0.034*xi*xi;
  }

  if (fI>3) {
    Double_t xi=0.5*(fLayers[3].GetR()+fLayers[4].GetR());
    d+=0.039*xi*xi;
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
    //if (c->IsUsed()) continue;
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



