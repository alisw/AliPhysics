//-------------------------------------------------------------------------
//               Implementation of the ITS tracker class
//    The pattern recongintion based on the "cooked covariance" approach
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSUTrackerCooked.h"
#include "AliITSUTrackCooked.h" 

using std::cout;
using std::endl;

ClassImp(AliITSUTrackerCooked)

//************************************************
// Constants hardcoded for the moment:
//************************************************
// radial positions of layers: default contructor
const 
Double_t klRadius[7]={2.34, 3.15, 3.93, 19.61, 24.55, 34.39, 39.34}; //tdr6
// seed "windows" in z and phi: MakeSeeds
const Double_t kzWin=0.33, kpWin=3.14/4;
// Maximal accepted impact parameters for the seeds 
const Double_t kmaxDCAxy=3.;
const Double_t kmaxDCAz= 3.;
// Layers for the seeding
const Int_t kSeedingLayer1=6, kSeedingLayer2=4, kSeedingLayer3=5;
// Space point resolution
const Double_t kSigma2=0.0005*0.0005;
// Max accepted chi2 per cluster
const Double_t kmaxChi2PerCluster=77.;

//************************************************
// TODO:
//************************************************
// Seeding:
// Precalculate cylidnrical (r,phi) for the clusters;
// use exact r's for the clusters


AliITSUTrackerCooked::AliITSUlayer
              AliITSUTrackerCooked::fgLayers[AliITSUTrackerCooked::kNLayers];

AliITSUTrackerCooked::AliITSUTrackerCooked(): 
AliTracker(),
fSeeds(0),
fI(kNLayers-1),
fBestTrack(0), 
fTrackToFollow(0) 
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------

  AliITSUGeomTGeo *gm  = new AliITSUGeomTGeo(kTRUE,kTRUE);
  AliITSUClusterPix::SetGeom(gm);

  for (Int_t i=0; i<kNLayers; i++) fgLayers[i].SetR(klRadius[i]);

  // Some default primary vertex
  Double_t xyz[]={0.,0.,0.};
  Double_t ers[]={2.,2.,2.};

  SetVertex(xyz,ers);

}

void AliITSUTrackerCooked::ResetTrackToFollow(const AliITSUTrackCooked &t) {
  //--------------------------------------------------------------------
  // Prepare to follow a new track seed
  //--------------------------------------------------------------------
     delete fTrackToFollow;
     fTrackToFollow = new AliITSUTrackCooked(t);
}
  
void AliITSUTrackerCooked::ResetBestTrack() {
  //--------------------------------------------------------------------
  // Replace the best track branch
  //--------------------------------------------------------------------
     delete fBestTrack;
     fBestTrack = new AliITSUTrackCooked(*fTrackToFollow);
}
  
AliITSUTrackerCooked::~AliITSUTrackerCooked() 
{
  //--------------------------------------------------------------------
  // Virtual destructor
  //--------------------------------------------------------------------

  if (fSeeds) fSeeds->Delete(); delete fSeeds; 
  delete fBestTrack;
  delete fTrackToFollow;

}

static Double_t 
f1(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3)
{
    //-----------------------------------------------------------------
    // Initial approximation of the track curvature
    //-----------------------------------------------------------------
    Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
    Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                    (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
    Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                    (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));
    
    Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
    
    return -xr*yr/sqrt(xr*xr+yr*yr);
}

static Double_t 
f2(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3)
{
    //-----------------------------------------------------------------
    // Initial approximation of the track curvature times center of curvature
    //-----------------------------------------------------------------
    Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
    Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                    (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
    Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                    (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));
    
    Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
    
    return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

static Double_t 
f3(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t z1, Double_t z2)
{
    //-----------------------------------------------------------------
    // Initial approximation of the tangent of the track dip angle
    //-----------------------------------------------------------------
    return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

Bool_t AliITSUTrackerCooked::
AddCookedSeed(const Float_t r1[3], Int_t l1, Int_t i1, 
              const Float_t r2[3], Int_t l2, Int_t i2,
              const AliCluster *c3,Int_t l3, Int_t i3) 
{
    //--------------------------------------------------------------------
    // This is the main cooking function.
    // Creates seed parameters out of provided clusters.
    //--------------------------------------------------------------------
    Float_t x,a;
    if (!c3->GetXAlphaRefPlane(x,a)) return kFALSE;

    Double_t ca=TMath::Cos(a), sa=TMath::Sin(a);
    Double_t x1 = r1[0]*ca + r1[1]*sa,
             y1 =-r1[0]*sa + r1[1]*ca, z1 = r1[2];
    Double_t x2 = r2[0]*ca + r2[1]*sa,
             y2 =-r2[0]*sa + r2[1]*ca, z2 = r2[2];
    Double_t x3 = x,  y3 = c3->GetY(), z3 = c3->GetZ();

    Double_t par[5];
    par[0]=y3;
    par[1]=z3;
    Double_t crv=f1(x1, y1, x2, y2, x3, y3); //curvature
    Double_t cx0=f2(x1, y1, x2, y2, x3, y3); //curvature*x0
    Double_t tgl12=f3(x1, y1, x2, y2, z1, z2);
    Double_t tgl23=f3(x2, y2, x3, y3, z2, z3);

    Double_t sf=x*crv - cx0;
    if (TMath::Abs(sf) >= kAlmost1) return kFALSE;
    par[2]=sf;

    par[3]=0.5*(tgl12 + tgl23);
    Double_t bz=GetBz();
    par[4]=(TMath::Abs(bz) < kAlmost0Field) ? kAlmost0 : crv/(bz*kB2C);

    Double_t cov[15];
    for (Int_t i=0; i<15; i++) cov[i]=0.;
    cov[0] =kSigma2*10;
    cov[2] =kSigma2*10;
    cov[5] =0.007*0.007*10;   //FIXME all these lines
    cov[9] =0.007*0.007*10;
    cov[14]=0.1*0.1*10;

    AliITSUTrackCooked *seed=new AliITSUTrackCooked();
    seed->Set(Double_t(x), Double_t(a), par, cov);

    Float_t dz[2]; 
    seed->GetDZ(GetX(),GetY(),GetZ(),GetBz(),dz);
    if (TMath::Abs(dz[0]) > kmaxDCAxy) {delete seed; return kFALSE;} 
    if (TMath::Abs(dz[1]) > kmaxDCAz ) {delete seed; return kFALSE;} 

    seed->SetClusterIndex(l1,i1);
    seed->SetClusterIndex(l2,i2);
    seed->SetClusterIndex(l3,i3);

    fSeeds->AddLast(seed);

    return kTRUE;
}

Int_t AliITSUTrackerCooked::MakeSeeds() {
  //--------------------------------------------------------------------
  // This is the main pattern recongition function.
  // Creates seeds out of two clusters and another point.
  //--------------------------------------------------------------------
   if (fSeeds) {fSeeds->Delete(); delete fSeeds;}
   fSeeds=new TObjArray(77777);

   Double_t zv=GetZ();

   AliITSUlayer &layer1=fgLayers[kSeedingLayer1];
   AliITSUlayer &layer2=fgLayers[kSeedingLayer2];
   AliITSUlayer &layer3=fgLayers[kSeedingLayer3];
   Double_t r1=layer1.GetR();
   Double_t r2=layer2.GetR();
   Double_t r3=layer3.GetR();

   Int_t nClusters1=layer1.GetNumberOfClusters();
   Int_t nClusters2=layer2.GetNumberOfClusters();
   Int_t nClusters3=layer3.GetNumberOfClusters();
   for (Int_t n1=0; n1<nClusters1; n1++) {
     AliCluster *c1=layer1.GetCluster(n1);
     //
     //Int_t lab=c1->GetLabel(0);
     //
     Double_t z1=c1->GetZ();
     Float_t xyz1[3]; c1->GetGlobalXYZ(xyz1);
     Double_t phi1=TMath::ATan2(xyz1[1],xyz1[0]);
     Double_t zr2=zv + r2/r1*(z1-zv);
     Int_t start2=layer2.FindClusterIndex(zr2-kzWin);
     for (Int_t n2=start2; n2<nClusters2; n2++) {
         AliCluster *c2=layer2.GetCluster(n2);
         //
         //if (c2->GetLabel(0)!=lab) continue;
	 //
         Double_t z2=c2->GetZ();
         
         if (z2 > (zr2+kzWin)) break;  //check in Z
         Float_t xyz2[3]; c2->GetGlobalXYZ(xyz2);
         Double_t phi2=TMath::ATan2(xyz2[1],xyz2[0]);
         if (TMath::Abs(phi2-phi1) > kpWin) continue;  //check in Phi
 
         Double_t zr3=z1 + (r3-r1)/(r2-r1)*(z2-z1);
         Double_t d13=r1-r3, d32=r3-r2, d=r1-r2;
         Double_t phir3=d32/d*phi1 + d13/d*phi2;  // FIXME
         Int_t start3=layer3.FindClusterIndex(zr3-kzWin/2);
         for (Int_t n3=start3; n3<nClusters3; n3++) {
             AliCluster *c3=layer3.GetCluster(n3);
             //
             //if (c3->GetLabel(0)!=lab) continue;
             //
             Double_t z3=c3->GetZ();
         
             if (z3 > (zr3+kzWin/2)) break;  //check in Z
             Float_t xyz3[3]; c3->GetGlobalXYZ(xyz3);
             Double_t phi3=TMath::ATan2(xyz3[1],xyz3[0]);
             if (TMath::Abs(phir3-phi3) > kpWin/100) continue;  //check in Phi

             AliITSUClusterPix cc(*((AliITSUClusterPix*)c2));
             cc.GoToFrameTrk();
             AddCookedSeed(xyz1, kSeedingLayer1, n1,
                           xyz3, kSeedingLayer3, n3, 
                           &cc,  kSeedingLayer2, n2);

	 }
     }
   }
   fSeeds->Sort();
   return fSeeds->GetEntriesFast();
}

Int_t AliITSUTrackerCooked::Clusters2Tracks(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This is the main tracking function
  // The clusters must already be loaded
  //--------------------------------------------------------------------

  Int_t nSeeds=MakeSeeds();

  // Possibly, icrement the seeds with additional clusters (Kalman)

  // Possibly, (re)fit the found tracks 

  Int_t nClusters1=fgLayers[kSeedingLayer1].GetNumberOfClusters();
  Int_t nClusters2=fgLayers[kSeedingLayer2].GetNumberOfClusters();
  cout<<nClusters1<<' '<<nClusters2<<endl;

    Int_t ngood=0;
    for (Int_t s=0; s<nSeeds; s++) {
        const AliITSUTrackCooked *track=(AliITSUTrackCooked*)fSeeds->At(s);
        ResetTrackToFollow(*track);
        ResetBestTrack();
        fI=kSeedingLayer2;
        fgLayers[fI].ResetTrack(*track); 

        for (FollowProlongation(); fI<kSeedingLayer2; fI++) {
            while (TakeNextProlongation()) FollowProlongation();
        }

        if (fBestTrack->GetNumberOfClusters() < 3) continue;

        CookLabel(fBestTrack,0.);
        Int_t label=fBestTrack->GetLabel();
        if (label>0) {ngood++;
            //cout<<track->label<<endl;
        }
        AliESDtrack iotrack;
        iotrack.UpdateTrackParams(fBestTrack,AliESDtrack::kITSin);
        iotrack.SetLabel(label);
        event->AddTrack(&iotrack);
        UseClusters(fBestTrack);
    }
    cout<<"Good tracks "<<ngood<<endl;
    cout<<"Good tracks/seeds "<<Float_t(ngood)/nSeeds<<endl;

  if (fSeeds) {fSeeds->Delete(); delete fSeeds;}
  fSeeds=0;
    
  return 0;
}

void AliITSUTrackerCooked::FollowProlongation() {
  //--------------------------------------------------------------------
  // Push this track tree branch towards the primary vertex
  //--------------------------------------------------------------------
  while (fI) {
    fI--;
    AliITSUlayer &layer = fgLayers[fI]; //fI is the number of the next layer
    Double_t r=layer.GetR();

    //Find intersection (fTrackToFollow is still at the previous layer)
    Double_t phi,z;  
    if (!fTrackToFollow->GetPhiZat(r,phi,z)) {
      //Warning("FollowProlongation","failed to estimate track !\n");
      return;
    }

    //Calculate the road at the next layer.  FIXME
    Double_t dr = fgLayers[fI+1].GetR() - r;
    Double_t xx0 = (fI > 2) ? 0.008 : 0.003;  //rough layer thickness
    Double_t p = fTrackToFollow->P();
    Double_t m = fTrackToFollow->GetMass();
    Double_t beta2 = p*p/(p*p + m*m); 
    Double_t scat2 = dr*dr*0.14*0.14*xx0/(beta2*p*p); 
    Double_t dz=7*TMath::Sqrt((scat2 + kSigma2) + kSigma2);
    Double_t dy=dz;
    //if (TMath::Abs(fTrackToFollow.GetZ()-GetZ()) > r+dz) return;
    Double_t zMin = z - dz; 
    Double_t zMax = z + dz;
    Double_t phiMin = phi - dy/r;
    Double_t phiMax = phi + dy/r;
    if (layer.SelectClusters(zMin, zMax, phiMin, phiMax)==0) return;  

    if (!TakeNextProlongation()) return;

  } 

  //deal with the best track
  Int_t ncl=fTrackToFollow->GetNumberOfClusters();
  Int_t nclb=fBestTrack->GetNumberOfClusters();
  if (ncl)
  if (ncl >= nclb) {
     Double_t chi2=fTrackToFollow->GetChi2();
     if (chi2/ncl < kmaxChi2PerCluster) {        
        if (ncl > nclb || chi2 < fBestTrack->GetChi2()) {
	   ResetBestTrack();
        }
     }
  }

}

Int_t AliITSUTrackerCooked::TakeNextProlongation() {
  //--------------------------------------------------------------------
  // Switch to the next track tree branch
  //--------------------------------------------------------------------
  AliITSUlayer &layer=fgLayers[fI];

  const AliCluster *c=0; Int_t ci=-1;
  const AliCluster *cc=0; Int_t cci=-1;
  UShort_t volId=-1;
  Double_t x, alpha;
  Double_t z, dz, y, dy, chi2; 
  while ((c=layer.GetNextCluster(ci))!=0) {
    if (c->IsClusterUsed()) continue;
    Int_t id=c->GetVolumeId();
    if (id != volId) {
       volId=id;
       Float_t xr,ar; c->GetXAlphaRefPlane(xr, ar);
       x=xr; alpha=ar;
       const AliITSUTrackCooked *t = fgLayers[fI+1].GetTrack();
       ResetTrackToFollow(*t);
       if (!fTrackToFollow->Propagate(alpha, x, GetBz())) {
         //Warning("TakeNextProlongation","propagation failed !\n");
          continue;
       }
       dz=7*TMath::Sqrt(fTrackToFollow->GetSigmaZ2() + kSigma2);
       dy=7*TMath::Sqrt(fTrackToFollow->GetSigmaY2() + kSigma2);
       z=fTrackToFollow->GetZ();
       y=fTrackToFollow->GetY();
    }

    //if (TMath::Abs(fTrackToFollow.GetZ()-GetZ())>layer.GetR()+dz) continue;

    if (TMath::Abs(z - c->GetZ()) > dz) continue;
    if (TMath::Abs(y - c->GetY()) > dy) continue;

    Double_t ch2=fTrackToFollow->GetPredictedChi2(c); 
    if (ch2 > kmaxChi2PerCluster) continue;
    chi2=ch2;
    cc=c; cci=ci;
    break;
  }

  if (!cc) return 0;

  if (!fTrackToFollow->Update(cc,chi2,(fI<<28)+cci)) {
     //Warning("TakeNextProlongation","filtering failed !\n");
     return 0;
  }
  Double_t xx0 = (fI > 2) ? 0.008 : 0.003;  // Rough layer thickness
  Double_t x0  = 9.36; // Radiation length of Si [cm]
  Double_t rho = 2.33; // Density of Si [g/cm^3] 
  Double_t mass = fTrackToFollow->GetMass();
  fTrackToFollow->CorrectForMeanMaterial(xx0, xx0*x0*rho, mass, kTRUE);
  layer.ResetTrack(*fTrackToFollow); 

  return 1;
}

Int_t AliITSUTrackerCooked::PropagateBack(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // Here, we implement the Kalman smoother ?
  // The clusters must already be loaded
  //--------------------------------------------------------------------
    Int_t n=event->GetNumberOfTracks();
    for (Int_t i=0; i<n; i++) {
        AliESDtrack *esdTrack=event->GetTrack(i);
        AliITSUTrackCooked track(*esdTrack);
        esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSout);
    }
    
  return 0;
}

Int_t AliITSUTrackerCooked::RefitInward(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // Some final refit, after the outliers get removed by the smoother ?  
  // The clusters must be loaded
  //--------------------------------------------------------------------
    Int_t n=event->GetNumberOfTracks();
    for (Int_t i=0; i<n; i++) {
        AliESDtrack *esdTrack=event->GetTrack(i);
        AliITSUTrackCooked track(*esdTrack);
        esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSrefit);
    }
    
  return 0;
}

Int_t AliITSUTrackerCooked::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  // This function reads the ITSU clusters from the tree,
  // sort them, distribute over the internal tracker arrays, etc
  //--------------------------------------------------------------------
  if (!cTree) {
     AliFatal("No cluster tree !");
     return 1;
  }

  //This TClonesArray is not the owner of the clusters
  TClonesArray dummy("AliITSUClusterPix",kMaxClusterPerLayer), *clusters=&dummy;

  for (Int_t i=0; i<kNLayers; i++) {
      TBranch *br = cTree->GetBranch(Form("ITSRecPoints%d",i));
      if (!br) AliFatal(Form("No cluster branch for layer %d",i));
      br->SetAddress(&clusters);
      br->GetEvent(0);
      switch (i) {
      case kSeedingLayer1: 
      case kSeedingLayer2: 
      case kSeedingLayer3: 
         fgLayers[i].InsertClusters(clusters,kTRUE);
         break;
      default:
         fgLayers[i].InsertClusters(clusters,kFALSE);
         break;
      }
      clusters->Delete();
  }

  return 0;
}

void AliITSUTrackerCooked::UnloadClusters() {
  //--------------------------------------------------------------------
  // This function unloads ITSU clusters from the RAM
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kNLayers; i++) fgLayers[i].DeleteClusters();
}

AliCluster *AliITSUTrackerCooked::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
    Int_t l=(index & 0xf0000000) >> 28;
    Int_t c=(index & 0x0fffffff) >> 00;
    return fgLayers[l].GetCluster(c);
}

AliITSUTrackerCooked::AliITSUlayer::AliITSUlayer():
  fR(0),
  fN(0),
  fNsel(0),
  fTrack(0) 
{
  //--------------------------------------------------------------------
  // This default constructor needs to be provided
  //--------------------------------------------------------------------
  for (Int_t i=0; i<kMaxClusterPerLayer; i++) fClusters[i]=0;
  for (Int_t i=0; i<kMaxSelected; i++) fIndex[i]=-1;
}

AliITSUTrackerCooked::AliITSUlayer::~AliITSUlayer()
{
  //--------------------------------------------------------------------
  // Simple destructor
  //--------------------------------------------------------------------
  DeleteClusters();
  delete fTrack;
}

void 
AliITSUTrackerCooked::AliITSUlayer::ResetTrack(const AliITSUTrackCooked &t) {
  //--------------------------------------------------------------------
  // Replace the track estimate at this layer
  //--------------------------------------------------------------------
   delete fTrack;
   fTrack=new AliITSUTrackCooked(t);
}

void AliITSUTrackerCooked::
  AliITSUlayer::InsertClusters(TClonesArray *clusters, Bool_t seedingLayer)
{
  //--------------------------------------------------------------------
  // Load clusters to this layer
  //--------------------------------------------------------------------
  Int_t ncl=clusters->GetEntriesFast();
 
  while (ncl--) {
     AliITSUClusterPix *c=(AliITSUClusterPix*)clusters->UncheckedAt(ncl);
     (seedingLayer) ? c->GoToFrameGlo() : c->GoToFrameTrk();
     //if (!c->Misalign()) AliWarning("Can't misalign this cluster !");
     InsertCluster(new AliITSUClusterPix(*c));
  }

}

void AliITSUTrackerCooked::AliITSUlayer::DeleteClusters()
{
  //--------------------------------------------------------------------
  // Load clusters to this layer
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) {delete fClusters[i]; fClusters[i]=0;}
  fN=0;
}

Int_t 
AliITSUTrackerCooked::AliITSUlayer::InsertCluster(AliCluster *c) {
  //--------------------------------------------------------------------
  // This function inserts a cluster to this layer in increasing
  // order of the cluster's fZ
  //--------------------------------------------------------------------
  if (fN>=kMaxClusterPerLayer) {
     ::Error("InsertCluster","Too many clusters !\n");
     return 1;
  }
  if (fN==0) fClusters[0]=c;
  else {
     Int_t i=FindClusterIndex(c->GetZ());
     Int_t k=fN-i;
     memmove(fClusters+i+1 ,fClusters+i,k*sizeof(AliCluster*));
     fClusters[i]=c;
  }
  fN++;
  return 0;
}

Int_t 
AliITSUTrackerCooked::AliITSUlayer::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the first 
  // with its fZ >= "z". 
  //--------------------------------------------------------------------
  if (fN==0) return 0;

  Int_t b=0;
  if (z <= fClusters[b]->GetZ()) return b;

  Int_t e=b+fN-1;
  if (z > fClusters[e]->GetZ()) return e+1;

  Int_t m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}

Int_t AliITSUTrackerCooked::AliITSUlayer::
SelectClusters(Float_t zMin,Float_t zMax,Float_t phiMin, Float_t phiMax) {
  //--------------------------------------------------------------------
  // This function selects clusters within the "road"
  //--------------------------------------------------------------------
  UShort_t volId=-1;
  Float_t x, alpha;
  for (Int_t i=FindClusterIndex(zMin); i<fN; i++) {
      AliCluster *c=fClusters[i];
      UShort_t id=c->GetVolumeId();
      if (id != volId) {
	volId=id;
        c->GetXAlphaRefPlane(x,alpha); //FIXME
      }
      Double_t cPhi=alpha + c->GetY()/fR;
      if (c->IsClusterUsed()) continue;
      if (c->GetZ() > zMax) break;
      if (cPhi <= phiMin) continue;
      if (cPhi >  phiMax) continue;
      fIndex[fNsel++]=i;
      if (fNsel==kMaxSelected) break;
  }
  return fNsel;
}

const AliCluster *AliITSUTrackerCooked::AliITSUlayer::GetNextCluster(Int_t &ci){
  //--------------------------------------------------------------------
  // This function returns clusters within the "road" 
  //--------------------------------------------------------------------
  AliCluster *c=0;
  ci=-1;
  if (fNsel) {
     fNsel--;
     ci=fIndex[fNsel]; 
     c=fClusters[ci];
  }
  return c; 
}

