// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************


#include "AliHLTITSLayer.h"
#include <algorithm>


//------------------------------------------------------------------------
AliHLTITSLayer::AliHLTITSLayer():
fR(0),
fPhiOffset(0),
fNladders(0),
fZOffset(0),
fNdetectors(0),
fDetectors(0),
fN(0),
fDy5(0),
fDy10(0),
fDy20(0),
fClustersCs(0),
fClusterIndexCs(0),
fYcs(0),
fZcs(0),
fNcs(0),
fCurrentSlice(-1),
fZmax(0),
fYmin(0),
fYmax(0),
fI(0),
fImax(0),
fSkip(0),
fAccepted(0),
fRoad(0)
{
  //--------------------------------------------------------------------
  //default AliHLTITSLayer constructor
  //--------------------------------------------------------------------
  for (Int_t i=0;i<6;i++)  fN5[i] =0;  
  for (Int_t i=0;i<11;i++) fN10[i]=0;  
  for (Int_t i=0;i<21;i++) fN20[i]=0;
  for (Int_t i=0;i<AliITSRecoParam::fgkMaxClusterPerLayer;i++){
    fClusters[i] = 0;
    fClusterIndex[i] = 0;
    fZ[i]        = 0;
    fY[i] = 0;
  }
  fYB[0] = 0;
  fYB[1] = 0;
}

//------------------------------------------------------------------------
AliHLTITSLayer::
AliHLTITSLayer(Double_t r,Double_t p,Double_t z,Int_t nl,Int_t nd):
fR(r),
fPhiOffset(p),
fNladders(nl),
fZOffset(z),
fNdetectors(nd),
fDetectors(0),
fN(0),
fDy5(0),
fDy10(0),
fDy20(0),
fClustersCs(0),
fClusterIndexCs(0),
fYcs(0),
fZcs(0),
fNcs(0),
fCurrentSlice(-1),
fZmax(0),
fYmin(0),
fYmax(0),
fI(0),
fImax(0),
fSkip(0),
fAccepted(0),
fRoad(0) 
{
  //--------------------------------------------------------------------
  //main AliHLTITSLayer constructor
  //--------------------------------------------------------------------
  fDetectors=new AliHLTITSDetector[fNladders*fNdetectors];
  fRoad=2*fR*TMath::Sqrt(TMath::Pi()/1.);//assuming that there's only one cluster

  for (Int_t i=0;i<6;i++)  fN5[i] =0;  
  for (Int_t i=0;i<11;i++) fN10[i]=0;  
  for (Int_t i=0;i<21;i++) fN20[i]=0;
  for (Int_t i=0;i<AliITSRecoParam::fgkMaxClusterPerLayer;i++){
    fClusters[i] = 0;
    fClusterIndex[i] = 0;
    fZ[i]        = 0;
    fY[i] = 0;
  }
  fYB[0] = 0;
  fYB[1] = 0;
}

//------------------------------------------------------------------------
AliHLTITSLayer::AliHLTITSLayer(const AliHLTITSLayer& layer):
fR(layer.fR),
fPhiOffset(layer.fPhiOffset),
fNladders(layer.fNladders),
fZOffset(layer.fZOffset),
fNdetectors(layer.fNdetectors),
fDetectors(layer.fDetectors),
fN(layer.fN),
fDy5(layer.fDy5),
fDy10(layer.fDy10),
fDy20(layer.fDy20),
fClustersCs(layer.fClustersCs),
fClusterIndexCs(layer.fClusterIndexCs),
fYcs(layer.fYcs),
fZcs(layer.fZcs),
fNcs(layer.fNcs),
fCurrentSlice(layer.fCurrentSlice),
fZmax(layer.fZmax),
fYmin(layer.fYmin),
fYmax(layer.fYmax),
fI(layer.fI),
fImax(layer.fImax),
fSkip(layer.fSkip),
fAccepted(layer.fAccepted),
fRoad(layer.fRoad)
{
  //Copy constructor
  for (Int_t i=0;i<6;i++)  fN5[i] =0;  
  for (Int_t i=0;i<11;i++) fN10[i]=0;  
  for (Int_t i=0;i<21;i++) fN20[i]=0;
  for (Int_t i=0;i<AliITSRecoParam::fgkMaxClusterPerLayer;i++){
    fClusters[i] = 0;
    fClusterIndex[i] = 0;
    fZ[i]        = 0;
    fY[i] = 0;
  }
  fYB[0] = 0;
  fYB[1] = 0;
}

//------------------------------------------------------------------------
AliHLTITSLayer::~AliHLTITSLayer() {
  //--------------------------------------------------------------------
  // AliHLTITSLayer destructor
  //--------------------------------------------------------------------
  delete [] fDetectors;
}
//------------------------------------------------------------------------
void AliHLTITSLayer::ResetClusters() {
  //--------------------------------------------------------------------
  // This function removes loaded clusters
  //--------------------------------------------------------------------
  
  fN=0;
  fI=0;
}


//------------------------------------------------------------------------
void AliHLTITSLayer::ResetRoad() {
  //--------------------------------------------------------------------
  // This function calculates the road defined by the cluster density
  //--------------------------------------------------------------------
  Int_t n=0;
  for (Int_t i=0; i<fN; i++) {
     if (TMath::Abs(fClusters[i]->GetZ())<fR) n++;
  }
  if (n>1) fRoad=2*fR*TMath::Sqrt(TMath::Pi()/n);
}
//------------------------------------------------------------------------
Int_t AliHLTITSLayer::InsertCluster(AliITSRecPoint *cl) {
  //--------------------------------------------------------------------
  //This function adds a cluster to this layer
  //--------------------------------------------------------------------
  if (fN==AliITSRecoParam::GetMaxClusterPerLayer()) {
    ::Error("InsertCluster","Too many clusters !\n");
    return 1;
  }
  fCurrentSlice=-1;
  fClusters[fN]=cl;
  fN++;
  AliHLTITSDetector &det=GetDetector(cl->GetDetectorIndex());    
  if (cl->GetY()<det.GetYmin()) det.SetYmin(cl->GetY());
  if (cl->GetY()>det.GetYmax()) det.SetYmax(cl->GetY());
  if (cl->GetZ()<det.GetZmin()) det.SetZmin(cl->GetZ());
  if (cl->GetZ()>det.GetZmax()) det.SetZmax(cl->GetZ());
			     
  return 0;
}


//------------------------------------------------------------------------
void  AliHLTITSLayer::SortClusters()
{
  //
  //sort clusters
  //
 
  AliITSRecPoint **clusters = new AliITSRecPoint*[fN];
  Float_t *z                = new Float_t[fN];
  Int_t   * index           = new Int_t[fN];
  //
  for (Int_t i=0;i<fN;i++){
    z[i] = fClusters[i]->GetZ();
  }
  TMath::Sort(fN,z,index,kFALSE);
  
  for (Int_t i=0;i<fN;i++){
    clusters[i] = fClusters[index[i]];
  }
  //
  for (Int_t i=0;i<fN;i++){
    fClusters[i] = clusters[i];
    fZ[i]        = fClusters[i]->GetZ();
    AliHLTITSDetector &det=GetDetector(fClusters[i]->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + fClusters[i]->GetY();
    if (y>2.*fR*TMath::Pi()) y -= 2.*fR*TMath::Pi();
    fY[i] = y;
  }
  delete[] index;
  delete[] z;
  delete[] clusters;
  //

  fYB[0]=10000000;
  fYB[1]=-10000000;
  for (Int_t i=0;i<fN;i++){
    if (fY[i]<fYB[0]) fYB[0]=fY[i];
    if (fY[i]>fYB[1]) fYB[1]=fY[i];
    fClusterIndex[i] = i;
  }
  //
  // fill slices
  fDy5 = (fYB[1]-fYB[0])/5.;
  fDy10 = (fYB[1]-fYB[0])/10.;
  fDy20 = (fYB[1]-fYB[0])/20.;
  for (Int_t i=0;i<6;i++)  fN5[i] =0;  
  for (Int_t i=0;i<11;i++) fN10[i]=0;  
  for (Int_t i=0;i<21;i++) fN20[i]=0;
  //  
  for (Int_t i=0;i<6;i++) {fBy5[i][0] =  fYB[0]+(i-0.75)*fDy5; fBy5[i][1] =  fYB[0]+(i+0.75)*fDy5;}
  for (Int_t i=0;i<11;i++) {fBy10[i][0] =  fYB[0]+(i-0.75)*fDy10; fBy10[i][1] =  fYB[0]+(i+0.75)*fDy10;} 
  for (Int_t i=0;i<21;i++) {fBy20[i][0] =  fYB[0]+(i-0.75)*fDy20; fBy20[i][1] =  fYB[0]+(i+0.75)*fDy20;}
  //
  //
  for (Int_t i=0;i<fN;i++)
    for (Int_t irot=-1;irot<=1;irot++){
      Float_t curY = fY[i]+irot*TMath::TwoPi()*fR; 
      // slice 5
      for (Int_t slice=0; slice<6;slice++){
	if (fBy5[slice][0]<curY && curY<fBy5[slice][1]&&fN5[slice]<AliITSRecoParam::GetMaxClusterPerLayer5()){
	  fClusters5[slice][fN5[slice]] = fClusters[i];
	  fY5[slice][fN5[slice]] = curY;
	  fZ5[slice][fN5[slice]] = fZ[i];
	  fClusterIndex5[slice][fN5[slice]]=i;
	  fN5[slice]++;
	}
      }
      // slice 10
      for (Int_t slice=0; slice<11;slice++){
	if (fBy10[slice][0]<curY && curY<fBy10[slice][1]&&fN10[slice]<AliITSRecoParam::GetMaxClusterPerLayer10()){
	  fClusters10[slice][fN10[slice]] = fClusters[i];
	  fY10[slice][fN10[slice]] = curY;
	  fZ10[slice][fN10[slice]] = fZ[i];
	  fClusterIndex10[slice][fN10[slice]]=i;
	  fN10[slice]++;
	}
      }
      // slice 20
      for (Int_t slice=0; slice<21;slice++){
	if (fBy20[slice][0]<curY && curY<fBy20[slice][1]&&fN20[slice]<AliITSRecoParam::GetMaxClusterPerLayer20()){
	  fClusters20[slice][fN20[slice]] = fClusters[i];
	  fY20[slice][fN20[slice]] = curY;
	  fZ20[slice][fN20[slice]] = fZ[i];
	  fClusterIndex20[slice][fN20[slice]]=i;
	  fN20[slice]++;
	}
      }      
    }

  //
  // consistency check
  //
  for (Int_t i=0;i<fN-1;i++){
    if (fZ[i]>fZ[i+1]){
      printf("Bug\n");
    }
  }
  //
  for (Int_t slice=0;slice<21;slice++)
  for (Int_t i=0;i<fN20[slice]-1;i++){
    if (fZ20[slice][i]>fZ20[slice][i+1]){
      printf("Bug\n");
    }
  }


}
//------------------------------------------------------------------------
Int_t AliHLTITSLayer::FindClusterIndex(Float_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  Int_t ncl=0;
  const Float_t *zcl;  
  if (fCurrentSlice<0) {
    ncl = fN;
    zcl   = fZ;
  }
  else{
    ncl   = fNcs;
    zcl   = fZcs;;
  }
  
  if (ncl==0) return 0;
  Int_t b=0, e=ncl-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    //    if (z > fClusters[m]->GetZ()) b=m+1;
    if (z > zcl[m]) b=m+1;
    else e=m; 
  }
  return m;
}

//------------------------------------------------------------------------
void AliHLTITSLayer::
SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin, Double_t ymax) {
  //--------------------------------------------------------------------
  // This function sets the "window"
  //--------------------------------------------------------------------
 
  Double_t circle=2*TMath::Pi()*fR;
  fYmin = ymin; fYmax =ymax;
  Float_t ymiddle = (fYmin+fYmax)*0.5;
  if (ymiddle<fYB[0]) {
    fYmin+=circle; fYmax+=circle; ymiddle+=circle;
  } else if (ymiddle>fYB[1]) {
    fYmin-=circle; fYmax-=circle; ymiddle-=circle;
  }
  
  //
  fCurrentSlice =-1;
  // defualt take all
  fClustersCs = fClusters;
  fClusterIndexCs = fClusterIndex;
  fYcs  = fY;
  fZcs  = fZ;
  fNcs  = fN;
  //
  //is in 20 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy20){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy20);
    if (slice<0) slice=0;
    if (slice>20) slice=20;
    Bool_t isOK = (fYmin>fBy20[slice][0]&&fYmax<fBy20[slice][1]);
    if (isOK) {
      fCurrentSlice=slice;
      fClustersCs = fClusters20[fCurrentSlice];
      fClusterIndexCs = fClusterIndex20[fCurrentSlice];
      fYcs  = fY20[fCurrentSlice];
      fZcs  = fZ20[fCurrentSlice];
      fNcs  = fN20[fCurrentSlice];
    }
  }  
  //
  //is in 10 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy10){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy10);
    if (slice<0) slice=0;
    if (slice>10) slice=10;
    Bool_t isOK = (fYmin>fBy10[slice][0]&&fYmax<fBy10[slice][1]);
    if (isOK) {
      fCurrentSlice=slice;
      fClustersCs = fClusters10[fCurrentSlice];
      fClusterIndexCs = fClusterIndex10[fCurrentSlice];
      fYcs  = fY10[fCurrentSlice];
      fZcs  = fZ10[fCurrentSlice];
      fNcs  = fN10[fCurrentSlice];
    }
  }  
  //
  //is in 5 slice?
  if (fCurrentSlice<0&&TMath::Abs(fYmax-fYmin)<1.49*fDy5){
    Int_t slice = int(0.5+(ymiddle-fYB[0])/fDy5);
    if (slice<0) slice=0;
    if (slice>5) slice=5;
    Bool_t isOK = (fYmin>fBy5[slice][0]&&fYmax<fBy5[slice][1]);
    if (isOK) {
      fCurrentSlice=slice;
      fClustersCs = fClusters5[fCurrentSlice];
      fClusterIndexCs = fClusterIndex5[fCurrentSlice];
      fYcs  = fY5[fCurrentSlice];
      fZcs  = fZ5[fCurrentSlice];
      fNcs  = fN5[fCurrentSlice];
    }
  }  
  //  
  fI=FindClusterIndex(zmin); fZmax=zmax;
  fImax = TMath::Min(FindClusterIndex(zmax)+1,fNcs);
  fSkip = 0;
  fAccepted =0;

  return;
}
//------------------------------------------------------------------------
Int_t AliHLTITSLayer::
FindDetectorIndex(Double_t phi, Double_t z) const {
  //--------------------------------------------------------------------
  //This function finds the detector crossed by the track
  //--------------------------------------------------------------------
  Double_t dphi;
  if (fZOffset<0)            // old geometry
    dphi = -(phi-fPhiOffset);
  else                       // new geometry
    dphi = phi-fPhiOffset;


  if      (dphi <  0) dphi += 2*TMath::Pi();
  else if (dphi >= 2*TMath::Pi()) dphi -= 2*TMath::Pi();
  Int_t np=Int_t(dphi*fNladders*0.5/TMath::Pi()+0.5);
  if (np>=fNladders) np-=fNladders;
  if (np<0)          np+=fNladders;


  Double_t dz=fZOffset-z;
  Double_t nnz = dz*(fNdetectors-1)*0.5/fZOffset+0.5;
  Int_t nz = (nnz<0 ? -1 : (Int_t)nnz);
  if (nz>=fNdetectors) return -1;
  if (nz<0)            return -1;

  // ad hoc correction for 3rd ladder of SDD inner layer,
  // which is reversed (rotated by pi around local y)
  // this correction is OK only from AliITSv11Hybrid onwards
  if (GetR()>12. && GetR()<20.) { // SDD inner
    if(np==2) { // 3rd ladder
      nz = (fNdetectors-1) - nz;
    } 
  }
  //printf("ndet %d phi %f z %f  np %d nz %d\n",fNdetectors,phi,z,np,nz);


  return np*fNdetectors + nz;
}
//------------------------------------------------------------------------
const AliITSRecPoint *AliHLTITSLayer::GetNextCluster(Int_t &ci,Bool_t test)
{
  //--------------------------------------------------------------------
  // This function returns clusters within the "window" 
  //--------------------------------------------------------------------

  if (fCurrentSlice<0) {
    Double_t rpi2 = 2.*fR*TMath::Pi();
    for (Int_t i=fI; i<fImax; i++) {
      Double_t y = fY[i];
      if (fYmax<y) y -= rpi2;
      if (fYmin>y) y += rpi2;
      if (y<fYmin) continue;
      if (y>fYmax) continue;
      if (fClusters[i]->GetQ()==0&&fSkip==2) continue;
      ci=i;
      if (!test) fI=i+1;
      return fClusters[i];
    }
  } else {
    for (Int_t i=fI; i<fImax; i++) {
      if (fYcs[i]<fYmin) continue;
      if (fYcs[i]>fYmax) continue;
      if (fClustersCs[i]->GetQ()==0&&fSkip==2) continue;
      ci=fClusterIndexCs[i];
      if (!test) fI=i+1;
      return fClustersCs[i];
    }
  }
  return 0;
}
//------------------------------------------------------------------------
Double_t AliHLTITSLayer::GetThickness(Double_t y,Double_t z,Double_t &x0)
const {
  //--------------------------------------------------------------------
  // This function returns the layer thickness at this point (units X0)
  //--------------------------------------------------------------------
  Double_t d=0.0085;
  x0=AliITSRecoParam::GetX0Air();
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
     if (TMath::Abs(y-3.03)<0.10) d+=0.014;
  } else
  if (3<fR&&fR<5) {   //SPD1
     Double_t dd=0.0063; x0=21.5;
     d=dd;
     if (TMath::Abs(y+0.21)>0.6) d+=dd;
     if (TMath::Abs(y+0.10)<0.10) d+=0.014;
  }

  return d;
}

//------------------------------------------------------------------------
Int_t AliHLTITSLayer::InRoad() const {
  //-------------------------------------------------------------------
  // This function returns number of clusters within the "window" 
  //--------------------------------------------------------------------
  Int_t ncl=0;
  for (Int_t i=fI; i<fN; i++) {
    const AliITSRecPoint *c=fClusters[i];
    if (c->GetZ() > fZmax) break;
    if (c->IsUsed()) continue;
    const AliHLTITSDetector &det=GetDetector(c->GetDetectorIndex());    
    Double_t y=fR*det.GetPhi() + c->GetY();

    if (y>2.*fR*TMath::Pi()) y -= 2*fR*TMath::Pi();
    if (y>1.*fR*TMath::Pi() && fYmax<y) y -= 2*fR*TMath::Pi();

    if (y<fYmin) continue;
    if (y>fYmax) continue;
    ncl++;
  }
  return ncl;
}

