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


//-------------------------------------------------------
//          Implementation of the TPC tracker helper clasess
//  AliTPCtrackerRow
//  AliTPCtrackerSector
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//  AliTPCtrakerMI -  parallel tracker helper clases
//

/* $Id: AliTPCtrackerSector.cxx 25837 2008-05-16 16:39:00Z marian $ */

#include "Riostream.h"
#include <TClonesArray.h>
#include "AliLog.h"
#include "AliComplexCluster.h"
#include "AliTPCcluster.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "AliTPCReconstructor.h"
#include "AliTPCreco.h"
//
#include "AliTPCtrackerSector.h"
#include "TStopwatch.h"
#include "TTreeStream.h"

//

ClassImp(AliTPCtrackerRow)
ClassImp(AliTPCtrackerSector)



AliTPCtrackerRow::AliTPCtrackerRow():
  fDeadZone(0.),
  fClusters1(0),
  fN1(0),
  fClusters2(0),
  fN2(0),
  fFastCluster(),
  fN(0),
  fClusters(),
  fIndex(),
  fX(0.)  
{
  //
  // default constructor
  //
}

AliTPCtrackerRow::~AliTPCtrackerRow(){
  //
  for (Int_t i = 0; i < fN1; i++)
    fClusters1[i].~AliTPCclusterMI();
  delete [] fClusters1;
  for (Int_t i = 0; i < fN2; i++)
    fClusters2[i].~AliTPCclusterMI();
  delete [] fClusters2;
}



//_________________________________________________________________________
void 
AliTPCtrackerRow::InsertCluster(const AliTPCclusterMI* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMaxClusterPerRow) {
    //AliInfo("AliTPCtrackerRow::InsertCluster(): Too many clusters"); 
    return;
  }
  if (fN>=fN1+fN2) {
    //AliInfo("AliTPCtrackerRow::InsertCluster(): Too many clusters !");
  }

  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetZ());
  if (i>=0 && i<=kMaxClusterPerRow-2) {
    memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTPCclusterMI*));
    memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  }
  fIndex[i]=index; fClusters[i]=c; fN++;
}

void AliTPCtrackerRow::ResetClusters() {
   //
   // reset clusters
   // MvL: Need to call destructors for AliTPCclusterMI, to delete fInfo
   for (Int_t i = 0; i < fN1; i++)
     fClusters1[i].~AliTPCclusterMI();
   delete [] fClusters1;  fClusters1=0;
   for (Int_t i = 0; i < fN2; i++)
     fClusters2[i].~AliTPCclusterMI();
   delete [] fClusters2;  fClusters2=0;

   fN  = 0; 
   fN1 = 0;
   fN2 = 0;
   //delete[] fClusterArray; 

   //fClusterArray=0;
}


//___________________________________________________________________
Int_t AliTPCtrackerRow::Find(Double_t z) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster 
  //-----------------------------------------------------------------------
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



//___________________________________________________________________
AliTPCclusterMI * AliTPCtrackerRow::FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;

  AliTPCclusterMI *cl =0;
  for (Int_t i=Find(z-roadz); i<fN; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( (c->GetY()-y) >  roady ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
      }
  }
  return cl;      
}

AliTPCclusterMI * AliTPCtrackerRow::FindNearest2(Double_t y, Double_t z, Double_t roady, Double_t roadz,UInt_t & index) const 
{
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = roady*roady + roadz*roadz;
  AliTPCclusterMI *cl =0;

  //PH Check boundaries. 510 is the size of fFastCluster
  Int_t iz1 = Int_t(z-roadz+254.5);
  if (iz1<0 || iz1>=510) return cl;
  iz1 = TMath::Max(GetFastCluster(iz1)-1,0);
  Int_t iz2 = Int_t(z+roadz+255.5);
  if (iz2<0 || iz2>=510) return cl;
  iz2 = TMath::Min(GetFastCluster(iz2)+1,fN);
  Bool_t skipUsed = !(AliTPCReconstructor::GetRecoParam()->GetClusterSharing());
  //FindNearest3(y,z,roady,roadz,index);
  //  for (Int_t i=Find(z-roadz); i<fN; i++) {
  for (Int_t i=iz1; i<iz2; i++) {
      AliTPCclusterMI *c=(AliTPCclusterMI*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      if ( c->GetY()-y >  roady ) continue;
      if ( y-c->GetY() >  roady ) continue;
      if (skipUsed && c->IsUsed(11)) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)+(c->GetY()-y)*(c->GetY()-y);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
	index =i;
	//roady = TMath::Sqrt(maxdistance);
      }
  }
  return cl;      
}


void AliTPCtrackerRow::SetFastCluster(Int_t i, Short_t cl){
  //
  // Set cluster info for fast navigation
  //
  if (i>=510|| i<0){
  }else{
    fFastCluster[i]=cl;
  }
}


Int_t  AliTPCtrackerSector::GetRowNumber(Double_t x) const 
{
  //return pad row number for this x
  Double_t r;
  if (fN < 64){
    r=fRow[fN-1].GetX();
    if (x > r) return fN;
    r=fRow[0].GetX();
    if (x < r) return -1;
    return Int_t((x-r)/fPadPitchLength + 0.5);}
  else{    
    r=fRow[fN-1].GetX();
    if (x > r) return fN;
    r=fRow[0].GetX();
    if (x < r) return -1;
    Double_t r1=fRow[64].GetX();
    if(x<r1){       
      return Int_t((x-r)/f1PadPitchLength + 0.5);}
    else{
      return (Int_t((x-r1)/f2PadPitchLength + 0.5)+64);} 
  }
}

//_________________________________________________________________________
void AliTPCtrackerSector::Setup(const AliTPCParam *par, Int_t f) {
  //-----------------------------------------------------------------------
  // Setup inner sector
  //-----------------------------------------------------------------------
  if (f==0) {
     fAlpha=par->GetInnerAngle();
     fAlphaShift=par->GetInnerAngleShift();
     fPadPitchWidth=par->GetInnerPadPitchWidth();
     fPadPitchLength=par->GetInnerPadPitchLength();
     fN=par->GetNRowLow();
     if(fRow)delete [] fRow;fRow = 0;
     fRow=new AliTPCtrackerRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiLow(i));
       fRow[i].SetDeadZone(1.5);  //1.5 cm of dead zone
     }
  } else {
     fAlpha=par->GetOuterAngle();
     fAlphaShift=par->GetOuterAngleShift();
     fPadPitchWidth  = par->GetOuterPadPitchWidth();
     fPadPitchLength = par->GetOuter1PadPitchLength();
     f1PadPitchLength = par->GetOuter1PadPitchLength();
     f2PadPitchLength = par->GetOuter2PadPitchLength();
     fN=par->GetNRowUp();
     if(fRow)delete [] fRow;fRow = 0;
     fRow=new AliTPCtrackerRow[fN];
     for (Int_t i=0; i<fN; i++) {
       fRow[i].SetX(par->GetPadRowRadiiUp(i)); 
       fRow[i].SetDeadZone(1.5);  // 1.5 cm of dead zone
     }
  } 
}

//_________________________________________________________________________
void AliTPCtrackerSector::InsertCluster(AliTPCclusterMI *cl, Int_t size, const AliTPCParam *par) {
  //-----------------------------------------------------------------------
  // Insert cluster to the sector
  //-----------------------------------------------------------------------

  if(!cl) return; 

  const Int_t fkNIS = par->GetNInnerSector()/2;
  const Int_t fkNOS = par->GetNOuterSector()/2;
  Int_t row = cl->GetRow();
  Int_t sec = cl->GetDetector();

  // add cluster to the corresponding pad row
  AliTPCtrackerRow *tpcrow = 0x0;

  Int_t left=0;
  if (sec<fkNIS*2){
    left = sec/fkNIS;
  }
  else{
    left = (sec-fkNIS*2)/fkNOS;
  }
  //
  if (left ==0){
    tpcrow = fRow+row;
    if(!tpcrow->GetClusters1()) {
       tpcrow->SetClusters1(new AliTPCclusterMI[size]); 
       tpcrow->SetN1(0);
    }
    if(size < kMaxClusterPerRow) {
      tpcrow->SetCluster1(tpcrow->GetN1(), *cl);
      //printf("inner: size %d, tpcrow->GetN1() %d  sec %d row %d tpcrow %p cl %p\n", size, tpcrow->GetN1(), sec, row, tpcrow, cl);

      tpcrow->IncrementN1();
    }
  }
  if (left ==1){
    tpcrow = fRow+row;
    if(!tpcrow->GetClusters2()) { 
      tpcrow->SetClusters2(new AliTPCclusterMI[size]); 
      tpcrow->SetN2(0);
    }
    if(size < kMaxClusterPerRow)  { 
      tpcrow->SetCluster2(tpcrow->GetN2(), *cl);
      //printf("outer: size %d, tpcrow->GetN2() %d  sec %d row %d tpcrow %p cl %p\n", size, tpcrow->GetN2(), sec, row, tpcrow, cl);

      tpcrow->IncrementN2();
    }
  }
}



