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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD propagation layer                                             //
//                                                                        //
//  Authors:                                                              //
//    Marian Ivanov <M.Ivanov@gsi.de>                                     //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "string.h"

#include "TMath.h"

#include "AliTRDpropagationLayer.h"
//#include "AliTRDtracker.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"

//_____________________________________________________________________________
AliTRDpropagationLayer::AliTRDpropagationLayer()
  :TObject()
  ,fN(0)
  ,fSec(0)
  ,fClusters(NULL)
  ,fIndex(NULL)
  ,fX(0.)
  ,fdX(0.)
  ,fRho(0.)
  ,fX0(0.)
  ,fTimeBinIndex(0)
  ,fPlane(0)
  ,fYmax(0)
  ,fYmaxSensitive(0)
  ,fHole(kFALSE)
  ,fHoleZc(0)
  ,fHoleZmax(0)
  ,fHoleYc(0)
  ,fHoleYmax(0)
  ,fHoleRho(0)
  ,fHoleX0(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDpropagationLayer::AliTRDpropagationLayer(Double_t x, Double_t dx, Double_t rho
                                             , Double_t radLength, Int_t tbIndex, Int_t plane)
  :TObject()
  ,fN(0)
  ,fSec(0)
  ,fClusters(NULL)
  ,fIndex(NULL)
  ,fX(x)
  ,fdX(dx)
  ,fRho(rho)
  ,fX0(radLength)
  ,fTimeBinIndex(tbIndex)
  ,fPlane(plane)
  ,fYmax(0)
  ,fYmaxSensitive(0)
  ,fHole(kFALSE)
  ,fHoleZc(0)
  ,fHoleZmax(0)
  ,fHoleYc(0)
  ,fHoleYmax(0)
  ,fHoleRho(0)
  ,fHoleX0(0)
{ 
  //
  // AliTRDpropagationLayer constructor
  //

  for (Int_t i = 0; i < (Int_t)kZones; i++) {
    fZc[i]   = 0; 
    fZmax[i] = 0;
  }

  if (fTimeBinIndex >= 0) { 
    fClusters = new AliTRDcluster*[kMaxClusterPerTimeBin];
    fIndex    = new UInt_t[kMaxClusterPerTimeBin];
  }

  for (Int_t i = 0; i < 5; i++) {
    fIsHole[i] = kFALSE;
  }

}

//_____________________________________________________________________________
AliTRDpropagationLayer::AliTRDpropagationLayer(const AliTRDpropagationLayer &p)
  :TObject((TObject&)p)
  ,fN(p.fN)
  ,fSec(p.fSec)
  ,fClusters(0x0)
  ,fIndex(0x0)
  ,fX(p.fX)
  ,fdX(p.fdX)
  ,fRho(p.fRho)
  ,fX0(p.fX0)
  ,fTimeBinIndex(p.fTimeBinIndex)
  ,fPlane(p.fPlane)
  ,fYmax(p.fYmax)
  ,fYmaxSensitive(p.fYmaxSensitive)
  ,fHole(p.fHole)
  ,fHoleZc(p.fHoleZc)
  ,fHoleZmax(p.fHoleZmax)
  ,fHoleYc(p.fHoleYc)
  ,fHoleYmax(p.fHoleYmax)
  ,fHoleRho(p.fHoleRho)
  ,fHoleX0(p.fHoleX0)
{
  //
  // AliTRDpropagationLayer copy constructor
  //

  for (Int_t i = 0; i < (Int_t)kZones; i++) {
    fZc[i]   = p.fZc[i]; 
    fZmax[i] = p.fZmax[i];
		fIsHole[i] = p.fIsHole[i];
		fZmaxSensitive[i] = p.fZmaxSensitive[i];  
	}

	// Make a deep copy of the Clusters array and the Index array unless they are needed in class AliTRDstackLayer
	Int_t arrsize = (Int_t)kMaxClusterPerTimeBin;
	 if (fTimeBinIndex >= 0) { 
    fClusters = new AliTRDcluster*[arrsize];
    fIndex    = new UInt_t[arrsize];
  }
	memset(fIndex, 0, sizeof(UInt_t)*arrsize);
	memset(fClusters, 0, sizeof(AliTRDcluster *)*arrsize);
	for(Int_t i = 0; i < arrsize; i++){
		fClusters[i] = p.fClusters[i];
		fIndex[i] = p.fIndex[i];
	}
}
 
//_____________________________________________________________________________
AliTRDpropagationLayer::~AliTRDpropagationLayer()
{
  //
  // Destructor
  //

  if (fTimeBinIndex >= 0) { 
    delete[] fClusters;
    delete[] fIndex;
  }

}

//_____________________________________________________________________________
void AliTRDpropagationLayer::Copy(TObject &o) const 
{
  //
  // Copy function
  //

  AliTRDpropagationLayer &p = (AliTRDpropagationLayer &)o; 
  p.fN   = fN;
  p.fSec = fSec;
  p.fX = fX;
  p.fdX = fdX;
  p.fRho = fRho;
  p.fX0  = fX0;
  p.fTimeBinIndex = fTimeBinIndex;
  p.fPlane  = fPlane;
  p.fYmax = fYmax;
  p.fYmaxSensitive  = fYmaxSensitive;
  p.fHole = fHole;
  p.fHoleZc = fHoleZc;
  p.fHoleZmax = fHoleZmax;
  p.fHoleYc = fHoleYc;
  p.fHoleYmax = fHoleYmax;
  p.fHoleRho = fHoleRho;
  p.fHoleX0 = fHoleX0;

  for (Int_t i = 0; i < (Int_t)kZones; i++) {
    p.fZc[i]   = fZc[i]; 
    p.fZmax[i] = fZmax[i];
    p.fIsHole[i] = fIsHole[i];
    p.fZmaxSensitive[i] = fZmaxSensitive[i];  
  }
	
  // Make a deep copy of the Clusters array and the Index array
  // unless they are needed in class AliTRDstackLayer
  if (fTimeBinIndex >= 0) { 
    if (!p.fClusters) 
      p.fClusters = new AliTRDcluster*[(Int_t)kMaxClusterPerTimeBin];
    if (!p.fIndex)
      p.fIndex    = new UInt_t[(Int_t)kMaxClusterPerTimeBin];
  }
  for (Int_t i = 0; i < (Int_t)kMaxClusterPerTimeBin; i++){
    //overwrite
    p.fClusters[i] = fClusters[i];
    p.fIndex[i] = fIndex[i];
  }

}

//_____________________________________________________________________________
void AliTRDpropagationLayer::SetZ(Double_t *center, Double_t *w, Double_t *wsensitive )
{
  //
  // Set centers and the width of sectors
  //

  for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
    fZc[icham]            = center[icham];  
    fZmax[icham]          = w[icham];
    fZmaxSensitive[icham] = wsensitive[icham];
  }  

}

//_____________________________________________________________________________
void AliTRDpropagationLayer::SetHoles(Bool_t *holes)
{
  //
  // Set centers and the width of sectors
  //

  fHole = kFALSE;

  for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
    fIsHole[icham] = holes[icham]; 
    if (holes[icham]) {
      fHole = kTRUE;
    }
  }  

}

//_____________________________________________________________________________
void AliTRDpropagationLayer::InsertCluster(AliTRDcluster *c, UInt_t index) 
{
  //
  // Insert cluster in cluster array.
  // Clusters are sorted according to Y coordinate.  
  //

  if (fTimeBinIndex < 0) { 
    //AliWarning("Attempt to insert cluster into non-sensitive time bin!\n");
    return;
  }

  if (fN == (Int_t) kMaxClusterPerTimeBin) {
    //AliWarning("Too many clusters !\n"); 
    return;
  }

  if (fN == 0) {
    fIndex[0]       = index; 
    fClusters[fN++] = c; 
    return;
  }

  Int_t i = Find(c->GetY());
  memmove(fClusters+i+1,fClusters+i,(fN-i)*sizeof(AliTRDcluster*));
  memmove(fIndex   +i+1,fIndex   +i,(fN-i)*sizeof(UInt_t)); 
  fIndex[i]    = index; 
  fClusters[i] = c; 
  fN++;

}  

//_____________________________________________________________________________
Int_t AliTRDpropagationLayer::Find(Float_t y) const
{
  //
  // Returns index of the cluster nearest in Y    
  //

  if (fN <= 0) {
    return 0;
  }
  if (y <= fClusters[0]->GetY()) {
    return 0;
  }
  if (y >  fClusters[fN-1]->GetY()) {
    return fN;
  }

  Int_t b = 0;
  Int_t e = fN - 1;
  Int_t m = (b + e) / 2;

  for ( ; b < e; m = (b + e) / 2) {
    if (y > fClusters[m]->GetY()) {
      b = m + 1;
    }
    else {
      e = m;
    }
  }

  return m;

}    

//_____________________________________________________________________________
Int_t AliTRDpropagationLayer::FindNearestCluster(Float_t y, Float_t z
                                               , Float_t maxroad
                                               , Float_t maxroadz) const 
{
  //
  // Returns index of the cluster nearest to the given y,z
  //

  Int_t   index   = -1;
  Int_t   maxn    = fN;
  Float_t mindist = maxroad;			

  for (Int_t i = Find(y-maxroad); i < maxn; i++) {
    AliTRDcluster *c = (AliTRDcluster *) (fClusters[i]);
    Float_t ycl = c->GetY();
    if (ycl > (y + maxroad)) {
      break;
    }
    if (TMath::Abs(c->GetZ() - z) > maxroadz) {
      continue;
    }
    if (TMath::Abs(ycl - y)       < mindist) {
      mindist = TMath::Abs(ycl - y);
      index   = fIndex[i];
    }
  }						

  return index;

}             

