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

/* $Id: AliTRDchamberTimeBin.cxx 23313 2008-01-11 14:56:43Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Organization of clusters at the level of 1 TRD chamber.                  //
//  The data structure is used for tracking at the stack level.              //
//                                                                           //
//  Functionalities:                                                         //
//   1. cluster organization and sorting                                     //
//   2. fast data navigation                                                 //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>
#include <TTreeStream.h>

#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDchamberTimeBin.h"
#include "AliTRDrecoParam.h"
#include "AliTRDReconstructor.h"

ClassImp(AliTRDchamberTimeBin)

//_____________________________________________________________________________
AliTRDchamberTimeBin::AliTRDchamberTimeBin(Int_t plane, Int_t stack, Int_t sector, Double_t z0, Double_t zLength)
  :TObject()
  ,fkReconstructor(NULL)
  ,fPlane(plane)
  ,fStack(stack)
  ,fSector(sector)
  ,fNRows(kMaxRows)
  ,fN(0)
  ,fX(0.)
  ,fZ0(z0)
  ,fZLength(zLength)
{
  //
  // Default constructor (Only provided to use AliTRDchamberTimeBin with arrays)
  //
  SetBit(kT0, kFALSE);
  SetBit(kOwner, kFALSE);
  memset(fPositions, 1, kMaxRows*sizeof(UChar_t));
  memset(fClusters, 0, kMaxClustersLayer*sizeof(AliTRDcluster*));
  memset(fIndex, 1, kMaxClustersLayer*sizeof(UInt_t));
}


//_____________________________________________________________________________
AliTRDchamberTimeBin::AliTRDchamberTimeBin(const AliTRDchamberTimeBin &layer):
  TObject()
  ,fkReconstructor(layer.fkReconstructor)
  ,fPlane(layer.fPlane)
  ,fStack(layer.fStack)
  ,fSector(layer.fSector)
  ,fNRows(layer.fNRows)
  ,fN(layer.fN)
  ,fX(layer.fX)
  ,fZ0(layer.fZ0)
  ,fZLength(layer.fZLength)
{
// Copy Constructor 
  
  SetBit(kT0, layer.IsT0());
  SetBit(kOwner, kFALSE);
  for(int i=0; i<kMaxRows; i++) fPositions[i] = layer.fPositions[i];
  memcpy(&fClusters[0], &layer.fClusters[0], kMaxClustersLayer*sizeof(AliTRDcluster*));
  memcpy(&fIndex[0], &layer.fIndex[0], kMaxClustersLayer*sizeof(UInt_t));


// 	BuildIndices();
}

//_____________________________________________________________________________
AliTRDchamberTimeBin &AliTRDchamberTimeBin::operator=(const AliTRDchamberTimeBin &layer)
{
// Assignment operator

  if (this != &layer) layer.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::Clear(const Option_t *) 
{ 
  // Reset the Chamber Timebin
  if(IsOwner())
    for(Int_t it = 0; it<kMaxClustersLayer; it++)
      delete fClusters[it];
  memset(fClusters,0,kMaxClustersLayer*sizeof(fClusters[0]));
  fN = 0; 
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::Copy(TObject &o) const
{
// Copy method. Performs a deep copy of all data from this object to object o.
  
  AliTRDchamberTimeBin &layer = (AliTRDchamberTimeBin &)o;
  layer.fkReconstructor = fkReconstructor;
  layer.fPlane       = fPlane;
  layer.fStack       = fStack;
  layer.fSector      = fSector;
  layer.fNRows       = fNRows;
  layer.fN           = fN;
  layer.fX           = fX;
  layer.fZ0          = fZ0;
  layer.fZLength     = fZLength;
  layer.SetT0(IsT0());
  layer.SetBit(kOwner, kFALSE);
  
  for(int i=0; i<kMaxRows; i++) layer.fPositions[i] = fPositions[i];
  memcpy(&layer.fClusters[0], &fClusters[0], kMaxClustersLayer*sizeof(AliTRDcluster*));
  memcpy(&layer.fIndex[0], &fIndex[0], kMaxClustersLayer*sizeof(UInt_t));
  
  TObject::Copy(layer); // copies everything into layer
  
// 	layer.BuildIndices();
}

//_____________________________________________________________________________
AliTRDchamberTimeBin::~AliTRDchamberTimeBin()
{
// Destructor
  if(IsOwner()){ 
    for(AliTRDcluster **cit = &fClusters[0]; (*cit); cit++) delete (*cit);
  }
}

//_____________________________________________________________________________
void  AliTRDchamberTimeBin::SetOwner(Bool_t copy) 
{
// Sets the ownership of the clusters to this 
// If option "copy" is kTRUE [default] the clusters 
// are also copied otherwise only the ownership bit 
// is flipped.

  SetBit(kOwner, kTRUE);
  if(!copy) return;
  for(AliTRDcluster **cit = &fClusters[0]; (*cit); ++cit){
    (*cit) = new AliTRDcluster(*(*cit)); 
  }
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::SetRange(Float_t z0, Float_t zLength)
{
// Sets the range in z-direction
//
// Parameters:
//   z0      : starting position of layer in the z direction
//   zLength : length of layer in the z direction 

  fZ0 = (z0 <= z0 + zLength) ? z0 : z0 + zLength;
  fZLength = TMath::Abs(zLength);
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::InsertCluster(AliTRDcluster *c, UInt_t index) 
{
  //
  // Insert cluster in cluster array.
  // Clusters are sorted according to Y coordinate.  
  //

  //if (fTimeBinIndex < 0) { 
    //AliWarning("Attempt to insert cluster into non-sensitive time bin!\n");
    //return;
  //}

  if (fN == (Int_t) kMaxClustersLayer) {
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

//___________________________________________________
void AliTRDchamberTimeBin::Bootstrap(const AliTRDReconstructor *rec, Int_t det)
{
// Reinitialize all data members from the clusters array
// It has to be used after reading from disk

  fkReconstructor = rec;
  fPlane  = AliTRDgeometry::GetLayer(det);
  fStack  = AliTRDgeometry::GetStack(det);
  fSector = AliTRDgeometry::GetSector(det);
  AliTRDgeometry g;
  fNRows  = g.GetPadPlane(fPlane, fStack)->GetNrows();
  fN = 0;
  for(AliTRDcluster **cit = &fClusters[0]; (*cit); cit++) fN++;
  BuildIndices();
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::BuildIndices(Int_t iter)
{
// Rearrangement of the clusters belonging to the propagation layer for the stack.
//
// Detailed description
//
// The array indices of all clusters in one PropagationLayer are stored in
// array. The array is divided into several bins.
// The clusters are sorted in increasing order of their y coordinate.
//
// Sorting algorithm: TreeSearch
//

  if(!fN) return;

  // Select clusters that belong to the Stack
  Int_t nClStack = 0;					// Internal counter
  for(Int_t i = 0; i < fN; i++){
    if(fClusters[i]->IsUsed() || fClusters[i]->IsShared()){
      fClusters[i] = NULL;
      fIndex[i] = 0xffff;
    } else nClStack++;
  }
  if(nClStack > kMaxClustersLayer) AliInfo(Form("Number of clusters in stack %d exceed buffer size %d. Truncating.", nClStack, kMaxClustersLayer));
    
  // Nothing in this time bin. Reset indexes 
  if(!nClStack){
    fN = 0;
    memset(&fPositions[0], 0, sizeof(UChar_t) * kMaxRows);
    memset(&fClusters[0], 0, sizeof(AliTRDcluster*) * kMaxClustersLayer);
    memset(&fIndex[0], 0, sizeof(UInt_t) * kMaxClustersLayer);
    return;
  }
  
  // Make a copy
  AliTRDcluster *helpCL[kMaxClustersLayer];
  UInt_t helpInd[kMaxClustersLayer];
  nClStack = 0;
  // Defining iterators
  AliTRDcluster **fcliter = &fClusters[0], **hcliter = &helpCL[0]; UInt_t *finditer = &fIndex[0], *hinditer = &helpInd[0];
  AliTRDcluster *tmpcl = NULL;
  for(Int_t i = 0; i < TMath::Min(fN, kMaxClustersLayer); i++){
    if(!(tmpcl = *(fcliter++))){
    	finditer++;
    	continue;
    }
    *(hcliter++)  = tmpcl;
    *(hinditer++) = *finditer;
    tmpcl = NULL;
    *(finditer++) = 0xffff;
    nClStack++;
  }
  
  // do clusters arrangement
  fX = 0.;
  fN =  nClStack;
  nClStack = 0;
  // Reset Positions array
  memset(fPositions, 0, sizeof(UChar_t)*kMaxRows);
  AliTRDcluster **cliter = &helpCL[0]; // Declare iterator running over the whole array
  const AliTRDrecoParam* const recoParam = fkReconstructor->GetRecoParam(); //the dynamic cast in GetRecoParam is slow, so caching the pointer to it
  Int_t tb(-1);
  for(Int_t i = 0; i < fN; i++){
    // boundary check
    AliTRDcluster *cl = *(cliter++);
    UChar_t rowIndex = cl->GetPadRow();
    if(tb<0) tb=cl->GetLocalTimeBin();
    // Insert Leaf
    Int_t pos = FindYPosition(cl->GetY(), rowIndex, nClStack);
    if(pos == -2) continue;   // zbin error;
    else if(pos == -1) {    // zbin is empty;
      Int_t upper = (rowIndex == fNRows - 1) ? nClStack : fPositions[rowIndex + 1];
      memmove(fClusters + upper + 1, fClusters + upper, (sizeof(AliTRDcluster *))*(nClStack-upper));
      memmove(fIndex + upper + 1, fIndex + upper, (sizeof(UInt_t))*(nClStack-upper));
      fClusters[upper] = cl;
      fIndex[upper] = helpInd[i]; 
      // Move All pointer one position back
      for(UChar_t j = rowIndex + 1; j < fNRows; j++) fPositions[j]++;
      nClStack++;
    } else {		// zbin not empty
      memmove(fClusters + pos + 2, fClusters + pos+1, (sizeof(AliTRDcluster *))*(nClStack-(pos+1)));
      memmove(fIndex + pos + 2, fIndex + pos+1, (sizeof(UInt_t))*(nClStack-(pos+1)));
      fClusters[pos + 1] = cl;	//fIndex[i];
      fIndex[pos + 1] = helpInd[i];
      // Move All pointer one position back
      for(UChar_t j = rowIndex + 1; j < fNRows; j++) fPositions[j]++;	
      nClStack++;
    }

    // calculate mean x
    fX += cl->GetX(); 
    
    // Debug Streaming
    if(recoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 3){
      AliTRDcluster dcl(*cl);
      TTreeSRedirector &cstream = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
      cstream << "BuildIndices"
      << "Plane="    << fPlane
      << "Stack="    << fStack
      << "Sector="   << fSector
      << "Iter="     << iter
      << "C.="       << &dcl
      << "rowIndex=" << rowIndex
      << "\n";
    }
  }

// 	AliInfo("Positions");
// 	for(int ir=0; ir<fNRows; ir++) printf("pos[%d] %d\n", ir, fPositions[ir]);
  if(nClStack < fN){
    AliWarning(Form("Found %d out of %d clusters outside in ChamberTimeBin[%02d_%d_%d|%2d]", fN-nClStack, fN, fSector, fStack, fPlane, tb));
    fN =  nClStack;
    if(!fN){ // Nothing left in this time bin. Reset indexes
      memset(&fPositions[0], 0, sizeof(UChar_t) * kMaxRows);
      memset(&fClusters[0], 0, sizeof(AliTRDcluster*) * kMaxClustersLayer);
      memset(&fIndex[0], 0, sizeof(UInt_t) * kMaxClustersLayer);
      return;
    }
  }
  if(fN) fX /= fN;
}

//_____________________________________________________________________________
Int_t AliTRDchamberTimeBin::Find(Float_t y) const
{
  //
  // Returns index of the cluster nearest in Y    
  //

  if (fN <= 0) return 0;
  
  if (y <= fClusters[0]->GetY()) return 0;
  
  if (y >  fClusters[fN-1]->GetY()) return fN;
  

  Int_t b = 0;
  Int_t e = fN - 1;
  Int_t m = (b + e) / 2;

  for ( ; b < e; m = (b + e) / 2) {
    if (y > fClusters[m]->GetY()) b = m + 1;
    else e = m;
  }

  return m;
}    

//_____________________________________________________________________________
Int_t AliTRDchamberTimeBin::FindYPosition(Double_t y, UChar_t z, Int_t nClusters) const
{
//
// Tree search Algorithm to find the nearest left cluster for a given
// y-position in a certain z-bin (in fact AVL-tree). 
// Making use of the fact that clusters are sorted in y-direction.
//
// Parameters:
//   y : y position of the reference point in tracking coordinates
//   z : z reference bin.
//   nClusters : 
//
// Output :
// Index of the nearest left cluster in the StackLayer indexing (-1 if no clusters are found)
//

  if(z>=fNRows){ // check pad row of cluster
    AliDebug(1, Form("Row[%2d] outside range [0 %2d] in %02d_%d_%d.", z, fNRows, fSector, fStack, fPlane));
    return -2;
  }
  Int_t start = fPositions[z]; 		// starting Position of the bin
  Int_t upper = (Int_t)((z != fNRows - 1) ? fPositions[z+1] : nClusters); // ending Position of the bin 
  Int_t end = upper - 1; // ending Position of the bin 
  if(end < start) return -1; // Bin is empty
  Int_t middle = static_cast<Int_t>((start + end)/2);
  // 1st Part: climb down the tree: get the next cluster BEFORE ypos
  while(start + 1 < end){
    if(y >= fClusters[middle]->GetY()) start = middle;
    else end = middle;
    middle = static_cast<Int_t>((start + end)/2);
  }
  if(y > fClusters[end]->GetY()) return end;
  return start;
}

//_____________________________________________________________________________
Int_t AliTRDchamberTimeBin::FindNearestYCluster(Double_t y, UChar_t z) const
{
//
// Tree search Algorithm to find the nearest cluster for a given
// y-position in a certain z-bin (in fact AVL-tree). 
// Making use of the fact that clusters are sorted in y-direction.
//
// Parameters:
//   y : y position of the reference point in tracking coordinates
//   z : z reference bin.
//
// Output 
// Index of the nearest cluster in the StackLayer indexing (-1 if no clusters are found)
//

  Int_t position = FindYPosition(y, z, fN);
  if(position == -2 || position == -1) return position; // bin empty
  // FindYPosition always returns the left Neighbor. We don't know if the left or the right Neighbor is nearest
  // to the Reference y-position, so test both
  Int_t upper = (Int_t)((z < fNRows-1) ? fPositions[z+1] : fN); // ending Position of the bin
  if((position + 1) < (upper)){
    if(TMath::Abs(y - fClusters[position + 1]->GetY()) < TMath::Abs(y - fClusters[position]->GetY())) return position + 1;
    else return position;
  }
  return position;
}

//_____________________________________________________________________________
Int_t AliTRDchamberTimeBin::SearchNearestCluster(Double_t y, Double_t z, Double_t maxroady, Double_t maxroadz) const
{
//
// Finds the nearest cluster from a given point in a defined range.
// Distance is determined in a 2D space by the 2-Norm.
//
// Parameters :
//   y : y position of the reference point in tracking coordinates
//   z : z reference bin.
//   maxroady : maximum searching distance in y direction
//   maxroadz : maximum searching distance in z direction
//
// Output 
// Index of the nearest cluster in the StackLayer indexing (-1 if no cluster is found).
// Cluster can be accessed with the operator[] or GetCluster(Int_t index)
//
// Detail description
//
// The following steps are perfomed:
// 1. Get the expected z bins inside maxroadz.
// 2. For each z bin find nearest y cluster.
// 3. Select best candidate
//
  Int_t   index   = -1;
  // initial minimal distance will be represented as ellipse: semi-major = z-direction
  // later 2-Norm will be used  
// 	Float_t nExcentricity = TMath::Sqrt(maxroadz*maxroadz - maxroad*maxroad)/maxroadz;
  Float_t mindist = maxroadz;
  
  // not very nice but unfortunately neccessarry: we have ho check the neighbors in both directions (+ and -) too. How 
  // much neighbors depends on the Quotient maxroadz/fZLength   
  UChar_t maxRows = 3;
  UChar_t zpos[kMaxRows];
  // Float_t mindist = TMath::Sqrt(maxroad*maxroad + maxroadz*maxroadz);
// 	UChar_t myZbin = FindTreePosition(z, fZ0 + fZLength/2, fZLength/4, 8, 8, kFALSE);
  UChar_t myZbin = fNRows - 1 - (UChar_t)(TMath::Abs(fZ0 - z)/fZLength * fNRows);
  if(z < fZ0) myZbin = fNRows - 1;
  if(z > fZ0 + fZLength) myZbin = 0;
  //printf("\n%f < %f < %f [%d]\n", fZ0, z, fZ0 + fZLength, myZbin);
  //for(int ic=0; ic<fN; ic++) printf("%d z = %f row %d\n", ic, fClusters[ic]->GetZ(), fClusters[ic]->GetPadRow());

  UChar_t nNeighbors = 0;
  for(UChar_t i = 0; i < maxRows; i++){
    if((myZbin - 1 + i) < 0) continue;
    if((myZbin - 1 + i) > fNRows - 1) break;
    zpos[nNeighbors] = myZbin - 1 + i;
    nNeighbors++;
  }
  Float_t ycl = 0, zcl = 0;
  for(UChar_t neighbor = 0; neighbor < nNeighbors; neighbor++){	// Always test the neighbors too
    Int_t pos = FindNearestYCluster(y, zpos[neighbor]);
    if(pos == -1) continue;	// No cluster in bin
    AliTRDcluster *c = (AliTRDcluster *) (fClusters[pos]);
    if(c->IsUsed()) continue;		// we are only interested in unused clusters
    ycl = c->GetY();
    // Too far away in y-direction (Prearrangement)
    if (TMath::Abs(ycl - y) > maxroady){ 
      //printf("y[%f] ycl[%f] roady[%f]\n", y, ycl, maxroady);
      continue;
    }
    zcl = c->GetZ();
    // Too far away in z-Direction
    // (Prearrangement since we have not so many bins to test)
    if (TMath::Abs(zcl - z) > maxroadz) continue;
    
    Float_t dist;		// distance defined as 2-Norm	
    // if we havent found a Particle that is in the ellipse around (y,z) with maxroad as semi-minor and
    // maxroadz as semi-major, we take the radius of the ellipse concerning the cluster as mindist, later we 
    // take the 2-Norm when we found a cluster inside the ellipse (The value 10000 is taken because it is surely
    // large enough to be usable as an indicator whether we have found a nearer cluster or not)
// 		if(mindist > 10000.){
// 			Float_t phi = ((zcl - z) == 0) ? TMath::Pi()/2 : TMath::ATan((ycl - y)/(zcl - z));
// 			mindist = maxroad/TMath::Sqrt(1 - nExcentricity*nExcentricity * (TMath::Cos(phi))*(TMath::Cos(phi)));
// 		}
    dist = TMath::Max(TMath::Abs(y-ycl),TMath::Abs(z-zcl));	// infinity Norm
// 		dist = TMath::Sqrt((ycl - y)*(ycl - y) + (zcl - z)*(zcl - z));
    if((Int_t)(dist * 100000) < (Int_t)(mindist * 100000)){
    //if((dist = TMath::Sqrt((ycl - y)*(ycl - y) + (zcl - z)*(zcl - z))) < mindist){
      mindist = dist;
      index   = pos;
    }	
  }
  // This is the Array Position in fIndex2D of the Nearest cluster: if a
  // cluster is called, then the function has to retrieve the Information
  // which is Stored in the Array called, the function
  return index;
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::BuildCond(AliTRDcluster * const cl, Double_t *cond, UChar_t Layer, Double_t theta, Double_t phi)
{
// Helper function to calculate the area where to expect a cluster in THIS
// layer. 
//
// Parameters :
//   cl    : 
//   cond  :
//   Layer : 
//   theta : 
//   phi   : 
//
// Detail description
//
// Helper function to calculate the area where to expect a cluster in THIS
// layer. by using the information of a former cluster in another layer
// and the angle in theta- and phi-direction between layer 0 and layer 3.
// If the layer is zero, initial conditions are calculated. Otherwise a
// linear interpolation is performed. 
//Begin_Html
//<img src="gif/build_cond.gif">
//End_Html
//

  if(!fkReconstructor){
    AliError("Reconstructor not set.");
    return;
  }
  
  if(Layer == 0){
    cond[0] = cl->GetY();			// center: y-Direction
    cond[1] = cl->GetZ();			// center: z-Direction
    cond[2] = fkReconstructor->GetRecoParam()->GetMaxPhi()   * (cl->GetX() - GetX()) + 1.0;			// deviation: y-Direction
    cond[3] = fkReconstructor->GetRecoParam()->GetMaxTheta() * (cl->GetX() - GetX()) + 1.0;			// deviation: z-Direction
  } else {
    cond[0] = cl->GetY() + phi   * (GetX() - cl->GetX()); 
    cond[1] = cl->GetZ() + theta * (GetX() - cl->GetX());
    cond[2] = fkReconstructor->GetRecoParam()->GetRoad0y() + phi;
    cond[3] = fkReconstructor->GetRecoParam()->GetRoad0z();
  }
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::GetClusters(const Double_t * const cond, Int_t *index, Int_t& ncl, Int_t BufferSize)
{
// Finds all clusters situated in this layer inside a rectangle  given by the center an ranges.
//
// Parameters :
//   cond  :
//   index : 
//   ncl :
//   BufferSize   :
//
// Output :
//
// Detail description
//
// Function returs an array containing the indices in the stacklayer of
// the clusters found an  the number of found clusters in the stacklayer

  ncl = 0;
  memset(index, 0, BufferSize*sizeof(Int_t));
  if(fN == 0) return;
    
  //Boundary checks
  Double_t zvals[2];
  if(((cond[1] - cond[3]) >= (fZ0 + fZLength)) || (cond[1] + cond[3]) <= fZ0) return; // We are outside of the chamvber
  zvals[0] = ((cond[1] - cond[3]) < fZ0) ? fZ0 : (cond[1] - cond[3]);
  zvals[1] = ((cond[1] + cond[3]) < fZ0 + fZLength) ? (cond[1] + cond[3]) : fZ0 + fZLength - 1.E-3;

  UChar_t zhi = fNRows - 1 - (UChar_t)(TMath::Abs(fZ0 - zvals[0])/fZLength * fNRows);
  UChar_t zlo = fNRows - 1 - (UChar_t)(TMath::Abs(fZ0 - zvals[1])/fZLength * fNRows);

/*	AliInfo(Form("yc[%f] zc[%f] dy[%f] dz[%f]", cond[0], cond[1], cond[2], cond[3]));
  PrintClusters();
  AliInfo(Form("zlo[%f] zhi[%f]", zvals[0], zvals[1]));
  AliInfo(Form("zlo[%d] zhi[%d]", zlo, zhi));*/

  Double_t ylo = cond[0] - cond[2],
           yhi = cond[0] + cond[2];
  //printf("CTB : ylo[%f] yhi[%f]\n", ylo, yhi);
  //Preordering in Direction z saves a lot of loops (boundary checked)
  for(UChar_t z = zlo; z <= zhi; z++){
    UInt_t upper = (z < fNRows-1) ? fPositions[z+1] : fN;
    //AliInfo(Form("z[%d] y [%d %d]", z, fPositions[z], upper));
    for(Int_t y = fPositions[z]; y < (Int_t)upper; y++){
      if(ncl == BufferSize){
        AliDebug(1, Form("Buffer size [%d] riched. Some clusters may be lost.", BufferSize));
        return;	//Buffer filled
      }
      
      if(fClusters[y]->GetY() > yhi) break;			// Abbortion conditions!!!
      if(fClusters[y]->GetY() < ylo) continue;	// Too small
      if(((Int_t)((fClusters[y]->GetZ())*1000) < (Int_t)(zvals[0]*1000)) || ((Int_t)((fClusters[y]->GetZ())*1000) > (Int_t)(zvals[1]*1000))){/*printf("exit z\n"); TODO*/ continue;}
      index[ncl] = y;
      ncl++;
    }
  }
  if(ncl>fN) AliError(Form("Clusters found %d > %d (clusters in layer)", ncl, fN));
}

//_____________________________________________________________________________
AliTRDcluster *AliTRDchamberTimeBin::GetNearestCluster(Double_t *cond)
{
// Function returning a pointer to the nearest cluster (nullpointer if not successfull).
//
// Parameters :
//   cond  :
//
// Output :
//   pointer to the nearest cluster (nullpointer if not successfull).
// 
// Detail description
//
// returns a pointer to the nearest cluster (nullpointer if not
// successfull) by the help of the method FindNearestCluster
  
  
  Double_t maxroad  = fkReconstructor->GetRecoParam()->GetRoad2y();
  Double_t maxroadz = fkReconstructor->GetRecoParam()->GetRoad2z();
  
  Int_t index = SearchNearestCluster(cond[0],cond[1],maxroad,maxroadz);
  AliTRDcluster *returnCluster = NULL;
  if(index != -1) returnCluster = (AliTRDcluster *) fClusters[index];
  return returnCluster;
}

//_____________________________________________________________________________
void AliTRDchamberTimeBin::Print(Option_t *) const
{
// Prints the position of each cluster in the stacklayer on the stdout
//
  if(!fN) return;
  AliInfo(Form("Layer[%d] Stack[%d] Sector[%2d] nRows[%2d]", fPlane, fStack, fSector, fNRows));
  AliInfo(Form("Z0[%7.3f] Z1[%7.3f]", fZ0, fZ0+fZLength));
  AliTRDcluster* const* icl = &fClusters[0];
  for(Int_t jcl = 0; jcl < fN; jcl++, icl++){
    AliInfo(Form("%2d X[%7.3f] Y[%7.3f] Z[%7.3f] tb[%2d] col[%3d] row[%2d] used[%s]", jcl,  (*icl)->GetX(), (*icl)->GetY(), (*icl)->GetZ(), (*icl)->GetLocalTimeBin(), (*icl)->GetPadCol(), (*icl)->GetPadRow(),
    (*icl)->IsUsed() ? "y" : "n"));
  }
  Int_t irow = 0;
  for(UChar_t const* pos = &fPositions[0]; irow<fNRows; pos++, irow++){ 
    if((*pos) != 0xff) AliInfo(Form("r[%2d] pos[%d]", irow, (*pos)));
  }
}
