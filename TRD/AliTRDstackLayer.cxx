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
#include <TROOT.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTreeStream.h>

#include "AliLog.h"
#include "AliTRDcluster.h"

#include "AliTRDstackLayer.h"
#include "AliTRDrecoParam.h"
#include "AliTRDReconstructor.h"

#define DEBUG

ClassImp(AliTRDstackLayer)

//_____________________________________________________________________________
AliTRDstackLayer::AliTRDstackLayer(Double_t z0, Double_t zLength
                                 , UChar_t stackNr, AliTRDrecoParam *p)
  :AliTRDpropagationLayer()
  ,fOwner(kFALSE)
  ,fStackNr(stackNr)
  ,fNRows(kMaxRows)
  ,fZ0(z0)
  ,fZLength(zLength)
  ,fRecoParam(p)
  ,fDebugStream(0x0)
{
  //
  // Default constructor (Only provided to use AliTRDstackLayer with arrays)
  //

	for(int i=0; i<kMaxRows; i++) fPositions[i] = 0;
}

//_____________________________________________________________________________
AliTRDstackLayer::AliTRDstackLayer(const AliTRDpropagationLayer &layer, Double_t
z0, Double_t zLength, UChar_t stackNr, AliTRDrecoParam *p):
	AliTRDpropagationLayer(layer)
	,fOwner(kFALSE)
	,fStackNr(stackNr)
	,fNRows(kMaxRows)
	,fZ0(z0)
	,fZLength(zLength)
	,fRecoParam(p)
	,fDebugStream(0x0)
{
// Standard constructor.
// Initialize also the underlying AliTRDpropagationLayer using the copy constructor.

	for(int i=0; i<kMaxRows; i++) fPositions[i] = 0;
}

//_____________________________________________________________________________
AliTRDstackLayer::AliTRDstackLayer(const AliTRDpropagationLayer &layer):
	AliTRDpropagationLayer(layer)
	,fOwner(kFALSE)
	,fStackNr(0)
	,fNRows(kMaxRows)
	,fZ0(0)
	,fZLength(0)
	,fRecoParam(0x0)
	,fDebugStream(0x0)
{
// Standard constructor using only AliTRDpropagationLayer.
	
	for(int i=0; i<kMaxRows; i++) fPositions[i] = 0;
}

//_____________________________________________________________________________
AliTRDstackLayer::AliTRDstackLayer(const AliTRDstackLayer &layer):
	AliTRDpropagationLayer(layer)
	,fOwner(layer.fOwner)
	,fStackNr(layer.fStackNr)
	,fNRows(layer.fNRows)
	,fZ0(layer.fZ0)
	,fZLength(layer.fZLength)
	,fRecoParam(layer.fRecoParam)
	,fDebugStream(layer.fDebugStream)
{
// Copy Constructor (performs a deep copy)
	
	for(Int_t i = 0; i < kMaxRows; i++) fPositions[i] = layer.fPositions[i];
// 	BuildIndices();
}

//_____________________________________________________________________________
AliTRDstackLayer &AliTRDstackLayer::operator=(const AliTRDpropagationLayer &layer)
{
// Assignment operator from an AliTRDpropagationLayer

	if (this != &layer) layer.Copy(*this);
	return *this;
}

//_____________________________________________________________________________
AliTRDstackLayer &AliTRDstackLayer::operator=(const AliTRDstackLayer &layer)
{
// Assignment operator

	if (this != &layer) layer.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliTRDstackLayer::Copy(TObject &o) const
{
// Copy method. Performs a deep copy of all data from this object to object o.
	
	AliTRDstackLayer &layer = (AliTRDstackLayer &)o;
	layer.fZ0          = fZ0;
	layer.fOwner       = kFALSE;
	layer.fNRows       = fNRows;
	layer.fZLength     = fZLength;
	layer.fStackNr     = fStackNr;
	layer.fDebugStream = fDebugStream;
	layer.fRecoParam   = fRecoParam;
	
	AliTRDpropagationLayer::Copy(layer); // copies everything into layer
	for(UChar_t i = 0; i < kMaxRows; i++) layer.fPositions[i] = 0;
// 	layer.BuildIndices();
}

//_____________________________________________________________________________
AliTRDstackLayer::~AliTRDstackLayer()
{
// Destructor
	if(fOwner) for(int ic=0; ic<fN; ic++) delete fClusters[ic];
}

//_____________________________________________________________________________
void AliTRDstackLayer::SetRange(const Float_t z0, const Float_t zLength)
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
void AliTRDstackLayer::BuildIndices(Int_t iter)
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
		Double_t zval = fClusters[i]->GetZ();
		if(zval < fZ0 || zval > fZ0 + fZLength || fClusters[i]->IsUsed()){
			fClusters[i] = 0x0;
			fIndex[i] = 9999;
		} else nClStack++;
	}
	if(nClStack > kMaxClustersLayer) AliWarning(Form("Number of clusters in stack %d exceed buffer size %d", nClStack, kMaxClustersLayer));
	
	// Nothing in this Stack 
	if(!nClStack){
		delete fClusters;
		delete fIndex;
		fClusters = 0x0;
		fIndex = 0x0;
		fN = 0;
		memset(fPositions, 0, sizeof(UChar_t) * 16);
		return;
	}
	
	// Make a copy
	AliTRDcluster *helpCL[kMaxClustersLayer];
	Int_t helpInd[kMaxClustersLayer];
	nClStack = 0;
	for(Int_t i = 0; i < TMath::Min(fN, kMaxClustersLayer); i++){
		if(!fClusters[i]) continue;
		helpCL[nClStack]  = fClusters[i];
		helpInd[nClStack] = fIndex[i];
		fClusters[i]      = 0x0;
		fIndex[i]         = 9999;
		nClStack++;
	}
	
	// do clusters arrangement
	fN =  nClStack;
	nClStack = 0;
	memset(fPositions,0,sizeof(UChar_t)*16);				// Reset Positions array
	for(Int_t i = 0; i < fN; i++){
		// boundarie check
		AliTRDcluster *cl = helpCL[i];
		Double_t zval = cl->GetZ();
		UChar_t treeIndex = (UChar_t)(TMath::Abs(fZ0 - zval)/fZLength * fNRows);
		if(treeIndex > fNRows - 1) treeIndex = fNRows - 1;
		// Insert Leaf
		Int_t pos = FindYPosition(cl->GetY(), treeIndex, i);
		if(pos == -1){		// zbin is empty;
			Int_t upper = (treeIndex == fNRows - 1) ? nClStack : fPositions[treeIndex + 1];
			memmove(fClusters + upper + 1, fClusters + upper, (sizeof(AliTRDcluster *))*(nClStack-upper));
			memmove(fIndex + upper + 1, fIndex + upper, (sizeof(UInt_t))*(nClStack-upper));
			fClusters[upper] = cl;
			fIndex[upper] = helpInd[i]; 
			// Move All pointer one position back
			for(UChar_t j = treeIndex + 1; j < fNRows; j++) fPositions[j]++;
			nClStack++;
		} else {		// zbin not empty
			memmove(fClusters + pos + 2, fClusters + pos+1, (sizeof(AliTRDcluster *))*(nClStack-(pos+1)));
			memmove(fIndex + pos + 2, fIndex + pos+1, (sizeof(UInt_t))*(nClStack-(pos+1)));
			fClusters[pos + 1] = cl;	//fIndex[i];
			fIndex[pos + 1] = helpInd[i];
			// Move All pointer one position back
			for(UChar_t j = treeIndex + 1; j < fNRows; j++) fPositions[j]++;	
			nClStack++;
		}

		
		// Debug Streaming
#ifdef DEBUG
		if(fDebugStream && AliTRDReconstructor::StreamLevel() >= 3){
			TTreeSRedirector &cstream = *fDebugStream;
			cstream << "BuildIndices"
			<< "StackNr="  << fStackNr
			<< "SectorNr=" << fSec
			<< "Iter="     << iter
			<< "C.="       << cl
			<< "TreeIdx="  << treeIndex
			<< "\n";
		}
#endif	
	}
}

//_____________________________________________________________________________
Int_t AliTRDstackLayer::FindYPosition(Double_t y, UChar_t z, Int_t nClusters) const
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
Int_t AliTRDstackLayer::FindNearestYCluster(Double_t y, UChar_t z) const
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
	if(position == -1) return position; // bin empty
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
Int_t AliTRDstackLayer::SearchNearestCluster(Double_t y, Double_t z, Double_t maxroady, Double_t maxroadz) const
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
	UChar_t myZbin = (UChar_t)(TMath::Abs(fZ0 - z)/fZLength * fNRows);
	if(z < fZ0) myZbin = 0;
	if(z > fZ0 + fZLength) myZbin = fNRows - 1;

	UChar_t nNeighbors = 0;
	for(UChar_t i = 0; i < maxRows; i++){
		if((myZbin - 1 + i) < 0) continue;
		if((myZbin - 1 + i) > fNRows - 1) break;
		zpos[nNeighbors] = myZbin - 1 + i;
		nNeighbors++;
	}
	Float_t ycl = 0, zcl = 0;
	for(UChar_t neighbor = 0; neighbor < nNeighbors; neighbor++){	// Always test the neighbors too
		Int_t pos = FindNearestYCluster(y,zpos[neighbor]);
		if(pos == -1) continue;	// No cluster in bin
		AliTRDcluster *c = (AliTRDcluster *) (fClusters[pos]);
		if(c->IsUsed()) continue;		// we are only interested in unused clusters
		ycl = c->GetY();
		// Too far away in y-direction (Prearrangement)
		if (TMath::Abs(ycl - y) > maxroady) continue;

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
void AliTRDstackLayer::BuildCond(AliTRDcluster *cl, Double_t *cond, UChar_t Layer, Double_t theta, Double_t phi)
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

	if(!fRecoParam){
		AliError("Reconstruction parameters not initialized.");
		return;
	}
	
	if(Layer == 0){
		cond[0] = cl->GetY();			// center: y-Direction
		cond[1] = cl->GetZ();			// center: z-Direction
		cond[2] = fRecoParam->GetMaxPhi()   * (cl->GetX() - GetX()) + 1.0;			// deviation: y-Direction
		cond[3] = fRecoParam->GetMaxTheta() * (cl->GetX() - GetX()) + 1.0;			// deviation: z-Direction
	} else {
		cond[0] = cl->GetY() + phi   * (GetX() - cl->GetX()); 
		cond[1] = cl->GetZ() + theta * (GetX() - cl->GetX());
		cond[2] = fRecoParam->GetRoad0y() + phi;
		cond[3] = fRecoParam->GetRoad0z();
	}
}

//_____________________________________________________________________________
void AliTRDstackLayer::GetClusters(Double_t *cond, Int_t *index, Int_t& ncl, Int_t BufferSize)
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
	zvals[0] = ((cond[1] - cond[3]) < fZ0) ? fZ0 : (cond[1] - cond[3]);
	zvals[1] = ((cond[1] + cond[3]) < fZ0 + fZLength) ? (cond[1] + cond[3]) : fZ0 + fZLength;

	UChar_t zlo = (fZ0>zvals[0]) ? 0 : (UChar_t)(TMath::Abs(fZ0 - zvals[0])/fZLength * fNRows);
	UChar_t zhi = (fZ0+fZLength<zvals[1]) ? fNRows - 1 : (UChar_t)(TMath::Abs(fZ0 - zvals[1])/fZLength * fNRows);

	//Preordering in Direction z saves a lot of loops (boundary checked)
	for(UChar_t z = zlo; z <= zhi; z++){
		UInt_t upper = (z != fNRows-1) ? fPositions[z+1] : fN;
		for(Int_t y = fPositions[z]; y < (Int_t)upper; y++){
			if(ncl == BufferSize){
				AliWarning("Buffer size riched. Some clusters may be lost.");
				return;	//Buffer filled
			}
			if(fClusters[y]->GetY() > (cond[0] + cond[2])) break;			// Abbortion conditions!!!
			if(fClusters[y]->GetY() < (cond[0] - cond[2])) continue;	// Too small
			if(((Int_t)((fClusters[y]->GetZ())*1000) < (Int_t)(zvals[0]*1000)) || ((Int_t)((fClusters[y]->GetZ())*1000) > (Int_t)(zvals[1]*1000))) continue;
			index[ncl] = y;
			ncl++;
		}
	}
	if(ncl>fN) AliError(Form("Clusters found %d > %d (clusters in layer)", ncl, fN));
}

//_____________________________________________________________________________
AliTRDcluster *AliTRDstackLayer::GetNearestCluster(Double_t *cond)
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
	
	
	Double_t maxroad  = fRecoParam->GetRoad2y();
	Double_t maxroadz = fRecoParam->GetRoad2z();
	
	Int_t index = SearchNearestCluster(cond[0],cond[1],maxroad,maxroadz);
	AliTRDcluster *returnCluster = 0x0;
	if(index != -1) returnCluster = (AliTRDcluster *) fClusters[index];
	return returnCluster;
}

//_____________________________________________________________________________
void AliTRDstackLayer::PrintClusters() const
{
// Prints the position of each cluster in the stacklayer on the stdout
//
	printf("fDebugStream = %#o\n", ((Int_t) fDebugStream));
	printf("fRecoParam   = %#o\n", ((Int_t) fRecoParam));
	
	for(Int_t i = 0; i < fN; i++){
		printf("AliTRDstackLayer: index=%i, Cluster: X = %3.3f, Y = %3.3f, Z = %3.3f\n", i,  fClusters[i]->GetX(),fClusters[i]->GetY(),fClusters[i]->GetZ());
		if(fClusters[i]->IsUsed()) printf("cluster allready used. rejected in search algorithm\n");
	}
}
