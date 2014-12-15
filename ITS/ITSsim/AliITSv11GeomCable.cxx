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



// General Root includes
//#include <Riostream.h>
//#include <TMath.h>
#include <TVectorD.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>

#include "AliITSv11GeomCable.h"



//*************************************************************************
//   Base class of cable classes
//
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

ClassImp(AliITSv11GeomCable)


//________________________________________________________________________
AliITSv11GeomCable::AliITSv11GeomCable(): TNamed(),
  fDebug(0),
  fPointArray(),
  fVolumeArray(),
  fCurrentVol(0),
  fInitialNode(0)
{ 
  // constructor
  fPointArray.SetOwner();
  for(Int_t i=0;i<fgkCableMaxNodeLevel;i++)fNodeInd[i]=0;
}

//________________________________________________________________________
AliITSv11GeomCable::AliITSv11GeomCable(const char* name): TNamed(name,""),
  fDebug(0),
  fPointArray(),
  fVolumeArray(),
  fCurrentVol(0),
  fInitialNode(0) { 
  // constructor
  fPointArray.SetOwner(); 
  for(Int_t i=0;i<fgkCableMaxNodeLevel;i++)fNodeInd[i]=0;
}


//________________________________________________________________________
AliITSv11GeomCable::~AliITSv11GeomCable() {
  fPointArray.Clear();
  fVolumeArray.Clear();
}

//________________________________________________________________________
void AliITSv11GeomCable::AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
					Double_t *coord)
{
  //
  // Add a check point and its volume container to the cable
  //

  if (iCheckPt>=fVolumeArray.GetEntriesFast()) {
    fVolumeArray.AddLast(vol);
    TVectorD *point = new TVectorD(3,coord);
    fPointArray.AddLast(point);

  } else if ((iCheckPt >= 0)&&(iCheckPt < fVolumeArray.GetEntriesFast())) {
    fVolumeArray.AddAt(vol, iCheckPt);
    TVectorD *point = new TVectorD(3,coord);
    fPointArray.AddAt(point, iCheckPt);
  };
}

//________________________________________________________________________
void AliITSv11GeomCable::ResetPoints() {
  //
  // Remove all points to the cable
  //
  fPointArray.Delete();
  fVolumeArray.Clear();
}


//________________________________________________________________________
Int_t AliITSv11GeomCable::GetPoint( Int_t iCheckPt,  Double_t *coord)
const {
  //
  // Get the check point #iCheckPt
  //
  TVectorD *coordVector =(TVectorD *)fPointArray.UncheckedAt(iCheckPt);
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
  CopyFrom(coord, coordVector->GetElements());
#else
  CopyFrom(coord, coordVector->GetMatrixArray());
#endif
  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSv11GeomCable::GetVect( Int_t iCheckPt,  Double_t *coord)
const {
  //
  //  Get the orientation vect. related to check point #iCheckPt
  //

  TVectorD *coordVector =(TVectorD *)fPointArray.UncheckedAt(iCheckPt);
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
  CopyFrom(coord, coordVector->GetElements());
#else
  CopyFrom(coord, coordVector->GetMatrixArray());
#endif
  return kTRUE;
}

//________________________________________________________________________
TGeoVolume *AliITSv11GeomCable::GetVolume( Int_t iCheckPt ) const {
  //
  // Get the volume of check point #iCheckPt
  //

  if (iCheckPt >= fVolumeArray.GetEntriesFast())
    return 0;
  else
    return (TGeoVolume *) fVolumeArray.UncheckedAt(iCheckPt);
}

//________________________________________________________________________
void AliITSv11GeomCable::SetInitialNode(TGeoVolume *vol) {
  //
  // Set the starting node, initializing the search for the volume
  // containing the cable check point
  //
  if (fInitialNode) delete fInitialNode;
  fInitialNode = new TGeoNodeMatrix(vol,0);
  fInitialNode->SetName("nodeInConstruction");
}

//________________________________________________________________________
void AliITSv11GeomCable::ResetInitialNode() {
  // Reset the initial node if it is set.
  if (fInitialNode) delete fInitialNode;
  fInitialNode = 0;
}

//________________________________________________________________________
bool AliITSv11GeomCable::CheckDaughter(const TGeoNode* node, Int_t i)
{
// Search where is the current volume in the tree of nodes
// stop each time it find the pointer of the current volume
// the path is recorded in fNodeInd[]
// node is the node where the search start.
// !!! recursive function !!!

  Int_t j = fNodeInd[i];
  if (node->GetVolume()==fCurrentVol) return kTRUE;
  TObjArray *array = node->GetNodes();
  if (array) {
    Int_t nDaughters = array->GetEntriesFast();
    if (j==-1) j++;
    while (j<nDaughters) {
      TGeoNode *subNode = (TGeoNode *) array->UncheckedAt(j);
      fNodeInd[i] = j;
      if (CheckDaughter(subNode, i+1)) return kTRUE;
      j++;
    };
    fNodeInd[i] = -1;
  };
  return kFALSE;
}

//________________________________________________________________________
Int_t AliITSv11GeomCable::
GetCheckPoint( Int_t iCheckPt, Int_t iOccur, Int_t motherLevel,
	       Double_t *coord ) {
// Get the coordinate of the check point number #iCheckPt, which is in the
// #iOccur occurrence of the containing volume in the node tree. Coordinates
// are given in the coordinate system of the #motherLevel mother level of
// this volume
  
  if (iCheckPt >= fVolumeArray.GetEntriesFast()) return kFALSE;

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  fCurrentVol = GetVolume(iCheckPt);
  ResetCheckDaughter();
  Int_t currentOccur = 0;

  // loop to get the volume position in the tree of nodes
  while ( CheckDaughter(mainNode) && (currentOccur!=iOccur) ) {
    currentOccur++;
    Int_t maxLevel = 0;
    while (fNodeInd[maxLevel]!=-1) maxLevel++;
    fNodeInd[maxLevel-1]++;
  };

  Int_t maxLevel = 0;
  while (fNodeInd[maxLevel]!=-1) maxLevel++;
  maxLevel--;
  if (maxLevel<-1) return kFALSE;

  TGeoNode *pathNode[fgkCableMaxNodeLevel];
  pathNode[0] = mainNode;
  for (Int_t i=0; i<=maxLevel; i++) {
    pathNode[i+1] = pathNode[i]->GetDaughter(fNodeInd[i]);
  };

  Double_t localCoord[3];
  GetPoint(iCheckPt, localCoord);
  CopyFrom(coord, localCoord);

  if (motherLevel>maxLevel+2) motherLevel = maxLevel+2;

  for (Int_t i=maxLevel; i>maxLevel-motherLevel; i--) {
    pathNode[i+1]->GetMatrix()->LocalToMaster(localCoord, coord);
    CopyFrom(localCoord, coord);
  };
  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSv11GeomCable::GetCheckVect( Int_t iCheckPt, Int_t iOccur,
					Int_t motherLevel, Double_t *coord)
{
// same as GetCheckPoint but with vectorial transformation ...

  if (iCheckPt >= fVolumeArray.GetEntriesFast()) return kFALSE;

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  fCurrentVol = GetVolume(iCheckPt);
  ResetCheckDaughter();
  Int_t currentOccur = 0;

  // loop to get the volume position in the tree of nodes
  while ( CheckDaughter(mainNode) && (currentOccur!=iOccur) ) {
    currentOccur++;
    Int_t maxLevel = 0;
    while (fNodeInd[maxLevel]!=-1) maxLevel++;
    fNodeInd[maxLevel-1]++;
  };

  Int_t maxLevel = 0;
  while (fNodeInd[maxLevel]!=-1) maxLevel++;
  maxLevel--;
  if (maxLevel<-1) return kFALSE;

  TGeoNode *pathNode[fgkCableMaxNodeLevel];
  pathNode[0] = mainNode;
  for (Int_t i=0; i<=maxLevel; i++) {
    pathNode[i+1] = pathNode[i]->GetDaughter(fNodeInd[i]);
  };

  Double_t localCoord[3];
  GetVect(iCheckPt, localCoord);
  CopyFrom(coord, localCoord);

  if (motherLevel>maxLevel+2) motherLevel = maxLevel+2;

  for (Int_t i=maxLevel; i>maxLevel-motherLevel; i--) {
    pathNode[i+1]->GetMatrix()->LocalToMasterVect(localCoord, coord);
    CopyFrom(localCoord, coord);
  };
  return kTRUE;
}


//________________________________________________________________________
Int_t AliITSv11GeomCable::GetCheckVect( const Double_t *localCoord,
					TGeoVolume *vol, Int_t iOccur,
					Int_t motherLevel, Double_t *coord)
{
  //
  // Get the global vect (in coord) correponding to the local vector (localCoord)
  // of the volume vol. Global at the level of #motherLevel level in the node tree
  //

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode; };

  fCurrentVol = vol;
  Int_t currentOccur = 0;

  // loop to get the volume position in the tree of nodes
  ResetCheckDaughter();
  while ( CheckDaughter(mainNode) && (currentOccur!=iOccur) ) {
    currentOccur++;
    Int_t maxLevel = 0;
    while (fNodeInd[maxLevel]!=-1) maxLevel++;
    fNodeInd[maxLevel-1]++;
  };

  Int_t maxLevel = 0;
  while (fNodeInd[maxLevel]!=-1) maxLevel++;
  maxLevel--;
  if (maxLevel<-1) return kFALSE;

  TGeoNode *pathNode[fgkCableMaxNodeLevel];
  pathNode[0] = mainNode;
  for (Int_t i=0; i<=maxLevel; i++) {
    pathNode[i+1] = pathNode[i]->GetDaughter(fNodeInd[i]);
  };

  if (motherLevel>maxLevel+2) motherLevel = maxLevel+2;

  Double_t tempCoord[3] = {localCoord[0], localCoord[1], localCoord[2]};
  CopyFrom(coord, tempCoord);
  for (Int_t i=maxLevel; i>maxLevel-motherLevel; i--) {
    pathNode[i+1]->GetMatrix()->LocalToMasterVect(tempCoord, coord);
    CopyFrom(tempCoord, coord);
  };
  return kTRUE;
}

