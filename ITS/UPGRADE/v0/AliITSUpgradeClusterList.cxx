/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////
// Author: A.Mastroserio, C.Terrevoli                               //
//         annalisa.mastroserio@cern.ch			  	    //
//         cristina.terrevoli@ba.infn.it			    //
// This class implements the use of a list of clusters.             //
//////////////////////////////////////////////////////////////////////  

/* $Id$ */

#include "AliITSUpgradeClusterList.h"
#include "AliITSUpgradeClusterListNode.h"

//______________________________________________________________________________
AliITSUpgradeClusterList::AliITSUpgradeClusterList():
  fNrEntries(0),
  fFirst(NULL),
  fLast(NULL),
  fFastAccess(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{}
//______________________________________________________________________________
AliITSUpgradeClusterList::AliITSUpgradeClusterList(AliITSUpgradeClusterListNode* first, UInt_t nrEntries):
  fNrEntries(nrEntries),
  fFirst(first),
  fLast(NULL),
  fFastAccess(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{
  SetLastNode();
}
//______________________________________________________________________________
AliITSUpgradeClusterList::AliITSUpgradeClusterList(const AliITSUpgradeClusterList& ilist):
  fNrEntries(0),
  fFirst(NULL),
  fLast(NULL),
  fFastAccess(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{
  // copy constructor
  *this = ilist;
}
//______________________________________________________________________________
AliITSUpgradeClusterList::~AliITSUpgradeClusterList() {
  Clear();
}
//______________________________________________________________________________
AliITSUpgradeClusterList& AliITSUpgradeClusterList::operator=(const AliITSUpgradeClusterList& ilist) {
  // assignment operator
  if (this!=&ilist) {
    this->Clear();
    fFirst = CloneNode(ilist.fFirst);
    SetLastNode();
    fFastAccess=kFALSE;
    fFastAccessArray=NULL;
    fDummyIndex=0;
  }
  return *this;
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::Clear() {
  // clear the whole list
  ClearFastAccess();
  ClearNode(fFirst);
  fLast=fFirst;
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::ClearNode(AliITSUpgradeClusterListNode* &node) {
  // clear this node and all children nodes
  if (node==NULL) return;
  ClearNode(node->Next());
  delete node;
  fNrEntries--;
  node = NULL;
  fFastAccess=kFALSE;
}
//______________________________________________________________________________
AliITSUpgradeClusterList* AliITSUpgradeClusterList::Clone() const {
  // returns a clone of the list
 AliITSUpgradeClusterListNode *newFirst;
  newFirst = CloneNode(fFirst);
  AliITSUpgradeClusterList* newList = new AliITSUpgradeClusterList(newFirst,fNrEntries);
  return newList;
}
//______________________________________________________________________________
AliITSUpgradeClusterListNode* AliITSUpgradeClusterList::CloneNode(AliITSUpgradeClusterListNode* node) const {
  if (node==NULL) return NULL;
  else return new AliITSUpgradeClusterListNode(node->Col(),node->Row(),node->Size(),node->WidthZ(),node->WidthPhi(),node->Type(),node->Charge(),CloneNode(node->Next()));
}
//______________________________________________________________________________
Bool_t AliITSUpgradeClusterList::Insert(Float_t col, Float_t row, UShort_t size, UShort_t widthZ, UShort_t widthPhi, UShort_t type, UShort_t charge, Int_t digLabels[12*kMaxLab]) {
  // insert a new node into the list (returns true if the node was not present before)
  fNrEntries++;
  AliITSUpgradeClusterListNode* node = new AliITSUpgradeClusterListNode(col,row,size,widthZ,widthPhi,type,charge,NULL);
  for(Int_t i=0; i< 12*kMaxLab; i++) node->AddDigitLabel(digLabels[i]); // adding digit label to the cluster
  if (fFirst==NULL) {
    fFirst = node;
  }
  else {
    fLast->Next() = node;
  }
  fLast = node;
  return kTRUE;
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::SetLastNode() {
  AliITSUpgradeClusterListNode* node = fFirst;
  if (node==NULL) {
    fLast = fFirst;
  }
  else {
    while (1) {
      if (node->Next()==NULL) {
	fLast = node;
	break;
      }
    }
  }
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::ClearFastAccess(){
  // clears the fast access array of pointers
  if (fFastAccessArray!=NULL) {
    delete [] fFastAccessArray;
    fFastAccessArray=NULL;
  }
  fFastAccess=kFALSE;
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::InitFastAccess(){
  // initializes the fast access array
  if (fFastAccess) return;
  ClearFastAccess();
  if (fNrEntries>0) {
    fFastAccessArray = new AliITSUpgradeClusterListNode*[fNrEntries];
    fDummyIndex=0;
    InitFastAccessNode(fFirst);
    fFastAccess=kTRUE;
  }
}
//______________________________________________________________________________
void AliITSUpgradeClusterList::InitFastAccessNode(AliITSUpgradeClusterListNode* node) {
  // initializes the fast access array starting from node (used recursively)
  if (node==NULL) return;
  fFastAccessArray[fDummyIndex++] = node;
  InitFastAccessNode(node->Next());
}
//______________________________________________________________________________
Float_t AliITSUpgradeClusterList::GetColIndex(UInt_t index) {
  // returns the col of the node at position 'index' in the list
  // returns -1 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->Col();
  }
  return -1;
}
//______________________________________________________________________________
Float_t AliITSUpgradeClusterList::GetRowIndex(UInt_t index) {
  // returns the row of the node at position 'index' in the list
  // returns -1 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->Row();
  }
  return -1;
}
//______________________________________________________________________________
UShort_t AliITSUpgradeClusterList::GetSizeIndex(UInt_t index) {
  // returns the size of the node at position 'index' in the list
  // returns 0 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->Size();
  }
  return 0;
}
//______________________________________________________________________________
UShort_t AliITSUpgradeClusterList::GetWidthZIndex(UInt_t index) {
  // returns the width z of the node at position 'index' in the list
  // returns 0 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->WidthZ();
  }
  return 0;
}
//______________________________________________________________________________
UShort_t AliITSUpgradeClusterList::GetWidthPhiIndex(UInt_t index) {
  // returns the width phi of the node at position 'index' in the list
  // returns 0 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->WidthPhi();
  }
  return 0;
}
//______________________________________________________________________________
UShort_t AliITSUpgradeClusterList::GetTypeIndex(UInt_t index) {
  // returns the type of the node at position 'index' in the list
  // returns 99 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess) InitFastAccess();
    return fFastAccessArray[index]->Type();
  }
  return 99;
}
//______________________________________________________________________________
UShort_t AliITSUpgradeClusterList::GetCharge(UInt_t index) {
  // returns the charge of the node at position 'index' in the list
  // returns 0 if out of bounds
  if (index<fNrEntries) {
  if (!fFastAccess) InitFastAccess();
   return fFastAccessArray[index]->Charge();
 }
  return 0;
}
//______________________________________________________________________________
Int_t * AliITSUpgradeClusterList::GetLabels(UInt_t index) {
  // returns the charge of the node at position 'index' in the list
  // returns 0 if out of bounds
  if (index<fNrEntries) {
  if (!fFastAccess) InitFastAccess();
   return fFastAccessArray[index]->GetLabels();
 }
  return 0;
}








