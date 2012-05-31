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

//***********************************************************************
//
// It consist of a TClonesArray of 
// AliITSpListItem objects
// This array can be accessed via 2 indexed
// it is used at digitization level by 
// all the 3 ITS subdetectors
//
// ***********************************************************************

#include "AliITSpList.h"
#include "AliITSpListItem.h"


//______________________________________________________________________

ClassImp(AliITSpList)
//______________________________________________________________________
AliITSpList::AliITSpList():
fNi(0),
fNj(0),
fa(0),
fEntries(0){
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A zeroed/empty AliITSpList class.

}
//______________________________________________________________________
AliITSpList::AliITSpList(Int_t imax,Int_t jmax):
fNi(imax),
fNj(jmax),
fa(new AliITSpListItem[imax*jmax]),
fEntries(0){
    // Standard constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A setup AliITSpList class.

}
//______________________________________________________________________
AliITSpList::~AliITSpList(){
    // Default destructor

  delete [] fa;
  fNi = 0;
  fNj = 0;
  fEntries = 0;
}

//______________________________________________________________________
void AliITSpList::ClearMap(){
    // Delete all AliITSpListItems and zero TClonesArray.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A zeroed AliITSpList class.
  for(Int_t i=0; i<fEntries; i++)(fa[i]).MarkUnused();
  fEntries = 0;
}
//______________________________________________________________________
void AliITSpList::DeleteHit(Int_t i,Int_t j){
    // Delete a particular AliITSpListItems.
    // Inputs:
    //    Int_t i   Row number
    //    Int_t j   Columns number
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t k = GetIndex(i,j);
    if((fa[k]).IsUsed())(fa[k]).MarkUnused();
    if(k==fEntries-1) fEntries--;
}
//______________________________________________________________________
AliITSpList& AliITSpList::operator=(const AliITSpList &source){
    // = operator
    // Inputs:
    //    const AliITSpList &source    A AliITSpList object.
    // Outputs:
    //    none.
    // Return:
    //    A copied AliITSpList object.

  this->~AliITSpList();
  new(this) AliITSpList(source);
  return *this;
}
//______________________________________________________________________
AliITSpList::AliITSpList(const AliITSpList &source) : AliITSMap(source),
fNi(source.fNi),
fNj(source.fNj),
fa(new AliITSpListItem[fNi*fNj]),
fEntries(source.fEntries){
    // Copy constructor
  for(Int_t i=0; i<fEntries; i++)(fa[i]).Build(source.fa[i]);
}
//______________________________________________________________________
void AliITSpList::AddItemTo(Int_t fileIndex, AliITSpListItem *pl) {
    // Adds the contents of pl to the list with track number off set given by
    // fileIndex.
    // Creates the AliITSpListItem if needed.
    // Inputs:
    //    Int_t fileIndex      track number offset value
    //    AliITSpListItem *pl  an AliITSpListItem to be added to this class.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t index = pl->GetIndex();
    AliITSpListItem &lit = fa[index];
    if(!lit.IsUsed())lit.Build(-2,-1,pl->GetModule(),index,0.);
    lit.AddTo(fileIndex,pl);
    if(index>=fEntries) fEntries = index +1;
}
//______________________________________________________________________
void AliITSpList::AddSignal(Int_t i,Int_t j,Int_t trk,Int_t ht,Int_t mod,
                       Double_t signal){
    // Adds a Signal value to the TClonesArray at i,j. 
    // Creates the AliITSpListItem
    // if needed.
    // Inputs:
    //    Int_t i         Row number for this signal
    //    Int_t j         Column number for this signal
    //    Int_t trk       Track number creating this signal
    //    Int_t ht        Hit number creating this signal
    //    Int_t mod       The module where this signal is in
    //    Double_t signal The signal (ionization)
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t index = GetIndex(i,j);
    if (index<0) return;
    AliITSpListItem &lit = fa[index];
    if(!lit.IsUsed()){
      lit.Build(trk,ht,mod,index,signal);
    }
    else {
      lit.AddSignal(trk,ht,mod,index,signal);
    }
    if(index>=fEntries) fEntries = index +1;
}
//______________________________________________________________________
void AliITSpList::AddNoise(Int_t i,Int_t j,Int_t mod,Double_t noise){
    // Adds a noise value to the TClonesArray at i,j. 
    // Creates the AliITSpListItem
    // if needed.
    // Inputs:
    //    Int_t i        Row number for this noise
    //    Int_t j        Column number for this noise
    //    Double_t noise The noise signal value.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t index = GetIndex(i,j);
    if (index<0) return;
    AliITSpListItem &lit = fa[index];
    if(!lit.IsUsed()){
      lit.Build(mod,index,noise);
    }
    else {
      lit.AddNoise(mod,index,noise);
    } 
    if(index>=fEntries) fEntries = index +1;
}
//______________________________________________________________________
void AliITSpList::GetCell(Int_t index,Int_t &i,Int_t &j) const {
  // returns the i,j index numbers from the linearized index computed
  // with GetIndex
  if(index<0 || index>=fNi*fNj){
    Warning("GetCell","Index out of range 0<=index=%d<%d",
	    index,fNi*fNj);
    i=-1;j=-1;
    return;
  } // end if
  i = index/fNj;
  j = index - fNj*i;
  return;
}
