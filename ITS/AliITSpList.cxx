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
// It consist of a TObjectArray of 
// AliITSpListItem objects
// This array can be accessed via 2 indexed
// it is used at digitization level by 
// all the 3 ITS subdetectors
//
// ***********************************************************************

#include <TObjArray.h>

#include "AliITSpList.h"
#include "AliITSpListItem.h"


//______________________________________________________________________

ClassImp(AliITSpList)
//______________________________________________________________________
AliITSpList::AliITSpList(){
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A zeroed/empty AliITSpList class.

    fNi = 0;
    fNj = 0;
    fa  = 0;
}
//______________________________________________________________________
AliITSpList::AliITSpList(Int_t imax,Int_t jmax){
    // Standard constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A setup AliITSpList class.

    fNi = imax;
    fNj = jmax;
    fEnteries = 0;
    fa  = new TObjArray(fNi*fNj); // elements are zeroed by 
                                  // TObjArray creator
}
//______________________________________________________________________
AliITSpList::~AliITSpList(){
    // Default destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    a properly destroyed class

    for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
        delete fa->At(i);
        fa->AddAt(0,i); // zero content
    } // end for i && if
    fNi = 0;
    fNj = 0;
    delete fa;
    fa  = 0;
    fEnteries = 0;
}

//______________________________________________________________________
void AliITSpList::ClearMap(){
    // Delete all AliITSpListItems and zero TObjArray.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A zeroed AliITSpList class.

    fa->Delete();
    /*
    for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
        delete fa->At(i);
        fa->AddAt(0,i); // zero content
    } // end for i && if
    */
    fEnteries = 0;
}
//______________________________________________________________________
void AliITSpList::DeleteHit(Int_t i,Int_t j){
    // Delete a particular AliITSpListItems and zero TObjArray.
    // Inputs:
    //    Int_t i   Row number
    //    Int_t j   Columns number
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t k = GetIndex(i,j);

    if(fa->At(k)!=0){
        delete fa->At(k);
        fa->AddAt(0,k); // zero content
    } // end for i && if
    if(k==fEnteries-1) fEnteries--;
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

    if(this == &source) return *this;

    if(this->fa!=0){ // if this->fa exists delete it first.
        for(Int_t i=0;i<GetMaxIndex();i++) if(fa->At(i)!=0){
            delete fa->At(i);
            fa->AddAt(0,i); // zero content
        } // end for i && if
        delete this->fa;
    } // end if this->fa!=0
    this->fNi = source.fNi;
    this->fNj = source.fNj;
    this->fa = new TObjArray(*(source.fa));
    this->fEnteries = source.fEnteries;

    return *this;
}
//______________________________________________________________________
AliITSpList::AliITSpList(const AliITSpList &source) : AliITSMap(source){
    // Copy constructor

  fNi = source.fNi;
  fNj = source.fNj;
  fa = new TObjArray(*(source.fa));
  fEnteries = source.fEnteries;
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

    if( fa->At( index ) == 0 ) { // most create AliITSpListItem
        fa->AddAt(new AliITSpListItem(-2,-1,pl->GetModule(),index,0.0),index);
    } // end if
 
    ((AliITSpListItem*)(fa->At(index)))->AddTo( fileIndex,pl);
    if(index>=fEnteries) fEnteries = index +1;
}
//______________________________________________________________________
void AliITSpList::AddSignal(Int_t i,Int_t j,Int_t trk,Int_t ht,Int_t mod,
                       Double_t signal){
    // Adds a Signal value to the TObjArray at i,j. Creates the AliITSpListItem
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

    if(GetpListItem(index)==0){ // most create AliITSpListItem
        fa->AddAt(new AliITSpListItem(trk,ht,mod,index,signal),index);
    }else{ // AliITSpListItem exists, just add signal to it.
        GetpListItem(index)->AddSignal(trk,ht,mod,index,signal);
    } // end if
    if(index>=fEnteries) fEnteries = index +1;
}
//______________________________________________________________________
void AliITSpList::AddNoise(Int_t i,Int_t j,Int_t mod,Double_t noise){
    // Adds a noise value to the TObjArray at i,j. Creates the AliITSpListItem
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

    if(GetpListItem(index)==0){ // most create AliITSpListItem
        fa->AddAt(new AliITSpListItem(mod,index,noise),index);
    }else{ // AliITSpListItem exists, just add signal to it.
        GetpListItem(index)->AddNoise(mod,index,noise);
    } // end if
    if(index>=fEnteries) fEnteries = index +1;
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
//______________________________________________________________________
Int_t AliITSpList::GetIndex(Int_t i, Int_t j) const {
 // returns the TObjArray index for a given set of map indexes.
  if(i<0||i>=fNi || j<0||j>=fNj){
    Warning("GetIndex","Index out of range 0<i=%d<%d and 0<0j=%d<%d",i,fNi,j,fNj);
    return -1;
  }
  else {
    return fNj*i+j;
  }
}
