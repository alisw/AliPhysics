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
//  Alice segment manager class                                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TTree.h>

#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdataArray.h"

ClassImp(AliTRDsegmentArray)

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray()
  :AliTRDsegmentArrayBase()
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray(const char *classname, Int_t n)
  :AliTRDsegmentArrayBase(classname,n)
{
  //
  // Constructor creating an array of AliTRDdataArray of size <n>
  //

  AliTRDdataArray *dataArray;  

  for (Int_t i = 0; i < n; i++) {
    dataArray = (AliTRDdataArray *) AddSegment(i);
  }

}

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray(AliTRDsegmentArray &a)
  :AliTRDsegmentArrayBase(a)
{
  //
  // AliTRDsegmentArray copy constructor
  //

  a.Copy(*this);

}

//_____________________________________________________________________________
AliTRDsegmentArray::~AliTRDsegmentArray()
{
  //
  // AliTRDsegmentArray destructor
  //

  Delete();

}

//_____________________________________________________________________________
void AliTRDsegmentArray::Copy(TObject &a) const
{
  //
  // Copy function
  //

  AliTRDsegmentArrayBase::Copy(a);

}

//_____________________________________________________________________________
void AliTRDsegmentArray::Delete()
{
  //
  // Deletes all detector segments from the array
  //

  for (Int_t iDet = 0; iDet < fNSegment; iDet++) {
    ClearSegment(iDet);
  }

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::LoadArray(const Char_t *branchname, TTree *tree)
{
  //
  // Loads all segments of the array from the branch <branchname> of
  // the digits tree <tree>
  //

  fTree = tree;

  if (!fTree) {
    AliError("Digits tree is not defined\n");
    return kFALSE;
  }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) {
    AliError(Form("Branch %s is not defined\n",branchname));
    return kFALSE;
  }

  // Loop through all segments and read them from the tree
  Bool_t status = kTRUE;
  for (Int_t iSegment = 0; iSegment < fNSegment; iSegment++) {
    AliTRDdataArray *dataArray = (AliTRDdataArray *) fSegment->At(iSegment);
    if (!dataArray) {
      status = kFALSE;
      break;    
    }
    fBranch->SetAddress(&dataArray);
    fBranch->GetEntry(iSegment);
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::StoreArray(const Char_t *branchname, TTree *tree)
{
  //
  // Stores all segments of the array in the branch <branchname> of 
  // the digits tree <tree>
  //

  fTree = tree;

  if (!fTree) {
    AliError("Digits tree is not defined\n");
    return kFALSE;
  }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) {
    AliError(Form("Branch %s is not defined\n",branchname));
    return kFALSE;
  }

  // Loop through all segments and fill them into the tree
  Bool_t status = kTRUE;
  for (Int_t iSegment = 0; iSegment < fNSegment; iSegment++) {
    const AliTRDdataArray *kDataArray = 
         (AliTRDdataArray *) AliTRDsegmentArrayBase::At(iSegment);
    if (!kDataArray) {
      status = kFALSE;
      break;
    }
    fBranch->SetAddress(&kDataArray);
    fBranch->Fill();
  }

  return status;

}

//_____________________________________________________________________________
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t det) const
{
  //
  // Returns the data array for a given detector
  //

  return ((AliTRDdataArray *) AliTRDsegmentArrayBase::At(det));

}

//_____________________________________________________________________________
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t pla
                                                , Int_t cha
                                                , Int_t sec) const
{
  //
  // Returns the data array for a given detector
  //

  Int_t det = AliTRDgeometry::GetDetector(pla,cha,sec);
  return GetDataArray(det);

}
