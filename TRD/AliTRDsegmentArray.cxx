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

/*
$Log$
Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:55:03  cblume
Bug fixes

Revision 1.1  2000/02/28 19:02:32  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Alice segment manager class                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRD.h"
#include "AliTRDgeometry.h"
#include "AliTRDsegmentArray.h"

ClassImp(AliTRDsegmentArray)

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray():AliTRDsegmentArrayBase()
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray(Text_t *classname, Int_t n)
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
}

//_____________________________________________________________________________
void AliTRDsegmentArray::Copy(TObject &a)
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
Bool_t AliTRDsegmentArray::LoadArray(const Char_t *branchname)
{
  //
  // Loads all segments of the array from the branch <branchname> of
  // the digits tree
  //

  // Connect the digits tree
  fTree = gAlice->TreeD();
  if (!fTree) return kFALSE;

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) return kFALSE;

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
Bool_t AliTRDsegmentArray::StoreArray(const Char_t *branchname)
{
  //
  // Stores all segments of the array in the branch <branchname> of 
  // the digits tree
  //

  // Connect the digits tree
  fTree = gAlice->TreeD();
  if (!fTree) return kFALSE;

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) return kFALSE;

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
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t det)
{
  //
  // Returns the data array for a given detector
  //

  return ((AliTRDdataArray *) AliTRDsegmentArrayBase::At(det));

}

//_____________________________________________________________________________
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t pla, Int_t cha, Int_t sec)
{
  //
  // Returns the data array for a given detector
  //

  if (gAlice) {

    AliTRDgeometry *geo = ((AliTRD*) gAlice->GetDetector("TRD"))->GetGeometry();  
    Int_t det = geo->GetDetector(pla,cha,sec);
    return GetDataArray(det);

  }
  else {

    printf("AliTRDsegmentArray::GetDigits -- ");
    printf("gAlice is not defined\n");
    return NULL;

  }

}
