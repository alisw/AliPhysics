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

////////////////////////////////////////////////////////////////////////
//
// AliTrackMap.cxx
// description: 
//   contains a relation between track label and it's index
//   in a TreeH.
//   The main method is At, it takes a particle label as an argument
//   and returns the correponding entry in TreeH (many particles
//   are usually stored in 1 TreeH entry, one has to check particle
//   labels for each hit). 'At' returns:
//       kNoEntry = -1 if the particle has no entry in TreeH
//       kOutOfBounds = -2 if the particle label is out of bounds
//                  
//  Author: Jiri Chudoba (CERN), 2002
//
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliLog.h"
#include "AliTrackMap.h"

ClassImp(AliTrackMap)

//_______________________________________________________________________
AliTrackMap::AliTrackMap():
  fSize(0),
  fArray(0)
{
  //
  // default ctor
  //
}

//_______________________________________________________________________
AliTrackMap::AliTrackMap(const AliTrackMap& trm):
  TNamed(trm),
  fSize(0),
  fArray(0)
{
  //
  // default ctor
  //
  trm.Copy(*this);
}

//_______________________________________________________________________
AliTrackMap::AliTrackMap(Int_t size, Int_t *array): 
  TNamed("AliTrackMap", "AliTrackMap"),
  fSize(size),
  fArray(new Int_t[fSize])
{
  //
  // ctor
  //
  for (Int_t i = 0; i < fSize; i++) fArray[i] = array[i];
}

//_______________________________________________________________________
void AliTrackMap::Copy(TObject& ) const
{
  AliFatal("Not implemented");
}

//_______________________________________________________________________
AliTrackMap::~AliTrackMap()
{
  //
  // dtor
  //
  delete [] fArray;
}

//_______________________________________________________________________
Int_t AliTrackMap::At(Int_t label) const
{
  //
  // returns entry number in the TreeH corresponding to particle with 
  // label label
  //
  if (label < 0 || label >= fSize) {
    AliError(Form("label %d out of range, fSize = %d", label, fSize));
    return kOutOfBounds;
  }
  return fArray[label];
}

//_______________________________________________________________________
void AliTrackMap::SetEventNr(Int_t eventNr) 
{
  //
  // map is identified by it's name, event number is part of it
  //
  char name[20];
  sprintf(name,"AliTrackMap_%5.5d",eventNr);
  SetName(name);  
}

//_______________________________________________________________________
void AliTrackMap::PrintValues() const
{
  //
  // method for debugging
  //
  cout<<this->GetName()<<" contains these values: "<<endl;
  for (Int_t i = 0; i < fSize; i++) {
    cout<<i<<"\t"<<fArray[i]<<"\n";
  }
}

