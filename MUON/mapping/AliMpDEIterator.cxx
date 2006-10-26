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

// $Id$ 
// $MpId: AliMpDEIterator.cxx,v 1.6 2006/05/24 13:58:34 ivana Exp $
// Category: management

// ------------------------
// Class AliMpDEIterator
// ------------------------
// The iterator over valid detection element IDs
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpFiles.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMpDEIterator)
/// \endcond

const  Int_t  AliMpDEIterator::fgkMaxNofDetElements = 250;
TArrayI  AliMpDEIterator::fgDetElemIds(fgkMaxNofDetElements);
Int_t    AliMpDEIterator::fgNofDetElemIds = 0;	

//
// static private methods
//

//______________________________________________________________________________
Bool_t AliMpDEIterator::ReadDEIds(AliMpStationType station)
{ 
/// Read det element ids from the file specified by name
/// and fill the map (the deNames are ignored)
/// Return true if the data were read ok

  // Open file
  TString filePath = AliMpFiles::DENamesFilePath(station);
  std::ifstream in(filePath);
  if (!in.good()) {
    AliErrorClassStream() << "Cannot open file " << filePath << endl;;
    return false;
  }
  
  // Skip plane types per cathods + empty lines
  //
  char line[80];
  in.getline(line, 80);
  in.getline(line, 80);
  in.getline(line, 80);
    
  // Read DE Ids
  //
  Int_t detElemId;
  TString word;
  in >> word;
  while ( ! in.eof() ) {
    if ( word[0] == '#' ) {
      in.getline(line, 80);
    }  
    else {  
      detElemId = word.Atoi();
      in.getline(line, 80);
      AliDebugClassStream(3) 
        << "Adding  " << fgNofDetElemIds << "  "  << detElemId << endl;
      fgDetElemIds.AddAt(detElemId, fgNofDetElemIds++);
    } 
    in >> word;
  }

  // Close file
  in.close();
  
  return true;
}

//______________________________________________________________________________
void AliMpDEIterator::ReadData()
{
/// Fill DE Ids array from DE names files
/// Return true if all data were read ok

  AliDebugClass(2,"");

  Bool_t result1 = ReadDEIds(kStation1);
  Bool_t result2 = ReadDEIds(kStation2);
  Bool_t result3 = ReadDEIds(kStation345);
  Bool_t result4 = ReadDEIds(kStationTrigger);
  
  Bool_t result = result1 && result2 && result3 && result4;
  if ( ! result ) {
    AliErrorClassStream() << "Error in reading DE names files" << endl;
  }  
}

//
// constructors, destructor
//

//______________________________________________________________________________
AliMpDEIterator::AliMpDEIterator()
    : TObject(),
      fIndex(-1),
      fChamberId(-1)
{  
/// Standard and default constructor

  if (! fgNofDetElemIds ) ReadData();
}

//______________________________________________________________________________
AliMpDEIterator::AliMpDEIterator(const AliMpDEIterator& rhs)
 : TObject(rhs),
   fIndex(rhs.fIndex),
   fChamberId(rhs.fChamberId)
{
/// Copy constructor
}

//______________________________________________________________________________

AliMpDEIterator::~AliMpDEIterator()
{
/// Destructor
}

//______________________________________________________________________________
AliMpDEIterator&  AliMpDEIterator::operator=(const AliMpDEIterator& rhs)
{
/// Assignement operator

  // check assignment to self
  if (this == &rhs) return *this;

  // base class assignment
  TObject::operator=(rhs);

  fIndex = rhs.fIndex;
  fChamberId = rhs.fChamberId;

  return *this;
} 

//
// public methods
//

//______________________________________________________________________________
void AliMpDEIterator::First()
{
/// Set iterator to the first DE Id defined 

  fIndex = 0;
  fChamberId = -1;
}  

//______________________________________________________________________________
void AliMpDEIterator::First(Int_t chamberId)
{
/// Reset the iterator, so that it points to the first DE
 
  fChamberId = -1;
  fIndex = -1;  
  if ( ! AliMpDEManager::IsValidChamberId(chamberId) ) {
    AliErrorStream() << "Invalid chamber Id " << chamberId << endl;
    return;
  }    

  Int_t i=0;
  while ( i < fgNofDetElemIds && fChamberId < 0 ) {
    Int_t detElemId = fgDetElemIds.At(i);
    if ( AliMpDEManager::GetChamberId(detElemId) == chamberId ) {
      fChamberId = chamberId;
      fIndex = i;
    } 
    i++; 
  }

  if ( fChamberId < 0 ) {
    AliErrorStream() 
      << "No DEs of Chamber Id " << chamberId << " found" << endl;
    return;
  }    

}

//______________________________________________________________________________
void AliMpDEIterator::Next()
{
/// Increment iterator to next DE

  fIndex++;

  // Invalidate if at the end
  if ( ( fIndex == fgNofDetElemIds ) ||
       ( fChamberId >= 0 &&    
         AliMpDEManager::GetChamberId(CurrentDE()) != fChamberId ) ) {
    fIndex = -1;
  }   
}

//______________________________________________________________________________
Bool_t AliMpDEIterator::IsDone() const
{
/// Is the iterator in the end?

  return ( fIndex < 0 );
}   

//______________________________________________________________________________
Int_t AliMpDEIterator::CurrentDE() const
{
/// Current DE Id

  if ( ! IsDone() )
    return fgDetElemIds.At(fIndex);
  else {   
    AliErrorStream()
      << "Not in valid position - returning invalid DE." << endl;
    return 0;
  }  
}
    
