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

//-----------------------------------------------------------------------------
// Class AliMUONAlignmentClusterRecord
// Class AliMUONAlignmentTrackRecord
//-----------------------------------------------------------------------------

#include "AliMUONAlignmentRecord.h"
#include "AliLog.h"

#include <iostream>

// class implementation in root context
ClassImp( AliMUONAlignmentClusterRecord )

// class implementation in root context
ClassImp( AliMUONAlignmentTrackRecord )

//__________________________________________________________________________
AliMUONAlignmentClusterRecord::AliMUONAlignmentClusterRecord( void ):
  fDetElemId(0),
  fDetElemNumber(0),
  fMeas(0),
  fSigma(0)
{

  // initialize derivatives
  for( Int_t index = 0; index < 4; ++index )
  {
    fLocalDerivatives[index] = 0;
    fGlobalDerivatives[index] = 0;
  }

}

//__________________________________________________________________________
void AliMUONAlignmentClusterRecord::Print( Option_t* ) const
{

  // detector id and measurement
  std::cout
    << " AliMUONAlignmentClusterRecord::Print - fDetElemId: " << fDetElemId
    << " fMeas: " << fMeas
    << " fSigma: " << fSigma;

  // local derivatives
  std::cout << " fLocalDerivatives: ";
  for( Int_t index = 0; index < 4; ++index )
  { std::cout << fLocalDerivatives[index] << " "; }

  // global derivatives
  std::cout << " fGlobalDerivatives: ";
  for( Int_t index = 0; index < 4; ++index )
  { std::cout << fGlobalDerivatives[index] << " "; }

  std::cout << std::endl;

}

//__________________________________________________________________________
AliMUONAlignmentTrackRecord::AliMUONAlignmentTrackRecord( void ):
  fClusterRecords (new TClonesArray( "AliMUONAlignmentClusterRecord", fSize ) ),
  fClusterCount( 0 )
  {}

//__________________________________________________________________________
AliMUONAlignmentTrackRecord::AliMUONAlignmentTrackRecord( const AliMUONAlignmentTrackRecord& other ):
  TObject(other),
  fClusterRecords (new TClonesArray( "AliMUONAlignmentClusterRecord", fSize ) ),
  fClusterCount( other.fClusterCount )
{

  // deep copy of records
  for( Int_t index = 0; index < other.GetNRecords(); ++index )
  { new ((*fClusterRecords)[index]) AliMUONAlignmentClusterRecord( *static_cast<const AliMUONAlignmentClusterRecord*>( other.fClusterRecords->UncheckedAt( index ) ) ); }

}

//__________________________________________________________________________
AliMUONAlignmentTrackRecord& AliMUONAlignmentTrackRecord::operator = ( const AliMUONAlignmentTrackRecord& other )
{

  if( this == &other ) return *this;
  Clear();

  TObject::operator =( other );

  // copy number of cluster
  fClusterCount = other.fClusterCount;

  // deep copy of records
  for( Int_t index = 0; index < other.GetNRecords(); ++index )
  { new ((*fClusterRecords)[index]) AliMUONAlignmentClusterRecord( *static_cast<const AliMUONAlignmentClusterRecord*>( other.fClusterRecords->UncheckedAt( index ) ) ); }

  return *this;

}


//__________________________________________________________________________
void AliMUONAlignmentTrackRecord::Print( Option_t* ) const
{

  std::cout << "AliMUONAlignmentTrackRecord::Print - fClusterCount: " << fClusterCount << std::endl;
  if( fClusterRecords )
  {
    for( Int_t index = 0; index < fClusterCount; ++index )
    { static_cast<AliMUONAlignmentClusterRecord*>( fClusterRecords->UncheckedAt( index ) )->Print(); }
  }

}

//__________________________________________________________________________
AliMUONAlignmentTrackRecord::~AliMUONAlignmentTrackRecord( void )
{

  // note: it is safe to delete a NULL pointer, so no need to check.
  delete fClusterRecords;

}

//__________________________________________________________________________
void AliMUONAlignmentTrackRecord::AddClusterRecord( const AliMUONAlignmentClusterRecord& record )
{

  // append copy to array
  new ((*fClusterRecords)[fClusterCount]) AliMUONAlignmentClusterRecord( record );
  fClusterCount++;

}

//__________________________________________________________________________
void AliMUONAlignmentTrackRecord::RemoveClusterRecord( AliMUONAlignmentClusterRecord* record )
{

  AliMUONAlignmentClusterRecord* local( static_cast<AliMUONAlignmentClusterRecord*>( fClusterRecords->Remove( record ) ) );

  if( local )
  {

    delete local;
    fClusterRecords->Compress();
    fClusterCount--;

  } else AliWarning("object to remove does not exist in array fClusterRecords");

}

//__________________________________________________________________________
void AliMUONAlignmentTrackRecord::Clear( Option_t* options )
{

  TObject::Clear( options );

  // only Clear, not delete, since fClusterRecords does not allocate memory
  fClusterRecords->Clear();

  // reset count
  fClusterCount = 0;

}
