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
// $MpId: AliMpSector.cxx,v 1.14 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpSector
// -----------------
// Class describing the sector of the MUON chamber of station 1.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpSector.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"
#include "AliMpMotifMap.h"
#include "AliMpEncodePair.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpSector)
/// \endcond

//_____________________________________________________________________________
AliMpSector::AliMpSector(const TString& id, Int_t nofZones, Int_t nofRows, 
                         AliMp::Direction direction, 
                         Double_t offsetx, Double_t offsety) 
  : TNamed("Sector", ""),
    fID(id),
    fOffsetX(offsetx),
    fOffsetY(offsety),
    fDimensionX(0.),
    fDimensionY(0.),
    fZones(),
    fRows(),
    fMotifMap(0),
    fDirection(direction),
    fMinPadDimensionX(1.e6),
    fMinPadDimensionY(1.e6),
    fMaxPadDimensionX(0.),
    fMaxPadDimensionY(0.),
    fLMaxPadIndices(0),
    fNofPads(0)
{
/// Standard constructor

  AliDebugStream(1) << "this = " << this << endl;

  fMotifMap = new AliMpMotifMap;

  for (Int_t izone = 0; izone<nofZones; izone++) 
    fZones.Add(new AliMpZone(izone+1));
    
  for (Int_t irow = 0; irow<nofRows; irow++) 
    fRows.Add(new AliMpRow(irow, fMotifMap));
    
}

//_____________________________________________________________________________
AliMpSector::AliMpSector() 
  : TNamed(),
    fID(""),    
    fOffsetX(0.),
    fOffsetY(0.),
    fDimensionX(0.),
    fDimensionY(0.),
    fZones(),
    fRows(),
    fMotifMap(0),
    fDirection(AliMp::kX),
    fMinPadDimensionX(0.),
    fMinPadDimensionY(0.),
    fMaxPadDimensionX(0.),
    fMaxPadDimensionY(0.),
    fLMaxPadIndices(0),
    fNofPads(0)
{
/// Default constructor

  AliDebugStream(1) << "this = " << this << endl;
}

//_____________________________________________________________________________
AliMpSector::~AliMpSector() 
{
/// Destructor 

  AliDebugStream(1) << "this = " << this << endl;

  // deletes 
  for (Int_t izone = 0; izone<GetNofZones(); izone++) 
    delete fZones[izone];
    
  for (Int_t irow = 0; irow<GetNofRows(); irow++) 
    delete fRows[irow];

  delete fMotifMap;
}

//
// private methods
//

//_____________________________________________________________________________
AliMpRow* AliMpSector::FindRow(Double_t y) const
{
/// Find the row for the specified y position.                              \n
/// If y is on border the lowest row is returned.
  
  for (Int_t i=0; i<GetNofRows(); i++) {
    if ( y >= ((AliMpRow*)fRows[i])->LowBorderY() && 
         y <= ((AliMpRow*)fRows[i])->UpperBorderY())
      return (AliMpRow*)fRows[i];
  }    
  
  return 0;
}

//_____________________________________________________________________________
Int_t AliMpSector::FindMotifPositionId(Double_t x, Double_t y) const
{
/// Find the motif position ID in the specified position.                   \n
/// Return 0 if no motif is found.
 
  // Find the row segment
  AliMpVRowSegment* rowSegment = FindRowSegment(x,y);
  
  if ( ! rowSegment ) return 0;
    
  // Find motif position ID
  return rowSegment->FindMotifPositionId(x, y);  
}

//_____________________________________________________________________________
AliMpVRowSegment* AliMpSector::FindRowSegment(Double_t x, Double_t y) const
{
/// Find the row segment in the specified position.                         \n
/// Return if no motif is found.
  
  // Find row
  AliMpRow* row = FindRow(y);
  
  if (!row) return 0;

  // Find the row segment and return its motif
  AliMpVRowSegment* rowSegment = row->FindRowSegment(x);
  
  return rowSegment;
}

//_____________________________________________________________________________
void  AliMpSector::SetRowOffsets()
{
/// For each row check consitency of the row segments
/// and calculate the row offset.

  Double_t offset = fOffsetY;
  
  for (Int_t irow=0; irow<GetNofRows(); irow++)
    offset = GetRow(irow)->SetOffsetY(offset);    
}

//_____________________________________________________________________________
void  AliMpSector::SetMotifPositions()
{
/// Create motif positions objects and fills them in the motif map.

  for (Int_t i=0; i<GetNofRows(); i++)
    GetRow(i)->SetMotifPositions();
}

//_____________________________________________________________________________
void  AliMpSector::SetGlobalIndices()
{
/// Set the indices limits to all indexed elements
/// (row, row segment, motif positions).

  AliMpRow* rowBefore=0;
  for (Int_t i=0; i<GetNofRows(); i++) {
    GetRow(i)->SetGlobalIndices(fDirection, rowBefore);
    rowBefore = GetRow(i);
  }
}

//_____________________________________________________________________________
void  AliMpSector::SetMinMaxPadDimensions()
{
/// Set the minimal pad dimensions.

  for (Int_t i=1; i<GetNofZones()+1; i++) {
    Double_t dx = GetZone(i)->GetPadDimensionX();
    Double_t dy = GetZone(i)->GetPadDimensionY();
    
    if ( ( fDirection == AliMp::kX && dy > 0. && dy < fMinPadDimensionY ) ||
         ( fDirection == AliMp::kY && dx > 0. && dx < fMinPadDimensionX ) ) {
      
      fMinPadDimensionX = dx;
      fMinPadDimensionY = dy;
    }  

    if ( ( fDirection == AliMp::kX && dy > 0. && dy > fMaxPadDimensionY ) ||
         ( fDirection == AliMp::kY && dx > 0. && dx > fMinPadDimensionX ) ) {
      
      fMaxPadDimensionX = dx;
      fMaxPadDimensionY = dy;
    }  
  }
}

//_____________________________________________________________________________
void  AliMpSector::SetMaxPadIndices()
{
/// Set maximum pad indices in x, y

  if ( fLMaxPadIndices != 0 ) return;
  
  Int_t maxIndexInX = 0;
  Int_t maxIndexInY = 0;
  for (Int_t i=0; i<GetNofRows(); i++) {

    Int_t ixh = GetRow(i)->GetHighLimitIx();
    if ( ixh > maxIndexInX ) maxIndexInX = ixh;

    Int_t iyh = GetRow(i)->GetHighLimitIy();
    if ( iyh > maxIndexInY ) maxIndexInY = iyh;
  }  
  
  fLMaxPadIndices = AliMp::Pair(maxIndexInX, maxIndexInY);
}


//_____________________________________________________________________________
void  AliMpSector::SetNofPads()
{
/// Set the total number of pads

  fNofPads = fMotifMap->CalculateNofPads();
}

//_____________________________________________________________________________
void  AliMpSector::SetDimensions()
{
/// Set the maximum halflengths in x, y.

  fDimensionX = 0.;
  fDimensionY = 0.;
  
  for (Int_t i=0; i<GetNofRows(); i++) {

    // take the largest x row dimension
    if ( ((AliMpRow*)fRows[i])->GetDimensionX() > fDimensionX ) 
      fDimensionX = ((AliMpRow*)fRows[i])->GetDimensionX();
      
    // add all rows y dimensions  
    fDimensionY += ((AliMpRow*)fRows[i])->GetDimensionY();
  }
}  

//
// public methods
//

//_____________________________________________________________________________
AliMpVPadIterator* AliMpSector::CreateIterator() const
{
/// Create sector pad iterator

  return new AliMpSectorPadIterator(this);
}


//_____________________________________________________________________________
void  AliMpSector::SetRowSegmentOffsets()
{
/// For all rows set the offset to all row segments.

  for (Int_t irow=0; irow<GetNofRows(); irow++)
    GetRow(irow)->SetRowSegmentOffsets(fOffsetX);    
}

//_____________________________________________________________________________
void AliMpSector::Initialize() 
{
/// Make needed settings after sector is read from
/// data files.

  SetRowOffsets();
  SetMotifPositions();
  SetGlobalIndices();
  SetMinMaxPadDimensions();
  SetMaxPadIndices();
  SetNofPads();
  SetDimensions();
}  

//_____________________________________________________________________________
void AliMpSector::PrintGeometry()  const
{
/// Print the positions of rows, rows segments

  for (Int_t i=0; i<GetNofRows(); i++) {
    AliMpRow* row = GetRow(i);
    
    cout << "ROW " << row->GetID() 
         << "  center Y " << row->GetPositionY() << endl;

    for (Int_t j=0; j<row->GetNofRowSegments(); j++) {
       AliMpVRowSegment* rowSegment = row->GetRowSegment(j);
	
       cout << "   ROW Segment " << j 
            << "  borders " 
            << rowSegment->LeftBorderX() << "  "
            << rowSegment->RightBorderX()
            << "  x-size " 
            << 2*rowSegment->GetDimensionX() << "  "
	    << endl;
    }
  }
}     
      	     

//_____________________________________________________________________________
AliMpRow* AliMpSector::FindRow(Int_t motifPositionId) const
{
/// Find the row with the the specified motif position.                     \n
/// Return 0 if no row is found.

  AliMpVRowSegment* segment = FindRowSegment(motifPositionId);
  
  if (segment) return segment->GetRow();
  
  return 0;  
}

//_____________________________________________________________________________
AliMpVRowSegment* AliMpSector::FindRowSegment(Int_t motifPositionId) const
{
/// Find the row segment with the the specified motif position.            \n
/// Return 0 if no row segment is found.

  for (Int_t irow=0; irow<GetNofRows(); irow++) {

    AliMpRow* row = (AliMpRow*)fRows[irow];

    for (Int_t iseg=0; iseg<row->GetNofRowSegments(); iseg++) {
      AliMpVRowSegment* segment = row->GetRowSegment(iseg); 
      if (segment->HasMotifPosition(motifPositionId)) return segment;
    }
  }
  
  return 0;    
}

//_____________________________________________________________________________
Double_t AliMpSector::GetPositionX() const
{
/// Return the sector offset.

  return fOffsetX;
}  


//_____________________________________________________________________________
Double_t AliMpSector::GetPositionY() const
{
/// Return the sector offset.

  return fOffsetY;
}  

//_____________________________________________________________________________
Double_t  AliMpSector::GetDimensionX() const
{
/// Return the maximum halflengths in x.

  return fDimensionX;
}  

//_____________________________________________________________________________
Double_t  AliMpSector::GetDimensionY() const
{
/// Return the maximum halflengths in y.

  return fDimensionY;
}  

//_____________________________________________________________________________
Int_t AliMpSector::GetNofZones() const
{    
/// Return the number of zones.

  return fZones.GetEntriesFast();
}  

//_____________________________________________________________________________
AliMpZone* AliMpSector::GetZone(Int_t zoneID) const
{
/// Return zone with specified ID.

  if (zoneID < 1 || zoneID > GetNofZones()) {
    AliWarningStream() << "Index outside range" << endl;
    return 0;
  }
  
  return (AliMpZone*)fZones[zoneID-1];
}  

//_____________________________________________________________________________
Int_t AliMpSector::GetNofRows() const
{
/// Return the number of rows.

  return fRows.GetEntriesFast();
}  

//_____________________________________________________________________________
AliMpRow* AliMpSector::GetRow(Int_t rowID) const
{
/// Return row with specified ID.

  if (rowID < 0 || rowID >= GetNofRows()) {
    AliWarningStream() << "Index outside range" << endl;
    return 0;
  }
  
  return (AliMpRow*)fRows[rowID];
}

//_____________________________________________________________________________
AliMp::PlaneType
AliMpSector::GetPlaneType() const
{
/// Return the plane type

  return GetDirection()==AliMp::kY ? AliMp::kBendingPlane : AliMp::kNonBendingPlane;
}

//_____________________________________________________________________________
Int_t 
AliMpSector::GetNofMotifPositions() const
{
  /// Return the number of manus
  
  return fMotifMap->GetNofMotifPositions();
}

//_____________________________________________________________________________
void 
AliMpSector::GetAllMotifPositionsIDs(TArrayI& ecn) const
{
/// Return the array of all motif positions IDs

  fMotifMap->GetAllMotifPositionsIDs(ecn);
}

//_____________________________________________________________________________
void
AliMpSector::Print(Option_t* opt) const
{
/// Print the map of motifs

  cout << "Sector," << PlaneTypeName(GetPlaneType()) << endl;
  fMotifMap->Print(opt);
}
