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
//
// Class AliMpSector
// -----------------
// Class describing the sector of the MUON chamber of station 1.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpSector.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"
#include "AliMpMotifMap.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpSector)
/// \endcond

//_____________________________________________________________________________
AliMpSector::AliMpSector(const TString& id, Int_t nofZones, Int_t nofRows, 
                         AliMp::Direction direction, const TVector2& offset) 
  : TNamed("Sector", ""),
    fID(id),
    fOffset(offset),
    fZones(),
    fRows(),
    fMotifMap(0),
    fDirection(direction),
    fMinPadDimensions(TVector2(1.e6, 1.e6)),
    fMaxPadIndices(AliMpIntPair::Invalid()),
    fNofPads(0)
{
/// Standard constructor

  AliDebugStream(1) << "this = " << this << endl;

  fMotifMap = new AliMpMotifMap(true);
  //fMotifMap = new AliMpMotifMap();

#ifdef WITH_STL
  for (Int_t izone = 0; izone<nofZones; izone++) 
    fZones.push_back(new AliMpZone(izone+1));
    
  for (Int_t irow = 0; irow<nofRows; irow++) 
    fRows.push_back(new AliMpRow(irow, fMotifMap));
#endif

#ifdef WITH_ROOT
  for (Int_t izone = 0; izone<nofZones; izone++) 
    fZones.Add(new AliMpZone(izone+1));
    
  for (Int_t irow = 0; irow<nofRows; irow++) 
    fRows.Add(new AliMpRow(irow, fMotifMap));
#endif
}

//_____________________________________________________________________________
AliMpSector::AliMpSector() 
  : TNamed(),
    fID(""),    
    fOffset(TVector2(0., 0.)),
    fZones(),
    fRows(),
    fMotifMap(0),
    fDirection(AliMp::kX),
    fMinPadDimensions(TVector2(0., 0.)),
    fMaxPadIndices(AliMpIntPair::Invalid()),
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
AliMpVPadIterator* AliMpSector::CreateIterator() const
{
/// Create sector pad iterator

  return new AliMpSectorPadIterator(this);
}


//_____________________________________________________________________________
AliMpVRowSegment* AliMpSector::FindRowSegment(const TVector2& position) const
{
/// Find the row segment in the specified position.                         \n
/// Return if no motif is found.
  
  // Find row
  AliMpRow* row = FindRow(position);
  
  if (!row) return 0;

  // Find the row segment and return its motif
  AliMpVRowSegment* rowSegment = row->FindRowSegment(position.X());
  
  return rowSegment;
}

//_____________________________________________________________________________
void  AliMpSector::SetRowOffsets()
{
/// For each row check consitency of the row segments
/// and calculate the row offset.

  Double_t offset = fOffset.Y();
  
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

  AliMpIntPair indices(0,0); 
  AliMpRow* rowBefore=0;
  for (Int_t i=0; i<GetNofRows(); i++) {
    GetRow(i)->SetGlobalIndices(fDirection, rowBefore);
    rowBefore = GetRow(i);
  }
}

//_____________________________________________________________________________
void  AliMpSector::SetMinPadDimensions()
{
/// Set the minimal pad dimensions.

  for (Int_t i=1; i<GetNofZones()+1; i++) {
    TVector2 padDimensions = GetZone(i)->GetPadDimensions();
    
    if ( fDirection == AliMp::kX &&  
         padDimensions.Y() > 0. && padDimensions.Y() < fMinPadDimensions.Y() ||
         fDirection == AliMp::kY && 
	 padDimensions.X() > 0. && padDimensions.X() < fMinPadDimensions.X())
      
      fMinPadDimensions = padDimensions;
  }
}

//_____________________________________________________________________________
void  AliMpSector::SetMaxPadIndices()
{
/// Set maximum pad indices in x, y

  if ( fMaxPadIndices != AliMpIntPair::Invalid() ) return;
  
  Int_t maxIndexInX = 0;
  Int_t maxIndexInY = 0;
  for (Int_t i=0; i<GetNofRows(); i++) {

    Int_t ixh = GetRow(i)->GetHighIndicesLimit().GetFirst();
    if ( ixh > maxIndexInX ) maxIndexInX = ixh;

    Int_t iyh = GetRow(i)->GetHighIndicesLimit().GetSecond();
    if ( iyh > maxIndexInY ) maxIndexInY = iyh;
  }  
  
  fMaxPadIndices = AliMpIntPair(maxIndexInX, maxIndexInY);
}


//_____________________________________________________________________________
void  AliMpSector::SetNofPads()
{
/// Set the total number of pads

  fNofPads = fMotifMap->CalculateNofPads();
}

//
// public methods
//

//_____________________________________________________________________________
void  AliMpSector::SetRowSegmentOffsets()
{
/// For all rows set the offset to all row segments.

  for (Int_t irow=0; irow<GetNofRows(); irow++)
    GetRow(irow)->SetRowSegmentOffsets(fOffset);    
}

//_____________________________________________________________________________
void AliMpSector::Initialize() 
{
/// Make needed settings after sector is read from
/// data files.

  SetRowOffsets();
  SetMotifPositions();
  SetGlobalIndices();
  SetMinPadDimensions();
  SetMaxPadIndices();
  SetNofPads();
}  

//_____________________________________________________________________________
void AliMpSector::PrintGeometry()  const
{
/// Print the positions of rows, rows segments

  for (Int_t i=0; i<GetNofRows(); i++) {
    AliMpRow* row = GetRow(i);
    
    cout << "ROW " << row->GetID() 
         << "  center Y " << row->Position().Y() << endl;

    for (Int_t j=0; j<row->GetNofRowSegments(); j++) {
       AliMpVRowSegment* rowSegment = row->GetRowSegment(j);
	
       cout << "   ROW Segment " << j 
            << "  borders " 
            << rowSegment->LeftBorderX() << "  "
            << rowSegment->RightBorderX()
            << "  x-size " 
            << 2*rowSegment->Dimensions().X() << "  "
	    << endl;
    }
  }
}     
      	     

//_____________________________________________________________________________
AliMpRow* AliMpSector::FindRow(const TVector2& position) const
{
/// Find the row for the specified y position.                              \n
/// If y is on border the lowest row is returned.
  
  Double_t y = position.Y();
  
#ifdef WITH_STL
  for (Int_t i=0; i<GetNofRows(); i++) {
    if ( y >= fRows[i]->LowBorderY() && y <= fRows[i]->UpperBorderY())
      return fRows[i];
  }    
#endif

#ifdef WITH_ROOT
  for (Int_t i=0; i<GetNofRows(); i++) {
    if ( y >= ((AliMpRow*)fRows[i])->LowBorderY() && 
         y <= ((AliMpRow*)fRows[i])->UpperBorderY())
      return (AliMpRow*)fRows[i];
  }    
#endif
  
  return 0;
}

//_____________________________________________________________________________
AliMpVMotif* AliMpSector::FindMotif(const TVector2& position) const
{
/// Find the motif in the specified position.                               \n
/// Return 0 if no motif is found.
  
  // Find the row segment
  AliMpVRowSegment* rowSegment = FindRowSegment(position);
  
  if (!rowSegment) return 0;
  
  // Find motif
  return rowSegment->FindMotif(position);  
}
//_____________________________________________________________________________
Int_t AliMpSector::FindMotifPositionId(const TVector2& position) const
{
/// Find the motif position ID in the specified position.                   \n
/// Return 0 if no motif is found.
 
  // Find the row segment
  AliMpVRowSegment* rowSegment = FindRowSegment(position);
  
  if (!rowSegment) return 0;
    
  // Find motif position ID
  return rowSegment->FindMotifPositionId(position);  
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

#ifdef WITH_STL
    AliMpRow* row = fRows[irow];
#endif
#ifdef WITH_ROOT
    AliMpRow* row = (AliMpRow*)fRows[irow];
#endif
   
    for (Int_t iseg=0; iseg<row->GetNofRowSegments(); iseg++) {
      AliMpVRowSegment* segment = row->GetRowSegment(iseg); 
      if (segment->HasMotifPosition(motifPositionId)) return segment;
    }
  }
  
  return 0;    
}

//_____________________________________________________________________________
TVector2  AliMpSector::FindPosition(Int_t motifPositionId) const
{
/// Find the position of the motif specified by its position Id.            \n
/// Return 0 if no row segment is found.

  AliMpVRowSegment* segment = FindRowSegment(motifPositionId);

  if (!segment) {
    AliWarningStream() << "Given motifPositionId not found." << endl;
    return TVector2();
  }   

  return segment->MotifCenter(motifPositionId);
}
   
//_____________________________________________________________________________
AliMpZone*  AliMpSector::FindZone(const TVector2& padDimensions) const
{
/// Find the zone with specified padDimensions.

  for (Int_t i=0; i<GetNofZones(); i++) {
    AliMpZone* zone = GetZone(i+1);
    if (AliMpConstants::IsEqual(padDimensions, zone->GetPadDimensions())) 
      return zone;
  }
  
  // Return 0 if not found
  return 0;  	 
}

//_____________________________________________________________________________
TVector2 AliMpSector::Position() const
{
/// Return the sector offset.

  return fOffset;
}  


//_____________________________________________________________________________
TVector2 AliMpSector::Dimensions() const
{
/// Return the maximum halflengths in x, y.

  Double_t x = 0.;
  Double_t y = 0.;
  for (Int_t i=0; i<GetNofRows(); i++) {

#ifdef WITH_STL
    // take the largest x row dimension
    if (fRows[i]->Dimensions().X() > x) 
      x = fRows[i]->Dimensions().X();
      
    // add all rows y dimensions  
    y += fRows[i]->Dimensions().Y();
#endif

#ifdef WITH_ROOT
    // take the largest x row dimension
    if ( ((AliMpRow*)fRows[i])->Dimensions().X() > x) 
      x = ((AliMpRow*)fRows[i])->Dimensions().X();
      
    // add all rows y dimensions  
    y += ((AliMpRow*)fRows[i])->Dimensions().Y();
#endif
  }
  
  return TVector2(x, y);  
}  

//_____________________________________________________________________________
Int_t AliMpSector::GetNofZones() const
{    
/// Return the number of zones.

#ifdef WITH_STL
  return fZones.size();
#endif

#ifdef WITH_ROOT
  return fZones.GetEntriesFast();
#endif
}  

//_____________________________________________________________________________
AliMpZone* AliMpSector::GetZone(Int_t zoneID) const
{
/// Return zone with specified ID.

  if (zoneID < 1 || zoneID > GetNofZones()) {
    AliWarningStream() << "Index outside range" << endl;
    return 0;
  }
  
#ifdef WITH_STL
  return fZones[zoneID-1];
#endif

#ifdef WITH_ROOT
  return (AliMpZone*)fZones[zoneID-1];
#endif
}  

//_____________________________________________________________________________
Int_t AliMpSector::GetNofRows() const
{
/// Return the number of rows.

#ifdef WITH_STL
  return fRows.size();
#endif

#ifdef WITH_ROOT
  return fRows.GetEntriesFast();
#endif
}  

//_____________________________________________________________________________
AliMpRow* AliMpSector::GetRow(Int_t rowID) const
{
/// Return row with specified ID.

  if (rowID < 0 || rowID >= GetNofRows()) {
    AliWarningStream() << "Index outside range" << endl;
    return 0;
  }
  
#ifdef WITH_STL
  return fRows[rowID];
#endif

#ifdef WITH_ROOT
  return (AliMpRow*)fRows[rowID];
#endif
}

//_____________________________________________________________________________
AliMp::PlaneType
AliMpSector::GetPlaneType() const
{
/// Return the plane type

  return GetDirection()==AliMp::kY ? AliMp::kBendingPlane : AliMp::kNonBendingPlane;
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
