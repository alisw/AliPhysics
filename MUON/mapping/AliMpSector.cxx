// $Id$
// Category: sector
//
// Class AliMpSector
// -----------------
// Class describing the sector of the MUON chamber of station 1.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>
#include <Riostream.h>

#include "AliMpSector.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"
#include "AliMpMotifMap.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"

ClassImp(AliMpSector)

//_____________________________________________________________________________
AliMpSector::AliMpSector(const TString& id, Int_t nofZones, Int_t nofRows, 
                         AliMpDirection direction) 
  : TObject(),
    fID(id),
    fOffset(TVector2(0., 0.)),
    fZones(),
    fRows(),
    fDirection(direction),
    fMinPadDimensions(TVector2(1.e6, 1.e6))
{
//
  fMotifMap = new AliMpMotifMap();

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
AliMpSector::AliMpSector(const AliMpSector& right) 
  : TObject(right) {
// 
  Fatal("AliMpSector", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpSector::AliMpSector() 
  : TObject(),
    fID(""),    
    fOffset(TVector2(0., 0.)),
    fZones(),
    fRows(),
    fMotifMap(0),
    fDirection(kX),
    fMinPadDimensions(TVector2(0., 0.))
{
//
}

//_____________________________________________________________________________
AliMpSector::~AliMpSector() {
// 
  // deletes 
  for (Int_t izone = 0; izone<GetNofZones(); izone++) 
    delete fZones[izone];
    
  for (Int_t irow = 0; irow<GetNofRows(); irow++) 
    delete fRows[irow];

  delete fMotifMap;
}

//
// operators
//

//_____________________________________________________________________________
AliMpSector& AliMpSector::operator=(const AliMpSector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// private methods
//

//_____________________________________________________________________________
AliMpVPadIterator* AliMpSector::CreateIterator() const
{
  return new AliMpSectorPadIterator(this);
}


//_____________________________________________________________________________
AliMpVRowSegment* AliMpSector::FindRowSegment(const TVector2& position) const
{
// Finds the row segment in the specified position.
// Returns 0 if no motif is found.
// ---
  
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
// For each row checks consitency of the row segments
// and calculates the row offset.
// ---

  Double_t offset = fOffset.Y();
  
  for (Int_t irow=0; irow<GetNofRows(); irow++)
    offset = GetRow(irow)->SetOffsetY(offset);    
}

//_____________________________________________________________________________
void  AliMpSector::SetMotifPositions()
{
// Creates motif positions objects and fills them in the motif map.
// ---

  for (Int_t i=0; i<GetNofRows(); i++)
    GetRow(i)->SetMotifPositions();
}

//_____________________________________________________________________________
void  AliMpSector::SetGlobalIndices()
{
// Set the indices limits to all indexed elements
// (row, row segment, motif positions).
// ---

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
// Sets the minimal pad dimensions.
// ---

  for (Int_t i=1; i<GetNofZones()+1; i++) {
    TVector2 padDimensions = GetZone(i)->GetPadDimensions();
    
    if ( fDirection == kX &&  
         padDimensions.Y() > 0. && padDimensions.Y() < fMinPadDimensions.Y() ||
         fDirection == kY && 
	 padDimensions.X() > 0. && padDimensions.X() < fMinPadDimensions.X())
      
      fMinPadDimensions = padDimensions;
  }
}

//
// public methods
//

//_____________________________________________________________________________
void  AliMpSector::SetRowSegmentOffsets()
{
// For all rows sets offset to all row segments.
// ---

  for (Int_t irow=0; irow<GetNofRows(); irow++)
    GetRow(irow)->SetRowSegmentOffsets(fOffset);    
}

//_____________________________________________________________________________
void AliMpSector::Initialize() 
{
// Makes needed settings after sector is read from
// data files.
// ---

  SetRowOffsets();
  SetMotifPositions();
  SetGlobalIndices();
  SetMinPadDimensions();
}  

//_____________________________________________________________________________
void AliMpSector::PrintGeometry()  const
{
// Prints the positions of rows, rows segments
// ---

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
// Finds the row for the specified y position.
// If y is on border the lowest row is returned.
// ---
  
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
// Finds the motif in the specified position.
// Returns 0 if no motif is found.
// ---
  
  // Find the row segment
  AliMpVRowSegment* rowSegment = FindRowSegment(position);
  
  if (!rowSegment) return 0;
  
  // Find motif
  return rowSegment->FindMotif(position);  
}

//_____________________________________________________________________________
Int_t AliMpSector::FindMotifPositionId(const TVector2& position) const
{
// Finds the motif position ID in the specified position.
// Returns 0 if no motif is found.
// ---
  
  // Find the row segment
  AliMpVRowSegment* rowSegment = FindRowSegment(position);
  
  if (!rowSegment) return 0;
    
  // Find motif position ID
  return rowSegment->FindMotifPositionId(position);  
}

//_____________________________________________________________________________
AliMpRow* AliMpSector::FindRow(Int_t motifPositionId) const
{
// Finds the row with the the specified motif position.
// Returns 0 if no row is found.
// ---

  AliMpVRowSegment* segment = FindRowSegment(motifPositionId);
  
  if (segment) return segment->GetRow();
  
  return 0;  
}

//_____________________________________________________________________________
AliMpVRowSegment* AliMpSector::FindRowSegment(Int_t motifPositionId) const
{
// Finds the row segment with the the specified motif position.
// Returns 0 if no row segment is found.
// ---

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
// Finds the position of the motif specified by its position Id.
// Returns 0 if no row segment is found.
// ---

  AliMpVRowSegment* segment = FindRowSegment(motifPositionId);

  if (!segment) {
    Warning("FindPosition", "Given motifPositionId not found.");
    return TVector2();
  }   

  return segment->MotifCenter(motifPositionId);
}
   
//_____________________________________________________________________________
AliMpZone*  AliMpSector::FindZone(const TVector2& padDimensions) const
{
// Finds the zone with specified padDimensions.
// ---

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
// Returns the offset.
// ---

  return fOffset;
}  


//_____________________________________________________________________________
TVector2 AliMpSector::Dimensions() const
{
// Returns the maximum halflength in x, y.
// ---

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
// Returns the number of zones.
// ---

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
// Returns zone with specified ID.
// ---

  if (zoneID < 1 || zoneID > GetNofZones()) {
    Warning("GetZone", "Index outside range");
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
// Returns the number of rows.
// ---

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
// Returns row with specified ID.
// ---

  if (rowID < 0 || rowID >= GetNofRows()) {
    Warning("GetRow", "Index outside range");
    return 0;
  }
  
#ifdef WITH_STL
  return fRows[rowID];
#endif

#ifdef WITH_ROOT
  return (AliMpRow*)fRows[rowID];
#endif
}
