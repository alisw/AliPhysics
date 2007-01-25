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
// $MpId: AliMpSectorSegmentation.cxx,v 1.15 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpSectorSegmentation
// -----------------------------
// Class describing the segmentation of the sector.        
// Provides methods related to pads:
// conversion between pad indices, pad location, pad position;
// finding pad neighbour.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpSectorSegmentation.h"
#include "AliMpSector.h"
#include "AliMpZone.h"
#include "AliMpSubZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpVMotif.h"
#include "AliMpMotifPosition.h"
#include "AliMpConnection.h"
#include "AliMpNeighboursPadIterator.h"
#include "AliMpSectorAreaHPadIterator.h"
#include "AliMpSectorAreaVPadIterator.h"
#include "AliMpIntPair.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMpSectorSegmentation)
/// \endcond

#ifdef WITH_ROOT
const Double_t AliMpSectorSegmentation::fgkS1 = 100000.;
const Double_t AliMpSectorSegmentation::fgkS2 = 1000.;
#endif

//______________________________________________________________________________
AliMpSectorSegmentation::AliMpSectorSegmentation(
                            const AliMpSector* sector, Bool_t own) 
  : AliMpVSegmentation(),
    fkSector(sector),
    fIsOwner(own),
    fPadBuffer(0),
    fPadDimensionsMap(),
    fMaxIndexInX(0),
    fMaxIndexInY(0)
{
/// Standard constructor

  AliDebugStream(1) << "this = " << this << endl;

  fPadBuffer = new AliMpPad(AliMpPad::Invalid());
  
  FillPadDimensionsMap();
}

//______________________________________________________________________________
AliMpSectorSegmentation::AliMpSectorSegmentation() 
  : AliMpVSegmentation(),
    fkSector(0),
    fIsOwner(false),
    fPadBuffer(0),
    fPadDimensionsMap(),      
    fMaxIndexInX(0),
    fMaxIndexInY(0)
{
/// Default constructor

  AliDebugStream(1) << "this = " << this << endl;
}

//______________________________________________________________________________
AliMpSectorSegmentation::~AliMpSectorSegmentation() 
{
/// Destructor 

  AliDebugStream(1) << "this = " << this << endl;

  if ( fIsOwner ) delete fkSector;

  delete fPadBuffer;
  
}

//
// private methods
//

//_____________________________________________________________________________
void 
AliMpSectorSegmentation::GetAllElectronicCardIDs(TArrayI& ecn) const
{
  GetSector()->GetAllMotifPositionsIDs(ecn);
}

#ifdef WITH_ROOT
//______________________________________________________________________________
Long_t AliMpSectorSegmentation::GetIndex(const TVector2& vector2) const
{
/// Convert the two vector to long.

  return Long_t(TMath::Floor((vector2.X()*fgkS1 + vector2.Y())*fgkS2));
}  

//______________________________________________________________________________
TVector2  AliMpSectorSegmentation::GetVector(Long_t index) const
{
/// Convert the long index to twovector.

  return TVector2( TMath::Floor(index/fgkS1)/fgkS2,
                   (index - TMath::Floor(index/fgkS1)*fgkS1)/fgkS2 );
}  
#endif

//______________________________________________________________________________
void AliMpSectorSegmentation::FillPadDimensionsMap()
{
/// Fill the maps between zone ids and pad dimensions.

  for (Int_t i=0; i<fkSector->GetNofZones(); i++) {
    AliMpZone* zone   = fkSector->GetZone(i+1);
    Int_t  zoneID = zone->GetID();
    
    if (!AliMpConstants::IsEqual(zone->GetPadDimensions(), TVector2())) {

      // regular zone
#ifdef WITH_STL
      fPadDimensionsMap[zoneID*10] = zone->GetPadDimensions();
#endif
#ifdef WITH_ROOT
     AliDebugStream(3)
       << "Filling fPadDimensions[" << zoneID*10 << "] = ("
       << zone->GetPadDimensions().X() << ", "
       << zone->GetPadDimensions().Y() << ")" << endl;

     fPadDimensionsMap.Add((Long_t)(zoneID*10), 
                            GetIndex(zone->GetPadDimensions()));
#endif
    }
    else {
      // special zone
      Int_t subIndex = 0;
      for (Int_t j=0; j<zone->GetNofSubZones(); j++) {
        AliMpSubZone* subZone = zone->GetSubZone(j);
	AliMpVMotif*  motif = subZone->GetMotif();
	
	for (Int_t k=0; k<motif->GetNofPadDimensions(); k++) {
	  Int_t index = zoneID*10 +  subIndex++;
#ifdef WITH_STL
          fPadDimensionsMap[index] = motif->GetPadDimensions(k);
#endif
#ifdef WITH_ROOT
          AliDebugStream(3)
            << "Filling fPadDimensions[" << index << "] = ("
            << motif->GetPadDimensions(k).X() << ", "
            << motif->GetPadDimensions(k).Y() << ") motif "
	    << motif->GetID().Data() << "-" << k << endl;

          fPadDimensionsMap.Add((Long_t)(index), 
                            GetIndex(motif->GetPadDimensions(k)));
#endif
	}
      }	  
    }	  
  }      
}

//______________________________________________________________________________
AliMpMotifPosition* 
AliMpSectorSegmentation::FindMotifPosition(const AliMpIntPair& indices) const
{
/// Find the motif position which contains the given pad indices
/// return 0 if not found

  switch (fkSector->GetDirection()) {
    case AliMp::kX : {
    // Case where all the pads have the same size along X direction

      for (Int_t irow=0; irow<fkSector->GetNofRows(); ++irow) {
        AliMpRow* row = fkSector->GetRow(irow);
        if (row->GetLowIndicesLimit().GetFirst()<=indices.GetFirst() &&
            row->GetHighIndicesLimit().GetFirst()>=indices.GetFirst()) {
            
           for (Int_t iseg=0;iseg<row->GetNofRowSegments();++iseg){
            AliMpVRowSegment* seg = row->GetRowSegment(iseg);
            if (seg->GetLowIndicesLimit().GetFirst()<=indices.GetFirst() &&
                seg->GetHighIndicesLimit().GetFirst()>=indices.GetFirst()) {

              AliMpMotifPosition* motifPos;
              for (Int_t imot=0;imot<seg->GetNofMotifs();++imot) {
                motifPos 
		  = fkSector->GetMotifMap()
		    ->FindMotifPosition(seg->GetMotifPositionId(imot));
                if (motifPos && motifPos->HasPad(indices)) return motifPos;
              }
            }
          }
        }
      }
      return 0;
    }
    break;
    ////////////////////////////////////////////////////////////////////////////////
    case AliMp::kY : {
      // Case where all the pads have the same size along Y direction   
      // look for the row which contains the indices
      AliMpRow* row=0;
      Int_t irow;
      for (irow=0; irow<fkSector->GetNofRows(); ++irow) {
        row = fkSector->GetRow(irow);
        AliMpVRowSegment* lastSeg = row->GetRowSegment(row->GetNofRowSegments()-1);
        if (lastSeg->GetLowIndicesLimit().GetSecond()<=indices.GetSecond() &&
            lastSeg->GetHighIndicesLimit().GetSecond()>=indices.GetSecond()) break;
        // NOTE : We use the last row segment in order to ensure that
        // we are not on a special motif
      }
      if (irow==fkSector->GetNofRows()) return 0;
      // look for the row segment, in the found row, which contains the indices
      AliMpVRowSegment* seg=0;
      Int_t iseg;
      for (iseg=0;iseg<row->GetNofRowSegments();++iseg){
        seg = row->GetRowSegment(iseg);
        if (seg->HasIndices(indices)) break;
      }
      if (iseg==row->GetNofRowSegments()) return 0;
  
      // look for the motif position which contains the indices
      AliMpMotifPosition* motifPos=0;
      Int_t imot=0;
      for (imot=0;imot<seg->GetNofMotifs();++imot) {
        motifPos 
	  = fkSector->GetMotifMap()
	    ->FindMotifPosition(seg->GetMotifPositionId(imot));
        if (motifPos && motifPos->HasPad(indices)) break;
      }      
      if (imot==seg->GetNofMotifs()) return 0;
   
      return motifPos;      
    }
    default: return 0;
  }
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByXDirection(const TVector2& startPosition, 
                                         Double_t maxX) const
{
/// Find the first valid pad from starting position in the
/// direction of pad lines up to distance dx.

  // Define step limits
  Double_t  stepX = fkSector->GetMinPadDimensions().X();
 
  // Search in X direction
  AliMpPad pad;
  TVector2 position(startPosition);    
  do {
    pad = PadByPosition(position, false);
    position += TVector2(stepX, 0.);
  }   
  while ( !pad.IsValid() && position.X() < maxX ); 
  
  // Invalidate pad if it is outside limits
  if ((pad.Position().X() - pad.Dimensions().X()) > maxX) 
    pad = AliMpPad::Invalid();

  return pad;
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByYDirection(const TVector2& startPosition, 
                                         Double_t maxY) const
{
/// Find the first valid pad from starting position in the
/// direction of pad columns up to distance dx.
  
  // Define step limits
  Double_t stepY = fkSector->GetMinPadDimensions().Y();
 
  // Search in Y direction
  AliMpPad pad;
  TVector2 position(startPosition);    
  do {
    pad = PadByPosition(position, false);
    position += TVector2(0., stepY);
  }   
  while ( !pad.IsValid() && position.Y() < maxY ); 
  
  // Invalidate pad if it is outside limits
  if ((pad.Position().Y() - pad.Dimensions().Y()) > maxY) 
    pad = AliMpPad::Invalid();

  return pad;
}

//
// public methods
//

//______________________________________________________________________________
AliMpVPadIterator* 
AliMpSectorSegmentation::CreateIterator(const AliMpArea& area) const
{
/// Create the area iterator. 

  switch (fkSector->GetDirection()) {
  
    case AliMp::kX: return new AliMpSectorAreaVPadIterator(this, area);
             ;;
    case AliMp::kY: return new AliMpSectorAreaHPadIterator(this, area);
             ;;
  }
  
  Fatal("CreateIterator", "Incomplete switch on Sector direction");
  return 0;  
}   
  
//______________________________________________________________________________
AliMpVPadIterator* 
AliMpSectorSegmentation::CreateIterator(const AliMpPad& centerPad,
                                        Bool_t includeCenter) const
{
/// Create the neighbours pad iterator.

  return new AliMpNeighboursPadIterator(this, centerPad, includeCenter);
}   
  
//______________________________________________________________________________
TVector2
AliMpSectorSegmentation::Dimensions() const
{
  return GetSector()->Dimensions();
}

//______________________________________________________________________________
AliMp::PlaneType
AliMpSectorSegmentation::PlaneType() const
{
  return GetSector()->GetPlaneType();
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByLocation(const AliMpIntPair& location, 
                                       Bool_t warning) const
{
/// Find the pad which corresponds to the given location
  
  if ((*fPadBuffer).GetLocation()==location) return (*fPadBuffer);
  
  AliMpMotifPosition* motifPos = 
    fkSector->GetMotifMap()->FindMotifPosition(location.GetFirst());
  if (!motifPos){
    if (warning) Warning("PadByLocation","The pad motif position ID doesn't exists");
    return AliMpPad::Invalid();
  }
  
  AliMpVMotif* motif = motifPos->GetMotif();
  AliMpIntPair localIndices = 
    motif->GetMotifType()->FindLocalIndicesByGassiNum(location.GetSecond());
  if (! localIndices.IsValid()) {
    if (warning) Warning("PadByLocation","The pad number doesn't exists");
    return AliMpPad::Invalid();
  }
  TVector2 delta = motif->PadPositionLocal(localIndices);
  return (*fPadBuffer) = AliMpPad(location,
              motifPos->GlobalIndices(localIndices),
              motifPos->Position()+delta,
              motif->GetPadDimensions(localIndices));

}
//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByIndices(const AliMpIntPair& indices,
                                      Bool_t warning ) const
{
/// Find the pad which corresponds to the given indices  

  if ((*fPadBuffer).GetIndices()==indices) return (*fPadBuffer);    
   
  AliMpMotifPosition* motifPos = FindMotifPosition(indices);
  if (!motifPos) {    
    if (warning) 
      Warning("PadByIndices","Pad indices not contained in any motif!");
    return AliMpPad::Invalid();
  }
  
  // retrieve the local indices in the found motif
  AliMpVMotif* motif = motifPos->GetMotif();
  AliMpIntPair localIndices = indices - motifPos->GetLowIndicesLimit();
  
  AliMpConnection* connection=
    motif->GetMotifType()->FindConnectionByLocalIndices(localIndices);
    
  if (!connection){
    if (warning) Warning("PadByIndices","No connection with the given indices!");
    return AliMpPad::Invalid();
  }

  TVector2 localPos = motif->PadPositionLocal(localIndices);

  return (*fPadBuffer) 
    = AliMpPad(AliMpIntPair(motifPos->GetID(),connection->GetGassiNum()),
               indices,
               motifPos->Position()+localPos,
               motif->GetPadDimensions(localIndices)); 

}
//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByPosition(const TVector2& position,
                                       Bool_t warning) const
{
/// Find the pad which corresponds to the given position

  if ((*fPadBuffer).Position().X()==position.X() && 
      (*fPadBuffer).Position().Y()==position.Y()) return (*fPadBuffer);  

  Int_t motifPosID = fkSector->FindMotifPositionId(position);
  AliMpMotifPosition* motifPos 
    = fkSector->GetMotifMap()
        ->FindMotifPosition(motifPosID);
    
  if (!motifPos){
    if (warning) Warning("PadByPosition","Position outside limits");
    return AliMpPad::Invalid();
  }

  AliMpVMotif* motif =  motifPos->GetMotif();  
  AliMpIntPair localIndices 
    = motif->PadIndicesLocal(position-motifPos->Position());
    
  AliMpConnection* connect = 
    motif->GetMotifType()->FindConnectionByLocalIndices(localIndices);

   if (!connect){
    if (warning) Warning("PadByPosition","Position outside motif limits");
    return AliMpPad::Invalid();
  }
  
  return (*fPadBuffer)
    = AliMpPad(AliMpIntPair(motifPosID,connect->GetGassiNum()),
               motifPos->GlobalIndices(localIndices),
               motifPos->Position()+motif->PadPositionLocal(localIndices),
               motif->GetPadDimensions(localIndices));

}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByDirection(const TVector2& startPosition, 
                                        Double_t distance) const
{
/// Find the first valid pad from starting position in the
/// direction of pad lines/columns up to the specified distance.
/// Pad lines are the lines of pads in the sector with constant pad y size,
/// pad columns are the columns of pads in the sector with constant pad x size. 

  switch (fkSector->GetDirection()) {
  
    case AliMp::kX: return PadByYDirection(startPosition, distance);
             ;;
    case AliMp::kY: return PadByXDirection(startPosition, distance);
             ;;
  }
  
  Fatal("PadByDirection", "Incomplete switch on Sector direction");
  return AliMpPad::Invalid();  
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::MaxPadIndexX() const
{
/// Return maximum pad index in x

  return fkSector->GetMaxPadIndices().GetFirst();
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::MaxPadIndexY() const
{
/// Return maximum pad index in y

  return fkSector->GetMaxPadIndices().GetSecond();
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::NofPads() const
{
/// Return number of pads defined in the sector

  return fkSector->GetNofPads();
}

//______________________________________________________________________________
Bool_t AliMpSectorSegmentation::HasPad(const AliMpIntPair& indices) const
{
/// Does the pad specified by \a indices exist ?

  return PadByIndices(indices,kFALSE) != AliMpPad::Invalid();
}

//______________________________________________________________________________
Bool_t AliMpSectorSegmentation::HasMotifPosition(Int_t motifPositionID) const
{
/// Does the motif position specified by motifPositionID exist ?

  return (fkSector->GetMotifMap()->FindMotifPosition(motifPositionID) != 0);
}

//______________________________________________________________________________
TVector2  AliMpSectorSegmentation::GetMinPadDimensions() const
{
/// Returne the dimensions of the smallest pad.

  return fkSector->GetMinPadDimensions();
}  

//______________________________________________________________________________
Int_t AliMpSectorSegmentation::Zone(const AliMpPad& pad, Bool_t warning) const
{
/// Return the zone index of the zone containing the specified pad.
/// This zone index is different from the zone ID,
/// as it is unique for each pad dimensions.
/// It is composed in this way:
///   zoneID*10 + specific index 
/// Specific index is present only for zones containing special motifs.

  if (!pad.IsValid()) {
    if (warning) Warning("Zone(AliMpPad)", "Invalid pad");
    return 0;
  }  

#ifdef WITH_STL
  PadDimensionsMapCIterator it;
  for (it = fPadDimensionsMap.begin(); it != fPadDimensionsMap.end(); ++it) {
    if (AliMpConstants::IsEqual(it->second, pad.Dimensions()))
      return it->first;
  }
#endif

#ifdef WITH_ROOT
  PadDimensionsMapCIterator it(&fPadDimensionsMap);
  Long_t key, value;
  while ( it.Next(key, value) ) {
    TVector2 dimensions =  GetVector(value);
    if (AliMpConstants::IsEqual(dimensions, pad.Dimensions()))
      return (Int_t)key;
  } 
  
  AliError(Form("fPadDimensionsMap size is %d",fPadDimensionsMap.GetSize()));
  
#endif

  // Should never happen
  AliErrorStream() 
    << "Zone(AliMpPad pad) not found, where pad is: " << pad << endl;
  return 0;
}  

//______________________________________________________________________________
TVector2 
AliMpSectorSegmentation::PadDimensions(Int_t zone, Bool_t warning) const
{
/// Return the pad dimensions for the zone with the specified zone index.

#ifdef WITH_STL
  PadDimensionsMapCIterator it = fPadDimensionsMap.find(zone);
  if (it != fPadDimensionsMap.end()) return it->second;
#endif

#ifdef WITH_ROOT
  Long_t value = fPadDimensionsMap.GetValue(zone);
  if (value) return GetVector(value);
#endif

  if (warning) Warning("PadDimensions(zone)", "not found");
  return TVector2();
}  

//______________________________________________________________________________
Bool_t AliMpSectorSegmentation::CircleTest(const AliMpIntPair& indices) const
{
/// Verify that all methods for retrieving pads are consistents between them.
/// Return true if the pad with specified indices was found and verified,
/// false otherwise.

  if (!HasPad(indices)) return false;

  // Verify the indice->location->position->indice way
  AliMpIntPair location = PadByIndices(indices).GetLocation();
  TVector2 position = PadByLocation(location).Position();
  AliMpIntPair retIndices = PadByPosition(position).GetIndices();
    
  if (retIndices != indices) {
    cout << "Pad " << indices << " lead to inconsistency" << endl;
    cout << "in indice->location->position->indice way..." << endl;
    cout << "starting from " << indices << "-->" << location << "-->" 
         << '(' << position.X() << ',' << position.Y() << ')'
         << " and retIndices: " << retIndices << endl;
  }
    
    
  // Verify the indice->position->location->indice way    
  position = PadByIndices(indices).Position();
  location = PadByPosition(position).GetLocation();
  retIndices = PadByLocation(location).GetIndices();

  if (retIndices != indices) {
    cout << "Pad " << indices << " lead to inconsistency" << endl;
    cout << "in indice->position->location->indice way..." <<endl;
    cout << "starting from " << indices 
         << " and retIndices: " << retIndices << endl;
  }
  
  return true;
}

//______________________________________________________________________________
void
AliMpSectorSegmentation::Print(Option_t* opt) const
{
  fkSector->Print(opt);
}

//______________________________________________________________________________
void AliMpSectorSegmentation::PrintZones() const
{
/// Print all zones and pads dimensions from the map.

  cout << "Zones: " << endl;

#ifdef WITH_STL
  PadDimensionsMapCIterator it;
  for (it = fPadDimensionsMap.begin(); it != fPadDimensionsMap.end(); ++it) {
    cout << "    zone: " <<   setw(4) << it->first;
    cout << "    pad dimensions: ( " 
         << it->second.X() << ", " << it->second.Y() << ")" << endl; 
  }
#endif

#ifdef WITH_ROOT
  PadDimensionsMapCIterator it(&fPadDimensionsMap);
  Long_t key, value;
  while ( it.Next(key, value) ) {
    //cout << "Iterating over: " << key << ", " << value << endl;
    TVector2 dimensions =  GetVector(value);

    cout << "    zone: " <<   setw(4) << key;
    cout << "    pad dimensions: ( " 
         << dimensions.X() << ", " << dimensions.Y() << ")" << endl; 
  }
#endif
}
  
