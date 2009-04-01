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

//-----------------------------------------------------------------------------
// Class AliMpSectorSegmentation
// -----------------------------
// Class describing the segmentation of the sector.        
// Provides methods related to pads:
// conversion between pad indices, pad location, pad position;
// finding pad neighbour.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

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
#include "AliMpSectorAreaHPadIterator.h"
#include "AliMpSectorAreaVPadIterator.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"
#include "AliMpEncodePair.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMpSectorSegmentation)
/// \endcond

//______________________________________________________________________________
AliMpSectorSegmentation::AliMpSectorSegmentation(
                            const AliMpSector* sector, Bool_t own) 
  : AliMpVSegmentation(),
    fkSector(sector),
    fIsOwner(own),
    fPadBuffer(0),
    fMaxIndexInX(0),
    fMaxIndexInY(0)
{
/// Standard constructor

  AliDebugStream(1) << "this = " << this << endl;

  fPadBuffer = new AliMpPad(AliMpPad::Invalid());
  
  //FillPadDimensionsMap();
}

//______________________________________________________________________________
AliMpSectorSegmentation::AliMpSectorSegmentation() 
  : AliMpVSegmentation(),
    fkSector(0),
    fIsOwner(false),
    fPadBuffer(0),
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

//______________________________________________________________________________
AliMpMotifPosition* 
AliMpSectorSegmentation::FindMotifPosition(Int_t ix, Int_t iy) const
{
/// Find the motif position which contains the given pad indices
/// return 0 if not found

  switch ( fkSector->GetDirection() ) {
    case AliMp::kX : {
    // Case where all the pads have the same size along X direction

      for ( Int_t irow=0; irow<fkSector->GetNofRows(); ++irow ) {
        AliMpRow* row = fkSector->GetRow(irow);
        if ( row->GetLowLimitIx() <= ix &&
             row->GetHighLimitIx()>= ix ) {
            
          for ( Int_t iseg=0;iseg<row->GetNofRowSegments();++iseg ) {
            AliMpVRowSegment* seg = row->GetRowSegment(iseg);
            if ( seg->GetLowLimitIx() <= ix &&
                 seg->GetHighLimitIx() >= ix ) {

              AliMpMotifPosition* motifPos;
              for ( Int_t imot=0;imot<seg->GetNofMotifs();++imot ) {
                motifPos 
		  = fkSector->GetMotifMap()
		    ->FindMotifPosition(seg->GetMotifPositionId(imot));
                if (motifPos && motifPos->HasPadByIndices(AliMp::Pair(ix,iy))) return motifPos;
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
      for ( irow=0; irow<fkSector->GetNofRows(); ++irow ) {
        row = fkSector->GetRow(irow);
        AliMpVRowSegment* lastSeg = row->GetRowSegment(row->GetNofRowSegments()-1);
        if ( lastSeg->GetLowLimitIy() <= iy &&
             lastSeg->GetHighLimitIy() >= iy ) break;
        // NOTE : We use the last row segment in order to ensure that
        // we are not on a special motif
      }
      if ( irow==fkSector->GetNofRows() ) return 0;
      // look for the row segment, in the found row, which contains the indices
      AliMpVRowSegment* seg=0;
      Int_t iseg;
      for ( iseg=0;iseg<row->GetNofRowSegments();++iseg ) {
        seg = row->GetRowSegment(iseg);
        if (seg->HasIndices(AliMp::Pair(ix, iy))) break;
      }
      if ( iseg==row->GetNofRowSegments() ) return 0;
  
      // look for the motif position which contains the indices
      AliMpMotifPosition* motifPos=0;
      Int_t imot=0;
      for ( imot=0;imot<seg->GetNofMotifs();++imot ) {
        motifPos 
	  = fkSector->GetMotifMap()
	    ->FindMotifPosition(seg->GetMotifPositionId(imot));
        if (motifPos && motifPos->HasPadByIndices(AliMp::Pair(ix, iy))) break;
      }      
      if (imot==seg->GetNofMotifs()) return 0;
   
      return motifPos;      
    }
    default: return 0;
  }
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByXDirection(Double_t startx, Double_t starty, 
                                         Double_t maxX) const
{
/// Find the first valid pad from starting position in the
/// direction of pad lines up to distance dx.

  // Define step limits
  Double_t  stepX = fkSector->GetMinPadDimensionX();
 
  // Search in X direction
  AliMpPad pad;
  Double_t posx = startx;
  do {
    pad = PadByPosition(posx, starty, false);
    posx += stepX;
  }   
  while ( ! pad.IsValid() && 
            posx - fkSector->GetMaxPadDimensionX() < maxX ); 
  
  // Invalidate pad if it is outside limits
  if ( ( pad.GetPositionX() - pad.GetDimensionX()) > maxX ) 
    pad = AliMpPad::Invalid();

  return pad;
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByYDirection(Double_t startx, Double_t starty, 
                                         Double_t maxY) const
{
/// Find the first valid pad from starting position in the
/// direction of pad columns up to distance dx.
  
  // Define step limits
  Double_t stepY = fkSector->GetMinPadDimensionY();
 
  // Search in Y direction
  AliMpPad pad;
  Double_t posy = starty;
  do {
    pad = PadByPosition(startx, posy, false);
    posy += stepY;
  }   
  while ( ! pad.IsValid() && 
            posy - fkSector->GetMaxPadDimensionY()< maxY ); 
  
  // Invalidate pad if it is outside limits
  if (( pad.GetPositionY() - pad.GetDimensionY()) > maxY ) 
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
AliMpSectorSegmentation::CreateIterator() const
{
/// Create the sector iterator

  return new AliMpSectorPadIterator(fkSector);
}

//______________________________________________________________________________
Int_t 
AliMpSectorSegmentation::GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                                       Bool_t includeSelf,
                                       Bool_t includeVoid) const
{
  /// Uses default implementation
  return AliMpVSegmentation::GetNeighbours(pad,neighbours,includeSelf,includeVoid);
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByLocation(Int_t manuId, Int_t manuChannel, 
                                       Bool_t warning) const
{
/// Find the pad which corresponds to the given location
 
  if ( fPadBuffer->GetManuId() == manuId &&
       fPadBuffer->GetManuChannel() == manuChannel ) return (*fPadBuffer);
  
  AliMpMotifPosition* motifPos = 
    fkSector->GetMotifMap()->FindMotifPosition(manuId);
  if (!motifPos){
    if (warning) Warning("PadByLocation","The pad motif position ID doesn't exists");
    return AliMpPad::Invalid();
  }
  
  AliMpVMotif* motif = motifPos->GetMotif();
  MpPair_t localIndices = 
    motif->GetMotifType()->FindLocalIndicesByGassiNum(manuChannel);
  if ( localIndices < 0 ) {
    if (warning) Warning("PadByLocation","The pad number doesn't exists");
    return AliMpPad::Invalid();
  }

  Double_t posx, posy;
  motif->PadPositionLocal(localIndices, posx, posy);
  posx += motifPos->GetPositionX();
  posy += motifPos->GetPositionY();

  Double_t dx, dy;
  motif->GetPadDimensionsByIndices(localIndices, dx, dy);

  return (*fPadBuffer) = AliMpPad(manuId, manuChannel,
              motifPos->GlobalIndices(localIndices),
              posx, posy, dx, dy);
}
//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByIndices(Int_t ix, Int_t iy, Bool_t warning ) const
{
/// Find the pad which corresponds to the given indices  

  if ( fPadBuffer->GetIx() == ix &&
       fPadBuffer->GetIy() == iy ) return (*fPadBuffer);    
       
  MpPair_t indices = AliMp::Pair(ix, iy);     
  AliMpMotifPosition* motifPos = FindMotifPosition(ix, iy);
  if (!motifPos) {    
    if (warning) 
      Warning("PadByIndices","Pad indices not contained in any motif!");
    return AliMpPad::Invalid();
  }
  
  // retrieve the local indices in the found motif
  AliMpVMotif* motif = motifPos->GetMotif();
  MpPair_t localIndices = indices - motifPos->GetLowIndicesLimit();
  
  AliMpConnection* connection=
    motif->GetMotifType()->FindConnectionByLocalIndices(localIndices);
    
  if (!connection){
    if (warning) Warning("PadByIndices","No connection with the given indices!");
    return AliMpPad::Invalid();
  }

  Double_t posx, posy;
  motif->PadPositionLocal(localIndices, posx, posy);
  posx += motifPos->GetPositionX();
  posy += motifPos->GetPositionY();

  Double_t dx, dy;
  motif->GetPadDimensionsByIndices(localIndices, dx, dy);

  return (*fPadBuffer) 
    = AliMpPad(motifPos->GetID(),connection->GetManuChannel(),
               ix, iy, posx, posy, dx, dy);
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByPosition(Double_t x, Double_t y,
                                       Bool_t warning) const
{
/// Find the pad which corresponds to the given position

  if (fPadBuffer->GetPositionX()==x && 
      fPadBuffer->GetPositionY()==y) return (*fPadBuffer);  

  Int_t motifPosID = fkSector->FindMotifPositionId(x,y);
  AliMpMotifPosition* motifPos 
    = fkSector->GetMotifMap()
        ->FindMotifPosition(motifPosID);
    
  if (!motifPos){
    if (warning) Warning("PadByPosition","Position outside limits");
    return AliMpPad::Invalid();
  }

  AliMpVMotif* motif =  motifPos->GetMotif();  
  MpPair_t localIndices 
    = motif->PadIndicesLocal(x-motifPos->GetPositionX(),
                             y-motifPos->GetPositionY());

  if ( localIndices < 0 ) {
    if (warning) Warning("PadByPosition","Position outside motif limits");
    return AliMpPad::Invalid();
  }
    
  AliMpConnection* connect = 
    motif->GetMotifType()->FindConnectionByLocalIndices(localIndices);

  if ( ! connect ) {
    if (warning) Warning("PadByPosition","Position outside motif limits");
    return AliMpPad::Invalid();
  }

  Double_t posx, posy;
  motif->PadPositionLocal(localIndices, posx, posy);
  posx += motifPos->GetPositionX();
  posy += motifPos->GetPositionY();

  Double_t dx, dy;
  motif->GetPadDimensionsByIndices(localIndices, dx, dy);
  
  return (*fPadBuffer)
    = AliMpPad(motifPosID, connect->GetManuChannel(),
               motifPos->GlobalIndices(localIndices),
               posx, posy, dx, dy);
}

//______________________________________________________________________________
AliMpPad 
AliMpSectorSegmentation::PadByDirection(Double_t startx, Double_t starty, 
                                        Double_t distance) const
{
/// Find the first valid pad from starting position in the
/// direction of pad lines/columns up to the specified distance.
/// Pad lines are the lines of pads in the sector with constant pad y size,
/// pad columns are the columns of pads in the sector with constant pad x size. 

  switch (fkSector->GetDirection()) {
  
    case AliMp::kX: return PadByYDirection(startx, starty, distance);
             ;;
    case AliMp::kY: return PadByXDirection(startx, starty, distance);
             ;;
  }
  
  Fatal("PadByDirection", "Incomplete switch on Sector direction");
  return AliMpPad::Invalid();  
}

//_____________________________________________________________________________
Bool_t 
AliMpSectorSegmentation::HasPadByIndices(Int_t ix, Int_t iy) const
{
  ///  Whether or not we have a pad at indices=(ix,iy) 
  
  MpPair_t indices = AliMp::Pair(ix, iy);

  AliMpMotifPosition* motifPos = FindMotifPosition(ix, iy);
  
  if (motifPos) return motifPos->HasPadByIndices(indices);
  
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t 
AliMpSectorSegmentation::HasPadByLocation(Int_t manuId, Int_t manuChannel) const
{
  /// Whether or not we have a pad at location=(manuId,manuChannel)
  
  AliMpMotifPosition* motifPos 
    = fkSector->GetMotifMap()->FindMotifPosition(manuId);
  
  if ( motifPos ) return motifPos->HasPadByManuChannel(manuChannel);
  
  return kFALSE;
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::MaxPadIndexX() const
{
/// Return maximum pad index in x

  return AliMp::PairFirst(fkSector->GetMaxPadIndices());
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::MaxPadIndexY() const
{
/// Return maximum pad index in y

  return AliMp::PairSecond(fkSector->GetMaxPadIndices());
}

//______________________________________________________________________________
Int_t  AliMpSectorSegmentation::NofPads() const
{
/// Return number of pads defined in the sector

  return fkSector->GetNofPads();
}

//_____________________________________________________________________________
void 
AliMpSectorSegmentation::GetAllElectronicCardIDs(TArrayI& ecn) const
{
  /// Fill the array ecn with all manuIds

  GetSector()->GetAllMotifPositionsIDs(ecn);
}

//_____________________________________________________________________________
Int_t 
AliMpSectorSegmentation::GetNofElectronicCards() const
{
  /// Get the number of manus of this sector
  
  return fkSector->GetNofMotifPositions();  
}

//_____________________________________________________________________________
Bool_t 
AliMpSectorSegmentation::HasMotifPosition(Int_t manuId) const
{
  /// Whether we get a given manu. Uses default implementation
  return (AliMpVSegmentation::HasMotifPosition(manuId) != 0x0);
}

//_____________________________________________________________________________
AliMpMotifPosition* 
AliMpSectorSegmentation::MotifPosition(Int_t manuId) const
{
  /// Return a given manu
  return fkSector->GetMotifMap()->FindMotifPosition(manuId);
}

//______________________________________________________________________________
AliMp::PlaneType
AliMpSectorSegmentation::PlaneType() const
{
  return GetSector()->GetPlaneType();
}

//_____________________________________________________________________________
Double_t  
AliMpSectorSegmentation::GetDimensionX() const
{
/// Return sector x dimensions
  return GetSector()->GetDimensionX();
}

//_____________________________________________________________________________
Double_t  
AliMpSectorSegmentation::GetDimensionY() const
{
/// Return sector y dimensions
  return GetSector()->GetDimensionY();
}

//_____________________________________________________________________________
Double_t  
AliMpSectorSegmentation::GetPositionX() const
{
/// Return x position 
  return 0.;
}

//_____________________________________________________________________________
Double_t  
AliMpSectorSegmentation::GetPositionY() const
{
/// Return y position 
  return 0.;
}

//______________________________________________________________________________
void
AliMpSectorSegmentation::Print(Option_t* opt) const
{
/// Print the sector

  fkSector->Print(opt);
}

//______________________________________________________________________________
Double_t AliMpSectorSegmentation::GetMinPadDimensionX() const
{
/// Return the x dimension of the smallest pad.

  return fkSector->GetMinPadDimensionX();
}  


//______________________________________________________________________________
Double_t AliMpSectorSegmentation::GetMinPadDimensionY() const
{
/// Return the y dimension of the smallest pad.

  return fkSector->GetMinPadDimensionY();
}  


//______________________________________________________________________________
Bool_t AliMpSectorSegmentation::CircleTest(Int_t ix, Int_t iy) const
{
/// Verify that all methods for retrieving pads are consistents between them.
/// Return true if the pad with specified indices was found and verified,
/// false otherwise.

  if ( ! HasPadByIndices(ix, iy) ) return false;

  // Verify the indice->location->position->indice way
  AliMpPad pad1 = PadByIndices(ix, iy);  
  AliMpPad pad2 = PadByLocation(pad1.GetManuId(), pad1.GetManuChannel());
  AliMpPad pad3 = PadByPosition(pad2.GetPositionX(),pad2.GetPositionY());
                                      
  MpPair_t retIndices = pad3.GetIndices();
    
  if ( retIndices != AliMp::Pair(ix, iy) ) {
    cout << "Pad (" << ix << ',' << iy << ") lead to inconsistency" << endl;
    cout << "in indice->location->position->indice way..." << endl;
    cout << "starting from indices " << pad1 << endl
         << "--> location " << pad2 << endl
         << "--> position " 
         << '(' << pad3.GetPositionX() << ',' << pad3.GetPositionY() << ')'
         <<  endl << endl;
  }
    
  // Verify the indice->position->location->indice way  
  AliMpPad pad2bis = PadByPosition(pad1.GetPositionX(),pad1.GetPositionY());
  AliMpPad pad3bis = PadByLocation(pad2bis.GetManuId(), pad2bis.GetManuChannel());
  
  retIndices = pad3bis.GetIndices();

  if ( retIndices != AliMp::Pair(ix, iy) ) {
    cout << "Pad (" << ix << ',' << iy << ") lead to inconsistency" << endl;
    cout << "in indice->position->location->indice way..." << endl;
    cout << "starting from indices " << pad1 << endl
         << "--> position " 
         << '(' << pad2bis.GetPositionX() << ',' << pad2bis.GetPositionY() << ')' << endl
         << "--> location " << pad3bis
         << endl << endl;
  }
  
  return true;
}
