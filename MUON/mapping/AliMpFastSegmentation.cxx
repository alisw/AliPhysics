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

/// \class AliMpFastSegmentation
/// An implementation of AliMpVSegmentation, which uses
/// some internal maps to speed up the (Has)PadByIndices and PadByLocation
/// methods.
/// 
/// L. Aphecetche, Subatech
///

#include "AliMpFastSegmentation.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpPad.h"
#include "AliMpSector.h"
#include "AliMpSlat.h"
#include "AliMpVPadIterator.h"
#include <TArrayI.h>

/// \cond CLASSIMP
ClassImp(AliMpFastSegmentation)
/// \endcond

//#define CHECK

#ifdef CHECK
#include <cassert>
#endif

namespace
{
  ///
  /// The values in Encode and Encode2 are not exactly random.
  /// They are the result of a few "try and see" efforts to optimize the
  /// timing of the TExMap::GetValue (you should note that the TExMap implementation
  /// speed depends on the "non-uniformity" of the keys).
  ///
  /// So don't change those values w/o at least testing a bit the implications...
  /// But feel free to experiment though, in order to optimizer further ;-)
  ///
  Int_t Encode(Int_t a, Int_t b)
  {
    return a*1009 + b;
  }

  Int_t Encode2(Int_t a)
  {
    /// Ideally this method should be different for sectors and slats, as we have
    /// much less manus per DE for slats, and hence the "non-uniformity" is less...
    return ( a ^ (1<<10) ) << 16 ;
  }
}

//_____________________________________________________________________________
AliMpFastSegmentation::AliMpFastSegmentation(AliMpVSegmentation* vseg)
: AliMpVSegmentation(),
fHelper(vseg),
fMotifPositions(),
fIxIy(),
fManuId(),
fPosition()
{
  /// Ctor. We adopt vseg.
  
  AliCodeTimerAuto(vseg->ClassName());
  
  if (!vseg) 
  {
    AliError("Will get a hard time working with a NULL vseg !");
    return;
  }
  
  fPosition = vseg->Position();
  
  TArrayI manus;
  
  vseg->GetAllElectronicCardIDs(manus);
  
  for ( Int_t i = 0; i < manus.GetSize(); ++i ) 
  {
    Int_t manuId = manus[i];
    
    AliMpMotifPosition* mp = vseg->MotifPosition(manuId);
    
    // Should never happen
    if ( ! mp ) {
      AliFatal("AliMpMotifPosition not found.");
    }  
    
    Int_t index = 1 + fMotifPositions.GetLast();

    fMotifPositions.AddLast(mp);

    fManuId.Add(Encode2(manuId),1+index);

    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      if ( vseg->HasPadByLocation(AliMpIntPair(manuId,manuChannel)) )
      {
        AliMpPad pad = vseg->PadByLocation(AliMpIntPair(manuId,manuChannel));
        
        fIxIy.Add(Encode(pad.GetIndices().GetFirst(),pad.GetIndices().GetSecond()),1+index);
      }
    }
  }
}

//_____________________________________________________________________________
AliMpFastSegmentation::~AliMpFastSegmentation()
{
  /// dtor
  delete fHelper;
}

//_____________________________________________________________________________
AliMpVPadIterator* 
AliMpFastSegmentation::CreateIterator(const AliMpArea& area) const
{
  /// Forward to our helper
  return fHelper->CreateIterator(area);
}

//_____________________________________________________________________________
AliMpVPadIterator* 
AliMpFastSegmentation::CreateIterator() const
{
  /// Forward to our helper
  return fHelper->CreateIterator();
}

//_____________________________________________________________________________
Int_t 
AliMpFastSegmentation::GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                                     Bool_t includeSelf,
                                     Bool_t includeVoid) const
{
  /// Use default implementation
  return AliMpVSegmentation::GetNeighbours(pad,neighbours,includeSelf,includeVoid);
}

//_____________________________________________________________________________
AliMpPad 
AliMpFastSegmentation::PadByLocation(const AliMpIntPair& location, Bool_t warning) const
{
  /// Get the pad by location, using the manuid map.
  
  Int_t index = fManuId.GetValue(Encode2(location.GetFirst()));
  
  if (!index) 
  {
    if (warning)
		{
			AliWarning(Form("Manu ID %d not found",location.GetFirst()));
      Print();
    }
    return AliMpPad::Invalid();
  }
  
  AliMpMotifPosition* motifPos = InternalMotifPosition(index);
  
  if (!motifPos)
  {
    AliError(Form("InternalMotifPosition(%d) failed",index));
    Print();
    return AliMpPad::Invalid();
  }
  
  AliMpVMotif* motif = motifPos->GetMotif();
  AliMpIntPair localIndices = 
  motif->GetMotifType()->FindLocalIndicesByGassiNum(location.GetSecond());
	
  if (!localIndices.IsValid()) 
	{
		if (warning) 
		{
			AliWarning(Form("The pad number %d doesn't exists",
			                location.GetSecond()));
      Print();
		}
		return AliMpPad::Invalid();
	}
	
#ifdef CHECK
  AliMpPad pad1 = AliMpPad(location,
                           motifPos->GlobalIndices(localIndices),
                           motifPos->Position() 
                           + motif->PadPositionLocal(localIndices) 
                           - fPosition,
                           motif->GetPadDimensions(localIndices));  
  AliMpPad pad2 = fHelper->PadByLocation(location,warning);
  if ( pad1 != pad2 ) 
  {
    Print();
    pad1.Print();
    pad2.Print();
    assert(pad1==pad2);
  }
#endif
  
  return AliMpPad(location,
                  motifPos->GlobalIndices(localIndices),
                  motifPos->Position() 
                  + motif->PadPositionLocal(localIndices) 
                  - fPosition,
                  motif->GetPadDimensions(localIndices));  
}

//_____________________________________________________________________________
AliMpMotifPosition*
AliMpFastSegmentation::InternalMotifPosition(Int_t index) const
{
  /// Get the internal manu from the index
  return static_cast<AliMpMotifPosition*>(fMotifPositions.UncheckedAt(index-1));
}

//_____________________________________________________________________________
AliMpPad 
AliMpFastSegmentation::PadByIndices (const AliMpIntPair& indices, Bool_t warning) const
{
  /// Get pad by indices
  
  Int_t index = fIxIy.GetValue(Encode(indices.GetFirst(),indices.GetSecond()));
  
  if ( !index )
  {
    if (warning)
    {
      AliWarning(Form("ManuID not found for pad indices (%d,%d)",
                      indices.GetFirst(),indices.GetSecond()));	  
      Print();
    }
    return AliMpPad::Invalid();
  }
  
  AliMpMotifPosition* motifPos = InternalMotifPosition(index);

  if (!motifPos)
	{
    AliError(Form("InternalMotifPosition(%d) failed",index));
    Print();
		return AliMpPad::Invalid();
	}
	
  AliMpVMotif* motif = motifPos->GetMotif();
  AliMpMotifType* motifType = motif->GetMotifType();
  AliMpIntPair localIndices(indices-motifPos->GetLowIndicesLimit());
  AliMpConnection* connection = motifType->FindConnectionByLocalIndices(localIndices);
  
  if (!connection)
	{
		if ( warning )
		{
			AliWarning(Form("No connection for pad indices (%d,%d)",
			                indices.GetFirst(),indices.GetSecond()));
    }
    return AliMpPad::Invalid();
	}
	
#ifdef CHECK
  AliMpPad pad2 = fHelper->PadByIndices(indices,warning);
  AliMpPad pad1 = AliMpPad(AliMpIntPair(motifPos->GetID(),connection->GetManuChannel()),
                          indices,
                          motifPos->Position()
                          + motif->PadPositionLocal(localIndices)
                          - fPosition,
                          motif->GetPadDimensions(localIndices));
  
  
  
  assert(pad1==pad2);
#endif
  
  return AliMpPad(AliMpIntPair(motifPos->GetID(),connection->GetManuChannel()),
                  indices,
                  motifPos->Position()
                  + motif->PadPositionLocal(localIndices)
                  - fPosition,
                  motif->GetPadDimensions(localIndices));
  
}

//_____________________________________________________________________________
AliMpPad 
AliMpFastSegmentation::PadByPosition(const TVector2& position, Bool_t warning ) const
{
  /// Forward to our helper
  return fHelper->PadByPosition(position,warning);
}

//_____________________________________________________________________________
Int_t 
AliMpFastSegmentation::MaxPadIndexX() const
{
  /// Forward to our helper
  return fHelper->MaxPadIndexX();
}

//_____________________________________________________________________________
Int_t  
AliMpFastSegmentation::MaxPadIndexY() const
{
  /// Forward to our helper
  return fHelper->MaxPadIndexY();
}

//_____________________________________________________________________________
Int_t  
AliMpFastSegmentation::NofPads() const
{
  /// Forward to our helper
  return fHelper->NofPads();
}

//_____________________________________________________________________________
Int_t
AliMpFastSegmentation::GetNofElectronicCards() const
{
  /// Forward to our helper
  return fHelper->GetNofElectronicCards();
}

//_____________________________________________________________________________
void
AliMpFastSegmentation::GetAllElectronicCardIDs(TArrayI& ecn) const
{
  /// Forward to our helper
  fHelper->GetAllElectronicCardIDs(ecn);
}

//_____________________________________________________________________________
Bool_t 
AliMpFastSegmentation::HasPadByIndices(const AliMpIntPair& indices) const
{
  /// Whether there is a pad at the given indices
  Int_t index = fIxIy.GetValue(Encode(indices.GetFirst(),indices.GetSecond()));
  
  if ( !index ) return kFALSE;
  
  AliMpMotifPosition* mp = InternalMotifPosition(index);
  
  Bool_t r1 = mp->HasPadByIndices(indices);
#ifdef CHECK
  Bool_t r2 = fHelper->HasPadByIndices(indices);
  
  assert(r1==r2);
#endif
  return r1;
}

//_____________________________________________________________________________
Bool_t 
AliMpFastSegmentation::HasPadByLocation(const AliMpIntPair& location) const
{
  /// Whether there is a pad at the given location (de,manuid)
  
  Int_t index = fManuId.GetValue(Encode2(location.GetFirst()));
  
  if (!index) return kFALSE;
  
  AliMpMotifPosition* mp = InternalMotifPosition(index);
  
  Bool_t r1 = mp->HasPadByManuChannel(location.GetSecond());
#ifdef CHECK
  Bool_t r2 = fHelper->HasPadByLocation(location);
  
  assert(r1==r2);
#endif
  return r1;
}

//_____________________________________________________________________________
void
AliMpFastSegmentation::Print(Option_t* opt) const
{
  /// Forward to our helper
  fHelper->Print(opt);
}

//_____________________________________________________________________________
AliMp::PlaneType
AliMpFastSegmentation::PlaneType() const
{
  /// Forward to our helper
  return fHelper->PlaneType();
}

//_____________________________________________________________________________
TVector2 
AliMpFastSegmentation::Dimensions() const
{
  /// Forward to our helper
  return fHelper->Dimensions();
}

//_____________________________________________________________________________
TVector2 
AliMpFastSegmentation::Position() const
{
  /// Forward to our helper
  return fHelper->Position();
}

//_____________________________________________________________________________
Bool_t
AliMpFastSegmentation::HasMotifPosition(Int_t manuId) const
{
  /// Whether or not we have a given manu
  return ( fManuId.GetValue(Encode2(manuId)) != 0);
}

//_____________________________________________________________________________
AliMpMotifPosition*
AliMpFastSegmentation::MotifPosition(Int_t manuId) const
{
  /// Get the motifPosition object of a given manu
  Int_t index = fManuId.GetValue(Encode2(manuId));

  if (!index) 
  {
    AliMpVPadIterator* it = CreateIterator();
    it->First();
    AliMpPad pad = it->CurrentItem();
    delete it;
    AliWarning(Form("DE %04d Manu ID %04d not found",pad.GetLocation().GetFirst(),manuId));
    return 0x0;
  }

  return InternalMotifPosition(index);
}

