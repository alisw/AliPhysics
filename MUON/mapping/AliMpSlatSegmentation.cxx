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
// $MpId: AliMpSlatSegmentation.cxx,v 1.12 2006/05/24 13:58:50 ivana Exp $

// Caution !!
// Implementation note.
// The position(s) used in the interface are supposed to be relative
// to the slat center (AliMpSlat::Position()), whereas internally
// the x,y are relative to bottom-left corner.

#include "AliMpSlatSegmentation.h"

#include "AliLog.h"
#include "AliMpArea.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliLog.h"
#include "AliMpMotif.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpSlat.h"
#include "AliMpSlatPadIterator.h"

/// \cond CLASSIMP
ClassImp(AliMpSlatSegmentation)
/// \endcond

//_____________________________________________________________________________
AliMpSlatSegmentation::AliMpSlatSegmentation() 
: AliMpVSegmentation(),
  fkSlat(0),
  fIsOwner(false)
{
  //
  // Default ctor. Not to be used really.
  //
  AliDebug(1,Form("this=%p Empty ctor",this));
}

//_____________________________________________________________________________
AliMpSlatSegmentation::AliMpSlatSegmentation(const AliMpSlat* slat, Bool_t own) 
: AliMpVSegmentation(), 
  fkSlat(slat),
  fIsOwner(own)
{
  //
  // Normal ctor.
  //
  AliDebug(1,Form("this=%p Normal ctor slat=%p",this,slat));
}

//_____________________________________________________________________________
AliMpSlatSegmentation::~AliMpSlatSegmentation()
{
  //
  // Dtor (empty).
  //
 
  if ( fIsOwner ) delete fkSlat;
 
  // Int_t i(0);//just to be able to put a breakpoint in gdb
  AliDebug(1,Form("this=%p",this));			
}

//_____________________________________________________________________________
AliMpVPadIterator*
AliMpSlatSegmentation::CreateIterator(const AliMpArea& area) const
{
  //
  // Returns an iterator to loop over the pad contained within given area.
  //
  AliMpArea a(area.Position()+fkSlat->Position(),area.Dimensions());
  AliDebug(3,Form("Converted input area wrt to slat center : "
                  "%7.2f,%7.2f->%7.2f,%7.2f to wrt slat lower-left : "
                  "%7.2f,%7.2f->%7.2f,%7.2f ",
                  area.LeftBorder(),area.DownBorder(),
                  area.RightBorder(),area.UpBorder(),
                  a.LeftBorder(),a.DownBorder(),
                  a.RightBorder(),a.UpBorder()));
                  
  return new AliMpSlatPadIterator(fkSlat,a);
}

//_____________________________________________________________________________
AliMpVPadIterator*
AliMpSlatSegmentation::CreateIterator() const
{
  /// Returns an iterator to loop over all pads of that segmentation
  ///
  /// FIXME: we currently just forward this to the other CreateIterator,
  /// with the proper region. Might be more efficient to write a dedicated
  /// iterator ? Test that idea.
  
  AliMpArea area(TVector2(0.0,0.0),fkSlat->Dimensions());
  return CreateIterator(area);
}

//_____________________________________________________________________________
Int_t 
AliMpSlatSegmentation::GetNeighbours(const AliMpPad& pad, 
                                     TObjArray& neighbours,
                                     Bool_t includeSelf,
                                     Bool_t includeVoid) const
{
  /// Uses default implementation
  return AliMpVSegmentation::GetNeighbours(pad,neighbours,includeSelf,includeVoid);
}

//_____________________________________________________________________________
TVector2
AliMpSlatSegmentation::Dimensions() const
{
  return Slat()->Dimensions();
}

//_____________________________________________________________________________
void 
AliMpSlatSegmentation::GetAllElectronicCardIDs(TArrayI& ecn) const
{
  Slat()->GetAllMotifPositionsIDs(ecn);
}

//_____________________________________________________________________________
const char*
AliMpSlatSegmentation::GetName() const
{
  // The name of this segmentation is "SlatSegmentation"+slatName

  TString name("SlatSegmentation");
  if ( fkSlat) 
  {
    name += ".";
    name += fkSlat->GetName();
  }
  return name.Data();
}

//_____________________________________________________________________________
Bool_t
AliMpSlatSegmentation::HasPad(const AliMpIntPair& indices) const
{
  //
  // Test if this slat has a pad located at the position referenced
  // by the integer indices.
  //
  
  return PadByIndices(indices,kFALSE) != AliMpPad::Invalid();
}

//_____________________________________________________________________________
Int_t 
AliMpSlatSegmentation::MaxPadIndexX() const
{
  //
  // Returns the value of the largest pad index in x-direction.
  //
  
  return fkSlat->GetMaxPadIndexX();
}

//_____________________________________________________________________________
Int_t 
AliMpSlatSegmentation::MaxPadIndexY() const
{
  //
  // Returns the value of the largest pad index in y-direction.
  //
  
  return fkSlat->GetMaxNofPadsY()-1;
}

//_____________________________________________________________________________
Int_t 
AliMpSlatSegmentation::NofPads() const
{
/// Return number of pads defined in the slat
  
  return fkSlat->NofPads();
}

//_____________________________________________________________________________
AliMpPad
AliMpSlatSegmentation::PadByLocation(const AliMpIntPair& location, 
                                     Bool_t warning) const
{
  //
  // Returns the pad specified by its location, where location is the 
  // pair (ManuID,ManuChannel).
  // If warning=kTRUE and the pad does not exist, a warning message is 
  // printed.
  //
  // AliMpPad::Invalid() is returned if there's no pad at the given location.
  //
  Int_t manuID = location.GetFirst();
	
  AliMpMotifPosition* motifPos = fkSlat->FindMotifPosition(manuID);
	
  if (!motifPos)
	{
		if (warning)
		{
			AliWarning(Form("Manu ID %d not found in slat %s",
			                 manuID, fkSlat->GetID()));
    }
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
		}
		return AliMpPad::Invalid();
	}
	
  return AliMpPad(location,
                  motifPos->GlobalIndices(localIndices),
                  motifPos->Position() 
                  + motif->PadPositionLocal(localIndices) 
                  - fkSlat->Position(),
                  motif->GetPadDimensions(localIndices));  
}

//_____________________________________________________________________________
AliMpPad
AliMpSlatSegmentation::PadByIndices(const AliMpIntPair& indices, 
                                    Bool_t warning) const
{
  //
  // Returns the pad specified by its integer indices.
  // If warning=kTRUE and the pad does not exist, a warning message is 
  // printed.
  //
  // AliMpPad::Invalid() is returned if there's no pad at the given location.
  //
  //  
  // FIXME: except for the FindMotifPosition below, this method
  // is exactly as the one in AliMpSectorSegmentation.
  // See if we can merge them somehow.
	
  AliMpMotifPosition* motifPos = fkSlat->FindMotifPosition(indices.GetFirst(),
																									 indices.GetSecond());
  if (!motifPos)
	{
		if ( warning ) 
		{
			AliWarning(Form("No motif found containing pad location (%d,%d)",
			                 indices.GetFirst(),indices.GetSecond()));	  
		}
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
			AliWarning(Form("No connection for pad location (%d,%d)",
			                indices.GetFirst(),indices.GetSecond()));
    }
    return AliMpPad::Invalid();
	}
	
  return AliMpPad(AliMpIntPair(motifPos->GetID(),connection->GetGassiNum()),
                  indices,
                  motifPos->Position()
                  + motif->PadPositionLocal(localIndices)
                  - fkSlat->Position(),
                  motif->GetPadDimensions(localIndices));
}

//_____________________________________________________________________________
AliMpPad
AliMpSlatSegmentation::PadByPosition(const TVector2& position, 
                                     Bool_t warning) const
{
  //
  // Returns the pad specified by its (floating point) position.
  // If warning=kTRUE and the pad does not exist, a warning message is 
  // printed.
  //
  // AliMpPad::Invalid() is returned if there's no pad at the given location.
  //
  
  TVector2 blPos(position);
  
  blPos += fkSlat->Position(); // position relative to bottom-left of the slat.
  
  AliMpMotifPosition* motifPos = fkSlat->FindMotifPosition(blPos.X(),blPos.Y());
	
  if (!motifPos)
	{
		if (warning) 
		{
			AliWarning(Form("Slat %s Position (%e,%e)/center (%e,%e)/bottom-left cm "
                      " outside limits",fkSlat->GetID(),
                      position.X(),position.Y(),
                      blPos.X(),blPos.Y()));
		}
		return AliMpPad::Invalid();
	}
	
  AliMpVMotif* motif =  motifPos->GetMotif();  
  blPos -= motifPos->Position();
  AliMpIntPair localIndices 
    = motif->PadIndicesLocal(blPos);
	
  AliMpConnection* connect = 
    motif->GetMotifType()->FindConnectionByLocalIndices(localIndices);
	
  if (!connect)
	{
		if (warning) 
		{
			AliWarning(Form("Slat %s localIndices (%d,%d) outside motif %s limits",
                      fkSlat->GetID(),localIndices.GetFirst(),
                      localIndices.GetSecond(),motif->GetID().Data()));
		}
		return AliMpPad::Invalid();
	}
  
  return AliMpPad(AliMpIntPair(motifPos->GetID(),connect->GetGassiNum()),
                  motifPos->GlobalIndices(localIndices),
                  motifPos->Position()
                  + motif->PadPositionLocal(localIndices)
                  - fkSlat->Position(),
                  motif->GetPadDimensions(localIndices));  
}

//_____________________________________________________________________________
AliMp::PlaneType
AliMpSlatSegmentation::PlaneType() const
{
  return Slat()->PlaneType();
}

//_____________________________________________________________________________
void
AliMpSlatSegmentation::Print(Option_t* opt) const
{
  fkSlat->Print(opt);
}

//_____________________________________________________________________________
const AliMpSlat* 
AliMpSlatSegmentation::Slat() const
{
  //
  // Returns the pointer to the referenced slat.
  //
  
  return fkSlat;
}
