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
// $MpId$

#include "AliMpTriggerSegmentation.h"

#include "AliLog.h"
#include "AliMpConnection.h"
#include "AliMpMotif.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpTrigger.h"

ClassImp(AliMpTriggerSegmentation)

//_____________________________________________________________________________
AliMpTriggerSegmentation::AliMpTriggerSegmentation() 
: AliMpVSegmentation(),
fkSlat(0)
{
  //
  // Default ctor. Not to be used really.
  //
  AliDebug(1,Form("this=%p Empty ctor",this));
}

//_____________________________________________________________________________
AliMpTriggerSegmentation::AliMpTriggerSegmentation(const AliMpTrigger* slat) 
: AliMpVSegmentation(), 
fkSlat(slat)
{
  //
  // Normal ctor.
  //
  AliDebug(1,Form("this=%p Normal ctor slat=%p",this,slat));
}

//_____________________________________________________________________________
AliMpTriggerSegmentation::~AliMpTriggerSegmentation()
{
  //
  // Dtor (empty).
  //
  AliDebug(1,Form("this=%p",this));			
}

//_____________________________________________________________________________
AliMpVPadIterator*
AliMpTriggerSegmentation::CreateIterator(const AliMpArea&) const
{
  //
  // Returns an iterator to loop over the pad contained within given area.
  // Not implemented for trigger.
  
  return 0;
}

//_____________________________________________________________________________
const char*
AliMpTriggerSegmentation::GetName() const
{
  TString name("TriggerSegmentation");
  if ( fkSlat) 
  {
    name += ".";
    name += fkSlat->GetName();
  }
  return name.Data();
}

//_____________________________________________________________________________
Bool_t
AliMpTriggerSegmentation::HasPad(const AliMpIntPair& indices) const
{
  //
  // Test if this slat has a pad located at the position referenced
  // by the integer indices.
  //
  
  return PadByIndices(indices,kFALSE) != AliMpPad::Invalid();
}

//_____________________________________________________________________________
Int_t 
AliMpTriggerSegmentation::MaxPadIndexX()
{
  //
  // Returns the value of the largest pad index in x-direction.
  //
  
  return fkSlat->GetNofPadsX()-1;
}

//_____________________________________________________________________________
Int_t 
AliMpTriggerSegmentation::MaxPadIndexY()
{
  //
  // Returns the value of the largest pad index in y-direction.
  //
  
  return fkSlat->GetMaxNofPadsY()-1;
}

//_____________________________________________________________________________
AliMpPad
AliMpTriggerSegmentation::PadByLocation(const AliMpIntPair& location, 
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
  AliMpPad pad;
  AliMpIntPair invloc;
  
  for ( Int_t i = 0; i < fkSlat->GetSize(); ++i )
  {
    const AliMpSlat* slat = fkSlat->GetLayer(i);
    AliMpSlatSegmentation seg(slat);
    AliMpPad p_i = seg.PadByLocation(location,kFALSE);
    if ( p_i.IsValid() ) 
    {
      if ( !pad.IsValid() )
      {
        pad = AliMpPad(invloc,p_i.GetIndices(),p_i.Position(),p_i.Dimensions());
        pad.AddLocation(p_i.GetLocation());
      }
      else
      {
        pad.AddLocation(p_i.GetLocation());
      }  
    }
  }
  if ( warning && !pad.IsValid()  )
  {
    AliWarning(Form("No pad found at location (%d,%d)",location.GetFirst(),
                    location.GetSecond()));
  }
  return pad;
}

//_____________________________________________________________________________
AliMpPad
AliMpTriggerSegmentation::PadByIndices(const AliMpIntPair& indices, 
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
 
  AliMpPad pad;
  AliMpIntPair invloc;
  
  for ( Int_t i = 0; i < fkSlat->GetSize(); ++i )
  {
    const AliMpSlat* slat = fkSlat->GetLayer(i);
    AliMpSlatSegmentation seg(slat);
    AliMpPad p_i = seg.PadByIndices(indices,kFALSE);
    if ( p_i.IsValid() ) 
    {      
      if ( !pad.IsValid() )
      {
        pad = AliMpPad(invloc,p_i.GetIndices(),p_i.Position(),p_i.Dimensions());
        pad.AddLocation(p_i.GetLocation());
      }
      else
      {
        pad.AddLocation(p_i.GetLocation());
      }  
    }
  }
  if ( warning && !pad.IsValid()  )
  {
    AliWarning(Form("No pad found at indices (%d,%d)",indices.GetFirst(),
                    indices.GetSecond()));
  }
  
  return pad;
}

//_____________________________________________________________________________
AliMpPad
AliMpTriggerSegmentation::PadByPosition(const TVector2& position, 
                                     Bool_t warning) const
{
  //
  // Returns the pad specified by its (floating point) position.
  // If warning=kTRUE and the pad does not exist, a warning message is 
  // printed.
  //
  // AliMpPad::Invalid() is returned if there's no pad at the given location.
  //
  AliMpPad pad;
  AliMpIntPair invloc;
  
  for ( Int_t i = 0; i < fkSlat->GetSize(); ++i )
  {
    const AliMpSlat* slat = fkSlat->GetLayer(i);
    AliMpSlatSegmentation seg(slat);
    AliMpPad p_i = seg.PadByPosition(position,kFALSE);
    if ( p_i.IsValid() ) 
    {
      if ( !pad.IsValid() )
      {
        pad = AliMpPad(invloc,p_i.GetIndices(),p_i.Position(),p_i.Dimensions());
        pad.AddLocation(p_i.GetLocation());
      }
      else
      {
        pad.AddLocation(p_i.GetLocation());
      }  
    }
  }
  if ( warning && !pad.IsValid()  )
  {
    AliWarning(Form("No pad found at position (%e,%e)",position.X(),
                    position.Y()));
  }
  
  return pad;  
}

//_____________________________________________________________________________
const AliMpTrigger* 
AliMpTriggerSegmentation::Slat() const
{
  //
  // Returns the pointer to the referenced slat.
  //
  
  return fkSlat;
}
