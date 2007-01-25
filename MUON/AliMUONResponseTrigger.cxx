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

// -------------------------------
// Class AliMUONResponseTrigger
// -------------------------------
// Implementation 
// of RPC response


#include "AliMUONResponseTrigger.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONHit.h"
#include "AliMUONSegmentation.h"
#include "AliMUONTriggerSegmentation.h"
#include "AliMUONConstants.h"

#include "AliMpPad.h"
#include "AliMpCathodType.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"

#include "AliRun.h"
#include "AliLog.h"
#include "TList.h"

/// \cond CLASSIMP
ClassImp(AliMUONResponseTrigger)
/// \endcond

namespace
{
  AliMUON* muon()
  {
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  }

  void Global2Local(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg,
                  Double_t& xl, Double_t& yl, Double_t& zl)
  {  
  // ideally should be : 
  // Double_t x,y,z;
  // AliMUONGeometry::Global2Local(detElemId,xg,yg,zg,x,y,z);
  // but while waiting for this geometry singleton, let's go through
  // AliMUON still.
  
    const AliMUONGeometryTransformer* transformer = muon()->GetGeometryTransformer();
    transformer->Global2Local(detElemId,xg,yg,zg,xl,yl,zl);
  }

  AliMUONSegmentation* Segmentation()
  {
    static AliMUONSegmentation* segmentation = muon()->GetSegmentation();
    return segmentation;
  }
}

//------------------------------------------------------------------   
AliMUONResponseTrigger::AliMUONResponseTrigger()
  : AliMUONResponse()
{
/// Default constructor
}

//------------------------------------------------------------------   
AliMUONResponseTrigger::~AliMUONResponseTrigger()
{
/// Destructor
}

//_____________________________________________________________________________
void 
AliMUONResponseTrigger::DisIntegrate(const AliMUONHit& hit, TList& digits)
{
  /// Generate 2 digits (one on each cathode) from 1 hit, i.e. no cluster-size
  /// generation (simplest response case).
  
  digits.Clear();
  
  Float_t xhit = hit.X();
  Float_t yhit = hit.Y();
  Float_t zhit = 0; // FIXME : should it be hit.Z() ?
  Int_t detElemId = hit.DetElemId();  
  
  Double_t x,y,z;
  Global2Local(detElemId,xhit,yhit,zhit,x,y,z);
  
  Float_t tof = hit.Age();
  Int_t twentyNano(100);
  if (tof<AliMUONConstants::TriggerTofLimit())
  {
    twentyNano=1;
  }
  
  for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
  {
    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
    
    AliMpPad pad = seg->PadByPosition(TVector2(x,y),kFALSE);
    Int_t ix = pad.GetIndices().GetFirst();
    Int_t iy = pad.GetIndices().GetSecond();
    
    AliDebug(1,Form("xhit,yhit=%e,%e lx,ly,lz=%e,%e,%e ix,iy=%d,%d",
                    xhit,yhit,x,y,z,ix,iy));
    
    if ( !pad.IsValid() )
    {
      AliWarning(Form("hit w/o strip %d-%d xhit,yhit=%e,%e local x,y,z "
                      "%e,%e,%e ix,iy=%d,%d",detElemId,
                      cath,
                      xhit,yhit,x,y,z,ix,iy));
      continue;
    }
    AliMUONDigit* d = new AliMUONDigit;
    d->SetDetElemId(detElemId);
/* pc 09/02/06 no need for that anymore : trigger is in local numbering

    //FIXME: >> the following code to get the ixGlo and iyGlo is a bad hack 
    // because trigger has not yet switched to local numbering of its indices !
    // We should be able to use directly the (local) ix,iy from the pad !
    const AliMUONTriggerSegmentationV2* old = 
      dynamic_cast<const AliMUONTriggerSegmentationV2*>
        (Segmentation()->GetDESegmentation(detElemId,cath));
    if ( !old )
    {
      AliFatal("Got a wrong TriggerSegmentation object! Check that!");
    }
    Int_t ixGlo;
    Int_t iyGlo;
    old->ILoc2IGlo(ix,iy,ixGlo,iyGlo);
    if ( xhit < 0 ) ixGlo = -ixGlo;
    // << end of bad hack.
    d->SetPadX(ixGlo);
    d->SetPadY(iyGlo);
*/
    d->SetPadX(ix);
    d->SetPadY(iy);

    d->SetSignal(twentyNano);
    d->AddPhysicsSignal(d->Signal());
    d->SetCathode(cath);
    digits.Add(d);   
 //   AliDebug(1,Form("Adding digit DE %d Cathode %d (%d,%d) signal %d",
//                    detElemId,cath,ixGlo,iyGlo,twentyNano));
  }
  
//  StdoutToAliDebug(1,digits.Print();); 
//  AliDebug(1,Form("Number of digits for detelem %d track %d : %d",
//                  hit.DetElemId(),hit.Track(),digits.GetSize()));
//   
}





