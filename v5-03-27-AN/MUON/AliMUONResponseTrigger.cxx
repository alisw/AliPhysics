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
// Class AliMUONResponseTrigger
// -------------------------------
// Implementation 
// of RPC response
//-----------------------------------------------------------------------------


#include "AliMUONResponseTrigger.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONHit.h"
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
AliMUONResponseTrigger::DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t /*timeDif*/)
{
  /// Generate 2 digits (one on each cathode) from 1 hit, i.e. no cluster-size
  /// generation (simplest response case).
  
  digits.Clear();
  
  Float_t xhit = hit.X();
  Float_t yhit = hit.Y();
  Float_t zhit = hit.Z();
  Int_t detElemId = hit.DetElemId();  
  
  Double_t x,y,z;
  Global2Local(detElemId,xhit,yhit,zhit,x,y,z);
  
  Float_t tof = hit.Age();
  Int_t twentyNano(100);
  if (tof<AliMUONConstants::TriggerTofLimit())
  {
    twentyNano=1;
  }
  
  Int_t nboard=0;

  for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
  {
    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
    
    AliMpPad pad = seg->PadByPosition(x,y,kFALSE);
    Int_t ix = pad.GetIx();
    Int_t iy = pad.GetIy();
    
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
    
    if ( cath == AliMp::kCath0 ) nboard = pad.GetLocalBoardId(0);
        
    AliMUONDigit* d = new AliMUONDigit(detElemId,nboard,
                                       pad.GetLocalBoardChannel(0),cath);
    d->SetPadXY(ix,iy);

    //FIXME : a trigger digit can have several locations. 
    //this is not currently supported by the digit class. Change that or not ?
    d->SetCharge(twentyNano);


    digits.Add(d);   
  }
}
