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
// $MpId: AliMpSlatPainter.cxx,v 1.10 2006/05/24 13:58:32 ivana Exp $

///
/// \class AliMpSlatPainter
/// 
/// Class for drawing a slat into canvas
///
/// \author Laurent Aphecetche

#include "AliMpSlatPainter.h"

#include "AliLog.h"
#include "AliMpGraphContext.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"

#include "TVirtualX.h"
#include "TPad.h"
#include <iostream>

/// \cond CLASSIMP
ClassImp(AliMpSlatPainter)
/// \endcond

//_____________________________________________________________________________
AliMpSlatPainter::AliMpSlatPainter()
 : AliMpVPainter(),
   fkSlat(0)
{
  //
  // Empty ctor.
  //
}

//_____________________________________________________________________________
AliMpSlatPainter::AliMpSlatPainter(const AliMpSlat* slat)
 : AliMpVPainter(),
   fkSlat(slat)
{
    //
    // Normal ctor.
    //
}

//_____________________________________________________________________________
AliMpSlatPainter::~AliMpSlatPainter()
{
  //
  // Dtor.
  //
}

//_____________________________________________________________________________
TVector2
AliMpSlatPainter::GetDimensions() const
{
  //
  // Returns the half-sizes of the slat.
  //
  return TVector2(fkSlat->DX(),fkSlat->DY());
}

//_____________________________________________________________________________
TVector2
AliMpSlatPainter::GetPosition() const
{
  //
  // Returns the (x,y) position of the slat.
  //
  return TVector2(fkSlat->DX(),fkSlat->DY());
}

//_____________________________________________________________________________
void
AliMpSlatPainter::Draw(Option_t* option)
{
  //
  // Draws the slat.
  //
  // If option[0] is 'P' then PCB are drawn too.
  //
  AliMpGraphContext *gr = AliMpGraphContext::Instance();

  gr->Push();
  InitGraphContext();

//   GetPosition().Print();
//   GetDimensions().Print();

  switch (option[0])
    {
    case 'P':
      for ( AliMpSlat::Size_t i = 0; i < fkSlat->GetSize(); ++i )
	{
	  AliMpPCB* pcb = fkSlat->GetPCB(i);
	  
	  gr->Push();

	  gr->SetPadPosForReal(TVector2(pcb->X(),pcb->Y()),
			       TVector2(pcb->DX(),pcb->DY()));
	  gr->SetColor(i+2);
	  
	  DrawObject(pcb,option+1);
	  
	  gr->Pop();
	}
      break;
    default:
      AppendPad(option);
    }

  gr->Pop();
}

//_____________________________________________________________________________
void
AliMpSlatPainter::Paint(Option_t* /*option*/)
{
  //
  // Paint the object.
  //
  AliMpGraphContext* gr = AliMpGraphContext::Instance();

  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();

  PaintWholeBox(kTRUE);
  
  gr->Pop();
  gVirtualX->SetFillColor(col);
}


