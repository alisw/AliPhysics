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
// $MpId: AliMpPCBPainter.cxx,v 1.4 2005/09/19 18:56:46 ivana Exp $
// Category: graphics
//
// Class AliMpSlatPainter
// -----------------------
// Class for drawing a PCB into canvas
// Authors: Laurent Aphecetche

#include "AliMpPCBPainter.h"

#include "AliMpGraphContext.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"

#include "AliLog.h"

#include "TVirtualX.h"
#include "TPad.h"

#include <iostream>

ClassImp(AliMpPCBPainter)

//_____________________________________________________________________________
AliMpPCBPainter::AliMpPCBPainter(AliMpPCB* pcb)
  : AliMpVPainter(), fPCB(pcb)
{
    //
    // Default ctor.
    //
}

//______________________________________________________________________________
AliMpPCBPainter::AliMpPCBPainter(const AliMpPCBPainter& right) 
  : AliMpVPainter(right) 
{  
/// Protected copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpPCBPainter::~AliMpPCBPainter()
{
  //
  // Dtor.
  //
}

//______________________________________________________________________________
AliMpPCBPainter& 
AliMpPCBPainter::operator=(const AliMpPCBPainter& right)
{
/// Protected assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//_____________________________________________________________________________
TVector2
AliMpPCBPainter::GetDimensions() const
{
  //
  // Returns the half-sizes of the PCB.
  //
  return TVector2(fPCB->DX(),fPCB->DY());
}

//_____________________________________________________________________________
TVector2
AliMpPCBPainter::GetPosition() const
{
  //
  // Returns the (x,y) position of the PCB.
  //
  return TVector2(fPCB->X(),fPCB->Y());
}

//_____________________________________________________________________________
void
AliMpPCBPainter::Draw(Option_t* option)
{
  //
  // Draws the PCB.
  //
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fPCB) return;

  gr->Push();
  InitGraphContext();

  switch (option[0])
  {
    case 'M':
      for ( AliMpPCB::Size_t i = 0; i < fPCB->GetSize(); ++i )
      {
        AliMpMotifPosition* pos = fPCB->GetMotifPosition(i);
        
        gr->Push();
        gr->SetPadPosForReal(pos->Position(),
                             pos->Dimensions());
        gr->SetColor(gr->GetColor()+i);
        
        DrawObject(pos,option+1);
        
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
AliMpPCBPainter::Paint(Option_t* /*option*/)
{
  //
  // Paint the object.
  //
  AliMpGraphContext* gr = AliMpGraphContext::Instance();
  if (!fPCB) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();

  PaintWholeBox(kTRUE);
  
  gr->Pop();
  gVirtualX->SetFillColor(col);
}
