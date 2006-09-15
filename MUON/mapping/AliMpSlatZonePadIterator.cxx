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
// $MpId: AliMpSlatZonePadIterator.cxx,v 1.9 2006/05/24 13:58:50 ivana Exp $

#include "AliMpSlatZonePadIterator.h"

#include "AliLog.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"

#include "Riostream.h"
#include <algorithm>
#include <limits>

/// 
/// \class AliMpSlatZonePadIterator
/// 
/// Iterates over slat pads within a region of constant pad size.
/// 
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMpSlatZonePadIterator)
/// \endcond

const Double_t AliMpSlatZonePadIterator::fgkEpsilon = 1E-4; // cm
const Double_t AliMpSlatZonePadIterator::fgkDmax = std::numeric_limits<Double_t>::max();

//_____________________________________________________________________________
AliMpSlatZonePadIterator::AliMpSlatZonePadIterator(const AliMpSlat* slat,
						   const AliMpArea& area)
: AliMpVPadIterator(),
fkSlat(slat),
fSlatSegmentation(new AliMpSlatSegmentation(slat)),
fArea(area),
fOffset(0.0,0.0),
fStep(0.0,0.0),
fCurrentPad(),
fIsDone(kTRUE)
{
  //
  // Normal ctor.
  // Iteration will be done on the slat, over the crop of (area,slat_area)
  //
  if (!CropArea()) 
    {
      AliError(Form("Could not crop area : (x,y)min=(%e,%e) ; max=(%e,%e) for slat %s",
		    area.LeftBorder(),area.DownBorder(),
		    area.RightBorder(),area.UpBorder(),fkSlat->GetID()));
    }
  Invalidate();
}

//_____________________________________________________________________________
AliMpSlatZonePadIterator::~AliMpSlatZonePadIterator()
{
  //
  // Dtor.
  //
  delete fSlatSegmentation;
}

//_____________________________________________________________________________
Bool_t
AliMpSlatZonePadIterator::CropArea()
{
  //
  // Checks the area is correct, and truncate it
  // if it goes outside the slat.

  AliDebug(3,Form("Input area (%7.2f,%7.2f)->(%7.2f,%7.2f)",
		  fArea.LeftBorder(),fArea.DownBorder(),
		  fArea.RightBorder(),fArea.UpBorder()));

  // Left and right x-limits have to come from first and last pcbs
  // to deal with short and rounded pcbs cases.
  AliMpPCB* first = fkSlat->FindPCB(fArea.LeftBorder(),fArea.DownBorder());
  AliMpPCB* last = fkSlat->FindPCB(fArea.RightBorder()-fgkEpsilon,
                                 fArea.DownBorder());

  AliDebug(3,Form("First PCB %s Ixmin %2d Last PCB %s Ixmax %2d",
		  first->GetID(),first->Ixmin(),
		  last->GetID(),last->Ixmax()));

  Double_t xleft = first->ActiveXmin();
  Double_t xright = last->ActiveXmax() - fgkEpsilon;

  AliDebug(3,Form("xleft,xright=%e,%e",xleft,xright));

  Double_t xmin = std::max(fArea.LeftBorder(),xleft);
  Double_t xmax = std::min(fArea.RightBorder(),xright);
  Double_t ymin = std::max(fArea.DownBorder(),0.0);
  Double_t ymax = std::min(fArea.UpBorder(),fkSlat->DY()*2.0-fgkEpsilon);

  AliDebug(3,Form("Cropped area (%e,%e)->(%e,%e)",
		  xmin,ymin,xmax,ymax));
  
  // At this point (xmin,ymin)->(xmax,ymax) should be a zone completely included
  // inside the slat.
  // But there's so far no guarantee that it is "filling" an integer number
  // of pads. The following lines solve this, by expanding the area up to 
  // the bottomLeft and topRight limits of the pads sitting at (xmin,ymin)
  // and (xmax,ymax).

  AliMpPad bottomLeft 
    = fSlatSegmentation->PadByPosition(TVector2(xmin,ymin)-fkSlat->Position(),
                                       kFALSE);

  StdoutToAliDebug(3,
                   cout << "bottomLeft=" << endl;
                   bottomLeft.Print();
                   cout << bottomLeft.Position().X()+fkSlat->Position().X()
                   << "," << bottomLeft.Position().Y()+fkSlat->Position().Y()
                   << endl;
                   );
  
  if ( bottomLeft.IsValid() )
    {
      xmin = std::min(xmin,fkSlat->DX() + 
		      bottomLeft.Position().X() - bottomLeft.Dimensions().X());
      ymin = std::min(ymin,fkSlat->DY() + 
		      bottomLeft.Position().Y() - bottomLeft.Dimensions().Y());
    }

  AliMpPad topRight 
    = fSlatSegmentation->PadByPosition(TVector2(xmax,ymax)-fkSlat->Position(),
                                       kFALSE);
  StdoutToAliDebug(3,
                   cout << "topRight=" << endl;
                   topRight.Print();
                   cout << topRight.Position().X()+fkSlat->Position().X()
                   << "," << topRight.Position().Y()+fkSlat->Position().Y()
                   << endl;
                   
                   );

  if ( topRight.IsValid() )
    {
      xmax = std::max(xmax,fkSlat->DX() + 
		      topRight.Position().X() + topRight.Dimensions().X());
      ymax = std::max(ymax,fkSlat->DY() + 
		      topRight.Position().Y() + topRight.Dimensions().Y());
    }

  fArea = AliMpArea(TVector2((xmin+xmax)/2.0,(ymin+ymax)/2.0),
		    TVector2((xmax-xmin)/2.0,(ymax-ymin)/2.0));
 
  AliDebug(3,Form("Paddified cropped area (%7.2f,%7.2f)->(%7.2f,%7.2f)",
		  fArea.LeftBorder(),fArea.DownBorder(),
		  fArea.RightBorder(),fArea.UpBorder())); 

  // Finally set the step sizes equal to the smallest pad sizes (which is
  // hereby assumed to be that of the first pcb).
  fStep.Set(first->PadSizeX(),first->PadSizeY());

  AliDebug(3,Form("Step Sizes (%7.2f,%7.2f)",fStep.X(),fStep.Y()));

  return fArea.IsValid();
}

//_____________________________________________________________________________
AliMpPad
AliMpSlatZonePadIterator::CurrentItem() const
{
  //
  // Returns the current iteration position (i.e. a pad)
  //
  return fCurrentPad;
}

//_____________________________________________________________________________
Bool_t
AliMpSlatZonePadIterator::GetNextPosition(Double_t& x, Double_t& y)
{
  // Get the next iteration position. 
  // On input, fOffset must be a valid position (i.e. within iteration
  // area already).

  AliDebug(3,Form("input (x,y)=(%7.2f,%7.2f)",x,y));

  x += fStep.X();

  if ( x > fArea.Dimensions().X() ) 
    {
      AliDebug(3,"Going back left and one step upper");
      // Go back leftmost position...
      x = -1.0*fArea.Dimensions().X();
      // ... and up
      y += fStep.Y();
      // Update y offset
      fOffset.Set(fOffset.X(),y);
      if ( y > fArea.Dimensions().Y() )
      {
        return false;
      }
    }
  AliDebug(3,Form("output (x,y)=(%7.2f,%7.2f",x,y));
  return true;
}

//_____________________________________________________________________________
void
AliMpSlatZonePadIterator::First()
{
  //
  // (re)Starts the iteration.
  //
  
  fOffset = fArea.Dimensions()*(-1.0);
  fIsDone = kFALSE;
  SetPad(fCurrentPad,fArea.Position()+fOffset);
  if (!fCurrentPad.IsValid()) Next();
  AliDebug(3,Form("fOffset after Next=%7.2f,%7.2f",fOffset.X(),fOffset.Y()));
  if ( !fCurrentPad.IsValid() ) 
    {
      // did not find any valid pad in there, bailing out.
      fIsDone = kTRUE;
      AliError(Form("Could not initiate iterator for slat %s. "
                    " Please check the area you gave : %e,%e +- %e,%e",
                    fkSlat->GetName(),
                    fArea.Position().X(),
                    fArea.Position().Y(),
                    fArea.Dimensions().X(),
                    fArea.Dimensions().Y()));
      return;
    }
  else
    {
      // Reposition to pad center (both in x and y).
      // Please note that repositionning y is valid here, and only here
      // (i.e. do *not* do this in Next() for instance).
      fOffset.Set(fCurrentPad.Position().X()+fkSlat->DX()-fArea.Position().X(),
		  fCurrentPad.Position().Y()+fkSlat->DY()-fArea.Position().Y());
      fIsDone = kFALSE;
    }
  AliDebug(3,Form("fOffset repositionned=%7.2f,%7.2f",fOffset.X(),fOffset.Y()));
}

//_____________________________________________________________________________
void
AliMpSlatZonePadIterator::Invalidate()
{
  //
  // Invalidate the iterator.
  //
  
  fOffset = TVector2(fgkDmax,fgkDmax);
  fCurrentPad = AliMpPad::Invalid();
  fIsDone = kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMpSlatZonePadIterator::IsDone() const
{
  //
  // Whether the iteration is finished or not.
  //
  return fIsDone;
}

//_____________________________________________________________________________
void
AliMpSlatZonePadIterator::Next()
{
  // This one is the meat of the class.
  // We're iterating in x-direction mainly, starting from 
  // lower-left of the iteration area, and proceeding right,
  // until we reach right border, in which case we increment y
  // and go back to leftmost position.
  // Put otherwise, here's basically how it should work:
  // try to do x+=xstep. If outside the area, get back to xmin
  // and try y+=ystep. If outside of the window end-of-game.
  // When new x,y known, get the corresponding pad.
  // If same pad as before (should not happen if step sizes are
  // well chosen) or not valid pad (might happen for e.g. rounded pcbs), 
  // restart.
  // End of iteration occurs when both x and y are outside the iteration
  // window.
  
  if (IsDone()) return;

  AliMpPad pad(fCurrentPad);
  int n = 0;
  Double_t x(fOffset.X());
  Double_t y(fOffset.Y());

  while ( ( pad == fCurrentPad || !pad.IsValid() ) && n<100 )
  {
    ++n;
    if (GetNextPosition(x,y)==kFALSE) 
    {
      Invalidate();
      return;
    } 
    SetPad(pad,fArea.Position()+TVector2(x,y));
  }
  fCurrentPad = pad;
}

//_____________________________________________________________________________
void
AliMpSlatZonePadIterator::SetPad(AliMpPad& pad, const TVector2& pos)
{
  //
  // Sets the current pad.
  //
  pad = fSlatSegmentation->PadByPosition(pos-fkSlat->Position(),kFALSE);
  if ( pad.IsValid() )
    {
      // Reposition fOffset to pad center (only in x-direction).
      fOffset.Set(pad.Position().X()+fkSlat->DX()-fArea.Position().X(),
                  fOffset.Y());
    }
  else
  {
    AliDebug(3,Form("No pad at pos=%e,%e",pos.X(),pos.Y()));
  }
}
