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

#include "AliMpIteratorPainter.h"

#include "AliMpVPadIterator.h"
#include "AliMpPad.h"
#include "AliMpGraphContext.h"

#include "AliLog.h"

#include "TObjArray.h"
#include "TVirtualX.h"
#include "TVirtualPad.h"
#include "TVector2.h"
#include "TMath.h"

//-----------------------------------------------------------------------------
/// \class AliMpIteratorPainter
///
/// A painter for a group of pads, which is defined by an iterator
///
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMpIteratorPainter)
/// \endcond

//_____________________________________________________________________________
AliMpIteratorPainter::AliMpIteratorPainter(AliMpVPadIterator *it)
: AliMpVPainter(), fPads(new TObjArray), fPosition(), fDimensions()
{
  /// Ctor. Iterator must not be null.
  
  //if (!it) throw;
  if (!it) {
    AliFatal("Iterator must not be null.");;
  }  
 
  Double_t xmin(1E9), xmax(-1E9);
  Double_t ymin(1E9), ymax(-1E9);
  
  fPads->SetOwner(kTRUE);
  it->First();
  while ( !it->IsDone() ) 
  {
    AliMpPad pad = it->CurrentItem();
    fPads->AddLast(new AliMpPad(pad));
    TVector2 lowerLeft(pad.Position()-pad.Dimensions());
    TVector2 upperRight(pad.Position()+pad.Dimensions());
    xmin = TMath::Min(lowerLeft.X(),xmin);
    ymin = TMath::Min(lowerLeft.Y(),ymin);
    xmax = TMath::Max(upperRight.X(),xmax);
    ymax = TMath::Max(upperRight.Y(),ymax);
    it->Next();
  }
  fPosition = TVector2((xmin+xmax)/2.0,(ymin+ymax)/2.0);
  fDimensions = TVector2((xmax-xmin)/2.0,(ymax-ymin)/2.0);
}

//_____________________________________________________________________________
AliMpIteratorPainter::~AliMpIteratorPainter()
{
  /// dtor
  delete fPads;
}

//_____________________________________________________________________________
void
AliMpIteratorPainter::Draw(Option_t* option)
{
  /// Append ourselves to the current graphic pad
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  gr->Push();
  InitGraphContext();
  gr->SetPadPosForReal(GetPosition(),GetDimensions());
  AppendPad(option);
  gr->Pop();
}

//_____________________________________________________________________________
void
AliMpIteratorPainter::Paint(Option_t*)
{
  /// Actual drawing method
  
  AliMpGraphContext* gr = AliMpGraphContext::Instance();
  gr->Push();
  InitGraphContext();
  gPad->Range(0.,0.,1.,1.);
  
  TIter next(fPads);
  AliMpPad* pad;
  
  while ( ( pad = static_cast<AliMpPad*>(next()) ) )
  {
    TVector2 padPadPos;
    TVector2 padPadDim;
    
    gr->RealToPad(pad->Position(),
                  pad->Dimensions(),
                  padPadPos,
                  padPadDim);

    TVector2 bl = padPadPos - padPadDim;
    TVector2 ur = padPadPos + padPadDim;
    
    Int_t manuId = pad->GetManuId();

    Style_t sty = gVirtualX->GetFillStyle();

    gVirtualX->SetFillStyle(1);
    if (manuId % 5 == 0) 
	gVirtualX->SetFillColor(0);
    if (manuId % 5 == 1) 
	gVirtualX->SetFillColor(38);
    if (manuId % 5 == 2) 
	gVirtualX->SetFillColor(33);
    if (manuId % 5 == 3) 
	gVirtualX->SetFillColor(16);
    if (manuId % 5 == 4) 
	gVirtualX->SetFillColor(44);

    gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
    gVirtualX->SetFillStyle(0);
    gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
    gVirtualX->SetFillStyle(sty);    

    Float_t textSize =   gVirtualX->GetTextSize();
    gVirtualX->SetTextSize(10);
    gVirtualX->SetTextAlign(22);
    gPad->PaintText((bl.X()+ur.X())/2.0,(bl.Y()+ur.Y())/2.0,
                    Form("%d",pad->GetManuChannel()));
    
    gVirtualX->SetTextSize(textSize);
    
  }
}


