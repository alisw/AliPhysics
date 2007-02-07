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
// $MpId: AliMpMotifPainter.cxx,v 1.9 2006/05/24 13:58:32 ivana Exp $
// Category: graphics
//
// Class AliMpMotifPainter
// -----------------------
// Class for drawing a motif into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay

#include "AliMpMotifPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpMotif.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"
#include "AliLog.h"

#include <TVirtualX.h>
#include <TPad.h>
 
/// \cond CLASSIMP
ClassImp(AliMpMotifPainter)
/// \endcond

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter()
  : AliMpVPainter(),
    fMotifPos(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(AliMpMotifPosition *motifPos)
  : AliMpVPainter(),
    fMotifPos(motifPos)
{
  /// Standard constructor 
      AliDebug(1,"Default ctor");
}

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(AliMpMotifType* motifType)
: AliMpVPainter(),
fMotifPos(0x0)
{
  /// Constructor from a motif Type. We hereby create a MotifPosition
  /// object from it, using arbitrary pad sizes, as this is just a way
  /// to visualize the *shape* of the motif.
  
  AliDebug(1,"Ctor from motifType");
  
  const Double_t kdx = 5;
  const Double_t kdy = 5; // cm but arbitrary anyway
  
  AliMpVMotif* motif = new AliMpMotif(motifType->GetID(),
                                      motifType,
                                      TVector2(kdx,kdy));

  fMotifPos = new AliMpMotifPosition(-1,motif,motif->Dimensions());
}

//_______________________________________________________________________
AliMpMotifPainter::~AliMpMotifPainter()
{
  /// Default constructor
}

//_______________________________________________________________________
void AliMpMotifPainter::DumpObject()
{
/// Dump the owned object

  fMotifPos->Dump();
}

//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetPosition() const
{
/// Get the owned object's position

  return fMotifPos->Position();
}

//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetDimensions() const
{
/// Get the owned object's dimensions

  return fMotifPos->Dimensions();
}

//_______________________________________________________________________
void AliMpMotifPainter::Paint(Option_t *option)
{
/// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fMotifPos) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();

  gVirtualX->SetLineWidth(1);
  
  switch (option[0]){
    case 'T':
    case 'I':
    case 'X':
    {
      PaintWholeBox();
      Float_t textSize =   gVirtualX->GetTextSize();
      gVirtualX->SetTextSize(10);
      TString str;
      switch (option[0]) {
        case 'T' : 
          str = Form("%d",fMotifPos->GetID());
          break;
        case 'I':{
          switch (option[1]){
            case '+' :
              str = Form("(%d,%d)",fMotifPos->GetHighIndicesLimit().GetFirst(),
                         fMotifPos->GetHighIndicesLimit().GetSecond());
              break;
            default:
              str = Form("(%d,%d)",fMotifPos->GetLowIndicesLimit().GetFirst(),
                         fMotifPos->GetLowIndicesLimit().GetSecond());
      	  }
        }
          break;
        case 'X' :
          str = Form("(%f,%f)",(GetPosition()-GetDimensions()).X(),
                     (GetPosition()-GetDimensions()).Y());
          break;
      }
      gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,str);
      
      gVirtualX->SetTextSize(textSize);
    }
      break;
    case 'P':
    case 'Z':
    {
      //PaintWholeBox(kFALSE);
      AliMpMotifType *motifType = fMotifPos->GetMotif()->GetMotifType();
      for (Int_t j=motifType->GetNofPadsY()-1;j>=0;j--){
        for (Int_t i=0;i<motifType->GetNofPadsX();i++){
          AliMpIntPair indices(i,j);
          AliMpConnection* connect = 
            motifType->FindConnectionByLocalIndices(indices);
          if (connect){
            TVector2 realPadPos = 
            GetPosition()+fMotifPos->GetMotif()->PadPositionLocal(indices);
            TVector2 padPadPos,padPadDim;
            gr->RealToPad(realPadPos,
                          fMotifPos->GetMotif()->GetPadDimensions(indices),
                          padPadPos,padPadDim);
            TVector2 bl = padPadPos - padPadDim;
            TVector2 ur = padPadPos + padPadDim;

            Style_t sty = gVirtualX->GetFillStyle();
            gVirtualX->SetFillStyle(1);
            gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
            gVirtualX->SetFillStyle(0);
            gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
            gVirtualX->SetFillStyle(sty);

            if (option[1]=='T'){
              Float_t textSize =   gVirtualX->GetTextSize();
              gVirtualX->SetTextSize(10);
              gVirtualX->SetTextAlign(22);
              gPad->PaintText((bl.X()+ur.X())/2.0,(bl.Y()+ur.Y())/2.0,
                              Form("%d",connect->GetGassiNum()));
              
              gVirtualX->SetTextSize(textSize);
                 }
	    }
          }
      }
      if ( option[0]=='Z' )
      {
        PaintContour(option,kFALSE);
      }
    }
      break;
      
    case 'C':
      PaintContour(option,kTRUE);
      break;
      
    default:
      PaintWholeBox(kFALSE);
  }
  gr->Pop();
  gVirtualX->SetFillColor(col);
}

//_______________________________________________________________________
void AliMpMotifPainter::PaintContour(Option_t* option, Bool_t fill)
{
/// Drawing real motif (not envelop) the real contour

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  
  Float_t xl = 0;
  Float_t yl = 0;
  Int_t manuId = 0;
  Int_t searchMotif = -1;
  TVector2 bl0 = TVector2(999, 999);
  TVector2 ur0 = TVector2(0,0); 
  TVector2 padPadPos,padPadDim;
  
  AliMpMotifType *motifType = fMotifPos->GetMotif()->GetMotifType();
  manuId = fMotifPos->GetID();
  
  if ( fill )
  {
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
  }
  
  Width_t lineW = gPad->GetLineWidth();
  Width_t lw = lineW*3;
  Double_t xlw = gPad->PixeltoX(lw/2);
  

  if (option[1] == 'I' && option[2] == ':')
      searchMotif = atoi(&option[3]);

  gVirtualX->SetLineWidth(lw);
    
    for (Int_t i = 0; i < motifType->GetNofPadsX(); i++){
      
      for (Int_t j = 0; j < motifType->GetNofPadsY(); j++){
        
        AliMpIntPair indices = AliMpIntPair(i,j);
        AliMpConnection* connect =  motifType->FindConnectionByLocalIndices(indices);

        if (connect){
          TVector2 realPadPos = 
          GetPosition()+fMotifPos->GetMotif()->PadPositionLocal(indices);
          gr->RealToPad(realPadPos, fMotifPos->GetMotif()->GetPadDimensions(indices),
                        padPadPos, padPadDim);
          
          TVector2 bl = padPadPos - padPadDim;
          TVector2 ur = padPadPos + padPadDim;

          if (bl0.X() > bl.X())
            bl0 = bl;
          
          if (ur0.Y() < ur.Y())
            ur0 = ur;
          
          if ( fill )
          {
	    Style_t csty = gVirtualX->GetFillColor();
	    Style_t sty = gVirtualX->GetFillStyle();
	    gVirtualX->SetFillStyle(1);
	    if (manuId == searchMotif) 
		gVirtualX->SetFillColor(5); // yellow
	    gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
            gVirtualX->SetFillStyle(sty);
            gVirtualX->SetFillColor(csty);
          } 

	  if (!motifType->FindConnectionByLocalIndices(AliMpIntPair(i,j-1)))
          {
            gPad->PaintLine(bl.X()-xlw, bl.Y(), bl.X()+ padPadDim.X()*2 + xlw, bl.Y());
          }
          
          if (!motifType->FindConnectionByLocalIndices(AliMpIntPair(i,j+1)))
          {
            gPad->PaintLine(bl.X()-xlw, bl.Y() + padPadDim.Y()*2, bl.X()+ padPadDim.X()*2+xlw, bl.Y() +  padPadDim.Y()*2);
          }
          if (!motifType->FindConnectionByLocalIndices(AliMpIntPair(i-1,j)))
          {
            gPad->PaintLine(bl.X(), bl.Y(), bl.X(), bl.Y()+ padPadDim.Y()*2);                  
          }          

          if (!motifType->FindConnectionByLocalIndices(AliMpIntPair(i+1,j)))
          {
            gPad->PaintLine(bl.X()+padPadDim.X()*2, bl.Y(), bl.X()+padPadDim.X()*2, bl.Y()+ padPadDim.Y()*2);  
          } 
        }

      }
    }
    
    switch (option[1]) {
      // add manudId indexes
      case 'I' :
        xl = bl0.X()+ padPadDim.X()/2.;
        
        yl = bl0.Y() + 1.5*padPadDim.Y();
        
        Float_t textSize =   gVirtualX->GetTextSize();
        gVirtualX->SetTextSize(12);
        gVirtualX->SetTextAlign(13);
        gVirtualX->SetTextAngle(90.);
        
        gPad->PaintText(xl, yl, Form("%d", manuId));
        
        gVirtualX->SetTextAngle(0.);
        gVirtualX->SetTextSize(textSize);
        break;
    }

    gVirtualX->SetLineWidth(lineW);

}
