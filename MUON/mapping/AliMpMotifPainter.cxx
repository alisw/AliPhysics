// $Id$
// Category: graphics
//
// Class AliMpMotifPainter
// -----------------------
// Class for drawing a motif into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay

#include <TVirtualX.h>
#include <TPad.h>
#include <TError.h>
 
#include "AliMpMotifPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"

ClassImp(AliMpMotifPainter)

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter()
  : AliMpVPainter(),
    fMotifPos(0)
{
  // default dummy constructor
}

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(AliMpMotifPosition *motifPos)
  : AliMpVPainter(),
    fMotifPos(motifPos)
{
  // normal constructor 

}

//_____________________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(const AliMpMotifPainter& right) 
  : AliMpVPainter(right) {
// 
  Fatal("AliMpMotifPainter", "Copy constructor not provided.");
}

//_______________________________________________________________________
AliMpMotifPainter::~AliMpMotifPainter()
{
  // default dummy constructor
}

//_____________________________________________________________________________
AliMpMotifPainter& 
AliMpMotifPainter::operator=(const AliMpMotifPainter& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//_______________________________________________________________________
void AliMpMotifPainter::DumpObject()
{
// Draw the owned object
  fMotifPos->Dump();

}

//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetPosition() const
{
// Get the owned object's position
  return fMotifPos->Position();

}
//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetDimensions() const
{
// Get the owned object's dimensions
  return fMotifPos->Dimensions();

}

//_______________________________________________________________________
void AliMpMotifPainter::Paint(Option_t *option)
{
// Paint the object
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fMotifPos) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();

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
    {
      //PaintWholeBox(kFALSE);
      AliMpMotifType *motifType = fMotifPos->GetMotif()->GetMotifType();
      for (Int_t j=motifType->GetNofPadsY()-1;j>=0;j--){
	   for (Int_t i=0;i<motifType->GetNofPadsX();i++){
	     AliMpIntPair indices = AliMpIntPair(i,j);
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
	         gPad->PaintText(padPadPos.X()-0.01,padPadPos.Y()-0.01,
			      Form("%d",connect->GetGassiNum()));
	      
	         gVirtualX->SetTextSize(textSize);
                 }
	    }
          }
      }
    }
    break;
  default:
    PaintWholeBox(kFALSE);
  }
  gr->Pop();
  gVirtualX->SetFillColor(col);
}
