// $Id$
// Category: graphics
//
// Class AliMpRowPainter
// ---------------------
// Class for drawing a row into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
  
#include <TVirtualX.h>
#include <TPad.h>
#include <TError.h>
 
#include "AliMpRowPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpRow.h"
#include "AliMpRowSegment.h"

ClassImp(AliMpRowPainter)

//_______________________________________________________________________
AliMpRowPainter::AliMpRowPainter()
  : AliMpVPainter(),
    fRow(0)
{
  // default dummy constructor
}

//_______________________________________________________________________
AliMpRowPainter::AliMpRowPainter(AliMpRow *row)
  : AliMpVPainter(),
    fRow(row)
{
  // normal constructor 
}

//_____________________________________________________________________________
AliMpRowPainter::AliMpRowPainter(const AliMpRowPainter& right) 
  : AliMpVPainter(right) {
// 
  // copy constructor (not implemented)

  Fatal("AliMpRowPainter", "Copy constructor not provided.");
}

//_______________________________________________________________________
AliMpRowPainter::~AliMpRowPainter()
{
  // destructor
}

//_____________________________________________________________________________
AliMpRowPainter& AliMpRowPainter::operator=(const AliMpRowPainter& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//_______________________________________________________________________
void AliMpRowPainter::DumpObject()
{
// Draw the owned object
  fRow->Dump();

}

//_______________________________________________________________________
TVector2 AliMpRowPainter::GetPosition() const
{
// Get the owned object's position
  return fRow->Position();

}

//_______________________________________________________________________
TVector2 AliMpRowPainter::GetDimensions() const
{
// Get the owned object's dimensions
  return fRow->Dimensions();

}

//_______________________________________________________________________
void AliMpRowPainter::Draw(Option_t *option)
{
// Draw the sector on the current pad
// The first letter of <option> is treated as follows:
// case "S" : each row segments are drawn separately
// case ""  : the whole row is drawn at once
// in both cases, the rest of the option is passed
// as argument to the Draw function of respectively
// zone or row objects.
// ---
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if( !fRow) return;

  gr->Push();
  InitGraphContext();

  switch (option[0]){
  case 'S':
    {
      for (Int_t iRowSeg=0;iRowSeg<fRow->GetNofRowSegments();++iRowSeg){
	AliMpVRowSegment *rowSegment = fRow->GetRowSegment(iRowSeg);
	gr->Push();

	gr->SetPadPosForReal(rowSegment->Position(),rowSegment->Dimensions());
	DrawObject(rowSegment,option+1);
      
	gr->Pop();
      }
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}

//_______________________________________________________________________
void AliMpRowPainter::Paint(Option_t *option)
{
  // Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if( !fRow) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();  
  PaintWholeBox(kTRUE);

  if (option[0]=='T'){
    Float_t textSize =   gVirtualX->GetTextSize();
    gVirtualX->SetTextSize(12);
    gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,
		    Form("%d",fRow->GetID()));
    gVirtualX->SetTextSize(textSize);
  }


  gr->Pop();
  gVirtualX->SetFillColor(col);
}
