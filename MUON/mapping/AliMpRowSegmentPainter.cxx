// $Id$
// Category: graphics
//
// Class AliMpRowSegmentPainter
// ----------------------------
// Class for drawing a motif into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
 
#include <TVirtualX.h>
#include <TPad.h>
#include <TError.h>

#include "AliMpRowSegmentPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpVRowSegment.h"
#include "AliMpRow.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"

ClassImp(AliMpRowSegmentPainter)

//_______________________________________________________________________
AliMpRowSegmentPainter::AliMpRowSegmentPainter()
  : AliMpVPainter(),
    fRowSegment(0)
{
  // default dummy constructor
}

//_______________________________________________________________________
AliMpRowSegmentPainter::AliMpRowSegmentPainter(AliMpVRowSegment *row)
  : AliMpVPainter(),
    fRowSegment(row)
{
  // normal constructor 

}

//_____________________________________________________________________________
AliMpRowSegmentPainter::AliMpRowSegmentPainter(
                                       const AliMpRowSegmentPainter& right) 
  : AliMpVPainter(right) 
{  
  // copy constructor (not implemented)

  Fatal("AliMpRowSegmentPainter", "Copy constructor not provided.");
}

//_______________________________________________________________________
AliMpRowSegmentPainter::~AliMpRowSegmentPainter()
{
  // destructor
}

//_____________________________________________________________________________
AliMpRowSegmentPainter& 
AliMpRowSegmentPainter::operator=(const AliMpRowSegmentPainter& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//_______________________________________________________________________
TVector2 AliMpRowSegmentPainter::GetPosition() const
{
// Get the owned object's position
  return fRowSegment->Position();

}
//_______________________________________________________________________
TVector2 AliMpRowSegmentPainter::GetDimensions() const
{
// Get the owned object's dimensions
  return fRowSegment->Dimensions();

}

//_______________________________________________________________________
void AliMpRowSegmentPainter::DumpObject()
{
// Draw the owned object
  fRowSegment->Dump();

}

//_______________________________________________________________________
void AliMpRowSegmentPainter::Draw(Option_t *option)
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
  if( !fRowSegment) return;

  gr->Push();
  InitGraphContext();
  switch (option[0]){
  case 'M':
    {

      for (Int_t iMotif=0;iMotif<fRowSegment->GetNofMotifs();++iMotif){
	gr->Push();

	Int_t motifPositionId = fRowSegment->GetMotifPositionId(iMotif);
	AliMpMotifPosition *motifPos = 
	    fRowSegment->GetRow()->GetMotifMap()
	      ->FindMotifPosition(motifPositionId);
				     
	gr->SetPadPosForReal(motifPos->Position(),motifPos->Dimensions());
      	gr->SetColor(GetColor());
	DrawObject(motifPos,option+1);

	gr->Pop();
      }
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}


//_______________________________________________________________________
void AliMpRowSegmentPainter::Paint(Option_t* /*option*/)
{
// Paint the object
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fRowSegment) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();
  PaintWholeBox(kTRUE);

//   Float_t textSize =   gVirtualX->GetTextSize();
//   if (option[0]=='T')
//     gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,
// 		     Form("%d",fRowSegment->GetMotif()->GetID()));
//   gVirtualX->SetTextSize(textSize);
  gr->Pop();
  gVirtualX->SetFillColor(col);
}
