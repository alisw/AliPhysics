// $Id$
// Category: basic
//
// Class AliMpTransformer
// ------------------------
// Class contains definition of transformation and
// provides functions for transforming pads.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpTransformer.h"

ClassImp(AliMpTransformer)

///_____________________________________________________________________________
AliMpTransformer::AliMpTransformer(const TVector2& offset, 
                                   const AliMpIntPair& scale) 
  : TObject(),
    fOffset(offset),
    fScale(scale) 
{
//
}

///_____________________________________________________________________________
AliMpTransformer::AliMpTransformer() 
  : TObject(),
    fOffset(TVector2()),
    fScale(AliMpIntPair(1,1)) 
{
//
}

//_____________________________________________________________________________
AliMpTransformer::~AliMpTransformer() {
// 
}

//
// private methods
//

//_____________________________________________________________________________
void AliMpTransformer::CutInterval(Double_t x1, Double_t x2, Double_t x0, 
                               Double_t s, Double_t& pos, Double_t& dx) const
{
// Cuts the interval <x1,x2> into <x1', x2'>
// by position x0 in direction of s*x.
// Returns the centre of the new interval (pos) and its half size.
// ---

  // Transform values

  Double_t sx0 = s * x0;
  Double_t sx1 = s * x1;
  Double_t sx2 = s * x2;

  if (s < 0) {
    // when inversion, replace x1 and x2
    Double_t tmp = sx1;
    sx1 = sx2;
    sx2 = tmp;
  }  

  if (sx0 > sx2) {
    // the interval outside region 
    pos = 0.; 
    dx = -1.;
  }
  else if (sx0 > sx1) {
    // x0 cuts the interval
    dx = (sx2 - sx0)/2.;
    pos = s * (sx0 + dx);
  }
  else {
    // the interval inside region
    dx = (sx2 - sx1)/2.;
    pos = s * ((sx2 + sx1)/2.);
  }  
}			       

//
// public methods
//

//_____________________________________________________________________________
AliMpIntPair AliMpTransformer::Scale(const AliMpIntPair& pair) const
{
// Returns the pair with values scaled by the given scale.
// ---

  return pair * fScale;
}  

//_____________________________________________________________________________
TVector2 AliMpTransformer::Scale(const TVector2& vector) const
{
// Returns vector with values scaled by the given scale.
// ---

  return TVector2(vector.X()*fScale.GetFirst(), vector.Y()*fScale.GetSecond());
}  

//_____________________________________________________________________________
AliMpIntPair AliMpTransformer::ScaleLocation(const AliMpIntPair& orig) const
{
// Returns location with values scaled by the given scale.
// ---

  return AliMpIntPair(orig.GetFirst() * fScale.GetSecond(), orig.GetSecond(), 
                  orig.IsValid());
}  

//_____________________________________________________________________________
AliMpPad AliMpTransformer::Scale(const AliMpPad& pad) const
{
// Returns pad with indices scaled by the given scale.
// ---

  return AliMpPad(ScaleLocation(pad.GetLocation()), 
              Scale(pad.GetIndices()),
	      Scale(pad.Position()),
	      pad.Dimensions(),
	      pad.IsValid());              
}

//_____________________________________________________________________________
TVector2 AliMpTransformer::Transform(const TVector2& vector) const
{
// Transforms given vector with scale and corresponding translation
// from sector (local) to plane (global). 
// ---

  return Scale(vector) + fOffset;
}  

//_____________________________________________________________________________
TVector2 AliMpTransformer::ITransform(const TVector2& vector) const
{
// Transforms given vector with scale and corresponding translation
// from plane (global) to sector (local).
// ---

  return Scale(vector - fOffset);
}  

//_____________________________________________________________________________
AliMpPad AliMpTransformer::Transform(const AliMpPad& pad) const
{
// Returns pad with characteristics transformed with given scale and 
// corresponding translation from sector (local) to plane (global). 
// ---

  return AliMpPad(ScaleLocation(pad.GetLocation()), 
              Scale(pad.GetIndices()),
	      Transform(pad.Position()),
	      pad.Dimensions(),
	      pad.IsValid());              
}


//_____________________________________________________________________________
AliMpPad AliMpTransformer::ITransform(const AliMpPad& pad) const
{
// Returns pad with characteristics transformed with given scale and 
// corresponding translation from plane (global) to sector (local).
// ---

  return AliMpPad(ScaleLocation(pad.GetLocation()), 
              Scale(pad.GetIndices()),
	      ITransform(pad.Position()),
	      pad.Dimensions(),
	      pad.IsValid());              
}

//_____________________________________________________________________________
AliMpArea AliMpTransformer::CutArea(const AliMpArea& area) const
{
// Cuts the area with its offset in the quadrant defined by scale
// and transforms its position into its local system.
// ---
 
  Double_t posx, dx;
  CutInterval(area.LeftBorder(), area.RightBorder(), fOffset.X(), 
              fScale.GetFirst(), posx, dx);
  
  Double_t posy, dy;
  CutInterval(area.DownBorder(), area.UpBorder(), fOffset.Y(), 
              fScale.GetSecond(), posy, dy);
	      
  return AliMpArea(ITransform(TVector2(posx, posy)), TVector2(dx, dy));	        
}
