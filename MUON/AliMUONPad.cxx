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

#include "AliMUONPad.h"

#include "AliLog.h"
#include "AliMpArea.h"

#include "Riostream.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TVector2.h"
#include "TMath.h"

/// \class AliMUONPad
///
/// Object gathering information about a hit pad.
/// 
/// Can be seen as a combination of a Digit (which brings the charge) 
/// and an MpPad (which brings location and position)
///
/// Also provided are some static functions to compute overlap and
/// get neighboring information.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONPad)
/// \endcond

namespace
{
  AliMpArea
  Intersect(const AliMpArea& a, const AliMpArea& b)
  { 
    //
    // Returns the common part of a and b.
    //
    Double_t xmin = TMath::Max(a.LeftBorder(),b.LeftBorder());
    Double_t xmax = TMath::Min(a.RightBorder(),b.RightBorder());
    Double_t ymin = TMath::Max(a.DownBorder(),b.DownBorder());
    Double_t ymax = TMath::Min(a.UpBorder(),b.UpBorder());
    AliMpArea c( TVector2( (xmin+xmax)/2.0, (ymin+ymax)/2.0 ),
                 TVector2( (xmax-xmin)/2.0, (ymax-ymin)/2.0 ) );
	
    return c;
  }
}

//_____________________________________________________________________________
AliMUONPad::AliMUONPad()
:
TObject(),
fIsSaturated(kFALSE),
fIsReal(kFALSE),
fClusterId(-1),
fCathode(-1),
fDetElemId(-1),
fDigitIndex(-1),
fIx(-1),
fIy(-1),
fStatus(0),
fDimensions(),
fPosition(),
fCharge(0.0),
fChargeBackup(0.0)
{
  /// Default ctor
  Init(-1,-1,-1,-1,TVector2(0,0),TVector2(0,0),0);
}

//_____________________________________________________________________________
AliMUONPad::AliMUONPad(Int_t detElemId, Int_t cathode,
                       Int_t ix, Int_t iy, Double_t x, Double_t y,
                       Double_t dx, Double_t dy, Double_t charge)
:
TObject(),
fIsSaturated(kFALSE),
fIsReal(kFALSE),
fClusterId(-1),
fCathode(-1),
fDetElemId(-1),
fDigitIndex(-1),
fIx(-1),
fIy(-1),
fStatus(0),
fDimensions(),
fPosition(),
fCharge(0.0),
fChargeBackup(0.0)

{
  /// Normal ctor, using full information
  Init(detElemId,cathode,ix,iy,TVector2(x,y),TVector2(dx,dy),charge);
}

//_____________________________________________________________________________
AliMUONPad::AliMUONPad(Double_t x, Double_t y,
                       Double_t dx, Double_t dy, Double_t charge)
: TObject(),
fIsSaturated(kFALSE),
fIsReal(kFALSE),
fClusterId(-1),
fCathode(-1),
fDetElemId(-1),
fDigitIndex(-1),
fIx(-1),
fIy(-1),
fStatus(0),
fDimensions(),
fPosition(),
fCharge(0.0),
fChargeBackup(0.0)
{
  /// Truncated constructor (w/o DE, cath, ix, iy)
  Init(-1,-1,-1,-1,TVector2(x,y),TVector2(dx,dy),charge);
}

//_____________________________________________________________________________
AliMUONPad::AliMUONPad(const TVector2& position, const TVector2& dimensions,
                       Double_t charge)
: TObject(),
fIsSaturated(kFALSE),
fIsReal(kFALSE),
fClusterId(-1),
fCathode(-1),
fDetElemId(-1),
fDigitIndex(-1),
fIx(-1),
fIy(-1),
fStatus(0),
fDimensions(),
fPosition(),
fCharge(0.0),
fChargeBackup(0.0)
{
  /// Alternate ctor
  Init(-1,-1,-1,-1,position,dimensions,charge);
}

//_____________________________________________________________________________
Bool_t 
AliMUONPad::AreNeighbours(const AliMUONPad& d1, const AliMUONPad& d2) 
{
  /// Whether 2 pads are neighbours or not
  if ( d1.DetElemId() != d2.DetElemId() || 
       d1.Cathode() != d2.Cathode() )
  {
    return kFALSE;
  }
  else
  {
    static Double_t precision = 1E-4; // cm
    static TVector2 precisionAdjustment(-precision,-precision);    
    return AreOverlapping(d1,d2,precisionAdjustment);
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONPad::AreOverlapping(const AliMUONPad& d1, const AliMUONPad& d2,
                           const TVector2& precision)
{
  /// Checks the overlap between 2 pads.
  /// The actual overlap is computed not on d1 and d2, but on d1 and d2
  /// "re-scaled" using the precision vector (conceptually equivalent to 
  /// d.Dimensions() += precision)
  ///
  /// So, if the elements (x,y) of precision are :
  ///
  /// - positive, the overlap is "computed" from smaller d1 and d2  
  ///   which is, modulo the precision, what you would think as a normal
  ///   overlap calculation
  /// - negative, overlap is from "bigger" d1 and d2, which is usefull to "tweek"
  ///   what we call an overlap, e.g. to consider 2 pads touching only by their
  ///   corners to be overlapping.
  
  AliMpArea a1(d1.Position(),d1.Dimensions());
  AliMpArea a2(d2.Position(),d2.Dimensions());
  
  if ( a1.LeftBorder() > a2.RightBorder() - precision.X() ||
       a1.RightBorder() < a2.LeftBorder() + precision.X() )
  {
    return kFALSE;
  }
  
  if ( a1.DownBorder() > a2.UpBorder() - precision.Y() ||
       a1.UpBorder() < a2.DownBorder() + precision.Y() )
  {
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONPad::AreOverlapping(const AliMUONPad& d1, const AliMUONPad& d2,
                           const TVector2& precision,
                           AliMpArea& overlapRegion) 
{
  /// Checks the overlap between 2 pads, and returns the overlap area
  /// 
  /// See comments on the other AreOverlapping method, too : in this
  /// method, the overlapRegion does *not* depend on the precision parameter,
  /// which is only used to decide whether the pads are overlapping, while
  /// the actual overlap region is computed w/o reference to precision.
  ///
  if ( AreOverlapping(d1,d2,precision) )
  {
    overlapRegion = Overlap(d1,d2);
    if ( !overlapRegion.IsValid() )
    {
      cerr << "Something is wrong : the 2 pads below are flagged as overlapping"
      << ", but the overlapRegion is not valid"
      << endl;
      d1.Print("corners");
      d2.Print("corners");
    }
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Int_t 
AliMUONPad::Compare(const TObject* obj) const
{
  /// Compare 2 pads. 
  /// Ordering is as complete as possible.
  
  const AliMUONPad* pad = static_cast<const AliMUONPad*>(obj);
    
  if ( DetElemId() > pad->DetElemId() )
  {
    return 1;
  }
  else if ( DetElemId() > pad->DetElemId() )
  {
    return 1;
  }
  else
  {
    // same DetElemId...
    if ( Cathode() > pad->Cathode() )
    {
      return 1;
    }
    else if ( Cathode() < pad->Cathode() ) 
    {
      return -1;
    }
    else
    {
      // same cathode...
      if ( Ix() > pad->Ix() )
      {
        return 1;
      }
      else if ( Ix() < pad->Ix() ) 
      {
        return -1;
      }
      else
      {
        // same ix....
        if ( Iy() > pad->Iy() )
        {
          return 1;
        }
        else if ( Iy() < pad->Iy() ) 
        {
          return -1;
        }
        else
        {
          // same iy....
          if ( X() > pad->X() )
          {
            return 1;
          }
          else if ( X() < pad->X() )
          {
            return -1;
          }
          else
          {
            // same X
            if ( Y() > pad->Y() )
            {
              return 1;
            }
            else if ( Y() < pad->Y() )
            {
              return -1;
            }
            else
            {
              // same Y
              return ( Charge() < pad->Charge() ) ? -1:1;              
            }            
          }          
        }        
      }
    }
  }
}

//_____________________________________________________________________________
Double_t
AliMUONPad::Coord(Int_t ixy) const
{
  /// To be friendly and backward compatible with AZ code, which 
  /// used that kind of coordinate accessing.
  
  if ( ixy == 0 ) 
  {
    return X();
  }
  else if ( ixy == 1 )
  {
    return Y();
  }
  AliError(Form("Incorrect coordinates index %d (only 0,1 are valid)",ixy));
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONPad::Init(Int_t detElemId, Int_t cathode,
                 Int_t ix, Int_t iy,
                 const TVector2& position,
                 const TVector2& dimensions,
                 Double_t charge)
{
  /// Called by all the ctors
  fIsSaturated = kFALSE;
  fIsReal = kTRUE;
  fDetElemId = detElemId;
  fCathode = cathode;
  fIx = ix;
  fIy = iy;
  fPosition = position;
  fDimensions = dimensions;
  fCharge = charge;
  fChargeBackup = fCharge;
  
  fClusterId = -1;
  fDigitIndex = -1;

  fStatus = 0;
}

//_____________________________________________________________________________
AliMpArea
AliMUONPad::Overlap(const AliMUONPad& d1, const AliMUONPad& d2)
{  
  /// Return the overlap region between two pads
  AliMpArea a1(d1.Position(),d1.Dimensions());
  AliMpArea a2(d2.Position(),d2.Dimensions());
  return Intersect(a1,a2);
}


//_____________________________________________________________________________
void
AliMUONPad::Paint(Option_t*)
{
  /// Paint pad on screen
  TVector2 ll = Position() - Dimensions();
  TVector2 ur = Position() + Dimensions();

  gPad->PaintBox(ll.X(),ll.Y(),ur.X(),ur.Y());
}

//_____________________________________________________________________________
void
AliMUONPad::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  sopt.ToLower();
  
  ios_base::fmtflags oldflags = cout.flags();
  if ( Cathode() >= 0 )
  {
    cout << "DetEle " << setw(5) << DetElemId()
    << " Cath " << setw(2) << Cathode()
    << " (Ix,Iy)=(" << setw(3) << Ix() << "," << setw(3) << Iy() << ") ";
  }
  cout.setf(ios::fixed);
  cout.precision(6);
  cout << " (x,y)=(" << setw(9) << X() << "," << setw(9) << Y() << ") "
  << " (dx,dy)=(" << setw(9) << DX() << "," << setw(9) << DY() << ") "
  << " Charge=";
  cout.precision(2);
  cout << setw(7) << Charge();
  if ( sopt.Contains("full") )
  {
    cout 
    << " Used=" << (IsUsed()?Form("YES (ClusterId %d)",fClusterId):"NO")
    << (IsSaturated()?"(S)":"   ")
    << (IsReal()?"   ":"(V)")
    << " Status=" << setw(4) << Status()
    << " ChargeBackup=" << ChargeBackup();
  }
  if ( sopt.Contains("corners") )
  {
    cout << Form(" (xmin,xmax)=(%e,%e) (ymin,ymax)=(%e,%e)",
                 X()-DX(),X()+DX(),
                 Y()-DY(),Y()+DY()) << endl;
  }
  cout << endl;
  cout.flags(oldflags);
}

//_____________________________________________________________________________
void 
AliMUONPad::SetCoord(Int_t ixy, Double_t coord)
{
  /// Set the pad coordinate (ixy=0 means x, ixy=1 means y)
  if ( ixy == 0 ) 
  {
    fPosition.Set(coord,Y());
  }
  else if ( ixy == 1 )
  {
    fPosition.Set(X(),coord);
  }
  else
  {
    AliError(Form("Incorrect coordinates index %d (only 0,1 are valid)",ixy));
  }
}

//_____________________________________________________________________________
void 
AliMUONPad::SetSize(Int_t ixy, Double_t size)
{
  /// Set the pad half size (ixy=0 means x half size, ixy=1 means y half size)
  if ( ixy == 0 ) 
  {
    fDimensions.Set(size,DY());
  }
  else if ( ixy == 1 )
  {
    fDimensions.Set(DX(),size);
  }
  else
  {
    AliError(Form("Incorrect coordinates index %d (only 0,1 are valid)",ixy));
  }
}

//_____________________________________________________________________________
void 
AliMUONPad::Shift(Int_t ixy, Double_t shift)
{
  /// Shift the position by "shift"
  SetCoord(ixy,Coord(ixy)+shift);
}

//_____________________________________________________________________________
Double_t
AliMUONPad::Size(Int_t ixy) const
{
  /// Returns the half size along a direction, given by ixy
  /// (see SetSize for ixy meaning)
  
  if ( ixy == 0 ) 
  {
    return DX();
  }
  else if ( ixy == 1 )
  {
    return DY();
  }
  AliError(Form("Incorrect coordinates index %d (only 0,1 are valid)",ixy));
  return 0;
}
