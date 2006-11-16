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

#include "AliMUONCluster.h"

#include "AliLog.h"
#include "AliMUONPad.h"
#include "TObjArray.h"
#include "Riostream.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"

/// \class AliMUONCluster
///
/// A group of adjacent pads
///
/// Besides holding an internal array of AliMUONPads, this object
/// also computes some global characteristics for that pad sets.
///
/// \author Laurent Aphecetche
///

/// \cond CLASSIMP
ClassImp(AliMUONCluster)
/// \endcond

namespace
{
  //___________________________________________________________________________
  Bool_t
  ShouldUsePad(const AliMUONPad& pad, 
               Int_t cathode, Int_t statusMask, Bool_t matchMask)
  {

      // FIXME : we should only use >=0 status, so we can fully
      // use masking possibility ?
    if ( pad.Status() < 0 ) return kFALSE;
    
    if ( pad.Cathode() == cathode && pad.IsReal() && !pad.IsSaturated() )
    {
      Bool_t test = ( ( pad.Status() & statusMask ) != 0 );
      if ( !statusMask ) 
      {
        test = ( pad.Status() == 0 );
      }
      if ( ( test && matchMask ) || ( !test && !matchMask ) )
      {
        return kTRUE;
      }
    }
    return kFALSE;
  }

  //___________________________________________________________________________
  Int_t Unique(Int_t n, Double_t* array, Double_t precision)
  {
    /// Return the number of *different* elements in array 
    /// where different is up to precision
    /// Note that we assume that n is >= 1
    
    Int_t count(1);
    
    Int_t* index = new Int_t[n];
    
    TMath::Sort(n,array,index);
        
    for ( Int_t i = 1; i < n; ++i )
    {
      if ( array[index[i]] - array[index[i-1]] < -precision ) ++count;
    }
    
    delete[] index;
        
    return count;
  }
}

//_____________________________________________________________________________
AliMUONCluster::AliMUONCluster() 
: TObject(), 
fPads(0x0),
fHasPosition(kFALSE),
fPosition(1E9,1E9),
fPositionError(1E9,1E9),
fHasCharge(kFALSE),
fChi2(0),
fIsSorted(kFALSE)
{
  /// ctor
  fMultiplicity[0]=fMultiplicity[1]=0;
  fRawCharge[0]=fRawCharge[1]=0;
  fCharge[0]=fCharge[1]=0;
  fIsSaturated[0]=fIsSaturated[1]=kFALSE;
}

//_____________________________________________________________________________
AliMUONCluster::AliMUONCluster(const AliMUONCluster& src)
: TObject(src),
fPads(0x0),
fHasPosition(kFALSE),
fPosition(1E9,1E9),
fPositionError(1E9,1E9),
fHasCharge(kFALSE),
fChi2(0),
fIsSorted(kFALSE)
{
  /// copy ctor
  src.Copy(*this);
}

//_____________________________________________________________________________
AliMUONCluster&
AliMUONCluster::operator=(const AliMUONCluster& src)
{
  /// assignement operator
  AliMUONCluster c(src);
  c.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONCluster::~AliMUONCluster()
{
  /// dtor : note that we're owner of our pads
  delete fPads;
}

//_____________________________________________________________________________
void
AliMUONCluster::AddPad(const AliMUONPad& pad)
{
  /// Add a pad to our pad array, and update some internal information
  /// accordingly.
  /// If pad array was sorted prior to this call, we re-sort it after
  /// actual addition.
  
  if (!fPads) 
  {
    fPads = new TObjArray(10);
    fPads->SetOwner(kTRUE);
  }
  AliMUONPad* p = new AliMUONPad(pad);
  fPads->AddLast(p);
  p->SetClusterId(GetUniqueID());
  Int_t cathode = p->Cathode();
  ++(fMultiplicity[cathode]);
  fRawCharge[cathode] += p->Charge();
  if ( p->IsSaturated() )
  {
    fIsSaturated[p->Cathode()]=kTRUE;
  }
  if ( fIsSorted ) Sort();
}

//_____________________________________________________________________________
AliMpArea
AliMUONCluster::Area() const
{
  /// Return the geometrical area covered by this cluster
  
  // Start by finding the (x,y) limits of this cluster
  TVector2 lowerLeft(1E9,1E9);
  TVector2 upperRight(-1E9,-1E9);
  
  for ( Int_t i = 0; i < Multiplicity(); ++i )
  {
    AliMUONPad* pad = Pad(i);
    TVector2 ll = pad->Position() - pad->Dimensions();
    TVector2 ur = pad->Position() + pad->Dimensions();
    lowerLeft.Set( TMath::Min(ll.X(),lowerLeft.X()),
                   TMath::Min(ll.Y(),lowerLeft.Y()) );
    upperRight.Set( TMath::Max(ur.X(),upperRight.X()),
                    TMath::Max(ur.Y(),upperRight.Y()) );
  }

  // then construct the area from those limits
  return AliMpArea((lowerLeft+upperRight)/2,
                   (upperRight-lowerLeft)/2);
}

//_____________________________________________________________________________
void
AliMUONCluster::Copy(TObject& obj) const
{
  //
  // Copy this cluster to (cluster&)obj
  //
  TObject::Copy(obj);
  AliMUONCluster& dest = static_cast<AliMUONCluster&>(obj);
  dest.fPads = static_cast<TObjArray*>(fPads->Clone());
  dest.fHasPosition = fHasPosition;
  dest.fPosition = fPosition;
  dest.fPositionError = fPositionError;
  dest.fHasCharge = fHasCharge;
  dest.fIsSorted = fIsSorted;
  dest.fChi2 = fChi2;
  for ( Int_t i = 0; i < 2; ++i )
  {
    dest.fRawCharge[i] = fRawCharge[i];
    dest.fCharge[i] = fCharge[i];
    dest.fMultiplicity[i] = fMultiplicity[i];
    dest.fIsSaturated[i] = fIsSaturated[i];
  }
}

//_____________________________________________________________________________
Float_t 
AliMUONCluster::Charge() const
{
  /// Return the average charge over both cathodes
  return (Charge(0)+Charge(1))/2.0;
}

//_____________________________________________________________________________
Float_t
AliMUONCluster::Charge(Int_t cathode) const
{
  /// Returns the charge of a given cathode
  if ( !fHasCharge ) return RawCharge(cathode);
  
  if ( cathode == 0 || cathode == 1 )
  {
    return fCharge[cathode];
  }
  return 0;
}

//_____________________________________________________________________________
Float_t
AliMUONCluster::ChargeAsymmetry() const
{
  /// Returns the charge asymmetry
  if ( Charge() > 0 )
  {
    return TMath::Abs(Charge(0)-Charge(1))/Charge();
  }
  return 0;
}

//_____________________________________________________________________________
TVector2
AliMUONCluster::MinPadDimensions(Int_t statusMask, Bool_t matchMask) const
{
  /// Returns the minimum pad dimensions (half sizes), only considering
  /// pads matching (or not, depending matchMask) a given mask
  
  TVector2 cath0(MinPadDimensions(0,statusMask,matchMask)); 
  TVector2 cath1(MinPadDimensions(1,statusMask,matchMask)); 
  
  return TVector2( TMath::Min(cath0.X(),cath1.X()),
                   TMath::Min(cath0.Y(),cath1.Y()) );
}

//_____________________________________________________________________________
TVector2
AliMUONCluster::MinPadDimensions(Int_t cathode, 
                                 Int_t statusMask, Bool_t matchMask) const
{
  /// Returns the minimum pad dimensions (half sizes), only considering
  /// pads matching (or not, depending matchMask) a given mask, within a
  /// given cathode
  
  Double_t xmin(1E9);
  Double_t ymin(1E9);
    
  for ( Int_t i = 0; i < Multiplicity(); ++i )
  {
    AliMUONPad* pad = Pad(i);
    if ( ShouldUsePad(*pad,cathode,statusMask,matchMask) )
    {
      xmin = TMath::Min(xmin,pad->DX());
      ymin = TMath::Min(ymin,pad->DY());
    }
  }
  return TVector2(xmin,ymin);
}

//_____________________________________________________________________________
Int_t 
AliMUONCluster::Multiplicity() const
{
  /// Returns the total number of pads in this cluster
  return Multiplicity(0)+Multiplicity(1);
}

//_____________________________________________________________________________
Int_t
AliMUONCluster::Multiplicity(Int_t cathode) const
{
  /// Returns the number of pads in this cluster, in the given cathode
  if ( cathode == 0 || cathode == 1 )
  {
    return fMultiplicity[cathode];
  }
  return 0;
}

//_____________________________________________________________________________
AliMpIntPair
AliMUONCluster::NofPads(Int_t statusMask, Bool_t matchMask) const
{
  /// Number of pads satisfying (or not, depending matchMask) a
  /// given mask
  
  Int_t nx, ny;
  
  TVector2 dim0(MinPadDimensions(0,statusMask,matchMask));
  TVector2 dim1(MinPadDimensions(1,statusMask,matchMask));
  
  AliMpIntPair npad0(NofPads(0,statusMask,matchMask));
  AliMpIntPair npad1(NofPads(1,statusMask,matchMask));
  
  if ( TMath::Abs( (dim0-dim1).X() ) < 1E-3 )
  {
    nx = TMath::Max( npad0.GetFirst(), npad1.GetFirst() );
  }
  else
  {
    nx = dim0.X() < dim1.X() ? npad0.GetFirst() : npad1.GetFirst();
  }
  
  if ( TMath::Abs( (dim0-dim1).Y() ) < 1E-3 )
  {
    ny = TMath::Max( npad0.GetSecond(), npad1.GetSecond() );
  }
  else
  {
    ny = dim0.Y() < dim1.Y() ? npad0.GetSecond() : npad1.GetSecond();
  }
  
  return AliMpIntPair(nx,ny);
}

//_____________________________________________________________________________
AliMpIntPair
AliMUONCluster::NofPads(Int_t cathode,
                        Int_t statusMask, Bool_t matchMask) const
{
  /// Number of pads of a given cathode, satisfying (or not, 
  /// depending matchMask) a given mask

  Int_t n = Multiplicity(cathode);
  if (!n) 
  {
    return AliMpIntPair(0,0);
  }
  Double_t* x = new Double_t[n];
  Double_t* y = new Double_t[n];
  Int_t np(0);
  
  for ( Int_t i = 0; i < Multiplicity(); ++i )
  {
    AliMUONPad* pad = Pad(i);
    if ( ShouldUsePad(*pad,cathode,statusMask,matchMask) )
    {
      x[np] = pad->X();
      y[np] = pad->Y();
      ++np;
    }
  }
  
  Int_t cx = Unique(np,x,0.01);
  Int_t cy = Unique(np,y,0.01);
  
  delete[] x;
  delete[] y;
  
  return AliMpIntPair(cx,cy);
}

//_____________________________________________________________________________
AliMUONPad*
AliMUONCluster::Pad(Int_t index) const
{
  /// Returns the index-th pad
  
  if (!fPads) return 0x0;
  if ( index < fPads->GetLast()+1 )
  {
    return static_cast<AliMUONPad*>(fPads->At(index));
  }
  else
  {
    AliError(Form("Requesting index %d out of bounds (%d)",index,fPads->GetLast()));
  }
  return 0x0;
}


//_____________________________________________________________________________
void
AliMUONCluster::Paint(Option_t*)
{
  /// Paint this cluster   
  if (!Multiplicity()) return;
  
  AliMpArea area(Area());
  
  gPad->Range(area.LeftBorder(),area.DownBorder(),area.RightBorder(),area.UpBorder());
      
  gVirtualX->SetFillStyle(0);
  
  gVirtualX->SetLineColor(2);
  gVirtualX->SetLineWidth(4);  
  for ( Int_t i = 0; i < Multiplicity(); ++i)
  {
    AliMUONPad* pad = Pad(i);
    if ( pad->Cathode() == 0 ) pad->Paint();
  }

  gVirtualX->SetLineColor(4);
  gVirtualX->SetLineWidth(2);  
  for ( Int_t i = 0; i < Multiplicity(); ++i)
  {
    AliMUONPad* pad = Pad(i);
    if ( pad->Cathode() == 1 ) pad->Paint();
  }
  
}

//_____________________________________________________________________________
void
AliMUONCluster::Print(Option_t* opt) const
{
  /// printout
  cout << "Cluster Id " << GetUniqueID() << " npads=" << Multiplicity() 
  << "(" << Multiplicity(0) << "," << Multiplicity(1) << ") RawCharge=" 
  << RawCharge() << " (" << RawCharge(0) << "," << RawCharge(1)
  << " Charge=(" << Charge(0) << "," << Charge(1) <<")";
  if ( HasPosition() )
  {
    cout << " (x,y)=(" << Position().X() << "," << Position().Y() << ")";
    cout << " (errX,errY)=(" << PositionError().X() << "," << PositionError().Y() << ")";
  }
  AliMpArea a(Area());
  cout << Form(" Area=(%e,%e,%e,%e)",a.LeftBorder(),a.RightBorder(),
               a.DownBorder(),a.UpBorder());
  cout << endl;
  if (fPads) 
  {
    fPads->Print("",opt);
  }
}

//_____________________________________________________________________________
void
AliMUONCluster::Sort()
{
  /// Sort the pad array
  fPads->Sort();
  fIsSorted = kTRUE;
}

//_____________________________________________________________________________
void
AliMUONCluster::RemovePad(AliMUONPad* pad)
{
  /// Remove a pad. 
  /// As a consequence, some internal information must be updated
  
  fPads->Remove(pad);
  fPads->Compress();
  // update cluster's data
  fIsSaturated[0]=fIsSaturated[1]=kFALSE;
  fMultiplicity[0]=fMultiplicity[1]=0;
  fRawCharge[0]=fRawCharge[1]=0;
  for ( Int_t i = 0; i <= fPads->GetLast(); ++i )
  {
    AliMUONPad* p = Pad(i);
    if ( p->IsSaturated() ) 
    {
      fIsSaturated[p->Cathode()] = kTRUE;
    }
    ++fMultiplicity[p->Cathode()];
    fRawCharge[p->Cathode()] += p->Charge();
  }
  if (fIsSorted) Sort();
}

//_____________________________________________________________________________
Float_t 
AliMUONCluster::RawCharge() const
{
  /// Returns the raw average charge
  return (RawCharge(0)+RawCharge(1))/2.0;
}

//_____________________________________________________________________________
Float_t
AliMUONCluster::RawCharge(Int_t cathode) const
{
  /// Returns the average charge of a given cathode
  if ( cathode == 0 || cathode == 1 )
  {
    return fRawCharge[cathode];
  }
  return 0;
}

//_____________________________________________________________________________
Float_t
AliMUONCluster::RawChargeAsymmetry() const
{
  /// Returns the raw charge asymmetry
   if ( RawCharge() > 0 )
   {
     return TMath::Abs(RawCharge(0)-RawCharge(1))/RawCharge();
   }
  return 0;
}
