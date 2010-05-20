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

#include <Riostream.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TVirtualPad.h>
#include <TVirtualX.h>

#include "AliMUONCluster.h"
#include "AliMUONPad.h"

#include "AliMpEncodePair.h"

#include "AliLog.h"

//-----------------------------------------------------------------------------
/// \class AliMUONCluster
///
/// A group of adjacent pads
///
/// Besides holding an internal array of AliMUONPads, this object
/// also computes some global characteristics for that pad sets.
///
/// \author Laurent Aphecetche
///
//-----------------------------------------------------------------------------

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
fPads(),
fHasPosition(kFALSE),
fPosition(1E9,1E9),
fPositionError(1E9,1E9),
fHasCharge(kFALSE),
fChi2(0)
{
  /// ctor
  fMultiplicity[0]=fMultiplicity[1]=0;
  fRawCharge[0]=fRawCharge[1]=0;
  fCharge[0]=fCharge[1]=0;
  fIsSaturated[0]=fIsSaturated[1]=kFALSE;
  fPads.SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONCluster::AliMUONCluster(const AliMUONCluster& src)
: TObject(src),
fPads(),
fHasPosition(kFALSE),
fPosition(1E9,1E9),
fPositionError(1E9,1E9),
fHasCharge(kFALSE),
fChi2(0)
{
  /// copy ctor
  fPads.SetOwner(kTRUE);
  src.Copy(*this);
}

//_____________________________________________________________________________
AliMUONCluster&
AliMUONCluster::operator=(const AliMUONCluster& src)
{
  /// assignement operator
  if ( this != &src ) 
  {
    src.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONCluster::~AliMUONCluster()
{
  /// dtor : note that we're owner of our pads
//   fPads.Delete();
}

//_____________________________________________________________________________
void
AliMUONCluster::Clear(Option_t*)
{
  /// Clear our pad array
  fPads.Clear();
//  fPads.Delete();
}

//_____________________________________________________________________________
Bool_t
AliMUONCluster::Contains(const AliMUONPad& pad) const
{
  /// Whether this cluster contains the pad
  if (fPads.IsEmpty()) return kFALSE;
  
  for ( Int_t i = 0; i < Multiplicity(); ++i ) 
  {
    AliMUONPad* p = Pad(i);
    if ( pad.Compare(p) == 0 ) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void
AliMUONCluster::AddCluster(const AliMUONCluster& cluster)
{
  /// Add all the pads for cluster to this one
  for ( Int_t i = 0; i < cluster.Multiplicity(); ++i )
  {
    AliMUONPad* p = cluster.Pad(i);
    if ( Contains(*p) ) 
    {
      AliError("I already got this pad : ");
      StdoutToAliError(p->Print(););
      AliFatal("");
    }
    AddPad(*p);
  }
  
}

//_____________________________________________________________________________
AliMUONPad*
AliMUONCluster::AddPad(const AliMUONPad& pad)
{
  /// Add a pad to our pad array, and update some internal information
  /// accordingly.

  AliMUONPad* p = new AliMUONPad(pad);
  fPads.AddLast(p);
  p->SetClusterId(GetUniqueID());
  Int_t cathode = p->Cathode();
  ++(fMultiplicity[cathode]);
  fRawCharge[cathode] += p->Charge();
  if ( p->IsSaturated() )
  {
    fIsSaturated[p->Cathode()]=kTRUE;
  }
  return p;
}

//___________________________________________________________________________
TString
AliMUONCluster::AsString() const
{
  /// Return a string containing a compact form of the pad list
  TString s(Form("NPADS(%d,%d)",Multiplicity(0),Multiplicity(1)));
  
  for (Int_t i = 0; i < Multiplicity(); ++i ) 
  {
    AliMUONPad* p = Pad(i);
    s += Form(" (%d,%d,%d) ",p->Cathode(),p->Ix(),p->Iy());
  }
  return s;
}


//___________________________________________________________________________
Bool_t
AliMUONCluster::AreOverlapping(const AliMUONCluster& c1, const AliMUONCluster& c2)
{
  /// Whether the two clusters overlap
  
  static Double_t precision = 1E-4; // cm
  static TVector2 precisionAdjustment(precision,precision);
    
  for ( Int_t i1 = 0; i1 < c1.Multiplicity(); ++i1 )
  {
    AliMUONPad* p1 = c1.Pad(i1);
    
    for ( Int_t i2 = 0; i2 < c2.Multiplicity(); ++i2 )
    {
      AliMUONPad* p2 = c2.Pad(i2);
      // Note: we use negative precision numbers, meaning
      // the area of the pads will be *increased* by these small numbers
      // prior to check the overlap by the AreOverlapping method,
      // so pads touching only by the corners will be considered as
      // overlapping.    
      if ( AliMUONPad::AreOverlapping(*p1,*p2,precisionAdjustment) )
      {
        return kTRUE;
      }
    }
  }
  return kFALSE;
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
  return AliMpArea((lowerLeft+upperRight).X()/2,(lowerLeft+upperRight).Y()/2, 
                   (upperRight-lowerLeft).X()/2, (upperRight-lowerLeft).Y()/2);
}

//_____________________________________________________________________________
AliMpArea
AliMUONCluster::Area(Int_t cathode) const
{
  /// Return the geometrical area covered by this cluster's pads on 
  /// a given cathode
  
  // Start by finding the (x,y) limits of this cluster
  TVector2 lowerLeft(1E9,1E9);
  TVector2 upperRight(-1E9,-1E9);
  
  for ( Int_t i = 0; i < Multiplicity(); ++i )
  {
    AliMUONPad* pad = Pad(i);
    if ( pad->Cathode() == cathode ) 
    {
      TVector2 ll = pad->Position() - pad->Dimensions();
      TVector2 ur = pad->Position() + pad->Dimensions();
      lowerLeft.Set( TMath::Min(ll.X(),lowerLeft.X()),
                     TMath::Min(ll.Y(),lowerLeft.Y()) );
      upperRight.Set( TMath::Max(ur.X(),upperRight.X()),
                      TMath::Max(ur.Y(),upperRight.Y()) );
    }
  }
  
  // then construct the area from those limits
  return AliMpArea((lowerLeft+upperRight).X()/2,(lowerLeft+upperRight).Y()/2,
                   (upperRight-lowerLeft).X()/2, (upperRight-lowerLeft).Y()/2);
}

//_____________________________________________________________________________
Bool_t
AliMUONCluster::IsMonoCathode() const
{
  /// Whether we have signals only in one of the two cathodes
  return (Cathode()<2);
}

//_____________________________________________________________________________
Int_t
AliMUONCluster::Cathode() const
{
  /// Return the cathode "number" of this cluster : 
  /// 0 if all its pads are on cathode 0
  /// 1 if all its pads are on cathode 1
  /// 2 if some pads on cath 0 and some on cath 1
  
  Int_t cathode(-1);
  if (Multiplicity(0)>0 && Multiplicity(1)>0) 
  {
    cathode=2;
  }
  else if (Multiplicity(0)>0) 
  {
    cathode=0;
  }
  else if (Multiplicity(1)>0) 
  {
    cathode=1;
  }
  
  return cathode;
}

//_____________________________________________________________________________
void
AliMUONCluster::Copy(TObject& obj) const
{
  ///
  /// Copy this cluster to (cluster&)obj
  ///
  TObject::Copy(obj);
  AliMUONCluster& dest = static_cast<AliMUONCluster&>(obj);

//  dest.fPads.Delete();
  dest.fPads.Clear();
  
  for ( Int_t i = 0; i <= fPads.GetLast(); ++i ) 
  {
    AliMUONPad* p = static_cast<AliMUONPad*>(fPads.UncheckedAt(i));
    dest.fPads.AddLast(new AliMUONPad(*p));
  }
  dest.fHasPosition = fHasPosition;
  dest.fPosition = fPosition;
  dest.fPositionError = fPositionError;
  dest.fHasCharge = fHasCharge;
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
  
  if ( Multiplicity(0) && Multiplicity(1) )
  {
    return (Charge(0)+Charge(1))/2.0;
  }
  else if ( Multiplicity(0) ) 
  {
    return Charge(0);
  }
  else if ( Multiplicity(1) ) 
  {
    return Charge(1);
  }
  AliError("Should not be here ?!");
  return -1.0;
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
AliMUONCluster::MaxPadDimensions(Int_t statusMask, Bool_t matchMask) const
{
  /// Returns the maximum pad dimensions (half sizes), only considering
  /// pads matching (or not, depending matchMask) a given mask
  
  TVector2 cath0(MaxPadDimensions(0,statusMask,matchMask)); 
  TVector2 cath1(MaxPadDimensions(1,statusMask,matchMask)); 
  
  return TVector2( TMath::Max(cath0.X(),cath1.X()),
                   TMath::Max(cath0.Y(),cath1.Y()) );
}

//_____________________________________________________________________________
TVector2
AliMUONCluster::MaxPadDimensions(Int_t cathode, 
                                 Int_t statusMask, Bool_t matchMask) const
{
  /// Returns the maximum pad dimensions (half sizes), only considering
  /// pads matching (or not, depending matchMask) a given mask, within a
  /// given cathode
  
  Double_t xmax(0);
  Double_t ymax(0);
  
  for ( Int_t i = 0; i < Multiplicity(); ++i )
  {
    AliMUONPad* pad = Pad(i);
    if ( ShouldUsePad(*pad,cathode,statusMask,matchMask) )
    {
      xmax = TMath::Max(xmax,pad->DX());
      ymax = TMath::Max(ymax,pad->DY());
    }
  }
  return TVector2(xmax,ymax);
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
Long_t
AliMUONCluster::NofPads(Int_t statusMask, Bool_t matchMask) const
{
  /// Number of pads satisfying (or not, depending matchMask) a
  /// given mask 
  
  Int_t nx, ny;
  
  TVector2 dim0(MinPadDimensions(0,statusMask,matchMask));
  TVector2 dim1(MinPadDimensions(1,statusMask,matchMask));
  
  Long_t npad0(NofPads(0,statusMask,matchMask));
  Long_t npad1(NofPads(1,statusMask,matchMask));
  
  if ( TMath::Abs( (dim0-dim1).X() ) < 1E-3 )
  {
    nx = TMath::Max( AliMp::PairFirst(npad0), AliMp::PairFirst(npad1) );
  }
  else
  {
    nx = dim0.X() < dim1.X() ? AliMp::PairFirst(npad0) : AliMp::PairFirst(npad1);
  }
  
  if ( TMath::Abs( (dim0-dim1).Y() ) < 1E-3 )
  {
    ny = TMath::Max( AliMp::PairSecond(npad0), AliMp::PairSecond(npad1) );
  }
  else
  {
    ny = dim0.Y() < dim1.Y() ? AliMp::PairSecond(npad0) : AliMp::PairSecond(npad1);
  }
  
  return AliMp::Pair(nx,ny);
}

//_____________________________________________________________________________
Long_t
AliMUONCluster::NofPads(Int_t cathode,
                        Int_t statusMask, Bool_t matchMask) const
{
  /// Number of pads of a given cathode, satisfying (or not, 
  /// depending matchMask) a given mask

  Int_t n = Multiplicity(cathode);
  if (!n) 
  {
    return 0;
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
  
  return AliMp::Pair(cx,cy);
}

//_____________________________________________________________________________
AliMUONPad*
AliMUONCluster::Pad(Int_t index) const
{
  /// Returns the index-th pad
  
  if (fPads.IsEmpty()) return 0x0;
  if ( index < fPads.GetLast()+1 )
  {
    return static_cast<AliMUONPad*>(fPads.At(index));
  }
  else
  {
    AliError(Form("Requested index %d out of bounds (%d) Mult is %d",index,
                  fPads.GetLast(),Multiplicity()));
    DumpMe();
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
AliMUONCluster::DumpMe() const
{
  /// printout
  cout << "Cluster Id " << GetUniqueID() << " npads=" << Multiplicity() 
  << "(" << Multiplicity(0) << "," << Multiplicity(1) << ") RawCharge=" 
  << RawCharge() << " (" << RawCharge(0) << "," << RawCharge(1)
  << ") Charge=(" << Charge(0) << "," << Charge(1) <<")";
  if ( HasPosition() )
  {
    cout << " (x,y)=(" << Position().X() << "," << Position().Y() << ")";
    cout << " (errX,errY)=(" << PositionError().X() << "," << PositionError().Y() << ")";
  }
  cout << endl;
//  cout << " " << Area() << endl;
  for (Int_t i = 0; i < fPads.GetSize(); ++i) 
  {
    cout << Form("fPads[%d]=%x",i,fPads.At(i)) << endl;
    if ( fPads.At(i) ) fPads.At(i)->Print();
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
  << ") Charge=(" << Charge(0) << "," << Charge(1) <<")";
  if ( HasPosition() )
  {
    cout << " (x,y)=(" << Position().X() << "," << Position().Y() << ")";
    cout << " (errX,errY)=(" << PositionError().X() << "," << PositionError().Y() << ")";
  }
  cout << " " << Area();

  TObjArray* a = static_cast<TObjArray*>(fPads.Clone());
  a->Sort();
  a->Print("",opt);
  delete a;
}

//_____________________________________________________________________________
//Bool_t
//AliMUONCluster::IsEqual(const TObject* obj) const
//{
//  const AliMUONCluster* c = static_cast<const AliMUONCluster*>(obj);
//  if ( c->Multiplicity() != Multiplicity() ) return kFALSE;
//  
//  for ( Int_t i = 0; i < c->Multiplicity(); ++i ) 
//  {
//    AliMUONPad* p = c->Pad(i);
//    if ( p->Compare(Pad(i)) ) return kFALSE;
//  }
//  return kTRUE;
//}

//_____________________________________________________________________________
Int_t 
AliMUONCluster::Compare(const TObject* obj) const
{
  /// Compare two clusters. Comparison is made on position and rawcharge only.
  
  const AliMUONCluster* cluster = static_cast<const AliMUONCluster*>(obj);
  
  AliMpArea carea(cluster->Area());
  AliMpArea area(Area());

  if ( carea.GetPositionX() > area.GetPositionX() ) 
  {
    return 1;
  }
  else if ( carea.GetPositionX() < area.GetPositionX() ) 
  {
    return -1;
  }
  else 
  {
    if ( carea.GetPositionY() > area.GetPositionY() ) 
    {
      return 1;
    }
    else if ( carea.GetPositionY() < area.GetPositionY() ) 
    {
      return -1;
    }
    else
    {
      if ( cluster->RawCharge() > RawCharge() ) 
      {
        return 1;
      }
      else if ( cluster->RawCharge() < RawCharge() )
      {
        return -1;
      }
    }
  }
  return 0;
}

//_____________________________________________________________________________
void
AliMUONCluster::RemovePad(AliMUONPad* pad)
{
  /// Remove a pad. 
  /// As a consequence, some internal information must be updated
  
  fPads.Remove(pad);
  fPads.Compress();
  delete pad;
  // update cluster's data
  fIsSaturated[0]=fIsSaturated[1]=kFALSE;
  fMultiplicity[0]=fMultiplicity[1]=0;
  fRawCharge[0]=fRawCharge[1]=0;
  for ( Int_t i = 0; i <= fPads.GetLast(); ++i )
  {
    AliMUONPad* p = Pad(i);
    if ( p->IsSaturated() ) 
    {
      fIsSaturated[p->Cathode()] = kTRUE;
    }
    ++fMultiplicity[p->Cathode()];
    fRawCharge[p->Cathode()] += p->Charge();
  }
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
