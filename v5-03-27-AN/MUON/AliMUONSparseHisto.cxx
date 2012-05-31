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

#include "AliMUONSparseHisto.h"

#include "AliLog.h"
#include <Riostream.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>

/// \class AliMUONSparseHisto
///
/// Tiny histogram-like class to hold some distributions of tracker data.
/// Only intent of this class is to minimize memory consumption, in
/// order to fit a maximum number of channel histograms into memory.
/// The rest is not supported ;-)
///
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONSparseHisto)
/// \endcond

//______________________________________________________________________________
AliMUONSparseHisto::AliMUONSparseHisto(Double_t xmin, Double_t xmax)
: TObject(),
fNbins(0),
fArray(0x0),
fXmin(xmin),
fXmax(xmax),
fFactor((1<<Nbits())/(xmax-xmin))
{
  /// ctor
  SetBit(kOverflow,0);
  SetBit(kUnderflow,0);
}

//______________________________________________________________________________
AliMUONSparseHisto::AliMUONSparseHisto(const AliMUONSparseHisto& rhs)
: TObject(rhs),
fNbins(0),
fArray(0x0),
fXmin(0.0),
fXmax(0.0),
fFactor(0.0)
{
  /// copy ctor
  rhs.Copy(*this);
}

//______________________________________________________________________________
AliMUONSparseHisto&
AliMUONSparseHisto::operator=(const AliMUONSparseHisto& rhs)
{
  /// assignment operator
  if ( this != &rhs )
  {
    rhs.Copy(*this);
  }
  return *this;
}

//______________________________________________________________________________
AliMUONSparseHisto::~AliMUONSparseHisto()
{
  /// dtor
  delete[] fArray;
}

//______________________________________________________________________________
Bool_t 
AliMUONSparseHisto::Add(const AliMUONSparseHisto& h)
{
  /// Add h to this

  if ( fXmin != h.Xmin() || fXmax != h.Xmax() ) 
  {
    AliError("Cannot add sparse histograms with different limits !");
    return kFALSE;
  }
  
  for ( Int_t i = 0; i < h.GetNbins(); ++i ) 
  {
    Fill(h.GetBinContent(i));
  }  
  
  return kTRUE;
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Clear(Option_t*)
{
  /// Reset the content
  delete[] fArray;
  fArray = 0x0;
  fNbins = 0;
}  

//______________________________________________________________________________
void 
AliMUONSparseHisto::Copy(TObject& object) const
{
  /// Copy this to *object
  TObject::Copy(object);
  AliMUONSparseHisto& h = static_cast<AliMUONSparseHisto&>(object);
  delete[] h.fArray;
  h.fArray = 0x0;
  h.fNbins = GetNbins();
  h.fXmin = Xmin();
  h.fXmax = Xmax();
  h.fFactor = Factor();
  
  if ( GetNbins() > 0 )
  {
    h.fArray = new UInt_t[GetNbins()];
    for ( Int_t i = 0; i < GetNbins(); ++i ) 
    {
      h.fArray[i] = GetBin(i);
    }
  }
}

//______________________________________________________________________________
Double_t 
AliMUONSparseHisto::DecodeValue(Int_t value) const
{
  /// From internal integer to "original" double
  return value/Factor() + Xmin();
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::EncodeValue(Double_t value) const
{
  /// From original double value to internal integer
  return TMath::Nint(Factor()*(value-Xmin()));
}

//______________________________________________________________________________
Int_t
AliMUONSparseHisto::BinCenter(UInt_t x) const
{
  /// Extract binCenter part from x
  
  return ( x & 0xFFF00000 ) >> 20;
}

//______________________________________________________________________________
Int_t
AliMUONSparseHisto::BinContent(UInt_t x) const
{
  /// Extract binContent part from x
  
  return (x & 0xFFFFF);
}

//______________________________________________________________________________
UInt_t 
AliMUONSparseHisto::Encode(Int_t binCenter, Int_t binContent) const
{
  /// Convert (binCenter,binContent) into a single value
  
  return ( ( binCenter & 0xFFF ) ) << 20 | ( ( binContent & 0xFFFFF ) );
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Expand()
{
  /// Make fArray of size n
  if (!fArray || !fNbins) 
  {
    delete[] fArray;
    fArray = new UInt_t[1];
    fNbins = 1;
  }
  else
  {
    UInt_t* tmp = new UInt_t[fNbins+1];
    for ( Int_t i = 0; i < fNbins; ++i ) 
    {
      tmp[i] = fArray[i];
    }
    delete[] fArray;
    fArray = tmp;
    ++fNbins;
  }
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::Fill(Double_t value)
{
  /// Fill
  
  if ( value < Xmin() ) 
  {
    SetBit(kUnderflow,1);
    return -1;
  }
  
  if ( value > Xmax() ) 
  {
    SetBit(kOverflow,1);
    return -1;
  }
  
  Int_t ivalue = EncodeValue(value);
  
  Int_t i = Find(ivalue);
  
  if ( i < 0 ) 
  {
    Int_t n = fNbins;
    Expand();
    fArray[n] = Encode(ivalue,1);
    i = n;
  }
  else
  {
    Int_t bc = GetBinContent(i);
    if ( bc < 0xFFFFF ) 
    {
      fArray[i] = Encode(ivalue,bc+1);
    }
  }
    
  return i;
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::Find(Int_t binCenter) const
{
  /// Return the index in fArray of value, or -1 if not found
  
  for ( Int_t i = 0; i < GetNbins(); ++i ) 
  {
    if ( binCenter == GetBinCenter(i) ) return i;
  }
  return -1;
}

//______________________________________________________________________________
UInt_t 
AliMUONSparseHisto::GetBin(Int_t bin) const
{
  /// Get bin, which is a compacted form of two integers : (binCenter,binContent)
  /// where binCenter itself might be an integer-fied double value.
  return fArray[bin];
}

//______________________________________________________________________________
Double_t 
AliMUONSparseHisto::GetBinCenter(Int_t bin) const
{
  /// Get bin center
  if ( bin < 0 ) return -FLT_MAX;
  if ( bin >= GetNbins() ) return FLT_MAX;
  
  UInt_t i = GetBin(bin);
  
  return DecodeValue(BinCenter(i));
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::GetBinContent(Int_t bin) const
{
  /// Get bin content
  
  if ( bin < 0 || bin >= GetNbins() ) return 0xFFFFFFFF;

  UInt_t i = GetBin(bin);

  return BinContent(i);
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Print(Option_t* opt) const
{
  /// Printout
  Int_t id1 = ( GetUniqueID() & 0xFFFF0000 ) >> 16;
  Int_t id2 = GetUniqueID() & 0xFFFF;
  
  cout << "ID=(" << id1 << "," << id2 << ") n bins = " << GetNbins();
  if ( HasUnderflow() ) cout << " has underflow(s)";
  if ( HasOverflow() ) cout << " has overflow(s)";
  cout << endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("FULL") )
  {
    for ( Int_t i = 0; i < GetNbins(); ++i ) 
    {
      cout << Form("Bin (%10u) %e = %6d",GetBin(i),GetBinCenter(i),GetBinContent(i)) << endl;
    }
  }
}
