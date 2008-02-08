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
#include <TString.h>
#include <Riostream.h>
#include <TH1.h>

/// \class AliMUONSparseHisto
///
/// Tiny histogram-like class to hold adc distributions of tracker data.
/// Only intent of this class is to minimize memory consumption, in
/// order to fit a maximum number of channel histograms into memory.
/// The rest is not supported ;-)
///
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONSparseHisto)
/// \endcond

//______________________________________________________________________________
AliMUONSparseHisto::AliMUONSparseHisto()
: TObject(),
fNbins(0),
fArray(0x0)
{
  /// ctor
}

//______________________________________________________________________________
AliMUONSparseHisto::AliMUONSparseHisto(const AliMUONSparseHisto& rhs)
: TObject(rhs),
fNbins(0),
fArray(0x0)
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
  if ( GetNbins() > 0 )
  {
    h.fArray = new Int_t[GetNbins()];
    for ( Int_t i = 0; i < GetNbins(); ++i ) 
    {
      h.fArray[i] = GetBinContent(i);
    }
  }
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Decode(Int_t value, Int_t& adc, Int_t& count) const
{
  /// Convert value into (adc,count) pair
  
  adc   = ( value & 0xFFF00000 ) >> 20;
  count = ( value & 0x000FFFFF );
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::Encode(Int_t adc, Int_t count) const
{
  /// Convert (adc,count) into a single value
  return ( ( adc & 0xFFF ) ) << 20 | ( count & 0xFFFFF );
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Expand()
{
  /// Make fArray of size n
  if (!fArray || !fNbins) 
  {
    delete[] fArray;
    fArray = new Int_t[1];
    fNbins = 1;
  }
  else
  {
    Int_t* tmp = new Int_t[fNbins+1];
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
AliMUONSparseHisto::Fill(Int_t adc)
{
  /// Fill
  
  if ( adc < 0 || adc > 4095 ) return -1;
  
  Int_t i = Find(adc);
  
  if ( i < 0 ) 
  {
    Int_t n = fNbins;
    Expand();
    fArray[n] = Encode(adc,1);
    i = n;
  }
  else
  {
    Int_t iadc,icontent;
    Decode(fArray[i],iadc,icontent);
    fArray[i] = Encode(adc,icontent+1);
  }
    
  return i;
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::Find(Int_t adc) const
{
  /// Return the index in fArray of adc, or -1 if not found
  for ( Int_t i = 0; i < GetNbins(); ++i ) 
  {
    Int_t content = GetBinContent(i);
    Int_t iadc,value;
    Decode(content,iadc,value);
    if ( iadc == adc ) return i;
  }
  return -1;
}

//______________________________________________________________________________
Int_t 
AliMUONSparseHisto::GetBinContent(Int_t bin) const
{
  /// Get bin content. Note that the content is compacted, so you must
  /// use Decode() method to get (adc,count) values.
  if ( bin >= 0 && bin < GetNbins() ) return fArray[bin];
  return 0;
}

//______________________________________________________________________________
void 
AliMUONSparseHisto::Print(Option_t* opt) const
{
  /// Printout
  Int_t id1 = ( GetUniqueID() & 0xFFFF0000 ) >> 16;
  Int_t id2 = GetUniqueID() & 0xFFFF;
  
  cout << "(" << id1 << "," << id2 << ") n bins = " << GetNbins() << endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("FULL") )
  {
    for ( Int_t i = 0; i < GetNbins(); ++i ) 
    {
      Int_t content = GetBinContent(i);
      Int_t adc,value;
      Decode(content,adc,value);
      cout << Form("i %4d content %10x adc %4d value %6d",i,content,adc,value)
        << endl;
    }
  }
}
