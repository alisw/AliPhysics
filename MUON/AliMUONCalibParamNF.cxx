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

#include "AliMUONCalibParamNF.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TMath.h"
#include "TString.h"

#include <limits.h>

//-----------------------------------------------------------------------------
/// \class AliMUONCalibParamNF
///
/// Handle the case of N floating point parameters per channel.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONCalibParamNF)
/// \endcond

//_____________________________________________________________________________
AliMUONCalibParamNF::AliMUONCalibParamNF() 
: AliMUONVCalibParam(),
  fDimension(0),
  fSize(0),
  fN(0),
  fValues(0x0)
{
/// Default constructor.
}

//_____________________________________________________________________________
AliMUONCalibParamNF::AliMUONCalibParamNF(Int_t dimension, Int_t theSize, 
                                         Int_t id0, Int_t id1,
                                         Float_t fillWithValue) 
: AliMUONVCalibParam(id0,id1),
  fDimension(dimension),
  fSize(theSize),
  fN(fSize*fDimension),
  fValues(0x0)
{
/// Normal constructor, where theSize specifies the number of channels handled
/// by this object, and fillWithValue the default value assigned to each
/// channel.

  if ( fN > 0 )
  {
    fValues = new Float_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      fValues[i] = fillWithValue;
    }
  }
}


//_____________________________________________________________________________
AliMUONCalibParamNF::AliMUONCalibParamNF(const AliMUONCalibParamNF& other) 
: AliMUONVCalibParam(),
fDimension(0),
fSize(0),
fN(0),
fValues(0x0)
{
/// Copy constructor.

  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUONCalibParamNF&
AliMUONCalibParamNF::operator=(const AliMUONCalibParamNF& other) 
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONCalibParamNF::~AliMUONCalibParamNF()
{
/// Destructor

  delete[] fValues;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::CopyTo(AliMUONCalibParamNF& destination) const
{
/// Copy *this to destination

  const TObject& o = static_cast<const TObject&>(*this);
  o.Copy(destination);
  
  delete[] destination.fValues;
  destination.fN = fN;
  destination.fSize = fSize;
  destination.fDimension = fDimension;
  
  if ( fN > 0 )
  {
    destination.fValues = new Float_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      destination.fValues[i] = fValues[i];
    }
  }
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNF::Index(Int_t i, Int_t j) const
{
/// Compute the 1D index of the internal storage from the pair (i,j)
/// Returns -1 if the (i,j) pair is invalid

  if ( i >= 0 && i < Size() && j >= 0 && j < Dimension() )
  {
    return IndexFast(i,j);
  }
  return -1;
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNF::IndexFast(Int_t i, Int_t j) const
{
  /// Compute the 1D index of the internal storage from the pair (i,j)
  
  return i + Size()*j;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::Print(Option_t* opt) const
{
/// Output this object to stdout.
/// If opt=="full", then all channels are printed, 
/// if opt=="mean#", only the mean and sigma value are printed for j-th dimension
/// otherwise only the general characteristics are printed.

  TString sopt(opt);
  sopt.ToUpper();
  cout << Form("AliMUONCalibParamNF Id=(%d,%d) Size=%d Dimension=%d",ID0(),
               ID1(),Size(),Dimension()) << endl;
  if ( sopt.Contains("FULL") )
  {
    for ( Int_t i = 0; i < Size(); ++i )
    {
      cout << Form("CH %3d",i);
      for ( Int_t j = 0; j < Dimension(); ++j )
      {
        cout << Form(" %e",ValueAsFloat(i,j));
      }
      cout << endl;
    }
  }  
  if ( sopt.Contains("MEAN") )
  {
    Int_t j(0);
    sscanf(sopt.Data(),"MEAN%d",&j);
    
    Float_t mean(0);
    Float_t v2(0);
    
    Int_t n = Size();
    
    for ( Int_t i = 0; i < Size(); ++i )
    {
      Float_t v = ValueAsFloat(i,j);
      mean += v;
      v2 += v*v;
    }
    mean /= n;
    float sigma = 0;
    if ( n > 1 ) sigma = TMath::Sqrt( (v2-n*mean*mean)/(n-1) );
    cout << Form(" Mean(j=%d)=%e Sigma(j=%d)=%e",j,mean,j,sigma) << endl;
  }
  
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::SetValueAsFloat(Int_t i, Int_t j, Float_t value)
{
/// Set one value as a float, after checking that the indices are correct.

  Int_t ix = Index(i,j);
  
  if ( ix < 0 )
  {
    AliError(Form("Invalid (i,j)=(%d,%d) max allowed is (%d,%d)",
                  i,j,Size()-1,Dimension()-1));
  }
  else
  {
    fValues[ix]=value;
  }
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::SetValueAsFloatFast(Int_t i, Int_t j, Float_t value)
{
  /// Set one value as a float, w/o checking that the indices are correct.
  
  fValues[IndexFast(i,j)] = value;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::SetValueAsInt(Int_t i, Int_t j, Int_t value)
{
/// Set one value as an int.

  SetValueAsFloat(i,j,static_cast<Float_t>(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParamNF::SetValueAsIntFast(Int_t i, Int_t j, Int_t value)
{
  /// Set one value as an int.
  
  SetValueAsFloatFast(i,j,static_cast<Float_t>(value));
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamNF::ValueAsFloat(Int_t i, Int_t j) const
{
/// Return the value as a float (which it is), after checking indices.

  Int_t ix = Index(i,j);
  
  if ( ix < 0 )
  {
    AliError(Form("Invalid (i,j)=(%d,%d) max allowed is (%d,%d)",
                  i,j,Size()-1,Dimension()-1));
    return 0.0;
  }
  else
  {
    return fValues[ix];
  }
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamNF::ValueAsFloatFast(Int_t i, Int_t j) const
{
  /// Return the value as a float (which it is), after checking indices.
  
  return fValues[IndexFast(i,j)];
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNF::ValueAsInt(Int_t i, Int_t j) const
{
/// Return the value as an int, by rounding the internal float value.

  Float_t v = ValueAsFloat(i,j);
  
  if ( v >= Float_t(INT_MAX) ) {
    AliErrorStream() 
      << "Cannot convert value " << v << " to Int_t." << endl;
    return 0;
  }  
  
  return TMath::Nint(v);
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNF::ValueAsIntFast(Int_t i, Int_t j) const
{
  /// Return the value as an int, by rounding the internal float value.
  
  Float_t v = ValueAsFloatFast(i,j);
  return TMath::Nint(v);
}
