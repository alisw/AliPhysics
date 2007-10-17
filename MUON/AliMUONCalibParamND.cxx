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

#include "AliMUONCalibParamND.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TMath.h"
#include "TString.h"

//-----------------------------------------------------------------------------
/// \class AliMUONCalibParamND
///
/// Handle the case of N floating point (double precision) parameters per channel.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONCalibParamND)
/// \endcond

//_____________________________________________________________________________
AliMUONCalibParamND::AliMUONCalibParamND() 
: AliMUONVCalibParam(),
  fDimension(0),
  fSize(0),
  fN(0),
  fValues(0x0)
{
/// Default constructor.
}

//_____________________________________________________________________________
AliMUONCalibParamND::AliMUONCalibParamND(Int_t dimension, Int_t theSize, 
                                         Int_t id0, Int_t id1,
                                         Double_t fillWithValue) 
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
    fValues = new Double_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      fValues[i] = fillWithValue;
    }
  }
}


//_____________________________________________________________________________
AliMUONCalibParamND::AliMUONCalibParamND(const AliMUONCalibParamND& other) 
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
AliMUONCalibParamND&
AliMUONCalibParamND::operator=(const AliMUONCalibParamND& other) 
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONCalibParamND::~AliMUONCalibParamND()
{
/// Destructor

  delete[] fValues;
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::CopyTo(AliMUONCalibParamND& destination) const
{
/// Copy *this to destination

  TObject::Copy(destination);
  
  delete[] destination.fValues;
  destination.fN = fN;
  destination.fSize = fSize;
  destination.fDimension = fDimension;

  if ( fN > 0 )
  {
    destination.fValues = new Double_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      destination.fValues[i] = fValues[i];
    }
  }
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamND::Index(Int_t i, Int_t j) const
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
AliMUONCalibParamND::IndexFast(Int_t i, Int_t j) const
{
  /// Compute the 1D index of the internal storage from the pair (i,j)
  
  return i + Size()*j;
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::Print(Option_t* opt) const
{
/// Output this object to stdout.
/// If opt=="full", then all channels are printed, 
/// if opt=="mean#", only the mean and sigma value are printed for j-th dimension
/// otherwise only the general characteristics are printed.

  TString sopt(opt);
  sopt.ToUpper();
  cout << Form("AliMUONCalibParamND Id=(%d,%d) Size=%d Dimension=%d",ID0(),
               ID1(),Size(),Dimension()) << endl;
  if ( sopt.Contains("FULL") )
  {
    for ( Int_t i = 0; i < Size(); ++i )
    {
      cout << Form("CH %3d",i);
      for ( Int_t j = 0; j < Dimension(); ++j )
      {
        cout << Form(" %g",ValueAsDouble(i,j));
      }
      cout << endl;
    }
  }  
  if ( sopt.Contains("MEAN") )
  {
    Int_t j(0);
    sscanf(sopt.Data(),"MEAN%d",&j);
    
    Double_t mean(0);
    Double_t v2(0);
    
    Int_t n = Size();
    
    for ( Int_t i = 0; i < Size(); ++i )
    {
      Float_t v = ValueAsDouble(i,j);
      mean += v;
      v2 += v*v;
    }
    mean /= n;
    float sigma = 0;
    if ( n > 1 ) sigma = TMath::Sqrt( (v2-n*mean*mean)/(n-1) );
    cout << Form(" Mean(j=%d)=%g Sigma(j=%d)=%g",j,mean,j,sigma) << endl;
  }
  
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::SetValueAsDouble(Int_t i, Int_t j, Double_t value)
{
  /// Set one value as a double, after checking that the indices are correct.
  
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
AliMUONCalibParamND::SetValueAsDoubleFast(Int_t i, Int_t j, Double_t value)
{
  /// Set one value as a double, w/o checking that the indices are correct.
  
  fValues[IndexFast(i,j)] = value;
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::SetValueAsFloat(Int_t i, Int_t j, Float_t value)
{
  /// Set one value as a float, after checking that the indices are correct.
  SetValueAsDouble(i,j,static_cast<Double_t>(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::SetValueAsFloatFast(Int_t i, Int_t j, Float_t value)
{
  /// Set one value as a float, after checking that the indices are correct.
  SetValueAsDoubleFast(i,j,static_cast<Double_t>(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::SetValueAsInt(Int_t i, Int_t j, Int_t value)
{
/// Set one value as an int.

  SetValueAsFloat(i,j,static_cast<Float_t>(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParamND::SetValueAsIntFast(Int_t i, Int_t j, Int_t value)
{
  /// Set one value as an int.
  
  SetValueAsFloatFast(i,j,static_cast<Float_t>(value));
}

//_____________________________________________________________________________
Double_t
AliMUONCalibParamND::ValueAsDouble(Int_t i, Int_t j) const
{
  /// Return the value as a double (which it is), after checking indices.
  
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
Double_t
AliMUONCalibParamND::ValueAsDoubleFast(Int_t i, Int_t j) const
{
  /// Return the value as a double (which it is), w/o checking indices.
  
  return fValues[IndexFast(i,j)];
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamND::ValueAsFloat(Int_t i, Int_t j) const
{
/// Return the value as a float 
  return static_cast<Float_t>(ValueAsDouble(i,j));
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamND::ValueAsFloatFast(Int_t i, Int_t j) const
{
  /// Return the value as a float 
  return static_cast<Float_t>(ValueAsDoubleFast(i,j));
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamND::ValueAsInt(Int_t i, Int_t j) const
{
/// Return the value as an int, by rounding the internal float value.

  Float_t v = ValueAsFloat(i,j);
  return TMath::Nint(v);
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamND::ValueAsIntFast(Int_t i, Int_t j) const
{
  /// Return the value as an int, by rounding the internal float value.
  
  Float_t v = ValueAsFloatFast(i,j);
  return TMath::Nint(v);
}
