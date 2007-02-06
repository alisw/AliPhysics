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

#include "AliMUONCalibParam1I.h"

#include "AliLog.h"
#include "Riostream.h"
#include "TMath.h"
#include "TString.h"

/// \class AliMUONCalibParam1I
///
/// This class is implementing the AliMUONVCalibParam interface.
/// 
/// It stores a given number of integers.
/// 
/// Those integers can also be retrieved as floats if really needed 
/// (this is to comply with the base class).
/// 
/// You might consider just as it is, namely a C++ wrapper to 
/// a plain int[] array.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONCalibParam1I)
/// \endcond

//_____________________________________________________________________________
AliMUONCalibParam1I::AliMUONCalibParam1I() 
: AliMUONVCalibParam(),
  fSize(0),
  fValues(0x0)
{
/// Default constructor.
}

//_____________________________________________________________________________
AliMUONCalibParam1I::AliMUONCalibParam1I(Int_t theSize, Int_t fillWithValue) 
: AliMUONVCalibParam(),
  fSize(theSize),
  fValues(0x0)
{
/// Normal constructor, where theSize specifies the number of channels handled
/// by this object, and fillWithValue the default value assigned to each
/// channel.

  if ( fSize > 0 )
  {
    fValues = new Int_t[fSize];
    for ( Int_t i = 0; i < fSize; ++i )
    {
      fValues[i] = fillWithValue;
    }
  }
}

//_____________________________________________________________________________
AliMUONCalibParam1I::AliMUONCalibParam1I(const AliMUONCalibParam1I& other) 
: AliMUONVCalibParam(other),
fSize(0),
fValues(0x0)
{
/// Copy constructor

  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUONCalibParam1I&
AliMUONCalibParam1I::operator=(const AliMUONCalibParam1I& other) 
{
/// Assignment operator

  AliMUONVCalibParam::operator=(other);
  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONCalibParam1I::~AliMUONCalibParam1I()
{
/// Destructor.

  delete[] fValues;
}

//_____________________________________________________________________________
void
AliMUONCalibParam1I::CopyTo(AliMUONCalibParam1I& destination) const
{
/// Copy this into destination.

  delete[] destination.fValues;
  destination.fSize = fSize;
  if ( fSize > 0 )
  {
    destination.fValues = new Int_t[fSize];
    for ( Int_t i = 0; i < fSize; ++i )
    {
      destination.fValues[i] = fValues[i];
    }
  }
}

//_____________________________________________________________________________
void
AliMUONCalibParam1I::Print(Option_t* opt) const
{
/// Output this object to stdout.
/// If opt=="full", then all channels are printed, 
/// if opt=="mean", only the mean and sigma value are printed,
/// otherwise only the general characteristics are printed.

  TString sopt(opt);
  sopt.ToUpper();
  cout << "AliMUONCalibParam1I - Size=" << Size()
    << " Dimension=" << Dimension();
  
  if ( sopt.Contains("FULL") )
  {
    cout << endl;
    for ( Int_t i = 0; i < Size(); ++i )
    {
      cout << Form("CH %3d %10d",i,ValueAsInt(i)) << endl;
    }
  }
  if ( sopt.Contains("MEAN") )
  {
    Int_t mean(0);
    Int_t v2(0);
    
    Int_t N = Size();
    
    for ( Int_t i = 0; i < Size(); ++i )
    {
      Int_t v = ValueAsInt(i);
      mean += v;
      v2 += v*v;
    }
    mean /= N;
    float sigma = 0;
    if ( N > 1 ) sigma = TMath::Sqrt( (v2-N*mean*mean)/(N-1) );
    cout << Form(" Mean=%d Sigma=%f",mean,sigma) << endl;
  }
}

//_____________________________________________________________________________
void
AliMUONCalibParam1I::SetValueAsFloat(Int_t i, Int_t j, Float_t value)
{
/// Set the value as a float, which is casted to an int prior to storage.

  SetValueAsInt(i,j,TMath::Nint(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParam1I::SetValueAsInt(Int_t i, Int_t j, Int_t value)
{
/// Set the value for a given channel.
/// (i,j) are checked for correctness before use.

  if ( j != 0 || i >= fSize || i < 0 )
  {
    AliError(Form("Invalid (i,j)=(%d,%d) max allowed is (%d,%d)",
                  i,j,Size()-1,Dimension()-1));
  }
  else
  {
    fValues[i]=value;
  }
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParam1I::ValueAsFloat(Int_t i, Int_t j) const
{
/// Return one value as a float.

  return 1.0*ValueAsInt(i,j);
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParam1I::ValueAsInt(Int_t i, Int_t j) const
{
/// Return one value as an integer, after checking that (i,j)
/// are valid indices.

  if ( j != 0 || i >= fSize || i < 0 )
  {
    AliError(Form("Invalid (i,j)=(%d,%d) max allowed is (%d,%d)",
                  i,j,Size()-1,Dimension()-1));
    return 0;
  }
  else
  {
    return fValues[i];
  }
}
