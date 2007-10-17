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

#include "AliMUONCalibParamNI.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TMath.h"
#include "TString.h"

//-----------------------------------------------------------------------------
/// \class AliMUONCalibParamNI
///
/// Handle the case of N integer parameters per channel.
///
/// Almost the same class as AliMUONCalibParamNF, but for ints.
/// We could have played with NF to store both int and float (using casts),
/// but for the sake of simplicity (e.g. the Print method must know whether
/// it should print floats or ints), we decided to "duplicate" the class
/// and use the correct type.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONCalibParamNI)
/// \endcond

//_____________________________________________________________________________
AliMUONCalibParamNI::AliMUONCalibParamNI() 
: AliMUONVCalibParam(),
  fDimension(0),
  fSize(0),
  fN(0),
  fPackingFactor(0),
  fValues(0x0)
{
/// Default constructor.
    AliDebug(1,Form("this=%p default ctor",this));
}

//_____________________________________________________________________________
AliMUONCalibParamNI::AliMUONCalibParamNI(Int_t dimension, Int_t theSize,
                                         Int_t id0, Int_t id1,
                                         Int_t fillWithValue, Int_t packingFactor) 
: AliMUONVCalibParam(id0,id1),
  fDimension(dimension),
  fSize(theSize),
  fN(fSize*fDimension),
  fPackingFactor(packingFactor),
  fValues(0x0)
{
/// Normal constructor, where theSize specifies the number of channels handled
/// by this object, and fillWithValue the default value assigned to each
/// channel.

//    AliDebug(1,Form("this=%p dimension=%d theSize=%d fillWithValue=%d packingFactor=%d",
//                    this,dimension,theSize,fillWithValue,packingFactor));

  if ( fN > 0 )
  {
    fValues = new Int_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      fValues[i] = fillWithValue;
    }
  }
}


//_____________________________________________________________________________
AliMUONCalibParamNI::AliMUONCalibParamNI(const AliMUONCalibParamNI& other) 
: AliMUONVCalibParam(),
fDimension(0),
fSize(0),
fN(0),
fPackingFactor(0),
fValues(0x0)
{
/// Copy constructor.

  AliDebug(1,Form("this=%p copy ctor",this));
  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUONCalibParamNI&
AliMUONCalibParamNI::operator=(const AliMUONCalibParamNI& other) 
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONCalibParamNI::~AliMUONCalibParamNI()
{
/// Destructor

  AliDebug(1,Form("this=%p",this));
  delete[] fValues;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNI::CopyTo(AliMUONCalibParamNI& destination) const
{
/// Copy *this to destination

  const TObject& o = static_cast<const TObject&>(*this);
  o.Copy(destination);
  
  delete[] destination.fValues;
  destination.fN = fN;
  destination.fSize = fSize;
  destination.fDimension = fDimension;
  destination.fPackingFactor = fPackingFactor;
  
  if ( fN > 0 )
  {
    destination.fValues = new Int_t[fN];
    for ( Int_t i = 0; i < fN; ++i )
    {
      destination.fValues[i] = fValues[i];
    }
  }
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNI::Index(Int_t i, Int_t j) const
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
AliMUONCalibParamNI::IndexFast(Int_t i, Int_t j) const
{
  /// Compute the 1D index of the internal storage from the pair (i,j)
  
  return i + Size()*j;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNI::Print(Option_t* opt) const
{
/// Output this object to stdout.
/// If opt=="full" then all channels are printed, 
/// if opt=="mean#", only the mean and sigma value are printed for j-th dimension
/// otherwise only the general characteristics are printed.

  TString sopt(opt);
  sopt.ToUpper();
  cout << "AliMUONCalibParamNI - Size=" << Size()
    << " Dimension=" << Dimension()
    << " Id=(" << ID0() << "," << ID1() << ")";
  
  if ( IsPacked() ) 
  {
    cout << " Packing Factor=" << fPackingFactor;
  }
  cout << endl;
  
  if ( sopt.Contains("FULL") )
  {
    for ( Int_t i = 0; i < Size(); ++i )
    {
      cout << Form("CH %3d",i);
      for ( Int_t j = 0; j < Dimension(); ++j )
      {
        Int_t v = ValueAsInt(i,j);
        if ( IsPacked() )
        {
          Int_t m,c;
          UnpackValue(v,m,c);
          cout << Form(" (%d,%d) (0x%x,0x%x)",m,c,m,c);
        }
        else
        {
          cout << Form(" %d (0x%x)",v,v);
        }
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
AliMUONCalibParamNI::SetValueAsInt(Int_t i, Int_t j, Int_t value)
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
AliMUONCalibParamNI::SetValueAsIntFast(Int_t i, Int_t j, Int_t value)
{
  /// Set one value as a float, w/o checking that the indices are correct.
  
  fValues[IndexFast(i,j)] = value;
}

//_____________________________________________________________________________
void
AliMUONCalibParamNI::SetValueAsFloat(Int_t i, Int_t j, Float_t value)
{
  /// Set one value as an int.

  AliWarning("Float will be rounded to be stored...");
  SetValueAsInt(i,j,TMath::Nint(value));
}

//_____________________________________________________________________________
void
AliMUONCalibParamNI::SetValueAsFloatFast(Int_t i, Int_t j, Float_t value)
{
  /// Set one value as an int.
  
  AliWarning("Float will be rounded to be stored...");
  SetValueAsIntFast(i,j,TMath::Nint(value));
}


//_____________________________________________________________________________
Int_t
AliMUONCalibParamNI::ValueAsInt(Int_t i, Int_t j) const
{
/// Return the value as an int (which it is), after checking indices.

  Int_t ix = Index(i,j);
  
  if ( ix < 0 )
  {
    AliError(Form("Invalid (i,j)=(%d,%d) max allowed is (%d,%d)",
                  i,j,Size()-1,Dimension()-1));
    return 0;
  }
  else
  {
    return fValues[ix];
  }
}

//_____________________________________________________________________________
Int_t
AliMUONCalibParamNI::ValueAsIntFast(Int_t i, Int_t j) const
{
  /// Return the value as an int (which it is), w/o checking indices.
  
  return fValues[IndexFast(i,j)];
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamNI::ValueAsFloat(Int_t i, Int_t j) const
{
  /// Return the value as a float
  return static_cast<Float_t>(ValueAsInt(i,j));
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParamNI::ValueAsFloatFast(Int_t i, Int_t j) const
{
  /// Return the value as a float
  return static_cast<Float_t>(ValueAsIntFast(i,j));
}

//_____________________________________________________________________________
Bool_t 
AliMUONCalibParamNI::UnpackValue(Int_t value, Int_t& a, Int_t& b) const
{
  /// Unpack single value into a couple (a,b), using packingFactor
  /// Returns false if IsPacked()==false
  
  if ( !IsPacked() ) return kFALSE;
  a = value / fPackingFactor;
  b = value % fPackingFactor;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONCalibParamNI::PackValues(Int_t a, Int_t b, Int_t& packedValue) const
{
  /// Pack a couple (a,b) into a single value, using packingFactor
  /// Returns false if IsPacked()==false

  if ( !IsPacked() ) return kFALSE;
  packedValue = a*fPackingFactor + b;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONCalibParamNI::IsPacked() const
{
  /// Whether the values we store are packed or not.
  /// If false, Pack and Unpack methods will always return false and do nothing
  return (fPackingFactor != 0);
}

