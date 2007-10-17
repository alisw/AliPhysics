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

#include "AliMUONVCalibParam.h"

#include "AliLog.h"

//-----------------------------------------------------------------------------
/// \class AliMUONVCalibParam
///  
/// Defines an interface for a calibration container object.
///
/// Note that a VCalibParam object is identified by a pair (id0,id1), 
/// where each member of the pair is a 16 bits word, so that id0 and id1
/// can be packed into a single 32 bits word.
///
/// id1 might be left to zero if not required (e.g. for calibparam which 
/// can be identified by a single integer)
///
/// Note that the ValueAsXXX methods have 2 versions : with or without bound
/// checking. The latter is to be used in e.g. loops, where you know for
/// sure the indices are ok, in order to gain some time.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONVCalibParam)
/// \endcond

//_____________________________________________________________________________
AliMUONVCalibParam::AliMUONVCalibParam() : TObject()
{
  ///  Default constructor
}

//_____________________________________________________________________________
AliMUONVCalibParam::AliMUONVCalibParam(Int_t id0, Int_t id1) : TObject()
{
  ///  constructor for 2D
  SetUniqueID(BuildUniqueID(id0,id1));
}

//_____________________________________________________________________________
AliMUONVCalibParam::~AliMUONVCalibParam()
{
/// Destructor.
}

//_____________________________________________________________________________
UInt_t
AliMUONVCalibParam::BuildUniqueID(Int_t id0, Int_t id1)
{
  /// Build a single index from the pair (id0,id1)
  return ( id0 | ( id1 << 16 ) );
}

//_____________________________________________________________________________
void
AliMUONVCalibParam::DecodeUniqueID(UInt_t uniqueID, Int_t& id0, Int_t& id1)
{
  /// Convert single integer into a pair (i,j)
  id0 = ID0(uniqueID);
  id1 = ID1(uniqueID);
}

//_____________________________________________________________________________
Int_t
AliMUONVCalibParam::ID0(UInt_t uniqueID)
{
  /// Extract id0 from uniqueID
  return uniqueID & 0xFFFF;
}

//_____________________________________________________________________________
Int_t
AliMUONVCalibParam::ID0() const
{
  /// Extract first identifier
  return ID0(GetUniqueID());
}

//_____________________________________________________________________________
Int_t
AliMUONVCalibParam::ID1(UInt_t uniqueID)
{
  /// Extract id1 from uniqueID
  return ( uniqueID & 0xFFFF0000 ) >> 16;
}

//_____________________________________________________________________________
Int_t
AliMUONVCalibParam::ID1() const
{
  /// Extract second identifier
  return ID1(GetUniqueID());
}

//_____________________________________________________________________________
const char* 
AliMUONVCalibParam::GetName() const
{
  /// Build a name for this object
  return Form("I=%d,J=%d",ID0(),ID1());
}

//_____________________________________________________________________________
void 
AliMUONVCalibParam::SetValueAsDouble(Int_t, Int_t, Double_t)
{
  /// By default, this one does not exist
  AliFatal("Not implemented");
}

//_____________________________________________________________________________
void 
AliMUONVCalibParam::SetValueAsDoubleFast(Int_t, Int_t, Double_t)
{
  /// By default, this one does not exist
  AliFatal("Not implemented");
}

//_____________________________________________________________________________
Double_t 
AliMUONVCalibParam::ValueAsDouble(Int_t , Int_t ) const
{
  /// By default, this one does not exist
  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
Double_t 
AliMUONVCalibParam::ValueAsDoubleFast(Int_t , Int_t ) const
{
  /// By default, this one does not exist
  AliFatal("Not implemented");
  return 0;
}

