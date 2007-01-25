/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpCathodType.h,v 1.8 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \enum AliMpCathodType
/// Enumeration for refering to cath0 and cath1.
///
/// \author Ivana Hrivnacova; IPN Orsay
 
#include "AliMpCathodType.h"

#include "AliLog.h"

//_____________________________________________________________________________
AliMp::CathodType AliMp::GetCathodType(Int_t cathodNumber)
{
  switch ( cathodNumber ) {
    case kCath0:  return kCath0;  break;
    case kCath1:  return kCath1;  break;
    default:  
      // Should reach this line
      AliErrorGeneral("AliMpCathodType.h", "Wrong cathod number"); 
      return kCath0;
  }
  
  // Should reach this line
  AliErrorGeneral("AliMpCathodType.h", "Wrong cathod number"); 
  return kCath0;
}       

//_____________________________________________________________________________
TString AliMp::CathodTypeName(AliMp::CathodType cathodType)
{
  switch ( cathodType ) {
    case kCath0:  return "cath0";  break;
    case kCath1:  return "cath1"; break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpCathodType.h", "Unknown cathod type"); 
  return "invalidCathod";
}       

//_____________________________________________________________________________
AliMp::CathodType AliMp::OtherCathodType(AliMp::CathodType cathodType)
{
  switch ( cathodType ) {
    case kCath0: return kCath1;  break;
    case kCath1: return kCath0;  break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpCathodType.h", "Unknown cathod type"); 
  return kCath0;
}       
