/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

//-----------------------------------------------------------------------------
// Enum AliMq::Station12Type
// Enumeration for refering to a MUON station12 type.
//
// Author: Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------
 
#include "AliMpStation12Type.h"

#include "AliLog.h" 

//_____________________________________________________________________________
TString AliMq::Station12TypeName(AliMq::Station12Type station12Type)
{
  switch ( station12Type ) {
    case kStation1:         return "st1";     break;
    case kStation2:         return "st2";     break;
    case kNotSt12:          return "";        break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpStation12Type.h", "Unknown sector station type"); 
  return "invalidStation12";
}
