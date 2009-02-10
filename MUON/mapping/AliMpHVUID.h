#ifndef ALIMPHVUID_H
#define ALIMPHVUID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpHVUID
/// \brief Unique ID class for HV channels
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMpHVUID : public TObject
{
public:
  AliMpHVUID();
  virtual ~AliMpHVUID();
  
  static UInt_t BuildUniqueID(Int_t detElemId, Int_t index);
  
  static Int_t DetElemId(UInt_t uniqueId);
  
  static Int_t Index(UInt_t uniqueID);
                     
  ClassDef(AliMpHVUID,1) // Unique ID class for HV channels
};

#endif
