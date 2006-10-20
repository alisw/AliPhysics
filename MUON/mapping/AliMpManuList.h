#ifndef ALIMPMANULIST_H
#define ALIMPMANULIST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpManuList
/// \brief Cache of often used information
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TList;

class AliMpManuList : public TObject
{
public:
  AliMpManuList();
  virtual ~AliMpManuList();
  
  /// return a list of AliMpIntPair(detElemId, manuID). Returned value must be deleted.
  static TList* ManuList();

  /// whether a given (de,id,ch) exists
  static Bool_t DoesChannelExist(Int_t detElemId, Int_t manuID, Int_t manuChannel);
  
  /// number of manu in a given DE
  static Int_t NumberOfManus(Int_t detElemId);
  
  /// number of channels in a given manu (<=64)
  static Int_t NumberOfChannels(Int_t detElemId, Int_t manuId);
  
  ClassDef(AliMpManuList,1) // 
};

#endif
