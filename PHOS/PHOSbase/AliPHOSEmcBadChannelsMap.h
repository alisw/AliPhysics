#ifndef ALIPHOSEMCBADCHANNELSMAP
#define ALIPHOSEMCBADCHANNELSMAP
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// This class keeps the EMC bad channels map 
// (bad means dead or noisy).

#include "TObject.h"

class AliPHOSEmcBadChannelsMap : public TObject {

public:

  AliPHOSEmcBadChannelsMap();
  AliPHOSEmcBadChannelsMap(const AliPHOSEmcBadChannelsMap &map);
  AliPHOSEmcBadChannelsMap& operator= (const AliPHOSEmcBadChannelsMap &map);
  ~AliPHOSEmcBadChannelsMap() {}

  void  SetBadChannel(Int_t module, Int_t col, Int_t row);
  Bool_t IsBadChannel(Int_t module, Int_t col, Int_t row) const { return fBadChannelEmc[module-1][col-1][row-1]; }
  Int_t GetNumOfBadChannels() const {  return fBads; }
  void BadChannelIds(Int_t *badIds=0);
  void Reset();

private:
  
  Bool_t fBadChannelEmc[5][56][64]; //[mod][col][row]
  Int_t fBads;

  ClassDef(AliPHOSEmcBadChannelsMap,2)

};

#endif
