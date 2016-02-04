#ifndef ALIPHOSCPVBADCHANNELSMAP
#define ALIPHOSCPVBADCHANNELSMAP
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// This class keeps the CPV bad channels map 
// (bad means dead or noisy).

#include "TObject.h"
#include "AliPHOSCpvParam.h"

class AliPHOSCpvBadChannelsMap : public TObject {

public:

  AliPHOSCpvBadChannelsMap();
  AliPHOSCpvBadChannelsMap(const AliPHOSCpvBadChannelsMap &map);
  AliPHOSCpvBadChannelsMap& operator= (const AliPHOSCpvBadChannelsMap &map);
  ~AliPHOSCpvBadChannelsMap() {}

  void  SetBadChannel(Int_t module, Int_t col, Int_t row);
  Bool_t IsBadChannel(Int_t module, Int_t col, Int_t row) const { return fBadChannelCpv[module-1][col-1][row-1]; }
  Int_t GetNumOfBadChannels() const {  return fBads; }
  void BadChannelIds(Int_t *badIds=0);
  void Reset();//reset all channels as good
  void Reset(Int_t module);//reset all channels in module as good

private:
  
  Bool_t fBadChannelCpv[AliPHOSCpvParam::kNModules][AliPHOSCpvParam::kPadPcX][AliPHOSCpvParam::kPadPcY]; //[mod][col][row]
  Int_t fBads;

  ClassDef(AliPHOSCpvBadChannelsMap,1)

};

#endif
