// @(#) $Id$

#ifndef ALIL3OFFLINEDATACOMPRESSOR_H
#define ALIL3OFFLINEDATACOMPRESSOR_H

#include "AliL3RootTypes.h"
#include "AliL3DataCompressor.h"

class AliTracker;

class AliL3OfflineDataCompressor : public AliL3DataCompressor {
  
 public:
  AliL3OfflineDataCompressor();
  AliL3OfflineDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape,Bool_t MI=kFALSE);
  virtual ~AliL3OfflineDataCompressor();
  
  void LoadData(Int_t event,Bool_t sp=kTRUE);
  void FillData(Int_t /*minhits*/,Bool_t /*expand*/) {return;};
  void WriteRemaining(Bool_t select);

 private:
  AliL3OfflineDataCompressor(const AliL3OfflineDataCompressor& /*ac*/) : AliL3DataCompressor() {;}
  AliL3OfflineDataCompressor& operator=(const AliL3OfflineDataCompressor& /*ac*/){return *this;}

  Bool_t fMarian;        // is Marian TPC tracking used
  AliTracker *fTracker;  //!
  
  void SelectRemainingClusters();

  ClassDef(AliL3OfflineDataCompressor,1) 

};

#endif
