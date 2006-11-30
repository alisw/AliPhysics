// @(#) $Id$

#ifndef ALIL3OFFLINEDATACOMPRESSOR_H
#define ALIL3OFFLINEDATACOMPRESSOR_H

#include "AliHLTRootTypes.h"
#include "AliHLTDataCompressor.h"

class AliTracker;

class AliHLTOfflineDataCompressor : public AliHLTDataCompressor {
  
 public:
  AliHLTOfflineDataCompressor();
  AliHLTOfflineDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape,Bool_t MI=kFALSE);
  virtual ~AliHLTOfflineDataCompressor();
  
  void LoadData(Int_t event,Bool_t sp=kTRUE);
  void FillData(Int_t /*minhits*/,Bool_t /*expand*/) {return;};
  void WriteRemaining(Bool_t select);

 private:
  AliHLTOfflineDataCompressor(const AliHLTOfflineDataCompressor& /*ac*/) : AliHLTDataCompressor() {;}
  AliHLTOfflineDataCompressor& operator=(const AliHLTOfflineDataCompressor& /*ac*/){return *this;}

  Bool_t fMarian;        // is Marian TPC tracking used
  AliTracker *fTracker;  //!
  
  void SelectRemainingClusters();

  ClassDef(AliHLTOfflineDataCompressor,1) 

};

typedef AliHLTOfflineDataCompressor AliL3OfflineDataCompressor; // for backward compatibility

#endif
