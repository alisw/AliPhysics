// @(#) $Id$

#ifndef AliL3_OfflineDataCompressor
#define AliL3_OfflineDataCompressor

#include "AliL3RootTypes.h"
#include "AliL3DataCompressor.h"

class AliTracker;

class AliL3OfflineDataCompressor : public AliL3DataCompressor {
  
 private:
  Bool_t fMarian;
  AliTracker *fTracker;  //!
  
  void SelectRemainingClusters();
  
 public:
  AliL3OfflineDataCompressor();
  AliL3OfflineDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape,Bool_t MI=kFALSE);
  virtual ~AliL3OfflineDataCompressor();
  
  void LoadData(Int_t event,Bool_t sp=kTRUE);
  void FillData(Int_t /*minhits*/,Bool_t /*expand*/) {return;};
  void WriteRemaining(Bool_t select);

  ClassDef(AliL3OfflineDataCompressor,1) 

};

#endif
