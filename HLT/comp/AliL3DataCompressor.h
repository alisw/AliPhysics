#ifndef AliL3_DataCompressor
#define AliL3_DataCompressor

#include "AliL3RootTypes.h"

class AliL3MemHandler;
class AliL3Benchmark;

class AliL3DataCompressor {
  
 private:
  AliL3MemHandler *fMemHandler;  //!
  AliL3Benchmark *fBenchmark;    //!
  Int_t fMinSlice;
  Int_t fMaxSlice;
  Char_t fPath[1024]; //!
  
  Int_t FindRemaining(Int_t slice,Int_t patch);

 public:
  AliL3DataCompressor();
  AliL3DataCompressor(Char_t *path,Int_t minslice,Int_t maxslice);
  virtual ~AliL3DataCompressor();
  
  void ProcessData(Char_t *trackpath,Int_t padoverlap,Int_t timeoverlap,Int_t padsearch,Int_t timesearch);
  void CompressAndExpand(Int_t bitspad,Int_t bitstime,Int_t bitscharge,Int_t bitsshape);
  void WriteRemainingDigits();
  void FindOfflineClusters(Bool_t remains);
  void RestoreData();
  void DoBench(Char_t *fname="benchmark");

  ClassDef(AliL3DataCompressor,1) 

};

#endif
