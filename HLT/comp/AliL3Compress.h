#ifndef AliL3_Compress
#define AliL3_Compress

#include "AliL3RootTypes.h"
#include "AliL3Models.h"

class AliL3TrackArray;

class AliL3Compress {
  
 private:
  
  
  
 public:
  AliL3Compress();
  virtual ~AliL3Compress();
  
  void Write2File(AliL3TrackArray *tracks);
  void ReadFile();
  void CompressFile();
  void ExpandFile();

  ClassDef(AliL3Compress,1) 

};

#endif
