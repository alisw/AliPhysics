#ifndef AliL3_Compress
#define AliL3_Compress

#include "AliL3RootTypes.h"
#include "AliL3Models.h"

class AliL3TrackArray;

class AliL3Compress {
  
 private:
  AliL3TrackArray *fTracks; //!
  
  
 public:
  AliL3Compress();
  virtual ~AliL3Compress();
  
  void WriteFile(AliL3TrackArray *tracks,Char_t *filename);
  void ReadFile(Char_t *filename);
  void CompressFile(Char_t *infile,Char_t *outfile);
  void ExpandFile(Char_t *infile,Char_t *outfile);
  
  AliL3TrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliL3Compress,1) 

};

#endif
