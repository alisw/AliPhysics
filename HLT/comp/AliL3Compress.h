// @(#) $Id$

#ifndef AliL3_Compress
#define AliL3_Compress

#include "AliL3RootTypes.h"
#include "AliL3DataCompressor.h"

class AliL3Compress {
  
 private:
  
 protected:
  AliL3TrackArray *fTracks; //!
  Int_t fSlice;
  Int_t fPatch;
  Char_t fPath[100];
  Bool_t fWriteShape;
  Int_t fEvent;

  
 public:
  AliL3Compress();
  AliL3Compress(Int_t slice,Int_t patch,Char_t *path="./",Bool_t writeshape=kFALSE,Int_t event=-1);
  virtual ~AliL3Compress();
  
  Bool_t WriteFile(AliL3TrackArray *tracks,Char_t *filename=0);
  Bool_t ReadFile(Char_t which,Char_t *filename=0);
  virtual Bool_t CompressFile();
  virtual Bool_t ExpandFile();
  void CompressRemaining(AliL3SpacePointData *points[36][6],UInt_t npoints[36][6]);
  void ExpandRemaining(TempCluster **clusters,Int_t *ncl,Int_t maxclusters);
  virtual void PrintCompRatio(STDOF *outfile=0);
  Int_t GetEntropy(Float_t &pad_entropy,Float_t &time_entropy,Float_t &charge_entropy);
  
  AliL3TrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliL3Compress,1) 

};

#endif
