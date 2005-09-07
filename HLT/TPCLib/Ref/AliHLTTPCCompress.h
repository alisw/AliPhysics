// @(#) $Id$

#ifndef AliHLTTPC_Compress
#define AliHLTTPC_Compress

class AliHLTTPCTrackArray;

class AliHLTTPCCompress {
  
 private:
  AliHLTTPCTrackArray *fTracks; //!

  Int_t fSlice;
  Int_t fPatch;
  Char_t fPath[100];
  Bool_t fWriteShape;
  Int_t fEvent;

 public:
  AliHLTTPCCompress();
  AliHLTTPCCompress(Int_t slice,Int_t patch,Char_t *path="./",Bool_t writeshape=kFALSE,Int_t event=-1);
  virtual ~AliHLTTPCCompress();
  
  Bool_t WriteFile(AliHLTTPCTrackArray *tracks,Char_t *filename=0);
  Bool_t ReadFile(Char_t which,Char_t *filename=0);
  Bool_t CompressFile();
  Bool_t ExpandFile();
  void PrintCompRatio(ofstream *outfile=0);
  
  AliHLTTPCTrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliHLTTPCCompress,1) 

};

#endif
