// @(#) $Id$

#ifndef AliHLT_Compress
#define AliHLT_Compress

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"
#include "AliHLTDataCompressor.h"

class AliHLTCompress {
  
 public:
  AliHLTCompress();
  AliHLTCompress(Int_t slice,Int_t patch,Char_t *path="./",Bool_t writeshape=kFALSE,Int_t event=-1);
  virtual ~AliHLTCompress();
  
  Bool_t WriteFile(AliHLTTrackArray *tracks,Char_t *filename=0);
  Bool_t ReadFile(Char_t which,Char_t *filename=0);
  virtual Bool_t CompressFile();
  virtual Bool_t ExpandFile();
  void CompressRemaining(AliHLTSpacePointData *points[36][6],UInt_t npoints[36][6]);
  void ExpandRemaining(TempCluster **clusters,Int_t *ncl,Int_t maxclusters);
  virtual void PrintCompRatio(STDOF *outfile=0);
  Int_t GetEntropy(Float_t &padEntropy,Float_t &timeEntropy,Float_t &chargeEntropy);
  
  AliHLTTrackArray *GetTracks() {return fTracks;}
  
 protected:
  AliHLTTrackArray *fTracks; //! Array of tracks
  Int_t fSlice; // Slice
  Int_t fPatch; // Patch
  Char_t fPath[100]; // Path to the files
  Bool_t fWriteShape; // Flag to write the shape
  Int_t fEvent; // Current event

  
  ClassDef(AliHLTCompress,1) 

};

typedef AliHLTCompress AliL3Compress; // for backward compatibility

#endif
