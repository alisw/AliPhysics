// @(#) $Id$

#ifndef AliL3_Trigger
#define AliL3_Trigger

#include "AliL3RootTypes.h"

class AliL3TrackSegmentData;
class AliL3DigitRowData;
class AliL3TrackArray;
class AliL3Vertex;

class AliL3Trigger {
 
 private:
  AliL3TrackArray *fTracks; //!
  AliL3DigitRowData *fDigitRowData; //!
  AliL3DigitRowData *fOutput; //!
  AliL3Vertex *fVertex; //!
  Int_t fDataSize;

  Float_t fZcut;
  Int_t fTimeMatch;
  Int_t fPadMatch;
  Int_t fSlice;
  Int_t fPatch;

 public:
  AliL3Trigger();
  virtual ~AliL3Trigger();
  
  void InitTrigger();
  void InitPatch(Int_t slice,Int_t patch);
  void FillTracks(Int_t ntracks,AliL3TrackSegmentData *tr);
  void FillData(AliL3DigitRowData *data);
  void SetOutputData(AliL3DigitRowData *ptr);
  void SetVertex(AliL3Vertex *vertex) {fVertex = vertex;}
  void SetParameters(Float_t zcut,Int_t timematch,Int_t padmatch);
  void RemovePileupTracks();
  void RemovePileupData();
  
  Int_t GetDataSize() {return fDataSize;}
  
  ClassDef(AliL3Trigger,1) 

};

#endif
