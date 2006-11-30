// @(#) $Id$

#ifndef AliHLT_Trigger
#define AliHLT_Trigger

#include "AliHLTRootTypes.h"

class AliHLTTrackSegmentData;
class AliHLTDigitRowData;
class AliHLTTrackArray;
class AliHLTVertex;

class AliHLTTrigger {
 
 private:
  AliHLTTrackArray *fTracks; //!
  AliHLTDigitRowData *fDigitRowData; //!
  AliHLTDigitRowData *fOutput; //!
  AliHLTVertex *fVertex; //!
  Int_t fDataSize;

  Float_t fZcut;
  Int_t fTimeMatch;
  Int_t fPadMatch;
  Int_t fSlice;
  Int_t fPatch;

 public:
  AliHLTTrigger();
  virtual ~AliHLTTrigger();
  
  void InitTrigger();
  void InitPatch(Int_t slice,Int_t patch);
  void FillTracks(Int_t ntracks,AliHLTTrackSegmentData *tr);
  void FillData(AliHLTDigitRowData *data);
  void SetOutputData(AliHLTDigitRowData *ptr);
  void SetVertex(AliHLTVertex *vertex) {fVertex = vertex;}
  void SetParameters(Float_t zcut,Int_t timematch,Int_t padmatch);
  void RemovePileupTracks();
  void RemovePileupData();
  
  Int_t GetDataSize() {return fDataSize;}
  
  ClassDef(AliHLTTrigger,1) 

};

typedef AliHLTTrigger AliL3Trigger; // for backward compatibility

#endif
