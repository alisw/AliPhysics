#ifndef ALIHLTTPCCLUSTERS
#define ALIHLTTPCCLUSTERS

#include "AliHLTLogging.h"

class AliHLTTPCClusters : public AliHLTLogging {

 public:
  AliHLTTPCClusters();
  AliHLTTPCClusters(const AliHLTTPCClusters& src);
  AliHLTTPCClusters& operator=(const AliHLTTPCClusters&);

  UInt_t fTotalCharge;   //tot charge of cluster
  UInt_t fPad;           //pad value
  UInt_t fTime;          //time value
  ULong64_t fPad2;       //for error in XY direction
  ULong64_t fTime2;      //for error in Z  direction
  UInt_t fMean;          //mean in time
  UInt_t fFlags;         //different flags
  UInt_t fChargeFalling; //for deconvolution
  UInt_t fLastCharge;    //for deconvolution
  UInt_t fLastMergedPad; //dont merge twice per pad
  UInt_t fRowNumber;
  Int_t fFirstPad;
  UInt_t fLastPad;
  ClassDef(AliHLTTPCClusters,0) //Fast cluster finder
    };
#endif
