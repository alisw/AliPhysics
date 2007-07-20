#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include "AliHLTTPCPadArray.h"
#include "AliHLTTPCPad.h"
#include "AliHLTStdIncludes.h"
#include "AliHLTTPCTransform.h"
#include "AliTPCRawStream.h"
#include "AliRawReaderMemory.h"
#include "AliHLTTPCDigitReader.h"
#include <vector>
#include "AliHLTTPCClusters.h"
/** ROOT macro for the implementation of ROOT specific class methods */
//ClassImp(AliHLTTPCClusters)

AliHLTTPCClusters::AliHLTTPCClusters()
  :
  fTotalCharge(0),
  fPad(0),
  fTime(0),
  fPad2(0),
  fTime2(0),
  fMean(0),
  fFlags(1),
  fChargeFalling(0),
  fLastCharge(0),
  fLastMergedPad(0),
  fFirstPad(0),
  fLastPad(0),
  fRowNumber(0)
{

}
AliHLTTPCClusters::AliHLTTPCClusters(const AliHLTTPCClusters& src)
  :
  fTotalCharge(src.fTotalCharge),
  fPad(src.fPad),
  fTime(src.fTime),
  fPad2(src.fPad2),
  fTime2(src.fTime2),
  fMean(src.fMean),
  fFlags(src.fFlags),
  fChargeFalling(src.fChargeFalling),
  fLastCharge(src.fLastCharge),
  fLastMergedPad(src.fLastMergedPad),
  fFirstPad(src.fFirstPad),
  fLastPad(src.fLastPad),
  fRowNumber(src.fRowNumber)
{
  //HLTInfo("Copy constructor called");
}
AliHLTTPCClusters& AliHLTTPCClusters::operator=(const AliHLTTPCClusters& src){
  fTotalCharge=src.fTotalCharge;
  fPad = src.fPad;
  fTime = src.fTime;
  fPad2 = src.fPad2;
  fTime2 = src.fTime2;
  fMean = src.fMean;
  fFlags = src.fFlags;
  fChargeFalling = src.fChargeFalling;
  fLastCharge = src.fLastCharge;
  fLastMergedPad = src.fLastMergedPad;
  fFirstPad = src.fFirstPad;
  fRowNumber= src.fRowNumber;
  return (*this);
}
