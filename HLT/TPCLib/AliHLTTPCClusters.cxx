#if __GNUC__>= 3
using namespace std;
#endif

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include <cerrno>
//#include "AliHLTTPCPadArray.h"
//#include "AliHLTTPCPad.h"
//#include "AliHLTStdIncludes.h"
//#include "AliHLTTPCTransform.h"
//#include "AliTPCRawStream.h"
//#include "AliRawReaderMemory.h"
//#include "AliHLTTPCDigitReader.h"
//#include <vector>
#include "AliHLTTPCClusters.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusters)

AliHLTTPCClusters::AliHLTTPCClusters() :
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
  fRowNumber(0),
  fFirstPad(0),
  fLastPad(0),
  fQMax(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusters::AliHLTTPCClusters(const AliHLTTPCClusters& src) :
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
  fRowNumber(src.fRowNumber),
  fFirstPad(src.fFirstPad),
  fLastPad(src.fLastPad),
  fQMax(src.fQMax)
{
  // see header file for class documentation
}

AliHLTTPCClusters& AliHLTTPCClusters::operator=(const AliHLTTPCClusters& src)
{
  // see header file for class documentation
  if (this==&src) return *this;
  fTotalCharge=src.fTotalCharge;
  fPad = src.fPad;
  fTime = src.fTime;
  fPad2 = src.fPad2;
  fTime2 = src.fTime2;
  fMean = src.fMean;
  fFlags = src.fFlags;
  fChargeFalling = src.fChargeFalling;
  fLastCharge = src.fLastCharge;
  fRowNumber= src.fRowNumber;
  fLastMergedPad = src.fLastMergedPad;
  fFirstPad = src.fFirstPad;
  fQMax = src.fQMax;
  return (*this);
}

AliHLTTPCClusters::~AliHLTTPCClusters()
{
  // Default destructor.
}

