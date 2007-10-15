//insert copyright

#ifndef ALIHLTPHOSDIGITMAKER_H
#define ALIHLTPHOSDIGITMAKER_H

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"


class AliHLTPHOSDigit;
class TClonesArray;
class TTree;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSDigitContainerDataStruct;
class AliHLTPHOSDigitDataStruct;
       
using namespace PhosHLTConst;

class AliHLTPHOSDigitMaker : public AliHLTPHOSBase
{
public:
  AliHLTPHOSDigitMaker();
  ~AliHLTPHOSDigitMaker();
 
  // void SetValidCellData(AliHLTPHOSValidCellDataStruct *data) { fCellDataPtr = data; }
  //  void SetDigitContainerStruct(AliHLTPHOSDigitContainerStruct *container) 
  //{ fDigitContainerStructPtr = container; }
                                
  void SetDigitContainerStruct(AliHLTPHOSDigitContainerDataStruct *container) 
  { fDigitContainerStructPtr = container; }
  
  void SetDigitArray(TClonesArray *array) { fDigitArrayPtr = array; }
  void ResetDigitCount() { fDigitCount = 0; }
  void SetDigitThreshold(Int_t threshold) { fDigitThreshold = threshold; }
  void SetNrPresamples(Int_t n) { fNrPresamples = n; }
  Int_t MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct*);
  void Reset();

private:

  AliHLTPHOSValidCellDataStruct *fCellDataPtr;
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerStructPtr;
  TClonesArray *fDigitArrayPtr;
  AliHLTPHOSDigit *fDigitPtr;
  //AliHLTPHOSDigitDataStruct *fDigitStructPtr;
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;
  Int_t fDigitCount; 
  Int_t fNrPresamples;

  Float_t fDigitThreshold;

  //  ClassDef(AliHLTPHOSDigitMaker, 1); 
};


#endif
 
