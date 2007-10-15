//insert copyright

#ifndef ALIHLTPHOSRCUDIGITMAKER_H
#define ALIHLTPHOSRCUDIGITMAKER_H

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"
//#include "AliHLTPHOSRcuProcessor.h"


class AliHLTPHOSDigit;
class TClonesArray;
class TTree;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuDigitContainerDataStruct;
class AliHLTPHOSDigitDataStruct;
       
using namespace PhosHLTConst;

class AliHLTPHOSRcuDigitMaker : public AliHLTPHOSBase
{
public:
  AliHLTPHOSRcuDigitMaker();
  ~AliHLTPHOSRcuDigitMaker();
 
  // void SetValidCellData(AliHLTPHOSValidCellDataStruct *data) { fCellDataPtr = data; }
  //  void SetDigitContainerStruct(AliHLTPHOSDigitContainerStruct *container) 
  //{ fDigitContainerStructPtr = container; }
                                
  void SetDigitContainerStruct(AliHLTPHOSRcuDigitContainerDataStruct *container) 
  { fDigitContainerStructPtr = container; }
  
  void SetDigitArray(TClonesArray *array) { fDigitArrayPtr = array; }
  void ResetDigitCount() { fDigitCount = 0; }
  void SetDigitThreshold(Int_t threshold) { fDigitThreshold = threshold; }
  void SetNrPresamples(Int_t n) { fNrPresamples = n; }
  Int_t MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct*);
  void Reset();
  //  virtual void De

private:

  AliHLTPHOSValidCellDataStruct *fCellDataPtr;
  AliHLTPHOSRcuDigitContainerDataStruct *fDigitContainerStructPtr;
  TClonesArray *fDigitArrayPtr;
  AliHLTPHOSDigit *fDigitPtr;
  //AliHLTPHOSDigitDataStruct *fDigitStructPtr;
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;
  Int_t fDigitCount; 
  Int_t fNrPresamples;

  Float_t fDigitThreshold;

  //  ClassDef(AliHLTPHOSRcuDigitMaker, 1); 
};


#endif
 
