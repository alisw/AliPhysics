//insert copyright

#ifndef ALIHLTPHOSRCUDIGITMAKERCOMPONENT_H
#define ALIHLTPHOSRCUDIGITMAKERCOMPONENT_H

#include "AliHLTPHOSRcuProcessor.h"
//#include "AliHLTPHOSDigitMaker.h"
//#include "TTree.h"
//#include "TClonesArray.h"


class AliHLTPHOSRcuDigitMaker;
class TTree;
class TClonesArray;
class AliHLTPHOSRcuDigitContainerDataStruct;



class AliHLTPHOSRcuDigitMakerComponent : public AliHLTPHOSRcuProcessor
{
public:
  AliHLTPHOSRcuDigitMakerComponent();
  ~AliHLTPHOSRcuDigitMakerComponent();

  const char* GetComponentID();

  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
	      AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
	      std::vector<AliHLTComponentBlockData>&);
  
  AliHLTComponent* Spawn();
  
protected:
  int DoInit(int argc, const char** argv);
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
private:
  AliHLTPHOSRcuDigitMaker *fDigitMakerPtr;
  AliHLTPHOSRcuDigitContainerDataStruct *fDigitContainerPtr;
  UInt_t fEventCount;
  Int_t fRunNb;

  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
 
