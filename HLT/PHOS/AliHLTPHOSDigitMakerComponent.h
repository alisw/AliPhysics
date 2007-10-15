//insert copyright

#ifndef ALIHLTPHOSDIGITMAKERCOMPONENT_H
#define ALIHLTPHOSDIGITMAKERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"
//#include "AliHLTPHOSDigitMaker.h"
//#include "TTree.h"
//#include "TClonesArray.h"


class AliHLTPHOSDigitMaker;
class TTree;
class TClonesArray;
class AliHLTPHOSDigitContainerDataStruct;



class AliHLTPHOSDigitMakerComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSDigitMakerComponent();
  ~AliHLTPHOSDigitMakerComponent();

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
  AliHLTPHOSDigitMaker *fDigitMakerPtr;
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerPtr;
  UInt_t fEventCount;
  Int_t fRunNb;

  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
 
