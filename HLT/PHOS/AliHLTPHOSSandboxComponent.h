#ifndef ALIHLTPHOSSANDBOXCOMPONENT
#define ALIHLTPHOSSANDBOXCOMPONENT

#include "AliHLTPHOSChannelCounter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSSandboxComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSSandboxComponent();
  ~AliHLTPHOSSandboxComponent();


  /*  AliHLTPHOSSandboxComponent(const AliHLTPHOSSandboxComponent &);
  AliHLTPHOSSandboxComponent & operator = (const AliHLTPHOSSandboxComponent &)
  {
    return *this;
    }*/

  const char* GetComponentID();

  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /*
  int DoProcessing(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
	      AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		   std::vector<AliHLTComponentBlockData>&, AliHLTComponentEventDoneData *&);
  */
  int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
	      AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		   std::vector<AliHLTComponentBlockData>&);

  AliHLTComponent* Spawn();
  
protected:
  int DoInit(int argc, const char** argv);
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor

  virtual int DoDeinit();
  
private:

  Int_t fEventCount;
  AliHLTPHOSChannelCounter *fChannelCounterPtr;
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
};
#endif
