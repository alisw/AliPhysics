#ifndef ALIHLTPHOSCLUSTERIZERCOMPONENT
#define ALIHLTPHOSCLUSTERIZERCOMPONENT

#include "AliHLTProcessor.h"

#include "AliHLTPHOSClusterizer.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSClusterDataStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "Rtypes.h"


class AliHLTPHOSClusterizerComponent: public AliHLTProcessor
{
 public:

  AliHLTPHOSClusterizerComponent();
  ~AliHLTPHOSClusterizerComponent();
  AliHLTPHOSClusterizerComponent(const AliHLTPHOSClusterizerComponent &);
  AliHLTPHOSClusterizerComponent & operator = (const AliHLTPHOSClusterizerComponent &)
    {
      return *this;
    }
  const char* GetComponentID();
  void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  Int_t DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
		AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		std::vector<AliHLTComponentBlockData>&);

  AliHLTComponent* Spawn();

 protected:

  Int_t DoInit(int argc, const char** argv);
  Int_t Deinit();
  Int_t DoDeinit();

 private:
  AliHLTPHOSClusterizer* fClusterizerPtr;
  AliHLTPHOSClusterDataStruct* fOutPtr;
  AliHLTPHOSRecPointDataStruct* fRecPointStructArrayPtr;
  AliHLTPHOSRecPointListDataStruct* fRecPointListPtr;
  static const AliHLTComponentDataType inputDataTypes[];
  static int fEventCount;

};

#endif
