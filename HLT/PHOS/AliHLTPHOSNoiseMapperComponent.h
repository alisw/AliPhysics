//insert copyright

#ifndef ALIHLTPHOSNOISEMAPPERCOMPONENT_H
#define ALIHLTPHOSNOISEMAPPERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSNoiseMapper;
class TH2I;


class AliHLTPHOSNoiseMapperComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSNoiseMapperComponent();
  ~AliHLTPHOSNoiseMapperComponent(); 

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
  void FillHistograms();
  
private:
  
  void CalculateAll();

  AliHLTPHOSNoiseMapper *fNoiseMapperPtr;
  UInt_t fEventCount;
  UInt_t fWriteInterval;
  char* fDirectory;
  char* fFilename;
  Int_t fRunNb;
  Int_t *fChannelArrayPtr;
  TH2I *fNoiseCountLowGainHistPtr;
  TH2I *fNoiseCountHighGainHistPtr;
  TH2I *fNoiseMapLowGainHistPtr;     
  TH2I *fNoiseMapHighGainHistPtr;
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
