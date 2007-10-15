//insert copyright

#ifndef ALIHLTPHOSBASELINEANALYZERCOMPONENT_H
#define ALIHLTPHOSBASELINEANALYZERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"

struct AliHLTPHOSBaselineAnalyzer;
struct TTree;

class AliHLTPHOSBaselineAnalyzerComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSBaselineAnalyzerComponent();
  ~AliHLTPHOSBaselineAnalyzerComponent();

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
  
  void CalculateAll();

  AliHLTPHOSBaselineAnalyzer *fBaselineAnalyzerPtr;
  TTree *fTreePtr;
  TClonesArray *fBaselineArrayPtr;
  UInt_t fEventCount;
  UInt_t fWriteInterval;
  UInt_t fFillInterval;
  char *fFilename;
  char* fDirectory;
  char* fHistPath;
  Int_t fRunNb;
  Bool_t fCalculateAll;

  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
