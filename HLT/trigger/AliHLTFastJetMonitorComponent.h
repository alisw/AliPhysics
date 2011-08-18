#ifndef ALIHLTFASTJETMONITORCOMPONENT_H
#define ALIHLTFASTJETMONITORCOMPONENT_H

#include "AliHLTProcessor.h"

class AliHLTFastJetMonitor;

class AliHLTFastJetMonitorComponent : public AliHLTProcessor
{
 public:
  // Constructor
  AliHLTFastJetMonitorComponent();

  // Destructor
  virtual ~AliHLTFastJetMonitorComponent();

  // interface function
  const char* GetComponentID();

  // interface function
  void GetInputDataTypes(std::vector<AliHLTComponentDataType> &list);

  // interface function
  AliHLTComponentDataType GetOutputDataType();

  // interface function
  void GetOutputDataSize(unsigned long &constBase, double &inputMultiplier);

  // interface function
  int DoEvent(const AliHLTComponentEventData &evtData, const AliHLTComponentBlockData *blocks,
	      AliHLTComponentTriggerData &/*trigData*/, AliHLTUInt8_t */*outputPtr*/, AliHLTUInt32_t &/*size*/,
	      std::vector<AliHLTComponentBlockData> &/*outputBlocks*/);
  
  // interface function
  AliHLTComponent* Spawn();

 protected:
  
  // interface function
  int DoInit(int argc, const char **argv);
  int DoDeinit() {return 0;};
  
  using AliHLTProcessor::DoEvent;

  // interface function
  virtual int Deinit();


 private:
  TString fRootFileName;
  int     fPushFraction;
  int     fLocalEventCount;
  int     fVerbose;

  // pointer to the histo maker itself
  AliHLTFastJetMonitor *fHistoPtr;

  AliHLTFastJetMonitorComponent(const AliHLTFastJetMonitorComponent &);
  AliHLTFastJetMonitorComponent & operator = (const AliHLTFastJetMonitorComponent &);

};

#endif
