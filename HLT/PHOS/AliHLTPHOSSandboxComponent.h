/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSSANDBOXCOMPONENT
#define ALIHLTPHOSSANDBOXCOMPONENT

//#include "AliHLTPHOSChannelCounter.h"
//#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSChannelCounter;
class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSSandboxComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSSandboxComponent();
  virtual ~AliHLTPHOSSandboxComponent();


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
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
		   std::vector<AliHLTComponentBlockData>& outputBlocks);

  AliHLTComponent* Spawn();
  
protected:
  int DoInit(int argc, const char** argv);
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor

  virtual int DoDeinit();
  
private:

  Int_t fEvtCnt; //comment
  AliHLTPHOSChannelCounter *fChannelCounterPtr; //comment
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
};
#endif
