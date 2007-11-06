/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTPHOSClusterizerComponent.h
    @author Ãystein Djuvsland
    @date   
    @brief  A clusterizer component for PHOS HLT
*/


#ifndef ALIHLTPHOSCLUSTERIZERCOMPONENT_H
#define ALIHLTPHOSCLUSTERIZERCOMPONENT_H



#include "AliHLTPHOSProcessor.h"

//#include "AliHLTPHOSBase.h"
//#include "AliHLTPHOSDefinitions.h"
//#include "AliHLTProcessor.h"

class AliHLTPHOSClusterizer;

//class Rtypes;

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSClusterDataStruct;
class AliHLTPHOSRecPointDataStruct;
class AliHLTPHOSRecPointContainerStruct;
class AliHLTPHOSRecPointListDataStruct;
class AliHLTPHOSDigitContainerDataStruct;



// PTH class AliHLTPHOSClusterizerComponent:  public AliHLTPHOSBase, public AliHLTProcessor
class AliHLTPHOSClusterizerComponent: public AliHLTPHOSProcessor
//class AliHLTPHOSClusterizerComponent:  public AliHLTPHOSBase, public AliHLTProcessor
{
 public:

  AliHLTPHOSClusterizerComponent();
  virtual ~AliHLTPHOSClusterizerComponent();

  //  AliHLTPHOSClusterizerComponent(const AliHLTPHOSClusterizerComponent &);
  //  AliHLTPHOSClusterizerComponent & operator = (const AliHLTPHOSClusterizerComponent &)
  //   {

  //     return *this;
  //   }



  const char* GetComponentID();
  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks);

  AliHLTComponent* Spawn();
  
 // void SetNoCrazyness(Bool_t val);

 protected:


  int DoInit(int argc, const char** argv);
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  int DoDeinit();

 private:
  AliHLTPHOSDigitContainerDataStruct *fAllDigitsPtr;            //comment
  AliHLTPHOSClusterizer* fClusterizerPtr;                       //Pointer to the clusterizer
  AliHLTPHOSRecPointContainerStruct* fOutPtr;                         //Pointer to the struct of output clusters
  AliHLTPHOSRecPointDataStruct* fRecPointStructArrayPtr;        //Pointer to the struct of output recpoints
  AliHLTPHOSRecPointListDataStruct* fRecPointListPtr;           //Pointer to the struct of list of output recpoints
  Int_t fDigitCount;                                            //comment
  Bool_t fNoCrazyness;                                          //comment
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
};

#endif
