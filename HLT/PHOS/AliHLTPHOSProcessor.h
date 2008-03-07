#ifndef ALIHLTPHOSPROCESSOR_H
#define ALIHLTPHOSPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHLTProcessor.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTDataTypes.h"

using namespace PhosHLTConst;

class AliHLTPHOSProcessor:public AliHLTProcessor, public AliHLTPHOSBase
{

 public:
  AliHLTPHOSProcessor();
  virtual ~AliHLTPHOSProcessor();
  virtual int DoInit(int argc, const char** argv) = 0;
  virtual int Deinit() = 0;
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list) =0;
  virtual AliHLTComponentDataType GetOutputDataType() =0;
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) =0;
  virtual AliHLTComponent* Spawn() = 0; 
 
  using  AliHLTProcessor::DoEvent;

 protected:
  void ScanRunNumberFromFile();
  virtual int ScanArguments(int argc, const char** argv);
  int fPhosEventCount;                  /**<Global event counter for this component*/
  AliHLTUInt8_t  fModuleID;             /**<ID of the module this component read data from (0-4)*/
  Bool_t fPrintInfo;                    /**<wether or not to print debugg info to std out*/
  int fPrintInfoFrequncy;               /**<Defines the update frequency for information printet to std out*/
  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/
  int fRunNumber;
 private:
  AliHLTPHOSProcessor(const AliHLTPHOSProcessor & );
  AliHLTPHOSProcessor & operator = (const AliHLTPHOSProcessor &);

};


#endif
