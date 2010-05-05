//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSPROCESSOR_H
#define ALIHLTPHOSPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHLTProcessor.h"
//#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTDataTypes.h"

// #include "AliHLTPHOSConstant.h"

// using namespace PhosHLTConst;

//class AliHLTPHOSProcessor:public AliHLTProcessor, public AliHLTPHOSBase
class AliHLTPHOSProcessor:public AliHLTProcessor
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


  char lineNumber[256];
  const char *IntToChar(int number);
  /*
   *Check file and write messages to AliLog system
   */
  bool CheckFileLog(const char *origin, const char *filename, const char *opt);
  void DoneWritingLog(const char *origin, const char *filename);

  using  AliHLTProcessor::DoEvent;

 protected:
  void ScanRunNumberFromFile();
  virtual int ScanArgumentsModule(int argc, const char** argv);
  int fPhosEventCount;                  /**<Global event counter for this component*/
  AliHLTUInt8_t  fModuleID;             /**<ID of the module this component read data from (0-4)*/

  Bool_t fPrintInfoModule;                    /**<wether or not to print debugg info to std out*/
  int fPrintInfoFrequncyModule;               /**<Defines the update frequency for information printet to std out*/

  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/
  int fRunNumber;
  char fFilepath[1024];
  char fMessage[1024];

 private:
  AliHLTPHOSProcessor(const AliHLTPHOSProcessor & );
  AliHLTPHOSProcessor & operator = (const AliHLTPHOSProcessor &);

};


#endif
