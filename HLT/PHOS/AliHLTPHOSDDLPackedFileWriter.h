#ifndef ALIHLTPHOSDDLPACKEDFILEWRITER_H
#define ALIHLTPHOSDDLPACKEDFILEWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTPHOSFileWriter.h"
#include <string>
#include "AliHLTDataTypes.h"

using std::string;

class AliHLTPHOSDDLPackedFileWriter: public AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSDDLPackedFileWriter();
  virtual ~AliHLTPHOSDDLPackedFileWriter();

  const virtual int WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, int evntCnt) const;
  

};


#endif
