#ifndef ALIHLTPHOSDDLPACKEDFILEWRITER_H
#define ALIHLTPHOSDDLPACKEDFILEWRITER_H

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
