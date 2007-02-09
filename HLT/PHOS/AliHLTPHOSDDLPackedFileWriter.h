#ifndef ALIHLTPHOSDDLPACKEDFILEWRITER_H
#define ALIHLTPHOSDDLPACKEDFILEWRITER_H

#include "AliHLTPHOSFileWriterDescriptorStruct.h"
#include "AliHLTPHOSFileWriter.h"
#include <string>
#include "AliHLTDataTypes.h"


using std::string;

class AliHLTPHOSDDLPackedFileWriter: public AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSDDLPackedFileWriter();
  ~AliHLTPHOSDDLPackedFileWriter();

virtual int WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, int evntCnt);


};


#endif
