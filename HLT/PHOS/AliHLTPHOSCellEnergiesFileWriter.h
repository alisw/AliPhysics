#ifndef ALIHLTPHOSCELLENERGIESFILEWRITER_H
#define ALIHLTPHOSCELLENERGIESFILEWRITER_H

#include "AliHLTPHOSFileWriterDescriptorStruct.h"
#include "AliHLTPHOSFileWriter.h"
#include <string>
#include "AliHLTDataTypes.h"


using std::string;

                                           
class AliHLTPHOSCellEnergiesFileWriter: public AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSCellEnergiesFileWriter();
  ~AliHLTPHOSCellEnergiesFileWriter();

  virtual int WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, int evntCnt);

 private:
  //  int fCurrentEvntCnt;
  //  FILE *fCurrentFile;
  //  char fCurrentFileName[256];
};


#endif
