#ifndef ALIHLTPHOSCELLENERGIESFILEWRITER_H
#define ALIHLTPHOSCELLENERGIESFILEWRITER_H

#include "AliHLTPHOSFileWriter.h"
#include <string>
#include "AliHLTDataTypes.h"
#include <iostream>


using std::string;

                                           
class AliHLTPHOSCellEnergiesFileWriter: public AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSCellEnergiesFileWriter();
  ~AliHLTPHOSCellEnergiesFileWriter();

  //  virtual int WriteFile(const AliHLTComponentEventData& evtData, 
  //			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData,AliHLTComponentDataType dataType, int evntCnt);

 private:
  int fCurrentEvntCnt;
};


#endif
