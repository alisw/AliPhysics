#ifndef ALIHLTPHOSFILEWRITER_H
#define ALIHLTPHOSFILEWRITER_H

//#include "AliHLTPHOSFileWriterDescriptorStruct.h"
#include "AliHLTPHOSFileWriter.h"
#include <string>
#include "AliHLTDataTypes.h"
#include <iostream>

using std::string;

class AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSFileWriter();
  virtual ~AliHLTPHOSFileWriter();
  virtual int WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, int evntCnt) = 0;
  void  MakeFilename(int eventNr, const AliHLTComponentDataType& dataType);
  void SetDirectory(string& directory); 

 protected:
  int fCurrentEvntCnt;
  FILE *fCurrentFile;
  string fDirectory;
  string fCurrentFilename;

 private:
  AliHLTPHOSFileWriter(const AliHLTPHOSFileWriter & );           /**<Never to be called*/
  AliHLTPHOSFileWriter & operator = (const AliHLTPHOSFileWriter &) /**<Never to be called*/
    {
      return *this;
    };
};


#endif
