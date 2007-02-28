#ifndef ALIHLTPHOSFILEWRITERCOMPONENT_H
#define ALIHLTPHOSFILEWRITERCOMPONENT_H

#include "AliHLTDataSink.h"
#include "AliHLTPHOSFileWriterComponent.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSDefinitions.h"
#include <string>
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSFileWriter.h"
#include "AliHLTPHOSCellEnergiesFileWriter.h"
#include "AliHLTPHOSDDLPackedFileWriter.h"  
#include "Rtypes.h"
#include "AliHLTFileWriter.h"

using std::string;


class AliHLTPHOSRcuCellEnergyDataStruct;
 

class AliHLTPHOSFileWriterComponent:public AliHLTFileWriter
{
 public:
  AliHLTPHOSFileWriterComponent();
  virtual ~AliHLTPHOSFileWriterComponent();
  int AddDataType(string dataType);
  virtual int Deinit();
  virtual int DoDeinit();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();

 protected:

 private:
  Bool_t IsRegisteredDataType(const AliHLTComponentDataType&);
  AliHLTPHOSFileWriterComponent(const AliHLTPHOSFileWriterComponent & );  
  AliHLTPHOSFileWriterComponent & operator = (const AliHLTPHOSFileWriterComponent)
    {
      return *this;
    };
  AliHLTPHOSCellEnergiesFileWriter *fCellEnergiesFileWriterPtr;
  AliHLTPHOSDDLPackedFileWriter    *fDDLPackedFileWriterPtr ;
  string  fDirectory; /**<target directory for files*/
  string  fFilename;  /**<the basename of the output file*/
  AliHLTComponentDataType fDataTypesToFile[N_DATATYPES];
  int fEventCount;
  static const AliHLTComponentDataType fInputDataTypes[];
};
#endif
