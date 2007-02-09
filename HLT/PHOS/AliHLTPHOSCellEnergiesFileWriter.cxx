#include  "AliHLTPHOSCellEnergiesFileWriter.h"



AliHLTPHOSCellEnergiesFileWriter::AliHLTPHOSCellEnergiesFileWriter()
{

}


AliHLTPHOSCellEnergiesFileWriter::~AliHLTPHOSCellEnergiesFileWriter()
{

}


int 
AliHLTPHOSCellEnergiesFileWriter::WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, int evntCnt)
{
  cout <<"AliHLTPHOSCellEnergiesFileWriter::WriteFile" << endl;

  if(fCurrentFile != 0)
    {
      
    }

  return 0;

}
