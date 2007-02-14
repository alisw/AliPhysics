#include  "AliHLTPHOSCellEnergiesFileWriter.h"



AliHLTPHOSCellEnergiesFileWriter::AliHLTPHOSCellEnergiesFileWriter():fCurrentEvntCnt(0)
{

}


AliHLTPHOSCellEnergiesFileWriter::~AliHLTPHOSCellEnergiesFileWriter()
{

}


/*
int 
AliHLTPHOSCellEnergiesFileWriter::WriteFile(const AliHLTComponentEventData& evtData, 
			const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTComponentDataTyp dataType, int evntCnt)
{
  cout <<"AliHLTPHOSCellEnergiesFileWriter::WriteFile" << endl;

  if(evntCnt != fCurrentEvntCnt)
    {
      if(fCurrentFile != 0)
	{
	  
	}
    }

  return 0;

}
*/
