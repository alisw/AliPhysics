#include  "AliHLTPHOSFileWriter.h"
#include <string>
#include <cstdlib>

using std::string;


AliHLTPHOSFileWriter::AliHLTPHOSFileWriter():fCurrentEvntCnt(0), fCurrentFile(0), fDirectory(0), fCurrentFilename(0)
{

}


AliHLTPHOSFileWriter::~AliHLTPHOSFileWriter()
{

}


void
AliHLTPHOSFileWriter::SetDirectory(string& directory)
{

}

void
AliHLTPHOSFileWriter::MakeFilename(int eventNr, const AliHLTComponentDataType& dataType)
{

  int charPos =fDirectory.size() +1;
  cout <<"charPos.size() = "<< charPos << endl;

  fCurrentFilename.erase(charPos);

  char tmpOr[kAliHLTComponentDataTypefOriginSize+1];
  char tmpID[kAliHLTComponentDataTypefIDsize+1];
  char tmpEvntNr[30];

  for(int i = 0; i< kAliHLTComponentDataTypefOriginSize; i++)
    {
      tmpOr[i] = dataType.fOrigin[i]; 
    }
  tmpOr[kAliHLTComponentDataTypefOriginSize] = '\0';

  for(int j = 0; j< kAliHLTComponentDataTypefIDsize; j++)
    {
      tmpID[j] = dataType.fID[j]; 
    }
  tmpID[kAliHLTComponentDataTypefIDsize] = '\0';

  fCurrentFilename.insert(charPos, tmpOr);
  charPos+= kAliHLTComponentDataTypefOriginSize;
  fCurrentFilename.insert(charPos, tmpID);
  charPos = fCurrentFilename.size();
  sprintf(tmpEvntNr,"%.16d", eventNr);
  fCurrentFilename.insert(charPos, tmpEvntNr);
  cout <<"AliHLTPHOSFileWriterComponent::MakeFilename, filename = " << fCurrentFilename <<endl;
  
}
