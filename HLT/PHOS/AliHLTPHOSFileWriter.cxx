#include  "AliHLTPHOSFileWriter.h"
#include <cstdlib>

using namespace std;


//_________________________________________________________________________________________________
AliHLTPHOSFileWriter::AliHLTPHOSFileWriter():fCurrentEvntCnt(0), fCurrentFile(0), fDirectory(""), fCurrentFilename("")
{

}


//_________________________________________________________________________________________________
AliHLTPHOSFileWriter::~AliHLTPHOSFileWriter()
{

}


//_________________________________________________________________________________________________
void
AliHLTPHOSFileWriter::SetDirectory(string& /*directory*/)
{

}


//_________________________________________________________________________________________________
void
AliHLTPHOSFileWriter::MakeFilename(int eventNr, const AliHLTComponentDataType& dataType)
{
  int charPos =fDirectory.size() +1;

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
  sprintf(tmpEvntNr,"_%.16d", eventNr);
  fCurrentFilename.insert(charPos, tmpEvntNr);
  cout <<"AliHLTPHOSFileWriterComponent::MakeFilename, filename = " << fCurrentFilename <<endl;
  
}
