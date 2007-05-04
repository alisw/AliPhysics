#ifndef ALIHLTPHOSFILEWRITER_H
#define ALIHLTPHOSFILEWRITER_H

#include "AliHLTDataTypes.h"
#include <iostream>
using std::string;


class AliHLTPHOSFileWriter
{
 public:
  AliHLTPHOSFileWriter();
  virtual ~AliHLTPHOSFileWriter();
  void  MakeFilename(int eventNr =0, const AliHLTComponentDataType& dataType =  kAliHLTVoidDataType);
  void SetDirectory(string& directory); 

 protected:
  int fCurrentEvntCnt;
  FILE *fCurrentFile;      /**<Flepointer to current file*/
  string fDirectory;       /**<Ouput directory for files produced by this component*/
  string fCurrentFilename; /**<Name of file for writng current data to file*/

 private:
  AliHLTPHOSFileWriter(const AliHLTPHOSFileWriter &);           /**<Never to be called*/
  AliHLTPHOSFileWriter & operator = (const AliHLTPHOSFileWriter &) /**<Never to be called*/
    {
      return *this;
    };
};


#endif
