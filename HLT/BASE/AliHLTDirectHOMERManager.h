//-*- Mode: C++ -*-
#ifndef ALIHLTDIRECTHOMERMANAGER_H
#define ALIHLTDIRECTHOMERMANAGER_H


#include "AliHLTHOMERReader.h"
#include "TString.h"
#include "TList.h"
#include <list>


class AliHLTDirectHOMERManager {
public:

  AliHLTDirectHOMERManager(int argc, char** argv);  
  ~AliHLTDirectHOMERManager();

  TList* GetNextBlocks();

private:
  AliHLTDirectHOMERManager();

  TString MakeOrigin(homer_uint32 origin);
  TString MakeDataType(homer_uint64 dataType);


  typedef std::pair<TString, unsigned short> connectionID;
  std::list<connectionID> fconnectionList;

  TList* fBlockList;

};

#endif //ALIHLTDIRECTHOMERMANAGER_H
