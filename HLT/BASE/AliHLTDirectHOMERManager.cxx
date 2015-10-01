//-*- Mode: C++ -*-
#include "AliHLTDirectHOMERManager.h"
#include <AliHLTHOMERReader.h>
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTOUT.h"
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <iostream>
  
AliHLTDirectHOMERManager::AliHLTDirectHOMERManager(int argc, char** argv) : fconnectionList(), fBlockList(0) {
  for (int i = 1; i<argc; ++i){
    TString addr = argv[i];
    TObjArray* tokens = addr.Tokenize(":");
    if(tokens->GetEntriesFast() == 2){
      TObjString* hostname = dynamic_cast<TObjString*>(tokens->At(0));
      TObjString* port = dynamic_cast<TObjString*>(tokens->At(1));
      std::cout << "Adding host \'" << hostname->String().Data() << "\', port \'" << port->String().Data() << "\'" << std::endl;
      fconnectionList.push_back( connectionID(hostname->String(), (unsigned short) port->String().Atoi()) );
    }
    delete tokens;
  }
}

AliHLTDirectHOMERManager::~AliHLTDirectHOMERManager() {
  if (fBlockList) delete fBlockList;
  fBlockList = NULL;
}


TList* AliHLTDirectHOMERManager::GetNextBlocks() {
  if (fBlockList) delete fBlockList;
  fBlockList = new TList();
  fBlockList->SetOwner(kTRUE);
  std::list<connectionID>::iterator iter;
  int result=0;
  unsigned long timeout = 100000000;
  for ( iter = fconnectionList.begin(); iter != fconnectionList.end(); ++iter ){
    AliHLTHOMERReader h(iter->first, iter->second);
    result = h.ReadNextEvent( timeout );
    if (result) {
      std::cout << "Error getting blocks from " << iter->first << ":" << iter->second << ": "
		<< strerror(result) << " (" << result << ")." << std::endl;
      continue;
    }
    for(unsigned long idx = 0; idx < h.GetBlockCnt(); ++idx){
      AliHLTHOMERBlockDesc * block = new AliHLTHOMERBlockDesc();
      
      block->SetBlock(const_cast<void*>( h.GetBlockData(idx)), h.GetBlockDataLength(idx), 
		      MakeOrigin(h.GetBlockDataOrigin(idx)),
		      MakeDataType(h.GetBlockDataType(idx)), 
		      (ULong_t) h.GetBlockDataSpec(idx) 
		      );
      fBlockList->Add(block);
    }
  }

  return fBlockList;
}


TString AliHLTDirectHOMERManager::MakeOrigin(homer_uint32 origin) {
  TString ret;
  union{
    AliHLTUInt32_t data;
    Char_t array[4];
  } reverseOrigin;

  reverseOrigin.data = AliHLTOUT::ByteSwap32(origin);
  for(int i=0; i<4; ++i)
    ret.Append(reverseOrigin.array[i]);

  ret.Remove( TString::kTrailing, ' ' );
  return ret;
}

TString AliHLTDirectHOMERManager::MakeDataType(homer_uint64 type) {
  TString ret;
  union{
    AliHLTUInt64_t data;
    Char_t array[8];
  } reverseType;

  reverseType.data = AliHLTOUT::ByteSwap64(type);
  for(int i=0; i<8; ++i)
    ret.Append(reverseType.array[i]);

  ret.Remove( TString::kTrailing, ' ' );
  return ret;
}
