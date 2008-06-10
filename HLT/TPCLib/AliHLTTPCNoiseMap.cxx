// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCNoiseMap.cxx
    @author Kalliopi Kanaki
    @date   06.06.2008
    @brief  Class for reading the noise map from HCDB.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCNoiseMap.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBPath.h"

#include "TObjString.h"

ClassImp(AliHLTTPCNoiseMap)

AliHLTTPCNoiseMap::AliHLTTPCNoiseMap(){ 
// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //
}

AliHLTTPCNoiseMap* AliHLTTPCNoiseMap::pNoiseMapInstance = 0; // initialize pointer

AliHLTTPCNoiseMap* AliHLTTPCNoiseMap::Instance(){
// see header file for class documentation

  if (pNoiseMapInstance == 0){  
      pNoiseMapInstance = new AliHLTTPCNoiseMap; // create sole instance
  }
  return pNoiseMapInstance; // address of sole instance
}

AliTPCCalPad* AliHLTTPCNoiseMap::ReadNoiseMap(){
// see header file for class documentation
  
  const char* pathNoiseMap = "TPC/Calib/PadNoise";
  AliTPCCalPad *padNoise = NULL;

  if(pathNoiseMap){    
    
    //AliCDBPath path(pathNoiseMap);
    //AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path); // read from the default storage 
    // for local testing purposes
    
    AliCDBStorage *stor   = AliCDBManager::Instance()->GetStorage(pathNoiseMap);
    AliCDBEntry   *pEntry = stor->Get(pathNoiseMap,AliCDBManager::Instance()->GetRun());
    
    if(pEntry){ 
       padNoise = (AliTPCCalPad*)pEntry->GetObject();
       TObjString *pString = dynamic_cast<TObjString*>(pEntry->GetObject());
       if(pString){ HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data()); }
    } // end if pEntry 
    else { HLTError("cannot fetch object \"%s\" from CDB", pathNoiseMap); } 
  } // end if pathNoiseMap
  return padNoise;
}
