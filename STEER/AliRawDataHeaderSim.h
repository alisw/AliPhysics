#ifndef ALIRAWDATAHEADERSIM_H
#define ALIRAWDATAHEADERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TError.h>

#include "AliRawDataHeader.h"
#include "AliRunLoader.h"
#include "AliCentralTrigger.h"

class AliRawDataHeaderSim : public AliRawDataHeader {

public:
  AliRawDataHeaderSim() : AliRawDataHeader() {
    // Takes the trigger mask and
    // stores it in the data header
    AliRunLoader *runloader = AliRunLoader::Instance();
    if (runloader) {
      if(!runloader->GetTrigger()) runloader->LoadTrigger();
      if (AliCentralTrigger *aCTP = runloader->GetTrigger()) {
	ULong64_t mask = aCTP->GetClassMask();
	SetTriggerClass(mask);
      }
      else
	Warning("SetTriggerClass","No trigger can be loaded! Putting empty trigger class into the raw data header !");
    }
    else
      Error("SetTriggerClass","No run loader is available! Putting empty trigger class into the raw data header !");
  }

};

#endif
