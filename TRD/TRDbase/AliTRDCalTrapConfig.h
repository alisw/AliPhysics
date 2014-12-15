/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALITRDCALTRAPCONFIG
#define ALITRDCALTRAPCONFIG

#include "TObject.h"
#include "TList.h"

#include "AliTRDtrapConfig.h"

class AliTRDCalTrapConfig : public TObject
{
public:
  AliTRDCalTrapConfig();
  ~AliTRDCalTrapConfig();

  void Add(AliTRDtrapConfig *cfg) { fConfigList.Add(cfg); }

  virtual void Print(Option_t *option = "") const;

  AliTRDtrapConfig* Get(const TString &name) const;

protected:
  TList fConfigList;

  ClassDef(AliTRDCalTrapConfig, 1);
};

#endif
