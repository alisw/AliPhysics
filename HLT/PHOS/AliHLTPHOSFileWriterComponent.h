#ifndef ALIHLTPHOSFILEWRITER_H
#define ALIHLTPHOSFILEWRITER_H

#include "AliHLTComponent.h"
#include "AliHLTFileWriter.h" 

/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliHLTPHOSFileWriterComponent : public AliHLTFileWriter
{
 public:
  AliHLTPHOSFileWriterComponent();
  ~AliHLTPHOSFileWriterComponent();
  //  virtual AliHLTComponent* Spawn();
  AliHLTComponent* Spawn();

  int DoInit( int argc, const char** argv );

  const char* GetComponentID();
};

#endif
