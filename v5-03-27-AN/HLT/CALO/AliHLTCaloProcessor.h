//-*- Mode: C++ -*-
// $Id: AliHLTCaloProcessor.h 35107 2009-09-30 01:45:06Z phille $

#ifndef ALIHLTCALOPROCESSOR_H
#define ALIHLTCALOPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHLTProcessor.h"


class AliHLTCaloProcessor:public AliHLTProcessor
{

 public:
  AliHLTCaloProcessor();
  virtual ~AliHLTCaloProcessor();
  virtual int DoInit(int argc, const char** argv) = 0;
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list) =0;
  virtual AliHLTComponentDataType GetOutputDataType() =0;
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) =0;
  virtual AliHLTComponent* Spawn() = 0; 
  char lineNumber[256];
  const char *IntToChar(int number);

 protected:
  int fCaloEventCount;                  /**<Global event counter for this component*/
  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/

 private:
  AliHLTCaloProcessor(const AliHLTCaloProcessor & );
  AliHLTCaloProcessor & operator = (const AliHLTCaloProcessor &);

};


#endif
