//-*- Mode: C++ -*-
// $Id: AliHLTPHOSSharedMemoryInterfacev2.h 35107 2009-09-30 01:45:06Z phille $

#ifndef ALIHLTCALOSHAREDMEMORYINTERFACEV2_H
#define ALIHLTCALOSHAREDMEMORYINTERFACEV2_H


/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "Rtypes.h"
//#include "AliHLTPHOSBase.h"

#include "AliHLTCaloChannelRawDataStruct.h"
#include "AliHLTDataTypes.h"
#include "AliHLTCaloConstantsHandler.h"



class AliHLTCaloChannelDataHeaderStruct;
class AliHLTCaloChannelDataStruct;
class AliHLTCaloCoordinate;
class AliHLTCaloMapper;

//class AliHLTCaloChannelRawDataStruct;

class  AliHLTCaloSharedMemoryInterfacev2 : public AliHLTCaloConstantsHandler
{
 public:
  AliHLTCaloSharedMemoryInterfacev2(TString det);
  virtual ~AliHLTCaloSharedMemoryInterfacev2();
  AliHLTCaloChannelDataStruct*   NextChannel();
  void  NextRawChannel();
  // void SetMemory(AliHLTCaloChannelDataHeaderStruct* channelDataHeaderPtr, const unsigned long specification);
 
  void SetMemory(AliHLTCaloChannelDataHeaderStruct* channelDataHeaderPtr );
 
  void Reset();
  const AliHLTCaloChannelRawDataStruct & GetRawData() { return  fRawData; };
 
protected:
  AliHLTCaloMapper  *fMapperPtr[32];
  
 private:
  AliHLTCaloSharedMemoryInterfacev2();
  AliHLTCaloSharedMemoryInterfacev2(const  AliHLTCaloSharedMemoryInterfacev2 & );
  AliHLTCaloSharedMemoryInterfacev2 & operator = (const  AliHLTCaloSharedMemoryInterfacev2 &);
  void Reset(AliHLTCaloChannelRawDataStruct &str);
  AliHLTCaloChannelDataStruct* fCurrentChannel;
  AliHLTUInt8_t* fChannelDataPtr;
  bool fIsSetMemory;
  bool fHasRawData;
  int fMaxCnt;
  int fCurrentCnt; 
  UShort_t *fRawDataPtr;
  AliHLTCaloChannelRawDataStruct fRawData;
  // unsigned long fSpecification;

  ClassDef(AliHLTCaloSharedMemoryInterfacev2, 1);
};

#endif
