//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSSHAREDMEMORYINTERFACE_H
#define ALIHLTPHOSSHAREDMEMORYINTERFACE_H

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

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSValidCellDataStruct;

class  AliHLTPHOSSharedMemoryInterface
{
 public:
  AliHLTPHOSSharedMemoryInterface();
  virtual ~AliHLTPHOSSharedMemoryInterface();
  AliHLTPHOSValidCellDataStruct*   NextChannel();
  void SetMemory(AliHLTPHOSRcuCellEnergyDataStruct *rcuCeelEnergyPtr);
  Int_t* GetRawData(Int_t& nSamples); //added by OD
  void Reset();

 private:
   AliHLTPHOSSharedMemoryInterface(const  AliHLTPHOSSharedMemoryInterface & );
   AliHLTPHOSSharedMemoryInterface & operator = (const  AliHLTPHOSSharedMemoryInterface &);

  void PingPongPointer();
  AliHLTPHOSValidCellDataStruct *fCurrentChannel;
  AliHLTPHOSRcuCellEnergyDataStruct *fCellEnergiesPtr ;
  bool fIsSetMemory;
  bool fHasRawData;
  int fMaxCnt;
  int fCurrentCnt; 
  int fCurrentX;   //added by OD
  int fCurrentZ;  //added by OD
  int fCurrentGain;//added by OD
  Int_t fCharDataOffset;
  char  *fCharPtr;
  Int_t *fIntPtr;
  Int_t *fRawDataPtr;
  //  Int_t *rawDataBufferPos = (Int_t *)outputPtr; 
};

#endif
