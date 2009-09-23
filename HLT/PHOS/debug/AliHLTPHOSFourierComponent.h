//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSFOURIERCOMPONENT_H
#define ALIHLTPHOSFOURIERCOMPONENT_H

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


#include "AliHLTPHOSRcuProcessor.h"

class  AliHLTPHOSFourier;
class  AliHLTPHOSSharedMemoryInterface;
class AliHLTPHOSRcuFFTDataStruct;

class  AliHLTPHOSFourierComponent :  public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSFourierComponent();
  virtual ~AliHLTPHOSFourierComponent();
  int DoInit(int argc, const char** argv);
  virtual int Deinit();
  int DoEvent(const AliHLTComponentEventData& evtData,
	      const AliHLTComponentBlockData* blocks, 
	      AliHLTComponentTriggerData& trigData,
	      AliHLTUInt8_t* outputPtr, 
	      AliHLTUInt32_t& size,
	      AliHLTComponentBlockDataList& outputBlocks );

  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn();
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
 
  //virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);

protected:
  //  int DoInit(int argc, const char** argv);
  //  int DoDeinit();
 
private:
  AliHLTPHOSFourierComponent  (const  AliHLTPHOSFourierComponent  & ): AliHLTPHOSRcuProcessor(),
    fFourierPtr(0),
    fShmPtr(0),
    fOutPtr(0)
      {
      };

  AliHLTPHOSFourierComponent  & operator = (const AliHLTPHOSFourierComponent)
  {
    return *this;
  };

  AliHLTPHOSFourier *fFourierPtr;
  AliHLTPHOSSharedMemoryInterface *fShmPtr; // Interface to read altro channel data from shared memory
  AliHLTPHOSRcuFFTDataStruct* fOutPtr;  //comment
};

#endif
