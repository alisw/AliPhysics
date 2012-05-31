// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTestProcessor.cxx
    @author Matthias Richter
    @date   
    @brief  A processor used for testing the external wrapper interface
 */

#include "AliHLTTestProcessor.h"
#include "TMath.h"
#include "AliTracker.h"

AliHLTTestProcessor::AliHLTTestProcessor()
  : AliHLTProcessor()
  , fState(kCreated)
  , fProcessedTypes() 
{
}

AliHLTTestProcessor::~AliHLTTestProcessor()
{
}

int AliHLTTestProcessor::DoInit( int argc, const char** argv ) {
  if (fState!=kCreated) {
    HLTError("wrong state (%d): component already initialized", fState);
    return -EBUSY;
  }

  if (argc>0 && argv) {
  }

  fState=kInitialized;
  return 0;
}
int AliHLTTestProcessor::DoDeinit() {
  if (fState!=kInitialized) {
    HLTError("wrong state (%d): required %d kInitialized", fState, kInitialized);
    return -ENODEV;
  }

  return 0;
}

int AliHLTTestProcessor::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/) {
  int iResult=0;
  if (fState!=kInitialized) {
    HLTError("wrong state (%d): component not initialized", fState);
    return -ENODEV;
  }

  if (!CheckMagneticField(-5)) {
    HLTError("initialization of magnetic field failed");
    return -EFAULT;
  }  

  HLTInfo("processing event");
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
	 pBlock!=NULL; 
	 pBlock=GetNextInputBlock()) {
    //gInputDt=pBlock->fDataType;
    AliHLTComponentDataType dt;
    SetDataType(dt, "-OUTPUT-", "MYCO");
    iResult=PushBack(pBlock->fPtr, pBlock->fSize/2, dt, ~pBlock->fSpecification);
    Forward();
  }

  return iResult;
}

bool AliHLTTestProcessor::CheckRunNo(unsigned runNo) const
{
  return runNo==GetRunNo();
}

bool AliHLTTestProcessor::CheckChainId(const char* chainId) const
{
  return strcmp(chainId, GetChainId())==0;
}

bool AliHLTTestProcessor::CheckDataType(const char* id, const char* origin) const
{
  AliHLTComponentDataType dt;
  SetDataType(dt, id, origin);
  for (unsigned i=0; i<fProcessedTypes.size(); i++) {
    if (MatchExactly(dt, fProcessedTypes[i])) return true;
  }
  return false;
}

bool AliHLTTestProcessor::CheckMagneticField(float bz)
{
  if (TMath::Abs(AliTracker::GetBz()-bz)<0.01) return true;
  HLTError("mismatch in magnetic field: %f, expecting %f", AliTracker::GetBz(), bz);
  return false;
}
