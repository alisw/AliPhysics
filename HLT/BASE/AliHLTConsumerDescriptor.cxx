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

/** @file   AliHLTConsumerDescriptor.cxx
    @author Matthias Richter
    @date   
    @brief  Helper class to describe a consumer component.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTConsumerDescriptor.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTConsumerDescriptor)

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor()
  :
  fpConsumer(NULL),
  fSegments()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fSegments.clear();
}

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor(AliHLTComponent* pConsumer)
  :
  fpConsumer(pConsumer),
  fSegments()
{
  // see header file for function documentation
  fSegments.clear();
}

AliHLTConsumerDescriptor::~AliHLTConsumerDescriptor()
{
  // see header file for function documentation
  if (fSegments.size()>0) {
    //HLTWarning("unreleased data segments found");
  }
}

int AliHLTConsumerDescriptor::SetActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  // see header file for function documentation
  int iResult=0;
  AliHLTDataBuffer::AliHLTDataSegment segment(offset, size);
  fSegments.push_back(segment);
  //HLTDebug("set active segment (%d:%d) for consumer %p", offset, size, this);
  return iResult;
}

int AliHLTConsumerDescriptor::CheckActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  // see header file for function documentation
  int iResult=0;
  if (fSegments.size()>0) {
    vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      if ((iResult=((*segment).fSegmentOffset==offset && (*segment).fSegmentSize==size))>0) {
	break;
      }
      segment++;
    }
  } else {
    //HLTWarning("no data segment active for consumer %p", this);
    iResult=-ENODATA;
  }
  return iResult;
}

int AliHLTConsumerDescriptor::ReleaseActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  // see header file for function documentation
  int iResult=0;
  if (fSegments.size()>0) {
    vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      if ((iResult=((*segment).fSegmentOffset==offset && (*segment).fSegmentSize==size))>0) {
	fSegments.erase(segment);
	break;
      }
      segment++;
    }
    if (iResult==0) {
      //HLTWarning("no data segment (%d:%d) active for consumer %p", offset, size, this);
      iResult=-ENOENT;
    }
  } else {
    //HLTWarning("no data segment active for consumer %p", this);
    iResult=-ENODATA;
  }
  return iResult;
}
