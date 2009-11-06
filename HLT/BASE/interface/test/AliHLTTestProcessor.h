//-*- Mode: C++ -*-
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

/** @file   AliHLTTestProcessor.h
    @author Matthias Richter
    @date   
    @brief  A processor used for testing the external wrapper interface
 */

#include <vector>
#include "AliHLTProcessor.h"
#include "AliHLTTest.h"

class AliHLTTestProcessor : public AliHLTProcessor, public AliHLTTest
{
public:
  AliHLTTestProcessor();
  ~AliHLTTestProcessor();

  const char* GetComponentID() {return "TestProcessor";}

  void GetInputDataTypes( vector<AliHLTComponentDataType>& list)
  {list.push_back(kAliHLTAnyDataType);}

  AliHLTComponentDataType GetOutputDataType() {return kAliHLTAnyDataType;}

  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
  {constBase=100; inputMultiplier=2.0;}
  
  AliHLTComponent* Spawn() {return new AliHLTTestProcessor;}

  bool CheckRunNo(unsigned runNo) const;
  bool CheckChainId(const char* chainId) const;
  bool CheckDataType(const char* id, const char* origin) const;
  bool CheckMagneticField(float bz);

private:
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/);

  enum {
    kCreated =0,
    kInitialized,
    kProcessing,
  };

  int fState; //!transient

  vector<AliHLTComponentDataType> fProcessedTypes; //!transient
};
