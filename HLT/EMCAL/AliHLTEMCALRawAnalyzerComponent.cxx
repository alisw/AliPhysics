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

// Base class fro anlyzing EMCAL raww data
// Further documentation found in base class
// --------------
// --------------
// --------------
// --------------


#include "AliHLTCaloMapper.h"
#include "AliHLTEMCALRawAnalyzerComponent.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALTRURawDigitMaker.h"

#include "AliCaloConstants.h"

using namespace Algo;

AliHLTEMCALRawAnalyzerComponent::AliHLTEMCALRawAnalyzerComponent( fitAlgorithm algo ) : AliHLTCaloRawAnalyzerComponentv3("EMCAL", algo)
{
  //fDebug = true;
  fDebug = false;
  for(Int_t iddl = 0; iddl < 44; iddl++) { fMapperArray[iddl] = NULL; }
}


AliHLTEMCALRawAnalyzerComponent::~AliHLTEMCALRawAnalyzerComponent()
{

}


void 
AliHLTEMCALRawAnalyzerComponent::GetInputDataTypes( vector <AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType   | kAliHLTDataOriginEMCAL );
}


AliHLTComponentDataType
AliHLTEMCALRawAnalyzerComponent::GetOutputDataType()
{
  //comment
  return AliHLTEMCALDefinitions::fgkChannelDataType;
}

int
AliHLTEMCALRawAnalyzerComponent::DoDeinit()
{
  //comment
  for(Int_t iddl = 0; iddl < 44; iddl++){
    if(fMapperArray[iddl]) delete fMapperArray[iddl];
  }

  return AliHLTCaloRawAnalyzerComponentv3::DoDeinit();
}


void 
AliHLTEMCALRawAnalyzerComponent::InitMapping( const int specification )
{
  // Comment
  if(specification >= 44) return; // Do not handle STU DDLs
  if ( fMapperArray[specification] == 0 )
    {
      fMapperArray[specification] = new AliHLTEMCALMapper( specification );
      fMapperPtr = fMapperArray[specification];
    } else {
      fMapperPtr = fMapperArray[specification];
    }

  if(fMapperPtr->GetIsInitializedMapping() == false )
    {
      HLTError("%d:%d, ERROR, mapping not initialized ", __FILE__, __LINE__ );
      exit(-2);
    }
}

