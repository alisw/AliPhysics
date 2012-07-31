// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt, Sergey Gorbunov                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTITSDigitPublisherComponent.cxx
    @author Kenneth Aamodt, Sergey Gorbunov
    @date   
    @brief  Component to run offline clusterfinders
*/

#include <TSystem.h>
#include <TROOT.h>
#include "AliHLTITSDigitPublisherComponent.h" 
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliGeomManager.h"
#include "AliITSInitGeometry.h"
#include "AliITSLoader.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSDigitPublisherComponent);

AliHLTITSDigitPublisherComponent::AliHLTITSDigitPublisherComponent()
  :
  fRunLoader(NULL),
  fITSLoader(NULL),
  fNumberOfEvents(0),
  fEventNumber(0),
  tD(NULL)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTITSDigitPublisherComponent::~AliHLTITSDigitPublisherComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSDigitPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "ITSDigitPublisher";
}

void AliHLTITSDigitPublisherComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  //list.push_back( ???? | ???? );
}

AliHLTComponentDataType AliHLTITSDigitPublisherComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeAliTreeD|kAliHLTDataOriginITS;
}

void AliHLTITSDigitPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation
  constBase = 20000;
  inputMultiplier = 1000;
}

AliHLTComponent* AliHLTITSDigitPublisherComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSDigitPublisherComponent();
}
	
Int_t AliHLTITSDigitPublisherComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation
  
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }

  fRunLoader = GetRunLoader();//AliRunLoader::Open("galice.root");
  if(!fRunLoader){
    HLTFatal("No RunLoader found");
    return -1;
  }
  fITSLoader = (AliITSLoader *)(fRunLoader->GetLoader("ITSLoader"));
  if(!fITSLoader){
    HLTFatal("No ITS RunLoader found");
    return -1;
  }
  fNumberOfEvents = fRunLoader->GetNumberOfEvents();
 
  return 0;
}

Int_t AliHLTITSDigitPublisherComponent::DoDeinit() {
  // see header file for class documentation

  return 0;
}

Int_t AliHLTITSDigitPublisherComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/,AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  if (!IsDataEvent()) return 0;

  fRunLoader->GetEvent(fEventNumber);
  fITSLoader->LoadDigits("read");
  tD = fITSLoader->TreeD();
  if(!tD){
    HLTFatal("No Digit Tree found");
    return -1;
  } 
  //tD->GetEntry(fEventNumber);
  
  PushBack((TObject*)tD,kAliHLTDataTypeAliTreeD|kAliHLTDataOriginITS,0x00000000);
  fEventNumber++;
  return 0;
}

