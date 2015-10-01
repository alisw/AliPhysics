/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTZeroComponent.cxx
@author  David Rohr (drohr@cern.ch)
*/

#include "AliHLTZeroComponent.h"
#include "AliHLTErrorGuard.h"

ClassImp(AliHLTZeroComponent)

AliHLTZeroComponent::AliHLTZeroComponent() {}

void AliHLTZeroComponent::GetOCDBObjectDescription(TMap*) {}

AliHLTZeroComponent::~AliHLTZeroComponent() {
}

const Char_t* AliHLTZeroComponent::GetComponentID() { 
	return "Zero";
}

void AliHLTZeroComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
	list.push_back(kAliHLTAllDataTypes|kAliHLTDataOriginAny);
	list.push_back(kAliHLTDataTypeEOR|kAliHLTDataOriginAny);
}

AliHLTComponentDataType AliHLTZeroComponent::GetOutputDataType() {
	return kAliHLTDataTypeCustomTrigger|kAliHLTDataOriginHLT;
}

void AliHLTZeroComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
	constBase = 1;
	inputMultiplier = 0.;
}

AliHLTComponent* AliHLTZeroComponent::Spawn() {
	return new AliHLTZeroComponent;
}

// #################################################################################
Int_t AliHLTZeroComponent::DoInit(Int_t argc, const Char_t** argv)
{
	return 0;
}


// #################################################################################
Int_t AliHLTZeroComponent::DoDeinit() {
	return 0;
}

// #################################################################################
Int_t AliHLTZeroComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {

	return 0;
}


// #################################################################################
Int_t AliHLTZeroComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
	return(0);
}

// #################################################################################
Int_t AliHLTZeroComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
	ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
	return 0;
}

