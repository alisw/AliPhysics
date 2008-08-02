#ifndef ALIHLTMUONROOTIFIERCOMPONENT_H
#define ALIHLTMUONROOTIFIERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONRootifierComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Component for converting dHLT raw data into ROOT objects.
///

#include "AliHLTMUONProcessor.h"

/**
 * Converts dHLT raw data blocks into ROOT objects.
 */
class AliHLTMUONRootifierComponent : public AliHLTMUONProcessor
{
public:

	AliHLTMUONRootifierComponent();
	virtual ~AliHLTMUONRootifierComponent();
	
	virtual const char* GetComponentID();

	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	virtual AliHLTComponent* Spawn();

protected:

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	virtual bool IgnoreArgument(const char* arg) const;
	using AliHLTProcessor::DoEvent;
	
private:

	// Prevent copying of these objects.
	AliHLTMUONRootifierComponent(const AliHLTMUONRootifierComponent& /*object*/);
	AliHLTMUONRootifierComponent& operator = (const AliHLTMUONRootifierComponent& /*object*/);
	
	bool fWarnForUnexpecedBlock;  /// Flag indicating if we should log a warning if we got a block of an unexpected type.

	ClassDef(AliHLTMUONRootifierComponent, 0); // Converter component of dHLT raw data.
};

#endif // ALIHLTMUONROOTIFIERCOMPONENT_H
