//-*- Mode: C++ -*-
// $Id: AliHLTZeroComponent $

#ifndef ALIHLTZEROCOMPONENT_H
#define ALIHLTZEROCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTZeroComponent.h
@author  David Rohr (drohr@cern.ch)
*/

#include "AliHLTProcessor.h"
#include "AliHLTAsyncMemberProcessor.h"

class TList;

class AliESDVZERO;
class AliESDtrackCuts;
class AliHLTCTPData;
class AliHLTMultiplicityCorrelations;
class AliHLTGlobalTriggerDecision;
class AliAnalysisManager;
class AliHLTTestInputHandler;

/**
* @class AliHLTZeroComponent
*/
class AliHLTZeroComponent : public AliHLTProcessor {
public:

	/*
	* ---------------------------------------------------------------------------------
	*                            Constructor / Destructor
	* ---------------------------------------------------------------------------------
	*/

	/** constructor */
	AliHLTZeroComponent();

	/** destructor */
	virtual ~AliHLTZeroComponent();

	/*
	* ---------------------------------------------------------------------------------
	* Public functions to implement AliHLTComponent's interface.
	* These functions are required for the registration process
	* ---------------------------------------------------------------------------------
	*/

	/** interface function, see @ref AliHLTComponent for description */
	const Char_t* GetComponentID();

	/** interface function, see @ref AliHLTComponent for description */
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

	/** interface function, see @ref AliHLTComponent for description */
	AliHLTComponentDataType GetOutputDataType();

	/** interface function, see @ref AliHLTComponent for description */
	void GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier );

	/** interface function, see @ref AliHLTComponent for description */
	void GetOCDBObjectDescription( TMap* const targetMap);

	/** interface function, see @ref AliHLTComponent for description */
	AliHLTComponent* Spawn();

protected:

	/*
	* ---------------------------------------------------------------------------------
	* Protected functions to implement AliHLTComponent's interface.
	* These functions provide initialization as well as the actual processing
	* capabilities of the component. 
	* ---------------------------------------------------------------------------------
	*/

	// AliHLTComponent interface functions

	/** interface function, see @ref AliHLTComponent for description */
	Int_t DoInit( Int_t /*argc*/, const Char_t** /*argv*/ );

	/** interface function, see @ref AliHLTComponent for description */
	Int_t DoDeinit();

	/** interface function, see @ref AliHLTComponent for description */
	Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

	using AliHLTProcessor::DoEvent;


	/** interface function, see @ref AliHLTComponent for description */
	Int_t Reconfigure(const Char_t* cdbEntry, const Char_t* chainId);

	/** interface function, see @ref AliHLTComponent for description */
	Int_t ReadPreprocessorValues(const Char_t* modules);

	///////////////////////////////////////////////////////////////////////////////////

private:

	/*
	* ---------------------------------------------------------------------------------
	* Private functions to implement AliHLTComponent's interface.
	* These functions provide initialization as well as the actual processing
	* capabilities of the component. 
	* ---------------------------------------------------------------------------------
	*/

	/** copy constructor prohibited */
	AliHLTZeroComponent(const AliHLTZeroComponent&);
	/** assignment operator prohibited */
	AliHLTZeroComponent& operator=(const AliHLTZeroComponent&);

	ClassDef(AliHLTZeroComponent, 0)
};
#endif
