/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///
/// @file   AliHLTMUONAgent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  Implementation of the AliHLTMUONAgent class.
///

#include "AliHLTMUONAgent.h"
#include "AliHLTMUONRecHitsSource.h"
#include "AliHLTMUONTriggerRecordsSource.h"
#include "AliHLTMUONRootifierComponent.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONMansoTrackerFSMComponent.h"
#include "AliRunLoader.h"

// The single global instance of the dimuon HLT agent.
AliHLTMUONAgent AliHLTMUONAgent::fgkInstance;

ClassImp(AliHLTMUONAgent);


AliHLTMUONAgent::AliHLTMUONAgent() : AliHLTModuleAgent()
{
	///
	/// Default constructor.
	///
}

AliHLTMUONAgent::~AliHLTMUONAgent()
{
	///
	/// Default destructor.
	///
}

const char* AliHLTMUONAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						     AliRunLoader* /*runloader*/
	) const
{
	///
	/// Inherited from AliHLTModuleAgent.
	/// Returns the top processing chain configurations for local event
	/// reconstruction.
	/// @param rawReader  [in] AliRoot rawreader instance.
	/// @param runloader  [in] AliRoot runloader
	/// @return string containing the top configurations separated by blanks.
	///
	
	return "dhlt-simhits";
}

const char* AliHLTMUONAgent::GetRequiredComponentLibraries() const
{
	///
	/// Inherited from AliHLTModuleAgent.
	/// Returns a list of libraries which the configurations registered by
	/// this module agent depend on.
	/// @return list of component libraries as a blank-separated string.
	///
	
	return "libGeom.so libMinuit.so libEG.so libTreePlayer.so libXMLIO.so "
		"libVMC.so libESD.so libSTEER.so libGui.so libMUONraw.so libMUONgeometry.so "
		"libMUONmapping.so libMUONcalib.so libMUONbase.so libMUONsim.so libAliHLTMUON.so";
}


int AliHLTMUONAgent::CreateConfigurations(
		AliHLTConfigurationHandler* handler,
		AliRawReader* /*rawReader*/,
		AliRunLoader* /*runloader*/
	) const
{
	///
	/// Register all processing configurations belonging to the dimuon HLT
	/// library with the AliHLTConfigurationHandler.
	/// @param handler      the configuration handler
	/// @param rawReader  [in] AliRoot rawreader instance.
	/// @param runloader    AliRoot runloader
	/// @return Zero on success and error code if failed.
	///
	
	if (handler == NULL) return 0;
	handler->CreateConfiguration("dhlt-simhits", "DimuoRecHitsSource", NULL, "-simdata");
	return 0;
}


int AliHLTMUONAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
	///
	/// Registers all available components of this module.
	/// @param pHandler  [in] instance of the component handler.
	///
	
	if (pHandler == NULL) return -EINVAL;
	pHandler->AddComponent(new AliHLTMUONRecHitsSource);
	pHandler->AddComponent(new AliHLTMUONTriggerRecordsSource);
	pHandler->AddComponent(new AliHLTMUONRootifierComponent);
	pHandler->AddComponent(new AliHLTMUONHitReconstructorComponent);
	pHandler->AddComponent(new AliHLTMUONTriggerReconstructorComponent);
	pHandler->AddComponent(new AliHLTMUONMansoTrackerFSMComponent);
	return 0;
}

