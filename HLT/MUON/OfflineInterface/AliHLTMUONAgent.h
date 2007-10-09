#ifndef ALIHLTMUONAGENT_H
#define ALIHLTMUONAGENT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONAgent.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  The HLT module agent for libAliHLTMUON.so which interfaces HLT
 *         components with offline.
 */

#include "AliHLTModuleAgent.h"
class AliRunLoader;

/**
 * This module agent handles dimuon HLT module registration and configurations
 * within the AliRoot framework.
 */
class AliHLTMUONAgent : public AliHLTModuleAgent
{
public:
	AliHLTMUONAgent();
	
	virtual ~AliHLTMUONAgent();

	/**
	 * Register all processing configurations belonging to the dimuon HLT
	 * library with the AliHLTConfigurationHandler.
	 * @param handler      the configuration handler
	 * @param runloader    AliRoot runloader
	 * @return Zero on success and error code if failed.
	 */
	int CreateConfigurations(
			AliHLTConfigurationHandler* handler,
			AliRawReader* rawReader=NULL,
			AliRunLoader* runloader = NULL
		) const;

	/**
	 * Returns the top processing chain configurations for local event
	 * reconstruction.
	 * @param runloader  [in] AliRoot runloader
	 * @return string containing the top configurations separated by blanks.
	 */
	const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
					    AliRunLoader* runloader = NULL) const;

	/**
	 * Component libraries which the configurations of this agent depend on.
	 * @return list of component libraries as a blank-separated string.
	 */
	const char* GetRequiredComponentLibraries() const;

	ClassDef(AliHLTMUONAgent, 1); // Dimuon HLT module agent which handles processing configurations.
};

#endif // ALIHLTMUONAGENT_H
