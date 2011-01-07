#ifndef ALIHLTMUONAGENT_H
#define ALIHLTMUONAGENT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// @file   AliHLTMUONAgent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 May 2007
/// @brief  The HLT module agent for libAliHLTMUON.so which interfaces HLT components with offline.
///

#include "AliHLTModuleAgent.h"
class AliRunLoader;
class AliHLTOUTHandlerChain;
class AliHLTOUTHandlerIgnore;

/**
 * This module agent handles dimuon HLT module registration and configurations
 * within the AliRoot framework.
 */
class AliHLTMUONAgent : public AliHLTModuleAgent
{
public:
	AliHLTMUONAgent();
	
	virtual ~AliHLTMUONAgent();

	// The following methods are all inherited from AliHLTModuleAgent:
	
	/**
	 * Register all processing configurations belonging to the dimuon HLT
	 * library with the AliHLTConfigurationHandler.
	 * @param handler      the configuration handler
	 * @param runloader    AliRoot runloader
	 * @return Zero on success and error code if failed.
	 */
	virtual int CreateConfigurations(
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
	virtual const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
					    AliRunLoader* runloader = NULL) const;

	/**
	 * Component libraries which the configurations of this agent depend on.
	 * @return list of component libraries as a blank-separated string.
	 */
	virtual const char* GetRequiredComponentLibraries() const;
	
	/**
	 * Registers all available components of this module.
	 * @param pHandler  [in] instance of the component handler.
	 */
	virtual int RegisterComponents(AliHLTComponentHandler* pHandler) const;
	
	/**
	 * Get handler decription for dHLT data in the HLTOUT data stream.
	 * @param dt        [in] data type of the block
	 * @param spec      [in] specification of the block
	 * @param desc      [out] handler description
	 * @return 1 if the agent can provide a handler, 0 if not.
	 */
	virtual int GetHandlerDescription(
			AliHLTComponentDataType dt,
			AliHLTUInt32_t spec,
			AliHLTOUTHandlerDesc& desc
		) const;
	
	/**
	 * Get specific handler for dHLT data in the HLTOUT data stream.
	 * @param dt        [in] data type of the block
	 * @param spec      [in] specification of the block
	 * @return pointer to handler
	 */
	virtual AliHLTOUTHandler* GetOutputHandler(
			AliHLTComponentDataType dt, AliHLTUInt32_t spec
		);
	
	/**
	 * Delete an HLTOUT handler.
	 * @param pInstance      pointer to handler
	 */
	virtual int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

	/**
	 * \returns true if a MUON module was added to gAlice.
	 */
	static bool IsMuonModuleLoaded();
	
	/**
	 * \returns the flag indicating if the dHLT-rootify-and-dump chain should
	 *    be run as a HLTOUT handler during reconstruction. The default is not
	 *    to run this chain. The chain is used to convert HLT raw data blocks
	 *    into ROOT objects, typically useful for testing and debugging.
	 */
	static bool RunRootifyChain() { return fgRunRootifyChain; }
	
	/**
	 * Sets the flag indicating if the dHLT-rootify-and-dump chain should be
	 * run during reconstruction.
	 */
	static void RunRootifyChain(bool value) { fgRunRootifyChain = value; }
	
private:
	// The following instance is used for automatic agent and component registration.
	static AliHLTMUONAgent fgkInstance;  ///< The single global instance of the dimuon HLT agent.
	
	static AliHLTOUTHandlerChain  fgkESDMakerChain;  ///< Chain handler for converting dHLT raw data to ESD format.
	static AliHLTOUTHandlerChain  fgkRootifyDumpChain;  ///< Chain handler for converting dHLT raw data to ROOT objects and dumping to file.
	static AliHLTOUTHandlerIgnore fgkDataIgnoreHandler;  ///< HLTOUT handler for ignoring data blocks.

	static Int_t fgMuonModuleLoaded; ///< Cached flag for indicating if the MUON module was loaded for a simulation.
	static bool fgRunRootifyChain; // Indicates if the dHLT-rootify-and-dump chain should be run.

	ClassDef(AliHLTMUONAgent, 0); // Dimuon HLT module agent which handles processing configurations.
};

#endif // ALIHLTMUONAGENT_H
