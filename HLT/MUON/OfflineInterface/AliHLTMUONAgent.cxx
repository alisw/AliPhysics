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

/**
 * @file   AliHLTMUONAgent.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONAgent class.
 */

#include "AliHLTMUONAgent.h"
#include "AliRunLoader.h"

namespace
{
	// The single global instance of the dimuon HLT agent.
	AliHLTMUONAgent gAliHLTMUONAgent;
}

ClassImp(AliHLTMUONAgent);


AliHLTMUONAgent::AliHLTMUONAgent() : AliHLTModuleAgent()
{
}

AliHLTMUONAgent::~AliHLTMUONAgent()
{
}

const char* AliHLTMUONAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						     AliRunLoader* /*runloader*/
	) const
{
	return "dhlt-simhits";
}

const char* AliHLTMUONAgent::GetRequiredComponentLibraries() const
{
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
	if (handler == NULL) return 0;
	handler->CreateConfiguration("dhlt-simhits", "DimuoRecHitsSource", NULL, "-simdata");
	return 0;
}
