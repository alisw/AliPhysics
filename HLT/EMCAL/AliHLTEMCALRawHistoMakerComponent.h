/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 * INFN, Laboratori Nazionali di Frascati                                 *
 * Primary Authors: Federico Ronchetti                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTEMCALRAWHISTOMAKERCOMPONENT_H
#define ALIHLTEMCALRAWHISTOMAKERCOMPONENT_H

/** @file   AliHLTEMCALRawHistoMakerComponent.h
    @author Federico Ronchetti
    @date   
    @brief  A histo maker component for EMCAL HLT
 */



#include "AliHLTCaloProcessor.h"
#include  "AliHLTEMCALDefinitions.h"
#include "TString.h"

class AliHLTEMCALRawHistoMaker;


class AliHLTEMCALRawHistoMakerComponent : public AliHLTCaloProcessor
{
public:

	/** Constructor */
	AliHLTEMCALRawHistoMakerComponent();

	/** Destructor */
	virtual ~AliHLTEMCALRawHistoMakerComponent();


	/** Assignment */
	AliHLTEMCALRawHistoMakerComponent & operator = (const AliHLTEMCALRawHistoMakerComponent)
	{
		//Assignment
		return *this;
	}

	/** interface function, see @ref AliHLTComponent for description */
	const char* GetComponentID();

	/** interface function, see @ref AliHLTComponent for description */
	void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

	/** interface function, see @ref AliHLTComponent for description */
	AliHLTComponentDataType GetOutputDataType();

	/** interface function, see @ref AliHLTComponent for description */
	void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	/** interface function, see @ref AliHLTComponent for description */
	int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
			std::vector<AliHLTComponentBlockData>& outputBlocks);

	/** interface function, see @ref AliHLTComponent for description */
	AliHLTComponent* Spawn();

protected:

	/** interface function, see @ref AliHLTComponent for description */

	int DoInit(int argc, const char** argv);
	int DoDeinit() {
		return(0);
	}

	using AliHLTCaloProcessor::DoEvent;

	/** interface function, see @ref AliHLTComponent for description */
	virtual int Deinit(); ////////// PTH WARNING

private:
	TString fRootFileName;
	int fPushFraction;
	int fLocalEventCount;
	int fBeVerbose;

	/** Pointer to the histo maker itself */
	AliHLTEMCALRawHistoMaker *fRawHistoMakerPtr;                    //! transient

	AliHLTEMCALRawHistoMakerComponent(const AliHLTEMCALRawHistoMakerComponent & );

};

#endif

