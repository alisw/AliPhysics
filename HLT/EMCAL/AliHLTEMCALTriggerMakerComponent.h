/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTEMCALTRIGGERMAKERCOMPONENT_H
#define ALIHLTEMCALTRIGGERMAKERCOMPONENT_H

/**
 * @file   AliHLTEMCALTriggerMakerComponent.h
 * @author Markus Fasel
 * @date   Oct. 30, 2015
 * @brief  A trigger maker component for EMCAL HLT
 */

#include "AliHLTCaloProcessor.h"

class AliHLTEMCALGeometry;
class AliHLTEMCALTriggerMaker;
struct AliHLTCaloTriggerPatchContainerStruct;

/**
 * @class AliHLTEMCALTriggerMakerComponent
 * @brief HLT component for EMCAL/DCAL trigger patch finder
 */
class AliHLTEMCALTriggerMakerComponent : public AliHLTCaloProcessor {
public:
	AliHLTEMCALTriggerMakerComponent();
	virtual ~AliHLTEMCALTriggerMakerComponent();

	/**
	 *
	 * @param evtData
	 * @param blocks
	 * @param trigData
	 * @param outputPtr
	 * @param size
	 * @param outputBlocks
	 * @return
	 */
  int DoEvent ( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                std::vector<AliHLTComponentBlockData>& outputBlocks );

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize ( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();


protected:
  /** interface function, see @ref AliHLTComponent for description */
   int DoInit ( int argc, const char** argv );

   /** interface function, see @ref AliHLTComponent for description */
   virtual int Deinit();

   /** Initialising geometry for the HLT Trigger maker */
   void InitialiseGeometry();

   bool CheckInputDataType(const AliHLTComponentDataType &datatype);

private:
   /** Pointer to the trigger maker from cells */
   AliHLTEMCALTriggerMaker                    *fTriggerMakerPtrCells;       //! Transient
    /** Pointer to trigger maker from fastors */
   AliHLTEMCALTriggerMaker                    *fTriggerMakerPtrFastor;      //! Transient
   /** EMCAL geometry data */
   AliHLTEMCALGeometry                        *fGeometry;                   //! Transient

   /** Copy constructor,  not implemented */
   AliHLTEMCALTriggerMakerComponent(const AliHLTEMCALTriggerMakerComponent &);

   /** Assignment operator, not implemented */
   AliHLTEMCALTriggerMakerComponent & operator = (const AliHLTEMCALTriggerMakerComponent);

	ClassDef(AliHLTEMCALTriggerMakerComponent, 1);
};

#endif
