//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSMONITORTRIGGERCOMPONENT_H
#define ALIHLTPHOSMONITORTRIGGERCOMPONENT_H

/**
 * Monitor component
 *
 * @file   AliHLTPHOSMonitorTriggerComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A monitor trigger component
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloProcessor.h"

class AliHLTCaloClusterHeaderStruct;

/**
 * @class AliHLTPHOSMonitorTriggerComponent
 *
 * @ingroup alihlt_phos
 */

class AliHLTPHOSMonitorTriggerComponent: public AliHLTCaloProcessor
{
 public:

  /** Constructor */
  AliHLTPHOSMonitorTriggerComponent();

  /** Destructor */
  virtual ~AliHLTPHOSMonitorTriggerComponent();

  /** Copy constructor */  
  AliHLTPHOSMonitorTriggerComponent(const AliHLTPHOSMonitorTriggerComponent &) : 
    AliHLTCaloProcessor(),
    fCheckClusterEnergy(false),
    fCheckClusterMultiplicities(false),
    fClusterEnergyThreshold(1),
    fMultiplicityThreshold(5),
    fMultEnergyThreshold(0.5),
    fDigitMultiplicityThreshold(16),
    fMultDigitMultiplicityThreshold(9),
    fLowerCentrality(0),
    fUpperCentrality(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSMonitorTriggerComponent & operator = (const AliHLTPHOSMonitorTriggerComponent&)
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
  
  using  AliHLTCaloProcessor::DoEvent;
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks);
  // Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();
  
protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int Deinit();
  
  Bool_t CheckClusters(AliHLTCaloClusterHeaderStruct* clusterHeader);

private:
  
  Bool_t fCheckClusterEnergy; //COMMENT
  Bool_t fCheckClusterMultiplicities; //COMMENT

  Float_t fClusterEnergyThreshold; //COMMENT
  UInt_t fMultiplicityThreshold; //COMMENT
  Float_t fMultEnergyThreshold; //COMMENT
  UInt_t fDigitMultiplicityThreshold; //COMMENT
  UInt_t fMultDigitMultiplicityThreshold; //COMMENT

  Float_t fLowerCentrality; //COMMENT
  Float_t fUpperCentrality; //COMMENT
  

  /** interface function, see @ref AliHLTComponent for description */
  static const AliHLTComponentDataType fgkInputDataTypes[];     //COMMENT  
  

};

#endif
