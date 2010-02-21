//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland, Svein Lindal                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALOHISTOCOMPONENT_H
#define ALIHLTCALOHISTOCOMPONENT_H

/**
 * 
 *
 * @file   AliHLTCaloHistoComponent.cxx
 * @author Svein Lindal
 * @date   
 * @brief  
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

#include "Rtypes.h"

class AliHLTPHOSPhysicsHistogramProducer;
class AliHLTCaloHistoCellEnergy;
class AliHLTCaloHistoClusterEnergy;
class AliHLTCaloHistoInvMass;
class AliHLTCaloHistoMatchedTracks;
class TRefArray;
class AliHLTCaloClusterReader;
/**
 * @class AliHLTPHOSHistogramProducerComponent
 *
 * 
 * @ingroup alihlt_phos
 */

class AliHLTCaloHistoComponent : public AliHLTProcessor
{

 public:

  /** Constructor */
  AliHLTCaloHistoComponent();
  /** Destructor */
  virtual ~AliHLTCaloHistoComponent();

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();


protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);
  /** interface function, see @ref AliHLTComponent for description */
  int DoDeinit();

  using AliHLTProcessor::DoEvent;

 private:

  /** Copy constructor prohibited*/  
  AliHLTCaloHistoComponent(const AliHLTCaloHistoComponent & );
  /** asssignment operator prohibited */
  AliHLTCaloHistoComponent& operator=(const AliHLTCaloHistoComponent&);

  TRefArray * fEmcalClustersArray;  //Array to contain EMCAL Clusters
  TRefArray * fPhosClustersArray;  //Array to contain PHOS Clusters

  Bool_t fDoPhos;       // Process PHOS data?
  Bool_t fDoEmcal;      // Process EMCAL data?
  Bool_t fDoCellEnergy; // make the cell energy histograms?
  Bool_t fDoClusterEnergy; // make the cluster energy histograms?
  Bool_t fDoInvariantMass; // make the invariant mass histograms?
  Bool_t fDoMatchedTracks; // make the matched tracks histograms?
  
  AliHLTCaloHistoCellEnergy *fPhosCellEnergyHistProducer; // cell energy histogram producer
  AliHLTCaloHistoCellEnergy *fEmcalCellEnergyHistProducer; // cell energy histogram producer

  AliHLTCaloHistoClusterEnergy *fPhosClusterEnergyHistProducer; // PHOS cluster energy histogram producer
  AliHLTCaloHistoClusterEnergy *fEmcalClusterEnergyHistProducer; // EMCAL cluster energy histogram producer

  AliHLTCaloHistoInvMass *fPhosInvariantMassHistProducer; // PHOS invariant mass histogram producer
  AliHLTCaloHistoInvMass *fEmcalInvariantMassHistProducer; // EMCAL insvariant mass histogram producer

  AliHLTCaloHistoMatchedTracks *fPhosMatchedTracksHistProducer; // PHOS matched tracks histogram producer
  AliHLTCaloHistoMatchedTracks *fEmcalMatchedTracksHistProducer; // EMCAL matched tracks histogram producer

  AliHLTCaloClusterReader * fClusterReader; //Class to read cluster data structs

  ClassDef(AliHLTCaloHistoComponent, 1);

};

#endif
