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
class TObjArray;
/**
 * @class AliHLTPHOSHistogramProducerComponent
 *
 * 
 * @ingroup alihlt_phos
 */

class AliHLTCaloHistoComponent : public AliHLTProcessor {

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

  /** function to get data content from blocks and pass it on to histogram produsers*/
  Int_t ProcessBlocks(const AliHLTComponentBlockData * pBlock, TObjArray * histoArray);
  
  AliHLTCaloClusterReader * fClusterReader; //!transient Class to read cluster data structs
  
  TRefArray * fEmcalClustersArray;  //!transient Array to contain EMCAL Clusters
  TRefArray * fPhosClustersArray;  //!transient Array to contain PHOS Clusters

  TObjArray * fPhosProducerArray; //!transient 
  TObjArray * fEmcalProducerArray; //!transient 
  
  TObjArray * fPhosHistogramArray; //!transient 
  TObjArray * fEmcalHistogramArray; //!transient 
  
  
  Bool_t fDoEmcal;  //Fill EMCAL histos?
  Bool_t fDoPhos;   //Fill PHOS histos?
  
  Bool_t fCutOnCentrality; // Cut on centrality on cluters with high energy
  Float_t fCentralityCut; //How large fraction of the energy do we want in the central tower to make the cut?
  Float_t fCentralityCutEnergy; //The minimum energy of the cluster to make the cut.

  ClassDef(AliHLTCaloHistoComponent, 0);

};

#endif
