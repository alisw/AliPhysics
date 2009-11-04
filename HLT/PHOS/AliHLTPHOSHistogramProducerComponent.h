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

#ifndef ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H



/**
 * 
 *
 * @file   AliHLTPHOSHistogramProducerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSPhysicsHistogramProducer;

class AliHLTPHOSHistoProdCellEnergy;

class AliHLTPHOSHistoProdClusterEnergy;

class AliHLTPHOSHistoProdInvMass;

class AliHLTPHOSHistoProdMatchedTracks;

/**
 * @class AliHLTPHOSHistogramProducerComponent
 *
 * 
 * @ingroup alihlt_phos
 */
class AliHLTPHOSHistogramProducerComponent: public AliHLTPHOSProcessor
{
 public:

  /** Constructor */
  AliHLTPHOSHistogramProducerComponent();

  /** Destructor */
  virtual ~AliHLTPHOSHistogramProducerComponent();

  /** Copy constructor */  
  AliHLTPHOSHistogramProducerComponent(const AliHLTPHOSHistogramProducerComponent & ) : 
    AliHLTPHOSProcessor(),
    fPhysicsHistogramProducerPtr(0),
    fPushModulo(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSHistogramProducerComponent & operator = (const AliHLTPHOSHistogramProducerComponent)
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
  
  using  AliHLTPHOSProcessor::DoEvent;

  int DoEvent(const AliHLTComponentEventData& evtData, 
	      AliHLTComponentTriggerData& trigData);
    
  // Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();
  
protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int Deinit();

 private:

  AliHLTPHOSPhysicsHistogramProducer* fPhysicsHistogramProducerPtr;
  UInt_t fPushModulo;
  
  
  Bool_t fCellEnergy; // make the cell energy histograms?
  Bool_t fClusterEnergy; // make the cluster energy histograms?
  Bool_t fInvariantMass; // make the invariant mass histograms?
  Bool_t fMatchedTracks; // make the matched tracks histograms?
  
  AliHLTPHOSHistoProdCellEnergy *fCellEnergyHistProducer; // cell energy histogram producer

  AliHLTPHOSHistoProdClusterEnergy *fClusterEnergyHistProducer; // cluster energy histogram producer

  AliHLTPHOSHistoProdInvMass *fInvariantMassHistProducer; // invariant mass histogram producer

  AliHLTPHOSHistoProdMatchedTracks *fMatchedTracksHistProducer; // matched tracks histogram producer

};

#endif
