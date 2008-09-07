
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

class TH1D;
class TNtuple;
class AliHLTPHOSHistogramProducer;
class AliHLTPHOSCaloClusterContainerStruct;

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
  AliHLTPHOSHistogramProducerComponent(const AliHLTPHOSHistogramProducerComponent &) : 
    AliHLTPHOSProcessor(),
    fClusterEnergiesHistPtr(0),
    fMultiplicitiesHistPtr(0),
    fClusterNtuplePtr(0),
    fDoFillClusterEnergies(false),
    fDoFillMultiplicities(false),
    fDoFillNtuple(false),
    fHistogramProducerPtr(0)
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

  TH1D* fClusterEnergiesHistPtr; 
  TH1D* fMultiplicitiesHistPtr; 
  TNtuple* fClusterNtuplePtr; 

  bool fDoFillClusterEnergies;
  bool fDoFillMultiplicities;
  bool fDoFillNtuple;

  AliHLTPHOSHistogramProducer* fHistogramProducerPtr;

  /** interface function, see @ref AliHLTComponent for description */
  static const AliHLTComponentDataType fgkInputDataTypes[];     //COMMENT
};

#endif
