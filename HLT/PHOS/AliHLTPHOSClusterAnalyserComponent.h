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

#ifndef ALIHLTPHOSCLUSTERANALYSERCOMPONENT_H
#define ALIHLTPHOSCLUSTERANALYSERCOMPONENT_H



/**
 * Cluster analyser component for PHOS HLT
 *
 * @file   AliHLTPHOSClusterAnalyserComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A cluster analyser component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSClusterAnalyser;

/**
 * @class AliHLTPHOSClusterAnalyserComponent
 *
 * Class for running cluster analysis for PHOS in HLT. It takes
 * reconstruction points as input, analyses them, and outputs the
 * rec points with updated information
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSClusterAnalyserComponent : public AliHLTPHOSProcessor
{
 public:

  /** Constructor */
  AliHLTPHOSClusterAnalyserComponent();

  /** Destructor */
  virtual ~AliHLTPHOSClusterAnalyserComponent();

  /** Copy constructor */
  /// AliHLTPHOSClusterAnalyserComponent(const AliHLTPHOSClusterAnalyserComponent &);
  
  /** Assignment */
  //AliHLTPHOSClusterAnalyserComponent & operator = (const AliHLTPHOSClusterAnalyserComponent);
  
  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  using AliHLTPHOSProcessor::DoEvent;
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();
  
protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int Deinit();

private:
  /** Copy constructor */
 AliHLTPHOSClusterAnalyserComponent(const AliHLTPHOSClusterAnalyserComponent &);
 /** Assignment */
 AliHLTPHOSClusterAnalyserComponent & operator = (const AliHLTPHOSClusterAnalyserComponent);
 
 AliHLTPHOSClusterAnalyser* fClusterAnalyserPtr;
 
 Bool_t fDoDeconvolution;
 Bool_t fDoCalculateMoments;
 
 static const AliHLTComponentDataType fgkInputDataTypes[];
   
};

#endif
