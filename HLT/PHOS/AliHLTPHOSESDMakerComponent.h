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

#ifndef ALIHLTPHOSESDMAKERCOMPONENT_H
#define ALIHLTPHOSESDMAKERCOMPONENT_H




/**
 * ESD maker component for PHOS HLT
 *
 * @file   AliHLTPHOSESDMakerComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  An ESD maker component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSESDMaker;
class AliHLTPHOSCaloClusterContainerStruct;
class AliESDEvent;

/**
 * @class AliHLTPHOSESDMakerComponent
 *
 * HLT component for making AliESDEvent from AliHLTPHOSCaloClusterDataStructs 
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSESDMakerComponent: public AliHLTPHOSProcessor
{
 public:

  /** Constructor */

  AliHLTPHOSESDMakerComponent();

  /** Destructor */
  virtual ~AliHLTPHOSESDMakerComponent();

  /** Copy constructor */  
 AliHLTPHOSESDMakerComponent(const AliHLTPHOSESDMakerComponent &) : 
  AliHLTPHOSProcessor(),
    fESDMakerPtr(0),
    fESDEventPtr(0)
    
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSESDMakerComponent & operator = (const AliHLTPHOSESDMakerComponent)
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
  using AliHLTPHOSProcessor::DoEvent;
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();
  
protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int Deinit();
 
private:

  /** Pointer to the ESD maker it self */
  AliHLTPHOSESDMaker* fESDMakerPtr;

  AliESDEvent* fESDEventPtr;

  /** interface function, see @ref AliHLTComponent for description */
  static const AliHLTComponentDataType fgkInputDataTypes[];
};


#endif
