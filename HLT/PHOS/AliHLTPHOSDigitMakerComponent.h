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
#ifndef ALIHLTPHOSDIGITMAKERCOMPONENT_H
#define ALIHLTPHOSDIGITMAKERCOMPONENT_H

/** @file   AliHLTPHOSDigitMakerComponent.h
    @author Oystein Djuvsland
    @date   
    @brief  A digit maker component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSDigitMaker;
class TTree;
class TClonesArray;
class AliHLTPHOSDigitContainerDataStruct;


/**
 * @class AliHLTPHOSDigitMakerComponent
 *
 * Class runs AliHLTPHOSDigitMaker, creating digits from "raw data"
 *
 * The component has the following component arguments:
 * -threshold              threshold for creating a digit, gives software zero suppression
 * -presamples             number of presamples (not really necessary)
 *
 * @ingroup alihlt_phos
 */

class AliHLTPHOSDigitMakerComponent : public AliHLTPHOSProcessor
{
public:

  /** Constructor */
  AliHLTPHOSDigitMakerComponent();

  /** Destructor */ 
  virtual ~AliHLTPHOSDigitMakerComponent();

  /** Copy constructor */  
  AliHLTPHOSDigitMakerComponent(const AliHLTPHOSDigitMakerComponent &) : 
    AliHLTPHOSProcessor(),
    fDigitMakerPtr(0),
    fDigitContainerPtr(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSDigitMakerComponent & operator = (const AliHLTPHOSDigitMakerComponent)
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

  using AliHLTPHOSProcessor::DoEvent;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
private:

  /** Pointer to the digit maker it self */
  AliHLTPHOSDigitMaker *fDigitMakerPtr;                    //! transient

  /** The output of the component, digits in a container */
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerPtr;  //! transient

  /** Event count */
  //  UInt_t fEvtCnt; 
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
 
