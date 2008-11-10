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



#ifndef ALIHLTPHOSTREEMAKERCOMPONENT_H
#define ALIHLTPHOSTREEMAKERCOMPONENT_H

/**
 * Tree maker component for PHOS HLT
 *
 * @file   AliHLTPHOSTreeMakerComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A tree maker component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSTreeMaker;
class TTree;


/**
 * @class AliHLTPHOSTreeMakerComponent
 *
 * Class for making trees of digits in PHOS HLT
 * Takes as input structs of type AliHLTPHOSDigitContainerStruct, and makes a
 * tree of TClonesArrays of objects of type AliHLTPHOSDigit 
 *
 * The component has the following component arguments:
 * -path                   File name base
 * -writeinterval          How often to write to disk, in events
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSTreeMakerComponent : public AliHLTPHOSProcessor
{
 public:

  /** Constructor */
  AliHLTPHOSTreeMakerComponent();

  /** Destructor */
  virtual ~AliHLTPHOSTreeMakerComponent();
    
  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();
  
  /** interface function, see @ref AliHLTComponent for description */
  void  GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
  
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  
  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputmultiplier);

/*
  int DoEvent(const AliHLTComponentEventData&,
	      AliHLTComponentTriggerData&);
  */

  /** interface function, see @ref AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t&size,
       std::vector<AliHLTComponentBlockData>& outputBlocks);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

  /** Write the tree */
  void Write(); 

  /** Reset the trees */
  void ResetTrees();
   

 protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  using AliHLTPHOSProcessor::DoEvent;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
 private:  

  /** Pointer to a AliHLTPHOSTreeMaker object */
  AliHLTPHOSTreeMaker *fTreeMakerPtr;                           //! transient

  /** Pointer to the tree */
  TTree *fDigitTreePtr;                                         //! transient

  /** Event count */
  UInt_t fEvtCount;                                             //COMMENT

  /** Write interval */
  UInt_t fWriteInterval;                                        //COMMENT

  /** Run number */
  //UInt_t fRunNb;                                                //COMMENT

  /** Path to the directory where to write the tree */
  char *fDirectory;                                             //! transient

  /** Framework variable */
  static const AliHLTComponentDataType fgkInputDataTypes[];     //COMMENT
  
};
#endif
  
