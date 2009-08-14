 
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

#ifndef ALIHLTPHOSESDCALOCLUSTERWRITERCOMPONENT_H
#define ALIHLTPHOSESDCALOCLUSTERWRITERCOMPONENT_H

/**
 * ESD Calo cluster writer component for PHOS HLT
 *
 * @file   AliHLTPHOSESDCaloClusterWriterComponent.h
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

class TClonesArray;
class TTree;
class TFile;
/**
 * @class AliHLTPHOSESDCaloClusterWriterComponent
 *
 * HLT component for making AliESDEvent from AliHLTPHOSCaloClusterDataStructs 
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSESDCaloClusterWriterComponent: public AliHLTPHOSProcessor
{
 public:

  /** Constructor */

  AliHLTPHOSESDCaloClusterWriterComponent();

  /** Destructor */
  virtual ~AliHLTPHOSESDCaloClusterWriterComponent();

  /** Copy constructor */  
  AliHLTPHOSESDCaloClusterWriterComponent(const AliHLTPHOSESDCaloClusterWriterComponent &) : 
    AliHLTPHOSProcessor(),
    fOutfile(0),
    fOutfileName(0),
    fWriteModulo(1000),
    fESDCaloClusterTreePtr(0),
    fESDCaloClustersPtr(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSESDCaloClusterWriterComponent & operator = (const AliHLTPHOSESDCaloClusterWriterComponent)
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
  
  /** Function for writing the tree containing the clusters */
  int WriteTree();

protected:

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int Deinit();
 
private:


  /** The file to which we will write */
  TFile* fOutfile;

  /** The filename */
  char* fOutfileName;

  /** Write modulo */
  UInt_t fWriteModulo;

  /** Pointer to the ESD calo cluster tree*/
  TTree* fESDCaloClusterTreePtr; //! transient

  /** Pointer to the array of clusters */
  TClonesArray* fESDCaloClustersPtr; //! transient


};


#endif
