//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTPHOSBASELINEANALYZERCOMPONENT_H
#define ALIHLTPHOSBASELINEANALYZERCOMPONENT_H


/** 
 * Class does baseline analysis
 * 
 * @file   AliHLTPHOSBaselineAnalyzerComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A baseline analyzer for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSBaselineAnalyzer;
class TTree;


/**
 * @class AliHLTPHOSBaselineAnalyzerComponent
 * 
 * Class for running baseline analysis for PHOS in HLT. Takes raw data as input and creates 
 * objects of class AliHLTPHOSBaseline. Also it fills histograms of baselines as well as RMS
 * of the channels
 */ 

class AliHLTPHOSBaselineAnalyzerComponent : public AliHLTPHOSProcessor
{
public:

  /** Constructor */ 
  AliHLTPHOSBaselineAnalyzerComponent();

  /** Destructor */ 
  virtual ~AliHLTPHOSBaselineAnalyzerComponent();

  /** Copy constructor */
  AliHLTPHOSBaselineAnalyzerComponent(const AliHLTPHOSBaselineAnalyzerComponent &) : 
    AliHLTPHOSProcessor(),
    fBaselineAnalyzerPtr(0),
    fTreePtr(0),
    fBaselineArrayPtr(0),
    fEvCnt(0),
    fWriteInterval(100),
    fFillInterval(100),
    fFilename(0),
    fDirectory(0),
    fHistPath(0),
    fRunNb(0),
    fCalculateAll(false)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSBaselineAnalyzerComponent & operator = (const AliHLTPHOSBaselineAnalyzerComponent)
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
  using AliHLTPHOSProcessor::DoEvent;

  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
private:
  
  /** Calculate all values baselines, RMS etc. */
  void CalculateAll();

  /** Pointer to a baseline analyzer */ 
  AliHLTPHOSBaselineAnalyzer *fBaselineAnalyzerPtr; //! transient

  /** Pointer the a tree containing the TClonesArrays of AliHLTPHOSBaseline objects */ 
  TTree *fTreePtr;                                  //! transient

  /** TClonesArray of AliHLTPHOSBaseline objects */ 
  TClonesArray *fBaselineArrayPtr;                  //! transient

  /** Event count */
  UInt_t fEvCnt;                                    //COMMENT

  /** Number of events between each writing of files */
  UInt_t fWriteInterval;                            //COMMENT

  /** Number of events between each filling of the tree */ 
  UInt_t fFillInterval;                             //COMMENT

  /** The filename base */
  char *fFilename;                                  //! transient

  /** The directory to write the files */
  char* fDirectory;                                 //! transient

  /** The path for the histograms */
  char* fHistPath;                                  //! transient

  /** The run number */
  Int_t fRunNb;                                     //COMMENT

  /** If every value should be calculated, not only RMS */ 
  Bool_t fCalculateAll;                             //COMMMENT

  /** Interface variable */
  static const AliHLTComponentDataType fgkInputDataTypes[];     //COMMENT

};
#endif
