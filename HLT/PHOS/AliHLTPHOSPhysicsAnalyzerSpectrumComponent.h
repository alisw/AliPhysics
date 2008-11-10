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

#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT_H
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT_H

/**
 * Invariant mass spectrum component for PHOS HLT
 *
 * @file   AliHLTPHOSPhysicsAnalyzerSpectrumComponent.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  An invariant mass spectrum component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

// removed  PTH#include "AliHLTProcessor.h"
#include "AliHLTPHOSProcessor.h" // added by PTH
#include "AliHLTPHOSBase.h"

class TH1F;
class AliHLTPHOSPhysicsAnalyzerSpectrum;
class AliHLTPHOSPhysicsAnalyzerPeakFitter;
class AliHLTPHOSDefinitions;
class TFile;
class  AliHLTPHOSRecPointDataStruct;

/**
 * @class AliHLTPHOSClusterizerComponent
 *
 * Class for making an invariant mass spectrum histogram from 2 gammas.
 * Takes as input AliHLTPHOSRecPointContainerStruct
 * 
 *                         
 * @ingroup alihlt_phos
 */
// PTH class AliHLTPHOSPhysicsAnalyzerSpectrumComponent: public AliHLTPHOSBase, public AliHLTProcessor
class AliHLTPHOSPhysicsAnalyzerSpectrumComponent: public AliHLTPHOSProcessor // added by PTH
{
 public:

  /** Constructor */
  AliHLTPHOSPhysicsAnalyzerSpectrumComponent();

  /** Destructor */
  ~AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
  
  /** Copy constructor */
  AliHLTPHOSPhysicsAnalyzerSpectrumComponent(const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &);

  /** Assignment */
  AliHLTPHOSPhysicsAnalyzerSpectrumComponent & operator = (const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &)
  {
    // Assignment
    return *this;
  }

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */ 
  Int_t DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
		AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		std::vector<AliHLTComponentBlockData>&);
  /** interface function, see @ref AliHLTComponent for description */
  /*
  int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  */

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
  using AliHLTProcessor::DoEvent;
  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  Int_t Deinit();

  //  using AliHLTPHOSProcessor::DoEvent;
  
 private:
  
  /** Pointer to spectrum analyzer */
  AliHLTPHOSPhysicsAnalyzerSpectrum* fAnalyzerPtr;                //! transient 

  /** Pointer to peak fitter */
  AliHLTPHOSPhysicsAnalyzerPeakFitter* fPeakFitter;               //! transient 

  /** Pointer to histogram */
  TH1F* fRootHistPtr;                                             //! transient 

  /** Pointer to array of clusters */
  AliHLTPHOSRecPointDataStruct* fRecPointArrayPtr[10000];           //! transient 

  /** Interval for writing to disk */
  Int_t fWriteInterval;                                               //COMMENT

  /** Data types */                          
  static const AliHLTComponentDataType fgkInputDataTypes[];           // COMMENT

  /** Event count */
  static UInt_t fgCount;                                              //COMMENT
};

#endif
