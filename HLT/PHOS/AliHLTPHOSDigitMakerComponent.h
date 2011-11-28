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


#include "AliHLTPHOSDefinitions.h" // PTH
#include "AliHLTCaloProcessor.h"
#include "AliHLTCaloConstantsHandler.h"


class AliHLTCaloDigitMaker;
struct AliHLTCaloDigitContainerDataStruct;
class AliPHOSEmcBadChannelsMap;
class AliPHOSEmcCalibData;
/**
 * @class AliHLTPHOSDigitMakerComponent
 * 
 * Creates the digit used for the clusterizer. Digits are equivalent to the ones in 
 * offline reconstruction
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b PhosDigitMaker <br>
 * Library: \b libAliHLTPHOS.so     <br>
 * Input Data Types: @ref AliHLTPHOSDefinitions::fkgChannelDataType<br>
 * Output Data Types: @ref AliHLTPHOSDefinitions::fgkDigitDataType<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li No mandatory arguments for component                           <br>
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -lowgainfactor      <i> value </i> <br>
 *      sets a global low gain factor 
 * \li -highgainfactor <i> value </i> <br>
 *      sets a global high gain factor
 * \li -reverseorder <br>
 *      set if one expects the low gain channels to come before the high gain ones
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li No configuration arguments 
 *
 * <h2>Default CDB entries:</h2>
 * \li No CDB entry yet, will come.
 *
 * <h2>Performance:</h2>
 * Pretty good
 *
 * <h2>Memory consumption:</h2>
 * Pretty low
 *
 * <h2>Output size:</h2>
 * Depends on the event...
 *
 * More detailed description. (Soon)
 *
 * @ingroup alihlt_phos
 */ 

class AliHLTPHOSDigitMakerComponent : public AliHLTCaloProcessor, public AliHLTCaloConstantsHandler
{
public:

  /** Constructor */
  AliHLTPHOSDigitMakerComponent();

  /** Destructor */ 
  virtual ~AliHLTPHOSDigitMakerComponent();
  
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

  using AliHLTCaloProcessor::DoEvent;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor

  /** Get bad channel map from CDB */
  virtual int GetBCMFromCDB();
  
  /** Get the ADC <-> Energy (GeV) gain factors */
  virtual int GetGainsFromCDB();
  
  
private:
   
   /** Copy constructor, prohibited */  
  AliHLTPHOSDigitMakerComponent(const AliHLTPHOSDigitMakerComponent &); 
  
  /** Assignment operator, prohibited */
  AliHLTPHOSDigitMakerComponent & operator = (const AliHLTPHOSDigitMakerComponent);
  
  /** Pointer to the digit maker itself */
  AliHLTCaloDigitMaker *fDigitMakerPtr;                    //! transient

  /** The output of the component, digits in a container */
  AliHLTCaloDigitContainerDataStruct *fDigitContainerPtr;  //! transient
  
  /** Bad channel map */
  AliPHOSEmcBadChannelsMap *fBadChannelMap; //! transient
   
//  /** Temporary holder for bad channel map */
  //Bool_t ***fBadChannelMap; //! transient

  /** Calibration data */
  AliPHOSEmcCalibData *fCalibData; //! transient
   

  /** Is the bad map initialised? */
  Bool_t fBCMInitialised; //! transient

   /** Are the gains initialised? */
  Bool_t fGainsInitialised; //! transient
  
  ClassDef(AliHLTPHOSDigitMakerComponent, 0);

};
#endif
 
