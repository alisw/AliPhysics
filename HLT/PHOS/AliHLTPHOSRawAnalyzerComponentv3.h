
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille, Oystein Djuvsland                   *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTPHOSRAWANALYZERCOMPONENTV3_H
#define ALIHLTPHOSRAWANALYZERCOMPONENTV3_H


/**
 * Raw data analyzer component base class for PHOS HLT
 *
 * @file   AliHLTPHOSRawAnalyzerComponentv3.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A clusterizer component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <AliHLTCaloRawAnalyzerComponentv3.h>

class AliHLTPHOSRawAnalyzer;

/**
 * @class AliHLTPHOSRawAnalyzerComponentv3
 * This the new and fast version of the component taking care of the decoding and energy and timing 
 * extraction of the raw data from PHOS.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b PhosRawAnalyzerv3 <br>
 * Library: \b libAliHLTPHOS.so     <br>
 * Input Data Types: @ref <br>
 * Output Data Types: @ref AliHLTPHOSDefinitions::fgkChannelDataType<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li No mandatory arguments for component                           <br>
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -offset      <i> value </i> <br>
 *      gives the offset added to the data during zero suppression (default value: 0)
 * \li -bunchsizecut <i> value </i> <br>
 *      minimum number of samples a bunch must contain to be considered  (default value: 0)
 * \li -minpeakposition <i> value </i> <br>
 *      cut on minimum postion of the peak in the bunch (defaul value: 0)
 * \li -maxpeakposition <i> value </i> <br>
 *      cut on maximum postion of the peak in the bunch (defaul value: 100)
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li No configuration arguments 
 *
 * <h2>Default CDB entries:</h2>
 * \li No CDB entry yet, will come.
 *
 * <h2>Performance:</h2>
 * Pretty good (~ 3 kHz), depends on amount of data...
 *
 * <h2>Memory consumption:</h2>
 * Depends on the amount of data, but pretty godd
 *
 * <h2>Output size:</h2>
 * Depends on the amount of data...
 *
 * More detailed description. (Soon)
 *
 * @ingroup alihlt_phos
 */ 


class AliHLTPHOSRawAnalyzerComponentv3 : public AliHLTCaloRawAnalyzerComponentv3
{
 public:
  AliHLTPHOSRawAnalyzerComponentv3(); /** Standard constructor */
  virtual ~AliHLTPHOSRawAnalyzerComponentv3();
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);/** interface function, see @ref AliHLTComponent for description */
  virtual AliHLTComponentDataType GetOutputDataType();/** interface function, see @ref AliHLTComponent for description */
  
  virtual const char* GetComponentID() = 0; 
  virtual AliHLTComponent* Spawn() = 0; /** interface function, see @ref AliHLTComponent for description */
 
 protected:
  virtual void InitMapping(const int specification);
  
 private:
  AliHLTPHOSRawAnalyzerComponentv3(const AliHLTPHOSRawAnalyzerComponentv3 & );
  AliHLTPHOSRawAnalyzerComponentv3 & operator = (const AliHLTPHOSRawAnalyzerComponentv3 &);
  
};

#endif

