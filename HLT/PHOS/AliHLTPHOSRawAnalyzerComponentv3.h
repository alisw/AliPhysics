
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
 * @brief  A raw analyzer component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSRcuProcessor.h"



class AliHLTPHOSRawAnalyzer;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSMapper;
class AliHLTPHOSSanityInspector;
class AliHLTPHOSDigitMaker;
class AliHLTPHOSDigitContainerDataStruct;
class AliRawReaderMemory;
class AliAltroRawStreamV3;
class AliHLTPHOSChannelDataStruct;

//class RawDataWriter;


/**
 * @class AliHLTPHOSRawAnalyzerComponentv3
 * @ingroup alihlt_phos
 */ 


class AliHLTPHOSRawAnalyzerComponentv3 : public AliHLTPHOSRcuProcessor
{
 public:

  /** Standard constructor */
  AliHLTPHOSRawAnalyzerComponentv3();

  /** Destructor */
  virtual ~AliHLTPHOSRawAnalyzerComponentv3();

  /** interface function, see @ref AliHLTComponent for description */
  virtual int DoInit(int argc =0, const char** argv  = 0);

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit();

  /** interface function, see @ref AliHLTComponent for description */
  virtual const char* GetComponentID() = 0;

  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  virtual AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  virtual AliHLTComponent* Spawn() = 0; 

 protected:

  /** 
   * Check for correct input data type (raw data from PHOS) 
   * @datatype is the data type specifier
   * @return true if the data type is correct
   */
  virtual bool CheckInputDataType(const AliHLTComponentDataType &datatype);

 private:

  /** Keep the copy constructor private since it should not be used */
  AliHLTPHOSRawAnalyzerComponentv3(const AliHLTPHOSRawAnalyzerComponentv3 & );

  /** Keep the assignement operator private since it should not be used */
  AliHLTPHOSRawAnalyzerComponentv3 & operator = (const AliHLTPHOSRawAnalyzerComponentv3 &);
  
  /** 
   * Initialise the mapping according to the specification
   * @specification is the specification provided by the HLT framework
   */
  virtual void InitMapping(const int specification);


};

#endif

