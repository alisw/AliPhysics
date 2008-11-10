
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


#ifndef ALIHLTPHOSRAWANALYZERCOMPONENTV2_H
#define ALIHLTPHOSRAWANALYZERCOMPONENTV2_H


/**
 * Raw data analyzer component base class for PHOS HLT
 *
 * @file   AliHLTPHOSRawAnalyzerComponentv2.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A clusterizer component for PHOS HLT
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
class AliAltroDecoder;      
class AliAltroData;         
class AliAltroBunch;        
class AliHLTPHOSDigitMaker;
class AliHLTPHOSDigitContainerDataStruct;


/**
 * @class AliHLTPHOSRawAnalyzerComponentv2
 *
 * Class for decoding and extract the energy and timing information from the
 * channels in PHOS. 
 * 
 * The component has the following component arguments:
 * -offset [int]             offset applied before doing ZS
 * -bunchsizecut [int]       minimum size a bunch can have to be considered   
 * -minpeakposition [int]    minimum time of the peak for a signal to be considered
 * -maxpeakposition [int]    maximum time of the peak for a signal to be considered
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSRawAnalyzerComponentv2 : public AliHLTPHOSRcuProcessor
{
 public:

  /** Standard constructor */
  AliHLTPHOSRawAnalyzerComponentv2();

  /** Destructor */
  virtual ~AliHLTPHOSRawAnalyzerComponentv2();

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

  /** interface function, see @ref AliHLTComponent for description */
  using AliHLTPHOSRcuProcessor::DoEvent;

    /** interface function, see @ref AliHLTComponent for description */
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );  

  /** 
   * Do the real processing in the component 
   * @param iter is the pointer to the data blocks
   * @param outputPtr is the pointer to the output buffer
   * @param size is the available size for output
   * @param totSize is the total size used for output
   * @return the size output size used
   */
  virtual Int_t DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize); 


  /** Pointer to an analyzer object used for raw data anlysis */ 
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr;   //COMMENT

 private:

  /** Keep the copy constructor private since it should not be used */
  AliHLTPHOSRawAnalyzerComponentv2(const AliHLTPHOSRawAnalyzerComponentv2 & );

  /** Keep the assignement operator private since it should not be used */
  AliHLTPHOSRawAnalyzerComponentv2 & operator = (const AliHLTPHOSRawAnalyzerComponentv2 &);

  /** Mapping from harware address to geometrical address */
  AliHLTPHOSMapper *fMapperPtr;                       //COMMENT

  /** Pointer to object which may check the integrity of the data */
  AliHLTPHOSSanityInspector *fSanityInspectorPtr;     //COMMENT

  /** Pointer to the decoder used to unpack the altro payload */
  AliAltroDecoder *fDecoderPtr;                       //COMMENT

  /** Pointer to struct containing all data from one altro channel */
  AliAltroData    *fAltroDataPtr;                     //COMMENT

  /** Pointer to struct containing information about single bunches in one altro channel */
  AliAltroBunch   *fAltroBunchPtr;                    //COMMENT

  /** Describing which algorithm we are using */
  Short_t fAlgorithm;                                 //COMMENT

  /** The offset applied before ZS */
  Int_t fOffset;                                      //COMMENT

  /** The minimum length a bunch can have to be considered */
  Int_t fBunchSizeCut;                                //COMMENT

  /** The lowest position a peak can have to be considered */
  Int_t fMinPeakPosition;                             //COMMENT
  
  /** The maximum position a peak can have to be considered */
  Int_t fMaxPeakPosition;                             //COMMENT

};

#endif

