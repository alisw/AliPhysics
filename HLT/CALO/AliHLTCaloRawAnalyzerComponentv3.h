
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


#ifndef ALIHLTCALORAWANALYZERCOMPONENTV3_H
#define ALIHLTCALORAWANALYZERCOMPONENTV3_H


/**
 * Raw data analyzer component base class for PHOS HLT
 *
 * @file   AliHLTCaloRawAnalyzerComponentv3.h
 * @author Oystein Djuvsland
 * @date   
 * @brief  A clusterizer component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTCaloRcuProcessor.h"


class AliCaloRawAnalyzer;
class AliHLTCaloRcuCellEnergyDataStruct;
class AliHLTCaloMapper;
class AliHLTCaloSanityInspector;
class AliHLTCaloDigitMaker;
class AliHLTCaloDigitContainerDataStruct;
class AliRawReaderMemory;
class AliAltroRawStreamV3;


/**
 * @class AliHLTCaloRawAnalyzerComponentv3
 * This the new and fast version of the component taking care of the decoding and energy and timing 
 * extraction of the raw data from PHOS.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b PhosRawAnalyzerv3 <br>
 * Library: \b libAliHLTCalo.so     <br>
 * Input Data Types: @ref <br>
 * Output Data Types: @ref AliHLTCaloDefinitions::fgkChannelDataType<br>
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



#include "AliHLTProcessor.h"
#include "AliHLTCaloDefinitions.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloRcuProcessor.h"

//#include "TObject.h"

class AliHLTCaloMapper;


class AliHLTCaloRawAnalyzerComponentv3 : public AliHLTCaloConstantsHandler, public AliHLTCaloRcuProcessor
{
 public:

  /** Constructor must be initialized to specific calorimeter */
  AliHLTCaloRawAnalyzerComponentv3(TString det);
  
  /** Destructor */
  virtual ~AliHLTCaloRawAnalyzerComponentv3();

  // virtual bool CheckInputDataType(const AliHLTComponentDataType &datatype) = 0;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int DoInit(int argc =0, const char** argv  = 0) ;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int DoDeinit();

  /** interface function, see @ref AliHLTComponent for description */
  virtual const char* GetComponentID() = 0;

  /** interface function, see @ref AliHLTComponent for description */
  //  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list) = 0; 

  /** interface function, see @ref AliHLTComponent for description */
  //  virtual AliHLTComponentDataType GetOutputDataType();
  virtual AliHLTComponentDataType GetOutputDataType() = 0;

  /** interface function, see @ref AliHLTComponent for description */
  //  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) = 0 ;

  /** interface function, see @ref AliHLTComponent for description */
  virtual AliHLTComponent* Spawn() = 0; 

 protected:

  //virtual bool CheckInputDataType(const AliHLTComponentDataType &datatype) = 0;
  /** interface function, see @ref AliHLTComponent for description */

  using AliHLTCaloRcuProcessor::DoEvent;

  /** interface function, see @ref AliHLTComponent for description */
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) = 0;  

  /** 
   * Do the real processing in the component 
   * @param iter is the pointer to the data blocks
   * @param outputPtr is the pointer to the output buffer
   * @param size is the available size for output
   * @param totSize is the total size used for output
   * @return the size output size used
   */
  virtual Int_t DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize); 


  // unsigned long fCaloEventCount;

  /** Pointer to an analyzer object used for raw data anlysis */ 
  AliCaloRawAnalyzer *fAnalyzerPtr;   //COMMENT

  //** Pointer to a mapper opbject */
  AliHLTCaloMapper *fMapperPtr;          //COMMENT

  virtual void InitMapping(const int specification ) = 0;

 private:

/** Keep default constructor private since it should not be used */
  AliHLTCaloRawAnalyzerComponentv3();

  /** Keep the copy constructor private since it should not be used */
  AliHLTCaloRawAnalyzerComponentv3(const AliHLTCaloRawAnalyzerComponentv3 & );

  /** Keep the assignement operator private since it should not be used */
  AliHLTCaloRawAnalyzerComponentv3 & operator = (const AliHLTCaloRawAnalyzerComponentv3 &);

  //virtual void InitMapping(const int specification ) = 0;
  
  /** Mapping from harware address to geometrical address */
  //  AliHLTCaloMapper *fMapperPtr;                       //!transient 



  /** Pointer to object which may check the integrity of the data */
  AliHLTCaloSanityInspector *fSanityInspectorPtr;     //!transient

  /** Pointer to the raw data reader which reads from memory */
  AliRawReaderMemory* fRawReaderMemoryPtr;            //!transient
  
  /** Pointer to the raw stream */
  AliAltroRawStreamV3* fAltroRawStreamPtr;              //!transient

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
  
  // AliHLTCaloMapper *fMapperPtr;

  /** Should we push the raw data when the channel is crazy? */
  Bool_t fDoPushBadRawData;                             //COMMENT
      
  /** Should we push all raw data (using the raw data writer) */
  Bool_t fDoPushRawData;                              //COMMENT

  class RawDataWriter 
  {
  public:
    RawDataWriter(AliHLTCaloConstants* cConst);
    ~RawDataWriter();
    //   void WriteChannelId(const UShort_t channeldid );
    void NewChannel( );
    void WriteBunchData(const UShort_t *bunchdata,  const int length,   const UInt_t starttimebin );
    void ResetBuffer();
    void SetChannelId( const UShort_t channeldid );
    //void CopyBufferToSharedMemory(UShort_t *memPtr, const int sizetotal, const int sizeused );
    int CopyBufferToSharedMemory(UShort_t *memPtr, const int sizetotal, const int sizeused );
    void NewEvent();
   
  private:
    
    //Default constructor, should not be used. 
    RawDataWriter();    
    RawDataWriter (const RawDataWriter  & );
    RawDataWriter & operator = (const RawDataWriter &);
    void Init();
    //    bool fIsFirstChannel;
    UShort_t* fRawDataBuffer;
    int fCurrentChannelSize;
    int fBufferIndex;
    int fBufferSize;
    UShort_t *fCurrentChannelIdPtr;
    UShort_t *fCurrentChannelSizePtr; 
    UShort_t *fCurrentChannelDataPtr; 
    int fTotalSize;
  };

  RawDataWriter *fRawDataWriter; 

  ClassDef(AliHLTCaloRawAnalyzerComponentv3, 1)

};

#endif

