/*//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef AliHLTEMCALRawAnalyzerComponentTRU_H
#define AliHLTEMCALRawAnalyzerComponentTRU_H


/**
 * TRU Raw data analyzer component
 *
 * @file   AliHLTEMCALRawAnalyzerComponentTRU.h
 * @author Markus Fasel
 * @date
 * @brief Extraction of TRU digits for EMCAL HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


class AliCaloRawAnalyzer;
class AliHLTCaloRcuCellEnergyDataStruct;
class AliHLTCaloMapper;
class AliHLTCaloSanityInspector;
class AliHLTCaloDigitMaker;
class AliHLTCaloDigitContainerDataStruct;
class AliRawReaderMemory;
class AliCaloRawStreamV3;



#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloProcessor.h"
#include "AliCaloRawAnalyzer.h"


#include "AliCaloConstants.h"
using namespace Algo;

class AliHLTCaloMapper;
class AliHLTEMCALTRURawDigitMaker;


class AliHLTEMCALRawAnalyzerComponentTRU :  public AliHLTCaloProcessor, protected AliHLTCaloConstantsHandler
{
 public:

  /** Constructor must be initialized to specific calorimeter */
  AliHLTEMCALRawAnalyzerComponentTRU();
  virtual ~AliHLTEMCALRawAnalyzerComponentTRU();
  virtual int DoInit(int argc =0, const char** argv  = 0) ;
  virtual int DoDeinit();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();

 protected:
  bool CheckInputDataType(const AliHLTComponentDataType &datatype);
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
  virtual Int_t DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr,
         const AliHLTUInt32_t size, UInt_t& totSize);


 protected:
  void PrintDebugInfo();
  /** Pointer to L0 TRU handler */
  AliHLTEMCALTRURawDigitMaker *fTRUhandler;
  bool fDebug;    // Turn on to enable debug info

 private:
  AliHLTEMCALRawAnalyzerComponentTRU(const AliHLTEMCALRawAnalyzerComponentTRU & );
  AliHLTEMCALRawAnalyzerComponentTRU & operator = (const AliHLTEMCALRawAnalyzerComponentTRU &);

  /** Pointer to the raw data reader which reads from memory */
  //AliRawReaderMemory* fRawReaderMemoryPtr;            //!transient

  /** Pointer to the raw stream */
  //AliCaloRawStreamV3* fAltroRawStreamPtr;              //!transient


  ClassDef(AliHLTEMCALRawAnalyzerComponentTRU, 0)

};

#endif
