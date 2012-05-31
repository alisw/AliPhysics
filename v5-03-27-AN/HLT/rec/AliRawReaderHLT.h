//-*- Mode: C++ -*-
// $Id$

#ifndef ALIRAWREADERHLT_H
#define ALIRAWREADERHLT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliRawReaderHLT.h
/// @author Matthias Richter
/// @date   
/// @brief  AliRawReader implementation which replaces original input of
///         detectors with the appropriate HLT output.

#include "AliHLTDataTypes.h"
#include "AliRawReader.h"      // RAW, base class
#include "TString.h"
#include <vector>

class AliHLTOUT;
class AliHLTOUTHandler;
class AliHLTPluginBase;

/**
 * @class AliRawReaderHLT
 * A specific AliRawReader for detector input replacement by HLTOUT data blocks.
 *
 * HLT components can produce output data in the detector ddl raw format.
 * Data blocks of this format can be fed into the normal detector reconstruction
 * without changes in the actual reconstruction code by means of the
 * AliRawReaderHLT implementation of the AliRawReader.
 *
 * @section sec_alirawreaderhlt_concept Conceptual design
 * The AliRawReader provides an abstract interface to the ddl raw data. All
 * reconstruction code uses this interface to access the data.
 * HLT components can send their data in the original ddl raw format. The only
 * difference of such data blocks is the location since they are shipped as
 * part of the HLTOUT data stream. The AliRawReaderHLT provides redirection of
 * those data blocks.
 *
 * The AliRawReaderHLT needs the original AliRawReader in order to get the
 * data. Furthermore, a string containing the detector specification defines
 * which data should be read from the HLT stream and which from the original
 * reader.
 *
 * @note An HLTOUT handler must be available for the HLTOUT data blocks to
 * be redirected. Please read @ref sec_alirawreaderhlt_module carefully.
 *
 * @section sec_alirawreaderhlt_usage   Selection of the HLTOUT data stream
 * The input data of a detector can be replaced by the corresponding HLT
 * data by calling the <tt>AliReconstruction::SetUseHLTData("...")</tt>, e.g.
 * <pre>
 *    AliReconstruction rec;
 *    rec.SetUseHLTData("TPC TRD");
 * </pre>
 * will replace the input of TPC and TRD.
 *
 * The reader can be used directly. In order to avoid library dependencies
 * downwards, the methed AliRawHLTManager::CreateRawReaderHLT is available
 * in the RAW package.
 * <pre>
 * {
 *   AliRawReader* orgReader=AliRawReader::Create("raw.root");
 *   AliRawReader* rawreader=AliRawHLTManager::CreateRawReaderHLT(orgReader, "ITSSDD");
 *   rawreader->Select("ITSSDD");
 *   int count=0;
 *   while (rawreader->NextEvent()) {
 *     cout << "scanning event " << count++ << endl;
 *     UChar_t* pSrc=NULL;
 *     while (rawreader->ReadNextData(pSrc)) {
 *       cout << "  equipment: " << rawreader->GetEquipmentId() << endl;
 *     }
 *   }
 * }
 * </pre>
 *
 * @section sec_alirawreaderhlt_detectorids  Detector selection
 * The constructor gets a detector selection string as parameter and initializes
 * the redirection according to that. Detector Ids are according to AliDAQ.
 * Please note the special strings for for ITS and MUON sub-detectors, ITSSPD,
 * ITSSDD, ITSSSD, and MUONTRK and MUONTRG respectively.
 *
 * @section sec_alirawreaderhlt_module  Module implementation
 * In order to determine the equipment id for the data block, the HLT module
 * must implement an HLTOUT handler of class AliHLTOUTHandlerEquId which is of 
 * type @ref AliHLTModuleAgent::AliHLTOUTHandlerType ::kRawReader.
 * The handler must implement the method
 * <pre>
 *  // AliHLTOUTHandlerEquId::ProcessData(AliHLTOUT*)
 *  virtual int ProcessData(AliHLTOUT* pData);
 * </pre>
 * which returns the equipment id and eventually decodes data to be retrieved
 * by calling AliHLTOUTHandler::GetProcessedData(). If the equipment id of the
 * DDL has been sent as data specification of the block, the AliHLTOUTHandlerEquId
 * can be used directly.
 *
 * Secondly, the AliHLTModuleAgent implementation of the module has to create
 * the handler for the data blocks. Depending on the data type and specification,
 * the the following interface methods return handler description and handler. 
 * <pre>
 *   int AliHLTModuleAgent::GetHandlerDescription(AliHLTComponentDataType dt,
 *                                                AliHLTUInt32_t spec,
 *                                                AliHLTOUTHandlerDesc& desc) const;
 *
 *   AliHLTOUTHandler* AliHLTModuleAgent::GetOutputHandler(AliHLTComponentDataType dt, 
 *                                                         AliHLTUInt32_t spec);
 * </pre>
 * See section @ref tut_alirawreaderhlt for sample implementation.
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliRawReaderHLT : public AliRawReader {
 public:
  /** constructor */
  AliRawReaderHLT(AliRawReader* pParentReader, const char* options=NULL);
  /** destructor */
  virtual ~AliRawReaderHLT();

  // interface methods of AliRawReader
  void     Select(Int_t detectorID, 
		  Int_t minDDLID = -1, Int_t maxDDLID = -1);
//   void     Select(const char *detectorName, 
// 		  Int_t minDDLID = -1, Int_t maxDDLID = -1);

  using AliRawReader::Select;

  void     SelectEquipment(Int_t equipmentType, 
			   Int_t minEquipmentId = -1, 
			   Int_t maxEquipmentId = -1);
  void     SkipInvalid(Bool_t skip = kTRUE);
  //  void     SelectEvents(Int_t type);

  UInt_t   GetType() const;
  UInt_t   GetRunNumber() const;
  const UInt_t* GetEventId() const;
  const UInt_t* GetTriggerPattern() const;
  const UInt_t* GetDetectorPattern() const;
  const UInt_t* GetAttributes() const;
  const UInt_t* GetSubEventAttributes() const;
  UInt_t   GetLDCId() const;
  UInt_t   GetGDCId() const;
  UInt_t   GetTimestamp() const;

  const UInt_t* GetEquipmentAttributes() const;
  Int_t    GetEquipmentElementSize() const;
  Int_t    GetEquipmentHeaderSize() const;

  Int_t    GetEquipmentSize() const;
  Int_t    GetEquipmentType() const;
  Int_t    GetEquipmentId() const;
  Bool_t   ReadHeader();
  Bool_t   ReadNextData(UChar_t*& data);
  Bool_t   ReadNextInt(UInt_t& data);
  Bool_t   ReadNextShort(UShort_t& data);
  Bool_t   ReadNextChar(UChar_t& data);
  Bool_t   ReadNext(UChar_t* data, Int_t size);

  Bool_t   Reset();

  Bool_t   NextEvent();
  Bool_t   RewindEvents();

 protected:

 private:
  /** standard constructor prohibited */
  AliRawReaderHLT();
  /** copy constructor prohibited */
  AliRawReaderHLT(const AliRawReaderHLT&);
  /** assignment operator prohibited */
  AliRawReaderHLT& operator=(const AliRawReaderHLT&);

  /**
   * Scan the options.
   * Set the ids for the specified detectors in the detector
   * list. Currently, no other options are available.
   */
  int ScanOptions(const char* options);

  /**
   * Read the next data block from the HLT stream
   */
  Bool_t ReadNextHLTData();

  /**
   * Check if a ddlid is part of the ones which are selected for
   * input replacement.
   */
  Bool_t IsHLTInput(int ddlid);

  /**
   * Check if redirection is enabled for at least one detector in the
   * selected range.
   * Set the fbHaveHLTData variable
   * @return true if data has to be read from the HLT stream.
   */
  Bool_t EvaluateSelection();

  /**
   * Release the current HLT data.
   * Releases the current buffer of either the active HLTOUT data
   * block handler or the HLTOUT instance. The latter implies a
   * reset of the reader concerning the HLT data blocks.
   * @param bReleaseHLTOUT   release HLTOUT instance if \em true
   *                         only current data buffer if \em false
   * @return neg. error code if failed
   */
  int ReleaseHLTData(bool bReleaseHLTOUT=true);

  /**
   * Backbone of all Read functions.
   * Reads the next data into the internal buffer and switches to next
   * block if enabled.
   *
   * @param data             target to receive pointer
   * @param readHeader       kTRUE: switch to next block if no more data
   */
  Bool_t   ReadNextData(UChar_t*& data, Bool_t readHeader);

  /** the rawreader */
  AliRawReader* fpParentReader; //!transient

  /** options */
  TString fOptions; //!transient

  /** system options = options w/o detector strings */
  TString fSystemOptions; //!transient

  /** current data set, either extracted from the HLT stream or parent raw reader */
  const AliHLTUInt8_t* fpData; // !transient

  /** size of the current data set */
  int fDataSize; // !transient

  /** current stream offset for reading from input stream */
  int fOffset; // !transient

  /** current stream position for block input ReadNextData function */
  int fPosition; // !transient

  /** equipment id of the current data set, >0 indicates data set from HLT stream */
  int fEquipmentId; // !transient

  /** indicates the availibility of data from the HLT stream */
  bool fbHaveHLTData; // !transient

  /** list of detectors for which data will be taken from HLT stream */
  vector<int> fDetectors; // !transient

  /** instance of the HLTOUT handler */
  AliHLTOUT* fpHLTOUT; // !transient

  /** start reading HLTOUT from beginning */
  bool fbReadFirst; //!transient

  /** instance of the data handler providing the current data buffer */
  AliHLTOUTHandler* fpDataHandler; // !transient

  /** base class for AliRoot HLT plugins */
  AliHLTPluginBase* fpPluginBase;                                     //!transient

  ClassDef(AliRawReaderHLT, 0)
};

#define ALIHLTREC_LIBRARY                   "libHLTrec.so"
#define ALIHLTREC_LIBRARY_VERSION           0
#define ALIRAWREADERHLT_CREATE_INSTANCE     "AliRawReaderHLTCreateInstance"

#ifdef __cplusplus
extern "C" {
#endif
  typedef AliRawReader* (*AliRawReaderHLTCreateInstance_t)(AliRawReader* pParentReader, const char* options);

  /**
   * Create an instance of the AliRawReader class
   */
  AliRawReader* AliRawReaderHLTCreateInstance(AliRawReader* pParentReader, const char* options);
#ifdef __cplusplus
}
#endif
#endif
