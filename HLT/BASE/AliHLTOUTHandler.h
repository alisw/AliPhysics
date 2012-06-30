//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHANDLER_H
#define ALIHLTOUTHANDLER_H
///* This file is property of and copyright by the                          * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///* See cxx source for full Copyright notice                               *

/// @file   AliHLTOUTHandler.h
/// @author Matthias Richter
/// @date   
/// @brief  Base class declaration of HLTOUT handlers

#include "AliHLTLogging.h"

class AliHLTOUT;

/**
 * @class AliHLTOUTHandler
 * Base class declaration of HLT output handlers.
 * The library implementation of the AliHLTModuleAgent allows to generate
 * handlers for data blocks of the HLT output. This can be the output of
 * the real HLT coming from the HLTOUT nodes, or simulated HLT output.   <br>
 * \em Note: The created instance of AliHLTOUTHandler is
 * deleted by the framework.
 */
class AliHLTOUTHandler : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTOUTHandler();
  /** standard destructor */
  virtual ~AliHLTOUTHandler();

  /**
   * Process the data.
   * The data blocks can be selected by AliHLTOUT::SelectFirstDataBlock() and
   * AliHLTOUT::SelectNextDataBlock()
   *
   * Properties of the current data block can be retrieved by the following member
   * functions of AliHLTOUT:
   * - AliHLTOUT::GetDataBlockDescription(AliHLTComponentDataType& dt, AliHLTUInt32_t& spec)
   * - AliHLTOUT::GetDataBlockIndex()
   * - AliHLTOUT::GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size)
   * - AliHLTOUT::ReleaseDataBuffer(const AliHLTUInt8_t* pBuffer)
   *
   * The handler might decode the data block and produce new data as a
   * replacement, see GetProcessedData()
   * @param pData     instance of the AliHLTOUT data
   * @return depending on the overloaded function, neg. error code if failed
   */
  virtual int ProcessData(AliHLTOUT* pData) = 0;

  /**
   * Get the output data, if available.
   * Some of the handlers might produce data to replace the original data
   * block. The handler must ensure the internal storage of the buffer and
   * is also responsible for cleaning the buffer. The buffer must be valid
   * until the next call of ProcessData() or ReleaseProcessedData().
   *
   * The default implementation just returns a NULL pointer to indicate
   * 'no data'.
   * @param pData     target to receive data pointer
   * @return size of the buffer
   */
  virtual int GetProcessedData(const AliHLTUInt8_t* &pData);

  /**
   * Release the buffer of processed data.
   * The handler implementation can do cleanup here.
   * @param pData     pointer to buffer
   * @param size      size of the buffer
   * @return neg. error code if failed
   */
  virtual int ReleaseProcessedData(const AliHLTUInt8_t* pData, int size);

  /**
   * Cleanup the current event processing.
   */
  virtual int FinishEvent();

  enum {
    kHandlerUndefined = 0,
    kHandlerOK    = 0,
    kHandlerError = 0x1000
  };

  /**
   * Check state flag of the handler.
   * @return true if flag matches
   */
  bool CheckStatus(unsigned int flag) const {
    return (fState&flag)!=0;
  }

  /**
   * Reset the state flag.
   */
  void ResetState() {
    fState=kHandlerOK;
  }

 protected:
  void SetStatusFlag(unsigned int flag) {
    fState|=flag;
  }

  void ClearStatusFlag(unsigned int flag) {
    fState&=~flag;
  }

 private:
  /** copy constructor prohibited */
  AliHLTOUTHandler(const AliHLTOUTHandler&);
  /** assignment operator prohibited */
  AliHLTOUTHandler& operator=(const AliHLTOUTHandler&);

  /** internal state of the handler */
  int fState; //!transient

  ClassDef(AliHLTOUTHandler, 0)
};
#endif
