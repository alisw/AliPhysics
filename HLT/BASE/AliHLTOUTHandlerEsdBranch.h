//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHANDLERESDBRANCH_H
#define ALIHLTOUTHANDLERESDBRANCH_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOUTHandlerEsdBranch.h
/// @author Matthias Richter
/// @date   01.07.2010
/// @brief  HLTOUT handler of type kEsd to merge objects into the hltEsd.

#include "AliHLTOUTHandler.h"
#include "TString.h"

class TArrayC;
class AliESDEvent;

/**
 * @class AliHLTOUTHandlerEsdBranch
 * An HLTOUT handler of type kEsd to add objects to hltEsd branches.
 *
 * The handler extracts objects from HLTOUT data blocks or converts
 * data to objects to be added to hltEsd branches. The default implementation
 * covers the first case right away, the class can be used directly for single
 * objects streamed to the HLTOUT.
 *
 * <h2>Object conversion</h2>
 * The method ExtractAndAddObjects() has to loop over all input blocks and
 * provide an appropriate conversion. If the data block simply contains a
 * streamed object it just needs to be extracted and added to the ESD using
 * the function Add(). Thhis case is covered by the default implementation.
 * Child classes can overload ExtractAndAddObjects() if there is further
 * conversion/formatting required.
 *
 * <h2>Usage example:</h2>
 * An agent implementation must announce to ability to process a certain
 * data block by implementing AliHLTModuleAgent::GetHandlerDescription()
 * and AliHLTModuleAgent::GetOutputHandler(). See AliHLTModuleAgent for
 * more details.
 * <pre>
 *  int AliHLTMyAgent::GetHandlerDescription(AliHLTComponentDataType dt,
 *  					     AliHLTUInt32_t spec,
 *  					     AliHLTOUTHandlerDesc& desc) const
 *  {
 *    // add TObject data blocks of type {ROOTTOBJ:SMPL} to ESD
 *    if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginSample)) {
 *        desc=AliHLTOUTHandlerDesc(kEsd, dt, GetModuleId());
 *        return 1;
 *    }
 *  
 *    return 0;
 *  }
 *
 *  AliHLTOUTHandler* AliHLTMyAgent::GetOutputHandler(AliHLTComponentDataType dt,
 *                                                    AliHLTUInt32_t spec)
 *  {
 *   // merge my objects into the hltEsd
 *   if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginSample)) {
 *     static AliHLTOUTHandlerEsdBranch handler;
 *     return &handler;
 *   }
 *
 *   return NULL;
 *  }
 * </pre>
 *
 * <h2>Data output</h2>
 * The handler produces a partial ESD containing the data objects. The framework
 * merges all the different partial ESDs in the AliHLTEsdManager, respectively the
 * specific implementation AliHLTEsdManagerImplementation.
 *
 * HLTOUT processing sequence:
 * - first handlers of type kChain
 * - handlers of type kEsd
 * - handlers of type kProprietary
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTOUTHandlerEsdBranch : public AliHLTOUTHandler {
 public:
  /** constructor */
  AliHLTOUTHandlerEsdBranch(const char* branchname=NULL);
  /** standard destructor */
  virtual ~AliHLTOUTHandlerEsdBranch();

  /**
   * Process a data block.
   * @return 
   */
  int ProcessData(AliHLTOUT* pData);

  int GetProcessedData(const AliHLTUInt8_t* &pData);
  int ReleaseProcessedData(const AliHLTUInt8_t* pData, int size);

 protected:
  /**
   * Extract and add objects
   * Loop over input blocks and extract/format the objects. Child class
   * can implement specific conversion. The default implementation just
   * extracts and adds objects.
   */
  virtual int ExtractAndAddObjects(AliHLTOUT* pData);

  /**
   * Add object to internal ESD.
   * If there is an object with 'branchname' the object is copied, otherwise
   * added.
   */
  virtual int Add(TObject* object, const char* branchname);

 private:
  /** copy constructor prohibited */
  AliHLTOUTHandlerEsdBranch(const AliHLTOUTHandlerEsdBranch&);
  /** assignment operator prohibited */
  AliHLTOUTHandlerEsdBranch& operator=(const AliHLTOUTHandlerEsdBranch&);

  TString fBranch; //! transient
  AliESDEvent* fESD; //! transient
  TArrayC* fpData;  //! transient
  int fSize; //! transient

  ClassDef(AliHLTOUTHandlerEsdBranch, 1)
};
#endif
