//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTESDMANAGER_H
#define ALIHLTESDMANAGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTEsdManager.h
    @author Matthias Richter
    @date   
    @brief  Manager for merging and writing of HLT ESDs
*/

#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "TString.h"

class AliESDEvent;

/**
 * @class AliHLTEsdManager
 * Tool to write and merge HLT ESD objects.
 *
 * HLT components can produce ESD output. The ESD objects are sent via
 * a TMessage like mechanism as part of the HLTOUT data. This class retrieves
 * streamed AliESDEvent objects from an HLT output block. An ESD object can be
 * copied to a global ESD provided by the caller or to files. The name of the
 * ROOT files follows the scheme AliHLTDETESDs.root where DET denotes a detector.
 * E.g. the ESD from a data block of type {ESD_TREE,TPC} will be added to the
 * file AliHLTTPCESDs.root.
 *
 * For the sake of library (in)dependencies, the concrete implementation of
 * the AliHLTEsdManager is separated from the libHLTbase class as this would
 * introduce dependencies to AliRoot libraries. The user does not notice the
 * handling apart from the fact that an instance can not be created and
 * deleted directly. Instead 
 * <pre>
 *   AliHLTEsdManager* pManager=AliHLTEsdManager::New();
 *   // ....
 *   AliHLTEsdManager::Delete(pManager);
 * </pre>
 * must be used.
 *
 * @ingroup alihlt_aliroot_reconstruction
 */
class AliHLTEsdManager : public AliHLTLogging {
 public:
  /** create an instance of the manager */
  static AliHLTEsdManager* New();
  /** delete an instance of the manager */
  static void Delete(AliHLTEsdManager* instance);

  /**
   * Set the options for the ESD merging and writing
   * Takes a string of blank separated options.
   * \li -writelocal   use local file writing in the specified diectory
   *                   The file name is derived from the data origin
   * \li -directory=<> Set the directory path
   *                   This makes the SetDirectory method obsolete
   * @return neg error code if failed
   */
  virtual int SetOption(const char* option)=0;

  /**
   * Convert data buffer to ESD.
   * The buffer is supposed to describe a streamed AliESDEvent object.
   * If no target object is specified, the ESD is written to a file AliHLTdetESDs.root,
   * where 'det' is derived from the data type origin. Each time the function is invoked
   * a new event is created. Dummy events are added if the previous events did not contain
   *
   * @param [in]  pBuffer  the data buffer
   * @param [in]  size     data buffer size
   * @param [in]  dt       data type of the block
   * @param [out] tgtesd   optional target
   * @param [in]  eventno  optional event no
   */
  virtual int WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size, AliHLTComponentDataType dt,
	       AliESDEvent* tgtesd=NULL, int eventno=-1)=0;

  /**
   * Merge content of source ESD into the target ESD.
   * Merging is done on the level of objects in the ESD and for the
   * moment it's only implemented for the TClonesArrays. In that case it's
   * easy to detect whether the object is empty or not.
   *
   * \b Note: The function can not match entries of the same type, like e.g.
   * tracks from the 'Tracks' member.
   */
  virtual int Merge(AliESDEvent* pTgt, AliESDEvent* pSrc) const =0;

  /**
   * Align all ESD to the same number of events.
   * The function adds empty events to all ESD files if their event number
   * does not match the specified one.
   * @param eventno     the desired event no
   * @return neg. error code if failed
   */
  virtual int PadESDs(int eventno)=0;

  /**
   * Set the target directory for the ESD files.
   */
  virtual void SetDirectory(const char* directory)=0;

  /**
   * Get the list of the internally created files.
   * Returns a blank separated list of the file names.
   */
  virtual TString GetFileNames(AliHLTComponentDataType dt=kAliHLTAnyDataType) const = 0;

  /**
   * Create an AliESDEvent object.
   * The standard content can optionally be initialized.
   */
  virtual TObject* CreateEsdEvent(bool bCreateStdContent=false) const = 0;

  /**
   * Delete instance of AliESDEvent
   */
  virtual int DestroyEsdEvent(TObject* pESDInstance) const=0;

  /**
   * Add object to ESD event.
   * Checks the existance of the object under the name 'branchname'
   * Note: some of the objects have (branch-)names which differ from the object name
   * However, parameter branchname is never used when adding an object not yet existing.
   */
  virtual int AddObject(TObject* pESDInstance, const TObject* pObject, const char* branchname) const = 0;

  /**
   * Reset the specified object.
   * The purpose of this method is to disentangle library dependencies.
   * The actual implementation is outside the HLT base library in a
   * child class.
   */
  virtual int ResetEsdEvent(TObject* pESDInstance) const = 0;

 protected:
  /** constructor */
  AliHLTEsdManager();
  /** destructor */
  virtual ~AliHLTEsdManager();

 private:
  /** copy constructor prohibited */
  AliHLTEsdManager(const AliHLTEsdManager&);
  /** assignment operator prohibited */
  AliHLTEsdManager& operator=(const AliHLTEsdManager&);

  /** the name of the actual implementation */
  static const char* fgkImplName; //!

  /** the library of the implementation */
  static const char* fgkImplLibrary; //!

  ClassDef(AliHLTEsdManager, 0)
};

#endif
