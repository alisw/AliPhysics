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
#include <vector>

class AliESDEvent;
class TTree;
class TFile;

/**
 * @class AliHLTEsdManager
 */
class AliHLTEsdManager : public AliHLTLogging {
 public:
  /** constructor */
  AliHLTEsdManager();
  /** destructor */
  virtual ~AliHLTEsdManager();

  /**
   * Convert data buffer to ESD.
   * The buffer is supposed to describe a streamed AliESDEvent object.
   * If no target object is specified, the ESD is written to a file AliHLTdetESDs.root,
   * where 'det' is derived from the data type origin. Each time the function is invoked
   * a new event is created. Dummy events are added if the previous events did not contain
   *
   * @param pBuffer  [in] the data buffer
   * @param size     [in] data buffer size
   * @param dt       [in] data type of the block
   * @param tgtesd   [out] optional target
   * @param eventno  [in] optional event no
   */
  int WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size, AliHLTComponentDataType dt,
	       AliESDEvent* tgtesd=NULL, int eventno=-1);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTEsdManager(const AliHLTEsdManager&);
  /** assignment operator prohibited */
  AliHLTEsdManager& operator=(const AliHLTEsdManager&);

  class AliHLTEsdListEntry : public AliHLTLogging {
  public:
    /** constructor */
    AliHLTEsdListEntry(AliHLTComponentDataType dt);
    /** copy constructor */
    AliHLTEsdListEntry(const AliHLTEsdListEntry& src);
    /** assignment operator */
    AliHLTEsdListEntry& operator=(const AliHLTEsdListEntry& src);
    /** destructor */
    ~AliHLTEsdListEntry();

    /**
     * Write the ESD to the corresponding file.
     * The tree is first synchronized with the eventno and additional empty
     * events might be inserted if there was a gap. Since we are writing
     * several files in parallel, we have to make sure that those files contain
     * the same number of events.
     * @param pESD        ESD to write
     * @param eventno     optional event no for tree synchronization
     */
    int WriteESD(AliESDEvent* pESD, int eventno=-1);

    bool operator==(AliHLTComponentDataType dt) const;

  private:

    /** root file name */
    TString fName; //!transient
    /** the root file for this esd */
    TFile* fpFile; //!transient
    /** the tree for this esd */
    TTree* fpTree; //!transient
    /** the esd to fill into the tree */
    AliESDEvent* fpEsd; //!transient
    /** data type of the corresponding block */
    AliHLTComponentDataType fDt; //!transient
  };

  typedef vector<AliHLTEsdListEntry> AliHLTEsdList;

  /**
   * Find list entry for given data type
   */
  AliHLTEsdListEntry* Find(AliHLTComponentDataType dt) const;

  /** the list of the ESDs */
  AliHLTEsdList fESDs; //!transient

  ClassDef(AliHLTEsdManager, 0)
};
#endif
