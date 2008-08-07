//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCONTROLTASK_H
#define ALIHLTCONTROLTASK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTControlTask.h
    @author Matthias Richter
    @date   
    @brief  Special task to produce the control events.
*/

#include "AliHLTTask.h"
#include "AliHLTDataSource.h"

class AliHLTComponentHandler;
class AliHLTConfiguration;

/**
 * @class AliHLTControlTask
 * This task is automatically added to the beginning of each chain and
 * produces the special steering events. The first component in every
 * branch get the special events from the task.
 *
 * This task gets initialized with data type, specification and payload
 * of the control event to be sent. It produces the data block if data
 * type differs from fAliHLTVoidDataType. The guard class can be used to
 * set the parameters.
 * <pre>
 * AliHLTControlEventGuard(task, kAliHLTDataTypeSOR, 0, payload, size);
 * </pre>
 *
 * @ingroup alihlt_system
 */
class AliHLTControlTask : public AliHLTTask {
 public:
  /** constructor */
  AliHLTControlTask();
  /** standard destructor */
  virtual ~AliHLTControlTask();

  // AliHLTTask interface function
  int CreateComponent(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH, AliHLTComponent*& pComponent) const;

  class AliHLTControlEventGuard {
  public:
    AliHLTControlEventGuard(AliHLTControlTask* me, AliHLTComponentDataType dt, AliHLTUInt32_t spec, AliHLTUInt8_t* pData, AliHLTUInt32_t size) :
      fTask(me) {
      if (!fTask) return;
      fTask->fEvent=dt; 
      fTask->fSpecification=spec; 
      fTask->fpData=pData; 
      fTask->fSize=size;
    }
      ~AliHLTControlEventGuard() {
	if (!fTask) return;
	fTask->fEvent=kAliHLTVoidDataType;
	fTask->fSpecification=kAliHLTVoidDataSpec;
	fTask->fpData=NULL;
	fTask->fSize=0;
      }

  private:
      /** standard constructor prohibited */
      AliHLTControlEventGuard();
      /** copy constructor prohibited */
      AliHLTControlEventGuard(const AliHLTControlEventGuard&);
      /** assignment operator prohibited */
      AliHLTControlEventGuard& operator=(const AliHLTControlEventGuard&);

      /** by the guard controlled task */
      AliHLTControlTask* fTask; //! transient
  };

  /**
   * Source component producing the data blocks
   */
  class AliHLTControlEventComponent : public AliHLTDataSource {
  public:
    AliHLTControlEventComponent(const AliHLTControlTask* pParent);
    ~AliHLTControlEventComponent();

    // AliHLTComponent interface functions
    const char* GetComponentID() {return "__priv_AliHLTControlTask";}
    AliHLTComponentDataType GetOutputDataType();
    int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
    void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
    AliHLTComponent* Spawn() {return NULL;}

  private:
    /** standard constructor prohibited */
    AliHLTControlEventComponent();
    /** copy constructor prohibited */
    AliHLTControlEventComponent(const AliHLTControlEventComponent&);
    /** assignment operator prohibited */
    AliHLTControlEventComponent& operator=(const AliHLTControlEventComponent&);

    // AliHLTDataSource interface function
    int GetEvent(const AliHLTComponentEventData& evtData,
		 AliHLTComponentTriggerData& trigData,
		 AliHLTUInt8_t* outputPtr, 
		 AliHLTUInt32_t& size,
		 vector<AliHLTComponentBlockData>& outputBlocks );

    const AliHLTControlTask* fpParent; //! transient
  };

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTControlTask(const AliHLTControlTask&);
  /** assignment operator prohibited */
  AliHLTControlTask& operator=(const AliHLTControlTask&);

  /** data type of the control event */
  AliHLTComponentDataType fEvent; //! transient
  /** specification of the control evtent */
  AliHLTUInt32_t fSpecification; //! transient
  /** payload to be sent with the control event */
  AliHLTUInt8_t* fpData; //! transient
  /** payload size */
  AliHLTUInt32_t fSize; //!transient

  ClassDef(AliHLTControlTask, 0)
    };
#endif
