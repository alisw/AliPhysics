// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTMONITORINGRELAY_H
#define ALIHLTMONITORINGRELAY_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTMonitoringRelay.h
/// @author Matthias Richter
/// @date   2009-11-11
/// @brief  Relay components for monitoring objects.
///

#include "AliHLTProcessor.h"
#include "TString.h"
#include <vector>

using std::vector;

class TArrayC;
class TObject;

/**
 * @class AliHLTMonitoringRelay
 * A relay component for monitoring data objects.
 * It keeps a copy of the last block of every parent and forwards all
 * the blocks together. By that, the output of histograms (rarely to
 * be published but for every event filled.
 *
 * Input data blocks must be uniquely identified by the combination of
 * - data type (id and origin)
 * - block specification
 * - object name
 * - object title
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MonitoringRelay                                         <br>
 * Library: \b libAliHLTUtil.so						<br>
 * Input Data Types: kAliHLTAnyDataType					<br>
 * Output Data Types: according to input blocks		<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -verbose                                                         <br>
 *      print out some more info messages, mainly for the sake of tutorials
 * \li -check-object                                                    <br>
 *      unpack the object from the binary block and use also object name
 *      and title for indexing
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * Low profile: input objects are unpacked and binary copied, no streaming
 * of obejcts.
 *
 * <h2>Memory consnumption:</h2>
 * The component allocates memory of the maximum size for every input
 * object.
 *
 * <h2>Output size:</h2>
 * 
 * @ingroup alihlt_util_components
 */
class AliHLTMonitoringRelay : public AliHLTProcessor
{
 public:
  /// standard constructor
  AliHLTMonitoringRelay();
  /// destructor
  virtual ~AliHLTMonitoringRelay();

  /// inherited from AliHLTComponent, get component id
  virtual const char* GetComponentID() {return "MonitoringRelay";};

  /// inherited from AliHLTComponent, get the input data type
  void GetInputDataTypes( AliHLTComponentDataTypeList& );

  /// inherited from AliHLTComponent, get the output data type
  AliHLTComponentDataType GetOutputDataType();

  /// inherited from AliHLTComponent, get the output data size estimator
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /// inherited from AliHLTComponent, create a component
  virtual AliHLTComponent* Spawn() {return new AliHLTMonitoringRelay;}

  enum {
    kCheckObject = 0x1
  };

  /// descriptor of monitoring items
  class AliHLTMonitoringItem {
  public:
    // standard constructor
    AliHLTMonitoringItem();
    // constructor
    AliHLTMonitoringItem(const AliHLTComponentBlockData* pBlock, const TObject* pObject);
    // destructor
    ~AliHLTMonitoringItem();

    /// copy data from buffer
    int SetData(void* pBuffer, int size);

    /// get buffer pointer of the current data
    void* GetBuffer() const;
    /// get size of the current data
    unsigned GetSize() const;
    /// get data type
    const AliHLTComponentDataType& GetDataType() const;
    /// get specification
    AliHLTUInt32_t GetSpecification() const;
    /// get object name
    const TString& GetObjectName() const {return fName;}

    bool operator==(const AliHLTComponentBlockData& bd) const;
    bool operator!=(const AliHLTComponentBlockData& bd) const;
    bool operator==(const TObject& object) const;
    bool operator!=(const TObject& object) const;

  protected:
  private:
    /// copy constructor prohibited
    AliHLTMonitoringItem(const AliHLTMonitoringItem&);
    /// assignment operator prohibited
    AliHLTMonitoringItem& operator=(const AliHLTMonitoringItem&);

    AliHLTComponentDataType fDt;                //! transient
    AliHLTUInt32_t          fSpecification;     //! transient
    TString                 fName;              //! transient
    TString                 fTitle;             //! transient
    TArrayC*                fData;              //! transient
    int                     fDataSize;          //! transient
  };
  typedef vector<AliHLTMonitoringItem*>  AliHLTMonitoringItemPList;
 protected:

  /// inherited from AliHLTProcessor, data processing
  int DoEvent( const AliHLTComponentEventData& evtData,
	       AliHLTComponentTriggerData& trigData );
  
  using AliHLTProcessor::DoEvent;

  /// inherited from AliHLTComponent, component initialisation
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent, scan argument
  int ScanConfigurationArgument(int argc, const char** argv);

  /// inherited from AliHLTComponent, component cleanup.
  int DoDeinit();

 private:
  /// copy constructor prohibited
  AliHLTMonitoringRelay(const AliHLTMonitoringRelay&);
  /// assignment operator prohibited
  AliHLTMonitoringRelay& operator=(const AliHLTMonitoringRelay&);

  /// find a block of data type and specificaton
  AliHLTMonitoringItem* FindItem(const AliHLTComponentBlockData* pBlock, const TObject* pObject) const;

  void SetFlag(unsigned flag) {fFlags|=flag;}
  bool CheckFlag(unsigned flag) const {return (fFlags&flag)!=0;}

  /// the list of items
  AliHLTMonitoringItemPList  fItems; //! transient
  /// actual size of the data sample
  unsigned fOutputSize; //! transient
  /// operation flags
  unsigned fFlags; //! transient

  ClassDef(AliHLTMonitoringRelay, 0)
};
#endif
