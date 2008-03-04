// @(#) $Id$

#ifndef ALIHLTOUTCOMPONENT_H
#define ALIHLTOUTCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTComponent.h
    @author Matthias Richter
    @date   
    @brief  The HLTOUT data sink component similar to HLTOUT nodes
*/

// see class description below
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTOfflineDataSink.h"

class AliHLTHOMERLibManager;
class AliHLTMonitoringWriter;
typedef vector<AliHLTMonitoringWriter*> AliHLTMonitoringWriterPVector;

/**
 * @class AliHLTOUTComponent
 * The HLTOUT data sink component which models the behavior of the HLTOUT
 * nodes of the HLT cluster.
 * <h2>General properties:</h2>
 *
 * Component ID: \b HLTOUT <br>
 * Library: \b libHLTrec.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 * Output Data Types: none (offline data sink) <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -links      <i> n   </i> <br>
 *      number of output ddl links
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * none
 *
 * <h2>Default CDB entries:</h2>
 * none
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component is an offline sink component and has no output data.
 *
 * The component can be used to write data in the same format as
 * the HLTOUT on the real HLT. In case of AliRoot simulation, the
 * component is automatically added to the chain if the specified
 * chains have output data. By that means, the HLT output is added
 * to the simulation.
 *
 * @ingroup alihlt_rec
 */
class AliHLTOUTComponent : public AliHLTOfflineDataSink  {
 public:
  /** standard constructor */
  AliHLTOUTComponent();
  /** destructor */
  virtual ~AliHLTOUTComponent();

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponent* Spawn();

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Deinit method.
   */
  int DoDeinit();

  /**
   * Data processing method for the component.
   * The function can be overloaded by other file writer components.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );

  using AliHLTDataSink::DumpEvent;

  /**
   * Fill ESD for one event.
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

 private:
  /** copy constructor prohibited */
  AliHLTOUTComponent(const AliHLTOUTComponent&);
  /** assignment operator prohibited */
  AliHLTOUTComponent& operator=(const AliHLTOUTComponent&);

  int ShuffleWriters(AliHLTMonitoringWriterPVector &list, AliHLTUInt32_t size);

  /**
   * Fill the output buffer and allocate if neccessary.
   * Assemble ouput buffer with Common Data Header, HLT header and data from the
   * writer. Works on the same buffer witch is allocated once and eventually
   * grown in order to avoid frequent allocs/deallocs.   
   * @param eventNo    number of the event
   * @param pWriter    [IN]  the HOMER writer
   * @param pBuffer    [OUT] target to receive the pointer to buffer
   * @return size of the buffer
   */
  int FillOutputBuffer(int eventNo, AliHLTMonitoringWriter* pWriter, const AliHLTUInt8_t* &pBuffer);

  /**
   * Write the digits for one DDL
   * @param eventNo    number of the event
   * @param runLoader  AliRoot run loader instance
   * @param hltddl     Number of DDL link within the range of HLT
   * @param pBuffer    buffer to write
   * @param size       size of the buffer
   * @return neg. error if failed
   */
  int WriteDigits(int eventNo, AliRunLoader* runLoader, int hltddl, const AliHLTUInt8_t* pBuffer, unsigned int size);

  /**
   * Write the raw file for one DDL
   * @param eventNo    number of the event
   * @param runLoader  AliRoot run loader instance
   * @param hltddl     Number of DDL link within the range of HLT
   * @param pBuffer    buffer to write
   * @param size       size of the buffer
   * @return neg. error if failed
   */
  int WriteRawFile(int eventNo, AliRunLoader* runLoader, int hltddl, const AliHLTUInt8_t* pBuffer, unsigned int size);

  /** list of HOMER writers */
  AliHLTMonitoringWriterPVector fWriters; //!transient

  /** number of DDLs used*/
  int fNofDDLs; //!transient

  /** equipment ID of first HLT DDL */
  int fIdFirstDDL; //!transient

  /** write digits or not */
  Bool_t fWriteDigits; //!transient

  /** write raw file or not */
  Bool_t fWriteRaw; //!transient

  /** output buffer, allocated once in order to avoid frequent alloc/dealloc */
  vector<AliHLTUInt8_t> fBuffer; //!transient

  /** instance of the HOMER library manager */
  AliHLTHOMERLibManager* fpLibManager; // !transient

  ClassDef(AliHLTOUTComponent, 0)
};
#endif
