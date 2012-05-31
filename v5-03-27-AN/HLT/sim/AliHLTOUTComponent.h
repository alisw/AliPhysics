//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTCOMPONENT_H
#define ALIHLTOUTCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOUTComponent.h
/// @author Matthias Richter
/// @date   
/// @brief  The HLTOUT data sink component similar to HLTOUT nodes.
/// @note   Used in the AliRoot environment only.

#include "AliHLTOfflineDataSink.h"

class AliHLTHOMERLibManager;
class AliHLTMonitoringWriter;
class TFile;
class TTree;
typedef vector<AliHLTMonitoringWriter*> AliHLTMonitoringWriterPVector;

/**
 * @class AliHLTOUTComponent
 * The HLTOUT data sink component which models the behavior of the HLTOUT
 * nodes of the HLT cluster.
 * <h2>General properties:</h2>
 * The HLTOUT component is attached at the end of a chain. It stores all input
 * block in the HOMER format, distributed over a number of DDL link. The data
 * is stored in a digit file or in raw ddl files.
 *
 * Component ID: \b HLTOUT <br>
 * Library: \b libHLTrec.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 * Output Data Types: none (offline data sink) <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -links      <i> n   </i> <br>
 *      number of output ddl links
 * \li -digitfile  <i> name   </i> <br>
 *      name of the digit file to write (default HLT.Digits.root)
 * \li -rawout[=on,off]  <br>
 *      switch raw output on/off (default on)
 * \li -digitout[=on,off]  <br>
 *      switch digit output on/off (default on)
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
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
 * @ingroup alihlt_aliroot_simulation
 */
class AliHLTOUTComponent : public AliHLTOfflineDataSink  {
 public:
  /// type of the HLTOUT component
  enum EType {
    kGlobal = 0, // generate according to global flags
    kDigits = 1, // generate only digits: ID HLTOUTdigits
    kRaw    = 2  // generate only raw:    ID HLTOUTraw
  };
  /// constructor for different component types
  AliHLTOUTComponent(EType type=kGlobal);
  /** destructor */
  virtual ~AliHLTOUTComponent();

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponent* Spawn();

  /**
   * Enable global options valid for all instances of the component
   * @param options   bit field
   */
  static void SetGlobalOption(unsigned int options);

  /**
   * Disable global options valid for all instances of the component
   * @param options   bit field
   */
  static void ClearGlobalOption(unsigned int options);

  /**
   * Test one of the global options
   */
  static bool TestGlobalOption(unsigned int option);

  enum {
    /** write the raw files of the HLT links */
    kWriteRawFiles = 0x1,
    /** write the digit file */
    kWriteDigits = 0x2
  };

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /// inherited from AliHLTComponent,  component specific argument scan
  int ScanConfigurationArgument(int argc, const char** argv);

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
   * Empty now, data written in Write() at the end of DumpEvent()
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @return neg. error code if failed 
   */
  int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

  /**
   * Write the ecoded HLTOUT data to raw and digits files.
   * Originally data was written in the FillESD function of the
   * AliHLTOfflineInterface. Mainly for the sake of availability of the
   * AliLoader. This concept has not turned out to be succesful and the
   * development went a slightly different direction with the concept of
   * HLTOUT handlers.
   * 2010-04-14 change the original FillESD() to Write(), keep the body
   * of the function
   *
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @return neg. error code if failed 
   */
  int Write(int eventNo, AliRunLoader* runLoader);

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
   * Write data for a DDL link.
   * @param hltddl     Number of DDL link within the range of HLT
   * @param pBuffer    buffer to write
   * @param bufferSize size of the buffer
   */
  int WriteDigitArray(int hltddl, const AliHLTUInt8_t* pBuffer, unsigned int bufferSize);

  /**
   * Write the digits for one DDL
   * @param eventNo    number of the event
   * @param runLoader  AliRoot run loader instance
   * @return neg. error if failed
   */
  int WriteDigits(int eventNo, AliRunLoader* runLoader);

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

  /** output buffer, allocated once in order to avoid frequent alloc/dealloc */
  vector<AliHLTUInt8_t> fBuffer; //!transient

  /** instance of the HOMER library manager */
  AliHLTHOMERLibManager* fpLibManager; // !transient

  /** global options for all instances */
  static int fgOptions; //! transient

  /// component options set from component type or global options at DoInit
  int fOptions; //! transient

  /** digit file name */
  TString fDigitFileName; //! transient

  /** the root file for the HLT 'digit' output */
  TFile* fpDigitFile; //!transient

  /** the tree for the HLT 'digit' output */
  TTree* fpDigitTree; //!transient

  /** array of TArrayC output buffers and branches */
  TArrayC** fppDigitArrays; //!transient

  /** Id of HOMER writer kept from previous event */
  int fReservedWriter; //!transient

  /** Data size kept in the internal buffer */
  int fReservedData; //!transient

  /// type of the component
  EType fType; //! type of the component

  /// counter for round robin usage of HLTOUT links
  int fRoundRobinCounter; //! counter for round robin usage of HLTOUT links

  ClassDef(AliHLTOUTComponent, 0)
};
#endif
