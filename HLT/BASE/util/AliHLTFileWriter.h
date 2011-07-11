// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTFILEWRITER_H
#define ALIHLTFILEWRITER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTFileWriter.h
    @author Matthias Richter
    @date   
    @brief  An HLT file dump (data sink) component.
*/

#include "AliHLTDataSink.h"
#include <TString.h>

class AliHLTBlockDataCollection;

/**
 * @class AliHLTFileWriter
 * An HLT data sink component which writes data to file(s).
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b FileWriter      <br>
 * Library: \b libAliHLTUtil.so     <br>
 * Input Data Types: ::kAliHLTAllDataTypes <br>
 * Output Data Types: none <br>
 *
 * \b Note: ::kAliHLTAllDataTypes contains both ::kAliHLTAnyDataType and 
 * ::kAliHLTVoidDataType
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -datafile     <i> filename   </i> <br>
 *      file name base
 * \li -directory    <i> directory  </i> <br>
 *      target directory
 * \li -subdir[=pattern] <br>
 *      create sub dir for each event, the format pattern can contain printf
 *      specifiers to print the event no into the dir name, default is
 *      'event%%03lu' (-subdir w/o additional pattern). The format specifyer
 *      %%lu is automatically added if missing in the pattern. Please note the
 *      \b long int type of the event id                                   <br>
 *      \b note: the idfmt string is reset since the subdir contains the id
 * \li -idfmt[=pattern] <br>
 *      format specifier for the event id in the file name,                <br>
 *      default: on, default pattern: '_0x%%08x'
 * \li -specfmt[=pattern] <br>
 *      format specifier for the data specification in the file name       <br>
 *      default: off, default pattern: '_0x%%08x'
 * \li -blocknofmt[=pattern] <br>
 *      format specifier for the block no in the file name                 <br>
 *      default: on, default pattern: '_0x%%02x'
 * \li -skip-datatype <br>
 *      do not consider data type when building the file name.
 * \li -enumerate <br>
 *      don't use the event number but an event counter beginning from 0
 * \li -concatenate-blocks <br>
 *      concatenate all blocks of one event into one file, this skips
 *      the block no, and the block data type in the file name
 * \li -concatenate-events <br>
 *      concatenate all events into one file, this skips the event no,
 *      the block no, and the block data type in the file name. Currently,
 *      this implies the -concatenate-blocks option.
 * \li -publisher-conf <i>filename</i> <br>
 *      write configuration file for FilePublisher component (AliHLTFilePublisher) <br>
 *      one line per file: -datatype id origin -datafile filename           <br>
 *      events separated by -nextevent
 * \li -write-all-events <br>
 *      by default, the file writer ignores all steering events like the
 *      the SOR/EOR events, with this option, all events will be considered
 *      the beginning.
 * \li -write-all-blocks <br>
 *      by default, the file writer ignores all blocks of origin {PRIV}
 *      (::kAliHLTDataOriginPrivate), with this option, all blocks will
 *      be written. For SOR/EOR events, a short string will be added in
 *      the beginning.
 * \li -write-all <br>
 *      combines both -write-all-events and -write-all-blocks
 * \li -burst-buffer <size> <br>
 *      size of burst buffer, blocks are written to buffer until it is filled
 *      and written in one burst (though to different files according to conf)<br>
 *      \b Note: burst write is currently only supported for mode
 *      -concatenate-events AND -concatenate-blocks (both enabled).
 * \li -datatype     <i> id origin      </i>                            <br>
 *     data block selection by AliHLTBlockDataCollection
 * \li -origin  <i> origin  </i>                                        <br>
 *     data block selection by AliHLTBlockDataCollection
 * \li -typeid  <i> id      </i>                                        <br>
 *     data block selection by AliHLTBlockDataCollection
 * \li -dataspec     <i> specification </i>                             <br>
 *     data block selection by AliHLTBlockDataCollection
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * No data published (AliHLTDataSink).
 *
 *
 * By default, file name is built from the basename, the event number, the
 * block number and the data type in the format:
 * <pre>
 * basename_eventno_dt
 * </pre>
 * If the basename was not given, \em 'event' ist used instead. A file
 * extension after the last dot is separated from the basename and appended
 * to the final name.
 *
 * The naming rule can be changed by the -xxfmt options, which can contain
 * printf format specifiers in order to print the corresponding variable. E.g.
 * <pre>
 * -specfmt             append specification
 * -subdir=test         store in sub folders
 * -blocknofmt=_0x%%x   format block no in hex
 * -idfmt=_%%04d        print id in 4-digits decimal number
 * -idfmt=              print no id
 * </pre>
 *
 * The class can be used as a base class for file writers. Additional
 * argument scan can be implemented in @ref ScanArgument which is called
 * for each unknown argument.
 *
 * @ingroup alihlt_util_components
 */
class AliHLTFileWriter : public AliHLTDataSink  {
 public:
  /** standard constructor */
  AliHLTFileWriter();
  /** destructor */
  virtual ~AliHLTFileWriter();

  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  virtual AliHLTComponent* Spawn();

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
   * Init the writer.
   * The DoInit function is not available for child classes. InitWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  virtual int InitWriter();

  /**
   * Close the writer.
   * The DoDeinit function is not available for child classes. CloseWriter is the
   * corresponding function for classes derived from AliHLTFileWriter.
   */
  virtual int CloseWriter();

  /**
   * Data processing method for the component.
   * The function can be overloaded by other file writer components.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& trigData );

  using AliHLTDataSink::DumpEvent;

  /**
   * Scan one argument and adjacent parameters.
   * Can be overloaded by child classes in order to add additional arguments
   * beyond the standard arguments of the file publisher. The method is called
   * whenever a non-standard argument is recognized. Make sure to return 
   * <tt> -EPROTO </tt> if the argument is not recognized be the child.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual int ScanArgument(int argc, const char** argv);

  /**
   * Build file name from eventID data type and the specified directory and basename.
   * @param [in] eventID   the ID of the event
   * @param [in] blockID   the ID of the current block
   *                       no block string appended if -1
   * @param [in] dataType  the data type of the data block
   *                       no type string appanded if @ref kAliHLTVoidDataType
   * @param [in] specification  data specification of the block
   * @param [out] filename string to receive the file name
   */
  int BuildFileName(const AliHLTEventID_t eventID, const int blockID,
		    const AliHLTComponentDataType& dataType,
		    const AliHLTUInt32_t specification,
		    TString& filename);

  /**
   * Set a mode flag.
   * @return current mode flags
   */
  int SetMode(Short_t mode);
    
  /**
   * Clear a mode flag.
   * @return current mode flags
   */
  int ClearMode(Short_t mode);

  /**
   * Check a mode flag.
   * @return 1 if flag is set, 0 if not
   */
  int CheckMode(Short_t mode) const;

  /**
   * Get the currently set file extension.
   */
  TString GetExtension() {return fExtension;}

  /**
   * Set the file extension.
   */
  void SetExtension(const char* extension) {fExtension=extension!=NULL?extension:"";}

  /**
   * Get the target directory
   */
  TString GetDirectory() {return fDirectory;}

  /**
   * Working modes of the writer
   * @internal
   */
  enum TWriterMode {
    /**
     * flag to indicate whether to write each incoming block to separate files
     * or all blocks of one event to one file. set = concatenate (one file).
     */
    kConcatenateBlocks = 0x1,

    /**
     * flag to indicate whether to concatenate incoming blocks of the same type
     * for all events to one file. If also @ref kConcatenateBlocks is set,
     * or all blocks of all events are written to the same file.
     */
    kConcatenateEvents = 0x2,

    /** event enumeration flag */
    kEnumerate = 0x4,

    /** write all events including steering events */
    kWriteAllEvents = 0x8,

    /** write all blocks including private ones */
    kWriteAllBlocks = 0x10,

    /** skip the data type information when creating the file name */
    kSkipDataType = 0x20
  };

  /** argument scan concerning block descriptor selections */
  AliHLTBlockDataCollection* fpBlockDataCollection;                //!transient

 private:
  /** copy constructor prohibited */
  AliHLTFileWriter(const AliHLTFileWriter&);
  /** assignment operator prohibited */
  AliHLTFileWriter& operator=(const AliHLTFileWriter&);

  /**
   * Set defaults for all internal properties
   */
  int SetDefaults();

  /**
   * Schedule block for writing.
   * The block is written immediately unless burst mode is activated.
   * In burst mode, the block is buffered in the burst buffer until it is filled.
   * Content of the burst buffer is then written in one burst.
   *
   * In the first implementation, burst write is only foreseen for the base
   * file writer.
   */
  int ScheduleBlock(int blockno, const AliHLTEventID_t& eventID,
		    const AliHLTComponentBlockData* pDesc);

  /**
   * Flush burst buffer.
   */
  int BurstWrite();

  /**
   * Write data block;
   * Build file name from the block attributes and compare with the
   * lat file name in order to correctly append data or not.
   */
  int WriteBlock(int blockno, const AliHLTEventID_t& eventID,
		 const AliHLTComponentBlockData* pDesc);

  /** the basename of the output file */
  TString    fBaseName;                                            // see above
  /** the extension of the output file */
  TString    fExtension;                                           // see above
  /** target directory */
  TString    fDirectory;                                           // see above
  /** base name of the event sub directories */
  TString    fSubDirFormat;                                        // see above
  /** event id format string (when added to file name) */
  TString    fIdFormat;                                            // see above
  /** specification format string (when added to file name) */
  TString    fSpecFormat;                                          // see above
  /** format string for block no (when added to file name) */
  TString    fBlcknoFormat;                                        // see above
 protected:
  /** enumeration format string */
  TString    fCurrentFileName;                                     // see above
 private:

  /** mode specifier, see @ref TWriterMode */
  Short_t    fMode;                                                // see above

  /** burst buffer for postponed data write */
  AliHLTUInt8_t* fpBurstBuffer;                                    //!transient

  /** size of burst buffer */
  AliHLTUInt32_t fBurstBufferSize;                                 //!transient

  /** block descriptor list for postponed burst write*/
  AliHLTComponentBlockDataList fBurstBlocks;                       //!transient

  /** event ids for the burst blocks */
  vector<AliHLTEventID_t> fBurstBlockEvents;                       //!transient

  /// configuration file of FilePublisher component
  TString fPublisherConfName;                                      // see above
  /// current event for FilePublisher configuration
  int fPublisherConfEvent;                                         // see above
  
  ClassDef(AliHLTFileWriter, 0)
};
#endif
