// @(#) $Id$

#ifndef ALIHLTFILEPUBLISHER_H
#define ALIHLTFILEPUBLISHER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTFilePublisher.h
    @author Matthias Richter
    @date   
    @brief  An HLT file publishing (data source) component.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTDataSource.h"
#include <TList.h>

/**
 * @class AliHLTFilePublisher
 * An HLT data source component which publishes data from one or a sequence
 * of files.<br>
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -datafile     <i> filename      </i>
 * \li -datafilelist <i> file pattern  </i> <br>
 *      not yet implemented
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype CLUSTERS TPC </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 *
 * Optional arguments:<br>
 *
 * The component needs at least one argument \em -datafile or \em -datafilelist.
 * Both can occur multiple times. The \em -datatype and \em -dataspec
 * parameters are valid for all files until the next occurrence of
 * \em -datatype/spec
 * @ingroup alihlt_component
 */
class AliHLTFilePublisher : public AliHLTDataSource  {
 public:
  /** standard constructor */
  AliHLTFilePublisher();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTFilePublisher(const AliHLTFilePublisher&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTFilePublisher& operator=(const AliHLTFilePublisher&);
  /** destructor */
  virtual ~AliHLTFilePublisher();

  const char* GetComponentID();
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

  /**
   * Open all files.
   * Opens all files from the file name list @ref fFileNames and adds TFile
   * opjects to the TFiles list.
   */
  int OpenFiles();

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
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   * @return
   */
  int GetEvent( const AliHLTComponentEventData& evtData,
		        AliHLTComponentTriggerData& trigData,
		        AliHLTUInt8_t* outputPtr, 
		        AliHLTUInt32_t& size,
		        vector<AliHLTComponentBlockData>& outputBlocks );

  /**
   * Scan one argument and adjacent parameters.
   * Can be overloaded by child classes in order to add additional arguments
   * beyond the standard arguments of the file publisher. The method is called
   * whenever a non-standard argument is recognized.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual int ScanArgument(int argc, const char** argv);

 private:
  TList                   fFileNames;
  TList                   fFiles;
  TObjLink*               fpCurrent; //! transient value
  AliHLTComponentDataType fDataType;
  AliHLTUInt32_t          fSpecification;
  Int_t                   fMaxSize;

  ClassDef(AliHLTFilePublisher, 0)
};
#endif
