// -*- Mode: C++ -*-
// $Id: AliHLTESDMCEventPublisherComponent.h 27447 2008-07-19 21:59:56Z richterm $

#ifndef ALIHLTESDMCEVENTPUBLISHERCOMPONENT_H
#define ALIHLTESDMCEVENTPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTESDMCEventPublisherComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Component for publishing ESD and MC events.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTFilePublisher.h"

#include "TList.h"
#include "TTree.h"
#include "TString.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"

#include "AliHLTMCEvent.h"

/**
 * @class AliHLTESDMCEventPublisherComponent
 * An HLT data source component which publishes AliESDEvent and AliMCEvent objects
 * out of a series of datapaths.<br>
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ESDMCEventPublisher <br>
 * Library: \b libAliHLTUtil.so     <br>
 * Input Data Types: none <br>
 * Output Data Types: according to arguments <br>
 *  - AliESDEvent    -> kAliHLTDataTypeESDObject
 *     - HLTESD      -> kAliHLTDataOriginHLT
 *     - ESD         -> kAliHLTDataOriginOffline
 *
 *  - AliMCEvent     -> kAliHLTDataTypeMCObject
 *                   -> kAliHLTDataOriginOffline
 *
 *  - AliHLTMCEvent  -> kAliHLTDataTypeMCObject
 *                   -> kAliHLTDataOriginHLT
 * 
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -entrytype    <i> Type of events to publish   </i> <br>
 *      Can be one, all or some of :<br>
 *      - ESD<br>
 *      - HLTESD<br>
 *      - MC (publishes both AliHLTMCEvent and AliMCEvent) <br>
 *      - MCFAST (publishes both AliHLTMCEvent and AliMCEvent created from FastSim) <br>
 *
 * \li -datapath     <i> Path to list of data files     </i><br>
 *      - AliESDs.root<br>
 *      - Kinematics.root<br>
 *      - galice.root<br>
 *      - TrackRefs.root<br>
 *       
 * <h2>Optional arguments:</h2>
 * \li -dataspec     <i> Specification </i> <br>
 *      Data specification treated as decimal number or hex number if
 *      prepended by '0x'<br>
 *      If not given void spec ist used. Otherwise each Bit corresponds to
 *      the detectorID specified at ALICE-INT-2007-016, Table 1
 *
 * \li -applyParticleCuts  <i> Apply particle cuts before filling in AliHLTMCEvent </i> <br>
 *
 * \li -skip-esd-object  <i> object </i> <br>
 *      remove object from ESD content before publishing, can be a blank separated list of
 *      objects, but don't forget to enclose multiple names in quotes
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
 * According to the available data. The component is an AliHLTDataSource
 * and inteded to be used in the AliHLTSystem framework only. The component
 * implements the standard AliHLTSystem adaptive buffer allocation. 
 *
 * The component needs at least one argument \em -datapath. 
 * It can occur multiple times. The \em -entrytype and \em -dataspec
 * parameters are valid for all data paths
 *
 * All files are broken up and published in individual events. Then one data 
 * block is pulished for each entrytype.
 *
 * @ingroup alihlt_util_components
 */

class AliHLTESDMCEventPublisherComponent : public AliHLTFilePublisher  {
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** standard constructor */
  AliHLTESDMCEventPublisherComponent();

  /** destructor */
  virtual ~AliHLTESDMCEventPublisherComponent();
  
  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /**
   * Init method. 
   * Overwrites the AliHLTFilePublisher::DoInit() method.  
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  Int_t DoInit( int argc, const char** argv );


  /**
   * Deinit method. Calls also the one of AliHLTFilePublisher.
   */
  Int_t DoDeinit();

  /**
   * Data processing method for the component.
   * The component uses the @ref alihltcomponent-high-level-interface
   * to put serialized Root object into the output stream. Despite of that it
   * implements the low-level DumpEvent method in order to allow child classes
   * to use the low-level method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr     not used
   * @param size          not used
   * @param outputBlocks  not used
   * @return
   */
  Int_t GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks);

  using AliHLTFilePublisher::GetEvent;

 private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /** copy constructor prohibited */
  AliHLTESDMCEventPublisherComponent(const AliHLTESDMCEventPublisherComponent&);

  /** assignment operator prohibited */
  AliHLTESDMCEventPublisherComponent& operator=(const AliHLTESDMCEventPublisherComponent&);


  /*
   * ---------------------------------------------------------------------------------
   *                    Helper functions - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Add output datatypes according to the fPublish* flags
   *  - AliESDEvent -> kAliHLTDataTypeESDObject
   *     - HLTESD   -> kAliHLTDataOriginHLT
   *     - ESD      -> kAliHLTDataOriginOffline
   *
   *  - AliMCEvent  -> kAliHLTDataTypeMCObject
   *                -> kAliHLTDataOriginOffline
   */
  void AddDataTypesToOutputlist();

  /** Insert datafiles according to the fPublish* flags 
   *  into folders.
   *  @return negative number in error case
   */
  Int_t InsertFiles();

  /** Open all files for current folder. Get ESD tree's and TreeE.
   *  @return negative number in error case
   */
  Int_t OpenCurrentFileList();

  /** Close all files for current folder.
   *  @return negative number in error case
   */
  Int_t CloseCurrentFileList();

  /** clone an ESD by copying all objects but skip the specified ones
   */
  Int_t CopyESDObjects(AliESDEvent* pTgt, const AliESDEvent* pSrc, const char* skippedObjects) const;

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** The current folder, containing 1 set of files */
  TObjLink *fpCurrentFolder;                 //! transient

  /** List of files in current folder*/
  TList *fpCurrentFileList;                  //! transient
  
  /** Event in current folder ( inside files ) */
  UInt_t fCurrentEvent;                      //  see above

  /** Number of event in current folder ( inside files ) */
  UInt_t fNEventsInFolder;                   //  see above

  /** List containing TObjStrings 
   *  -> Contain paths to reconstructed data
   */
  TList fFolderList;                         //! see above

  /** Data specification */
  AliHLTUInt32_t fSpecification;             //  see above

  /** Publish class AliESDEvent, containing normal ESD */
  Bool_t fPublishESD;                        //  see above

  /** Publish class AliESDEvent, containing normal HLTESD */
  Bool_t fPublishHLTESD;                     //  see above

  /** Publish class AliMCEvent */
  Bool_t fPublishMC;                         //  see above

  /** Fill AliMCEvent without TrackRefs, out of Fast Simulation */
  Bool_t fFastMC;                            //  see above

  /** Pointer to ESD tree in current file */
  TTree* fpTreeESD;                          //! transient

  /** Pointer to HLT ESD tree in current file */
  TTree* fpTreeHLTESD;                       //! transient

  /** Pointer to TreeE tree in current galice file */
  TTree* fpTreeE;                            //! transient

  /** Pointer to TreeK tree in current kinematics file 
   *  - changes every event
   */
  TTree* fpTreeK;                            //! transient

  /** Pointer to TreeTR tree in current track refernce file 
   *  - changes every event
   */
  TTree* fpTreeTR;                           //! transient

  /* Ptr to current AliESDEvent, to be shipped out*/
  AliESDEvent* fpESD;                        //! transient

  /* Ptr to current HLT - AliESDEvent, to be shipped out*/
  AliESDEvent* fpHLTESD;                     //! transient

  /* Ptr for ESD with selected objects, to be shipped out*/
  AliESDEvent* fpESDClone;                   //! transient

  /* Ptr to current AliMCEvent, to be shipped out*/
  AliMCEvent* fpMC;                          //! transient
  
  /* Ptr to current AliHLTMCEvent, to be shipped out*/
  AliHLTMCEvent* fpHLTMC;                    //! transient

  /* Maximum required output size */
  UInt_t fOutputSize;                        //! transient

  /** Apply particle cuts, before filling in AliHLTMCEvent */
  Bool_t fApplyParticleCuts;                 // see above

  /// list of ESD objects to be skipped
  TString fSkippedEsdObjects;                //! transient

  ClassDef(AliHLTESDMCEventPublisherComponent, 0)
};
#endif
