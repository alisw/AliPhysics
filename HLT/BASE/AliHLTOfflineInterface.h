// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOFFLINEINTERFACE_H
#define ALIHLTOFFLINEINTERFACE_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOfflineInterface.h
/// @author Matthias Richter
/// @date   
/// @brief  the HLT interface to AliRoot
///

#include <TObject.h>
#include <TList.h>

class AliRunLoader;
class AliRawReader;
class AliESDEvent;
class TTree;

/******************************************************************************/

/**
 * @class AliHLTOfflineInterface
 * The class implements the basic interface to the AliRoot objects during
 * reconstructions.
 * It serves as a base class for offline source and sink interface components
 * and provides access methods for the AliRunLoader, AliRawReader and AliESDEvent
 * objects. The AliRunLoader and the AliRawReader are fixed during one run,
 * while the AliESDEvent object will be changed from event to event.<br>
 * \em Note: The digits and clusters trees are not available through this
 * interface class as they are completetly detector (AliLoader) dependend.
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTOfflineInterface : public TObject {
 public:
  /** standard constructor */
  AliHLTOfflineInterface();
  /** constructor
   *  @param pRunLoader   pointer to AliRoot run loader
   *  @param pRawReader   pointer to AliRoot raw reader
   */
  AliHLTOfflineInterface(AliRunLoader* pRunLoader, AliRawReader* pRawReader);
  /** destructor */
  virtual ~AliHLTOfflineInterface();

  /**
   * Get the AliRoot run loader.
   */
  AliRunLoader* GetRunLoader() const;

  /**
   * Get the AliRoot raw reader
   */
  AliRawReader* GetRawReader() const;

  /**
   * Set AliRoot ESD for the current event.
   */
  int SetESD(Int_t eventNo, AliESDEvent* pESD);

  /**
   * Get the AliRoot ESD
   */
  AliESDEvent* GetESD() const;

  /**
   * Set AliRoot external params.
   *
   * @param runLoader     the AliRoot runloader
   * @param rawReader     the AliRoot RawReader
   * @return neg. error code if failed 
   */
  int SetParams(AliRunLoader* runLoader, AliRawReader* rawReader);

  /**
   * Reset AliRoot internal params.
   */
  int Reset();

  /**
   * Set AliRoot external params.
   * This method works on the global list.
   * @param runLoader     the AliRoot runloader
   * @param rawReader     the AliRoot RawReader
   * @return neg. error code if failed 
   */
  static int SetParamsToComponents(AliRunLoader* runLoader, AliRawReader* rawReader);

  /**
   * Fill ESD for one event.
   * Fill the ESD with the previously reconstructed data. It must be implmented
   * by the child class.
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  virtual int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd)=0;

  /**
   * Fill ESD for one event.
   * The FillESD method of all active AliHLTOfflineDataSink's is called in
   * order to fill the ESD with the previously reconstructed data. This method
   * works on the global list.
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  static int FillComponentESDs(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

  /**
   * Reset AliRoot internal params of all active components.
   * This method works on the global list.
   */
  static int ResetComponents();

protected:
  /**
   * Register an OfflineInterface.
   * @param me            instance of AliHLTOfflineInterface
   * @return neg. error code if failed
   */
  static int Register(AliHLTOfflineInterface* me);

  /**
   * Unregister  an OfflineInterface.
   * @param me            instance of AliHLTOfflineInterface
   * @return neg. error code if failed
   */
  static int Unregister(AliHLTOfflineInterface* me);

 private:
  /** copy constructor prohibited */
  AliHLTOfflineInterface(const AliHLTOfflineInterface&);
  /** assignment operator prohibited */
  AliHLTOfflineInterface& operator=(const AliHLTOfflineInterface&);

  /** global AliRoot run loader instance (for all components) */
  static AliRunLoader* fgpRunLoader;                              //! transient
  /** global AliRoot raw reader instance (for all components) */
  static AliRawReader* fgpRawReader;                              //! transient

  /** private AliRoot run loader instance */
  AliRunLoader* fpRunLoader;                                      //! transient
  /** private AliRoot raw reader instance */
  AliRawReader* fpRawReader;                                      //! transient
  /** AliRoot HLT ESD instance */
  AliESDEvent* fpESD;                                                  //! transient

  /** the list of active interfaces */
  static AliHLTOfflineInterface* fgAnchor;                        //! transient

  /** next element in the list */
  AliHLTOfflineInterface* fpNext;                                 //! transient

  /** the current element */
  static AliHLTOfflineInterface* fgCurrent;                       //! transient

  /** number of interfaces */
  static int fgCount;                                             //! see above

  ClassDef(AliHLTOfflineInterface, 0);
};

#endif
