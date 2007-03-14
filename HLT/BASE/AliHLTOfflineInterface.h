// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOFFLINEINTERFACE_H
#define ALIHLTOFFLINEINTERFACE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOfflineInterface.h
    @author Matthias Richter
    @date   
    @brief  the HLT interface to AliRoot
*/

#include <TObject.h>
#include <TList.h>

class AliRunLoader;
class AliRawReader;
class AliESD;
class TTree;

/******************************************************************************/

/**
 * @class AliHLTOfflineInterface
 * The class implements the basic interface to the AliRoot objects during
 * reconstructions.
 * It serves as a base class for offline source and sink interface components
 * and provides access methods for the AliRunLoader, AliRawReader and AliESD
 * objects. The AliRunLoader and the AliRawReader are fixed during one run,
 * while the AliESD object will be changed from event to event.<br>
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
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTOfflineInterface(const AliHLTOfflineInterface&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTOfflineInterface& operator=(const AliHLTOfflineInterface&);
  /** destructor */
  virtual ~AliHLTOfflineInterface();

  /**
   * Get the AliRoot run loader.
   */
  const AliRunLoader* GetRunLoader() const;

  /**
   * Get the AliRoot raw reader
   */
  const AliRawReader* GetRawReader() const;

  /**
   * Set AliRoot ESD for the current event.
   */
  int SetESD(Int_t eventNo, AliESD* pESD);

  /**
   * Get the AliRoot ESD
   */
  AliESD* GetESD() const;

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
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESD instance
   * @return neg. error code if failed 
   */
  virtual int FillESD(AliRunLoader* runLoader, AliESD* esd)=0;

  /**
   * Fill ESD for one event.
   * The FillESD method of all active AliHLTOfflineDataSink's is called in
   * order to fill the ESD with the previously reconstructed data. This method
   * works on the global list.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESD instance
   * @return neg. error code if failed 
   */
  static int FillComponentESDs(AliRunLoader* runLoader, AliESD* esd);

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
  /** the list of active interfaces */
  static TList fgList;                                            // see above

  /** the current object link (list position) */
  static TObjLink* fgCurrentLnk;                                  // see above

  /** AliRoot run loader instance */
  AliRunLoader* fpRunLoader;                                      //! transient
  /** AliRoot raw reader instance */
  AliRawReader* fpRawReader;                                      //! transient
  /** AliRoot HLT ESD instance */
  AliESD* fpESD;                                                  //! transient

  ClassDef(AliHLTOfflineInterface, 1);
};

#endif
