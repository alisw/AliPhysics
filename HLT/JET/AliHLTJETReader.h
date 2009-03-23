//-*- Mode: C++ -*-

// $Id:  $

#ifndef ALIHLTJETREADER_H
#define ALIHLTJETREADER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETReader.h
    @author Jochen Thaeder
    @date   
    @brief  Reader for jet finder
*/

#include "AliHLTLogging.h"

#include "AliJetReader.h"
#include "AliJetReaderHeader.h"

#include "AliHLTJETReaderHeader.h"

#include "AliHLTMCEvent.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"


/**
 * @class  AliHLTJETReader
 * Reader for jet finder
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETReader : public AliJetReader, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETReader();

  /** destructor */
  virtual ~AliHLTJETReader();

  /*
   * ---------------------------------------------------------------------------------
   *                               Reader functionality
   * ---------------------------------------------------------------------------------
   */

  /** Fill tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArray();

  /** Fill MC tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayMC();

  /** Fill ESD tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayESD();

  /** Fill AOD tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayAOD();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set pointer to input event
   *  @param esd an AliESDEvent
   *  @param aod an AliAODEvent
   *  @param mc an AliHLTMCEvent
   */
  void SetInputEvent(TObject* esd, TObject* aod, TObject* mc);

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to AliHLTJETReaderHeader
   *  @return ptr to AliHLTJETReaderHeader
   */
  AliHLTJETReaderHeader* GetReaderHeader() {return dynamic_cast<AliHLTJETReaderHeader*>(fReaderHeader);}

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETReader (const AliHLTJETReader&);

  /** assignment operator prohibited */
  AliHLTJETReader& operator= (const AliHLTJETReader&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** ESD event */
  AliESDEvent*       fESD;                     //! transient

  /** MC event */
  AliHLTMCEvent*     fMC;                      //! transient

  /** AOD event */
  AliAODEvent*       fAOD;                     //! transient

  ClassDef(AliHLTJETReader, 1)

};
#endif

