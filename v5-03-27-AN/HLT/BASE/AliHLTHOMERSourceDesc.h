//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTHOMERSOURCEDESC_H
#define ALIHLTHOMERSOURCEDESC_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTHOMERSourceDesc.h
    @author Jochen Thaeder
    @date   
    @brief  Container for HOMER Sources
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

/**
 * @defgroup alihlt_homer HOMER handling for AliROOT
 * This section describes the handling of HOMER Sources, Blocks 
 * and the HOMER Reader inside the HLT and AliROOT
 */

#include "TString.h"
#include "TNamed.h"

/**
 * @class AliHLTHOMERSourceDesc
 * This class contains the information of 1 homer source: hostname, port for the HOMER 
 * interface as well as data specifications. It used in order to fill these sources in 
 * TLists. ( It has to inherit from TObject ). Further more it knows if this source was
 * selected for read from a user. Mainly used in the AliEVEHOMERManager, as source 
 * objects for AliEVE.
 * 
 * @ingroup alihlt_homer
 */
class AliHLTHOMERSourceDesc : public TNamed {

public:

  /**  constructor */
  AliHLTHOMERSourceDesc();       

  /** destructor */
  virtual ~AliHLTHOMERSourceDesc();

  /*
   * ---------------------------------------------------------------------------------
   *                        Selection - public
   * ---------------------------------------------------------------------------------
   */

  // -- SELECTION --

  /** Set selection state
   * @param state      state, either kTRUE or kFALSE
   */
  void SetState( Bool_t state ) { fSelected = state; }

  /** Checks if Source is selected to readout 
   * @return           returns state, either kTRUE or kFALSE
   */
  Bool_t IsSelected()           { return fSelected; }

  /** Select this source */
  void Select()                 { fSelected = kTRUE; }
  
  /** Deselect this source */
  void Deselect()               { fSelected = kFALSE; }

  /*
   * ---------------------------------------------------------------------------------
   *                        Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Set Service of this source 
   *  @param hostname  hostname of the source
   *  @param port      port of the source
   *  @param origin    HLT data origin
   *  @param type      HLT data type
   *  @param spec      HLT data specification
   */
  void SetService( TString hostname, Int_t port, TString origin, 
		   TString type, TString spec );

  /*
   * ---------------------------------------------------------------------------------
   *                        Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Get node name of this source 
   * @return           hostname
   */
  TString& GetHostname()       { return fHostname; }

  /** Get node name of this source 
   * @return           port
   */
  Int_t   GetPort()           { return fPort; }

  /** Get name of this source 
   * @return           name
   */
  TString& GetSourceName()     { return fSourceName; }

  /** Get detector of this source 
   * @return           detector
   */
  TString& GetDetector()       { return fDetector; }

  /** Get sub detector of this source 
   * @return           subdetector
   */
  Int_t   GetSubDetector()    { return fSubDetector; }

  /** Get sub sub detector of this source 
   * @return           subsubdetector
   */
  Int_t   GetSubSubDetector() { return fSubSubDetector; }

  /** Get HLT data type of this source
   * @return           HLT data type
   */
  TString& GetDataType()       { return fDataType; }

  /** Get HLT specification of this source
   * @return           HLT specification
   */
  ULong_t GetSpecification()  { return fSpecification; }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTHOMERSourceDesc(const AliHLTHOMERSourceDesc&);
  
  /** assignment operator prohibited */
  AliHLTHOMERSourceDesc& operator=(const AliHLTHOMERSourceDesc&);

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  /** is selected to read out */
  Bool_t fSelected;             // see above

  /** Name of Source */
  TString fSourceName;          // see above

  // -- Service Specifications --
  // ----------------------------

  /** Name of HOMER Node */
  TString fHostname;            // see above

  /** Name of HOMER port */
  Int_t fPort;                  // see above

  // -- Data Specifications --
  // -------------------------

  /** HLT DataType */
  TString fDataType;            // see above

  /** Detector Name, e.g. PHOS 
   *  corresponds to HLT origin
   */
  TString fDetector;            // see above

  /** HLT Specification */
  ULong_t fSpecification;       // see above

  /** SubDetector Name e.g. MODULE */
  Int_t   fSubDetector;         // see above

  /** SubSubDetector Name e.g. PARTITION */
  Int_t   fSubSubDetector;      // see above

  ClassDef( AliHLTHOMERSourceDesc, 0 )
};

#endif
