//-*- Mode: C++ -*-
// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVEHOMERSRCTRANSLATOR_H
#define ALIEVEHOMERSRCTRANSLATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliEveHOMERSrcTranslator.h
    @author Jochen Thaeder
    @date
    @brief  Src Translator of HomerManger
*/

#include "TString.h"
#include "TDOMParser.h"
#include "TXMLNode.h"
#include "TList.h"

#include "AliEveHOMERSourceList.h"
#include "AliHLTHOMERSourceDesc.h"

#include "TMap.h"

class AliEveHOMERSrcTranslator : public TObject
{
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliEveHOMERSrcTranslator( TString realm );

  /** destructor */
  virtual ~AliEveHOMERSrcTranslator();

  /*
   * ---------------------------------------------------------------------------------
   *                            Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Sets realm ( which can be ACR, GPN, HLT, KIP ) */
  void SetRealm( TString s ) { fRealm = s; } // Sets realm ( which can be ACR, GPN, HLT, KIP )

  /*
   * ---------------------------------------------------------------------------------
   *                            Translation - public
   * ---------------------------------------------------------------------------------
   */
  
  /** Resolve Information of nodename and port for source which has to be used by HOMER */
  Int_t Translate( TString Nodename, TString xmlPort, TString &hostname, Int_t &port );

  /** Apply corrections for differnt detectors and subdetectors */
  void ApplyDetectorCorrections( TString &detector, TString &subDetector);

  /** Fill SourceDesc with object Information */
  Int_t FillSourceDesc( AliHLTHOMERSourceDesc* source, TString name );

  ///////////////////////////////////////////////////////////////////////////////////

private:

  AliEveHOMERSrcTranslator(const AliEveHOMERSrcTranslator&);            // Not implemented.
  AliEveHOMERSrcTranslator& operator=(const AliEveHOMERSrcTranslator&); // Not implemented.

  /*
   * ---------------------------------------------------------------------------------
   *                            Source Resolving - private
   * ---------------------------------------------------------------------------------
   */

  /** Resolve hostname out of XML and realm */
  TString ResolveHostname( TString nodename );

  /** Resolve port out of XML and realm */
  Int_t ResolvePort( TString srcPort, TString srcHostname );

  /*
   * ---------------------------------------------------------------------------------
   *                             Setup - private
   * ---------------------------------------------------------------------------------
   */

  /** Setup the BasePort mapping */
  void SetupPortMap();

  /** Setup the Object mapping */
  void SetupObjectMap();

  /** Setup the Object mapping for TPC */
  void SetupObjectMapTPC();

  /** Setup the Object mapping for TRD */
  void SetupObjectMapTRD();

  /** Setup the Object mapping for PHOS */
  void SetupObjectMapPHOS();

  /** Setup the Object mapping for DIMUON */
  void SetupObjectMapDIMUON();

  /** Setup the common Object mappings */
  void SetupObjectMapCommon( TMap* objectMap );

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  TMap* fBasePortMap;                  //! Map of BasePorts on the gateways   

  TMap* fObjectMap;                    //! Map of Objects and DataTypes

  TString fRealm;                      // Indicates the realm where AliEve can connect to ( HLT, GPN, ACR, KIP );

  /*
   * ---------------------------------------------------------------------------------
   *                          Constants - private
   * ---------------------------------------------------------------------------------
   */

  static const Int_t kNodeBasePort = 49152;  // BasePort of TCP ports on the nodes

  ClassDef(AliEveHOMERSrcTranslator, 0); // Translate HLT data-sources.
};

#endif
