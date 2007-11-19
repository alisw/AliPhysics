//-*- Mode: C++ -*-
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

#include "TString.h"
#include "TObject.h"

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
class AliHLTHOMERSourceDesc : public TObject {

public:

  /** standard constructor */
  AliHLTHOMERSourceDesc();       

  /** constructor 
   * @param hostname   hostname of the source
   * @param port       port of the source
   */
  AliHLTHOMERSourceDesc( TString hostname, Int_t port );

  /** destructor */
  virtual ~AliHLTHOMERSourceDesc();

  // -- SELECTION --

  /** Set selection state
   * @param state      state, either kTRUE or kFALSE
   */
  void SetState( Bool_t state ) { fSelected = state; }

  /** Checks if Source is selected to readout 
   * @return           returns state, either kTRUE or kFALSE
  */
  Bool_t IsSelected() { return fSelected; }

  /** Select this source */
  void Select() { fSelected = kTRUE; }

  /** Deselect this source */
  void Deselect() { fSelected = kFALSE; }

  // -- SETTER --

  /** Set node name of this source 
   * @param hostname   hostname of the source
   * @param port       port of the source
   */
   void SetHostnamePort( TString hostname, Int_t port ) { fHostname = hostname; fPort = port; }

  /** Set name/obj name of this source 
   *  @param s1        source name
   *  @param s2        source class name, default is ""
   */
  void SetSourceName( TString s1, TString s2="" ) { fSourceName = s1; fClassName = s2; }

  /** Set detector, sub detector and sub sub detector of this source 
   *  @param s1        detector name
   *  @param s2        subdetector name
   *  @param s3        subsubdetector name
   */
  void SetDetectors( TString s1, TString s2, TString s3 ) { fDetector = s1; fSubDetector = s2; fSubSubDetector = s3; }

  /** Set HLT specification anf HLT data type of this source 
   *  @param ul        HLT specification
   *  @param s         HLT data type
   */
  void SetSourceType( ULong_t ul, TString s ) { fSpecification = ul, fDataType = s; }

  // -- GETTER --

  /** Get node name of this source 
   * @return           hostname
   */
  TString GetHostname() { return fHostname; }

  /** Get node name of this source 
   * @return           port
   */
  Int_t GetPort() { return fPort; }

  /** Get name of this source 
   * @return     #include "TString.h"      name
   */
  TString GetSourceName() { return fSourceName; }

  /** Get object name of this source
   * @return           class name
   */
  TString GetClassName() { return fClassName; }

  /** Get detector of this source 
   * @return           detector
   */
  TString GetDetector() { return fDetector; }

  /** Get sub detector of this source 
   * @return           subdetector
   */
  TString GetSubDetector() { return fSubDetector; }

  /** Get sub sub detector of this source 
   * @return           subsubdetector
   */
  TString GetSubSubDetector() { return fSubSubDetector; }

  /** Get HLT data type of this source
   * @return           HLT data type
   */
  TString GetDataType() { return fDataType; }

  /** Get HLT specification of this source
   * @return           HLT specification
   */
  ULong_t GetSpecification() { return fSpecification; }


private:
  /** copy constructor prohibited */
  AliHLTHOMERSourceDesc(const AliHLTHOMERSourceDesc&);
  
  /** assignment operator prohibited */
  AliHLTHOMERSourceDesc& operator=(const AliHLTHOMERSourceDesc&);

  /** is selected to read out */
  Bool_t fSelected;             // see above

  /** Name of HOMER Node */
  TString fHostname;            // see above

  /** Name of HOMER port */
  Int_t fPort;                  // see above

  /** Name of Source */
  TString fSourceName;          // see above

  /** Object Name of Source */
  TString fClassName;           // see above

  /** Detector Name, e.g. PHOS 
   *  corresponds to HLT origin
   */
  TString fDetector;            // see above

  /** SubDetector Name e.g. MODULE */
  TString fSubDetector;         // see above

  /** SubSubDetector Name e.g. PARTITION */
  TString fSubSubDetector;      // see above

  /** HLT Specification */
  ULong_t fSpecification;       // see above

  /** HLT DataType */
  TString fDataType;            // see above

  ClassDef( AliHLTHOMERSourceDesc, 0 )
};

#endif
