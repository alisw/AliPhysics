// $Id: 
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//-*- Mode: C++ -*-
#ifndef ALIEVEHOMERXMLHANDLER_H
#define ALIEVEHOMERXMLHANDLER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliEveHOMERXMLHandler.h
    @author Jochen Thaeder
    @date
    @brief  XML Handler for HomerManger
*/

#include "TString.h"
#include "TDOMParser.h"
#include "TXMLNode.h"
#include "TList.h"

#include "AliEveHOMERSrcTranslator.h"
#include "AliHLTHOMERSourceDesc.h"


class AliEveHOMERXMLHandler : public TObject
{
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliEveHOMERXMLHandler( TString xmlFile );

  /** destructor */
  virtual ~AliEveHOMERXMLHandler();

  /*
   * ---------------------------------------------------------------------------------
   *                             Source List - public
   * ---------------------------------------------------------------------------------
   */

  /** Fill's source list, with entries */
  Int_t FillSourceList(TList *srcList);

  /** Sets realm ( which can be ACR, GPN, HLT, KIP ) */ 
  void SetRealm( TString s ) { fSrcTranslator->SetRealm(s); } // Sets realm ( which can be ACR, GPN, HLT, KIP )


  ///////////////////////////////////////////////////////////////////////////////////

private:

  AliEveHOMERXMLHandler(const AliEveHOMERXMLHandler&);            // Not implemented.
  AliEveHOMERXMLHandler& operator=(const AliEveHOMERXMLHandler&); // Not implemented.

  /*
   * ---------------------------------------------------------------------------------
   *                            Source Resolving - private
   * ---------------------------------------------------------------------------------
   */

  /** Initialize the XML Handler, create the SrcTranslator */
  Int_t Initialize();

  /** Get Information out of a TDS process in XML file */
  Int_t AddSourceTDS( TXMLNode * xmlNode );

  /** Get xml nodename out of xml hostname */
  TString GetNodename( TString xmlHostname );

  /** Resolve information of source */
  Int_t FillSourceInformation( TString xmlParent, AliHLTHOMERSourceDesc * source );

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  // == XML parser ==
  TString     fXMLFile;                           //  XML input file
  TDOMParser* fXMLParser;                         //! XML parser into DOM model
  TXMLNode*   fRootNode;                          //! Root node of parsed config file

  // == Source Translator ==
  AliEveHOMERSrcTranslator* fSrcTranslator;       //! Translates HOMER sources to real hostname,port

  // == Sources List ==
  TList* fSourceList;                             //! List to HOMER sources

  ClassDef(AliEveHOMERXMLHandler, 0); // Handles HLT xml sources.
};

#endif
