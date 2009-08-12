//-*- Mode: C++ -*-

// $Id: AliHLTHOMERProxyHandler.h $

#ifndef ALIHLTHOMERPROXYHANDLER_H
#define ALIHLTHOMERPROXYHANDLER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliHLTHOMERProxyHandler.h
    @author Jochen Thaeder
    @date
    @brief  HOMER proxy handler for HomerManger
*/

#include "TString.h"
#include "TList.h"
#include "TXMLNode.h" 
// -- -- -- -- -- -- -- 
#include "AliHLTHOMERSourceDesc.h"
// -- -- -- -- -- -- -- 
#include "AliHLTLoggingVariadicFree.h"

/**
 * @class AliHLTHOMERProxyHandler
 * This Class should handle the communication with the proxy
 * and fill the source list.
 *
 * @ingroup alihlt_homer
 */

class AliHLTHOMERProxyHandler : public TObject, public AliHLTLogging
{
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTHOMERProxyHandler();

  /** destructor */
  virtual ~AliHLTHOMERProxyHandler();

  /** Initialize 
   *  @return 0 on success, <0 for failure
   */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                             Source List - public
   * ---------------------------------------------------------------------------------
   */

  /** Fill's source list, with entries 
   *  @return 0 on success, <0 for failure, 1 for no active service
   */
  Int_t FillSourceList(TList *srcList);

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTHOMERProxyHandler(const AliHLTHOMERProxyHandler&);            // Not implemented.

  /** assignment operator prohibited */
  AliHLTHOMERProxyHandler& operator=(const AliHLTHOMERProxyHandler&); // Not implemented.

  /*
   * ---------------------------------------------------------------------------------
   *                        Realms - private
   * ---------------------------------------------------------------------------------
   */

  /** Realms */
  enum HOMERRealms_t { 
    kHLT,            /**< HLT realm */
    kACR,            /**< ACR realm */
    kGPN,            /**< GPN realm */
    kKIP,            /**< KIP realm */
    kHOMERRealmsMax  /**< Number of enum entries */
  };

  /** Array of proxy nodes per realm */
  static const Char_t *fgkHOMERProxyNode[];     

  /** Indentifies the realm and sets it
   *  @return 0 on success, <0 for failure
   */
  void IdentifyRealm();
  
  /*
   * ---------------------------------------------------------------------------------
   *                        Proxy Communication - private
   * ---------------------------------------------------------------------------------
   */

  /** Get xmlrpc response from the proxy 
   *  @return 0 on success, <0 for failure
   */
  Int_t RequestXmlRpcResponse();

  /** process xmlrpc response and fill the source list
   *  @return 0 on success, <0 for failure, 1 for no active service
   */
  Int_t ProcessXmlRpcResponse();
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Source Resolving - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Add a new Service to list
   *  @param xmlNode   Ptr to service node
   *  @return          0 on sucess, <0 for failure
   */
  Int_t AddService(TXMLNode *innerNode);
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Realm, which can be ACR, GPN, HLT, KIP */
  Int_t       fRealm;                             // see above

  /** xmlRPC response */
  TString     fXmlRpcResponse;                    // see above

  /** List to HOMER sources */
  TList*      fSourceList;                        //! transient

  ClassDef(AliHLTHOMERProxyHandler, 0); // Handles HLT xml sources.
};

#endif
