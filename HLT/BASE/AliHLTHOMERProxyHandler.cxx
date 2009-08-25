//-*- Mode: C++ -*-
// $Id: AliHLTHOMERProxyHandler.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTHOMERProxyHandler.cxx
    @author Jochen Thaeder
    @date
    @brief  HOMER proxy handler for HomerManger
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TDOMParser.h"
#include "TSocket.h"
#include "TSystem.h"
// -- -- -- -- -- -- -- 
#include "AliHLTHOMERProxyHandler.h"
// -- -- -- -- -- -- -- 

ClassImp(AliHLTHOMERProxyHandler)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTHOMERProxyHandler::AliHLTHOMERProxyHandler() :
  fRealm(kHLT),
  fXmlRpcResponse(""),
  fSourceList(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliHLTHOMERProxyHandler::~AliHLTHOMERProxyHandler() {
  // see header file for class documentation

}

//##################################################################################
Int_t AliHLTHOMERProxyHandler::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0 ;

  IdentifyRealm();

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Source List - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERProxyHandler::FillSourceList(TList *srcList) {
  // see header file for class documentation

  Int_t iResult = 0;

  fSourceList = srcList;

  iResult = RequestXmlRpcResponse();

  if (!iResult)
    iResult = ProcessXmlRpcResponse();

  if (iResult < 0) {
    HLTError(Form("Filling SourceList failed."));
  }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                        Realms - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
const Char_t *AliHLTHOMERProxyHandler::fgkHOMERProxyNode[] = { 
  "portal-dcs0.internal", 
  "alihlt-dcs0.cern.ch",
  "alihlt-vobox0.cern.ch",
  "alihlt-gw0.kip.uni-heidelberg.de",
  "portal-dcs1.internal", 
  "alihlt-dcs1.cern.ch",
  "alihlt-vobox1.cern.ch",
  "alihlt-gw1.kip.uni-heidelberg.de"
};

//##################################################################################
void AliHLTHOMERProxyHandler::IdentifyRealm() {
  // see header file for class documentation

  TString hostIP(gSystem->GetHostByName(gSystem->HostName()).GetHostAddress());

  if ( hostIP.Contains("10.162.") )
    fRealm = kHLT;
  else if ( hostIP.Contains("10.160.") || hostIP.Contains("10.161.") )
    fRealm = kACR;
  else if ( hostIP.Contains("129.206.") )
    fRealm = kKIP;
  else 
    fRealm = kGPN;
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                        Proxy Communication - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERProxyHandler::RequestXmlRpcResponse() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- open socket
  // ----------------

  Int_t proxyPort = 19999;

  TSocket *socket = new TSocket(fgkHOMERProxyNode[fRealm], proxyPort);
  if ( ! socket->IsValid() ) {
    HLTWarning(Form("Failed to create socket to %s:%d,",fgkHOMERProxyNode[fRealm], proxyPort));
    HLTWarning(Form("trying %s:%d now.", fgkHOMERProxyNode[fRealm+kHOMERRealmsMax],proxyPort));

    socket = new TSocket(fgkHOMERProxyNode[fRealm+kHOMERRealmsMax], proxyPort);
    if ( ! socket->IsValid() ) {
      HLTError(Form("Failed to create socket to %s:%d and %s:%d.",
		    fgkHOMERProxyNode[fRealm], proxyPort,
		    fgkHOMERProxyNode[fRealm+kHOMERRealmsMax],proxyPort));

      fRealm = -1;
      return -1;
    }
    else
      fRealm = fRealm+kHOMERRealmsMax;
  }

  // -- send request
  // -----------------

  Char_t reqMsg[] = "PUT / HTTP/1.1\r\n\
User-Agent: curl/7.18.0 (x86_64-pc-linux-gnu) libcurl/7.18.0 OpenSSL/0.9.8g zlib/1.2.3.3 libidn/1.1\r\n\
Host: localhost:10000\r\n\
Accept: */*\r\n\
Content-type: text/xml\r\n\
Content-Length: 68\r\n\
\r\n<methodCall><methodName>getTcpDumpServices</methodName></methodCall>\r\n";

  iResult = socket->SendRaw( reqMsg, strlen(reqMsg) );
  if ( iResult < 1 || 
       iResult !=  static_cast<Int_t>(strlen(reqMsg))) {
    HLTError(Form("Error sending! -- send length %d  -- msg length %d.", iResult, static_cast<Int_t>(strlen(reqMsg)) ));
    socket->Close();
    return iResult;
  }
  
  // -- receive answer
  // -------------------

  const Int_t bufferSize = 1024;
  Char_t buffer[bufferSize];

  Bool_t isXmlRpc = kFALSE;

  fXmlRpcResponse = "";

  // -- loop for getting full xmlRPC response
  while(1) {

    Int_t bufferLength = 0;
    
    // -- loop for reading until end of line
    while (1) {

      iResult = socket->RecvRaw(&buffer[bufferLength], 1);
      if ( iResult < 0) {
	HLTError(Form("Error reading form socket."));
	socket->Close();
	return iResult;
      }
            
      // -- Checking for  end of line
      if ( buffer[bufferLength] == 10 ) {
	buffer[bufferLength] = 0;
	break;
      }
      
      ++bufferLength;
    }

    TString bufferString(buffer);

    // -- Checking for start of XML response
    if ( bufferString.BeginsWith("<?xml") )
      isXmlRpc = kTRUE;

    // -- Append the xml response
    if (isXmlRpc) {
      fXmlRpcResponse.Append(bufferString);
    }

    // -- Checking for end of XML response
    if( ! bufferString.CompareTo("</methodResponse>") ) 
      break;
  }
  
  // -- close socket
  socket->Close();

  return 0;
}

//##################################################################################
Int_t AliHLTHOMERProxyHandler::ProcessXmlRpcResponse() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Parse XML RPC Response
  // ---------------------------

  TDOMParser xmlParser;
  xmlParser.SetValidate(kFALSE);

  HLTDebug(Form("XMLResponse: %s",fXmlRpcResponse.Data()));
  
  iResult = xmlParser.ParseBuffer(fXmlRpcResponse.Data(), fXmlRpcResponse.Length());
  if ( iResult < 0 ) {
    HLTError(Form("Parsing buffer with error: %s", 
		  xmlParser.GetParseCodeMessage(xmlParser.GetParseCode()) ));
    

    return iResult;
  }

  TXMLNode * node = xmlParser.GetXMLDocument()->GetRootNode()->
    GetChildren()->GetChildren()->GetChildren()->GetChildren();
  
  if ( strcmp( node->GetNodeName(), "string" ) ) {
    HLTError(Form("No node 'string' in XmlRpcResponse."));
    return -1;
  }

  // -- Parse Content
  // ------------------

  // -- Get Content
  TString xmlContent(node->GetText() );

  HLTDebug(Form("XMLContent: %s",xmlContent.Data()));

  iResult = xmlParser.ParseBuffer(xmlContent.Data(), xmlContent.Length());
  if ( iResult < 0 ) {
    HLTError(Form("Parsing buffer with error: %s", 
		  xmlParser.GetParseCodeMessage(xmlParser.GetParseCode()) ));

    return iResult;
  }
  
  if ( !xmlParser.GetXMLDocument()->GetRootNode()->HasChildren() ) {
    HLTWarning(Form("No Services active."));
    return 1;
  }

  // -- Loop over all service nodes
  TXMLNode* serviceNode = xmlParser.GetXMLDocument()->GetRootNode()->GetChildren();
  TXMLNode* prevServiceNode = NULL;
  
  do {
    prevServiceNode = serviceNode;

    // -- Add service to list
    iResult = AddService( serviceNode->GetChildren() );
    if ( iResult > 0 ) {
      HLTWarning(Form("Incomplete Service not added."));
      iResult = 0;
    }
  } while ( ( serviceNode = prevServiceNode->GetNextNode() ) && !iResult );



  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Source Resolving - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERProxyHandler::AddService(TXMLNode *innerNode) {
  // see header file for class documentation

  Int_t iResult = 0;

  HLTInfo(Form(">> New service"));    

  TXMLNode* serviceNode = innerNode;

  // -- Loop over all service properties and 
  //    read them from the service tag
  // -----------------------------------------

  TString hostname          = "";
  Int_t   port              = 0;
  TString dataType          = "";
  TString dataOrigin        = "";
  TString dataSpecification = "";

  TXMLNode* prevInnerNode = NULL;

  // -- Retrieve hostname and port
  // -------------------------------

  do {
    prevInnerNode = innerNode;
    
    if ( ! strcmp(innerNode->GetNodeName(), "text" ) )
      continue;
        
    // -- hostname
    if ( ! strcmp( innerNode->GetNodeName(), "address") ) {
      HLTInfo(Form("  > %s ++ %s", innerNode->GetNodeName(), innerNode->GetText() ));
      hostname = innerNode->GetText();
    }
    else if ( ! strcmp( innerNode->GetNodeName(), "port") ) {
      HLTInfo(Form("  > %s ++ %s", innerNode->GetNodeName(), innerNode->GetText() ));
      TString portS(innerNode->GetText());
      if ( portS.IsDigit() )
	port = portS.Atoi();
      else {
	HLTError(Form("Port %s is not a digit.", portS.Data()));
	iResult = -1;
      }
    }
  } while ( ( innerNode = prevInnerNode->GetNextNode() ) && !iResult );


  // -- Change hostame from service with proxy, if outside HLT
  if ( fRealm != kHLT || fRealm != kHLT+kHOMERRealmsMax )
    hostname = fgkHOMERProxyNode[fRealm];


  // -- Get Data Specifications from blocks
  // ----------------------------------------

  do {
    prevInnerNode = serviceNode;

    if ( strcmp( serviceNode->GetNodeName(), "blocks") )
      continue;
 
    TXMLNode* blocks = serviceNode->GetChildren();

    if ( ! blocks ) {
      HLTError(Form("No blocks present"));
      return 1;
    }
      
    TXMLNode* blockNode = blocks->GetNextNode();
    TXMLNode* prevBlockNode = NULL;

    if ( ! blockNode ) {
      HLTError(Form("No block present in the blocks tag"));
      return 1;
    }
      
    // -- blocks loop 
    
    do {
      prevBlockNode = blockNode;
      
      if ( strcmp( blockNode->GetNodeName(), "block") )
	continue;

      TXMLNode *dataNode = blockNode->GetChildren();
      TXMLNode *prevDataNode = NULL;

      if ( ! dataNode ) {
	HLTError(Form("No data specification tags present in block tag."));
	return 1;
      }
      // -- data spec loop
      
      do {
	prevDataNode = dataNode;

	if ( ! strcmp(dataNode->GetNodeName(), "text" ) )
	  continue;

	HLTInfo(Form(" %s ++ %s", dataNode->GetNodeName(), dataNode->GetText() ));	

	if ( ! strcmp( dataNode->GetNodeName(), "dataorigin") ) {
	  dataOrigin = dataNode->GetText();
	}
	else if ( ! strcmp( dataNode->GetNodeName(), "datatype") ) {
	  dataType = dataNode->GetText();
	}
	else if ( ! strcmp( dataNode->GetNodeName(), "dataspecification") ) {
	  dataSpecification = dataNode->GetText();    
	}
      } while ( ( dataNode = prevDataNode->GetNextNode() ) && !iResult );
      
      // -- data spec loop

      // -- Check the service properties
      // ---------------------------------
      
      // -- Check for completeness of the source properties
      if ( hostname.IsNull() || !port || dataOrigin.IsNull() ||
	   dataType.IsNull() || dataSpecification.IsNull() ) {
	HLTWarning(Form("Service provides not all values:\n\thostname\t\t %s\n\tport\t\t\t %d\n\tdataorigin\t\t %s\n\tdatatype\t\t %s\n\tdataspecification\t 0x%08X", 
			hostname.Data(), port, dataOrigin.Data(), dataType.Data(), dataSpecification.Atoi()));
	
	return 1;
      }

      // -- Create new source
      // ----------------------
      
      AliHLTHOMERSourceDesc * source = new AliHLTHOMERSourceDesc();
      source->SetService( hostname, port, dataOrigin, dataType, dataSpecification );
      
      fSourceList->Add( source );
      
      HLTInfo(Form( "New Source added : %s", source->GetSourceName().Data()));

    } while ( ( blockNode = prevBlockNode->GetNextNode() ) && !iResult );


    // -- blocks loop
    
  } while ( ( serviceNode = prevInnerNode->GetNextNode() ) && !iResult );
  
  return iResult;
}


