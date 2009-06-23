// $Id: 
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
//         for The ALICE HLT Project.                                    *

//-*- Mode: C++ -*-

/** @file   AliEveHOMERXMLHandler.cxx
    @author Jochen Thaeder
    @date
    @brief  Src Translator of HomerManger
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#define use_aliroot
#define use_root
#define ROWHOUGHPARAMS
#define use_reconstruction
#define use_newio
#define ROOTVERSION    "unchecked"
#define ALIROOTVERSION "unchecked"
#define __ROOT__
#define USE_ALILOG
#define LINUX

#define EVE_DEBUG 1
// -- -- -- -- -- -- -- 
#include "AliEveHOMERSource.h"
#include "AliEveHOMERXMLHandler.h"
// -- -- -- -- -- -- -- 
#include "TString.h"
#include <TApplication.h>
#include "Riostream.h"
#include "TXMLAttr.h"
#include "TCollection.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"
// -- -- -- -- -- -- -- 
#include "AliLog.h"



//______________________________________________________________________________
//
// Manage connections to HLT data-sources.

ClassImp(AliEveHOMERXMLHandler)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliEveHOMERXMLHandler::AliEveHOMERXMLHandler( TString xmlFile ) :
    fXMLFile(xmlFile),
    fXMLParser(NULL),
    fRootNode(NULL),
    fSrcTranslator(NULL),
    fSourceList(NULL)
{
  // This Class should handle the HLT XML config file.
  // host the XML parser, and do all searching in the 
  // XML File
  
  Initialize();
}

//##################################################################################
AliEveHOMERXMLHandler::~AliEveHOMERXMLHandler() {
  // The destructor

  if ( fXMLParser )
    delete fXMLParser;
  fXMLParser = NULL;

  if ( fSrcTranslator != NULL )
    delete fSrcTranslator;
  fSrcTranslator = NULL;

}

//##################################################################################
Int_t AliEveHOMERXMLHandler::Initialize() {
  // Initialize the XML Parser, set the root node

  Int_t iResult = 0 ;

  // -- Initialize XML parser
  if ( fXMLParser != NULL )
    delete fXMLParser;
  fXMLParser = NULL;

  fXMLParser = new TDOMParser();
  fXMLParser->SetValidate( kFALSE );

  iResult = fXMLParser->ParseFile( fXMLFile );
  if ( iResult < 0 ) {
    iResult = 1;
    AliError( Form("Parsing file with error: %s", fXMLParser->GetParseCodeMessage( fXMLParser->GetParseCode() )) );
    return iResult;
  }
  
  // -- Set root node
  fRootNode = fXMLParser->GetXMLDocument()->GetRootNode();

  // -- Initialize Src Translator
  fSrcTranslator = new AliEveHOMERSrcTranslator( "GPN" );

  return iResult;
}

//##################################################################################
Int_t AliEveHOMERXMLHandler::FillSourceList(TList *srcList) {
  // Files the source list of HOMER sources

  fSourceList = srcList;

  Int_t iResult = 0;

  TXMLNode * node = NULL;
  TXMLNode * prevNode = fRootNode->GetChildren();

  // -- Loop over all nodes
  while ( ( node = prevNode->GetNextNode() ) ) {
    prevNode = node;

    // -- Find only "Process" nodes, otherwise continue to next node
    if ( strcmp( node->GetNodeName(), "Proc" ) != 0 )
      continue;

    // -- Get Attributes of current node
    TList *attrList = node->GetAttributes();
    TXMLAttr *attr = 0;
    TIter next(attrList);

    while ( ( attr = (TXMLAttr*)next() ) ) {

      // -- Find "ID" attribute, otherwise continue to next attribute
      if ( strcmp( attr->GetName(), "ID" ) != 0 )
	continue;

      TString nodeId( attr->GetValue() );

      // -- Find only TDS processes
      TObjArray * nodeIdTok = nodeId.Tokenize("_");

      for ( Int_t ii=0 ; ii < nodeIdTok->GetEntries() ; ii++ ) {

	if ( ! ( (TObjString*) nodeIdTok->At(ii) )->GetString().CompareTo("TDS") ) {

	  iResult = AddSourceTDS( node->GetChildren() );
	  if ( iResult ) {
	    AliError( Form("Error processing TDS process : %s", nodeId.Data()) );
	  }

	}

      } // for ( Int_t ii=0 ; ii < nodeIdTok->GetEntries() ; ii++ ) {

    } // while ( ( attr = (TXMLAttr*)next() ) ) {

  } // while ( ( node = prevNode->GetNextNode() ) ) {
  
  return iResult;
}

//##################################################################################
Int_t AliEveHOMERXMLHandler::AddSourceTDS( TXMLNode * xmlNode ) {
  // Get Information out of a TDS process in XML file
  // * param xmlNode   Pointer to childs of TDS node
  // * return          0 on sucess, > 0 on errorsee header file for class documentation

  Int_t iResult = 0;

  TXMLNode * attrNode = NULL;
  TXMLNode * prevNode = xmlNode;

  TString xmlHostname;
  TString xmlPort;

  TString hostname;
  Int_t port = 0;

  // -- Get hostname and port from TDS node out of XML
  while ( ( attrNode = prevNode->GetNextNode() ) ) {
    prevNode = attrNode;

    // -- Get port out of the commandline
    if ( strcmp( attrNode->GetNodeName(), "Cmd" ) == 0 ) {
      TString cmd( attrNode->GetText() );

      TObjArray * cmdTok = cmd.Tokenize(" ");
      xmlPort = ((TObjString*) cmdTok->At(2))->GetString();
    }
    // -- Get hostname
    else if ( strcmp( attrNode->GetNodeName(), "Node" ) == 0 )
      xmlHostname = attrNode->GetText();

  } // while ( ( attrNode = prevNode->GetNextNode() ) ) {

  TString xmlNodename = GetNodename( xmlHostname );

  // -- Resolve hostname and port information --
  iResult = fSrcTranslator->Translate( xmlNodename, xmlPort, hostname, port );
  if ( iResult ) {
    if ( iResult == 1 )
      { AliError( Form("Error resolving hostname : %s", xmlHostname.Data()) ); }
    else if ( iResult == 2 )
      { AliError( Form("Error resolving port : %s", xmlPort.Data()) ); }
    return iResult;
  }

  // -- Reset loop to TDS node
  prevNode = xmlNode;

  // -- Get Sources out of XML, resolve sources, add to sources List
  while ( ( attrNode = prevNode->GetNextNode() ) ) {
    prevNode = attrNode;

    // Find only "Parent" tags, otherwise continue to next tag
    if ( strcmp( attrNode->GetNodeName(), "Parent" ) != 0 )
      continue;

    TString xmlParent = attrNode->GetText();

    AliHLTHOMERSourceDesc * source = new AliHLTHOMERSourceDesc( hostname, port );

    if ( FillSourceInformation( xmlParent, source ) ) {
      AliError( Form("Error resolving source : %s", xmlParent.Data()) );
      iResult = 3;
      delete source;
    }
    else {
      fSourceList->Add( source );
#if EVE_DEBUG
      AliInfo( Form("New Source added : %s", xmlParent.Data()) );
#endif
    }

  } // while ( ( attrNode = prevNode->GetNextNode() ) ) {

  return iResult;
}

//##################################################################################
Int_t AliEveHOMERXMLHandler::FillSourceInformation( TString xmlParent, AliHLTHOMERSourceDesc *source ) {
  // Resolve information of source
  // * param xmlParent   ParentString out of the XML
  // * param source      Return the filled AliHLTHOMERSourceDesc object
  // * return            0 on sucess, 1 on errorsee header file for class documentation

  Int_t iResult = 0;

  if ( ! xmlParent.Contains( "_" ) ) {
    AliError( Form("Source %s could not be resolved", xmlParent.Data() ) );
    iResult = 1;

    return iResult;
  }

  // -- Get detector / subDetector / subSubDetector
  TObjArray * parentTokens = xmlParent.Tokenize("_");
  Int_t nParentTokens = parentTokens->GetEntries();

  TString detector = ((TObjString*) parentTokens->At(0) )->GetString();
  TString subDetector = "";
  TString subSubDetector = "";
  TString name = ((TObjString*) parentTokens->At(1) )->GetString();
  
  if ( nParentTokens == 3 )
    subDetector = ((TObjString*) parentTokens->At(2) )->GetString();
  else if ( nParentTokens == 4 ) {
    subDetector = ((TObjString*) parentTokens->At(2) )->GetString();
    subSubDetector = ((TObjString*) parentTokens->At(3) )->GetString();
  }

  // -- Apply detector corrections
  fSrcTranslator->ApplyDetectorCorrections( detector, subDetector );

  // -- Remove Leading '0' in sub detector and subsubdetector
  subDetector.Remove( TString::kLeading, '0' );
  subSubDetector.Remove( TString::kLeading, '0' );
  
  // -- Set detector / subDetector / subSubDetector
  source->SetDetectors( detector, subDetector, subSubDetector );  

  // -- Fill dataType / specification / className
  iResult = fSrcTranslator->FillSourceDesc( source, name );

#if EVE_DEBUG
  AliInfo( Form("Set Source %s , Type %s, ClassName %s .", name.Data(), 
		source->GetDataType().Data(), source->GetClassName().Data()) );
  AliInfo( Form("    Detector %s , SubDetector : %s, SubSubDetector %s .",
		detector.Data(), subDetector.Data(), subSubDetector.Data()) );
#endif

  return iResult;
}

//##################################################################################
TString AliEveHOMERXMLHandler::GetNodename( TString xmlHostname ) {
  // Get xml nodename out of xml hostname

  TString nodename;

  TXMLNode * node = NULL;
  TXMLNode * prevNode = fRootNode->GetChildren();

  while ((node = prevNode->GetNextNode()) != 0)
  {
    prevNode = node;

    // -- Find only "Node" nodes, otherwise continue
    if ( strcmp( node->GetNodeName(), "Node" ) != 0 )
      continue;

    // -- Get Attributes of current node
    TList *attrList = node->GetAttributes();
    TXMLAttr *attr = 0;
    TIter next(attrList);

    TString nodeId;

    // Get "nodeID" and "nodeName" of this "Node" node
    while ( ( attr = (TXMLAttr*)next() ) ) {
      if ( strcmp( attr->GetName(), "ID" ) == 0 )
	nodeId = attr->GetValue();
      else if ( strcmp( attr->GetName(), "hostname" ) == 0 )
	nodename = attr->GetValue();
    }

    // -- if this is not the correct "nodeID" continue
    if ( nodeId != xmlHostname )
      continue;

    break;

  } // while ( node = prevNode->GetNextNode() ) {

  return nodename;
}
