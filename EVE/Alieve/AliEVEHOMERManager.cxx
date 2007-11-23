/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE HLT Project.                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliEVEHOMERManager.cxx
    @author Jochen Thaeder
    @date   
    @brief  Manger for HOMER in offline
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliEVEHOMERManager.h"

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

#include "AliHLTHOMERLibManager.h"

#include "AliHLTHOMERSourceDesc.h"
#include "AliHLTHOMERBlockDesc.h"

#include "AliLog.h"

#include "TString.h"
#include <TApplication.h>
#include "Riostream.h"
#include "TXMLAttr.h"
#include "TCollection.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"

ClassImp(AliEVEHOMERManager)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor 
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
AliEVEHOMERManager::AliEVEHOMERManager( TString xmlFile ) :
  Reve::RenderElementList("AliEVEHOMERManager"),

  fLibManager(new AliHLTHOMERLibManager),
  fXMLFile(xmlFile),
  fXMLParser(NULL), 
  fRootNode(NULL),
  fSourceList(NULL),
  fReader(NULL),
  fBlockList(NULL),
  fNBlks(0),
  fEventID(0),
  fCurrentBlk(0),
  fConnected(kFALSE),
  fStateHasChanged(kTRUE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliEVEHOMERManager::AliEVEHOMERManager( const AliEVEHOMERManager& m) :
  Reve::RenderElementList(m.GetName(), m.GetTitle())
{
  // see header file for class documentation

  AliError( Form( "copy constructor to be tested." ) );
}

//##################################################################################
AliEVEHOMERManager& AliEVEHOMERManager::operator=( const AliEVEHOMERManager& ) {
  // see header file for class documentation

  AliError( Form( "assignment constructor to be tested." ) );
  return *this;
}

//##################################################################################
AliEVEHOMERManager::~AliEVEHOMERManager() {
  // see header file for class documentation

  if ( fLibManager ) {
    if ( fReader )
      fLibManager->DeleteReader(fReader);
    delete fLibManager;
    fLibManager = NULL;
    fReader = NULL;
  }

  if ( fXMLParser ) 
    delete fXMLParser;
  fXMLParser = NULL;

  if ( fSourceList != NULL )
    delete fSourceList;
  fSourceList = NULL;

  if ( fBlockList != NULL )
    delete fBlockList;
  fBlockList = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Source Handling
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
Int_t AliEVEHOMERManager::CreateHOMERSourcesList() {
  // see header file for class documentation

  // -- Initialize XML parser
  if ( fXMLParser != NULL )
    delete fXMLParser;
  fXMLParser = NULL;

  fXMLParser = new TDOMParser();
  fXMLParser->SetValidate( kFALSE );
  
  Int_t iResult = fXMLParser->ParseFile( fXMLFile );
  if ( iResult < 0 ) {
    iResult = 1;
    AliError( Form("Parsing file with error: %s", fXMLParser->GetParseCodeMessage( fXMLParser->GetParseCode() )) );
    return iResult;
  }

  // -- Initialize sources list
  DestroyElements();
  if ( fSourceList != NULL )
    delete fSourceList;
  fSourceList = NULL;
  
  fSourceList = new TList();
  fSourceList->SetOwner( kTRUE );

  // -- Set ROOT node
  fRootNode = fXMLParser->GetXMLDocument()->GetRootNode();
  
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
      if ( nodeId.BeginsWith( "TDS_" ) ) {
	iResult = GetTDSAttributes( node->GetChildren() );
	if ( iResult ) {
	  AliError( Form("Error processing TDS process : %s", nodeId.Data()) );
	}
      } 
    
    } // while ( ( attr = (TXMLAttr*)next() ) ) {

  } // while ( ( node = prevNode->GetNextNode() ) ) {

  // -- New SourceList has been created --> All Sources are new --> State has changed
  fStateHasChanged = kTRUE;

  TIter next(fSourceList);
  AliHLTHOMERSourceDesc* src = 0;
  while ((src = (AliHLTHOMERSourceDesc*) next())) {
    Reve::RenderElementObjPtr* re = new Reve::RenderElementObjPtr(src, kFALSE);
    re->SetRnrElNameTitle
      (Form("%s-%s-%s %s", src->GetDetector().Data(), src->GetSubDetector().Data(),
	    src->GetSubSubDetector().Data(), src->GetDataType().Data()),
       "Title");
    AddElement(re);
  }

  if ( iResult ) {
    AliWarning( Form("There have been errors, while creating the sources list.") );
  }
  else {
    AliInfo( Form("New sources list created.") );
  }

  return iResult;
}

//##################################################################################
void AliEVEHOMERManager::SetSourceState( AliHLTHOMERSourceDesc * source, Bool_t state ) {
  // see header file for class documentation

  if ( source->IsSelected() != state ) {
    source->SetState( state );
    fStateHasChanged = kTRUE;
  }

  return;
}

//##################################################################################
Int_t AliEVEHOMERManager::GetTDSAttributes( TXMLNode * xmlNode ) {
  // see header file for class documentation
  
  Int_t iResult = 0;

  TXMLNode * attrNode = NULL;
  TXMLNode * prevNode = xmlNode;

  TString xmlHostname = 0;
  TString xmlPort = 0;

  TString hostname = 0; 
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
  
  // -- Resolve hostname and port information
  iResult = ResolveHostPortInformation ( xmlHostname, xmlPort, hostname, port );
  if ( iResult == 1 ) {
    AliError( Form("Error resolving hostname : %s", xmlHostname.Data()) );
    return iResult;
  }
  else if ( iResult == 2 ) {AliInfo( Form("Connection established") );
    AliError( Form("Error resolving port : %s", xmlPort.Data()) );
    return iResult;
  }

  // -- Reset loop to TDS node
  prevNode = xmlNode;  

  // -- Get Sources out of XML, resolve sources, add to sources ListxmlHostname.Data()
  while ( ( attrNode = prevNode->GetNextNode() ) ) {
    prevNode = attrNode;
    
    // Find only "Parent" tags, otherwise continue to next tag
    if ( strcmp( attrNode->GetNodeName(), "Parent" ) != 0 ) 
      continue;

    TString xmlParent = attrNode->GetText();
      
    AliHLTHOMERSourceDesc * source = new AliHLTHOMERSourceDesc( hostname, port );
      
    if ( ResolveSourceInformation( xmlParent, source ) ) {
      iResult = 3;
      AliError( Form("Error resolving source : %s", xmlParent.Data()) );
      
      delete source;
    }
    else {
      fSourceList->Add( source );
      AliInfo( Form("New Source added : %s", xmlParent.Data()) );
    }
    
  } // while ( ( attrNode = prevNode->GetNextNode() ) ) {
  
  return iResult;
}

//##################################################################################
Int_t AliEVEHOMERManager::ResolveHostPortInformation ( TString xmlHostname, TString xmlPort, TString &hostname, Int_t &port ) {
  // see header file for class documentation

  Int_t iResult = 1;

  // *** Resolve hostname

  TXMLNode * node = NULL;
  TXMLNode * prevNode = fRootNode->GetChildren();

  while ( ( node = prevNode->GetNextNode() ) && iResult == 1 ) {
    prevNode = node;
    
    // -- Find only "Node" nodes, otherwise continue
    if ( strcmp( node->GetNodeName(), "Node" ) != 0 )
      continue;

    // -- Get Attributes of current node
    TList *attrList = node->GetAttributes();
    TXMLAttr *attr = 0;
    TIter next(attrList);
    
    TString nodeId = 0;
    TString nodeName = 0;
    
    // Get "nodeID" and "nodeName" of this "Node" node
    while ( ( attr = (TXMLAttr*)next() ) ) {
      if ( strcmp( attr->GetName(), "ID" ) == 0 )
	nodeId = attr->GetValue();
      else if ( strcmp( attr->GetName(), "hostname" ) == 0 )
	nodeName = attr->GetValue();
    } 

    // -- if this is not the correct nodeID continue
    if ( nodeId != xmlHostname )
      continue;

    // -- Set hostname
    hostname = nodeName;
    iResult = 0;

    break;

  } // while ( ( node = prevNode->GetNextNode() ) ) {

  if ( iResult ) {
    AliError( Form("Error resolving hostname : %s", xmlHostname.Data()) );
    return iResult;
  }

  // *** Resolve port

  if ( xmlPort.IsDigit() ) {
    port = xmlPort.Atoi();
  }
  else {
    AliError ( Form("Error resolving port : %s", xmlPort.Data()) );
    iResult = 2;
  }

  // *** Summary

  if ( !iResult ) {
    AliInfo( Form("%s:%i resolved out of %s:%s", hostname.Data(), port, xmlHostname.Data(), xmlPort.Data()) );    
  }

  return iResult;
}

//##################################################################################
Int_t AliEVEHOMERManager::ResolveSourceInformation( TString xmlParent, AliHLTHOMERSourceDesc *source ) {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( ! xmlParent.Contains( "_" ) ) {
    AliError( Form("Source %s could not be resolved", xmlParent.Data() ) );
    iResult = 1;

    return iResult;
  }
    
  TObjArray * parentTokens = xmlParent.Tokenize("_");
    
  Int_t nEntries = parentTokens->GetEntries();
  
  TString detector = ((TObjString*) parentTokens->At(0) )->GetString();
  TString subDetector = "";
  TString subSubDetector = "";
  TString dataType = "";
  ULong_t specification = 0;

  TString name = ((TObjString*) parentTokens->At(1) )->GetString();
  TString objName = "";
  
  if ( nEntries == 3 ) 
    subDetector = ((TObjString*) parentTokens->At(2) )->GetString(); 
  else if ( nEntries == 4 ) {
    subDetector = ((TObjString*) parentTokens->At(2) )->GetString(); 
    subSubDetector = ((TObjString*) parentTokens->At(3) )->GetString(); 
  }
  
  // -- Corecct TPC subdetector, because in we have somtimes "A","C"
  if ( ! detector.CompareTo("TPC") ) {
    if ( subDetector.BeginsWith('A') ) {
      subDetector.Remove( TString::kLeading, 'A' );
    }
    else if ( subDetector.BeginsWith('C') ) {
      subDetector.Remove( TString::kLeading, 'C' );
      Int_t tmp = subDetector.Atoi() + 18;
      subDetector = "";
      subDetector += tmp;
    }
  }

  // -- Remove Leading '0' in sub detector and subsubdetector
  subDetector.Remove( TString::kLeading, '0' );
  subSubDetector.Remove( TString::kLeading, '0' );

  // -- Set Object Names

  // **** General ****
  if ( name == "RP" || name == "FP" ) {
    objName = "";
    dataType = "DDL_RAW";
    specification = 0;
  } 

  // **** TPC ****
  else if ( detector == "TPC" ) {
    
    if ( name == "CalibPedestal" ) {
      objName = "AliTPCCalibPedestal";
      dataType = "HIS_CAL";
      specification = 0;
    }
    else if ( name == "CalibPulser" ) {
      objName = "AliTPCCalibPulser"; 
      dataType = "HIS_CAL";
      specification = 0;
    }
  } // if ( detector == "TPC" ) {

  // **** TRD ****
  else if ( detector == "TRD" ) {

    if ( name == "foo" ) {
      objName = "bar";
      dataType = "FOO_BAR";
      specification = 0;
    }
  } // else if ( detector == "TRD" ) {

  // **** PHOS ****
  else if ( detector == "PHOS" ) {

  } // else if ( detector == "PHOS" ) {

  // **** DIMU ****
  else if ( detector == "MUON" ) {

  } // else if ( detector == "MUON" ) {
  
  // -- Fill object
  source->SetSourceName( name, objName );
  source->SetSourceType( specification, dataType );
  source->SetDetectors( detector, subDetector, subSubDetector );


  AliInfo( Form("Set Source %s , Type %s, ClassName %s .", name.Data(), dataType.Data(), objName.Data()) );

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Connection Handling
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
Int_t AliEVEHOMERManager::ConnectHOMER(){
  // see header file for class documentation

  Int_t iResult = 0;    

  // -- Check if already connected and state has not changed
  if ( fStateHasChanged == kFALSE && IsConnected() ) {
    AliInfo( Form("No need for reconnection.") );
    return iResult;
  }
    
  // -- If already connected, disconnect before connect
  if ( IsConnected() ) 
    DisconnectHOMER();

  // *** Create the Readoutlist

  UShort_t* sourcePorts = new UShort_t [fSourceList->GetEntries()];
  const char ** sourceHostnames = new const char* [fSourceList->GetEntries()];
  UInt_t sourceCount = 0;

  CreateReadoutList( sourceHostnames, sourcePorts, sourceCount );
  
  if ( sourceCount == 0 ) {
    AliError(Form("No sources selected, aborting.") );
    return iResult;
  }

  // *** Connect to data sources
  
  if ( !fReader ) {
    if ( fLibManager ) 
      fReader = fLibManager->OpenReader( sourceCount, sourceHostnames, sourcePorts );
  }
  
  iResult = fReader->GetConnectionStatus();
  
  if ( iResult ) {	
    // -- Connection failed 
      
    UInt_t ndx = fReader->GetErrorConnectionNdx();
      
    if ( ndx < sourceCount ) {
      AliError( Form("Error : Error establishing connection to TCP source %s:%hu: %s (%d)", 
		     sourceHostnames[ndx], sourcePorts[ndx], strerror(iResult), iResult) );
    }
    else {
      AliError( Form("Error : Error establishing connection to unknown source with index %d: %s (%d)",
		     ndx, strerror(iResult), iResult) );
    }

    if ( fReader )
      fLibManager->DeleteReader( fReader );
    fReader = NULL;

  } 
  else {
    // -- Connection ok - set reader
    fConnected = kTRUE; 
    
    AliInfo( Form("Connection established") );
  }
  
  delete[] sourceHostnames;
  delete[] sourcePorts;
  

  // -- Get next event
  if ( ! iResult )
    NextEvent();
  
  return iResult;
}

//##################################################################################
void AliEVEHOMERManager::DisconnectHOMER(){
  // see header file for class documentation

  if ( ! IsConnected() )
    return;

  if ( fReader )
    fLibManager->DeleteReader( fReader );
  fReader = NULL;

  fStateHasChanged = kTRUE;
  fConnected = kFALSE;
  
  AliInfo( Form("Connection closed") );

  return;
}

//##################################################################################
Int_t AliEVEHOMERManager::ReconnectHOMER(){
  // see header file for class documentation

  Int_t iResult = 0;

  if ( IsConnected() )
    DisconnectHOMER();

  iResult = ConnectHOMER();
  if ( iResult ) {
    AliError( Form("Error connecting.") );
  }

  return iResult;
}


//##################################################################################
void AliEVEHOMERManager::CreateReadoutList( const char** sourceHostnames, UShort_t *sourcePorts, UInt_t &sourceCount ){
  // see header file for class documentation

  AliHLTHOMERSourceDesc * source= NULL;

  // -- Read all sources and check if they should be read out
  TIter next( fSourceList );
  while ( ( source = (AliHLTHOMERSourceDesc*)next() ) ) {

    if ( ! source->IsSelected() ) 
      continue;
    
    Bool_t exists = kFALSE;

    // -- Loop over existing entries and check if entry is already in readout list
    for ( UInt_t ii = 0; ii < sourceCount; ii++ ){
      if ( !strcmp( sourceHostnames[ii], source->GetHostname().Data() ) &&  sourcePorts[ii] == source->GetPort() ) {
	exists = kTRUE;
	break;
      }
    }
    
    // -- Add new entires to readout list
    if ( ! exists ) {
      sourcePorts[sourceCount] = source->GetPort();
      sourceHostnames[sourceCount] = source->GetHostname().Data();
      sourceCount++;
    }

  } // while ( ( source = (AliHLTHOMERSourceDesc*)next() ) ) {

  fStateHasChanged = kFALSE;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Event Handling
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
Int_t AliEVEHOMERManager::NextEvent(){
  // see header file for class documentation
  
  Int_t iResult = 0;
  
  if ( !fReader || ! IsConnected() ) {
    AliWarning( Form( "Not connected yet." ) );
    return 1;
  }
  
  // -- Read next event data and error handling for HOMER (error codes and empty blocks)
  while( 1 ) {
    iResult = fReader->ReadNextEvent( 5000000 /*timeout in us*/);
    
    if ( iResult == 111 || iResult == 32 || iResult == 6 ) {
      Int_t ndx = fReader->GetErrorConnectionNdx();
      AliError( Form("Error, No Connection to source %d: %s (%d)", ndx, strerror(iResult), iResult) );
 
     return 2;
    }
    else if ( iResult == 110 ) {
      Int_t ndx = fReader->GetErrorConnectionNdx();
      AliError( Form("Timout occured, reading event from source %d: %s (%d)", ndx, strerror(iResult), iResult) );
      return 3; 
    }
    else if ( iResult == 56) {
      Int_t ndx = fReader->GetErrorConnectionNdx();
      AliError( Form("Error reading event from source %d: %s (%d) -- IRESULTRY", ndx, strerror(iResult), iResult) );
      continue; 
    }
    else if ( iResult ) {
      Int_t ndx = fReader->GetErrorConnectionNdx();
      AliError( Form("General Error reading event from source %d: %s (%d)", ndx, strerror(iResult), iResult) );
      return 2;
    }
    else {
      break;
    }
  } // while( 1 ) {

  if ( iResult ) 
    return iResult;


  // -- Get blockCnt and eventID
  fNBlks = (ULong_t) fReader->GetBlockCnt();
  fEventID = (ULong64_t) fReader->GetEventID();
  fCurrentBlk = 0;

  AliInfo( Form("Event 0x%016LX (%Lu) with %lu blocks", fEventID, fEventID, fNBlks) );

#if 0
  /*
  // Loop for Debug only
  for ( ULong_t i = 0; i < fNBlks; i++ ) {
    Char_t tmp1[9], tmp2[5];
    memset( tmp1, 0, 9 );
    memset( tmp2, 0, 5 );
    void *tmp11 = tmp1;
    ULong64_t* tmp12 = (ULong64_t*)tmp11;
    *tmp12 = fReader->GetBlockDataType( i );
    void *tmp21 = tmp2;
    ULong_t* tmp22 = (ULong_t*)tmp21;
    *tmp22 = fReader->GetBlockDataOrigin( i );
    AliInfo( Form("Block %lu length: %lu - type: %s - origin: %s",i, fReader->GetBlockDataLength( i ), tmp1, tmp2) );
  } // end for ( ULong_t i = 0; i < fNBlks; i++ ) {
  */
#endif

  // -- Create BlockList
  CreateBlockList();

    return iResult;
}

//##################################################################################
Int_t AliEVEHOMERManager::CreateBlockList() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Initialize block list
  if ( fBlockList != NULL )
    delete fBlockList;
  fBlockList = NULL;

  fBlockList = new TList();
  fBlockList->SetOwner( kTRUE );

  void* iter = GetFirstBlk();

  // -- Fill block list
  while ( iter != NULL ){
    
    // -- Create new block
    AliHLTHOMERBlockDesc * block = new AliHLTHOMERBlockDesc( GetBlk(), GetBlkSize(), GetBlkOrigin(), 
							     GetBlkType(), GetBlkSpecification() );
    
    // -- Check sources list if block is requested
    if ( CheckIfRequested( block ) )
      fBlockList->Add( block );      
    else 
      delete block;
    
    iter = GetNextBlk();
    
  } // while ( iter != NULL ){
  
  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                            BlockHandling
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
void* AliEVEHOMERManager::GetBlk( Int_t ndx ) {
  // see header file for class documentation

  void* data = NULL;
  
  if ( !fReader || ! IsConnected() ) {
    AliError( Form("Not connected yet.") ); 
  }
  else {
    if ( ( ndx ) < (Int_t) fNBlks ) 
      data = (void*) fReader->GetBlockData( ndx );
  }

  return data;
}

//##################################################################################
ULong_t AliEVEHOMERManager::GetBlkSize( Int_t ndx ) {
  // see header file for class documentation

  ULong_t length = 0;

  if ( !fReader || ! IsConnected() ) {
    AliError( Form("Not connected yet.") ); 
  }
  else {
    if ( ( ndx ) < (Int_t) fNBlks ) 
      length = (ULong_t) fReader->GetBlockDataLength( ndx );
  }

  return length;
}

//##################################################################################
TString AliEVEHOMERManager::GetBlkOrigin( Int_t ndx ) {
  // see header file for class documentation

  TString origin = "";

  // -- Check for Connection
  if ( !fReader || ! IsConnected() ) {
    AliError( Form("Not connected yet.") ); 
    return origin;
  }

  // -- Check block index
  if ( ( ndx ) >= (Int_t) fNBlks ) {
    AliError( Form("Block index %d out of range.", ndx ) ); 
    return origin;  
  }

  // -- Get origin
  union{
    UInt_t data;
    Char_t array[4];
  } reverseOrigin;

  reverseOrigin.data = (UInt_t) fReader->GetBlockDataOrigin( ndx );

  // -- Reverse the order
  for (Int_t ii = 3; ii >= 0; ii-- )
    if ( reverseOrigin.array[ii] != ' ')
      origin.Append( reverseOrigin.array[ii] );

  return origin;
}

//##################################################################################
TString AliEVEHOMERManager:: GetBlkType( Int_t ndx ) {
  // see header file for class documentation

  TString type = "";

  // -- Check for Connection
  if ( !fReader || ! IsConnected() ) {
    AliError( Form("Not connected yet.") ); 
    return type;
  }

  // -- Check blockk index
  if ( ( ndx ) >= (Int_t) fNBlks ) {
    AliError( Form("Block index %d out of range.", ndx ) ); 
    return type;  
  }

  // -- Get type
  union{
    ULong64_t data;
    Char_t array[8];
  } reverseType;

  reverseType.data = (ULong64_t) fReader->GetBlockDataType( ndx );

  // -- Reverse the order
  for (Int_t ii = 7; ii >= 0; ii-- )
    if ( reverseType.array[ii] != ' ')
      type.Append( reverseType.array[ii] );

  return type;
}


//##################################################################################
ULong_t AliEVEHOMERManager:: GetBlkSpecification( Int_t ndx ) {
  // see header file for class documentation

  ULong_t spec = 0;


  // -- Check for Connection
  if ( !fReader || ! IsConnected() ) {
    AliError( Form("Not connected yet.") ); 
    return spec;
  }

  // -- Check blockk index
  if ( ( ndx ) >= (Int_t) fNBlks ) {
    AliError( Form("Block index %d out of range.", ndx ) ); 
    return spec;  
  }
  
  spec = (ULong_t) fReader->GetBlockDataSpec( ndx );

  return spec;
}

//##################################################################################
Bool_t AliEVEHOMERManager::CheckIfRequested( AliHLTHOMERBlockDesc * block ) {
  // see header file for class documentation

  Bool_t requested = kFALSE;
  
  AliHLTHOMERSourceDesc * source= NULL;
  
  // -- Read all sources and check if they should be read out
  TIter next( fSourceList );
  while ( ( source = (AliHLTHOMERSourceDesc*)next() ) ) {

    if ( ! source->IsSelected() ) 
      continue;

    if ( source->GetDetector().CompareTo( block->GetDetector() ) )
      continue;

    if ( source->GetDataType().CompareTo( block->GetDataType() ) )
      continue;

    if ( ! block->HasSubDetectorRange() ) {
      
      if ( source->GetSubDetector().CompareTo( block->GetSubDetector() ) )
	continue;
      
      if ( ! block->HasSubSubDetectorRange() ) {
	
	if ( source->GetSubSubDetector().CompareTo( block->GetSubSubDetector() ) )
	   continue;
	 
      } // if ( ! block->HasSubSubDetectorRange ) {
    } //  if ( ! block->HasSubDetectorRange ) {
    
    requested = kTRUE;
   
  } // while ( ( source = (AliHLTHOMERSourceDesc*)next() ) ) {
  
  if ( requested) {
    AliInfo( Form("Block requested : %s - %s : %s/%s -> %s ", block->GetDetector().Data(), block->GetDataType().Data(), 
		  block->GetSubDetector().Data(), block->GetSubSubDetector().Data(), block->GetClassName().Data() ) );
  }

  return requested;
}

//##################################################################################
void AliEVEHOMERManager::TestSelect() {
  // see header file for class documentation

  for (Int_t ii =0; ii < fSourceList->GetEntries() ; ii++ ) {
    if ( (ii%2) == 0 )
      ((AliHLTHOMERSourceDesc*) fSourceList->At(ii))->Select();
  }
}


