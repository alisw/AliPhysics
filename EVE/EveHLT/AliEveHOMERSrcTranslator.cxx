// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
//         for The ALICE HLT Project.                                    *

//-*- Mode: C++ -*-

/** @file   AliEveHOMERSrcTranslator.cxx
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
#include "AliEveHOMERSrcTranslator.h"
#include "AliEveHOMERSrcObject.h"
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
// Translate HLT data-sources.

ClassImp(AliEveHOMERSrcTranslator)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliEveHOMERSrcTranslator::AliEveHOMERSrcTranslator( TString realm ) :
    fBasePortMap(NULL),
    fObjectMap(NULL),
    fRealm(realm)
{
  // This Class should handle the translation of 
  // internal hostnames and ports to the ones used by
  // HOMER according to realm, where AliEVE is running in.
  
  SetupPortMap();
  SetupObjectMap();
}

//##################################################################################
AliEveHOMERSrcTranslator::~AliEveHOMERSrcTranslator() {
  // The destructor

  if ( fBasePortMap )
    delete fBasePortMap;
  fBasePortMap = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                           Translation - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliEveHOMERSrcTranslator::Translate( TString xmlNodename, TString xmlPort, 
					   TString &hostname, Int_t &port ) {
  // Translate hostname and port for source which has to be used by HOMER
  // ( due to port mapping inside the HLT )
  // * param xmlNodename  Nodename out of the XML
  // * param xmlPort      Port out of the XML
  // * param hostname     Return of the hostname
  // * param port         Return of the port
  // * return             0 on sucess, 1 if port couldn't be resolved,

  Int_t iResult = 0;

  // *** Resolve hostname
  hostname = ResolveHostname( xmlNodename );
  
  // *** Resolve port
  port = ResolvePort( xmlPort, xmlNodename );
  
  if ( port == -1 ) {
    AliError( Form("Error resolving port : %s", xmlPort.Data()) );
    iResult = 1;
  }

  // *** Summary
#if EVE_DEBUG
  if ( !iResult ) {
    AliInfo( Form("%s:%i resolved out of %s:%s", hostname.Data(), port, xmlNodename.Data(), xmlPort.Data()) );
  }
#endif

  return iResult;
}

//##################################################################################
void AliEveHOMERSrcTranslator::ApplyDetectorCorrections( TString &detector, TString &subDetector) {
  // Apply corrections for differnt detectors and subdetectors */

  // -- Correct TPC subdetector, because in we have somtimes "A","C"
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
  
  // -- Correct for MUON
  if ( ! detector.CompareTo("DIMU") ) {
    detector = "MUON";
    
    if ( ! subDetector.CompareTo("TRG") )
      subDetector = "1";
    else if ( ! subDetector.CompareTo("TRK") )
      subDetector = "2";
  }
}

//##################################################################################
Int_t AliEveHOMERSrcTranslator::FillSourceDesc( AliHLTHOMERSourceDesc* source, TString name ) {
  // Fill SourceDesc with object Information

  Int_t iResult = 0;

  TString detector = source->GetDetector();

  if ( ! fObjectMap->FindObject( detector ) ) {
    AliError( Form("Error mapping for detector not known : %s", detector.Data()) );
    AliError( Form("Error mapping for NAME : %s", name.Data()) );
    iResult = 1;

    return iResult;
  }


  TMap * objectMap = (TMap*) fObjectMap->GetValue( detector );



  if ( ! objectMap->FindObject( name ) ) {
  cout << "DET..." << detector.Data() << endl;    

    source->SetSourceName( name, "" );
    source->SetSourceType( 0, "*******" );
  }
  else {
    AliEveHOMERSrcObject* srcObject = (AliEveHOMERSrcObject*) objectMap->FindObject( name );
    source->SetSourceName( name, srcObject->GetClassName() );
    source->SetSourceType( srcObject->GetSpecification(), srcObject->GetDataType() );
  }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Source Resolving - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
TString AliEveHOMERSrcTranslator::ResolveHostname( TString nodename ){
  // resolves the hostname, out of the XML nodename, and the realm set

  TString hostname;

  // -- Set hostname according to realm
  if ( ! fRealm.CompareTo( "ACR" ) ) 
    hostname = "alihlt-dcs0.cern.ch";
  else if ( ! fRealm.CompareTo( "GPN" ) ) 
    hostname = "alihlt-vobox0.cern.ch";
  else if ( ! fRealm.CompareTo( "KIP" ) ) 
    hostname = "alihlt-gw0.kip.uni-heidelberg.de";
  else 
    hostname = nodename;

  return hostname;
}

//##################################################################################
Int_t AliEveHOMERSrcTranslator::ResolvePort( TString srcPort, TString nodename ) {
  // resolves the port, out of the XML port, and the realm set
  Int_t port = -1;

  if ( ! srcPort.IsDigit() )
    return port;

  if ( srcPort.Atoi() < kNodeBasePort )
    return port;
  
  if ( ! fBasePortMap->FindObject( nodename ) ) 
    return port;

  port = srcPort.Atoi();
  
  if ( ! fRealm.CompareTo( "HLT" ) )
    return port;

  Int_t offset = port - kNodeBasePort;

  port = ( ( (TObjString*) fBasePortMap->GetValue(nodename) )->GetString() ).Atoi();
  port += offset;

  return port;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Setup - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliEveHOMERSrcTranslator::SetupPortMap () {
  // Setup hostname to port mapping

  fBasePortMap = new TMap();	
  //fBasePortMap->SetOwnerKeyValue();

  fBasePortMap->Add( new TObjString("cntpca000"), new TObjString("49408"));
  fBasePortMap->Add( new TObjString("feptpcao00"), new TObjString("49436"));
  fBasePortMap->Add( new TObjString("feptpcai00"), new TObjString("49440"));
  fBasePortMap->Add( new TObjString("feptpcao01"), new TObjString("49444"));
  fBasePortMap->Add( new TObjString("feptpcao02"), new TObjString("49564"));
  fBasePortMap->Add( new TObjString("feptpcai02"), new TObjString("49568"));
  fBasePortMap->Add( new TObjString("feptpcao03"), new TObjString("49572"));
  fBasePortMap->Add( new TObjString("cntpca040"), new TObjString("49664"));
  fBasePortMap->Add( new TObjString("feptpcao04"), new TObjString("49692"));
  fBasePortMap->Add( new TObjString("feptpcai04"), new TObjString("49696"));
  fBasePortMap->Add( new TObjString("feptpcao05"), new TObjString("49700"));
  fBasePortMap->Add( new TObjString("feptpcao06"), new TObjString("49820"));
  fBasePortMap->Add( new TObjString("feptpcai06"), new TObjString("49824"));
  fBasePortMap->Add( new TObjString("feptpcao07"), new TObjString("49828"));
  fBasePortMap->Add( new TObjString("cntpca080"), new TObjString("49920"));
  fBasePortMap->Add( new TObjString("feptpcao08"), new TObjString("49948"));
  fBasePortMap->Add( new TObjString("feptpcai08"), new TObjString("49952"));
  fBasePortMap->Add( new TObjString("feptpcao09"), new TObjString("49956"));
  fBasePortMap->Add( new TObjString("feptpcao10"), new TObjString("50076"));
  fBasePortMap->Add( new TObjString("feptpcai10"), new TObjString("50080"));
  fBasePortMap->Add( new TObjString("feptpcao11"), new TObjString("50084"));
  fBasePortMap->Add( new TObjString("cntpca120"), new TObjString("50176"));
  fBasePortMap->Add( new TObjString("feptpcao12"), new TObjString("50204"));
  fBasePortMap->Add( new TObjString("feptpcai12"), new TObjString("50208"));
  fBasePortMap->Add( new TObjString("feptpcao13"), new TObjString("50212"));
  fBasePortMap->Add( new TObjString("feptpcao14"), new TObjString("50332"));
  fBasePortMap->Add( new TObjString("feptpcai14"), new TObjString("50336"));
  fBasePortMap->Add( new TObjString("feptpcao15"), new TObjString("50340"));
  fBasePortMap->Add( new TObjString("cntpca160"), new TObjString("50432"));
  fBasePortMap->Add( new TObjString("feptpcao16"), new TObjString("50460"));
  fBasePortMap->Add( new TObjString("feptpcai16"), new TObjString("50464"));
  fBasePortMap->Add( new TObjString("feptpcao17"), new TObjString("50468"));
  fBasePortMap->Add( new TObjString("cntrd0"), new TObjString("54144"));
  fBasePortMap->Add( new TObjString("feptrd00"), new TObjString("54168"));
  fBasePortMap->Add( new TObjString("feptrd04"), new TObjString("54172"));
  fBasePortMap->Add( new TObjString("feptrd08"), new TObjString("54176"));
  fBasePortMap->Add( new TObjString("feptrd10"), new TObjString("54180"));
  fBasePortMap->Add( new TObjString("feptrd14"), new TObjString("54184"));
  fBasePortMap->Add( new TObjString("feptpcco16"), new TObjString("54428"));
  fBasePortMap->Add( new TObjString("feptpcci16"), new TObjString("54432"));
  fBasePortMap->Add( new TObjString("feptpcco17"), new TObjString("54436"));
  fBasePortMap->Add( new TObjString("feptpcco14"), new TObjString("54556"));
  fBasePortMap->Add( new TObjString("feptpcci14"), new TObjString("54560"));
  fBasePortMap->Add( new TObjString("feptpcco15"), new TObjString("54564"));
  fBasePortMap->Add( new TObjString("cntpcc120"), new TObjString("54656"));
  fBasePortMap->Add( new TObjString("feptpcco12"), new TObjString("54684"));
  fBasePortMap->Add( new TObjString("feptpcci12"), new TObjString("54688"));
  fBasePortMap->Add( new TObjString("feptpcco13"), new TObjString("54692"));
  fBasePortMap->Add( new TObjString("feptpcco10"), new TObjString("54812"));
  fBasePortMap->Add( new TObjString("feptpcci10"), new TObjString("54816"));
  fBasePortMap->Add( new TObjString("feptpcco11"), new TObjString("54820"));
  fBasePortMap->Add( new TObjString("cntpcc080"), new TObjString("54912"));
  fBasePortMap->Add( new TObjString("feptpcco08"), new TObjString("54940"));
  fBasePortMap->Add( new TObjString("feptpcci08"), new TObjString("54944"));
  fBasePortMap->Add( new TObjString("feptpcco09"), new TObjString("54948"));
  fBasePortMap->Add( new TObjString("feptpcco06"), new TObjString("55068"));
  fBasePortMap->Add( new TObjString("feptpcci06"), new TObjString("55072"));
  fBasePortMap->Add( new TObjString("feptpcco07"), new TObjString("55076"));
  fBasePortMap->Add( new TObjString("cntpcc040"), new TObjString("55168"));
  fBasePortMap->Add( new TObjString("feptpcco04"), new TObjString("55196"));
  fBasePortMap->Add( new TObjString("feptpcci04"), new TObjString("55200"));
  fBasePortMap->Add( new TObjString("feptpcco05"), new TObjString("55204"));
  fBasePortMap->Add( new TObjString("feptpcco02"), new TObjString("55324"));
  fBasePortMap->Add( new TObjString("feptpcci02"), new TObjString("55328"));
  fBasePortMap->Add( new TObjString("feptpcco03"), new TObjString("55332"));
  fBasePortMap->Add( new TObjString("cntpcc000"), new TObjString("55424"));
  fBasePortMap->Add( new TObjString("feptpcco00"), new TObjString("55452"));
  fBasePortMap->Add( new TObjString("feptpcci00"), new TObjString("55456"));
  fBasePortMap->Add( new TObjString("feptpcco01"), new TObjString("55460"));
  fBasePortMap->Add( new TObjString("cnspd0"), new TObjString("57728"));
  fBasePortMap->Add( new TObjString("fepspd0"), new TObjString("57752"));
  fBasePortMap->Add( new TObjString("fepspd1"), new TObjString("57756"));
  fBasePortMap->Add( new TObjString("fepspd2"), new TObjString("57760"));
  fBasePortMap->Add( new TObjString("fepspd3"), new TObjString("57764"));
  fBasePortMap->Add( new TObjString("fepspd4"), new TObjString("57768"));
  fBasePortMap->Add( new TObjString("cnssd0"), new TObjString("57856"));
  fBasePortMap->Add( new TObjString("fepssd0"), new TObjString("57880"));
  fBasePortMap->Add( new TObjString("fepssd1"), new TObjString("57884"));
  fBasePortMap->Add( new TObjString("fepssd2"), new TObjString("57888"));
  fBasePortMap->Add( new TObjString("fepssd3"), new TObjString("57892"));
  fBasePortMap->Add( new TObjString("fepphos4"), new TObjString("57896"));
  fBasePortMap->Add( new TObjString("fepsdd5"), new TObjString("57896"));
  fBasePortMap->Add( new TObjString("fephmpid0"), new TObjString("58008"));
  fBasePortMap->Add( new TObjString("fepsdd0"), new TObjString("58008"));
  fBasePortMap->Add( new TObjString("fephmpid1"), new TObjString("58012"));
  fBasePortMap->Add( new TObjString("fepsdd1"), new TObjString("58012"));
  fBasePortMap->Add( new TObjString("fephmpid2"), new TObjString("58016"));
  fBasePortMap->Add( new TObjString("fepsdd2"), new TObjString("58016"));
  fBasePortMap->Add( new TObjString("fephmpid3"), new TObjString("58020"));
  fBasePortMap->Add( new TObjString("fepsdd3"), new TObjString("58020"));
  fBasePortMap->Add( new TObjString("fepphos1"), new TObjString("58024"));
  fBasePortMap->Add( new TObjString("fepsdd4"), new TObjString("58024"));
  fBasePortMap->Add( new TObjString("feptriggerdet"), new TObjString("58140"));
  fBasePortMap->Add( new TObjString("fepfmdaccorde"), new TObjString("58144"));
  fBasePortMap->Add( new TObjString("fephltout0"), new TObjString("58400"));
  fBasePortMap->Add( new TObjString("fephltout1"), new TObjString("58404"));
  fBasePortMap->Add( new TObjString("cnphos0"), new TObjString("58624"));
  fBasePortMap->Add( new TObjString("fepphos2"), new TObjString("58656"));
  fBasePortMap->Add( new TObjString("fepphos3"), new TObjString("58660"));
  fBasePortMap->Add( new TObjString("cndimutrg0"), new TObjString("58752"));
  fBasePortMap->Add( new TObjString("fepdimutrg"), new TObjString("58784"));
  fBasePortMap->Add( new TObjString("cndimutrk0"), new TObjString("58880"));
  fBasePortMap->Add( new TObjString("fepdimutrk1"), new TObjString("58904"));
  fBasePortMap->Add( new TObjString("fepdimutrk2"), new TObjString("58908"));
  fBasePortMap->Add( new TObjString("fepdimutrk3"), new TObjString("58912"));
  fBasePortMap->Add( new TObjString("fepdimutrk4"), new TObjString("58916"));
  fBasePortMap->Add( new TObjString("fepdimutrk5"), new TObjString("58920"));
 

//   fBasePortMap->Add( new TObjString("feptpcao00"), new TObjString("49436"));
//   fBasePortMap->Add( new TObjString("feptpcai00"), new TObjString("49440"));
//   fBasePortMap->Add( new TObjString("feptpcao01"), new TObjString("49444"));
//   fBasePortMap->Add( new TObjString("feptpcao02"), new TObjString("49564"));
//   fBasePortMap->Add( new TObjString("feptpcai02"), new TObjString("49568"));
//   fBasePortMap->Add( new TObjString("feptpcao03"), new TObjString("49572"));
//   fBasePortMap->Add( new TObjString("feptpcao04"), new TObjString("49692"));
//   fBasePortMap->Add( new TObjString("feptpcai04"), new TObjString("49696"));
//   fBasePortMap->Add( new TObjString("feptpcao05"), new TObjString("49700"));
//   fBasePortMap->Add( new TObjString("feptpcao06"), new TObjString("49820"));
//   fBasePortMap->Add( new TObjString("feptpcai06"), new TObjString("49824"));
//   fBasePortMap->Add( new TObjString("feptpcao07"), new TObjString("49828"));
//   fBasePortMap->Add( new TObjString("feptpcao08"), new TObjString("49948"));
//   fBasePortMap->Add( new TObjString("feptpcai08"), new TObjString("49952"));
//   fBasePortMap->Add( new TObjString("feptpcao09"), new TObjString("49956"));
//   fBasePortMap->Add( new TObjString("feptpcao10"), new TObjString("50076"));
//   fBasePortMap->Add( new TObjString("feptpcai10"), new TObjString("50080"));
//   fBasePortMap->Add( new TObjString("feptpcao11"), new TObjString("50084"));
//   fBasePortMap->Add( new TObjString("feptpcao12"), new TObjString("50204"));
//   fBasePortMap->Add( new TObjString("feptpcai12"), new TObjString("50208"));
//   fBasePortMap->Add( new TObjString("feptpcao13"), new TObjString("50212"));
//   fBasePortMap->Add( new TObjString("feptpcao14"), new TObjString("50332"));
//   fBasePortMap->Add( new TObjString("feptpcai14"), new TObjString("50336"));
//   fBasePortMap->Add( new TObjString("feptpcao15"), new TObjString("50340"));
//   fBasePortMap->Add( new TObjString("feptpcao16"), new TObjString("50460"));
//   fBasePortMap->Add( new TObjString("feptpcai16"), new TObjString("50464"));
//   fBasePortMap->Add( new TObjString("feptpcao17"), new TObjString("50468"));
//   fBasePortMap->Add( new TObjString("feptrd00"), new TObjString("54168"));
//   fBasePortMap->Add( new TObjString("feptrd04"), new TObjString("54172"));
//   fBasePortMap->Add( new TObjString("feptrd08"), new TObjString("54176"));
//   fBasePortMap->Add( new TObjString("feptrd10"), new TObjString("54180"));
//   fBasePortMap->Add( new TObjString("feptrd14"), new TObjString("54184"));
//   fBasePortMap->Add( new TObjString("feptpcco16"), new TObjString("54428"));
//   fBasePortMap->Add( new TObjString("feptpcci16"), new TObjString("54432"));
//   fBasePortMap->Add( new TObjString("feptpcco17"), new TObjString("54436"));
//   fBasePortMap->Add( new TObjString("feptpcco14"), new TObjString("54556"));
//   fBasePortMap->Add( new TObjString("feptpcci14"), new TObjString("54560"));
//   fBasePortMap->Add( new TObjString("feptpcco15"), new TObjString("54564"));
//   fBasePortMap->Add( new TObjString("feptpcco12"), new TObjString("54684"));
//   fBasePortMap->Add( new TObjString("feptpcci12"), new TObjString("54688"));
//   fBasePortMap->Add( new TObjString("feptpcco13"), new TObjString("54692"));
//   fBasePortMap->Add( new TObjString("feptpcco10"), new TObjString("54812"));
//   fBasePortMap->Add( new TObjString("feptpcci10"), new TObjString("54816"));
//   fBasePortMap->Add( new TObjString("feptpcco11"), new TObjString("54820"));
//   fBasePortMap->Add( new TObjString("feptpcco08"), new TObjString("54940"));
//   fBasePortMap->Add( new TObjString("feptpcci08"), new TObjString("54944"));
//   fBasePortMap->Add( new TObjString("feptpcco09"), new TObjString("54948"));
//   fBasePortMap->Add( new TObjString("feptpcco06"), new TObjString("55068"));
//   fBasePortMap->Add( new TObjString("feptpcci06"), new TObjString("55072"));
//   fBasePortMap->Add( new TObjString("feptpcco07"), new TObjString("55076"));
//   fBasePortMap->Add( new TObjString("feptpcco04"), new TObjString("55196"));
//   fBasePortMap->Add( new TObjString("feptpcci04"), new TObjString("55200"));
//   fBasePortMap->Add( new TObjString("feptpcco05"), new TObjString("55204"));
//   fBasePortMap->Add( new TObjString("feptpcco02"), new TObjString("55324"));
//   fBasePortMap->Add( new TObjString("feptpcci02"), new TObjString("55328"));
//   fBasePortMap->Add( new TObjString("feptpcco03"), new TObjString("55332"));
//   fBasePortMap->Add( new TObjString("feptpcco00"), new TObjString("55452"));
//   fBasePortMap->Add( new TObjString("feptpcci00"), new TObjString("55456"));
//   fBasePortMap->Add( new TObjString("feptpcco01"), new TObjString("55460"));
//   fBasePortMap->Add( new TObjString("fepspd0"), new TObjString("57752"));
//   fBasePortMap->Add( new TObjString("fepspd1"), new TObjString("57756"));
//   fBasePortMap->Add( new TObjString("fepspd2"), new TObjString("57760"));
//   fBasePortMap->Add( new TObjString("fepspd3"), new TObjString("57764"));
//   fBasePortMap->Add( new TObjString("fepspd4"), new TObjString("57768"));
//   fBasePortMap->Add( new TObjString("fepssd0"), new TObjString("57880"));
//   fBasePortMap->Add( new TObjString("fepssd1"), new TObjString("57884"));
//   fBasePortMap->Add( new TObjString("fepssd2"), new TObjString("57888"));
//   fBasePortMap->Add( new TObjString("fepssd3"), new TObjString("57892"));
//   fBasePortMap->Add( new TObjString("feptriggerdet"), new TObjString("58140"));
//   fBasePortMap->Add( new TObjString("fepfmdaccorde"), new TObjString("58144"));
//   fBasePortMap->Add( new TObjString("fephmpid0"), new TObjString("58264"));
//   fBasePortMap->Add( new TObjString("fephmpid1"), new TObjString("58268"));
//   fBasePortMap->Add( new TObjString("fephmpid2"), new TObjString("58272"));
//   fBasePortMap->Add( new TObjString("fephmpid3"), new TObjString("58276"));
//   fBasePortMap->Add( new TObjString("fephltout0"), new TObjString("58400"));
//   fBasePortMap->Add( new TObjString("fephltout1"), new TObjString("58404"));
//   fBasePortMap->Add( new TObjString("fepphos2"), new TObjString("58656"));
//   fBasePortMap->Add( new TObjString("fepphos3"), new TObjString("58660"));
//   fBasePortMap->Add( new TObjString("fepphos4"), new TObjString("58664"));
//   fBasePortMap->Add( new TObjString("fepdimutrg"), new TObjString("58784"));
//   fBasePortMap->Add( new TObjString("fepdimutrk1"), new TObjString("58904"));
//   fBasePortMap->Add( new TObjString("fepdimutrk2"), new TObjString("58908"));
//   fBasePortMap->Add( new TObjString("fepdimutrk3"), new TObjString("58912"));
//   fBasePortMap->Add( new TObjString("fepdimutrk4"), new TObjString("58916"));
//   fBasePortMap->Add( new TObjString("fepdimutrk5"), new TObjString("58920"));

}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMap () {
  // Setup hostname to port mapping

  fObjectMap = new TMap();	
  // fObjectMap->SetOwnerKeyValue();

  SetupObjectMapTPC();
  SetupObjectMapTRD();
  SetupObjectMapPHOS();
  SetupObjectMapDIMUON();
}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMapTPC() {
  //Setup the Object mapping for TPC

  TMap* objectMap =  new TMap();
  fObjectMap->Add( new TObjString("TPC"), objectMap );
  
  SetupObjectMapCommon( objectMap );

  objectMap->Add( new TObjString("CF"),            new AliEveHOMERSrcObject( "AliHLTTPCClusterDataFormat", "CLUSTERS", 0 ) );
  objectMap->Add( new TObjString("RelayCF"),       new AliEveHOMERSrcObject( "AliHLTTPCClusterDataFormat", "CLUSTERS", 0 ) );
  objectMap->Add( new TObjString("CalibPedestal"), new AliEveHOMERSrcObject( "AliTPCCalibPedestal", "HIS_CAL", 0 ) );
  objectMap->Add( new TObjString("CalibPulser"),   new AliEveHOMERSrcObject( "AliTPCCalibPulser", "HIS_CAL", 0 ) );
  objectMap->Add( new TObjString("ESDConv"),       new AliEveHOMERSrcObject( "TTree", "ESD_TREE", 0 ) );
  objectMap->Add( new TObjString("ESDCM"),         new AliEveHOMERSrcObject( "TTree", "ESD_TREE", 0 ) );
  objectMap->Add( new TObjString("ESDCA"),         new AliEveHOMERSrcObject( "TTree", "ESD_TREE", 0 ) );
  objectMap->Add( new TObjString("RelayESD"),      new AliEveHOMERSrcObject( "TTree", "ESD_TREE", 0 ) );
  objectMap->Add( new TObjString("KRCF"),          new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );
  objectMap->Add( new TObjString("RelayKR"),       new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );
  objectMap->Add( new TObjString("CLHI"),          new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );
  objectMap->Add( new TObjString("RelayCLHI"),     new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );
  objectMap->Add( new TObjString("NM"),            new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );
  objectMap->Add( new TObjString("HH"),            new AliEveHOMERSrcObject( "TH1F", "ROOTHIST", 0 ) );

}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMapTRD(){
  //Setup the Object mapping for TRD

  TMap* objectMap =  new TMap();
  fObjectMap->Add( new TObjString("TRD"), objectMap );

  SetupObjectMapCommon( objectMap );

}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMapPHOS(){
  //Setup the Object mapping for PHOS

  TMap* objectMap =  new TMap();
  fObjectMap->Add( new TObjString("PHOS"), objectMap );

  SetupObjectMapCommon( objectMap );


}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMapDIMUON(){
  //Setup the Object mapping for DIMUON

  TMap* objectMap =  new TMap();
  fObjectMap->Add( new TObjString("MUON"), objectMap );

  objectMap->Add( new TObjString("RECHITS"),        new AliEveHOMERSrcObject( "", "RECHITS", 0 ) );
  objectMap->Add( new TObjString("TRIGRECS"),       new AliEveHOMERSrcObject( "", "TRIGRECS", 0 ) );
  objectMap->Add( new TObjString("DECIDSIN"),       new AliEveHOMERSrcObject( "", "DECIDSIN", 0 ) );
  objectMap->Add( new TObjString("DECIDPAR"),       new AliEveHOMERSrcObject( "", "DECIDPAR", 0 ) );
  objectMap->Add( new TObjString("MANTRACK"),       new AliEveHOMERSrcObject( "", "MANTRACK", 0 ) );


  SetupObjectMapCommon( objectMap );
  

}

//##################################################################################
void AliEveHOMERSrcTranslator::SetupObjectMapCommon( TMap* objectMap) {
  // Setup the common Object mappings

  objectMap->Add( new TObjString("RP"),    new AliEveHOMERSrcObject( "", "DDL_RAW", 0 ) );
  objectMap->Add( new TObjString("FP"),    new AliEveHOMERSrcObject( "", "DDL_RAW", 0 ) );
  objectMap->Add( new TObjString("Relay"), new AliEveHOMERSrcObject( "", "DDL_RAW", 0 ) );
}



