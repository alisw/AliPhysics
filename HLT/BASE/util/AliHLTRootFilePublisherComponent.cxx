// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTRootFilePublisherComponent.cxx
    @author Matthias Richter, Jochen Thaeder
    @date   
    @brief  HLT file publisher component implementation. */

#include "AliHLTRootFilePublisherComponent.h"
//#include <TObjString.h>
//#include <TMath.h>
//#include <TFile.h>

#include "TList.h"
#include "TTree.h"
#include "TKey.h"
#include "TFile.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootFilePublisherComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTRootFilePublisherComponent::AliHLTRootFilePublisherComponent()
  :
  AliHLTFilePublisher(),
  fpCurrentEvent(NULL),
  fObjectName("") {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // Set file to ROOT-File
  SetIsRawFile( kFALSE );
}

// #################################################################################
AliHLTRootFilePublisherComponent::~AliHLTRootFilePublisherComponent() {
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const char* AliHLTRootFilePublisherComponent::GetComponentID() {
  // see header file for class documentation
  return "ROOTFilePublisher";
}

// #################################################################################
AliHLTComponent* AliHLTRootFilePublisherComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTRootFilePublisherComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTRootFilePublisherComponent::ScanArgument(Int_t argc, const char** argv) {
  // see header file for class documentation

  Int_t iResult = 0;

  TString argument = "";
  TString parameter = "";
  Int_t bMissingParam = 0;
  fObjectName = "";  // Reset this to the default: read all objects in the file.
  
  argument=argv[iResult];
  if (argument.IsNull()) return -EINVAL;

  // -objectname
  if ( !argument.CompareTo("-objectname") ) {
    if ( ! (bMissingParam=(++iResult>=argc)) ) {
      parameter = argv[iResult];
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      fObjectName = parameter;
    } 
  }
  else {
    HLTError("unknown argument %s", argument.Data());
    iResult = -EINVAL;    
  }
  
  if ( bMissingParam ) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult = -EPROTO;
  }

  return iResult;
}

 // #################################################################################
Int_t AliHLTRootFilePublisherComponent::GetEvent( const AliHLTComponentEventData& /*evtData*/,
						AliHLTComponentTriggerData& /*trigData*/,
						AliHLTUInt8_t* /*outputPtr*/, 
						AliHLTUInt32_t& size,
						vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation

  if ( !IsDataEvent() ) return 0;

  Int_t iResult=0;
  size=0;

  // -- Ptr to current event
  TObjLink *lnk = fpCurrentEvent;
  if ( lnk == NULL) {
    lnk = GetEventList()->FirstLink();
    fpCurrentEvent = lnk;
  }

  if ( lnk ) {
    EventFiles* pEventDesc = dynamic_cast<EventFiles*>( lnk->GetObject() );
    if (pEventDesc) {
    
      HLTDebug("publishing files for event %p", pEventDesc);
      TList& files=*pEventDesc; // type conversion operator defined
      TObjLink *flnk=files.FirstLink();

      while (flnk && iResult>=0) {

	FileDesc* pFileDesc=dynamic_cast<FileDesc*>(flnk->GetObject());

	if (not fOpenFilesAtStart) pFileDesc->OpenFile();
	TFile* pFile=NULL;

	if (pFileDesc && (pFile=*pFileDesc)!=NULL) {

	  for ( Int_t i = 0; i < pFile->GetListOfKeys()->GetEntries(); i++  ){
	    TKey * key= dynamic_cast<TKey*>( pFile->GetListOfKeys()->At(i) );

	    if ( fObjectName != "" ) {
	      if ( !( ((TString) key->GetName()).CompareTo(fObjectName) ) )
		PushBack( key->ReadObj(), *pFileDesc, *pFileDesc ); 
	    }
	    else 
	      PushBack( key->ReadObj(), *pFileDesc, *pFileDesc ); 	      
	      // above : type conversion operator defined for DataType and Spec
	  }

	  if (not fOpenFilesAtStart) pFileDesc->CloseFile();
	} else {
	  HLTError("no file available");
	  iResult=-EFAULT;
	}
	flnk = flnk->Next();
      }
    } else {
      HLTError("can not get event descriptor from list link");
      iResult=-EFAULT;
    }
  } else {
    iResult=-ENOENT;
  }
  if (iResult>=0 && fpCurrentEvent) fpCurrentEvent=fpCurrentEvent->Next();

  return iResult;
}
