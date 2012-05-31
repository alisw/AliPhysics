// $Id: AliHLTESDMCEventPublisherComponent.cxx 25452 2008-04-25 21:46:04Z richterm $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTESDMCEventPublisherComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Component for publishing ESD and MC events.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTESDMCEventPublisherComponent.h"
#include "AliHLTErrorGuard.h"

#include "TList.h"
#include "TTree.h"
#include "TKey.h"
#include "TFile.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTESDMCEventPublisherComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTESDMCEventPublisherComponent::AliHLTESDMCEventPublisherComponent()
  :
  AliHLTFilePublisher(),
  fpCurrentFolder(NULL),
  fpCurrentFileList(NULL),
  fCurrentEvent(0),
  fNEventsInFolder(0),
  fFolderList(),
  fSpecification(0),
  fPublishESD(kFALSE),
  fPublishHLTESD(kFALSE),
  fPublishMC(kFALSE),
  fFastMC(kFALSE),
  fpTreeESD(NULL),
  fpTreeHLTESD(NULL),
  fpTreeE(NULL),
  fpTreeK(NULL),
  fpTreeTR(NULL),
  fpESD(NULL),
  fpHLTESD(NULL),
  fpESDClone(NULL),
  fpMC(NULL),
  fpHLTMC(NULL),
  fOutputSize(0),
  fApplyParticleCuts(kFALSE),
  fSkippedEsdObjects() {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // Set file to ROOT-File
  SetIsRawFile( kFALSE );

  fFolderList.SetOwner();
}

// #################################################################################
AliHLTESDMCEventPublisherComponent::~AliHLTESDMCEventPublisherComponent() {
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
const char* AliHLTESDMCEventPublisherComponent::GetComponentID() {
  // see header file for class documentation
  return "ESDMCEventPublisher";
}

// #################################################################################
AliHLTComponent* AliHLTESDMCEventPublisherComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTESDMCEventPublisherComponent;
}

void AliHLTESDMCEventPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation
  constBase=fOutputSize;
  inputMultiplier=0.0;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::DoInit(int argc, const char** argv) {
  // see header file for class documentation

  Int_t iResult = 0;
  Int_t bMissingParam=0;
  
  TString argument="";

  fSpecification = kAliHLTVoidDataSpec;

  // -- Loop over all arguments
  for ( Int_t iter = 0; iter<argc && iResult>=0; iter++) {
    argument=argv[iter];
    if (argument.IsNull()) continue;
    
    // -- Check for which type to publish
    //    Can be one, all or some of : ESD, HLTESD and MC
    if ( !argument.CompareTo("-entrytype") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');
      
      if ( !parameter.CompareTo("ESD") )
	fPublishESD = kTRUE;
      else if ( !parameter.CompareTo("HLTESD") )
	fPublishHLTESD = kTRUE;
      else if ( !parameter.CompareTo("MC") )
	fPublishMC = kTRUE;
      else if ( !parameter.CompareTo("MCFAST") ) {
	fPublishMC = kTRUE;
	fFastMC = kTRUE;
      }
      else {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    } 
    
    // -- Data path to reconstructed data
    else if ( ! argument.CompareTo("-datapath") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');
      
      fFolderList.Add( new TObjString(parameter) );
    }

    // -- Data Specification 
    //    -- see header file !!!
    else if ( ! argument.CompareTo("-dataspec") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');
            
      if ( parameter.IsDigit() )
	fSpecification = (AliHLTUInt32_t) parameter.Atoi();
      else if ( parameter.BeginsWith("0x") && parameter.Replace(0,2,"",0).IsHex() )
	sscanf(parameter.Data(),"%x", &fSpecification);
      else {
	HLTError("Wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    }
    // -- applyParticleCuts
    else if ( !argument.CompareTo("-applyParticleCuts") ) {
      fApplyParticleCuts = kTRUE;
    } 

   else if ( !argument.CompareTo("-skip-esd-object") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      if (!fSkippedEsdObjects.IsNull()) fSkippedEsdObjects+=" ";
      fSkippedEsdObjects+=argv[iter];
    } 

    // -- Argument not known
    else {
      HLTError("Unknown argument %s.", argument.Data());
      iResult = -EINVAL;
    }
    
  } // for ( Int iter = 0; iter<argc && iResult>=0; iter++) {
  
  // -- Check if parameter is missing
  if ( bMissingParam ) {
    HLTError("Missing parameter for argument %s.", argument.Data());
     iResult=-EINVAL;
  }
  
  //  -- Check if at least on entrytype was specified
  if ( fPublishESD||fPublishHLTESD||fPublishMC ) 
    AddDataTypesToOutputlist();
  else {
    HLTError("At least one -entrytype argument needs to be specified.");
    iResult=-EINVAL;
  }
  
  // -- Check if at least one datapath was specified
  if ( !fFolderList.GetEntries() ) {
    HLTError("At least one -datapath argument needs to be specified.");
    iResult=-EINVAL;
  }
  
  // -- Add files according to entry type to event
  if ( iResult >= 0 ) {
    iResult = InsertFiles();
  }

  // -- Calculate event size
  if ( fPublishESD )
    fMaxSize += sizeof(AliESDEvent);
  if ( fPublishHLTESD )
    fMaxSize += sizeof(AliESDEvent);
  if ( fPublishMC )
    fMaxSize += sizeof(AliMCEvent) + 20000;
  
  if ( fEvents.GetSize() == 0) {
    HLTError("No Files have been added.");
    iResult=-EINVAL;
  }

  if ( iResult < 0 ) {
    fEvents.Clear();
  }

  if (!fSkippedEsdObjects.IsNull()) {
    fpESDClone=new AliESDEvent;
  }
  
  return iResult;
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::DoDeinit() {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( fpESD ) 
    delete fpESD;
  fpESD = NULL;

  if ( fpHLTESD ) 
    delete fpHLTESD;
  fpHLTESD = NULL;

  if ( fpESDClone ) 
    delete fpESDClone;
  fpESDClone = NULL;

  if ( fpMC ) 
    delete fpMC;
  fpMC = NULL;

  if ( fpHLTMC ) 
    delete fpHLTMC;
  fpHLTMC = NULL;


  CloseCurrentFileList();
  
  AliHLTFilePublisher::DoDeinit();

  return iResult;
}

// #################################################################################
void AliHLTESDMCEventPublisherComponent::AddDataTypesToOutputlist() {
  // see header file for class documentation

  // add all different data types to the list
  if ( fPublishESD  )
    fOutputDataTypes.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline);
  else if ( fPublishHLTESD )
    fOutputDataTypes.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT);
  else if( fPublishMC )
    fOutputDataTypes.push_back(kAliHLTDataTypeMCObject);
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::InsertFiles(){
  // see header file for class documentation

  Int_t iResult=0;

  EventFiles* pCurrentFolder = NULL;

  TObjString *objPath = 0;
  TIter next(&fFolderList);

  // Loop over all folders
  while ( ( objPath = (TObjString*) next() ) && iResult>=0 ) {

    TString pathAliESDs = objPath->GetString();      
    pathAliESDs.Append("/AliESDs.root");

    TString pathKinematics = objPath->GetString();      
    pathKinematics.Append("/Kinematics.root");

    TString pathgalice = objPath->GetString();      
    pathgalice.Append("/galice.root");

    TString pathTrackRefs = objPath->GetString();      
    pathTrackRefs.Append("/TrackRefs.root");


    // Add ESD file for this folder
    if ( fPublishESD || fPublishHLTESD ) {
      FileDesc* pDesc=new FileDesc( pathAliESDs.Data(), kAliHLTDataTypeESDObject, kAliHLTVoidDataSpec, kFALSE);
      if ( pDesc ) 
	iResult = InsertFile( pCurrentFolder, pDesc );
      else 
	iResult = -ENOMEM;
    }

    // Add MC files for this folder
    if ( fPublishMC ) {
      FileDesc* pDescKinematics=new FileDesc( pathKinematics.Data(), kAliHLTDataTypeMCObject, kAliHLTVoidDataSpec, kFALSE);
      if ( pDescKinematics )  
	iResult = InsertFile( pCurrentFolder, pDescKinematics );
      else 
	iResult = -ENOMEM;

      FileDesc* pDescgalice=new FileDesc( pathgalice.Data(), kAliHLTDataTypeMCObject, kAliHLTVoidDataSpec, kFALSE);
      if ( pDescgalice )  
	iResult = InsertFile( pCurrentFolder, pDescgalice );
      else 
	iResult = -ENOMEM;

      if ( ! fFastMC ) {
	FileDesc* pDescTrackRefs=new FileDesc( pathTrackRefs.Data(), kAliHLTDataTypeMCObject, kAliHLTVoidDataSpec, kFALSE);
	if ( pDescTrackRefs )  
	  iResult = InsertFile( pCurrentFolder, pDescTrackRefs );
	else 
	  iResult = -ENOMEM;
      }
    }
    
    // -- Add all files belonging to one eventfolder
    InsertEvent( pCurrentFolder );
    
  } // while ( ( path = (TObjString*)next() ) ) {
  
  return iResult;
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::GetEvent( const AliHLTComponentEventData& /*evtData*/,
						AliHLTComponentTriggerData& /*trigData*/,
						AliHLTUInt8_t* /*outputPtr*/, 
						AliHLTUInt32_t& size,
						vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation

  if ( !IsDataEvent() ) return 0;

  Int_t iResult=0;
  size=0;

  // -- Initialize folder list and 
  //    re-initialize after reaching the end in the event before
  //    
  if ( !fpCurrentFolder ) {
    fpCurrentFolder = GetEventList()->FirstLink();

    // -- Set file list of this folder and open the files
    if ( fpCurrentFolder ) {
      iResult = OpenCurrentFileList();
    }
    else {
      HLTError("Folder link can not be retrieved from folder list.");
      iResult=-ENOENT;
    }
  }
  
  // -- if file list was properly created
  if ( iResult >= 0 && fpCurrentFileList ) {
    
    // -- Publish ESDs
    if ( fPublishESD  && fpTreeESD ) {
      fpTreeESD->GetEntry(fCurrentEvent);
      AliESDEvent* pESDpublish=fpESD;
      if (fpESDClone && CopyESDObjects(fpESDClone, fpESD, fSkippedEsdObjects.Data())>=0) {
	pESDpublish=fpESDClone;
      }
      if ((iResult=PushBack( pESDpublish, kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline , fSpecification  ))==-ENOSPC) {
	fOutputSize+=GetLastObjectSize();
      }
    }

    // -- Publish HLTESDs
    if ( fPublishHLTESD && fpTreeHLTESD ) {
      fpTreeHLTESD->GetEntry(fCurrentEvent);
      AliESDEvent* pESDpublish=fpHLTESD;
      if (fpESDClone && CopyESDObjects(fpESDClone, fpHLTESD, fSkippedEsdObjects.Data())>=0) {
	pESDpublish=fpESDClone;
      }
      if ((iResult=PushBack(pESDpublish, kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT , fSpecification ))==-ENOSPC) {
	fOutputSize+=GetLastObjectSize();	  
      }
    }

    // -- Publish MC
    if ( fPublishMC ) {

      TObjLink *flnk = fpCurrentFileList->FirstLink();
      
      // -- Loop over files in the current list
      while (flnk && iResult>=0) {
	
	FileDesc* pFileDesc = dynamic_cast<FileDesc*>( flnk->GetObject() );
	if (!pFileDesc) {
	  ALIHLTERRORGUARD(1, "internal mismatch, object is not of type FileDesc");
	  break;
	}
	TFile* pFile = *pFileDesc;
	if ( !pFile ) {
	  HLTError("No pointer to file");
	  iResult=-EFAULT;
	  continue;
	}

	TString filename = pFile->GetName();
	filename.Remove(0,filename.Last('/')+1);
	
	TString foldername( Form("Event%d", fCurrentEvent) );

	// -- Open Kinematics
	if ( !filename.CompareTo("Kinematics.root") ){

	  if ( fNEventsInFolder != (UInt_t)(pFile->GetListOfKeys()->GetEntries()) ) {
	    HLTError("Not all files contain the same number of events : %u and %d", 
		     fNEventsInFolder, pFile->GetListOfKeys()->GetEntries() );
	    iResult=-EFAULT;
	  }
	  
	  // -- Get Tree K
	  else {
	    TDirectoryFile *dirK = dynamic_cast<TDirectoryFile*> (pFile->Get(foldername));
	    if ( dirK )
	      fpTreeK = dynamic_cast<TTree*> (dirK->Get("TreeK"));
	    else {
	      HLTError("Directory %s could not be found.",  foldername.Data() );
	      iResult=-EFAULT;
	    }
	  }
	} // if ( !filename.CompareTo("Kinematics.root") ){

	// -- Open TrackReferences
	else if ( !fFastMC && !filename.CompareTo("TrackRefs.root") ) {
	  
	  if ( fNEventsInFolder != (UInt_t)(pFile->GetListOfKeys()->GetEntries()) ) {
	    HLTError("Not all files contain the same number of events : %u and %d", 
		     fNEventsInFolder, pFile->GetListOfKeys()->GetEntries() );
	    iResult=-EFAULT;
	  }
	  
	  // -- Get Tree TR
	  else {
	    TDirectoryFile *dirTR = dynamic_cast<TDirectoryFile*> (pFile->Get(foldername));
	    if ( dirTR )
	      fpTreeTR = dynamic_cast<TTree*> (dirTR->Get("TreeTR"));
	    else {
	      HLTError("Directory %s could not be found.",  foldername.Data() );
	      iResult=-EFAULT;
	    }
	  }
	} // else if ( !filename.CompareTo("TrackRefs.root") ){
	
	flnk = flnk->Next();
	
      } // while (flnk && iResult>=0) {
      
      if ( fpMC ) delete fpMC;
      fpMC = new AliMCEvent();
      
      if ( fpTreeE && fpTreeK && fpTreeTR ) {
	fpMC->ConnectTreeE(fpTreeE);
	fpTreeE->GetEntry(fCurrentEvent);

	fpMC->ConnectTreeK(fpTreeK);
	fpMC->ConnectTreeTR(fpTreeTR);	
      }
      else if( fpTreeE && fpTreeK && fFastMC ) {
	fpMC->ConnectTreeE(fpTreeE);
	fpTreeE->GetEntry(fCurrentEvent);

	fpMC->ConnectTreeK(fpTreeK);
      }
      else {
	HLTError("Not all trees needed for MC Event have been created.");
	iResult=-EFAULT;
      }

      // -- Fill AliHLTMCEvent 
      if ( iResult>=0 && fpMC ) {
	if ((iResult=PushBack( fpMC, kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline , fSpecification ))==-ENOSPC) {
	  fOutputSize+=GetLastObjectSize();
	}

	fpHLTMC = new AliHLTMCEvent( fApplyParticleCuts );
	
	if ( iResult>=0 && !(iResult=fpHLTMC->FillMCEvent(fpMC)) )
	  if ((iResult=PushBack( fpHLTMC, kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT , fSpecification ))==-ENOSPC) {
	    fOutputSize+=GetLastObjectSize();	    
	  }
      }

    } // if ( fPublishMC ) {
    
    // -- Next Event in Folder, 
    //    if last event was already reached:
    //    -> Reser and close files
    if (iResult>=0 && ++fCurrentEvent >= fNEventsInFolder )
      iResult = CloseCurrentFileList();
    
  } // if ( iResult >= 0 && fCurrentFileList ) {

  // -- Next Folder 
  //   (fpCurrentEvent will be NULL*pFileDesc if end of events in file is reached)
  if ( iResult >= 0 && fpCurrentFolder && !fCurrentEvent ) {
    fpCurrentFolder=fpCurrentFolder->Next();

    // -- if next Folder present : Set current file list of next folder
    //    else : end of folder list reached, reset fpCurrentFileList
    if ( fpCurrentFolder )
      iResult = OpenCurrentFileList();

  }

  return iResult;
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::OpenCurrentFileList() {
  // see header file for class documentation

  Int_t iResult = 0;
  Int_t nEntries = -1;

  // -- Get Current Folder
  EventFiles* pCurrentFolderDesc = dynamic_cast<EventFiles*>( fpCurrentFolder->GetObject() );

  // -- Get List of files in current folder  
  if ( pCurrentFolderDesc )
    fpCurrentFileList = *pCurrentFolderDesc; // type conversion operator defined
  else {
    HLTError("Folder descriptor can not be retrieved from folder list link");
    iResult=-EFAULT;
  }
  
  if ( fpCurrentFileList && iResult >= 0 ) {
  
    TObjLink *flnk = fpCurrentFileList->FirstLink();
  
    // -- Loop over all files in folder
    // ----------------------------------
    while ( flnk && iResult>=0 ) {
      FileDesc* pFileDesc = dynamic_cast<FileDesc*>( flnk->GetObject() );
      if (!pFileDesc) {
	ALIHLTERRORGUARD(1, "internal mismatch, object is not of type FileDesc");
	break;
      }
      TFile* pFile = *pFileDesc;
    
      // If file is not opened yet, open it ... which it should not
      if (pFile != NULL ) {
	HLTError("File already open, this should not happen.");
	iResult=-EFAULT;
	continue;
      }
      
      pFileDesc->OpenFile();
      if ( (pFile = *pFileDesc ) == NULL ) {
	HLTError("Opening file failed." );
	iResult=-EFAULT;
	continue;
      }

      // -- Open trees
      // ---------------
      TString filename = pFile->GetName();
      filename.Remove(0,filename.Last('/')+1);
      
      // -- Open ESD tree
      if ( !filename.CompareTo("AliESDs.root") && fPublishESD ) {
	fpTreeESD = dynamic_cast<TTree*>(pFile->Get("esdTree") );

	if ( fpTreeESD ) {
	  if ( nEntries == -1 ) 
	    nEntries = (Int_t)(fpTreeESD->GetEntriesFast());
	  else if ( nEntries != (Int_t)(fpTreeESD->GetEntriesFast()) ) {
	    HLTError("Not all files contain the same number of events : %d and %ld", 
		     nEntries, fpTreeESD->GetEntriesFast());
	    iResult=-EFAULT;
	  }
	
	  // Setup AliESDEvent
	  if ( fpESD ) delete fpESD;
	  fpESD = new AliESDEvent();	    
	  fpESD->ReadFromTree(fpTreeESD); 

	} // if ( fpTreeESD ) {
	else {
	  HLTError("Getting ESD tree failed for file %s.", pFile->GetName() );
	  iResult=-EFAULT;
	}
      } // if ( !filename.CompareTo("AliESDs.root") && fPublishESD ){


      // -- Open HLTESD tree
      if ( !filename.CompareTo("AliESDs.root") && fPublishHLTESD ) {
	fpTreeHLTESD = dynamic_cast<TTree*>(pFile->Get("HLTesdTree") );

	if ( fpTreeHLTESD ) {
	  if ( nEntries == -1 ) 
	    nEntries = (Int_t)(fpTreeHLTESD->GetEntriesFast());
	  else if ( nEntries != (Int_t)(fpTreeHLTESD->GetEntriesFast()) ) {
	    HLTError("Not all files contain the same number of events : %d and %ld", 
		     nEntries, fpTreeHLTESD->GetEntriesFast());
	    iResult=-EFAULT;
	  }

	  // Setup AliESDEvent
	  if ( fpHLTESD ) delete fpHLTESD;
	  fpHLTESD = new AliESDEvent();	    
	  fpHLTESD->ReadFromTree(fpTreeHLTESD); 

	} // if ( fpTreeHLTESD ) {
	else {
	  HLTError("Getting HLTESD tree failed for file %s.", pFile->GetName() );
	  iResult=-EFAULT;
	}

      } // if ( !filename.CompareTo("AliESDs.root") && fPublishHLTESD ){

      // -- Open galice tree : TreeE
      else if ( !filename.CompareTo("galice.root") && fPublishMC ) {
	fpTreeE = dynamic_cast<TTree*>(pFile->Get("TE") );

	if ( fpTreeE ) {
	  if ( nEntries == -1 ) 
	    nEntries = (Int_t)(fpTreeE->GetEntriesFast());
	  else if ( nEntries != (Int_t)(fpTreeE->GetEntriesFast()) ) {
	    HLTError("Not all files contain the same number of events : %d and %ld", 
		     nEntries, fpTreeE->GetEntriesFast());
	    iResult=-EFAULT;
	  }
	} // if ( fpTreeE ) {
	else {
	  HLTError("Getting TreeE failed for file %s.", pFile->GetName() );
	  iResult=-EFAULT;
	}
      } // if ( !filename.CompareTo("galice.root") && fPublishMC ){

      flnk = flnk->Next();
      
    } // while ( flnk && iResult>=0 ) {

  } // if ( fpCurrentFileList && iResult >= 0 ) {
  else {
    HLTError("Current file list could not get created.");
    iResult=-EFAULT;
  }
  
  if ( iResult >= 0 )
    fNEventsInFolder = (UInt_t)(nEntries);

  return iResult;  
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::CloseCurrentFileList() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Reset 
  fCurrentEvent = 0;
  fNEventsInFolder = 0;

  // if not set, no files are open ...
  if ( fpCurrentFileList ) {
    
    TObjLink *flnk = fpCurrentFileList->FirstLink();
  
    // -- Loop over all files in folder
    while ( flnk && iResult>=0 ) {
      FileDesc* pFileDesc = dynamic_cast<FileDesc*>( flnk->GetObject() );
      if (!pFileDesc) {
	ALIHLTERRORGUARD(1, "internal mismatch, object is not of type FileDesc");
	break;
      }
      
      TFile* pFile = *pFileDesc;
    
      // If is opened, close it
      if ( pFile != NULL ) {
      	pFileDesc->CloseFile();
      }
      
      if ( (pFile = *pFileDesc ) != NULL ) {
	HLTError("Closing file %s failed.", pFile->GetName() );
	iResult=-EFAULT;
      }

      flnk = flnk->Next();
      
    } // while ( flnk && iResult>=0 ) {

    fpCurrentFileList = NULL;

  } // if ( fpCurrentFileList && iResult >= 0 ) {
  
  return iResult;  
}

// #################################################################################
Int_t AliHLTESDMCEventPublisherComponent::CopyESDObjects(AliESDEvent* pTgt,
							 const AliESDEvent* pSrc,
							 const char* skippedObjects) const {
  // clone an ESD by copying all objects but skip the specified ones

  Int_t iResult = 0;
  if (!pSrc || !pTgt) return -EINVAL;

  // copy the full ESD
  pTgt->Reset();
  *pTgt=*pSrc;

  // filter according to the list of objects to be skipped
  if (pTgt->GetList() && skippedObjects!=NULL) {
    TString skipObjects=skippedObjects;
    TObjArray* pTokens=skipObjects.Tokenize(" ");
    if (pTokens) {
      const char* id=NULL;
      TIter next(pTokens);
      TObject* pObject=NULL;
      while ((pObject=next())!=NULL) {
	id=pObject->GetName();
	TObject* pObj=pTgt->GetList()->FindObject(id);
	if (pObj) {
	  HLTDebug("removing object %s", id);
	  pTgt->GetList()->Remove(pObj);
	  delete pObj;
	} else {
	  HLTWarning("failed to remove object '%s' from ESD", id);
	}
      }
      pTgt->GetStdContent();
      delete pTokens;
    }
  }

  return iResult;
}
