/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
// class AliRunLoader                                               //
//                                                                  //
// This class aims to be the unque interface for managing data I/O. //
// It stores Loaders for all modules which, knows names             //
// of the files were data are to be stored.                         //
//                                                                  //
// It aims to substitud AliRun in automatic data managing           //
// thus there is no necessity of loading gAlice from file in order  //
// to get access to the data.                                       //
//                                                                  //
// Logical place to put the specific Loader to the given            //
// detector is detector  itself (i.e ITSLoader in ITS).             //
// But, to load detector object one need to load gAlice, and        //
// by the way all other detectors with their geometrieces and       //
// so on. So, if one need to open TPC clusters there is no          //
// principal need to read everything.                               //
//                                                                  //
//                                                                  //
// When RunLoader is read from the file it does not connect to      //
// the folder structure automatically. It must be connected         //
// (mounted) manualy. Default event folder is defined by            //
// AliConfig::GetDefaultEventFolderName()                           //
// but can be mounted elsewhere. Usefull specially in merging case, //
// when more than pone session needs to be loaded                   //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TBranch.h>
#include <TFile.h>
#include <TFolder.h>
#include <TGeometry.h>
#include <TObjArray.h>
#include <TString.h>
class TTask;
#include <TTree.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliLoader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliDetector.h"
#include "AliCDBManager.h"
#include "AliCDBLocal.h"
#include "AliCentralTrigger.h"

ClassImp(AliRunLoader)

AliRunLoader* AliRunLoader::fgRunLoader = 0x0;

const TString AliRunLoader::fgkRunLoaderName("RunLoader");
const TString AliRunLoader::fgkHeaderBranchName("Header");
const TString AliRunLoader::fgkTriggerBranchName("ClassMask");
const TString AliRunLoader::fgkHeaderContainerName("TE");
const TString AliRunLoader::fgkTriggerContainerName("TreeCT");
const TString AliRunLoader::fgkKineContainerName("TreeK");
const TString AliRunLoader::fgkTrackRefsContainerName("TreeTR");
const TString AliRunLoader::fgkKineBranchName("Particles");
const TString AliRunLoader::fgkDefaultKineFileName("Kinematics.root");
const TString AliRunLoader::fgkDefaultTrackRefsFileName("TrackRefs.root");
const TString AliRunLoader::fgkGAliceName("gAlice");
const TString AliRunLoader::fgkDefaultTriggerFileName("Trigger.root");
/**************************************************************************/

AliRunLoader::AliRunLoader():
 fLoaders(0x0),
 fEventFolder(0x0),
 fCurrentEvent(0),
 fGAFile(0x0),
 fHeader(0x0),
 fStack(0x0),
 fCTrigger(0x0),
 fKineDataLoader(0x0),
 fTrackRefsDataLoader(0x0),
 fNEventsPerFile(1),
 fUnixDirName(".")
{
  AliConfig::Instance();//force to build the folder structure
  if (!fgRunLoader) fgRunLoader = this;
}
/**************************************************************************/

AliRunLoader::AliRunLoader(const char* eventfoldername):
 TNamed(fgkRunLoaderName,fgkRunLoaderName),
 fLoaders(new TObjArray()),
 fEventFolder(0x0),
 fCurrentEvent(0),
 fGAFile(0x0),
 fHeader(0x0),
 fStack(0x0),
 fCTrigger(0x0),
 fKineDataLoader(new AliDataLoader(fgkDefaultKineFileName,fgkKineContainerName,"Kinematics")),
 fTrackRefsDataLoader(new AliDataLoader(fgkDefaultTrackRefsFileName,fgkTrackRefsContainerName,"Track References")),
 fNEventsPerFile(1),
 fUnixDirName(".")
{
//ctor
  SetEventFolderName(eventfoldername);
 if (!fgRunLoader) fgRunLoader = this;
}
/**************************************************************************/

AliRunLoader::AliRunLoader(const AliRunLoader &rl):
 TNamed(rl),
 fLoaders(0x0),
 fEventFolder(0x0),
 fCurrentEvent(0),
 fGAFile(0x0),
 fHeader(0x0),
 fStack(0x0),
 fCTrigger(0x0),
 fKineDataLoader(0x0),
 fTrackRefsDataLoader(0x0),
 fNEventsPerFile(0),
 fUnixDirName(".")
{
  //
  // Copy ctor
  //
  rl.Copy(*this);
}
/**************************************************************************/

AliRunLoader::~AliRunLoader()
{
//dtor
  if (fgRunLoader == this) fgRunLoader = 0x0;
  
  UnloadHeader();
  UnloadgAlice();
  
  if(fLoaders) {
    fLoaders->SetOwner();
    delete fLoaders;
  }
  
  delete fKineDataLoader;
  delete fTrackRefsDataLoader;
  
  
  RemoveEventFolder();
  
  //fEventFolder is deleted by the way of removing - TopAliceFolder owns it
  if( fCTrigger ) delete  fCTrigger;
  delete fHeader;
  delete fStack;
  delete fGAFile;
}
/**************************************************************************/

AliRunLoader::AliRunLoader(TFolder* topfolder):
 TNamed(fgkRunLoaderName,fgkRunLoaderName),
 fLoaders(new TObjArray()),
 fEventFolder(topfolder),
 fCurrentEvent(0),
 fGAFile(0x0),
 fHeader(0x0),
 fStack(0x0),
 fCTrigger(0x0),
 fKineDataLoader(new AliDataLoader(fgkDefaultKineFileName,fgkKineContainerName,"Kinematics")),
 fTrackRefsDataLoader(new AliDataLoader(fgkDefaultTrackRefsFileName,fgkTrackRefsContainerName,"Track References")),
 fNEventsPerFile(1),
 fUnixDirName(".")
{
//ctor
 if(topfolder == 0x0)
  {
    TString errmsg("Parameter is NULL");
    AliError(errmsg.Data());
    throw errmsg;
    return;
  }
 
 TObject* obj = fEventFolder->FindObject(fgkRunLoaderName);
 if (obj)
  { //if it is, then sth. is going wrong... exits aliroot session
    TString errmsg("In Event Folder Named ");
    errmsg+=fEventFolder->GetName();
    errmsg+=" object named "+fgkRunLoaderName+" already exists. I am confused ...";

    AliError(errmsg.Data());
    throw errmsg;
    return;//never reached
  }

 if (!fgRunLoader) fgRunLoader = this;
   
 fEventFolder->Add(this);//put myself to the folder to accessible for all
  
}
/**************************************************************************/

void AliRunLoader::Copy(TObject &) const 
{
  AliFatal("Not implemented");
}
/**************************************************************************/

Int_t AliRunLoader::GetEvent(Int_t evno)
{
//Gets event number evno
//Reloads all data properly
//PH  if (fCurrentEvent == evno) return 0;
  
  if (evno < 0)
   {
     AliError("Can not give the event with negative number");
     return 4;
   }

  if (evno >= GetNumberOfEvents())
   {
     AliError(Form("There is no event with number %d",evno));
     return 3;
   }
  
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
  AliDebug(1, Form("          GETTING EVENT  %d",evno));
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
  AliDebug(1, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
   
  fCurrentEvent = evno;

  Int_t retval;
  
  //Reload header (If header was loaded)
  if (GetHeader())
   {
     retval = TreeE()->GetEvent(fCurrentEvent);
     if ( retval == 0)
      {
        AliError(Form("Cannot find event: %d\n ",fCurrentEvent));
        return 5;
      }
   }
  //Reload stack (If header was loaded)
  if (TreeE()) fStack = GetHeader()->Stack();
  //Set event folder in stack (it does not mean that we read kinematics from file)
  if (fStack) 
   { 
     fStack->SetEventFolderName(fEventFolder->GetName());
   }
  else
   {
     AliWarning("Stack not found in header");
   }

   if( GetTrigger() && TreeCT() ) {
      retval = TreeCT()->GetEvent(fCurrentEvent);
      if ( retval < 0 )      {
         AliError(Form("Error occured while GetEvent for Trigger. Event %d",evno));
         return 2;
      }
   }
  
  retval = SetEvent();
  if (retval)
   {
     AliError(Form("Error occured while setting event %d",evno));
     return 1;
   }
   
  //Post Track References
  retval = fTrackRefsDataLoader->GetEvent();
  if (retval)
   {
     AliError(Form("Error occured while GetEvent for Track References. Event %d",evno));
     return 2;
   }

  //Read Kinematics if loaded
  fKineDataLoader->GetEvent();
  if (retval)
   {
     AliError(Form("Error occured while GetEvent for Kinematics. Event %d",evno));
     return 2;
   }

  if (fStack && fKineDataLoader->GetBaseLoader(0)->IsLoaded())
    {
      if (fStack->GetEvent() == kFALSE)
	{
	  AliError(Form("Error occured while GetEvent for Stack. Event %d",evno));
	  return 2;
	}
    }

  //Trigger data reloading in all loaders 
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
     retval = loader->GetEvent();
     if (retval)
      {
       AliError(Form("Error occured while getting event for %s. Event %d.",
		     loader->GetDetectorName().Data(), evno));
       return 3;
      }
   }
  
  SetDetectorAddresses();
  
  return 0;
}
/**************************************************************************/
Int_t AliRunLoader::SetEvent()
{
//if kinematics was loaded Cleans folder data

  Int_t retval;
  
  retval = fKineDataLoader->SetEvent();
  if (retval)
   {
     AliError("SetEvent for Kinamtics Data Loader retutned error.");
     return retval;
   }
  retval = fTrackRefsDataLoader->SetEvent(); 
  if (retval)
   {
     AliError("SetEvent for Track References Data Loader retutned error.");
     return retval;
   }

  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
     retval = loader->SetEvent();
     if (retval)
      {
        AliError(Form("SetEvent for %s Data Loader retutned error.",loader->GetName()));
        return retval;
      }
   }

  return 0;
}
/**************************************************************************/

Int_t AliRunLoader::SetEventNumber(Int_t evno)
{
  //cleans folders and sets the root dirs in files 
  if (fCurrentEvent == evno) return 0;
  fCurrentEvent = evno;
  return SetEvent();
}

/**************************************************************************/
AliCDBEntry* AliRunLoader::GetCDBEntry(const char* name) const
{
//Get an AliCDBEntry from the run data storage

  if ( !(AliCDBManager::Instance()->IsDefaultStorageSet()) ) {
    AliError("No run data storage defined!");
    return 0x0;
  }
  return AliCDBManager::Instance()->GetDefaultStorage()->Get(name, GetHeader()->GetRun());

}

/**************************************************************************/
AliRunLoader* AliRunLoader::Open
  (const char* filename, const char* eventfoldername, Option_t* option)
{
//Opens a desired file 'filename'
//gets the the run-Loader and mounts it desired folder
//returns the pointer to run Loader which can be further used for accessing data 
//in case of error returns NULL
 
 static const TString kwebaddress("http://alisoft.cern.ch/people/skowron/codedoc/split/index.html");
 AliDebugClass(1,Form("\n\n\nNew I/O strcture: See more info:\n %s\n\n\n",kwebaddress.Data()));
 
 AliRunLoader* result = 0x0;
 
 /* ************************************************ */
 /* Chceck if folder with given name already exists  */
 /* ************************************************ */
 
 TObject* obj = AliConfig::Instance()->GetTopFolder()->FindObject(eventfoldername);
 if(obj)
  {
    TFolder* fold = dynamic_cast<TFolder*>(obj);
    if (fold == 0x0)
     {
      AliErrorClass("Such a obejct already exists in top alice folder and it is not a folder.");
      return 0x0;
     }
    
    //check if we can get RL from that folder
    result = AliRunLoader::GetRunLoader(eventfoldername);
    if (result == 0x0)
     {
       AliErrorClass(Form("Folder %s already exists, and can not find session there. Can not mount.",eventfoldername));
       return 0x0;
     }

    if (result->GetFileName().CompareTo(filename) != 0)
     {
       AliErrorClass("Other file is mounted in demanded folder. Can not mount.");
       return 0x0;
     }

    //check if now is demanded (re)creation 
    if ( AliLoader::TestFileOption(option) == kFALSE)
     {
       AliErrorClass(Form("Session already exists in folder %s and this session option is %s. Unable to proceed.",
			  eventfoldername,option));
       return 0x0;
     }
     
    //check if demanded option is update and existing one 
    TString tmpstr(option);
    if ( (tmpstr.CompareTo("update",TString::kIgnoreCase) == 0) && 
         (result->fGAFile->IsWritable() == kFALSE) )
     { 
       AliErrorClass(Form("Session already exists in folder %s and is not writable while this session option is %s. Unable to proceed.",
			  eventfoldername,option));
       return 0x0;
     }
     
    AliWarningClass("Session is already opened and mounted in demanded folder");	
    if (!fgRunLoader) fgRunLoader = result; //PH get access from any place
    return result;
  } //end of checking in case of existance of object named identically that folder session is being opened
 
 
 TFile * gAliceFile = TFile::Open(filename,option);//open a file
 if (!gAliceFile) 
  {//null pointer returned
    AliFatalClass(Form("Can not open file %s.",filename));
    return 0x0;
  }
  
 if (gAliceFile->IsOpen() == kFALSE)
  {//pointer to valid object returned but file is not opened
    AliErrorClass(Form("Can not open file %s.",filename));
    return 0x0;
  }
 
 //if file is "read" or "update" than we try to find AliRunLoader there - if not found cry and exit
 //else create new AliRunLoader
 if ( AliLoader::TestFileOption(option) )
  { 
    AliDebugClass(1, "Reading RL from file");
    
    result = dynamic_cast<AliRunLoader*>(gAliceFile->Get(fgkRunLoaderName));//get the run Loader from the file
    if (result == 0x0)
     {//didn't get
       AliErrorClass(Form("Can not find run-Loader in file %s.",filename));
       delete gAliceFile;//close the file
       return 0x0;
     }
    Int_t tmp = result->SetEventFolderName(eventfoldername);//mount a event folder   
    if (tmp)//if SetEvent  returned error
     {
       AliErrorClass(Form("Can not mount event in folder %s.",eventfoldername));
       delete result; //delete run-Loader
       delete gAliceFile;//close the file
       return 0x0;
     }
  }
 else
  {
    AliDebugClass(1, Form("Creating new AliRunLoader. Folder name is %s",eventfoldername));
    try
     {  
       result = new AliRunLoader(eventfoldername);
     }
    catch (TString& errmsg)
     {
       AliErrorClass(Form("AliRunLoader constrcutor has thrown exception: %s\n",errmsg.Data()));
       delete result;
       delete gAliceFile;//close the file
       return 0x0;
     }
  }
 
//procedure for extracting dir name from the file name 
 TString fname(filename);
 Int_t  nsl = fname.Last('#');//look for hash in file name
 TString dirname;
 if (nsl < 0) {//hash not found
   nsl = fname.Last('/');// look for slash
   if (nsl < 0) 
     nsl = fname.Last(':');// look for colon e.g. rfio:galice.root
 }

 if (nsl < 0) dirname = "./";      // no directory path, use "."
 else dirname = fname.Remove(nsl+1);// directory path
 
 AliDebugClass(1, Form("Dir name is : %s",dirname.Data()));
 
 result->SetDirName(dirname); 
 result->SetGAliceFile(gAliceFile);//set the pointer to gAliceFile
 if (!fgRunLoader) fgRunLoader = result; //PH get access from any place
 return result;
}
/**************************************************************************/
Int_t AliRunLoader::GetNumberOfEvents()
{
 //returns number of events in Run
 Int_t retval;
 if( TreeE() == 0x0 )
  {
    retval = LoadHeader();
    if (retval) 
     {
       AliError("Error occured while loading header");
       return -1;
     }
  }
 return (Int_t)TreeE()->GetEntries();
}
/**************************************************************************/

void AliRunLoader::MakeHeader()
{
 //Makes header and connects it to header tree (if it exists)
  AliDebug(1, "");
  if(fHeader == 0x0)
   {
     AliDebug(1, "Creating new Header Object");
     fHeader= new AliHeader();
   }
  TTree* tree = TreeE();
  if (tree)
   {
     AliDebug(1, "Got Tree from folder.");
     TBranch* branch = tree->GetBranch(fgkHeaderBranchName);
     if (branch == 0x0)
      {
        AliDebug(1, "Creating new branch");
        branch = tree->Branch(fgkHeaderBranchName, "AliHeader", &fHeader, 4000, 0);
        branch->SetAutoDelete(kFALSE);
      }
     else
      {
        AliDebug(1, "Got Branch from Tree");
        branch->SetAddress(&fHeader);
        tree->GetEvent(fCurrentEvent);
        fStack = fHeader->Stack(); //should be safe - if we created Stack, header returns pointer to the same object
        if (fStack)
         {
           fStack->SetEventFolderName(fEventFolder->GetName());
           if (TreeK()) fStack->GetEvent();
         }
        else
        {
          AliDebug(1, "Header does not have a stack.");
        }
      }
   } 
  AliDebug(1, "Exiting MakeHeader method");
}
/**************************************************************************/

void AliRunLoader::MakeStack()
{
//Creates the stack object -  do not connect the tree
  if(fStack == 0x0)
   { 
     fStack = new AliStack(10000);
     fStack->SetEventFolderName(fEventFolder->GetName());
   }
}
/**************************************************************************/

void AliRunLoader::MakeTrigger()
{
 // Makes trigger object and connects it to trigger tree (if it exists)
   AliDebug( 1, "" );
   if( fCTrigger == 0x0 ) {
      AliDebug( 1, "Creating new Trigger Object" );
      fCTrigger = new AliCentralTrigger();
   }
   TTree* tree = TreeCT();
   if( tree ) {
      fCTrigger->MakeBranch( fgkTriggerBranchName, tree );
      tree->GetEvent( fCurrentEvent );
   }

   AliDebug( 1, "Exiting MakeTrigger method" );
}
/**************************************************************************/

void AliRunLoader::MakeTree(Option_t *option)
{
//Creates trees
  const char *oK  = strstr(option,"K");  //Kine
  const char *oE  = strstr(option,"E");  //Header
  const char *oGG = strstr(option,"GG"); //Central TriGGer

  if(oK && !TreeK())
   { 
     if (fKineDataLoader->GetBaseLoader(0)->IsLoaded() == kFALSE)
      {
        AliError("Load Kinematics first");
      }
     else
      {
        fKineDataLoader->MakeTree();
        MakeStack();
        fStack->ConnectTree();
        WriteKinematics("OVERWRITE");
     }
   }
  
  if(oE && !TreeE())
   { 
     fGAFile->cd();
     TTree* tree = new TTree(fgkHeaderContainerName,"Tree with Headers");
     GetEventFolder()->Add(tree);
     MakeHeader();
     WriteHeader("OVERWRITE");
   }
  
   if(oGG && !TreeCT())
   {
      // create the CTP Trigger output file and tree
      TFile* file = gROOT->GetFile( fgkDefaultTriggerFileName );
      if( !file ) {
         file = TFile::Open( gSystem->ConcatFileName( fUnixDirName.Data(), fgkDefaultTriggerFileName.Data() ), "RECREATE" ) ;
      }

      file->cd();
      TTree* tree = new TTree( fgkTriggerContainerName, "Tree with Central Trigger Mask" );
      GetEventFolder()->Add(tree);
      MakeTrigger();
  //    WriteHeader("OVERWRITE");
   }

  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
    loader->MakeTree(option);
   }

}
/**************************************************************************/
    
Int_t AliRunLoader::LoadgAlice()
{
//Loads gAlice from file
 if (GetAliRun())
  {
    AliWarning("AliRun is already in folder. Unload first.");
    return 0;
  }
 AliRun* alirun = dynamic_cast<AliRun*>(fGAFile->Get(fgkGAliceName));
 if (alirun == 0x0)
  {
    AliError(Form("Can not find gAlice in file %s",fGAFile->GetName()));
    return 2;
  }
 alirun->SetRunLoader(this);
 if (gAlice)
  {
    AliWarning(Form("gAlice already exists. Putting retrived object in folder named %s",
		    GetEventFolder()->GetName()));
  }
 else
  {
    gAlice = alirun;
  }
 SetDetectorAddresses();//calls SetTreeAddress for all detectors
 return 0; 
}
/**************************************************************************/

Int_t AliRunLoader::LoadHeader()
{
//loads treeE and reads header object for current event
 if (TreeE())
  {
     AliWarning("Header is already loaded. Use ReloadHeader to force reload. Nothing done");
     return 0;
  }
 
 if (GetEventFolder() == 0x0)
  {
    AliError("Event folder not specified yet");
    return 1;
  }

 if (fGAFile == 0x0)
  {
    AliError("Session not opened. Use AliRunLoader::Open");
    return 2;
  }
 
 if (fGAFile->IsOpen() == kFALSE)
  {
    AliError("Session not opened. Use AliRunLoader::Open");
    return 2;
  }

 TTree* tree = dynamic_cast<TTree*>(fGAFile->Get(fgkHeaderContainerName));
 if (tree == 0x0)
  {
    AliError(Form("Can not find header tree named %s in file %s",
		  fgkHeaderContainerName.Data(),fGAFile->GetName()));
    return 2;
  }

 if (tree == TreeE()) return 0;

 CleanHeader();
 GetEventFolder()->Add(tree);
 MakeHeader();//creates header object and connects to tree
 return 0; 

}
/**************************************************************************/

Int_t AliRunLoader::LoadTrigger(Option_t* option)
{
   //Load treeCT

   if( TreeCT() ) {
      AliWarning("Trigger is already loaded. Nothing done");
      return 0;
   }
 
   if( GetEventFolder() == 0x0 ) {
      AliError("Event folder not specified yet");
      return 1;
   }
   // get the CTP Trigger output file and tree
   TString trgfile = gSystem->ConcatFileName( fUnixDirName.Data(),
                                              fgkDefaultTriggerFileName.Data() );
   TFile* file = gROOT->GetFile( trgfile );
   if( !file ) {
      file = TFile::Open( trgfile, option ) ;
      if (!file || file->IsOpen() == kFALSE ) {
         AliError( Form( "Can not open trigger file %s", trgfile.Data() ) );
         return 2;
      }
   }
   file->cd();

   TTree* tree = dynamic_cast<TTree*>(file->Get( fgkTriggerContainerName ));
   if( !tree ) {
      AliError( Form( "Can not find trigger tree named %s in file %s",
                      fgkTriggerContainerName.Data(), file->GetName() ) );
      return 2;
   }

   CleanTrigger();

   fCTrigger = dynamic_cast<AliCentralTrigger*>(file->Get( "AliCentralTrigger" ));
   GetEventFolder()->Add( tree );
   MakeTrigger();

   return 0;
}

/**************************************************************************/

Int_t AliRunLoader::LoadKinematics(Option_t* option)
{
//Loads the kinematics 
 Int_t retval = fKineDataLoader->GetBaseLoader(0)->Load(option);
 if (retval)
  {
    AliError("Error occured while loading kinamatics tree.");
    return retval;
  }
 if (fStack) 
  {
    retval = fStack->GetEvent();
    if ( retval == kFALSE)
     {
       AliError("Error occured while loading kinamatics tree.");
       return retval;
     }
    
  }
 return 0;
}
/**************************************************************************/

Int_t AliRunLoader::OpenDataFile(const TString& filename,TFile*& file,TDirectory*& dir,Option_t* opt,Int_t cl)
{
//Opens File with kinematics
 if (file)
  {
    if (file->IsOpen() == kFALSE)
     {//pointer is not null but file is not opened
       AliWarning("Pointer to file is not null, but file is not opened");//risky any way
       delete file;
       file = 0x0; //proceed with opening procedure
     }
    else
     { 
       AliWarning(Form("File  %s already opened",filename.Data()));
       return 0;
     }
  }
//try to find if that file is opened somewere else
 file = (TFile *)( gROOT->GetListOfFiles()->FindObject(filename) );
 if (file)
  {
   if(file->IsOpen() == kTRUE)
    {
     AliWarning(Form("File %s already opened by sombody else.",file->GetName()));
     return 0;
    }
  }

 file = TFile::Open(filename,opt);
 if (file == 0x0)
  {//file is null
    AliError(Form("Can not open file %s",filename.Data()));
    return 1;
  }
 if (file->IsOpen() == kFALSE)
  {//file is not opened
    AliError(Form("Can not open file %s",filename.Data()));
   return 1;
  }
  
 file->SetCompressionLevel(cl);
 
 dir = AliLoader::ChangeDir(file,fCurrentEvent);
 if (dir == 0x0)
  {
    AliError(Form("Can not change to root directory in file %s",filename.Data()));
    return 3;
  }
 return 0; 
}
/**************************************************************************/

TTree* AliRunLoader::TreeE() const
{
 //returns the tree from folder; shortcut method
 if (AliDebugLevel() > 10) fEventFolder->ls();
 TObject *obj = fEventFolder->FindObject(fgkHeaderContainerName);
 return (obj)?dynamic_cast<TTree*>(obj):0x0;
}
/**************************************************************************/

TTree* AliRunLoader::TreeCT() const
{
 //returns the tree from folder; shortcut method
   if (AliDebugLevel() > 10) fEventFolder->ls();
   TObject *obj = fEventFolder->FindObject(fgkTriggerContainerName);
   return (obj)?dynamic_cast<TTree*>(obj):0x0;
}
/**************************************************************************/

AliHeader* AliRunLoader::GetHeader() const
{
//returns pointer header object
 return fHeader;
}
/**************************************************************************/

AliCentralTrigger* AliRunLoader::GetTrigger() const
{
//returns pointer trigger object
   return fCTrigger;
}

/**************************************************************************/
 
TTree* AliRunLoader::TreeK() const
{
 //returns the tree from folder; shortcut method
 TObject *obj = GetEventFolder()->FindObject(fgkKineContainerName);
 return (obj)?dynamic_cast<TTree*>(obj):0x0;
}
/**************************************************************************/

TTree* AliRunLoader::TreeTR() const
{
 //returns the tree from folder; shortcut method
 TObject* obj = GetEventFolder()->FindObject(fgkTrackRefsContainerName);
 return (obj)?dynamic_cast<TTree*>(obj):0x0;
}
/**************************************************************************/

AliRun* AliRunLoader::GetAliRun() const
{
//returns AliRun which sits in the folder
 if (fEventFolder == 0x0) return 0x0;
 TObject *obj = fEventFolder->FindObject(fgkGAliceName);
 return (obj)?dynamic_cast<AliRun*>(obj):0x0;
}
/**************************************************************************/

Int_t AliRunLoader::WriteGeometry(Option_t* /*opt*/)
{
//writes geometry to the file
  fGAFile->cd();
  TGeometry* geo = GetAliRun()->GetGeometry();
  if (geo == 0x0)
   {
     AliError("Can not get geometry from gAlice");
     return 1;
   }
  geo->Write();
  return 0;
}
/**************************************************************************/

Int_t AliRunLoader::WriteHeader(Option_t* opt)
{
//writes treeE
  AliDebug(1, "WRITING HEADER");
  
  TTree* tree = TreeE();
  if ( tree == 0x0)
   {
     AliWarning("Can not find Header Tree in Folder");
     return 0;
   } 
  if (fGAFile->IsWritable() == kFALSE)
   {
     AliError(Form("File %s is not writable",fGAFile->GetName()));
     return 1;
   }

  TObject* obj = fGAFile->Get(fgkHeaderContainerName);
  if (obj)
   { //if they exist, see if option OVERWRITE is used
     TString tmp(opt);
     if(tmp.Contains("OVERWRITE",TString::kIgnoreCase) == 0)
      {//if it is not used -  give an error message and return an error code
        AliError("Tree already exisists. Use option \"OVERWRITE\" to overwrite previous data");
        return 3;
      }
   }
  fGAFile->cd();
  tree->SetDirectory(fGAFile);
  tree->Write(0,TObject::kOverwrite);

  AliDebug(1, "WRITTEN\n\n");
  
  return 0;
}

/**************************************************************************/

Int_t AliRunLoader::WriteTrigger(Option_t* opt)
{
   //writes TreeCT
   AliDebug( 1, "WRITING TRIGGER" );
  
   TTree* tree = TreeCT();
   if ( tree == 0x0) {
      AliWarning("Can not find Trigger Tree in Folder");
      return 0;
   }

   TFile* file = gROOT->GetFile( gSystem->ConcatFileName( fUnixDirName.Data(), fgkDefaultTriggerFileName.Data() ) ) ;
   if( !file || !file->IsOpen() ) {
      AliError( "can't write Trigger, file is not open" );
      return kFALSE;
   }

   TObject* obj = file->Get( fgkTriggerContainerName );
   if( obj ) { //if they exist, see if option OVERWRITE is used
      TString tmp(opt);
      if( tmp.Contains( "OVERWRITE", TString::kIgnoreCase ) == 0) {
         //if it is not used -  give an error message and return an error code
         AliError( "Tree already exisists. Use option \"OVERWRITE\" to overwrite previous data" );
         return 3;
      }
   }
   file->cd();
   fCTrigger->Write( 0, TObject::kOverwrite );
   tree->Write( 0, TObject::kOverwrite );
   file->Flush();

   AliDebug(1, "WRITTEN\n\n");
  
   return 0;
}
/**************************************************************************/

Int_t AliRunLoader::WriteAliRun(Option_t* /*opt*/)
{
//writes AliRun object to the file
  fGAFile->cd();
  if (GetAliRun()) GetAliRun()->Write();
  return 0;
}
/**************************************************************************/

Int_t AliRunLoader::WriteKinematics(Option_t* opt)
{
//writes Kinematics
  return fKineDataLoader->GetBaseLoader(0)->WriteData(opt);
}
/**************************************************************************/
Int_t AliRunLoader::WriteTrackRefs(Option_t* opt)
{
//writes Track References tree
  return fTrackRefsDataLoader->GetBaseLoader(0)->WriteData(opt);
}
/**************************************************************************/

Int_t AliRunLoader::WriteHits(Option_t* opt)
{
//Calls WriteHits for all loaders
  Int_t res;
  Int_t result = 0;
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     res = loader->WriteHits(opt);
     if (res)
      {
        AliError(Form("Failed to write hits for %s (%d)",loader->GetDetectorName().Data(),res));
        result = 1;
      }
   }
  return result;
}
/**************************************************************************/

Int_t AliRunLoader::WriteSDigits(Option_t* opt)
{
//Calls WriteSDigits for all loaders
  Int_t res;
  Int_t result = 0;
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     res = loader->WriteSDigits(opt);
     if (res)
      {
        AliError(Form("Failed to write summable digits for %s.",loader->GetDetectorName().Data()));
        result = 1;
      }
   }
  return result;
}
/**************************************************************************/

Int_t AliRunLoader::WriteDigits(Option_t* opt)
{
//Calls WriteDigits for all loaders
  Int_t res;
  Int_t result = 0;
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   { 
     res = loader->WriteDigits(opt);
     if (res)
      {
        AliError(Form("Failed to write digits for %s.",loader->GetDetectorName().Data()));
        result = 1;
      }
   }
  return result;
}
/**************************************************************************/

Int_t AliRunLoader::WriteRecPoints(Option_t* opt)
{
//Calls WriteRecPoints for all loaders
  Int_t res;
  Int_t result = 0;
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     res = loader->WriteRecPoints(opt);
     if (res)
      {
        AliError(Form("Failed to write Reconstructed Points for %s.",
		      loader->GetDetectorName().Data()));
        result = 1;
      }
   }
  return result;
}
/**************************************************************************/

Int_t AliRunLoader::WriteTracks(Option_t* opt)
{
//Calls WriteTracks for all loaders
  Int_t res;
  Int_t result = 0;
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     res = loader->WriteTracks(opt);
     if (res)
      {
        AliError(Form("Failed to write Tracks for %s.",
		      loader->GetDetectorName().Data()));
        result = 1;
      }
   }
  return result;
}
/**************************************************************************/

Int_t AliRunLoader::WriteRunLoader(Option_t* /*opt*/)
{
//Writes itself to the file
  CdGAFile();
  this->Write(0,TObject::kOverwrite);
  return 0;
}
/**************************************************************************/

Int_t AliRunLoader::SetEventFolderName(const TString& name)
{  
//sets top folder name for this run; of alread
  if (name.IsNull())
   {
     AliError("Name is empty");
     return 1;
   }
  
  //check if such a folder already exists - try to find it in alice top folder
  TObject* obj = AliConfig::Instance()->GetTopFolder()->FindObject(name);
  if(obj)
   {
     TFolder* fold = dynamic_cast<TFolder*>(obj);
     if (fold == 0x0)
      {
       AliError("Such a obejct already exists in top alice folder and it is not a folder.");
       return 2;
      }
     //folder which was found is our folder
     if (fEventFolder == fold)
      {
       return 0;
      }
     else
      {
       AliError("Such a folder already exists in top alice folder. Can not mount.");
       return 2;
      }
   }

  //event is alredy connected, just change name of the folder
  if (fEventFolder) 
   {
     fEventFolder->SetName(name);
     return 0;
   }

  if (fKineDataLoader == 0x0)
    fKineDataLoader = new AliDataLoader(fgkDefaultKineFileName,fgkKineContainerName,"Kinematics");

  if ( fTrackRefsDataLoader == 0x0)
    fTrackRefsDataLoader = new AliDataLoader(fgkDefaultTrackRefsFileName,fgkTrackRefsContainerName,"Track References");
   
  //build the event folder structure
  AliDebug(1, Form("Creating new event folder named %s",name.Data()));
  fEventFolder = AliConfig::Instance()->BuildEventFolder(name,"Event Folder");
  fEventFolder->Add(this);//put myself to the folder to accessible for all
  
  if (Stack()) Stack()->SetEventFolderName(fEventFolder->GetName());
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     loader->Register(fEventFolder);//build folder structure for this detector
   }
  
  fKineDataLoader->SetEventFolder(GetEventFolder());
  fTrackRefsDataLoader->SetEventFolder(GetEventFolder());
  fKineDataLoader->SetFolder(GetEventFolder());
  fTrackRefsDataLoader->SetFolder(GetEventFolder());
  
  fEventFolder->SetOwner();
  return 0;
}
/**************************************************************************/

void AliRunLoader::AddLoader(AliLoader* loader)
 {
 //Adds the Loader for given detector 
  if (loader == 0x0) //if null shout and exit
   {
     AliError("Parameter is NULL");
     return;
   }
  loader->SetDirName(fUnixDirName);
  if (fEventFolder) loader->SetEventFolder(fEventFolder); //if event folder is already defined, 
                                                          //pass information to the Loader
  fLoaders->Add(loader);//add the Loader to the array
 }
/**************************************************************************/

void AliRunLoader::AddLoader(AliDetector* det)
 {
//Asks module (detector) ro make a Loader and stores in the array
   if (det == 0x0) return;
   AliLoader* get = det->GetLoader();//try to get loader
   if (get == 0x0)  get = det->MakeLoader(fEventFolder->GetName());//if did not obtain, ask to make it

   if (get) 
    {
      AliDebug(1, Form("Detector: %s   Loader : %s",det->GetName(),get->GetName()));
      AddLoader(get);
    }
 }

/**************************************************************************/

AliLoader* AliRunLoader::GetLoader(const char* detname) const
{
//returns loader for given detector
//note that naming convention is TPCLoader not just TPC
  return (AliLoader*)fLoaders->FindObject(detname);
}

/**************************************************************************/

AliLoader* AliRunLoader::GetLoader(AliDetector* det) const
{
//get loader for detector det
 if(det == 0x0) return 0x0;
 TString getname(det->GetName());
 getname+="Loader";
 AliDebug(1, Form(" Loader name is %s",getname.Data()));
 return GetLoader(getname);
}

/**************************************************************************/

void AliRunLoader::CleanFolders()
{
//  fEventFolder->Add(this);//put myself to the folder to accessible for all

  CleanDetectors();
  CleanHeader();
  CleanKinematics();
  CleanTrigger();
}
/**************************************************************************/

void AliRunLoader::CleanDetectors()
{
//Calls CleanFolders for all detectors
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
     loader->CleanFolders();
   }
}
/**************************************************************************/

void AliRunLoader::RemoveEventFolder()
{
//remove all the tree of event 
//all the stuff changing EbE stays untached (PDGDB, tasks, etc.)

 if (fEventFolder == 0x0) return;
 fEventFolder->SetOwner(kFALSE);//don't we want to deleted while removing the folder that we are sitting in
 fEventFolder->Remove(this);//remove us drom folder
 
 AliConfig::Instance()->GetTopFolder()->SetOwner(); //brings ownership back for fEventFolder since it sits in top folder
 AliConfig::Instance()->GetTopFolder()->Remove(fEventFolder); //remove the event tree
 delete fEventFolder;
}
/**************************************************************************/

void AliRunLoader::SetGAliceFile(TFile* gafile)
{
//sets pointer to galice.root file
 fGAFile = gafile;
}

/**************************************************************************/

Int_t AliRunLoader::LoadHits(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  AliDebug(1, "Loading Hits");
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     AliDebug(1, "Option is All");
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer array
   }   

  AliDebug(1, Form("For detectors. Number of detectors chosen for loading %d",loaders->GetEntries()));
  
  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    AliDebug(1, Form("    Calling LoadHits(%s) for %s",opt,loader->GetName()));
    loader->LoadHits(opt);
   }
  AliDebug(1, "Done");
  return 0;
} 

/**************************************************************************/

Int_t AliRunLoader::LoadSDigits(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->LoadSDigits(opt);
   }
  return 0;
} 

/**************************************************************************/

Int_t AliRunLoader::LoadDigits(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->LoadDigits(opt);
   }
  return 0;
} 
/**************************************************************************/

Int_t AliRunLoader::LoadRecPoints(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->LoadRecPoints(opt);
   }
  return 0;
} 
/**************************************************************************/

Int_t AliRunLoader::LoadRecParticles(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->LoadRecParticles(opt);
   }
  return 0;
} 
/**************************************************************************/

Int_t AliRunLoader::LoadTracks(Option_t* detectors,Option_t* opt)
{
//LoadHits in selected detectors i.e. detectors="ITS TPC TRD" or "all"

  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->LoadTracks(opt);
   }
  return 0;
} 
/**************************************************************************/

void AliRunLoader::UnloadHits(Option_t* detectors)
{
  //unloads hits for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadHits();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadSDigits(Option_t* detectors)
{
  //unloads SDigits for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadSDigits();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadDigits(Option_t* detectors)
{
  //unloads Digits for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadDigits();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadRecPoints(Option_t* detectors)
{
  //unloads RecPoints for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadRecPoints();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadAll(Option_t* detectors)
{
  //calls UnloadAll for detectors names specified in parameter
  // option "all" passed can be passed
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadAll();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadTracks(Option_t* detectors)
{
  //unloads Tracks for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadTracks();
   }
}
/**************************************************************************/

void AliRunLoader::UnloadRecParticles(Option_t* detectors)
{
  //unloads Particles for detectors specified in parameter
  TObjArray* loaders;
  TObjArray arr;

  const char* oAll = strstr(detectors,"all");
  if (oAll)
   {
     loaders = fLoaders;
   }
  else
   {
     GetListOfDetectors(detectors,arr);//this method looks for all Loaders corresponding to names (many) specified in detectors option
     loaders = &arr;//get the pointer to array
   }   

  TIter next(loaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
    loader->UnloadRecParticles();
   }
}
/**************************************************************************/

AliRunLoader* AliRunLoader::GetRunLoader(const char* eventfoldername)
{
//returns RunLoader from folder named eventfoldername
  TFolder* evfold= dynamic_cast<TFolder*>(AliConfig::Instance()->GetTopFolder()->FindObject(eventfoldername));
  if (evfold == 0x0)
   {
     return 0x0;
   }
  AliRunLoader* runget = dynamic_cast<AliRunLoader*>(evfold->FindObject(AliRunLoader::fgkRunLoaderName));
  return runget;
  
}
/**************************************************************************/

AliLoader* AliRunLoader::GetDetectorLoader(const char* detname, const char* eventfoldername)
{
//get the loader of the detector with the given name from the global
//run loader object
  AliRunLoader* runLoader = GetRunLoader(eventfoldername);
  if (!runLoader) {
    AliErrorClass("No run loader found");
    return NULL;
  }
  return runLoader->GetDetectorLoader(detname);
}
/**************************************************************************/

AliLoader* AliRunLoader::GetDetectorLoader(const char* detname)
{
//get the loader of the detector with the given name from the global
//run loader object
  
  char loadername[256];
  sprintf(loadername, "%sLoader", detname);
  AliLoader* loader = GetLoader(loadername);
  if (!loader) {
    AliError(Form("No loader for %s found", detname));
    return NULL;
  }
  return loader;
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeH(const char* detname, Bool_t maketree, const char* eventfoldername)
{
//get the tree with hits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeH() && maketree) loader->MakeTree("H");
  return loader->TreeH();
}

/**************************************************************************/

TTree* AliRunLoader::GetTreeH(const char* detname, Bool_t maketree)
{
//get the tree with hits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeH() && maketree) loader->MakeTree("H");
  return loader->TreeH();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeS(const char* detname, Bool_t maketree,const char* eventfoldername)
{
//get the tree with summable digits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeS() && maketree) loader->MakeTree("S");
  return loader->TreeS();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeS(const char* detname, Bool_t maketree)
{
//get the tree with summable digits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeS() && maketree) loader->MakeTree("S");
  return loader->TreeS();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeD(const char* detname, Bool_t maketree,const char* eventfoldername)
{
//get the tree with digits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeD() && maketree) loader->MakeTree("D");
  return loader->TreeD();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeD(const char* detname, Bool_t maketree)
{
//get the tree with digits of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeD() && maketree) loader->MakeTree("D");
  return loader->TreeD();
}
/**************************************************************************/
TTree* AliRunLoader::GetTreeR(const char* detname, Bool_t maketree,const char* eventfoldername)
{
//get the tree with clusters of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeR() && maketree) loader->MakeTree("R");
  return loader->TreeR();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeR(const char* detname, Bool_t maketree)
{
//get the tree with clusters of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeR() && maketree) loader->MakeTree("R");
  return loader->TreeR();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeT(const char* detname, Bool_t maketree,const char* eventfoldername)
{
//get the tree with tracks of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeT() && maketree) loader->MakeTree("T");
  return loader->TreeT();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeT(const char* detname, Bool_t maketree)
{
//get the tree with tracks of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeT() && maketree) loader->MakeTree("T");
  return loader->TreeT();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeP(const char* detname, Bool_t maketree,const char* eventfoldername)
{
//get the tree with particles of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname,eventfoldername);
  if (!loader) return NULL;
  if (!loader->TreeP() && maketree) loader->MakeTree("P");
  return loader->TreeP();
}
/**************************************************************************/

TTree* AliRunLoader::GetTreeP(const char* detname, Bool_t maketree)
{
//get the tree with particles of the detector with the given name
//if maketree is true and the tree does not exist, the tree is created
  AliLoader* loader = GetDetectorLoader(detname);
  if (!loader) return NULL;
  if (!loader->TreeP() && maketree) loader->MakeTree("P");
  return loader->TreeP();
}

/**************************************************************************/

void AliRunLoader::CdGAFile()
{
//sets gDirectory to galice file
//work around 
  if(fGAFile) fGAFile->cd();
}
 
/**************************************************************************/

void AliRunLoader::GetListOfDetectors(const char * namelist,TObjArray& pointerarray) const
 {
//this method looks for all Loaders corresponding 
//to names (many) specified in namelist i.e. namelist ("ITS TPC TRD")
  
   char buff[10];
   char dets [200];
   strcpy(dets,namelist);//compiler cries when char* = const Option_t*;
   dets[strlen(dets)+1] = '\n';//set endl at the end of string 
   char* pdet = dets;
   Int_t tmp;
   for(;;)
    {
      tmp = sscanf(pdet,"%s",buff);//read the string from the input string pdet into buff
      if ( (buff[0] == 0) || (tmp == 0) ) break; //if not read
     
      pdet = strstr(pdet,buff) + strlen(buff);//move the input pointer about number of bytes (letters) read
      //I am aware that is a little complicated. I don't know the number of spaces between detector names
      //so I read the string, than I find where it starts (strstr) and move the pointer about length of a string
      //If there is a better way, please write me (Piotr.Skowronski@cern.ch)
      //construct the Loader name
      TString getname(buff);
      getname+="Loader";
      AliLoader* loader = GetLoader(getname);//get the Loader
      if (loader)
       {
        pointerarray.Add(loader);
       }
      else
       {
        AliError(Form("Can not find Loader for %s",buff));
       }
        
      buff[0] = 0;
    }
 }
/*****************************************************************************/ 

void AliRunLoader::Clean(const TString& name)
{
//removes object with given name from event folder and deletes it
  if (GetEventFolder() == 0x0) return;
  TObject* obj = GetEventFolder()->FindObject(name);
  if(obj)
   {
     AliDebug(1, Form("name=%s, cleaning %s.",GetName(),name.Data()));
     GetEventFolder()->Remove(obj);
     delete obj;
     obj = 0x0;
   }
}

/*****************************************************************************/ 

TTask* AliRunLoader::GetRunDigitizer()
{
//returns Run Digitizer from folder

 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetDigitizerTaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}
/*****************************************************************************/ 

TTask* AliRunLoader::GetRunSDigitizer()
{
//returns SDigitizer Task from folder

 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetSDigitizerTaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}
/*****************************************************************************/ 

TTask* AliRunLoader::GetRunReconstructioner()
{
//returns Reconstructioner Task from folder
 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetReconstructionerTaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}
/*****************************************************************************/ 

TTask* AliRunLoader::GetRunTracker()
{
//returns Tracker Task from folder
 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetTrackerTaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}
/*****************************************************************************/ 

TTask* AliRunLoader::GetRunPIDTask()
{
//returns Tracker Task from folder
 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetPIDTaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}
/*****************************************************************************/ 

TTask* AliRunLoader::GetRunQATask()
{
//returns Quality Assurance Task from folder
 TFolder* topf = AliConfig::Instance()->GetTaskFolder();
 if (topf == 0x0)
  {
    AliErrorClass("Can not get task folder from AliConfig");
    return 0x0;
  }
 TObject* obj = topf->FindObjectAny(AliConfig::Instance()->GetQATaskName());
 return (obj)?dynamic_cast<TTask*>(obj):0x0;
}

/*****************************************************************************/ 

void AliRunLoader::SetCompressionLevel(Int_t cl)
{
//Sets Compression Level in all files
 if (fGAFile) fGAFile->SetCompressionLevel(cl);
 SetKineComprLevel(cl);
 SetTrackRefsComprLevel(cl);
 TIter next(fLoaders);
 AliLoader *loader;
 while((loader = (AliLoader*)next()))
  {
   loader->SetCompressionLevel(cl);
  }
}
/**************************************************************************/

void AliRunLoader::SetKineComprLevel(Int_t cl)
{
//Sets comression level in Kine File
  fKineDataLoader->SetCompressionLevel(cl);
}
/**************************************************************************/

void AliRunLoader::SetTrackRefsComprLevel(Int_t cl)
{
//Sets comression level in Track Refences File
  fTrackRefsDataLoader->SetCompressionLevel(cl);
}
/**************************************************************************/

void AliRunLoader::UnloadHeader()
{
 //removes TreeE from folder and deletes it
 // as well as fHeader object
 CleanHeader();
 delete fHeader;
 fHeader = 0x0;
}
/**************************************************************************/

void AliRunLoader::UnloadTrigger()
{
 //removes TreeCT from folder and deletes it
 // as well as fHeader object
   CleanTrigger();
   delete fCTrigger;
   fCTrigger = 0x0;
}

/**************************************************************************/

void AliRunLoader::UnloadKinematics()
{
//Unloads Kinematics
 fKineDataLoader->GetBaseLoader(0)->Unload();
}
/**************************************************************************/

void AliRunLoader::UnloadTrackRefs()
{
//Unloads Track Refernces
 fTrackRefsDataLoader->GetBaseLoader(0)->Unload();
}
/**************************************************************************/

void AliRunLoader::UnloadgAlice()
{
//Unloads gAlice
 if (gAlice == GetAliRun())
  {
   AliDebug(1, "Set gAlice = 0x0");
   gAlice = 0x0;//if gAlice is the same that in folder (to be deleted by the way of folder)
  }
 AliRun* alirun = GetAliRun();
 if (GetEventFolder()) GetEventFolder()->Remove(alirun);
 delete alirun;
}
/**************************************************************************/

void  AliRunLoader::MakeTrackRefsContainer()
{
// Makes a tree for Track References
  fTrackRefsDataLoader->MakeTree();
}
/**************************************************************************/

Int_t AliRunLoader::LoadTrackRefs(Option_t* option)
{
//Load track references from file (opens file and posts tree to folder)

 return fTrackRefsDataLoader->GetBaseLoader(0)->Load(option);
}
/**************************************************************************/

void  AliRunLoader::SetDirName(TString& dirname)
{
//sets directory name 
  if (dirname.IsNull()) return;
  fUnixDirName = dirname;
  fKineDataLoader->SetDirName(dirname);
  fTrackRefsDataLoader->SetDirName(dirname);
  
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
    loader->SetDirName(dirname);
   }

}
/*****************************************************************************/ 

Int_t AliRunLoader::GetFileOffset() const
{
//returns the file number that is added to the file name for current event
  return Int_t(fCurrentEvent/fNEventsPerFile);
}

/*****************************************************************************/ 
const TString AliRunLoader::SetFileOffset(const TString& fname)
{
//adds the the number to the file name at proper place for current event
  Long_t offset = (Long_t)GetFileOffset();
  if (offset < 1) return fname;
  TString soffset;
  soffset += offset;//automatic conversion to string
  TString dotroot(".root");
  const TString& offfsetdotroot = offset + dotroot;
  TString out = fname;
  out = out.ReplaceAll(dotroot,offfsetdotroot);
  AliDebug(1, Form(" in=%s out=%s",fname.Data(),out.Data()));
  return out;
}
/*****************************************************************************/ 

void AliRunLoader::SetDigitsFileNameSuffix(const TString& suffix)
{
//adds the suffix before ".root", 
//e.g. TPC.Digits.root -> TPC.DigitsMerged.root
//made on Jiri Chudoba demand

  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next())) 
   {
     loader->SetDigitsFileNameSuffix(suffix);
   }
}
/*****************************************************************************/ 

TString AliRunLoader::GetFileName() const
{
//returns name of galice file
 TString result;
 if (fGAFile == 0x0) return result;
 result = fGAFile->GetName();
 return result;
}
/*****************************************************************************/ 

void AliRunLoader::SetDetectorAddresses()
{
 //calls SetTreeAddress for all detectors
  if (GetAliRun()==0x0) return;
  TIter next(GetAliRun()->Modules());
  AliModule* mod;
  while((mod = (AliModule*)next())) 
   {
     AliDetector* det = dynamic_cast<AliDetector*>(mod);
     if (det) det->SetTreeAddress();
   }
}
/*****************************************************************************/ 

void AliRunLoader::Synchronize()
{
  //synchrinizes all writtable files 
  TIter next(fLoaders);
  AliLoader *loader;
  while((loader = (AliLoader*)next()))
   {
     loader->Synchronize();
   }
  
  fKineDataLoader->Synchronize();
  fTrackRefsDataLoader->Synchronize();

  TFile* file = gROOT->GetFile( gSystem->ConcatFileName( fUnixDirName.Data(), fgkDefaultTriggerFileName.Data() ) ) ;
  if( file ) file->Flush();
  
  if (fGAFile) fGAFile->Flush();
}
/*****************************************************************************/ 
/*****************************************************************************/ 
