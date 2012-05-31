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

//  Description:
//  This class is responsible for creating folder structure
//  All stuff of aliroot sits in one folder with name defined by
//  fgkTopFolderName data wich do not very trough event to event are
//  sitting in directly in "top folder" all data which changes from
//  event to event are sitting in one folder (which has more subfolders)
//  Idea is to have more than one event in folder structure which allows
//  usage of standard procedures in merging
//  Add(AliDetector*) calls Add(AliModule*) as AliDetector is a AliModule
//  as well and should be listed in module list

#include <TDatabasePDG.h>
#include <TFolder.h>
#include <TInterpreter.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TVirtualMC.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliGenerator.h" 
#include "AliLoader.h"
#include "AliLog.h"

enum
 {
   kDetFolderData = 0,
   kDetFolderCalibration,
   kDetFolderAligmnet,
   kDetFolderLast
 };
ClassImp(AliConfig)

AliConfig* AliConfig::fgInstance = 0;

//0 level folder
const TString AliConfig::fgkTopFolderName("Folders");

//1st level folder
const TString AliConfig::fgkConstantsFolderName("Constants");
const TString AliConfig::fgkDefaultEventFolderName("Event");  //default folder for event, always used except merging

//2st level folder
//subfolder of event folder
const TString AliConfig::fgkDataFolderName("Data");//folder for data (hits, digits, points, tracks) grouped by detectors
const TString AliConfig::fgkModuleFolderName("Modules");//folder with modules objects
const TString AliConfig::fgkConditionsFolderName("Conditions");//folder with conditions (mag. field etc.)
const TString AliConfig::fgkConfigurationFolderName("Configuration");//folder with configuration (setup) of the detector
const TString AliConfig::fgkHeaderFolderName("Header");//folder with header and other MC information

//3rd level folder
//fgkConditionsFolderName subfolders
const TString AliConfig::fgkCalibrationFolderName("Calibration");
const TString AliConfig::fgkAligmentFolderName("Aligment");
  
//3rd level folder
//fgkConfigurationFolderName subfolders
const TString AliConfig::fgkFieldFolderName("Field");
const TString AliConfig::fgkGeneratorsFolderName("Generators");
const TString AliConfig::fgkVirtualMCFolderName("VirtualMC");


const TString AliConfig::fgkPDGFolderName("Constants/DatabasePDG");//folder with PDG Database
const TString AliConfig::fgkGeneratorFolderName("Configuration/Generators");//folder with generators
const TString AliConfig::fgkMCFolderName("Configuration/VirtualMC");

//____________________________________________________________________________
AliConfig* AliConfig::Instance ()
{
  //
  // Instance method for singleton class
  //
   if(fgInstance == 0) 
    {
     fgInstance = new AliConfig (fgkTopFolderName,"Alice data exchange board");
    }
   return fgInstance;
}
//____________________________________________________________________________

AliConfig::AliConfig(const char *name, const char *title): 
  TNamed(name,title), 
  fTopFolder(gROOT->GetRootFolder()->AddFolder(name,title)),
  fConstFolder(0x0),
  fDetectorFolder(new TString[kDetFolderLast+1])
{
// Constructor

  //Main AliRoot Folder
  if (fTopFolder == 0x0)
   {
     AliFatal("Can not create Top Alice Folder.");
     return;//never reached
   }
  fTopFolder->SetOwner();
  
  fDetectorFolder[kDetFolderData] = fgkDataFolderName;
  fDetectorFolder[kDetFolderCalibration] = fgkConditionsFolderName+"/"+fgkCalibrationFolderName;
  fDetectorFolder[kDetFolderAligmnet] = fgkConditionsFolderName+"/"+fgkAligmentFolderName;
  fDetectorFolder[kDetFolderLast] = "";
  
  gROOT->GetListOfBrowsables()->Add(fTopFolder, name);

  //Constants folder
  fConstFolder = fTopFolder->AddFolder (fgkConstantsFolderName, "Constant parameters");
  fConstFolder->AddFolder("DatabasePDG", "PDG database");
  
  fgInstance=this;
}

//____________________________________________________________________________
AliConfig::~AliConfig()
{ 
  // destructor
  delete [] fDetectorFolder ;  
  if (fTopFolder)
   {
    fTopFolder->SetOwner();
    delete fTopFolder; 
   }
}
//____________________________________________________________________________

void AliConfig::AddInFolder (const char *dir, TObject *obj)
{
  // Adds object "obj" to folder "dir"
  TFolder *folder = dynamic_cast<TFolder *>(fTopFolder->FindObject(dir));
  if (folder)
    folder->Add (static_cast<TObject *>(obj));
}

//____________________________________________________________________________
TObject* AliConfig::FindInFolder (const char *dir, const char *name)
{
  // Finds object with name "name" in folder "dir"
  if(!name) return(fTopFolder->FindObject(name));
  TFolder * folder = dynamic_cast<TFolder *>(fTopFolder->FindObject(dir));
  if (!folder) return (NULL);
  return(folder->FindObject(name));
}

//____________________________________________________________________________
void    AliConfig::Add (AliGenerator * obj,const char* eventfolder)
{
  // Adds generator "obj" to the event folder "eventfolder"
  TString path(eventfolder);
  path = path + "/" + fgkGeneratorsFolderName;
  AddInFolder(path,obj);
}

//____________________________________________________________________________
void AliConfig::Add (TVirtualMC * obj,const char* eventfolder)
{
  // Adds TVirtualMC object to the event folder
  TString path(eventfolder);
  path = path + "/" + fgkMCFolderName;
  AddInFolder(path, obj);
}

//____________________________________________________________________________
void  AliConfig::Add (TDatabasePDG * obj)
{
  // Adds TDataBase object
  AddInFolder(fgkPDGFolderName, obj);
}

//____________________________________________________________________________
void AliConfig::Add(AliModule* obj,const char* eventfolder)
{
  // Adds module to the event folder
  TString path(eventfolder);
  path = path + "/" + fgkModuleFolderName;
  AliDebug(1, Form("module name = %s, Ev. Fold. Name is %s.",
		   obj->GetName(),eventfolder));
  AddInFolder(path, obj);
}
//____________________________________________________________________________

Int_t AliConfig::AddDetector(TFolder* evntfolder, const char *name, const char* title)
{
//creates folders for the detector 'name'
 Int_t retval;//returned value
 retval = CreateDetectorFolders(evntfolder,name,title);
 if (retval)
  {
    AliError(Form("CreateDetectorFolders returned error for detector %s",name));
    return retval;
  }
 return 0; 
}
//____________________________________________________________________________

Int_t AliConfig::AddDetector(const char* evntfoldername,const char *name, const char* title)
{
//creates folders for the detector 'name'
 Int_t retval;//returned value
 retval = CreateDetectorFolders(evntfoldername,name,title);
 if (retval)
  {
    AliError(Form("CreateDetectorFolders returned error for detector %s",name));
    return retval;
  }
 return 0; 
}
//____________________________________________________________________________

void  AliConfig::Add(AliDetector * obj,const char* eventfolder)
{
  // Adds new AliDetector objest to the correspondent event folder
  AliDebug(1, Form("detector name = %s, Ev. Fold. Name is %s.",
		   obj->GetName(),eventfolder));

  TObject* foundobj = GetTopFolder()->FindObject(eventfolder);
  TFolder* evfolder = (foundobj)?dynamic_cast<TFolder*>(foundobj):0x0;
  if (evfolder == 0x0)
   {
     AliFatal(Form("Can not find folder %s while adding detector %s",eventfolder,obj->GetName()));
     return;
   } 
  CreateDetectorFolders(evfolder, obj->GetName(), obj->GetTitle());
  
}
//____________________________________________________________________________

Int_t  AliConfig::CreateDetectorFolders(const char* evntfoldername,const char *name, const char* title)
{
//creates a folders for detector named 'name' and titled 'title'
//in a event folder named 'evntfoldername'
//list of folder names where new folders are created is defined in fDetectorFolder array 
//detector folders are named 'name' and titled 'title' as well

 TFolder* evfolder = dynamic_cast<TFolder*>(GetTopFolder()->FindObject(evntfoldername));
 if (evfolder == 0x0)
  {
   AliError(Form("Can not find folder %s while adding detector %s",evntfoldername,name));
   return 1;
  }
 return CreateDetectorFolders(evfolder,name,title);
}
//____________________________________________________________________________
Int_t  AliConfig::CreateDetectorFolders(TFolder* evntfolder,const char *name, const char* title)
{
//creates a folders for detector named 'name' and titled 'title'
//in a event folder 'evntfolder'
//list of folder names where new folders are created is defined in fDetectorFolder array 
//detector folders are named 'name' and titled 'title' as well
//Here we add only detector not an modules
 
 Int_t tmp;
 Int_t i = 0;//iterator
 while(!fDetectorFolder[i].IsNull())
  {
    tmp = AddSubFolder(evntfolder,fDetectorFolder[i],name,title);
    if (tmp)
     {
      AliError(Form("Failed to create subfolder of %s for detector %s",fDetectorFolder[i].Data(),name));
      return 1;
     }
    i++;
  }
 return 0;
}

/*****************************************************************************/

TFolder* AliConfig::BuildEventFolder(const char* name,const char* title)
{
/*
 creates the folder structure for one event
 TopFolder
         |_
         | \
         |  Constants
         |_
         | \
         |  Event_
         |      | \
         |      |  Modules(detector objects)
         |      |_
         |      | \              
         |      |  Header
         |      |_
         |      | \              
         |      |  Data_
         |      |     | \ 
         |      |     |  TPC_
         |      |     |    | \
         |      |     |    |  Hits(object;e.g. tree)
         |      |     |    |_  
         |      |     |    | \ 
         |      |     |    |  SDigits(object)
         |      |     |    |_
         |      |     |    | \ 
         |      |     |    |  Digits(object)
         |      |     |    |_
         |      |     |    | \ 
         |      |     |    |  RecPoints(object)
         |      |     |    |_
         |      |     |      \ 
         |      |     |       Tracks(object)
         |      |     |_ 
         |      |       \
         |      |        ITS_
         |      |          | \
         |      |          |  Hits(object;e.g. tree)
         |      |          |_  
         |      |          | \ 
         |      |          |  SDigits(object)
         |      |          |_
         |      |          | \ 
         |      |          |  Digits(object)
         |      |          |_
         |      |          | \ 
         |      |          |  RecPoints(object)
         |      |          |_
         |      |            \ 
         |      |             Tracks(object)
         |      |_         
         |        \       
         |         Configuration
         |               
         |_
           \
            Event2_  (to be merged with event)
                |  \
                |   Modules(detector objects)
                |_
                | \              
                |  Header
                |_
                | \              
                |  Data_
                |     | \ 
                |     |  TPC_
                |     |    | \
                |     |    |  Hits(object;e.g. tree)
                |     |    |_  
                |     |    | \ 
                |     |    |  SDigits(object)
                |     |    |_
                |     |    | \ 
                |     |    |  Digits(object)
                |     |    |_
                |     |    | \ 
                |     |    |  RecPoints(object)
                |     |    |_
                |     |      \ 
                |     |       Tracks(object)
                |     |_ 
                |       \
                |        ITS_
                |          | \
                |          |  Hits(object;e.g. tree)
                |          |_  
                |          | \ 
                |          |  SDigits(object)
                |          |_
                |          | \ 
                |          |  Digits(object)
                |          |_
                |          | \ 
                |          |  RecPoints(object)
                |          |_
                |            \ 
                |             Tracks(object)
                |_         
                  \       
                   Configuration
                         
*/
  TFolder* eventfolder = fTopFolder->AddFolder(name,title); 
   
  //modules
  eventfolder->AddFolder(fgkModuleFolderName, "Detector objects");
  //event data
  eventfolder->AddFolder(fgkDataFolderName, "Detector data");
  
    //Conditions
  TFolder *conditions = eventfolder->AddFolder(fgkConditionsFolderName, "Run conditions");
  conditions->AddFolder(fgkCalibrationFolderName,"Detector calibration data");
  conditions->AddFolder(fgkAligmentFolderName,"Detector aligment");
  //Configuration
  TFolder *configuration = eventfolder->AddFolder(fgkConfigurationFolderName, "Run configuration");
  configuration->AddFolder(fgkFieldFolderName, "Magnetic field maps");
  configuration->AddFolder(fgkGeneratorsFolderName,"list of generator objects");
  //PH configuration->AddFolder(fgkVirtualMCFolderName,"the Virtual MC");

  eventfolder->AddFolder(fgkHeaderFolderName,"MonteCarlo event header");

  eventfolder->SetOwner();

  return eventfolder;
}

/*****************************************************************************/

const TString& AliConfig::GetDataFolderName() const
{
//returns name of data folder path relative to event folder
 return fgkDataFolderName;
}

/*****************************************************************************/

Int_t AliConfig::AddSubFolder(TFolder* topfolder, const char* infoler, 
                     const char* newfoldname, const char* newfoldtitle)
{
//helper method
//in topfolder looks for folder named 'infolder'
//and if it exist it creates inside a new folder named 'newfoldname' and titled 'newfoldtitle'

 if (topfolder == 0x0)//check if exists top folder
  {
   AliError("Parameter TFolder* points to NULL.");
   return 1;
  }
 
 TObject *obj;
 TFolder* folder;
 
 folder = dynamic_cast<TFolder*>(topfolder->FindObject(infoler));
 if (folder == 0x0) //check if we got inolder
  {
   AliError(Form("Can not find folder %s in folder %s.",infoler,topfolder->GetName()));
   return 1;
  }
 obj = folder->FindObject(newfoldname); //see if such a subfolder already exists
 if (obj == 0x0) //nope
  {
   TFolder *newfolder = folder->AddFolder(newfoldname,newfoldtitle);//add the desired subfolder
   if (newfolder == 0x0) //check if we managed
    {
     AliError(Form("Can not create folder %s in folder %s",newfoldname,infoler));
     return 2;
    }
   else return 0;//success
  }
 //such an object already exists
 TFolder* fol = dynamic_cast<TFolder*>(obj);
 if (fol == 0x0)
   {
     AliError(Form("Object named %s already exists in folder %s AND IT IS NOT A FOLDER",newfoldname,infoler));
     return 3;
   }
 AliWarning(Form("Folder named %s already exists in folder %s",newfoldname,infoler));
 return 0;
}
