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
#include <TTask.h>
#include <TVirtualMC.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliGenerator.h" 
#include "AliLoader.h"
#include "AliLog.h"

enum 
 {
  kDetTaskQA = 0,
  kDetTaskSDigitizer,
  kDetTaskDigitizer,
  kDetTaskRecontructioner,
  kDetTaskTracker,
  kDetTaskPID,
  kDetTaskLast
 };

enum
 {
   kDetFolderData = 0,
   kDetFolderCalibration,
   kDetFolderAligmnet,
   kDetFolderQA,
   kDetFolderLast
 };
ClassImp(AliConfig)

AliConfig* AliConfig::fgInstance = 0;

//0 level folder
const TString AliConfig::fgkTopFolderName("Folders");

//1st level folder
const TString AliConfig::fgkTasksFolderName("Tasks"); //task folder, commn
const TString AliConfig::fgkConstantsFolderName("Constants");
const TString AliConfig::fgkDefaultEventFolderName("Event");  //default folder for event, always used except merging

//2st level folder
//subfolder of event folder
const TString AliConfig::fgkDataFolderName("Data");//folder for data (hits, digits, points, tracks) grouped by detectors
const TString AliConfig::fgkModuleFolderName("Modules");//folder with modules objects
const TString AliConfig::fgkConditionsFolderName("Conditions");//folder with conditions (mag. field etc.)
const TString AliConfig::fgkConfigurationFolderName("Configuration");//folder with configuration (setup) of the detector
const TString AliConfig::fgkHeaderFolderName("Header");//folder with header and other MC information

//Tasks names, goes into fgkTasksFolderName folder
const TString AliConfig::fgkDigitizerTaskName("Digitizer");
const TString AliConfig::fgkSDigitizerTaskName("SDigitizer");
const TString AliConfig::fgkReconstructionerTaskName("Reconstructioner");
const TString AliConfig::fgkTrackerTaskName("Tracker");
const TString AliConfig::fgkPIDTaskName("PIDTask");//;=) PIDer???
const TString AliConfig::fgkQATaskName("QAtask");

//3rd level folder
//fgkConditionsFolderName subfolders
const TString AliConfig::fgkCalibrationFolderName("Calibration");
const TString AliConfig::fgkAligmentFolderName("Aligment");
const TString AliConfig::fgkQAFolderName("QAout");
  
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
AliConfig::AliConfig():
  fTopFolder(0x0),
  fTaskFolder(0x0),
  fConstFolder(0x0),
  fDetectorTask(0x0),
  fDetectorFolder(0x0)
{
  //
  // Default constructor, mainly to keep coding conventions
  //
  fgInstance=0;//never mind, its going to exit in next step
  AliFatal("Constructor should not be called for a singleton");
}
//____________________________________________________________________________

AliConfig::AliConfig(const AliConfig& conf):
  TNamed(conf),
  fTopFolder(0x0),
  fTaskFolder(0x0),
  fConstFolder(0x0),
  fDetectorTask(0x0),
  fDetectorFolder(0x0)
{
  //
  // Copy constructor, mainly to keep coding conventions
  //
  fgInstance=0;
    
  AliFatal("Copy constructor should not be called for a singleton");
}
//____________________________________________________________________________

AliConfig& AliConfig::operator=(const AliConfig& /*conf*/)
{
  //
  // Assignment, mainly to keep coding conventions
  //
  fgInstance=0;
    
  AliFatal("Assignment operator should not be called for a singleton");
  return *this;
}
//____________________________________________________________________________

AliConfig::AliConfig(const char *name, const char *title): 
  TNamed(name,title), 
  fTopFolder(gROOT->GetRootFolder()->AddFolder(name,title)),
  fTaskFolder(fTopFolder->AddFolder(fgkTasksFolderName, "ALICE Tasks")),
  fConstFolder(0x0),
  fDetectorTask(0x0),
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
  fDetectorFolder[kDetFolderQA] = fgkConditionsFolderName+"/"+fgkQAFolderName;
  fDetectorFolder[kDetFolderLast] = "";
  
  gROOT->GetListOfBrowsables()->Add(fTopFolder, name);

  //Constants folder
  fConstFolder = fTopFolder->AddFolder (fgkConstantsFolderName, "Constant parameters");
  fConstFolder->AddFolder("DatabasePDG", "PDG database");
  
  // Add the tasks to //Folders
  
  TTask * qa = new TTask(fgkQATaskName, "Alice QA tasks");
  fTaskFolder->Add(qa); 
  TTask * sd = new TTask(fgkSDigitizerTaskName, "Alice SDigitizer") ;
  fTaskFolder->Add(sd); 
  TTask * di = new TTask(fgkDigitizerTaskName, "Alice Digitizer") ;
  fTaskFolder->Add(di); 
  TTask * re = new TTask(fgkReconstructionerTaskName, "Alice Reconstructioner") ;
  fTaskFolder->Add(re); 
  TTask * tr = new TTask(fgkTrackerTaskName,"Alice Tracker");
  fTaskFolder->Add(tr);
  TTask * pid = new TTask(fgkPIDTaskName,"Alice Particle Identification Task");
  fTaskFolder->Add(pid);
  fDetectorTask    =  new TString[kDetTaskLast+1];
  
  fDetectorTask[kDetTaskQA] = fgkQATaskName;
  fDetectorTask[kDetTaskSDigitizer] = fgkSDigitizerTaskName;
  fDetectorTask[kDetTaskDigitizer] =  fgkDigitizerTaskName;
  fDetectorTask[kDetTaskRecontructioner] = fgkReconstructionerTaskName;
  fDetectorTask[kDetTaskTracker] = fgkTrackerTaskName;
  fDetectorTask[kDetTaskPID] = fgkPIDTaskName;
  fDetectorTask[kDetTaskLast] = "";

  fgInstance=this;
}

//____________________________________________________________________________
AliConfig::~AliConfig()
{ 
  // destructor
  delete [] fDetectorFolder ;  
  delete [] fDetectorTask;
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

Int_t AliConfig::AddSubTask(const char *taskname, const char* name,const char* title)
{
//Create new task named 'name' and titled 'title' 
//as a subtask of the task named 'taskname'

   AliDebug(1, Form("Try to get folder named %s",taskname));
   TObject* obj = fTopFolder->FindObject(taskname);
   TTask * task = (obj)?dynamic_cast<TTask*>(obj):0x0;
   if (task)
     {
      AliDebug(1, "          Got");
      TTask * subtask = static_cast<TTask*>(task->GetListOfTasks()->FindObject(name));
      if (!subtask) 
        {
          subtask = new TTask(name,title);
          task->Add(subtask);
        }
      else
       {
         AliWarning(Form("Task named \"%s\" already exists in Task %s",name,taskname));
       }
     }
   else
     {
       AliError(Form("Can not find task %s to put a new task in.",taskname));
       return 1;
     }
  return 0;
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
//creates folders and tasks for the detector 'name'
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
//creates folders and tasks for the detector 'name'
 Int_t retval;//returned value
 retval = CreateDetectorFolders(evntfoldername,name,title);
 if (retval)
  {
    AliError(Form("CreateDetectorFolders returned error for detector %s",name));
    return retval;
  }
// retval = CreateDetectorTasks(name,title);
// if (retval)
//  {
//    Error("AddDetector","CreateDetectorTasks returned error for detector %s",name);
//    return retval;
//  }
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
  
//  CreateDetectorTasks(obj->GetName(),obj->GetTitle());

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
//____________________________________________________________________________
Int_t AliConfig::CreateDetectorTasks(const char *name, const char* title)
{
  // Creates new detector's task "name"
   Int_t i = 0;
   Int_t tmp;
   while (i < kDetTaskLast)
    {
      tmp = AddSubTask(fgkTasksFolderName+"/"+fDetectorTask[i],
                       name+fDetectorTask[i],(fDetectorTask[i]+" for ")+title);
      if (tmp)
       {
         AliError(Form("Error occured while creating task for %s in %s.",
		       name,fDetectorTask[i-1].Data()));
         return 1;
       }
      i++;
    }
   return 0;
}

//____________________________________________________________________________
void    AliConfig::Add (char *list)
{
  // Adds space separated list of objects to the
  // configuration macro (?)
  char *path;
  
  const char   *confPath = gSystem->Getenv ("ALICE_CONFIG_PATH");
  if  (confPath) {
    path = new char[strlen (confPath)];
    strcpy (path, confPath);
  } else {
    const char   *alice = gSystem->Getenv ("ALICE_ROOT");
    path = new char[strlen (alice) + 32];
    
    strcpy (path, ".:");
    if (alice) {
      strcat (path, alice);
    }
    strcat (path, "/macros/config");
  }
  
  char   *token = strtok (path, ":");
  
  TList  *dirlist = new TList;
  
  while (token != NULL)	
    {
      dirlist->Add (new TObjString(token));
      token = strtok (NULL, ":");
    }
  
  token = strtok (list, " ");
  
  while (token != NULL)
    { 
      AliInfo(Form("Configuring token=%s",token));
      
      TObject *obj;
      TIter   next (dirlist);
      TString found = "\0";
      
      while ((obj = next ()))
        {
          TString dir(obj->GetName());
          TString apath  = dir + "/" + token;
          TString macro = apath + ".C";
          if (!gSystem->AccessPathName (macro.Data()))	
           {
            gInterpreter->ExecuteMacro (macro.Data());				   
            found = "(" + macro + ")";
            if (macro.Contains("/")) 
             {
               TString dirname = gSystem->DirName(macro.Data());
               TString macroConfigure = dirname + "/Configure.C";
               if (!gSystem->AccessPathName (macroConfigure.Data()))
                {
                  gInterpreter->ExecuteMacro (macroConfigure.Data());				    
                  found += " => Configured";
                }			      
             }
            break;
           } 
          else 
           {
            TString macroDefault = apath + "/Default.C";
            if (!gSystem->AccessPathName (macroDefault.Data()))
              {
                gInterpreter->ExecuteMacro (macroDefault.Data());
                found = "(" + macro + ")";
                TString macroConfigure = apath + "/Configure.C";
                if (!gSystem->AccessPathName (macroConfigure.Data()))	
                  {
                    gInterpreter->ExecuteMacro (macroConfigure.Data());				    
                    found += " => Configured";
                  }
                break; 				    
              }
           }
        }
      
      if (strlen(found.Data())) 
        {
          AliInfo(Form("found=%s  => OK",found.Data()));
        } 
      else 
        {
          AliError(" => FAILED.");
          exit(1); 
        }   	    
      
      token = strtok (NULL," ");
    }
  
  if (dirlist) delete dirlist;
  
}

/*****************************************************************************/

TFolder* AliConfig::BuildEventFolder(const char* name,const char* title)
{
/*
 creates the folder structure for one event
 TopFolder_
         | \
         |  Tasks
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
  conditions->AddFolder(fgkQAFolderName,"Quality Asurance Output"); //Folder with output of the QA task(s)
  //Configuration
  TFolder *configuration = eventfolder->AddFolder(fgkConfigurationFolderName, "Run configuration");
  configuration->AddFolder(fgkFieldFolderName, "Magnetic field maps");
  configuration->AddFolder(fgkGeneratorsFolderName,"list of generator objects");
  configuration->AddFolder(fgkVirtualMCFolderName,"the Virtual MC");

  eventfolder->AddFolder(fgkHeaderFolderName,"MonteCarlo event header");

  eventfolder->SetOwner();

  return eventfolder;
}

/*****************************************************************************/

TString AliConfig::GetQATaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskQA];
 }
/*****************************************************************************/
 
TString AliConfig::GetDigitizerTaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskDigitizer];
 }
/*****************************************************************************/
 
TString AliConfig::GetSDigitizerTaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskSDigitizer];
 }
/*****************************************************************************/

TString AliConfig::GetReconstructionerTaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskRecontructioner];
 }
/*****************************************************************************/

TString AliConfig::GetTrackerTaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskTracker];
 }
/*****************************************************************************/

TString AliConfig::GetPIDTaskName() const
 {
 //returns task name
  return fDetectorTask[kDetTaskPID];
 }
/*****************************************************************************/

const TString& AliConfig::GetQAFolderName() const
{
//returns pathname of folder with QA output relative to Top Alice Folder
  return fDetectorFolder[kDetFolderQA];
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
