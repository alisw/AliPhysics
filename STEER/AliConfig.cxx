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

/*
$Log$
Revision 1.6  2002/10/22 15:02:15  alibrary
Introducing Riostream.h

Revision 1.5  2002/10/14 14:57:32  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.3.8.1  2002/06/10 14:43:06  hristov
Merged with v3-08-02

Revision 1.4  2002/05/27 14:26:59  hristov
New folder for track references added

Revision 1.3  2001/10/04 15:30:56  hristov
Changes to accommodate the set of PHOS folders and tasks (Y.Schutz)

Revision 1.2  2001/05/21 17:22:50  buncic
Fixed problem with missing AliConfig while reading galice.root

Revision 1.1  2001/05/16 14:57:22  alibrary
New files for folders and Stack

*/

#include <Riostream.h>
#include <TDatabasePDG.h>
#include <TFolder.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliGenerator.h" 
#include "TObjString.h" 
#include "TString.h"
#include "TTask.h" 

ClassImp(AliConfig)

AliConfig* AliConfig::fInstance = 0;

//____________________________________________________________________________
AliConfig* AliConfig::Instance ()
{
  //
  // Instance method for singleton class
  //
  if(!fInstance) fInstance = new AliConfig ("Folders","Alice data exchange");

  return fInstance;
}

//____________________________________________________________________________
AliConfig::AliConfig():
  fTopFolder(0),
  fTasks(0),
  fPDGFolder(0),
  fGeneratorFolder(0),
  fMCFolder(0),
  fModuleFolder(0),
  fDetectorFolder(0),
  fDetectorTask(0)
{
  //
  // Default constructor, mainly to keep coding conventions
  //
  fInstance=0;
    
  Fatal("ctor",
   "Constructor should not be called for a singleton\n");
}

//____________________________________________________________________________
AliConfig::AliConfig(const AliConfig&):
  fTopFolder(0),
  fTasks(0),
  fPDGFolder(0),
  fGeneratorFolder(0),
  fMCFolder(0),
  fModuleFolder(0),
  fDetectorFolder(0),
  fDetectorTask(0)
{
  //
  // Copy constructor, mainly to keep coding conventions
  //
  fInstance=0;
    
  Fatal("copy ctor",
   "Copy constructor should not be called for a singleton\n");
}

//____________________________________________________________________________
AliConfig::AliConfig(const char *name, const char *title): 
  TNamed(name,title), 
  fTopFolder(gROOT->GetRootFolder ()->AddFolder (name,title)),
  fTasks(0),
  fPDGFolder("Constants/DatabasePDG"),
  fGeneratorFolder("RunMC/Configuration/Generators"),
  fMCFolder("RunMC/Configuration/VirtualMC"),
  fModuleFolder("Run/Configuration/Modules"),
  fDetectorFolder(new char*[kFolders]),
  fDetectorTask(new char*[kTasks])
{
  //
  // Constructor
  //
  fInstance=this;
  
  fDetectorFolder[0] = "Run/Conditions/Calibration" ;
  fDetectorFolder[1] = "Run/Event/Data" ;
  fDetectorFolder[2] = "Run/Event/RecData" ;
  fDetectorFolder[3] = "RunMC/Event/Data/Hits" ;
  fDetectorFolder[4] = "RunMC/Event/Data/SDigits" ;
  fDetectorFolder[5] = "Run/Conditions/QA" ;  
  fDetectorFolder[6] = "RunMC/Event/Data/TrackReferences" ;  
  fDetectorFolder[7] = 0 ;  
  fDetectorTask[0] = "Tasks/QA" ;  
  fDetectorTask[1] = "Tasks/SDigitizer" ;  
  fDetectorTask[2] = "Tasks/Digitizer" ;  
  fDetectorTask[3] = "Tasks/Reconstructioner" ;  
  fDetectorTask[4] = 0 ;  

  fTopFolder->SetOwner() ; 
  gROOT->GetListOfBrowsables ()->Add (fTopFolder, name);
  
  TFolder *subfolder;
  
  TFolder *constants =
    fTopFolder->AddFolder ("Constants", "Detector constants");
  
  subfolder = 
    constants->AddFolder ("DatabasePDG", "PDG database");
  
  TFolder *run = 
    fTopFolder->AddFolder ("Run", "Run dependent folders");
  
  TFolder *conditions = 
    run->AddFolder ("Conditions", "Run conditions");
  
  subfolder =
    conditions->AddFolder ("Calibration","Detector calibration data");
  
  subfolder =
    conditions->AddFolder ("Aligment", "Detector aligment");
  
  subfolder =
    conditions->AddFolder ("QA", "Detector QA");
  
  TFolder *configuration =
    run->AddFolder ("Configuration", "Run configuration");
  
  subfolder =
    configuration->AddFolder ("Modules", "Detector objects");
  
  subfolder =
    configuration->AddFolder ("Field", "Magnetic field maps");
  
  TFolder *event = 
    run->AddFolder ("Event", "Event folders");
  
  subfolder = 
    event->AddFolder ("Data", "Detector raw data");
  
  subfolder =
    event->AddFolder ("RecData", "Detectors reconstucted data");
  
  TFolder *run_mc =
    fTopFolder->AddFolder ("RunMC", "MonteCarlo run dependent folders");
  
  TFolder *configuration_mc =
    run_mc->AddFolder ("Configuration","MonteCarlo run configuration");
  
  subfolder =
    configuration_mc->AddFolder ("Generators","list of generator objects");
  
  subfolder =
    configuration_mc->AddFolder ("VirtualMC", "the Virtual MC");
  
  TFolder *event_mc =
    run_mc->AddFolder ("Event", "MonteCarlo event folders");
  
  subfolder =
    event_mc->AddFolder ("Header", "MonteCarlo event header");
  
  //		subfolder =
  //		    event_mc->AddFolder ("Kinematics", "MonteCarlo generated particles");
  
  TFolder *data_mc =
    event_mc->AddFolder ("Data", "MonteCarlo data");
  
  subfolder = 
    data_mc->AddFolder ("Hits", "MonteCarlo Hits") ; 
 
 subfolder = 
    data_mc->AddFolder ("SDigits", "MonteCarlo SDigits") ; 

 subfolder = 
    data_mc->AddFolder ("TrackReferences", "MonteCarlo track references") ; 


  
  // Add the tasks to //Folders
  
  TFolder * tasksfolder = fTopFolder->AddFolder("Tasks", "ALICE Tasks") ; 
  
  TTask * qa = new TTask("QA", "Alice QA tasks") ;
  tasksfolder->Add(qa); 
  
  TTask * sd = new TTask("SDigitizer", "Alice SDigitizer") ;
  tasksfolder->Add(sd); 

  TTask * di = new TTask("Digitizer", "Alice Digitizer") ;
  tasksfolder->Add(di); 

  TTask * re = new TTask("Reconstructioner", "Alice Reconstructioner") ;
  tasksfolder->Add(re); 

}

//____________________________________________________________________________
AliConfig::~AliConfig()
{ 
  delete [] fDetectorFolder ;  
  delete fDetectorTask ;  
  delete fTopFolder ; 
}

//____________________________________________________________________________
void AliConfig::AddInFolder (char *dir, TObject *obj)
{
  TFolder *folder =
    dynamic_cast<TFolder *>(fTopFolder->FindObject(dir));
  if (folder)
    folder->Add (static_cast<TObject *>(obj));
}

//____________________________________________________________________________
void    AliConfig::AddSubFolder(char * dir[], TObject *obj)
{
  int iDir = 0;
  
  while (dir[iDir]) 
    {
      TFolder * folder = dynamic_cast<TFolder *>(fTopFolder->FindObject (dir[iDir++]));
      if (folder) {
	TFolder * subfolder = dynamic_cast<TFolder *>(folder->FindObject (obj->GetName()));
	if (!subfolder)
	  subfolder = folder->AddFolder (obj->GetName(),obj->GetTitle());			   
      }
    }
}

//____________________________________________________________________________
void    AliConfig::AddSubTask(char * dir[], TObject *obj)
{
  int iDir = 0;
  
  while (dir[iDir]) 
    {
      TTask * task = dynamic_cast<TTask *>(fTopFolder->FindObject (dir[iDir++]));
      if (task) {
	TTask * subtask = static_cast<TTask*>(task->GetListOfTasks()->FindObject (obj->GetName()));
	if (!subtask) {
	  subtask = new TTask(obj->GetName(), obj->GetTitle()) ; 
	  task->Add(subtask);			   
	}
      }
    }
}

//____________________________________________________________________________
TObject* AliConfig::FindInFolder (char *dir, const char *name)
{
  if(!name) return(fTopFolder->FindObject(name));
  TFolder * folder = dynamic_cast<TFolder *>(fTopFolder->FindObject(dir));
  if (!folder) return (NULL);
  return(folder->FindObject(name));
}

//____________________________________________________________________________
void    AliConfig::Add (AliGenerator * obj)
{
  AddInFolder(fGeneratorFolder, obj);
}

//____________________________________________________________________________
void    AliConfig::Add (AliMC * obj)
{
  AddInFolder(fMCFolder, obj);
}

//____________________________________________________________________________
void    AliConfig::Add (TDatabasePDG * obj)
{
  AddInFolder(fPDGFolder, obj);
}

//____________________________________________________________________________
void    AliConfig::Add (AliModule* obj)
{
  AddInFolder(fModuleFolder, obj);
}

//____________________________________________________________________________
void    AliConfig::Add (AliDetector * obj)
{
  AddSubFolder(fDetectorFolder, obj); 
  AddSubTask(fDetectorTask, obj);
}


//____________________________________________________________________________
void    AliConfig::Add (char *list)
{
  char *path;
  
  const char   *conf_path = gSystem->Getenv ("ALICE_CONFIG_PATH");
  if  (conf_path) {
    path = new char[strlen (conf_path)];
    strcpy (path, conf_path);
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
      cout << "Configuring " << token << ": ";
      
      TObject *obj;
      TIter   next (dirlist);
      TString found = "\0";
      
      while ((obj = next ()))
	{
	  TString dir(obj->GetName());
	  TString path  = dir + "/" + token;
	  TString macro = path + ".C";
	  if (!gSystem->AccessPathName (macro.Data()))	{
	    gInterpreter->ExecuteMacro (macro.Data());				   
	    found = "(" + macro + ")";
	    if (macro.Contains("/")) {
	      TString dirname = gSystem->DirName(macro.Data());
	      TString macroConfigure = dirname + "/Configure.C";
	      if (!gSystem->AccessPathName (macroConfigure.Data()))	{
		gInterpreter->ExecuteMacro (macroConfigure.Data());				    
		found += " => Configured";
	      }			      
	    }
	    break;
	  } else {
	    TString macroDefault = path + "/Default.C";
	    if (!gSystem->AccessPathName (macroDefault.Data()))	{
	      gInterpreter->ExecuteMacro (macroDefault.Data());
	      found = "(" + macro + ")";
	      TString macroConfigure = path + "/Configure.C";
	      if (!gSystem->AccessPathName (macroConfigure.Data()))	{
		gInterpreter->ExecuteMacro (macroConfigure.Data());				    
		found += " => Configured";
	      }
	      break; 				    
	    }
	  }
	}
      
      if (strlen(found.Data())) {
	cout << found << " => OK" << endl;
      } else {
	cout << " => FAILED." << endl;
	exit(1); 
      }   	    
      
      token = strtok (NULL, " ");
    }
  
  if (dirlist) delete dirlist;
  
}


