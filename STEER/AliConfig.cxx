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
*/

#include "TString.h"

#include "AliConfig.h"
#include "AliDetector.h"

#include <iostream.h>

ClassImp(AliConfig)


static char* gPDGFolder= 
	    	"Constants/DatabasePDG";

static char* gGeneratorFolder =  
	    	"RunMC/Configuration/Generators";

static char* gMCFolder = 
	    	"RunMC/Configuration/VirtualMC";

static char* gModuleFolder =
	    	"Run/Configuration/Modules";

static char* gDetectorFolder[] = {
	        "Run/Conditions/Calibration",
	        "Run/Event/Data",
	        "Run/Event/RecData",
	        "RunMC/Event/Data",0};

AliConfig* AliConfig::fInstance = 0;

AliConfig* AliConfig::Instance ()
{
  //
  // Instance method for singleton class
  //
        if(fInstance == 0) {
		    fInstance = new AliConfig ();
		}
		return fInstance;
}


AliConfig::AliConfig(const char *name, const char *title)
{
  //
  // Default constructor
  //
  //
 	    fInstance=this;
	    
		fTopFolder = gROOT->GetRootFolder ()->AddFolder (name,title);

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

		subfolder =
		    event_mc->AddFolder ("Data", "MonteCarlo data");

}

AliConfig::~AliConfig()
{ 
}

void    AliConfig::AddInFolder (char *dir, TObject *obj)
{
	    TFolder *folder =
               (TFolder *) fTopFolder->FindObject (dir);
        if (folder)
               folder->Add ((TObject *)obj);
}

void    AliConfig::AddSubFolder(char *dir[], TObject *obj)
{
        int iDir = 0;
        
		while (dir[iDir]) 
		{
		    TFolder * folder = (TFolder *) fTopFolder->FindObject (dir[iDir++]);
		    if (folder) {
  		        TFolder * subfolder = (TFolder *) folder->FindObject (obj->GetName());
		       if (!subfolder)
			        subfolder = folder->AddFolder (obj->GetName(),obj->GetTitle());			   
            }
		}

}

TObject* AliConfig::FindInFolder (char *dir, const char *name)
{
      if(!name) return(fTopFolder->FindObject(name));
  	  TFolder * folder = (TFolder *) fTopFolder->FindObject (dir);
      if (!folder) return (NULL);
      return(folder->FindObject(name));
}

void    AliConfig::Add (AliGenerator * obj)
{
		AddInFolder(gGeneratorFolder, (TObject *) obj);
}

void    AliConfig::Add (AliMC * obj)
{
		AddInFolder(gMCFolder, (TObject *) obj);
}

void    AliConfig::Add (TDatabasePDG * obj)
{
		AddInFolder(gPDGFolder, (TObject *) obj);
}

void    AliConfig::Add (AliModule* obj)
{
		AddInFolder(gModuleFolder, (TObject *) obj);
}

void    AliConfig::Add (AliDetector * obj)
{
		AddSubFolder(gDetectorFolder, (TObject *) obj);
}

	
void    AliConfig::Add (const char *list)
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
		    dirlist->Add ((TObject *) token);
		    token = strtok (NULL, ":");
		}

		token = strtok ((char *)list, " ");

		while (token != NULL)
		{
		    cout << "Configuring " << token << ": ";

		    TObject *obj;
		    TIter   next (dirlist);
		    TString found = "\0";
		    
		    while ((obj = next ()))
		    {
		        TString dir   = (char *) obj;
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


