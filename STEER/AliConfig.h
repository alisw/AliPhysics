#ifndef ALICONFIG_H
#define ALICONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
class TDatabasePDG;
class TFolder;
class TString;
class TVirtualMC;

class AliConfig;
class AliDetector;
class AliGenerator;
class AliModule;
class AliTasks;

class AliConfig : public TNamed {
  
public:
  
  AliConfig();
  
  virtual ~ AliConfig (); 

  void       Add(TDatabasePDG *pdg);
  void       Add(char *list);
  
  void       Add(AliGenerator *generator,const char* eventfolder = fgkDefaultEventFolderName);
  void       Add (TVirtualMC *mc,const char* eventfolder = fgkDefaultEventFolderName);
  void       Add (AliModule *module,const char* eventfolder = fgkDefaultEventFolderName);
  void       Add (AliDetector *detector,const char* eventfolder = fgkDefaultEventFolderName);

  Int_t      AddDetector(const char* evntfoldername,const char *name, const char* title);
  Int_t      AddDetector(TFolder* evntfolder,const char *name, const char* title);
  
  Int_t      CreateDetectorFolders(const char* evntfoldername,const char *name, const char* title);//Used by AliRunGetter
  Int_t      CreateDetectorFolders(TFolder* evntfolder,const char *name, const char* title);//Used by AliRunGetter
  Int_t      CreateDetectorTasks(const char *name, const char* title);
  
  static     AliConfig* Instance();
  
private:
  AliConfig(const char * name, const char * title );
  AliConfig(const AliConfig&);

  void          AddInFolder (const char * dir, TObject *obj);
  Int_t         AddSubTask(const char *taskname, const char* name, const char* title);
  Int_t         AddSubFolder(TFolder* topfolder, const char* infoler, //helper method
                     const char* newfoldname, const char* newfoldtitle);
  TObject*      FindInFolder (const char *dir, const char *name);
  
  // folders
  TFolder*              fTopFolder;
  TFolder*              fTaskFolder;
  TFolder*              fConstFolder;

  static const TString  fgkPDGFolderName; 
  static const TString  fgkGeneratorFolderName; 
  static const TString  fgkMCFolderName;
  static const TString  fgkModuleFolderName;
  
  TString              *fDetectorTask;//!array with names for detector tasks
  TString              *fDetectorFolder;//!array with names for detector folders (where detector is going to be put)
  
  static AliConfig*     fInstance;
  
  
 public:

  TFolder*              BuildEventFolder(const char* name,const char* tilte);
  
  TFolder*              GetTopFolder(){return fTopFolder;}
  TFolder*              GetTaskFolder(){return fTaskFolder;}
  TFolder*              GetConstFolder(){return fConstFolder;}

  static const TString& GetModulesFolderName(){return fgkModuleFolderName;}
  
  TString               GetQATaskName() const; //returns path to QA tasks
  TString               GetDigitizerTaskName () const;
  TString               GetSDigitizerTaskName () const;
  TString               GetReconstructionerTaskName () const;
  TString               GetTrackerTaskName () const;
  TString               GetPIDTaskName () const;
  
  
  const TString&        GetQAFolderName() const; //returns path to folder with QA output
  
  const TString&        GetDataFolderName();//returns name of data folder
  
  static const TString  fgkTopFolderName; //name of top AliRoot folder
 
  static const TString  fgkDefaultEventFolderName; 
  static const TString  fgkTasksFolderName;
  static const TString  fgkConstantsFolderName;
  
  static const TString  fgkDataFolderName;  
  static const TString  fgkConditionsFolderName;
  static const TString  fgkConfigurationFolderName;
  static const TString  fgkHeaderFolderName;
  
  static const TString  fgkDigitizerTaskName;
  static const TString  fgkSDigitizerTaskName;
  static const TString  fgkQATaskName;
  static const TString  fgkReconstructionerTaskName;
  static const TString  fgkTrackerTaskName;
  static const TString  fgkPIDTaskName;
  
  static const TString  fgkCalibrationFolderName;
  static const TString  fgkAligmentFolderName;
  static const TString  fgkQAFolderName;
  
  static const TString  fgkFieldFolderName;
  static const TString  fgkGeneratorsFolderName;
  static const TString  fgkVirtualMCFolderName;
  
  ClassDef(AliConfig,2) //Configuration class for AliRun
};				// end class AliConfig

#endif
