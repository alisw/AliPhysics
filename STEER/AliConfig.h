#ifndef ALICONFIG_H
#define ALICONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


class TDatabasePDG;
class TFolder;
class TString;
class TVirtualMC;

class AliConfig;
class AliDetector;
class AliGenerator;
class AliModule;
class AliTasks;

#include <TNamed.h>

class AliConfig : public TNamed {
  
public:

  AliConfig();
  virtual ~ AliConfig ();
  
  void  Add (AliGenerator *generator);
  void  Add (TVirtualMC *mc);
  void  Add (TDatabasePDG *pdg);
  void  Add (AliModule *module);
  void  Add (AliDetector *detector);
  
  void  Add (char *list);
  
  static AliConfig* Instance();
  
private:

  enum {kFolders=8, kTasks=5};
  AliConfig(const char * name, const char * title);
  AliConfig(const AliConfig& conf);
  void  AddInFolder (const char * dir, TObject *obj);
  void  AddSubFolder(const char * dir[], TObject *obj);
  void  AddSubTask(const char * dir[], TObject *obj);
  TObject* FindInFolder (const char *dir, const char *name);
  AliConfig& operator = (const AliConfig&) {return *this;}
  
  TFolder  *fTopFolder;                // Pointer of the top folder
  AliTasks *fTasks;                    // Pointer for the tasks
  // folders
  const char*  fPDGFolder ;            // Names of the PDG folders
  const char*  fGeneratorFolder ;      // Names of the Generator folders
  const char*  fMCFolder ;             // Names of MC folders
  const char*  fModuleFolder ;         // Names of Module folders
  const char** fDetectorFolder ;       // Names of Detector folders
  const char** fDetectorTask ;         // Names of Detector Task folders

  static AliConfig*  fgInstance;       // Instance of the singleton
    
    ClassDef(AliConfig,1) //Configuration class for AliRun
};				// end class AliConfig

#endif
