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
  AliConfig(const AliConfig&);
  void  AddInFolder (const char * dir, TObject *obj);
  void  AddSubFolder(const char * dir[], TObject *obj);
  void  AddSubTask(const char * dir[], TObject *obj);
  TObject* FindInFolder (const char *dir, const char *name);
  AliConfig& operator = (const AliConfig&) {return *this;}
  
  TFolder  *fTopFolder;
  AliTasks *fTasks;
  // folders
  const char*  fPDGFolder ; 
  const char*  fGeneratorFolder ; 
  const char*  fMCFolder ; 
  const char*  fModuleFolder ; 
  const char** fDetectorFolder ; 
  const char** fDetectorTask ; 

  static AliConfig*  fInstance;
    
    ClassDef(AliConfig,1) //Configuration class for AliRun
};				// end class AliConfig

#endif
