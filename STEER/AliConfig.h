#ifndef ALICONFIG_H
#define ALICONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* 
 * $Log$
 * Revision 1.5  2002/10/14 14:57:32  hristov
 * Merging the VirtualMC branch to the main development branch (HEAD)
 *
 * Revision 1.4.8.1  2002/06/10 14:43:06  hristov
 * Merged with v3-08-02
 *
 * Revision 1.4  2001/10/05 12:11:40  hristov
 * iostream.h used instead of iostream (HP)
 *
 * Revision 1.3  2001/10/04 15:30:56  hristov
 * Changes to accommodate the set of PHOS folders and tasks (Y.Schutz)
 *
 * Revision 1.2  2001/05/21 17:22:51  buncic
 * Fixed problem with missing AliConfig while reading galice.root
 *
 * Revision 1.1  2001/05/16 14:57:22  alibrary
 * New files for folders and Stack
 * 
 */

#include "AliMC.h"

class TDatabasePDG;
class TFolder;
class TString;

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
  void  Add (AliMC *mc);
  void  Add (TDatabasePDG *pdg);
  void  Add (AliModule *module);
  void  Add (AliDetector *detector);
  
  void  Add (char *list);
  
  static AliConfig* Instance();
  
private:

  enum {kFolders=8, kTasks=5};
  AliConfig(const char * name, const char * title);
  AliConfig(const AliConfig&);
  void  AddInFolder (char * dir, TObject *obj);
  void  AddSubFolder(char * dir[], TObject *obj);
  void  AddSubTask(char * dir[], TObject *obj);
  TObject* FindInFolder (char *dir, const char *name);
  AliConfig& operator = (const AliConfig&) {return *this;}
  
  TFolder  *fTopFolder;
  AliTasks *fTasks;
  // folders
  char*  fPDGFolder ; 
  char*  fGeneratorFolder ; 
  char*  fMCFolder ; 
  char*  fModuleFolder ; 
  char** fDetectorFolder ; 
  char** fDetectorTask ; 

  static AliConfig*  fInstance;
    
    ClassDef(AliConfig,1) //Configuration class for AliRun
};				// end class AliConfig

#endif
