#ifndef ALICONFIG_H
#define ALICONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* 
 * $Log$
 * Revision 1.8  2002/10/29 14:59:45  alibrary
 * Some more code cleanup
 *
 * Revision 1.7  2002/10/29 14:26:49  hristov
 * Code clean-up (F.Carminati)
 *
 * Revision 1.6  2002/10/22 15:02:15  alibrary
 * Introducing Riostream.h
 *
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
