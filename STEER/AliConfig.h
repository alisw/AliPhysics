#ifndef ALICONFIG_H
#define ALICONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* 
 * $Log$ 
 */

#include <TFolder.h>
#include <TList.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TDatabasePDG.h>

class AliGenerator;
class AliModule;
class AliDetector;
class AliMC;
class AliConfig;
class AliTasks;

class AliConfig : public TNamed {

  public:

    AliConfig(const char *name="gAlice", 
    		  const char *title = "Alice simulation and reconstruction framework");
    virtual ~ AliConfig ();

    void  Add (AliGenerator *generator);
    void  Add (AliMC *mc);
    void  Add (TDatabasePDG *pdg);
    void  Add (AliModule *module);
    void  Add (AliDetector *detector);
   
    void  Add (const char *list);
    
    static AliConfig* Instance();

  private:
     void  AddInFolder (char *dir, TObject *obj);
     void  AddSubFolder(char *dir[], TObject *obj);
     TObject* FindInFolder (char *dir, const char *name);

    TFolder  *fTopFolder;
    AliTasks *fTasks;

    static AliConfig*  fInstance;
    
    ClassDef(AliConfig,1) 

};				// end class AliConfig

#endif
