#ifndef ALIJOB_H
#define ALIJOB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TFolder.h"
#include "TObjArray.h"

class AliJob : public TTask
{
 public :
  AliJob(const char* name="AliJob",const char* title=""); // Constructor
  virtual ~AliJob();                                      // Destructor
  void ListEnvironment();                                 // Provide listing of the job environment
  void ExecuteJob(Int_t mode=0);                          // Invokation of the top level processing
  void MakeFolder();                                      // Select creation of the folder structure 
  TFolder* GetFolder() const;                             // Provide pointer to the whiteboard folder 
  TObject* GetMainObject() const;                         // Provide pointer to the main object structure
  void AddObject(TObject* obj);                           // Add an object into the environment
  void AddObjects(TObjArray* arr);                        // Add all array objects into the environment
  void RemoveObject(TObject* obj);                        // Remove an object from the environment
  void RemoveObjects(const char* classname);              // Remove all objects inheriting from "classname"
  TObject* GetObject(const char* classname) const;        // Provide first stored object inheriting from "classname" 
  TObject* GetObject(Int_t j) const;                      // Provide j-th stored object
  TObjArray* GetObjects() const;                          // Provide all stored object pointers
  TObjArray* GetObjects(const char* classname);           // Provide all objects inheriting from "classname" 
  void ProcessObject(TObject* obj);                       // Process all sub-tasks for the specified object

 protected :
  Int_t fMakefolder;    // Flag to indicate creation of the folder structure
  TFolder* fFolder;     // Pointer to the folder which serves as the job's whiteboard
  TObject* fMainObject; // Pointer to the main processing object structure within the job
  TObjArray* fObjects;  // Pointers to the various user-added objects 
  TObjArray* fSelect;   //! Temp. array of pointers to user-selected stored objects 

  void SetMainObject(TObject* obj); // Store pointer to the main object structure

 ClassDef(AliJob,4) // Base class for top level job in a task based procedure 
};
#endif
