#ifndef ALIMODULEINFO_H
#define ALIMODULEINFO_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE MODULE INFORMATION CLASS                                      //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

class AliModuleInfo{
 public:
  AliModuleInfo(int n);
  AliModuleInfo(const AliModuleInfo& rh);
  AliModuleInfo& operator = (const AliModuleInfo& rh);
  
  virtual ~AliModuleInfo();

  void     SetId(const char* name,Int_t id);
  void     Add(const char * name,Int_t i);
  Int_t    Id(const char *name) const;
  const char*    Name(Int_t id) const;
  Bool_t   IsEnabled(Int_t id) const;
  Bool_t   IsEnabled(const char *name) const {return IsEnabled(Id(name));};
  void     Disable(Int_t id);
  void     Disable(const char* name){Disable(Id(name));};
  void     Enable(Int_t id);
  void     Enable(const char *name){Enable(Id(name));};
  void     Print() const;
  
 private:
  //The purposes of this class is to link each module to its Id
  char     **fName; // Array containing the names of the modules
  Int_t    *fId;    // Array of module's Ids
  Bool_t   *fEnabled; // Array of flags to enable/disable modules
  Int_t    fNb; // Number of modules

  ClassDef(AliModuleInfo,0);
};

#endif
