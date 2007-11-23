#ifndef AliAODPid_H
#define AliAODPid_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Pid object for additional pid information
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TObject.h>

class AliAODPid : public TObject {

 public:
  
  AliAODPid();
  virtual ~AliAODPid();
  AliAODPid(const AliAODPid& pid); 
  AliAODPid& operator=(const AliAODPid& pid);

 private :

  ClassDef(AliAODPid,1);
};

#endif
