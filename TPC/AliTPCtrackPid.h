#ifndef ALITPCTRACKPID_H
#define ALITPCTRACKPID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TMath.h>
#include <TObject.h>
#include "Riostream.h"

//_____________________________________________________________________________
class AliTPCtrackPid : public TObject {
public:
    AliTPCtrackPid();
    virtual ~AliTPCtrackPid(){}
public:
    Float_t fWpi,fWk,fWp,fSignal,fMom,fPhi,fLam;
    Int_t   fPcode,fLabel;
    
    Float_t fGSignal,fGMom,fGpx,fGpy,fGpz,fGx,fGy,fGz;
    Int_t   fGcode,fGlab;

  ClassDef(AliTPCtrackPid,1)  // TPC track PID
};

#endif


