#ifndef ALITPCTRACKPID_H
#define ALITPCTRACKPID_H

#include <TMath.h>
#include "Riostream.h"
#include "../ITS/AliITStrackV2.h"

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


