#ifndef ALIITSTRACKV2PID_H
#define ALIITSTRACKV2PID_H

#include <TMath.h>
#include <iostream.h>
#include "AliITStrackV2.h"

//_____________________________________________________________________________
class AliITStrackV2Pid : public TObject {
public:
    AliITStrackV2Pid();
    virtual ~AliITStrackV2Pid(){}
public:
    Float_t fSignal,fMom;
    Int_t   fPcode;
    Float_t fWpi,fWk,fWp;

  ClassDef(AliITStrackV2Pid,1)  // ITS trackV2 PID
};

#endif


