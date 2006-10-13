#ifndef ALIITSTRACKV2PID_H
#define ALIITSTRACKV2PID_H
/////////////////////////////////////////
// ITS trackV2 pid class               //
/////////////////////////////////////////                     
#include <TMath.h>
#include <Riostream.h>
#include "AliITStrackV2.h"

//_____________________________________________________________________________
class AliITStrackV2Pid : public TObject {
public:
    AliITStrackV2Pid();
    virtual ~AliITStrackV2Pid(){}
protected:
    Float_t fWpi,fWk,fWp,fSignal,fMom;
    Float_t fPhi,fLam;
    Int_t   fPcode;
    
    Float_t fGSignal,fGMom,fGpx,fGpy,fGpz,fGx,fGy,fGz;
    Int_t   fGcode,fGlab,fFound;

    Float_t fQ1,fQ2,fQ3,fQ4,fQ5,fQ6;
    Float_t fD,fZ;
    
  ClassDef(AliITStrackV2Pid,2)  // ITS trackV2 PID
};

#endif


