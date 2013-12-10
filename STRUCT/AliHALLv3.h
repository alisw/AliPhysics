#ifndef ALIHALLV3_H
#define ALIHALLV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: HALL          //
////////////////////////////////////////////////
 
#include "AliHALL.h"
 
 
class AliHALLv3 : public AliHALL {
 
public:
   AliHALLv3();
   AliHALLv3(const char *name, const char *title);
   virtual      ~AliHALLv3() {}
   virtual void  CreateGeometry();
   virtual void  StepManager();
   virtual void  Init();
   virtual void  SetNewShield24() {fNewShield24 = 1;}
   virtual void  SetScoring()     {fScoring     = 1;}
   virtual void  SetRackShield()  {fRackShield  = 1;}
private:
   Bool_t fNewShield24;   // Option for new shielding in PX24 and RB24
   Int_t  fRefVolumeId;   // Volume ID of scoring plane
   Bool_t fScoring;       // Scoring Option
   Bool_t fRackShield;    // Additional rack shielding
   ClassDef(AliHALLv3,2)  //Class for ALICE experimental hall
};

#endif
