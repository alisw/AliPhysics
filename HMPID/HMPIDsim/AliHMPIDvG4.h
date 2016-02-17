#ifndef AliHMPIDvG4_h
#define AliHMPIDvG4_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//.
//HMPID class for new geometry based on TGeo to test G4 production
//.
//.

#include "AliHMPIDv3.h"             //base class 
  
class AliHMPIDvG4 : public AliHMPIDv3 //TObject-TNamed-AliModule-AliDetector-AliHMPID-AliHMPIDv3
{
public:
 AliHMPIDvG4()                                   :AliHMPIDv3(          ),fMaxFeed(10) {;}          // default ctor
 AliHMPIDvG4(const char *name, const char *title):AliHMPIDv3(name,title),fMaxFeed(10) {;}          // named ctor
  virtual       ~AliHMPIDvG4()                                                         {;}            // dtor
          void    StepManager      (                                 );                               // from AliModule invoked from AliMC::Stepping()
          void    GenFee           (Float_t qtot                     );                               // generates feedback photons
	  void    SetMaxFeed(Int_t max) {fMaxFeed = max;}
 private:
	  Int_t fMaxFeed;
  ClassDef(AliHMPIDvG4,1)                                                                //HMPID full version for simulation
};

#endif
