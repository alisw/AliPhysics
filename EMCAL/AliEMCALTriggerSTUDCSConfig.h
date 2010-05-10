#ifndef ALIEMCALTRIGGERSTUDCSCONFIG_H
#define ALIEMCALTRIGGERSTUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*


EMCAL STU DCS parameters to be stored in OCDB
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "TObject.h"

class AliEMCALTriggerSTUDCSConfig : public TObject 
{
 public:

  AliEMCALTriggerSTUDCSConfig();
  virtual ~AliEMCALTriggerSTUDCSConfig() {};
	  
  void    SetGA(Int_t ga)        { fGA         = ga; }
  void    SetGB(Int_t gb)        { fGB         = gb; }
  void    SetGC(Int_t gc)        { fGC         = gc; }
  void    SetJA(Int_t ja)        { fJA         = ja; }
  void    SetJB(Int_t jb)        { fJB         = jb; }
  void    SetJC(Int_t jc)        { fJC         = jc; }
  void    SetRawData(Int_t rd)   { fGetRawData = rd; }
  void    SetRegion(Int_t rg)    { fRegion     = rg; }
  void    SetFw(Int_t fv)        { fFw         = fv; }
	
  Int_t   GetGA()        const { return fGA;         }
  Int_t   GetGB()        const { return fGB;         }
  Int_t   GetGC()        const { return fGC;         }
  Int_t   GetJA()        const { return fJA;         }
  Int_t   GetJB()        const { return fJB;         }
  Int_t   GetJC()        const { return fJC;         }
  Int_t   GetRawData()   const { return fGetRawData; }
  Int_t   GetRegion()    const { return fRegion;     }
  Int_t   GetFw()        const { return fFw;         }
	
protected:

	AliEMCALTriggerSTUDCSConfig(const AliEMCALTriggerSTUDCSConfig &cd);
	AliEMCALTriggerSTUDCSConfig &operator=(const AliEMCALTriggerSTUDCSConfig &cd);

private:
	
  Int_t   fGA;         //
  Int_t   fGB;         //
  Int_t   fGC;         //
  Int_t   fJA;         //
  Int_t   fJB;         //
  Int_t   fJC;         //
  Int_t   fGetRawData; //
  Int_t   fRegion;     //
  Int_t   fFw;         //
  
  ClassDef(AliEMCALTriggerSTUDCSConfig,1) //
};
#endif

