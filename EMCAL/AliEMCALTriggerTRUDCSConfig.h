#ifndef ALIEMCALTRIGGERTRUDCSCONFIG_H
#define ALIEMCALTRIGGERTRUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
 
 
 Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "TObject.h"

class AliEMCALTriggerTRUDCSConfig : public TObject 
{
public:

	AliEMCALTriggerTRUDCSConfig();
	virtual ~AliEMCALTriggerTRUDCSConfig() {}

	void    SetSELPF( Int_t pf)              { fSELPF  = pf; }  
	void    SetL0SEL( Int_t la)              { fL0SEL  = la; }  
	void    SetL0COSM(Int_t lc)              { fL0COSM = lc; }
	void    SetGTHRL0(Int_t lg)              { fGTHRL0 = lg; }
	void    SetMaskReg(Int_t arr[], Int_t n) { for (Int_t i=0;i<n;i++) fMaskReg[i] = arr[i]; }

	Int_t   GetSELPF()                       const { return fSELPF;  }
	Int_t   GetL0SEL()                       const { return fL0SEL;  }
	Int_t   GetL0COSM()                      const { return fL0COSM; }
	Int_t   GetGTHRL0()                      const { return fGTHRL0; }
	void    GetMaskReg(Int_t arr[], Int_t n) const { for (Int_t i=0;i<n;i++) arr[i] = fMaskReg[i]; }
	
protected:

	AliEMCALTriggerTRUDCSConfig(const AliEMCALTriggerTRUDCSConfig &cd);
	AliEMCALTriggerTRUDCSConfig &operator=(const AliEMCALTriggerTRUDCSConfig &cd);

private:
	
	Int_t   fSELPF;                          // 
	Int_t   fL0SEL;                          // 
	Int_t   fL0COSM;                         // 
	Int_t   fGTHRL0;                         // 
	Int_t   fMaskReg[5];                     //
	
	ClassDef(AliEMCALTriggerTRUDCSConfig,1)  //
};
#endif
