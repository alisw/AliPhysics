#ifndef ALIAD_H
#define ALIAD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : AD       //
//////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliADLoader.h"
#include "AliADDigitizer.h"
#include "AliADTrigger.h"

  
class AliAD : public AliDetector {
 
public:

                        AliAD();
                        AliAD(const char *name, const char *title);
  virtual              ~AliAD();
  virtual	void	CreateMaterials();
  virtual      Int_t     IsVersion() const { return -1;}
  virtual      TString   Version() { return TString("");}
  virtual      void     SetTreeAddress();  
  virtual      void MakeBranch(Option_t* opt = "");
  virtual      AliLoader*    MakeLoader(const char* topfoldername);
  AliDigitizer*  CreateDigitizer(AliDigitizationInput* digInput) const;
  virtual AliTriggerDetector* CreateTriggerDetector() const { return new AliADTrigger();}
  
  virtual      void     Digits2Raw();
  virtual    Bool_t     Raw2SDigits(AliRawReader*);
  virtual    void 	SetADAToInstalled(Bool_t b){fSetADAToInstalled = b;}
  virtual    void  	SetADCToInstalled(Bool_t b){fSetADCToInstalled = b;}
  virtual    Bool_t 	GetADAToInstalled() const {return fSetADAToInstalled;}
  virtual    Bool_t 	GetADCToInstalled() const {return fSetADCToInstalled;}


private:
                       AliAD(const AliAD&); 
                       AliAD& operator = (const AliAD&); 
  Bool_t	fSetADAToInstalled; 
  Bool_t	fSetADCToInstalled; 


  ClassDef(AliAD,1)  // Base Class for the AD detector
};

#endif
