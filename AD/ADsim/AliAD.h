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

class AliADCalibData;

  
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
  
  virtual    void       Hits2Digits();
  virtual    void       Hits2SDigits();
  virtual    void       Digits2Raw();
  virtual    Bool_t     Raw2SDigits(AliRawReader*);
  virtual    void 	SetADATwoInstalled(Bool_t b){fSetADATwoInstalled = b;} // ecv
  virtual    void  	SetADCTwoInstalled(Bool_t b){fSetADCTwoInstalled = b;} // ecv
  virtual    Bool_t 	GetADATwoInstalled() const {return fSetADATwoInstalled;}  // ecv
  virtual    Bool_t 	GetADCTwoInstalled() const {return fSetADCTwoInstalled;}  // ecv
  void                  GetCalibData();


private:
                       AliAD(const AliAD&); 
                       AliAD& operator = (const AliAD&);
  AliADCalibData *fCalibData;      //! Pointer to the calibration object 
  Bool_t	fSetADATwoInstalled; 
  Bool_t	fSetADCTwoInstalled; 


  ClassDef(AliAD,1)  // Base Class for the AD detector
};

#endif
