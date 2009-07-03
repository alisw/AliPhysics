#ifndef ALIESDZDCSCALERS_H
#define ALIESDZDCSCALERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//  This class is used to store data in events where scalers are read
//  Scaler format: 1(in)+1(out) floats from Benotto's card downstairs
//		   4(in)+4(out) floats from Benotto's card upstairs
//		   32 floats from VME scaler
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>


class AliESDZDCScalers: public TObject {

public:
  AliESDZDCScalers();
  AliESDZDCScalers(const AliESDZDCScalers& zdc);
  AliESDZDCScalers& operator=(const AliESDZDCScalers& zdc);
  virtual void Copy(TObject &obj) const;
  
  UInt_t GetScalerDown(Int_t i) const {return fScalerDown[i];}
  UInt_t GetScalerUp(Int_t i)   const {return fScalerUp[i];}
  UInt_t GetVMEScaler(Int_t i)  const {return fVMEScaler[i];}
  const UInt_t* GetScalerDown() const {return fScalerDown;}
  const UInt_t* GetScalerUp()   const {return fScalerUp;}
  const UInt_t* GetVMEScaler()  const {return fVMEScaler;}
  
  void SetScalerDown(UInt_t count[2]) 
  	{for(Int_t k=0; k<2; k++) fScalerDown[k] = count[k];}
  void SetScalerUp(UInt_t count[8]) 
  	{for(Int_t k=0; k<8; k++) fScalerUp[k] = count[k];}
  void SetVMEScaler(UInt_t count[32]) 
  	{for(Int_t k=0; k<32; k++) fVMEScaler[k] = count[k];}

  void    Reset();
  void    Print(const Option_t *opt=0) const;

private:
  UInt_t fScalerDown[2]; // counts from Benotto's card downstairs
  UInt_t fScalerUp[8];   // counts from Benotto's card upstairs
  UInt_t fVMEScaler[32]; // counts from VME scaler

  ClassDef(AliESDZDCScalers,1)
  
};

#endif
