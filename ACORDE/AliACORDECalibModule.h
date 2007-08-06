#ifndef AliACORDECalibModule_H
#define AliACORDECalibModule_H

////////////////////////////////////////////////
//  class for ACORDE calibration                 //
////////////////////////////////////////////////

#include "TObject.h"
#include "TH1F.h"

class AliACORDEDataModule;

class AliACORDECalibModule: public TObject 
{

 public:
  enum {kNmodules=60};
  AliACORDECalibModule();
  AliACORDECalibModule(const AliACORDECalibModule &calibdata);
  AliACORDECalibModule& operator=(const AliACORDECalibModule & calibdata);
  virtual ~AliACORDECalibModule();
  void Print_Module();
  AliACORDEDataModule* GetModule(Int_t module) const {return fModule[module];}
  void SetModuleRate(Int_t module,Float_t value);
  Float_t GetModuleRate(Int_t module);
  void SetModule(Int_t module, Float_t value, Bool_t status,const char* name);
  void Create_Histo();
  void Draw(const Option_t* /*option*/);
  

 protected:
  AliACORDEDataModule  *fModule[kNmodules];  
  TH1F   *fHRate;
  ClassDef(AliACORDECalibModule,1)    // ACORDE Calibration data
};

#endif
