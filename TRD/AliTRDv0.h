#ifndef TRDv0_H
#define TRDv0_H
////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 0    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"
 
class AliTRDv0 : public AliTRD {

public:
  AliTRDv0() {}
  AliTRDv0(const char *name, const char *title);
  virtual      ~AliTRDv0() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  StepManager();
  virtual void  Init();
  virtual void  DrawModule();
  
protected:
  Int_t        fIdSens1;    // 1st sensitive volume identifier
  Int_t        fIdSens2;    // 2nd sensitive volume identifier
  Int_t        fIdSens3;    // 3rd sensitive volume identifier
  
  ClassDef(AliTRDv0,1)      // Transition Radiation Detector version 0
};

#endif
