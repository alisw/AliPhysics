#ifndef TOFv1_H
#define TOFv1_H
///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 1  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv1 : public AliTOF {
 
protected:
  Int_t fIdFBT2; // First sensitive element identifier
  Int_t fIdFBT3; // Second sensitive element identifier
 
public:
  AliTOFv1();
  AliTOFv1(const char *name, const char *title);
  virtual       ~AliTOFv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   StepManager();
  virtual void   DrawModule();
  
  ClassDef(AliTOFv1,1)  // Time Of Flight version 1
};
 
#endif
