#ifndef TOFv0_H
#define TOFv0_H
///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 0  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv0 : public AliTOF {
  
protected:
  Int_t fIdFBT2; // Identifier of the first sensitive volume
  Int_t fIdFBT3; // Identifier of the second sensitive volume
  
public:
  AliTOFv0();
  AliTOFv0(const char *name, const char *title);
  virtual      ~AliTOFv0() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  StepManager();
  virtual void    DrawDetector();
  
  ClassDef(AliTOFv0,1)  // Time Of Flight version 0
};
 
#endif
