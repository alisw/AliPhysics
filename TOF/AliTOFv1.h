#ifndef TOFv1_H
#define TOFv1_H
///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 1  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv1 : public AliTOF {

private:
  Int_t fIdFTO2; // First sensitive volume identifier
  Int_t fIdFTO3; // Second sensitive volume identifier
  Int_t fIdFLT1; // Third sensitive volume identifier
  Int_t fIdFLT2; // Fourth sensitive volume identifier
  Int_t fIdFLT3; // Fifth sensitive volume identifier
 
public:
  AliTOFv1();
  AliTOFv1(const char *name, const char *title);
  virtual       ~AliTOFv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t);
  virtual void   StepManager();
  virtual void   DrawModule();
 
   ClassDef(AliTOFv1,1)  //Time Of Flight version 1
};
 
#endif
