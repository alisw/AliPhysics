#ifndef TOFv2_H
#define TOFv2_H
///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 2  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv2 : public AliTOF {

private:
  Int_t fIdFTO2; // First sensitive volume identifier
  Int_t fIdFTO3; // Second sensitive volume identifier
  Int_t fIdFLT1; // Third sensitive volume identifier
  Int_t fIdFLT2; // Fourth sensitive volume identifier
  Int_t fIdFLT3; // Fifth sensitive volume identifier
 
public:
  AliTOFv2();
  AliTOFv2(const char *name, const char *title);
  virtual       ~AliTOFv2() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 2;}
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t);
  virtual void   StepManager();
  virtual void   DrawDetector();
 
   ClassDef(AliTOFv2,1)  //Time Of Flight version 2
};
 
#endif
