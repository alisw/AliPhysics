#ifndef TRDv1_H
#define TRDv1_H
////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 1    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"
             
class AliTRDv1 : public AliTRD {

public:
  AliTRDv1() {}
  AliTRDv1(const char *name, const char *title);
  virtual      ~AliTRDv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  virtual void  Init();
  virtual void  DrawModule();

protected:
  Int_t        fIdSens1;    // 1st sensitive volume identifier
  Int_t        fIdSens2;    // 2nd sensitive volume identifier
  Int_t        fIdSens3;    // 3rd sensitive volume identifier
            
  ClassDef(AliTRDv1,1)      // Transition Radiation Detector version 1
};

#endif
