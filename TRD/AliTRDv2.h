#ifndef TRDv2_H
#define TRDv2_H
////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 2    //
////////////////////////////////////////////////////////
 
#include <TF1.h> 
#include "AliTRD.h"

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);

class AliTRDv2 : public AliTRD {

public:
  AliTRDv2() {}
  AliTRDv2(const char *name, const char *title);
  virtual      ~AliTRDv2();
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 2;}
  virtual void  StepManager();
  virtual void  Init();
  virtual void  DrawModule();

protected:
  Int_t        fIdSensI[ncham];  // Sensitive volume identifier (inner chambers)
  Int_t        fIdSensN[ncham];  // Sensitive volume identifier (neighbouring chambers)
  Int_t        fIdSensO[ncham];  // Sensitive volume identifier (outer chambers)

private:
  virtual Double_t BetheBloch(Double_t bg);

  TF1         *fDeltaE;          // Energy distribution of the delta-electrons
  
  ClassDef(AliTRDv2,1)           // Transition Radiation Detector version 2
};

#endif
