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
  virtual void  SetSensPlane(Int_t iplane = 0);
  virtual void  SetSensChamber(Int_t ichamber = 0);
  virtual void  SetSensSector(Int_t isector = 0);
  virtual void  Init();

protected:
  Int_t        fIdSens;          // Sensitive volume identifier

  Int_t        fIdSpace1;        // Spaceframe volume identifier
  Int_t        fIdSpace2;        // 
  Int_t        fIdSpace3;        // 

  Int_t        fIdChamber1;      // Driftchamber volume identifier
  Int_t        fIdChamber2;      // 
  Int_t        fIdChamber3;      // 

  Int_t        fSensSelect;      // Switch to select only parts of the detector
  Int_t        fSensPlane;       // Sensitive detector plane
  Int_t        fSensChamber;     // Sensitive detector chamber
  Int_t        fSensSector;      // Sensitive detector sector

private:
  virtual Double_t BetheBloch(Double_t bg);

  TF1         *fDeltaE;          // Energy distribution of the delta-electrons
  
  ClassDef(AliTRDv2,1)           // Transition Radiation Detector version 2

};

#endif
