#ifndef ALIZDC_H
#define ALIZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and classes for set ZDC           //
////////////////////////////////////////////////

#include "AliDetector.h"

class AliZDCCalibData;
 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual       ~AliZDC();
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t IsVersion() const =0;
  virtual Float_t ZMin() const;	// Minimum overall dimension of the ZDC
  virtual Float_t ZMax() const;	// Maximum overall dimension of the ZDC
  virtual void  SetTreeAddress();
  virtual void  MakeBranch(Option_t* opt, const char *file=0);
  virtual void  Hits2SDigits();
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  Digits2Raw();
  virtual void  StepManager() {}
    
  // Switching off the shower development in ZDCs
  void  NoShower(){fNoShower=1;}
  void  Shower()  {fNoShower=0;}


//Calibration methods (by Alberto Colla)
  void          CreateCalibData();
  void          WriteCalibData(Int_t option=TObject::kOverwrite);
  void          LoadCalibData();
  void          SetCalibData(AliZDCCalibData* data) {fCalibData = data;}
  AliZDCCalibData* GetCalibData() const  {return fCalibData;}
//Calibration methods (by Alberto Colla)


protected:

  Int_t        fNoShower;	// Flag to switch off the shower	

//Calibration methods (by Alberto Colla)
  AliZDCCalibData* fCalibData;	// Calibration data for ZDC
//Calibration methods (by Alberto Colla)
  
  ClassDef(AliZDC,4)  	// Zero Degree Calorimeter base class
};
 
R__EXTERN  AliZDC *gZDC;

#endif
