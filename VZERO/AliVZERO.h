#ifndef VZERO_H
#define VZERO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO    //
//////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
#include "TNamed.h"
#include "TTree.h"

class TDirectory;
R__EXTERN TDirectory *  gDirectory;
 
 
class AliVZERO : public AliDetector {
 
public:

  AliVZERO() {}
  AliVZERO(const char *name, const char *title);
  virtual       ~AliVZERO() {}
  virtual void   AddDigit( Int_t* tracks, Int_t* digits) = 0;
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const = 0;
  virtual void   Init();
  virtual void   DrawModule() {};
  virtual void   StepManager() {};
  virtual inline  void   SetThickness(Float_t thick)  {fThickness = thick;};
  virtual inline  void   SetThickness1(Float_t thick) {fThickness1 = thick;};
// Set Stepping Parameters
  virtual void   SetMaxStepQua(Float_t p1);
  virtual void   SetMaxStepAlu(Float_t p1);
  virtual void   SetMaxDestepQua(Float_t p1);
  virtual void   SetMaxDestepAlu(Float_t p1);

   Int_t         fNCerenkovs;    //  Number of cerenkovs (detected on photocathode)
   Int_t         fNGCerenkovs;   //  Number of cerenkovs (generated)
   Int_t         fMulti; // Multiplicity of charged particles

   Float_t      fThickness;
   Float_t      fThickness1;

protected:

   Int_t fIdSens1;
  
// Stepping Parameters
   Float_t fMaxStepQua;      // Maximum step size inside the quartz volumes
   Float_t fMaxStepAlu;      // Maximum step size inside the  aluminum volumes
   Float_t fMaxDestepQua;    // Maximum relative energy loss in quartz
   Float_t fMaxDestepAlu;    // Maximum relative energy loss in aluminum
  
  ClassDef(AliVZERO,1)  //Class for the VZERO detector
};

//____________________________________________________________

#endif
