#ifndef ALITRD_H
#define ALITRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////
 
#include "AliRun.h"
#include "AliDetector.h"
#include "AliTRDhit.h" 
#include "AliDigit.h"

#include "AliTRDconst.h"
#include "AliTRDgeometry.h"

//_____________________________________________________________________________
class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  AliTRD(AliTRD &trd);
  virtual           ~AliTRD();
  virtual void       AddHit(Int_t track, Int_t *det, Float_t *hits);
  virtual void       AddDigit(Int_t *digits, Int_t *amp);    
  virtual void       AddRecPoint(Float_t *pos, Int_t *digits
                               , Int_t det, Float_t amp);
  virtual void       BuildGeometry();
  virtual void       Copy(AliTRD &trd);
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual void       DrawModule();
  Int_t              DistancetoPrimitive(Int_t px, Int_t py);
  TObjArray         *RecPoints()           { return fRecPoints;   };
  virtual void       Init();
  virtual Int_t      IsVersion() const = 0;
  virtual void       MakeBranch(Option_t* option);     
  virtual void       ResetRecPoints();
  virtual void       StepManager() = 0; 
  virtual void       SetTreeAddress();

  virtual void       SetGasMix(Int_t imix = 0);
  virtual void       SetHits()             {};
  virtual void       SetPHOShole()         { fGeometry->SetPHOShole(); };
  virtual void       SetRICHhole()         { fGeometry->SetRICHhole(); };

  AliTRDgeometry    *GetGeometry()         { return fGeometry; };

  virtual void       SetSensChamber(Int_t ichamber)              = 0;
  virtual void       SetSensPlane(Int_t iplane)                  = 0;
  virtual void       SetSensSector(Int_t isector)                = 0;
  virtual void       SetSensSector(Int_t isector, Int_t nsector) = 0;

  virtual Int_t      GetSensChamber()     = 0;
  virtual Int_t      GetSensPlane()       = 0;
  virtual Int_t      GetSensSector()      = 0;
  virtual Int_t      GetSensSectorRange() = 0; 

  inline  AliTRD    &operator=(AliTRD &trd);

 protected:

  Int_t              fGasMix;            //  Gas mixture. 0: Xe/Isobutane 1: Xe/CO2

  AliTRDgeometry    *fGeometry;          //  The TRD geometry

  TObjArray         *fRecPoints;         //  Array of reconstructed points
  Int_t              fNRecPoints;        //! Number of reconstructed points

  ClassDef(AliTRD,1)                     //  Transition Radiation Detector base class

};

//_____________________________________________________________________________
AliTRD &AliTRD::operator=(AliTRD &trd)
{
  //
  // Assignment operator
  //

  if (this != &trd) trd.Copy(*this);
  return *this;

}

#endif
