#ifndef ALITRD_H
#define ALITRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////

#include "AliDetector.h"

class TFile;

class AliRun;
class AliDigit;

class AliTRDhit;
class AliTRDsim;
class AliTRDgeometry;

//_____________________________________________________________________________
class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  AliTRD(const AliTRD &trd);
  virtual           ~AliTRD();

          AliTRD    &operator=(const AliTRD &trd);

  virtual void       AddHit(Int_t track, Int_t det, Float_t *hits, Int_t q);
  virtual void       AddDigit(Int_t *digits, Int_t *amp);    
  virtual void       AddRecPoint(Float_t *pos, Int_t *digits
                               , Int_t det, Float_t amp, Int_t *tracks);
  virtual void       BuildGeometry();
  virtual void       Copy(TObject &trd);
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual void       DrawModule();
  Int_t              DistancetoPrimitive(Int_t px, Int_t py);
  virtual void       LoadPoints(Int_t track);    
  TObjArray         *RecPoints() const              { return fRecPoints;   };
  virtual void       Init();
  virtual Int_t      IsVersion() const = 0;
  virtual void       MakeBranch(Option_t* option, char *file=0);     
  virtual void       ResetRecPoints();
  virtual void       StepManager() = 0; 
  virtual void       SetTreeAddress();

  virtual void       SetGasMix(Int_t imix = 0);
  virtual void       SetHits()             {};
  virtual void       SetPHOShole();
  virtual void       SetRICHhole();
  virtual void       SetDrawTR(Int_t idraw = 1)     { fDrawTR      = idraw; };
  virtual void       SetDisplayType(Int_t type = 0) { fDisplayType = type;  };

  AliTRDgeometry    *GetGeometry() const            { return fGeometry; };

  virtual void       SetSensChamber(Int_t ichamber)              = 0;
  virtual void       SetSensPlane(Int_t iplane)                  = 0;
  virtual void       SetSensSector(Int_t isector)                = 0;
  virtual void       SetSensSector(Int_t isector, Int_t nsector) = 0;

  virtual Int_t      GetSensChamber() const     = 0;
  virtual Int_t      GetSensPlane() const       = 0;
  virtual Int_t      GetSensSector() const      = 0;
  virtual Int_t      GetSensSectorRange() const = 0; 

  virtual void       SDigits2Digits();

  virtual AliTRDsim *CreateTR()     = 0;
  virtual AliTRDsim *GetTR() const  = 0;

 protected:

  Int_t              fGasMix;            //  Gas mixture. 0: Xe/Isobutane 1: Xe/CO2

  AliTRDgeometry    *fGeometry;          //  The TRD geometry

  TObjArray         *fRecPoints;         //  Array of reconstructed points
  Int_t              fNRecPoints;        //! Number of reconstructed points

  Float_t            fGasDensity;        //  The density of the drift gas
  Float_t            fFoilDensity;       //  The density of the entrance window foil

  Int_t              fDrawTR;            //  Switches marking the TR photons in the display
  Int_t              fDisplayType;       //  Display type (0: normal, 1: detailed) 

  ClassDef(AliTRD,1)                     //  Transition Radiation Detector base class

};

#endif
