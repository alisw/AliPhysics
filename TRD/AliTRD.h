#ifndef ALITRD_H
#define ALITRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TRD     //
////////////////////////////////////////////////

#include <TLorentzVector.h>
#include "AliDetector.h"
#include <TVirtualMC.h>

class TFile;
class TLorentzVector;

class AliRun;
class AliDigit;

class AliTRDhit;
class AliTRDsim;
class AliTRDgeometry;
class AliTRDtrackHits;

//_____________________________________________________________________________
class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  AliTRD(const AliTRD &trd);
  virtual           ~AliTRD();

          AliTRD    &operator=(const AliTRD &trd);

  virtual void       AddHit(Int_t track, Int_t det, Float_t *hits, Int_t q, Bool_t inDrift); 
  virtual void       AddCluster(Float_t *pos, Int_t det, Float_t amp, Int_t *tracks
                              , Float_t *sig, Int_t iType);
  virtual void       BuildGeometry();
  virtual void       Copy(TObject &trd);
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual void       DrawModule() const;
  Int_t              DistancetoPrimitive(Int_t px, Int_t py) const;
  virtual void       LoadPoints(Int_t track);    
  TObjArray         *RecPoints() const              { return fRecPoints;   };
  virtual void       Init();
  virtual Int_t      IsVersion() const = 0;
  virtual void       MakeBranch(Option_t* option);
  virtual void       ResetDigits();     
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

  virtual void       SetHitTypeStandard()           { fHitType = 1; };
  virtual void       SetHitTypeCompressed()         { fHitType = 2; };
  virtual void       SetHitTypeBoth()               { fHitType = 3; };

  virtual Int_t      GetSensChamber() const     = 0;
  virtual Int_t      GetSensPlane() const       = 0;
  virtual Int_t      GetSensSector() const      = 0;
  virtual Int_t      GetSensSectorRange() const = 0; 
 
  virtual void       Hits2Digits();
  virtual void       Hits2SDigits();
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const; 
  virtual void       SDigits2Digits();

  virtual AliTRDsim *CreateTR()     = 0;
  virtual AliTRDsim *GetTR() const  = 0;

  virtual void       MakeBranch2(Option_t *opt=" ", const char *file=0); 
  virtual void       SetTreeAddress2();
  virtual void       AddHit2(Int_t track, Int_t det, Float_t *hits, Int_t q, Bool_t inDrift);
  virtual void       ResetHits();

  virtual AliHit*    FirstHit(Int_t track);
  virtual AliHit*    NextHit();
  virtual AliHit*    FirstHit2(Int_t track);
  virtual AliHit*    NextHit2();

  virtual void       FinishPrimary();
  virtual void       RemapTrackHitIDs(Int_t *map);

 protected:

  Int_t                fGasMix;             //  Gas mixture. 0: Xe/Isobutane 1: Xe/CO2

  AliTRDgeometry      *fGeometry;           //  The TRD geometry

  TObjArray           *fRecPoints;          //  Array of reconstructed points / cluster
  Int_t                fNRecPoints;         //! Number of reconstructed points / cluster

  Float_t              fGasDensity;         //  The density of the drift gas
  Float_t              fFoilDensity;        //  The density of the entrance window foil

  Int_t                fDrawTR;             //  Switches marking the TR photons in the display
  Int_t                fDisplayType;        //  Display type (0: normal, 1: detailed) 


  AliTRDtrackHits     *fTrackHits;          //! Compressed hits for a given track 
  Int_t                fHitType;            //  Switch for the different hit data structures

  ClassDef(AliTRD,6)                        //  Transition Radiation Detector base class

};

#endif
