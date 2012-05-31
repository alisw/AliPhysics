#ifndef ALIITSMIXTURE_H
#define ALIITSMIXTURE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <TGeoMaterial.h>

class AliITSMixture : public TGeoMixture{
 public:
    AliITSMixture(){};
    AliITSMixture(const char *name,Int_t N,Double_t *w,TObjArray *m,
		  Double_t rho=-1.,Double_t radlen=0.,Double_t intleng=0.);
    virtual ~AliITSMixture(){};
 private:
    ClassDef(AliITSMixture,1) // Extension of TGeoMixture class
}
;
#endif

#ifndef ALIITSGEOCABLE_H
#define ALIITSGEOCABLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <TGeoTube.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TVector3.h>

class AliITSGeoCable : public TObject{
  public:
    AliITSGeoCable();
    AliITSGeoCable(const char *name,const TObjArray *vect,const Double_t Rmin,
                   const Double_t Rmax,const TVector3 ns = TVector3(0.,0.,1.),
                   const TVector3 ne = TVector3(0.,0.,1.));
    virtual ~AliITSGeoCable();
    virtual Double_t GetRmin(){return fRmin;}
    virtual Double_t GetRmax(){return fRmax;}
    virtual TVector3& GetNormStart(){return fNs;}
    virtual TVector3& GetNormEnd(){return fNe;}
    virtual TObjArray* GetArrayOfTubes(){return fTubes;}
    virtual TObjArray* GetArrayOfCombiTrans(){return fTranRot;}
  private:
    Double_t  fRmin; // Minimum radius
    Double_t  fRmax; // Maximum radius
    TVector3  fNs;  // Starting normal vector
    TVector3  fNe;  // Ending normal vector
    TObjArray *fTubes; // Array of Ctub objects
    TObjArray *fTranRot;  // Array of Rotations
    ClassDef(AliITSGeoCable,1) // Extension of TGeoMixture class
        ;
};
#endif
