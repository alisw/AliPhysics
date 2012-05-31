#ifndef ALIITSV11GEOMBEAMPIPE_H
#define ALIITSV11GEOMBEAMPIPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*************************************************************************
// This class Defines the Geometry of the Beam Pipe
// for the ITS Upgrade using TGeo
//
//  Mario Sitta <sitta@to.infn.it>
//*************************************************************************


/*
  $Id: AliITSv11GeomBeamPipe.h
 */

#include "AliITSv11Geometry.h"
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoXtru.h>

class TGeoVolume;

class AliITSv11GeomBeamPipe : public AliITSv11Geometry {
  public:
    AliITSv11GeomBeamPipe();
    AliITSv11GeomBeamPipe(Double_t rmin, Double_t rmax, Double_t zlen,
			  Int_t debug);
    virtual ~AliITSv11GeomBeamPipe();
    //
    Double_t  GetRmin()     {return fBeamPipeRmin;};
    Double_t  GetRmax()     {return fBeamPipeRmax;};
    Double_t  GetHalfZLen() {return fBeamPipeZlen;};

    void      SetRmin(const Double_t r)     {fBeamPipeRmin = r;};
    void      SetRmax(const Double_t r)     {fBeamPipeRmax = r;};
    void      SetHalfZLen(const Double_t z) {fBeamPipeZlen = z;};

    virtual void CreateBeamPipe(TGeoVolume *moth,
				const TGeoManager *mgr=gGeoManager);

  private:
    AliITSv11GeomBeamPipe(const AliITSv11GeomBeamPipe &source); // copy constructor
    AliITSv11GeomBeamPipe& operator=(const AliITSv11GeomBeamPipe &source); // assignment operator

    Double_t  fBeamPipeRmin;   // Rmin of beam pipe
    Double_t  fBeamPipeRmax;   // Rmax of beam pipe
    Double_t  fBeamPipeZlen;   // Half Z length of beam pipe

  ClassDef(AliITSv11GeomBeamPipe,0) // Beam Pipe geometry for v11 Upgrade
};

#endif
