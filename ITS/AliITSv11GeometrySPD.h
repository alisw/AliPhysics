#ifndef ALIITSV11GEOMETRYSPD_H
#define ALIITSV11GEOMETRYSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <AliITSv11Geometry.h>
class TGeoVolume;

class AliITSv11GeometrySPD : public AliITSv11Geometry {
  public:
    AliITSv11GeometrySPD(){};
    AliITSv11GeometrySPD(Bool_t debug):AliITSv11Geometry(debug){};
    virtual ~AliITSv11GeometrySPD(){};
    //
    virtual void CarbonFiberSector(TGeoVolume *Moth);
    virtual void HalfStave(TGeoVolume *Moth);
    //
    // Create figures for the documentation of this class
    virtual void CreateFigure0(const Char_t *filepath="",
                               const Char_t *type="gif");
  private:
    void SPDsectorShape(Int_t n,const Double_t *xc,const Double_t *yc,
          const Double_t *r,const Double_t *ths,const Double_t *the,Int_t npr,
                        Int_t &m,Double_t **xp,Double_t **yp);

  private:
    ClassDef(AliITSv11GeometrySPD,1) // ITS v11 Centeral SPD geometry
};

#endif
