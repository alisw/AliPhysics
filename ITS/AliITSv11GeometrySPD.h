#ifndef ALIITSV11GEOMETRYSPD_H
#define ALIITSV11GEOMETRYSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Class which defines the SPD v11 centeral geometry, defines the
  materials/media used for this geometry and sets the related transport
  parameters (GEANT3 types for the moment.
 */


#include <TGeoManager.h>
#include <TVirtualMC.h>
#include <AliITSv11Geometry.h>

class TGeoVolume;

class AliITSv11GeometrySPD : public AliITSv11Geometry {
  public:
    // Default constructor
    AliITSv11GeometrySPD():AliITSv11Geometry(){};
    // Default constructor with debuging level switch
    AliITSv11GeometrySPD(Bool_t debug):AliITSv11Geometry(debug){};
    // Destructor
    virtual ~AliITSv11GeometrySPD(){};
    // Creates SPD Sector geometry
    virtual void SPDSector(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    // Creates SPD Carbon Fiber Sector only
    virtual void CarbonFiberSector(TGeoVolume *moth,Double_t &xAAtubeCenter0,
				   Double_t &yAAtubeCenter0,
				   TGeoManager *mgr=gGeoManager);
    // creates SPD Half Stave
    virtual void HalfStave(TGeoVolume *moth,Double_t &thicknessAA,
			   TGeoManager *mgr=gGeoManager);
    //
    // Creates standard figures for the documentation of this class
    virtual void CreateFigure0(const Char_t *filepath="",
                               const Char_t *type="gif",
			       TGeoManager *mgr=gGeoManager);
    // Defnes/creates SPD Centeral detector materials
    virtual Int_t CreateSPDCenteralMaterials(Int_t &Medffset,Int_t &matOffset);
    // Sets SPD Centeral, GEANT3 type, Tracking parameters
    virtual void InitSPDCenteral(Int_t offset,TVirtualMC *mc=gMC);
    //
  private:
    // Computes shape of the SPD Sector given specific inputs Internal use only
    void SPDsectorShape(Int_t n,const Double_t *xc,const Double_t *yc,
          const Double_t *r,const Double_t *ths,const Double_t *the,Int_t npr,
                        Int_t &m,Double_t **xp,Double_t **yp);

  private:
    ClassDef(AliITSv11GeometrySPD,1) // ITS v11 Centeral SPD geometry
};

#endif
