#ifndef ALIITSV11GEOMETRYSUPPORT_H
#define ALIITSV11GEOMETRYSUPPORT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


// This class Defines the Geometry for the ITS services and support cones
// outside of the ceneteral volume (except for the Ceneteral support 
// cylinders. Other classes define the rest of the ITS. Specificaly the ITS
// The SSD support cone,SSD Support centeral cylinder, SDD support cone,
// The SDD cupport centeral cylinder, the SPD Thermal Sheald, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC. 


/*
  $Id$
 */
#include "AliITSv11Geometry.h"
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoXtru.h>

class TGeoVolume;

class AliITSv11GeometrySupport : public AliITSv11Geometry {
  public:
    AliITSv11GeometrySupport(){};
    AliITSv11GeometrySupport(Int_t debug):AliITSv11Geometry(debug){};
    virtual ~AliITSv11GeometrySupport(){};
    //
    virtual void SPDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void SDDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void SSDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void ServicesCableSupport(TGeoVolume *moth,
                                      TGeoManager *mgr=gGeoManager);
    virtual void ServicesCableSupportSPD(TGeoVolume *moth,
					 TGeoManager *mgr=gGeoManager);
    virtual void ServicesCableSupportSDD(TGeoVolume *moth,
					 TGeoManager *mgr=gGeoManager);
    virtual void ServicesCableSupportSSD(TGeoVolume *moth,
					 TGeoManager *mgr=gGeoManager);


  private:
    void CreateSPDThermalShape(Double_t ina, Double_t inb, Double_t inr,
			       Double_t oua, Double_t oub, Double_t our,
			       Double_t   t, Double_t *x , Double_t *y );
    void CreateSPDOmegaShape(Double_t *xin, Double_t *yin, Double_t  d,
			     Double_t   *x, Double_t *y);
    void FillSPDXtruShape(Double_t a, Double_t  b, Double_t  r,
			  Double_t t, Double_t *x, Double_t *y);
    void PointFromParallelLines(Double_t x1, Double_t y1,
				Double_t x2, Double_t y2, Double_t d,
				Double_t &x, Double_t &y);

    void ReflectPoint(Double_t x1, Double_t y1, Double_t x2, Double_t y2,
		      Double_t x3, Double_t y3, Double_t &x, Double_t &y);

    void  TraySupportsSideA(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SPDCableTraysSideA(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SPDCableTraysSideC(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SDDCableTraysSideA(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SDDCableTraysSideC(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SSDCableTraysSideA(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
    void SSDCableTraysSideC(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);

    TGeoVolumeAssembly* CreateSDDForwardTraySideA(TGeoManager *mgr);

    TGeoCompositeShape* CreateTrayAForwardCover(const Double_t coverLen);
    TGeoCompositeShape* CreateTrayAExternalCover(const Double_t coverLen);
    void CreateTrayACoverHolesShape(const Double_t wide, const Double_t length,
				    const Double_t r10,
				    Double_t *x, Double_t *y);

    TGeoXtru* CreateSDDSSDTraysSideA(const Double_t trayLen,
				     const Double_t trayHi);

    TGeoVolumeAssembly* CreateSDDSSDTraysSideC(const char *trayName,
					       TGeoManager *mgr=gGeoManager);

    ClassDef(AliITSv11GeometrySupport,1) // ITS v11 Support geometry
};

#endif
