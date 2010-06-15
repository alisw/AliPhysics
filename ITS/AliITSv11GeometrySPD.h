#ifndef ALIITSV11GEOMETRYSPD_H
#define ALIITSV11GEOMETRYSPD_H

/*
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice.
 */

 // Implementation of the SPD v11 central geometry.
 // Contains also:
 //  - the materials/media used for its volumes;
 //  - settings for the related transport parameters
 //   (GEANT3 types for the moment).
 //

/*
 * $Id$
 */

#include <TGeoManager.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TArrayI.h>
#include <TPolyLine.h>
#include <TPolyMarker.h>
#include <AliITSv11Geometry.h>

class TGeoVolume;
class TGeoCompositeShape;

class AliITSv11GeometrySPD : public AliITSv11Geometry
{
 public:

    // Default constructor
    AliITSv11GeometrySPD(/*Double_t gap = 0.0075*/);
    // Standard Constructor
    AliITSv11GeometrySPD(Int_t debug/*, Double_t gap = 0.0075*/);
    // Copy constructor
    AliITSv11GeometrySPD(const AliITSv11GeometrySPD &s);
    // Assignment operator
    AliITSv11GeometrySPD& operator=(const AliITSv11GeometrySPD &s);
    // Destructor
    virtual ~AliITSv11GeometrySPD() {};

    /* Settings */

    // define/create materials
    virtual Int_t CreateSPDCentralMaterials(Int_t &medOffset,
                                            Int_t &matOffset) const;
    /* Monitoring */

    // creates standard figures for the documentation of this class
    virtual void CreateFigure0(const Char_t *path = "",
                               const Char_t *type = "gif",
                               TGeoManager *mgr = gGeoManager) const;
    // fill TPolylines with crossections of the SPD Carbon fiber sectors.
    Bool_t Make2DCrossSections(TPolyLine &a0, TPolyLine &a1, TPolyLine &b0,
                               TPolyLine &b1, TPolyMarker &p) const;

    /* Services */

    // get names
    virtual const char *GetSenstiveVolumeName1() const
        {return "ITSSPDlay1-sensor";}
    virtual const char *GetSenstiveVolumeName2() const
        {return "ITSSPDlay2-sensor";}
    virtual const char *GetSenstiveVolumeName(Int_t lay) const
        {return (lay==1) ? GetSenstiveVolumeName1():GetSenstiveVolumeName2();}
    // get medium
    virtual TGeoMedium* GetMedium(const char* mediumName,
                                  TGeoManager *mgr = gGeoManager) const;
    // retrieve the mounting location and rotation needed to mount an SPD stave
    virtual Bool_t GetSectorMountingPoints(Int_t index, Double_t &x0,
                               Double_t &y0, Double_t &x1, Double_t &y1) const;
    // displace the staves on the carbon fiber sector
    virtual void StavesInSector(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    // (debug purposes) define which staves to put in the sector
    virtual void SetAddStave(Bool_t *mask);
    // print class in ascii form to stream
    virtual void PrintAscii(ostream *os) const;
    // read in class in ascii form from stream
    virtual void ReadAscii(istream *is);

    /* Parts of the geometry */

    // a single ladder (= 1 detector + 5 chips)
    virtual TGeoVolume* CreateLadder(Int_t layer, TArrayD &sizes,
                                     TGeoManager *mgr = gGeoManager) const;
    // a clip on the central ladders
    virtual TGeoVolume* CreateClip(TArrayD &sizes,Bool_t isDummy,
                                   TGeoManager *mgr = gGeoManager) const;
    // the grounding foil (splitted in many components)
    //virtual TGeoVolumeAssembly* CreateGroundingFoilSingle(Int_t type,
    //                 TArrayD &sizes, TGeoManager *mgr = gGeoManager) const;
    virtual  TGeoCompositeShape* CreateGroundingFoilShape(Int_t itype,
        Double_t &length,Double_t &width,Double_t thickness,TArrayD &sizes);
    virtual TGeoVolumeAssembly* CreateGroundingFoil(Bool_t isRight, TArrayD &sizes,
                                        TGeoManager *mgr = gGeoManager);
    // the MCM (thin part + thick part with chips inside)
    virtual TGeoVolumeAssembly* CreateMCM(Bool_t isRight, TArrayD &sizes,
                                       TGeoManager *mgr = gGeoManager) const;
    // the pixel bus (flat part + pt1000s + large capacitors/resistors)
    virtual TGeoVolumeAssembly* CreatePixelBus(Bool_t isRight, Int_t layer, TArrayD &sizes,
                                        TGeoManager *mgr = gGeoManager) const;
    // the extender complicated geometry
    virtual TGeoVolume* CreateExtender(const Double_t *params,
                              const TGeoMedium *medium, TArrayD &sizes) const;
    // the Pixel Bus & extenders (old method which will be removed)
    virtual TGeoVolumeAssembly* CreatePixelBusAndExtensions(Bool_t zpos=kTRUE,
                                        TGeoManager *mgr = gGeoManager) const;

    virtual TList* CreateConeModule(const Double_t angle,
				    TGeoManager *mgr = gGeoManager) const;
    virtual void CreateCones(TGeoVolume *moth) const;
    // a half-stave (put together ladders + MCM + bus, and add clips
    // if requested)
    virtual TGeoVolumeAssembly* CreateHalfStave(Bool_t isRight, Int_t layer,
                 Int_t idxCentral, Int_t idxSide,TArrayD &sizes/*,
              Bool_t addClips = kFALSE*/, TGeoManager *mgr = gGeoManager);
    // the whole stave (2 half-staves of different orientation)
    virtual TGeoVolumeAssembly* CreateStave(Int_t layer, TArrayD &sizes,
             /*Bool_t addClips = kFALSE,*/TGeoManager *mgr = gGeoManager);
    // the complete Carbon Fiber sector (support + staves)
    virtual void CarbonFiberSector(TGeoVolume *moth, Double_t &xAAtubeCenter0,
                     Double_t &yAAtubeCenter0, TGeoManager *mgr = gGeoManager);
    // the whole SPD barrel (the 10 sectors at once)
    virtual void SPDSector(TGeoVolume *moth, TGeoManager *mgr = gGeoManager);
    // Returns the location of the SPD cooling tube ends. RB26 (muon absober
    // side) and RB24 (open side). Staves number 0,1 inner Staves, 2-5 outer
    // staves. Sectors numbers 0-9.
    virtual void GetSPDCoolingTubeRB26(Int_t sector,Int_t stave,
                                 Double_t &x,Double_t &y,Double_t &z)const{
                            x = fTubeEndSector[sector][1][stave][0];
                            y = fTubeEndSector[sector][1][stave][1];
                            z = fTubeEndSector[sector][1][stave][2];return;};
    virtual void GetSPDCoolingTubeRB24(Int_t sector,Int_t stave,
                                 Double_t &x,Double_t &y,Double_t &z)const{
                            x = fTubeEndSector[sector][0][stave][0];
                            y = fTubeEndSector[sector][0][stave][1];
                            z = fTubeEndSector[sector][0][stave][2];return;};
 private:
    // NOTE:
    // all of the member functions which define a component of the final SPD
    // will need to be defined as private once the design is fixed and
    // does not need any longer to be checked and debugged.

    /* Service methods for internal use only */

    // compute shape of the SPD Sector given specific inputs
    void SPDsectorShape(Int_t n,const Double_t *xc, const Double_t *yc,
                        const Double_t *r,const Double_t *ths,
                        const Double_t *the, Int_t npr,Int_t &m,
                        Double_t **xp, Double_t **yp) const;
    // compute a point o a line parallel to a given direction
    // and with a fixed distance from it
    void ParallelPosition(Double_t dist1, Double_t dist2, Double_t phi,
                          Double_t &x, Double_t &y) const;
    // comutes the radial translation of a sector to give the
    // proper distance between SPD detectors and the beam pipe.
    Double_t GetSPDSectorTranslation(Double_t x0,Double_t y0,Double_t x1,
                                     Double_t y1,Double_t r)const;
    Bool_t CFHolePoints(Double_t s,Double_t r1,Double_t r2,Double_t l,
                        Double_t &x,Double_t &y)const;

    /* Data members */

    static const Double_t fgkGapLadder;// thicknes of the empty (air) gap left
                               // between the ladder and the grounding
                               // foil for alignment
    static const Double_t fgkGapHalfStave;//thickness of the empty (air) gap
                                          // left between HS and Carbon Suport
    Bool_t  fAddStave[6];      // [DEBUG] must be TRUE for all staves
	                       // which will be mounted in the sector
                               // (used to check overlaps)
    TArrayD fSPDsectorX0;      // X of first edge of sector plane for stave
    TArrayD fSPDsectorY0;      // Y of first edge of sector plane for stave
    TArrayD fSPDsectorX1;      // X of second edge of sector plane for stave
    TArrayD fSPDsectorY1;      // Y of second edge of sector plane for stave
    //
    Double_t fTubeEndSector[10][2][6][3]; // Location of tube end in sector
    /* ROOT dictionary */

    ClassDef(AliITSv11GeometrySPD,2) // ITS v11 Central SPD geometry
};

// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os, const AliITSv11GeometrySPD &s);
istream &operator>>(istream &is, AliITSv11GeometrySPD &s);

#endif
