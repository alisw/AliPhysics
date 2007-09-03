#ifndef ALIITSV11GEOMETRYSPD_H
#define ALIITSV11GEOMETRYSPD_H

/* 
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice.
 */



/*
 * $Id$
 */

//
// Implementation of the SPD v11 central geometry.
// Contains also:
//  - the materials/media used for its volumes;
//  - settings for the related transport parameters (GEANT3 types for the moment).
//
#include <TGeoManager.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TArrayI.h>
#include <AliITSv11Geometry.h>

class TGeoVolume;
class TPolyLine;
class TPolyMarker;


class AliITSv11GeometrySPD : public AliITSv11Geometry
{
public:
        // Default Constructor, should not be used by default
        AliITSv11GeometrySPD() : AliITSv11Geometry(),fSPDsectorShapeName(),
                                 fSPDsectorPoints0(),fSPDsectorPoints1() {};
        // Standard Constructor, set explicitly debug level
        AliITSv11GeometrySPD(Int_t debug) : AliITSv11Geometry(debug),
              fSPDsectorShapeName("ITS SPD Carbon fiber support Sector A0"),
              fSPDsectorPoints0(6),fSPDsectorPoints1() {}; 
        virtual ~AliITSv11GeometrySPD() {}; // Destructor

	/* Settings */

	// define/create materials
	virtual Int_t CreateSPDCentralMaterials(Int_t &medOffset,
						Int_t &matOffset) const;
	// set SPD Central, GEANT3 type, tracking parameters
	virtual void InitSPDCentral(Int_t offset,TVirtualMC *mc=gMC) const;
	
	/* Monitoring */
	
	// creates standard figures for the documentation of this class
	virtual void CreateFigure0(const Char_t *filepath = "",
				   const Char_t *type = "gif",
				   TGeoManager *mgr=gGeoManager);
	
	/* Member functions which create pieces of the geometry */
	
	// a single ladder (= 1 detector + 5 chips)
	TGeoVolume* CreateLadder(Int_t layer, Double_t &width,
				 Double_t &height, Double_t &thickness,
				 TGeoManager *mgr = gGeoManager);
	// the grounding foil (splitted in two components)
	TGeoVolume* CreateGroundingFoilSingle(Bool_t kapLayer,Double_t &len, 
					      Double_t &wid, Double_t &thick,
					      TGeoManager *mgr = gGeoManager);
	TGeoVolume* CreateGroundingFoil(Double_t &thickness,
					TGeoManager *mgr = gGeoManager);
	// the MCM (incomplete: missing the internal chips)
	TGeoVolume* CreateMCMBase(TGeoManager *mgr = gGeoManager) const;
	TGeoVolume* CreateMCMCoverBorder(TGeoManager *mgr = gGeoManager);
	TGeoVolume* CreateMCMCoverTop(TGeoManager *mgr = gGeoManager);
	// the Pixel Bus & extenders
	TGeoVolumeAssembly* CreatePixelBusAndExtensions(Bool_t zpos = kTRUE,
					      TGeoManager *mgr = gGeoManager);
	// the thin part of a stave (grounding + ladders)
	TGeoVolume *CreateStaveBase(Int_t layer, Double_t &width,
				    Double_t &height, Double_t &thickness,
				    TGeoManager *mgr=gGeoManager);
	// the whole stave, including the thick parts (MCM cover, 
	// pixel bus & extensions)
	TGeoVolumeAssembly* CreateStave(Int_t layer, Double_t &thickness,
					TGeoManager *mgr);
	// displacement of staves on the carbon fiber sector
	virtual void StavesInSector(TGeoVolume *moth,
				    TGeoManager *mgr=gGeoManager);
	// the complete Carbon Fiber sector (support + staves)
	virtual void CarbonFiberSector(TGeoVolume *moth,
				       Double_t &xAAtubeCenter0,
				       Double_t &yAAtubeCenter0,
				       TGeoManager *mgr=gGeoManager);
	// the whole SPD barrel
	virtual void SPDSector(TGeoVolume *moth, TGeoManager *mgr=gGeoManager);
	// Returns the mounting location and rotation needed to mount
        // and SPD ladder.
	virtual Bool_t AliITSv11GeometrySPD::GetSectorMountingPoints(
                             Int_t index,Double_t &x0,Double_t &y0,
                             Double_t &x1,Double_t &y1,
			     TGeoManager *mgr=gGeoManager)const;
	// Functions to test aspects of this geometry.
        // Fill TPolylines with crossections of the SPD Carbon
	// fiber sectors.
	Bool_t Make2DcrossSections(TPolyLine &a0,TPolyLine &a1,
				   TPolyLine &b0,TPolyLine &b1,
				   TPolyMarker &p)const;

private:
	
	// NOTE:
	// all of the member functions which define a component of the 
	// final SPD will need to be defined as private once the design 
	// is fixed and does not
	// need any longer to be checked and debugged.
	
	// Computes shape of the SPD Sector given specific inputs 
	// (internal use only)
	void SPDsectorShape(Int_t n, const Double_t *xc, const Double_t *yc,
			    const Double_t *r,
	                    const Double_t *ths, const Double_t *the, Int_t npr,
	                    Int_t &m, Double_t **xp, Double_t **yp);
	
	// computes a point to a line parallel to a given direction
	// and with a fixed distance from it (internal use only)
	void ParallelPosition(Double_t dist1, Double_t dist2,
			      Double_t phi, Double_t &x, Double_t &y);

	// Some internal data
	TString fSPDsectorShapeName; // SPD Sector Shape name
	TArrayI fSPDsectorPoints0;    // Array of sector points index
                                     // used to mount SPD ladder to.
	TArrayI fSPDsectorPoints1;    // Array of sector points index
                                     // used to mount SPD ladder to.
	
	ClassDef(AliITSv11GeometrySPD,1) // ITS v11 Centeral SPD geometry
};

#endif

