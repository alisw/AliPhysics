#ifndef ALIITSUV2LAYER_H
#define ALIITSUV2LAYER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*************************************************************************
// This class Defines the Geometry for the ITS Upgrade using TGeo
// This is a work class used to study different configurations
// during the development of the new ITS structure.
//
//  Mario Sitta <sitta@to.infn.it>
//*************************************************************************


/*
  $Id: AliITSUv2Layer.h
 */

#include "AliITSUv2.h"
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoXtru.h>

class TGeoVolume;

class AliITSUv2Layer : public TObject {
  public:
  enum {kStave,kHalfStave,kModule,kChip,kNHLevels};

  public:
    AliITSUv2Layer();
    AliITSUv2Layer(Int_t debug);
    AliITSUv2Layer(Int_t lay, Int_t debug);
    AliITSUv2Layer(Int_t lay, Bool_t turbo, Int_t debug);
    AliITSUv2Layer(const AliITSUv2Layer &source);
    AliITSUv2Layer& operator=(const AliITSUv2Layer &source);
    virtual ~AliITSUv2Layer();
    //
    Bool_t    IsTurbo() const {return fIsTurbo;};

    Double_t  GetChipThick()   const {return fChipThick;};
    Double_t  GetStaveTilt()   const {return fStaveTilt;};
    Double_t  GetStaveWidth()  const {return fStaveWidth;};
    Double_t  GetSensorThick() const {return fSensorThick;};
    Double_t  GetNStaves()     const {return fNStaves;};
    Double_t  GetNChips()      const {return fNChips;};
    Double_t  GetRadius()      const {return fLayRadius;};
    Double_t  GetPhi0()        const {return fPhi0;};
    Double_t  GetZLength()     const {return fZLength;};
    Int_t     GetChipType()    const {return fChipTypeID;}
    //
    Int_t     GetNStavesPerParent()     const {return fHierarchy[kStave];}
    Int_t     GetNHalfStavesPerParent() const {return fHierarchy[kHalfStave];}
    Int_t     GetNModulesPerParent()    const {return fHierarchy[kModule];}
    Int_t     GetNChipsPerParent()      const {return fHierarchy[kChip];}
    //
    Int_t     GetBuildLevel()           const {return fBuildLevel;}
    AliITSUv2::AliITSUModel_t GetStaveModel() const {return fStaveModel;}
    //
    void      SetChipThick(Double_t t)       {fChipThick = t;};
    void      SetStaveTilt(Double_t t);
    void      SetStaveWidth(Double_t w);
    void      SetSensorThick(Double_t t)     {fSensorThick = t;};
    void      SetNStaves(Int_t n)            {fHierarchy[kStave] = fNStaves = n;};
    void      SetNUnits(Int_t u);
    void      SetRadius(Double_t r)          {fLayRadius = r;};
    void      SetPhi0(Double_t phi)          {fPhi0 = phi;}
    void      SetZLength(Double_t z)         {fZLength   = z;};
    void      SetChipType(Int_t tp)          {fChipTypeID = tp;}
    void      SetBuildLevel(Int_t buildLevel){fBuildLevel=buildLevel;}
    void      SetStaveModel(AliITSUv2::AliITSUModel_t model) {fStaveModel=model;}
    virtual void CreateLayer(TGeoVolume *moth);

    virtual void CreateBarrelLayer(TGeoVolume* moth);	// #pnamwong
	static void CreateBarrel(TGeoVolume* moth);		// #pnamwong

    // Define Trig functions for use with degrees (standerd TGeo angles).
    // Sine function
    Double_t SinD(Double_t deg)const{return TMath::Sin(deg*TMath::DegToRad());}
    // Cosine function
    Double_t CosD(Double_t deg)const{return TMath::Cos(deg*TMath::DegToRad());}
    // Tangent function
    Double_t TanD(Double_t deg)const{return TMath::Tan(deg*TMath::DegToRad());}	
	
 protected:

    // Units, Convert from k?? to cm,degree,GeV,seconds,
    static const Double_t fgkmicron; // Convert micron to TGeom's cm.
    static const Double_t fgkmm; // Convert mm to TGeom's cm.
    static const Double_t fgkcm; // Convert cm to TGeom's cm.

  private:
    void CreateLayerTurbo(TGeoVolume *moth);

    TGeoVolume* CreateStave(const TGeoManager *mgr=gGeoManager);
    //TGeoVolume* CreateChip(Double_t x, Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateModuleInnerB(Double_t x,Double_t y, Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateChipInnerB(Double_t x,Double_t y, Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateModuleOuterB(const TGeoManager *mgr=gGeoManager);


    TGeoVolume* CreateStaveInnerB(Double_t x, Double_t y, Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveStructInnerB(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerBDummy(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager) const;
    TGeoVolume* CreateStaveModelInnerB0(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerB1(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerB21(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerB22(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerB3(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelInnerB4(Double_t x,Double_t z, const TGeoManager *mgr=gGeoManager);
    void        CreateIBConnectors(const TGeoManager *mgr=gGeoManager);
    void        CreateIBConnectorsASide(const TGeoManager *mgr=gGeoManager);
    void        CreateIBConnectorsCSide(const TGeoManager *mgr=gGeoManager);


    TGeoVolume* CreateStaveOuterB(const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelOuterBDummy(const TGeoManager *mgr=gGeoManager) const;
    TGeoVolume* CreateStaveModelOuterB0(const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateStaveModelOuterB12(const TGeoManager *mgr=gGeoManager);
    void        CreateOBColdPlateConnectors();
    void        CreateOBColdPlateConnectorsASide();
    void        CreateOBColdPlateConnectorsCSide();
    TGeoVolume* CreateSpaceFrameOuterB(const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateSpaceFrameOuterBDummy(const TGeoManager *mgr=gGeoManager) const;
    TGeoVolume* CreateSpaceFrameOuterB1(const TGeoManager *mgr=gGeoManager);
    void        CreateOBSpaceFrameObjects(const TGeoManager *mgr=gGeoManager);

    TGeoXtru* CreateStaveSide(const char *name,
			       Double_t dz, Double_t alpha, Double_t beta,
			       Double_t L, Double_t H, Bool_t top);
    TGeoCombiTrans* CreateCombiTrans( const char *name,
				      Double_t dy, Double_t dz, Double_t dphi,
				      Bool_t planeSym=kFALSE);
    void AddTranslationToCombiTrans( TGeoCombiTrans* ct,
				     Double_t dx=0, Double_t dy=0,
				     Double_t dz=0) const;

    // ----------- Barrel -------------
    // #pnamwong
    // CreateBarrelLayer()
    TGeoVolume* CreateInnerBEWheelA2(const Color_t, const TGeoMedium*);

    TGeoVolume* CreateInnerBEWheelA3(const Color_t, const TGeoMedium*);
    TGeoVolume* CreateInnerBEWheelC3(const Color_t, const TGeoMedium*);

   	TGeoVolume* CreateInnerBEWheelC(const TGeoMedium*);
	TGeoVolume* CreateOuterBEWheelC(const TGeoManager *mgr=gGeoManager);
	TGeoVolume* CreateOuterBEWheelA(const TGeoManager *mgr=gGeoManager);	
	TGeoVolume* CreateInnerBHandleA(const char *name, const Double_t pr, Double_t cr, const TGeoMedium*);

	TGeoVolume* Create_ALC_0334_xxA(const TGeoMedium*);			// Inner Barrel End Wheel Handle (in each)
	TGeoVolume* Create_ALC_0334_xxC(const TGeoMedium*);			// Inner Barrel End Wheel Handle (in each)

	TGeoVolume* Create_ALC_0336_xxA_ins(const TGeoMedium*);					// Inner Barrel End Wheel Plate inside the warp volume (in half stave)
	TGeoVolume* Create_ALC_0336_xxA_ous(const Int_t, const TGeoMedium*);	// Inner Barrel End Wheel Plate outside the warp volume (in half stave)
	TGeoVolume* Create_ALC_0337_xxC_ins(const TGeoMedium*);					// Inner Barrel End Wheel Plate inside the warp volume (in half stave)
	TGeoVolume* Create_ALC_0337_xxC_ous(const Int_t, const TGeoMedium*);	// Inner Barrel End Wheel Plate inside the warp volume (in half stave)

	TGeoVolume* CreateInnerDCDC3(const TGeoMedium*);
	TGeoVolume* CreateInnerDCCNT3(const TGeoMedium*);
	
	TGeoVolume* CreateInnerBHandleC(const TGeoMedium*);

    TGeoVolume* CreateInnerBEWheelA0(const TGeoMedium*);
    TGeoVolume* CreateInnerBEWheelA1(const TGeoMedium*);
    
    TGeoVolume* CreateSamtecCable(const TGeoMedium*);
    
    // CreateBarrel()
    static TGeoVolume* CreateInnerBSupporterA(const TGeoManager *mgr=gGeoManager);
	static TGeoVolume* CreateInnerBSupporterC(const TGeoManager *mgr=gGeoManager);
    
    static TGeoVolume* CreateInnerBCShell(const TGeoManager *mgr=gGeoManager);
	static TGeoVolume* CreateOuterBCShell(const TGeoManager *mgr=gGeoManager);	

	static TGeoVolume* CreateInnerBServiceB(const TGeoManager *mgr=gGeoManager);
	static TGeoVolume* CreateInnerBServiceB2(const TGeoManager *mgr=gGeoManager);
	static TGeoVolume* CreateOuterBServiceB(const TGeoManager *mgr=gGeoManager);
	
	// Create DC Powersupply()
	TGeoVolume* CreateDCSet(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	
	TGeoVolume* CreateDCUnit(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	TGeoVolume* CreateDCBase(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	TGeoVolume* CreateRoundP(const char *nid = "", const TGeoManager *mgr=gGeoManager);

	TGeoVolume* CreateDCPSLayer(const char *nid = "", const TGeoManager *mgr=gGeoManager);

	TGeoVolume* CreateDCPS1(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	TGeoVolume* CreateDCBox1(const char *nid = "", const TGeoManager *mgr=gGeoManager);

	TGeoVolume* CreateDCPS2(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	
	TGeoVolume* CreateDCPipe(const char *nid = "", const TGeoManager *mgr=gGeoManager);
	TGeoVolume* CreateDCPipes(const char *nid = "", const TGeoManager *mgr=gGeoManager);

	TGeoVolume* CreateFPGAConnector(const char *nid = "", const TGeoManager *mgr=gGeoManager);


    // Basis Geometry
    TGeoVolume* CreateSTube45(TString, Double_t, Double_t, Double_t, Double_t, Double_t, 
		Bool_t, Int_t color, const TGeoMedium*);
	TGeoVolume* CreateTube(TString, Double_t, Double_t, Double_t,
		Bool_t, Int_t color, Double_t, const TGeoMedium*);
	TGeoVolume* CreateLTube(TString, Double_t, Double_t, Double_t, Double_t, Double_t, 
		Bool_t, Bool_t, Int_t color, const TGeoMedium*);
	
    // Outer Barrel Detector Barrel
	//vector<TString> CreateDBOuterHoles(TString, Double_t, Double_t, Double_t, Double_t, Double_t);
	//vector<TString> CreateDBOuterHandles(TString, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

    // --------------------------------	
    Int_t     fLayerNumber; // Current layer number
    Double_t  fPhi0;        // lab phi of 1st stave, in degrees!!!
    Double_t  fLayRadius;   // Inner radius of this layer
    Double_t  fZLength;     // Z length of this layer
    Double_t  fSensorThick; // Sensor thickness
    Double_t  fChipThick;   // Chip thickness
    Double_t  fStaveWidth;  // Stave width (for turbo layers only)
    Double_t  fStaveTilt;   // Stave tilt angle (for turbo layers only) in degrees
    Int_t     fNStaves;     // Number of staves in this layer
    Int_t     fNModules;    // Number of modules per container if defined (HalfStave, Stave, whatever is container)
    Int_t     fNChips;      // N. chips per container (module, HalfStave, Stave, whatever is container)
    Int_t     fHierarchy[kNHLevels]; // array to query number of staves, hstaves, modules, chips per its parent volume
    //    
    UInt_t    fChipTypeID;  // detector type id
    Bool_t    fIsTurbo;     // True if this layer is a "turbo" layer
    Int_t     fBuildLevel;  // Used for material studies

    AliITSUv2::AliITSUModel_t fStaveModel; // The stave model
    
    // ---------------------------------
    // #pnamwong
    Int_t fNHandleCreated;  // Number of handles which are created in hole geometry.
    Int_t fN_DCDC_Created;  // Number of DCDCset which are created in hole geometry.	
    Int_t fN_DCCNT_Created; // Number of DC Connector which are created in hole geometry.
    Int_t fN_SamtecCable_Created;
    // ---------------------------------    

    // Parameters for the Upgrade geometry

    // General Parameters
    static const Int_t    fgkNumberOfInnerLayers;// Number of IB Layers

    static const Double_t fgkDefaultSensorThick; // Default sensor thickness
    static const Double_t fgkDefaultChipThick;   // Default chip thickness

    // Inner Barrel Parameters
    static const Int_t    fgkIBChipsPerRow;      // IB chips per row in module
    static const Int_t    fgkIBNChipRows;        // IB chip rows in module

    static const Double_t fgkIBFlexCableAlThick; // Thickness of FPC Aluminum
    static const Double_t fgkIBFlexCableKapThick;// Thickness of FPC Kapton
    static const Double_t fgkIBGlueThick;        // IB glue thickness
    static const Double_t fgkIBCarbonFleeceThick;// IB carbon fleece thickness
    static const Double_t fgkIBCarbonPaperThick; // IB Carbon Paper Thickness
    static const Double_t fgkIBK13D2UThick;      // IB k13d2u prepreg thickness
    static const Double_t fgkIBCoolPipeInnerD;   // IB cooling inner diameter
    static const Double_t fgkIBCoolPipeThick;    // IB cooling pipe thickness
    static const Double_t fgkIBCoolPipeXDist;    // IB cooling pipe separation
    static const Double_t fgkIBTopVertexWidth1;  // IB TopVertex width
    static const Double_t fgkIBTopVertexWidth2;  // IB TopVertex width
    static const Double_t fgkIBTopVertexHeight;  // IB TopVertex height
    static const Double_t fgkIBTopVertexAngle ;  // IB TopVertex aperture angle
    static const Double_t fgkIBSideVertexWidth;  // IB SideVertex width
    static const Double_t fgkIBSideVertexHeight; // IB SideVertex height
    static const Double_t fgkIBTopFilamentLength;// IB TopFilament length
    static const Double_t fgkIBTopFilamentSide;  // IB TopFilament side
    static const Double_t fgkIBTopFilamentAlpha; // IB TopFilament angle
    static const Double_t fgkIBTopFilamentGamma; // IB TopFilament angle

    static const Double_t fgkIBConnectorXWidth;  // IB Connectors Width
    static const Double_t fgkIBConnectorYTot;    // IB Connectors total height
    static const Double_t fgkIBConnectBlockZLen; // IB Connector Block Z length
    static const Double_t fgkIBConnBodyYHeight;  // IB Connector Body Y height
    static const Double_t fgkIBConnTailYMid;     // IB Connector Tail Y mid pt
    static const Double_t fgkIBConnTailYShift;   // IB Connector Tail Y shift
    static const Double_t fgkIBConnTailZLen;     // IB Connector Tail Z length
    static const Double_t fgkIBConnTailOpenPhi;  // IB Connector Tail Angle
    static const Double_t fgkIBConnRoundHoleD;   // IB Connector Hole diameter
    static const Double_t fgkIBConnRoundHoleZ;   // IB Connector Hole Z pos
    static const Double_t fgkIBConnSquareHoleX;  // IB Connector Hole X len
    static const Double_t fgkIBConnSquareHoleZ;  // IB Connector Hole Z len
    static const Double_t fgkIBConnSquareHoleZPos;//IB Connector Hole Z pos
    static const Double_t fgkIBConnInsertHoleD;  // IB Connector Insert diam
    static const Double_t fgkIBConnInsertHoleZPos;//IB Connector Insert Z pos
    static const Double_t fgkIBConnTubeHole1D;   // IB Connector Tube1 diam
    static const Double_t fgkIBConnTubeHole1ZLen;// IB Connector Tube1 Z len
    static const Double_t fgkIBConnTubeHole2D;   // IB Connector Tube2 diam
    static const Double_t fgkIBConnTubeHole3XPos;// IB Connector Tube3 X pos
    static const Double_t fgkIBConnTubeHole3ZPos;// IB Connector Tube3 Z pos
    static const Double_t fgkIBConnTubesXDist;   // IB Connector Tubes X dist
    static const Double_t fgkIBConnTubesYPos;    // IB Connector Tubes Y pos
    static const Double_t fgkIBConnInsertInnerX; // IB Connector Insert X in
    static const Double_t fgkIBConnInsertZThick; // IB Connector Insert Z thick
    static const Double_t fgkIBConnInsertD;      // IB Connector Insert diam
    static const Double_t fgkIBConnInsertHeight; // IB Connector Insert height
    static const Double_t fgkIBConnectAFitExtD;  // IB ConnectorA Fitting ext D
    static const Double_t fgkIBConnectAFitIntD;  // IB ConnectorA Fitting int D
    static const Double_t fgkIBConnectAFitZLen;  // IB ConnectorA Fitting Z len
    static const Double_t fgkIBConnectAFitZOut;  // IB ConnectorA Fitting Z Out
    static const Double_t fgkIBConnPlugInnerD;   // IB Connector Plug int diam
    static const Double_t fgkIBConnPlugTotLen;   // IB Connector Plug tot le
    static const Double_t fgkIBConnPlugThick;    // IB Connector Plug thickness

    static const Double_t fgkIBStaveHeight;      // IB Stave Total Y Height

    // Outer Barrel Parameters
    static const Int_t    fgkOBChipsPerRow;      // OB chips per row in module
    static const Int_t    fgkOBNChipRows;        // OB chip rows in module

    static const Double_t fgkOBHalfStaveWidth;   // OB Half Stave Width
    static const Double_t fgkOBModuleWidth;      // OB Module Width
    static const Double_t fgkOBModuleGap;        // Gap between OB modules
    static const Double_t fgkOBChipXGap;         // Gap between OB chips on X
    static const Double_t fgkOBChipZGap;         // Gap between OB chips on Z
    static const Double_t fgkOBFlexCableAlThick; // Thickness of FPC Aluminum
    static const Double_t fgkOBFlexCableCuThick; // Thickness of FPC Copper
    static const Double_t fgkOBFlexCableKapThick1;// Thickness of FPC Kapton
    static const Double_t fgkOBFlexCableKapThick;// Thickness of FPC Kapton
    static const Double_t fgkOBBusCableAlThick;  // Thickness of Bus Aluminum
    static const Double_t fgkOBBusCableKapThick; // Thickness of Bus Kapton
    static const Double_t fgkOBCarbonPlateThick; // OB Carbon Plate Thickness
    static const Double_t fgkOBColdPlateThick;   // OB Cold Plate Thickness
    static const Double_t fgkOBGlueThickM1;      // OB Glue total Thickness
    static const Double_t fgkOBGlueThick;        // OB Glue Thickness in Model2
    static const Double_t fgkOBModuleZLength;    // OB Chip Length along Z
    static const Double_t fgkOBHalfStaveYPos;    // OB half staves Y position
    static const Double_t fgkOBHalfStaveYTrans;  // OB half staves Y transl.
    static const Double_t fgkOBHalfStaveXOverlap;// OB half staves X overlap
    static const Double_t fgkOBGraphiteFoilThick;// OB graphite foil thickness
    static const Double_t fgkOBCarbonFleeceThick;// OB carbon fleece thickness
    static const Double_t fgkOBCoolTubeInnerDM1; // OB cooling inner diameter
    static const Double_t fgkOBCoolTubeInnerD;   // OB cooling inner diameter
    static const Double_t fgkOBCoolTubeThick;    // OB cooling tube thickness
    static const Double_t fgkOBCoolTubeXDist;    // OB cooling tube separation

    static const Double_t fgkOBCPConnectorXWidth;// OB Cold Plate Connect Width
    static const Double_t fgkOBCPConnBlockZLen;  // OB CP Connect Block Z len
    static const Double_t fgkOBCPConnBlockYHei;  // OB CP Connect Block Z len
    static const Double_t fgkOBCPConnHollowZLen; // OB CP Connect Block Z len
    static const Double_t fgkOBCPConnHollowYHei; // OB CP Connect Block Z len
    static const Double_t fgkOBCPConnSquareHoleX;// OB Conn Square Hole X len
    static const Double_t fgkOBCPConnSquareHoleZ;// OB Conn Square Hole Z len
    static const Double_t fgkOBCPConnSqrHoleZPos;// OB Conn Square Hole Z pos
    static const Double_t fgkOBCPConnSqrInsertRZ;// OB Conn Square Insert RZpos
    static const Double_t fgkOBCPConnRoundHoleD; // OB Conn Round Hole diam
    static const Double_t fgkOBCPConnRndHoleZPos;// OB Conn Round Hole Z pos
    static const Double_t fgkOBCPConnTubesXDist; // OB Connector Tubes X dist
    static const Double_t fgkOBCPConnTubesYPos;  // OB Connector Tubes Y pos
    static const Double_t fgkOBCPConnTubeHole1D; // OB Connector Tube1 diam
    static const Double_t fgkOBCPConnTubeHole1Z; // OB Connector Tube1 Z len
    static const Double_t fgkOBCPConnTubeHole2D; // OB Connector Tube2 diam
    static const Double_t fgkOBCPConnFitHoleD;   // OB Connector Fit Hole diam
    static const Double_t fgkOBCPConnTubeHole3XP;// OB Connector Tube3 X pos
    static const Double_t fgkOBCPConnTubeHole3ZP;// OB Connector Tube3 Z pos
    static const Double_t fgkOBCPConnInstInnerX; // OB Connector Insert X in
    static const Double_t fgkOBCPConnInstInnerR; // OB Connector Insert R in
    static const Double_t fgkOBCPConnInstZThick; // OB Connector Insert height
    static const Double_t fgkOBCPConnInsertYHei; // OB Connector Insert height
    static const Double_t fgkOBCPConnInsertD;    // OB Connector Insert diam
    static const Double_t fgkOBCPConnAFitExtD;   // OB ConnectorA Fitting ext D
    static const Double_t fgkOBCPConnAFitThick;  // OB ConnectorA Fitting thick
    static const Double_t fgkOBCPConnAFitZLen;   // OB ConnectorA Fitting Z len
    static const Double_t fgkOBCPConnAFitZOut;   // OB ConnectorA Fitting Z Out
    static const Double_t fgkOBCPConnPlugInnerD; // OB Connector Plug int diam
    static const Double_t fgkOBCPConnPlugTotLen; // OB Connector Plug tot le
    static const Double_t fgkOBCPConnPlugThick;  // OB Connector Plug thickness

    static const Double_t fgkOBSpaceFrameZLen[2];// OB Space Frame Length
    static const Int_t    fgkOBSpaceFrameNUnits[2];//OB Number of SF Units
    static const Double_t fgkOBSpaceFrameUnitLen;// OB Space Frame Unit length
    static const Double_t fgkOBSpaceFrameWidth;  // OB Space Frame Width
    static const Double_t fgkOBSpaceFrameHeight; // OB Space Frame Height
    static const Double_t fgkOBSpaceFrameTopVL;  // Parameters defining...
    static const Double_t fgkOBSpaceFrameTopVH;  // ...the Top V shape
    static const Double_t fgkOBSpaceFrameSideVL; // Parameters defining...
    static const Double_t fgkOBSpaceFrameSideVH; // ...the Side V shape
    static const Double_t fgkOBSpaceFrameVAlpha; // Angles of aperture...
    static const Double_t fgkOBSpaceFrameVBeta;  // ...of the V shapes
    static const Double_t fgkOBSFrameBaseRibDiam;// OB SFrame Base Rib Diam
    static const Double_t fgkOBSFrameBaseRibPhi; // OB SF base beam angle
    static const Double_t fgkOBSFrameSideRibDiam;// OB SFrame Side Rib Diam
    static const Double_t fgkOBSFrameSideRibPhi; // OB SF side beam angle
    static const Double_t fgkOBSFrameULegLen;    // OB SF U-Leg length
    static const Double_t fgkOBSFrameULegWidth;  // OB SF U-Leg width
    static const Double_t fgkOBSFrameULegHeight1;// OB SF U-Leg height
    static const Double_t fgkOBSFrameULegHeight2;// OB SF U-Leg height
    static const Double_t fgkOBSFrameULegThick;  // OB SF U-Leg thickness
    static const Double_t fgkOBSFrameULegXPos;   // OB SF U-Leg X position

	// Barrel Geomatrix Parameters 
    // #pnamwong
//	static const Double_t AliITSUv2Layer::wrpZSpan[3];
	
	// Inner Barrel
	static const Double_t fgkIBStaveLength;
	static Double_t fgkIBConnOffset;
	static const Double_t fgkIBConnOffsetExt;

	// Inner Barrel - End Wheel 
	//static const TGeoMedium *fkIBEWheelMedium;
	static const Double_t fkALC_0334_RotateAngle[3];
	static const Double_t fkALC_0334_RotateOffset[3];
	
	static const Double_t fkALC_0334_PlateRadius[3];
	static const Double_t fkALC_0334_PlateThick;

	static const Double_t fkALC_0334_ContactDX[3];
	static const Double_t fkALC_0334_ContactDY[3];

	static const Double_t fkALC_0334_HoleDia;
	static const Double_t fkALC_0334_HoleFromEdge;		// center-of-hole to inside-edge
	static const Double_t fkALC_0334_HoleLength;		// approximate for hole object

	static const Double_t fkALC_0334_HndLength;
	static const Double_t fkALC_0334_HndWidth;
	static const Double_t fkALC_0334_HndThick;
	
	static const Double_t fkThermalConnOffsetZ[3];
	static const Double_t fkThermalConnOffsetY[3];
	
	static const Double_t fkDCDCOffsetZ[3];
	static const Double_t fkDCDCOffsetY[3];
	static const Double_t fkDCDCRotateBeta[3];

//---------------------------------------------------
	static const Double_t fgkInnerBarrelOffsetX;
	static const Double_t fgkInnerBarrelOffsetY;
	static const Double_t fgkInnerBarrelOffsetZ;

	static const Double_t fgkOuterBarrelOffsetX;
	static const Double_t fgkOuterBarrelOffsetY;
	static const Double_t fgkOuterBarrelOffsetZ;
	
	static const Double_t fgkBarrelRotateAlpha;
	static const Double_t fgkBarrelRotateBetha;
	static const Double_t fgkBarrelRotateGamma;


  ClassDef(AliITSUv2Layer,0) // ITS Upgrade v2 geometry
};

#endif
