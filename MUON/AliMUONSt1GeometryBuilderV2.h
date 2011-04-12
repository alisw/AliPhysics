/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
/// \ingroup sim
/// \class AliMUONSt1GeometryBuilderV2
/// \brief MUON Station1 detailed geometry construction class
///
/// \author David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay

#ifndef ALI_MUON_ST1_GEOMETRY_BUILDER_V2_H
#define ALI_MUON_ST1_GEOMETRY_BUILDER_V2_H


#include "AliMUONVGeometryBuilder.h"

#include <TExMap.h>

// typedef Float_t GReal_t; // for AliGeant3
typedef Double_t GReal_t;  // for VirtualMC

class AliMUON;
class AliMpSector;

class TTree;
class TVector2;
class TVector3;

class AliMUONSt1GeometryBuilderV2 : public AliMUONVGeometryBuilder 
{
  public:
    AliMUONSt1GeometryBuilderV2(AliMUON* muon);
    AliMUONSt1GeometryBuilderV2();
    virtual ~AliMUONSt1GeometryBuilderV2();

    virtual void CreateMaterials();
    virtual void CreateGeometry();
    virtual void SetVolumes();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
   
  protected:
 
  private:
    /// Not implemented
    AliMUONSt1GeometryBuilderV2(const AliMUONSt1GeometryBuilderV2& rMUON);
    /// Not implemented
    AliMUONSt1GeometryBuilderV2& operator = (const AliMUONSt1GeometryBuilderV2& rhs);    

    // Constants
    //
    static const GReal_t fgkHzPadPlane;    ///< Pad plane
    static const GReal_t fgkHzFoam;        ///< Foam of mechanicalplane
    static const GReal_t fgkHzFR4;         ///< FR4 of mechanical plane
    static const GReal_t fgkHzSnPb;        ///< Pad/Kapton connection (66 pt)
    static const GReal_t fgkHzKapton;      ///< Kapton
    static const GReal_t fgkHzBergPlastic; ///< Berg connector 
    static const GReal_t fgkHzBergCopper;  ///< Berg connector (80 pt)
    static const GReal_t fgkHzDaughter;    ///< Daughter board
    static const GReal_t fgkHzGas;         ///< ArCO2 Gas
        
    // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxHole;          ///< foam hole paremeter
    static const GReal_t fgkHyHole;          ///< foam hole paremeter
    static const GReal_t fgkHxBergPlastic;   ///< Berg connector parameter
    static const GReal_t fgkHyBergPlastic;   ///< Berg connector parameter
    static const GReal_t fgkHxBergCopper;    ///< Berg connector parameter
    static const GReal_t fgkHyBergCopper;    ///< Berg connector parameter
    static const GReal_t fgkHxKapton;        ///< Kapton parameter
    static const GReal_t fgkHyKapton;        ///< Kapton parameter
    static const GReal_t fgkHxDaughter;      ///< Electronics parameter
    static const GReal_t fgkHyDaughter;      ///< Electronics parameter
    static const GReal_t fgkOffsetX;         ///< Offset X
    static const GReal_t fgkOffsetY;         ///< Offset Y
    static const GReal_t fgkDeltaFilleEtamX; ///< Electronics parameter
    static const GReal_t fgkDeltaFilleEtamY; ///< Electronics parameter

    static const GReal_t fgkDeltaQuadLHC; ///< LHC Origin wrt Quadrant Origin
    static const GReal_t fgkFrameOffset;  ///< Frame offset

    // Pad planes offsets
    static const GReal_t fgkPadXOffsetBP; ///< Horizontal offset in bending plane  
    static const GReal_t fgkPadYOffsetBP; ///< Vertical offset in bending plane 

    // Quadrant Mother volume - TUBS1   
    static const GReal_t fgkMotherIR1;    ///< Middle Layer Rin
    static const GReal_t fgkMotherOR1;    ///< Middle Layer Rout
    static const GReal_t fgkMotherThick1; ///< Middle Layer Hz 
    static const GReal_t fgkMotherPhiL1;  ///< Middle Layer Sphi
    static const GReal_t fgkMotherPhiU1;  ///< Middle Layer Endphi

    // Quadrant Mother volume - TUBS2 (2 copies at different Z's)   
    static const GReal_t fgkMotherIR2;    ///< Near and Far Layer Rin
    static const GReal_t fgkMotherOR2;    ///< Near and Far Layer Rout
    static const GReal_t fgkMotherThick2; ///< Near and Far Layer Hz 
    static const GReal_t fgkMotherPhiL2;  ///< Near and Far Layer Sphi
    static const GReal_t fgkMotherPhiU2;  ///< Near and Far Layer Endphi  

    static const char* fgkHoleName;          ///< prefix for automatic volume naming
    static const char* fgkQuadrantEnvelopeName; ///< prefix for automatic volume naming
    static const char* fgkQuadrantMLayerName;///< prefix for automatic volume naming
    static const char* fgkQuadrantNLayerName;///< prefix for automatic volume naming
    static const char* fgkQuadrantFLayerName;///< prefix for automatic volume naming
    static const char* fgkQuadrantMFLayerName;    ///< prefix for automatic volume naming
    static const char* fgkDaughterName;      ///< prefix for automatic volume naming
    static const Int_t fgkFoamBoxNameOffset; ///< coefficient for automatic volume naming
    static const Int_t fgkFR4BoxNameOffset;  ///< coefficient for automatic volume naming
    static const Int_t fgkDaughterCopyNoOffset; ///< \brief copy number offset for daughter 
                                                /// boards positions in non-bending plane

    // Methods
    //
    void CreateHole();
    void CreateDaughterBoard();
    void CreateInnerLayers();
    void CreateSpacer0();
    void CreateSpacer();
    void CreateQuadrant(Int_t chamber);
    void CreateFoamBox(Int_t segNumber, const TVector2& dimensions);
    void CreatePlaneSegment(Int_t segNumber, const TVector2& dimensions,
                     Int_t nofHoles);
    void CreateQuadrantLayersAsVolumes(Int_t chamber);
    void CreateQuadrantLayersAsAssemblies(Int_t chamber);
    void CreateFrame(Int_t chamber);

    void PlaceInnerLayers(Int_t chamber);
    void PlaceSpacer0(Int_t chamber);
    void PlaceSector(const AliMpSector* sector, TExMap specialMap,
                     const TVector3& where, Bool_t reflectZ, Int_t chamber);
		     
    TString QuadrantEnvelopeName(Int_t chamber, Int_t quadrant) const;
    TString QuadrantMLayerName(Int_t chamber) const;
    TString QuadrantNLayerName(Int_t chamber) const;
    TString QuadrantFLayerName(Int_t chamber) const;
    TString QuadrantMFLayerName(Int_t chamber) const;
    TString PlaneSegmentName(Int_t segNumber) const;
    TString FoamBoxName(Int_t segNumber) const;
    TString FR4BoxName(Int_t segNumber) const;
    TString GasVolumeName(const TString& name, Int_t chamber) const;

    GReal_t TotalHzPlane() const ; 
    GReal_t TotalHzDaughter() const ;
    GReal_t TotalHz() const ;
       
    // Data members
    //
    //Float_t  fRadlCopper;  //! copper computed radiation length
    //Float_t  fRadlFoam;    //! foam   computed radiation length
    //Float_t  fRadlFR4;     //! FR4    computed radiation length
    AliMUON*  fMUON; ///< the MUON detector class 
    
  ClassDef(AliMUONSt1GeometryBuilderV2,1)  // MUON Detector base class
};

// inline functions

/// Return total mechanical plane half Size
inline GReal_t AliMUONSt1GeometryBuilderV2::TotalHzPlane() const 
//{ return fgkHzPadPlane + fgkHzFoam + fgkHzFR4; }
{ return fgkHzFoam + fgkHzFR4; }

/// Return total daughter plane half Size
inline GReal_t AliMUONSt1GeometryBuilderV2::TotalHzDaughter() const 
{ return fgkHzBergPlastic + fgkHzDaughter; }

/// Return total plane half Size
inline GReal_t AliMUONSt1GeometryBuilderV2::TotalHz() const 
{ return TotalHzPlane() + TotalHzDaughter(); }

/// Return middle quadrant layer name for chamber \a chamber
inline TString AliMUONSt1GeometryBuilderV2::QuadrantMLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantMLayerName,chamber); }

/// Return middle quadrant frame layer name for chamber \a chamber
inline TString AliMUONSt1GeometryBuilderV2::QuadrantMFLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantMFLayerName,chamber); }

/// Return nearer quadrant layer name for chamber \a chamber
inline TString AliMUONSt1GeometryBuilderV2::QuadrantNLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantNLayerName,chamber); }

/// Return farther quadrant layer name for chamber \a chamber
inline TString AliMUONSt1GeometryBuilderV2::QuadrantFLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantFLayerName,chamber); }

/// Return plane segment name for segment \a segNumber
inline TString AliMUONSt1GeometryBuilderV2::PlaneSegmentName(Int_t segNumber) const
{ return Form("S%.3d", segNumber); }

/// Return foam box name for segment \a segNumber
inline TString AliMUONSt1GeometryBuilderV2::FoamBoxName(Int_t segNumber) const
{ return Form("S%.3d", segNumber + fgkFoamBoxNameOffset); }

/// Return FR4 box name for segment \a segNumber
inline TString AliMUONSt1GeometryBuilderV2::FR4BoxName(Int_t segNumber) const
{ return Form("S%.3d", segNumber + fgkFR4BoxNameOffset); }

#endif //ALI_MUON_ST1_GEOMETRY_BUILDER_V2_H
