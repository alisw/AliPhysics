#ifndef ALI_MUON_V2_H
#define ALI_MUON_V2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONv2
// ---------------
// Inherits from AliMUONv1 but with a more detailed
// geometrical description of station 1 

#include <TVector2.h>
#include <TVector3.h>

#include "AliMUONv1.h"
#include "AliMUONSt1Types.h"
#include "AliMUONSt1SpecialMotif.h"

//typedef Float_t GReal_t; // for AliGeant3
typedef Double_t GReal_t; // for VirtualMC

class TTree;
class TArrayI;
class AliMpSector;

class AliMUONv2 : public  AliMUONv1 
{
  public:
    AliMUONv2();
    AliMUONv2(const char* name, const char* title);
    AliMUONv2(const AliMUONv2& rMUON);
    virtual ~AliMUONv2();

    virtual Int_t IsVersion() const;
    virtual void  CreateMaterials();
    virtual void  CreateGeometry();
    virtual void  Init();
   
  protected:
    // Copy Operator
    AliMUONv2& operator = (const AliMUONv2& rhs);    
    virtual Int_t  GetChamberId(Int_t volId) const;

  private:
    // Methods
    //
    void CreateHole();
    void CreateDaughterBoard();
    void CreateInnerLayers();
    void CreateQuadrant(Int_t chamber);
    void CreateFoamBox(const char* name,const TVector2& dimensions);
    void CreatePlaneSegment(const char* name,const TVector2& dimensions,
                     Int_t nofHoles);
    void CreateFrame(Int_t chamber);

    void PlaceInnerLayers(Int_t chamber);
    void PlaceSector(AliMpSector* sector, TSpecialMap specialMap,
                     const TVector3& where, Bool_t reflectZ, Int_t chamber);
		     
    TString QuadrantMLayerName(Int_t chamber) const;
    TString QuadrantNLayerName(Int_t chamber) const;
    TString QuadrantFLayerName(Int_t chamber) const;
    TString GasVolumeName(const TString& name, Int_t chamber) const;

    void   AddChamberGid(Int_t id,Int_t volName,Int_t idx);
    Bool_t IsInChamber(Int_t ich, Int_t volGid) const;   

    GReal_t TotalHzPlane() const ;         // Total mechanical plane half Size
    GReal_t TotalHzDaughter() const ;      // Total daughter plane half Size
    GReal_t TotalHz() const ;              // Total plane half Size
       
    // Constants
    //
    static const GReal_t fgkHzPadPlane;    // Pad plane
    static const GReal_t fgkHzFoam;        // Foam of mechanicalplane
    static const GReal_t fgkHzFR4;         // FR4 of mechanical plane
    static const GReal_t fgkHzSnPb;        // Pad/Kapton connection (66 pt)
    static const GReal_t fgkHzKapton;      // Kapton
    static const GReal_t fgkHzBergPlastic; // Berg connector 
    static const GReal_t fgkHzBergCopper;  // Berg connector (80 pt)
    static const GReal_t fgkHzDaughter;    // Daughter board
    static const GReal_t fgkHzGas;         // ArCO2 Gas
        
    // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxHole;   // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHyHole;   // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxBergPlastic; // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHyBergPlastic; // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxBergCopper;  // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHyBergCopper;  // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxKapton;      // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHyKapton;      // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHxDaughter;    // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkHyDaughter;    // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkOffsetX;       // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkOffsetY;       // Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkDeltaFilleEtamX;// Sensitive copper pads, foam layer, PCB and electronics model parameters
    static const GReal_t fgkDeltaFilleEtamY;// Sensitive copper pads, foam layer, PCB and electronics model parameters

    static const GReal_t fgkDeltaQuadLHC; //LHC Origin wrt Quadrant Origin
    static const GReal_t fgkFrameOffset;  

    // Quadrant Mother volume - TUBS1   
    static const GReal_t fgkMotherIR1;
    static const GReal_t fgkMotherOR1; 
    static const GReal_t fgkMotherThick1;  
    static const GReal_t fgkMotherPhiL1; 
    static const GReal_t fgkMotherPhiU1;

    // Quadrant Mother volume - TUBS2 (2 copies at different Z's)   
    static const GReal_t fgkMotherIR2;
    static const GReal_t fgkMotherOR2; 
    static const GReal_t fgkMotherThick2;  
    static const GReal_t fgkMotherPhiL2; 
    static const GReal_t fgkMotherPhiU2;    

    static const char* fgkHoleName;          // prefix for automatic volume naming
    static const char* fgkQuadrantMLayerName;// prefix for automatic volume naming
    static const char* fgkQuadrantNLayerName;// prefix for automatic volume naming
    static const char* fgkQuadrantFLayerName;// prefix for automatic volume naming
    static const char* fgkDaughterName;      // prefix for automatic volume naming
    static const char  fgkFoamLayerSuffix;   // suffix for automatic volume naming

    // Data members
    //
    Float_t  fRadlCopper;  //! copper computed radiation length
    Float_t  fRadlFoam;    //! foam   computed radiation length
    Float_t  fRadlFR4;     //! FR4    computed radiation length
    TArrayI* fChamberV2[2];// Sensitive volumes IDs    
    
  ClassDef(AliMUONv2,1)  // MUON Detector base class
};

// inline functions

inline Int_t AliMUONv2::IsVersion () const 
{ return 2; }

inline GReal_t AliMUONv2::TotalHzPlane() const 
//{ return fgkHzPadPlane + fgkHzFoam + fgkHzFR4; }
{ return fgkHzFoam + fgkHzFR4; }

inline GReal_t AliMUONv2::TotalHzDaughter() const 
{ return fgkHzBergPlastic + fgkHzDaughter; }

inline GReal_t AliMUONv2::TotalHz() const 
{ return TotalHzPlane() + TotalHzDaughter(); }

inline TString AliMUONv2::QuadrantMLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantMLayerName,chamber); }

inline TString AliMUONv2::QuadrantNLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantNLayerName,chamber); }

inline TString AliMUONv2::QuadrantFLayerName(Int_t chamber) const
{ return Form("%s%d",fgkQuadrantFLayerName,chamber); }

inline void AliMUONv2::AddChamberGid(Int_t id,Int_t volName,Int_t idx)
{ fChamberV2[id]->AddAt(volName,idx); }


#endif //ALI_MUON_V2_H
