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

#include <map>

#include "AliMUONv1.h"
#include "AliMUONSt1SpecialMotif.h"

#include <TVector2.h>
#include <TVector3.h>

typedef Float_t GReal_t; // for AliGeant3
//typedef Double_t GReal_t; // for VirtualMC

class TTree;
class MSector;
using std::map;


class AliMUONv2 : public  AliMUONv1 {
  public:
    AliMUONv2();
    AliMUONv2(const char* name, const char* title);
    AliMUONv2(const AliMUONv2& rMUON);
    virtual ~AliMUONv2();

    virtual Int_t IsVersion() const;
    virtual void  CreateGeometry();
    virtual void  CreateMaterials();

  protected:
    // Copy Operator
    AliMUONv2& operator = (const AliMUONv2& rhs);
    
  private:
    
    typedef map< Int_t , AliMUONSt1SpecialMotif > TSpecialMap;
    static const GReal_t fgkHzPadPlane; // Pad plane
    static const GReal_t fgkHzFoam;  // Foam of mechanicalplane
    static const GReal_t fgkHzFR4;  // FR4 of mechanical plane
    static const GReal_t fgkHzSnPb; //Pad/Kapton connection (66 pt)
    static const GReal_t fgkHzKapton; //Kapton
    static const GReal_t fgkHzBergPlastic; //Berg connector 
    static const GReal_t fgkHzBergCopper; //Berg connector (80 pt)
    static const GReal_t fgkHzDaughter; //Daughter board
    static const GReal_t fgkHzGas ; // ArCO2 Gas

    GReal_t totalHzPlane() const ; //Total mechanical plane half Size
    GReal_t totalHzDaughter() const ; //Total daughter plane half Size
    GReal_t totalHz() const ; //Total plane half Size

    static const GReal_t fgkHxHole;
    static const GReal_t fgkHyHole;
    static const GReal_t fgkHxBergPlastic;
    static const GReal_t fgkHyBergPlastic;
    static const GReal_t fgkHxBergCopper;
    static const GReal_t fgkHyBergCopper;
    static const GReal_t fgkHxKapton;
    static const GReal_t fgkHyKapton;
    static const GReal_t fgkHxDaughter;
    static const GReal_t fgkHyDaughter;
    static const GReal_t fgkOffsetX;
    static const GReal_t fgkOffsetY;
    static const GReal_t fgkDeltaFilleEtamX;
    static const GReal_t fgkDeltaFilleEtamY;
    static const GReal_t fgkHxQuadrant;
    static const GReal_t fgkHyQuadrant;
    static const GReal_t fgkMotherIR;
    static const GReal_t fgkMotherOR; 
    static const GReal_t fgkMotherThick;  
    static const GReal_t fgkMotherPhiL; 
    static const GReal_t fgkMotherPhiU;

    static const char* fgkHoleName;
    static const char* fgkQuadrantName;
    static const char* fgkDaughterName;
    static const char fgkFoamLayerSuffix;


    void CreateHole();
    void CreateDaughterBoard();
    void CreateFrame(Int_t chamber);
    void CreateQuadrant(Int_t chamber);
    void CreatePlaneBox(const char* name,const TVector2& dimensions);
    void CreatePlaneSegment(const char* name,const  TVector2& dimensions
                      ,Int_t nofHoles);
    void CreateDaughterSegment(const char* name,const  TVector2& dimensions
                      ,Int_t nofHoles);
    void PlaceSector(MSector* sector,TSpecialMap specialMap
                    ,const TVector3& where,Int_t chamber);
    TString QuadrantName(Int_t chamber);
    Int_t      fIdSens;  // Sensitive volume identifier
    

  ClassDef(AliMUONv2,1)  // MUON Detector base class
};

// inline functions

inline Int_t AliMUONv2::IsVersion () const 
 { return 2; }
inline GReal_t AliMUONv2::totalHzPlane() const 
 { return fgkHzPadPlane + fgkHzFoam + fgkHzFR4;}
inline GReal_t AliMUONv2::totalHzDaughter() const 
 { return fgkHzBergPlastic + fgkHzDaughter;}
inline GReal_t AliMUONv2::totalHz() const 
 { return totalHzPlane() + totalHzDaughter();}
inline TString AliMUONv2::QuadrantName(Int_t chamber)
{return Form("%s%d",fgkQuadrantName,chamber);}


#endif //ALI_MUON_V2_H
