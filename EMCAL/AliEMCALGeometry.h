#ifndef ALIEMCALGEOMETRY_H
#define ALIEMCALGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton
// EMCAL consists of a layers of scintillator, and lead.
//                  
//*-- Author: Sahal Yacoob (LBL / UCT)
//*--   and : Yves Schutz (Subatech)

#include <assert.h> 

// --- ROOT system ---
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h"
class TParticle ; 

// --- AliRoot header files ---

#include "AliGeometry.h"

class AliEMCALGeometry : public AliGeometry {
public:
  AliEMCALGeometry() {
    // default ctor,  must be kept public for root persistency purposes,
    // but should never be called by the outside world
  };
  AliEMCALGeometry(const AliEMCALGeometry & geom) {
    // cpy ctor requested by Coding Convention but not yet needed
    assert(0==1);
  };
  virtual ~AliEMCALGeometry(void) ; 
  static AliEMCALGeometry * GetInstance(const Text_t* name,
					const Text_t* title="") ; 
  static AliEMCALGeometry * GetInstance() ;
  AliEMCALGeometry & operator = (const AliEMCALGeometry  & rvalue) const {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *(GetInstance()) ; 
  };
  virtual void GetGlobal(const AliRecPoint *, TVector3 &, TMatrix &) const {}
  virtual void GetGlobal(const AliRecPoint *, TVector3 &) const {}
  virtual Bool_t Impact(const TParticle * particle) const {return kTRUE;}
  // General
  Bool_t  IsInitialized(void) const { return fgInit ; }
  // Return EMCA geometrical parameters
  // geometry
  const Float_t GetAlFrontThickness() const { return fAlFrontThick;}
  const Float_t GetArm1PhiMin() const { return fArm1PhiMin ; }
  const Float_t GetArm1PhiMax() const { return fArm1PhiMax ; }
  const Float_t GetArm1EtaMin() const { return fArm1EtaMin;}
  const Float_t GetArm1EtaMax() const { return fArm1EtaMax;}
  const Float_t GetIPDistance() const { return fIPDistance  ; } 
  const Float_t GetIP2PreShower()  const { return (GetIPDistance() + GetAlFrontThickness() +  GetGap2Active() ) ;}
  const Float_t GetIP2Tower() const { return ( GetIP2PreShower() + 2 * (  GetPreSintThick() +  GetPbRadThick() ) ) ; }  
  const Float_t GetEnvelop(Int_t index) const { return fEnvelop[index] ; }  
  const Float_t GetShellThickness() const { return fShellThickness ; }
  const Float_t GetZLength() const { return fZLength ; } 
  const Float_t GetGap2Active() const {return  fGap2Active ; }
  const Float_t GetDeltaEta() const {return (fArm1EtaMax-fArm1EtaMin)/
				       ((Float_t)fNZ);}
  const Float_t GetDeltaPhi() const {return (fArm1PhiMax-fArm1PhiMin)/
				       ((Float_t)fNPhi);}
  const Int_t   GetNECLayers() const {return fNECLayers ;}
  const Int_t   GetNHCLayers() const {return fNHCLayers ;}
  const Int_t   GetNPRLayers() const {return fNPRLayers;}
  const Int_t   GetNZ() const {return fNZ ;}
  const Int_t   GetNEta() const {return fNZ ;}
  const Int_t   GetNPhi() const {return fNPhi ;}
  const Int_t   GetNTowers() const {return fNPhi * fNZ ;}
  const Float_t GetPbRadThick()const {return fPbRadThickness;}
  const Float_t GetCuRadThick()const {return fCuRadThickness;}
  const Float_t GetFullSintThick() const { // returns Full tower sintilator
    // thickness in cm.
    return fFullShowerSintThick;
  }
  const Float_t GetPreSintThick() const { // returns PreShower tower sintilator
    // thickness in cm.
    return fPreShowerSintThick;
  }
  Float_t AngleFromEta(Float_t eta){ // returns angle in radians for a given
    // pseudorapidity.
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta){ // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }
  Int_t TowerIndex(Int_t iz,Int_t iphi,Int_t ipre) const; // returns tower index
  // returns tower indexs iz, iphi.
  void TowerIndexes(Int_t index,Int_t &iz,Int_t &iphi,Int_t &ipre) const;
  // for a given tower index it returns eta and phi of center of that tower.
  void EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi) const;
  // returns x, y, and z (cm) on the inner surface of a given EMCAL Cell specified by relid.
  void XYZFromIndex(const Int_t *relid,Float_t &x,Float_t &y, Float_t &z) const;
  // for a given eta and phi in the EMCAL it returns the tower index.
  Int_t TowerIndexFromEtaPhi(Float_t eta,Float_t phi) const;
  // for a given eta and phi in the EMCAL it returns the pretower index.
  Int_t PreTowerIndexFromEtaPhi(Float_t eta,Float_t phi) const;
  // Returns theta and phi (degree) for a given EMCAL cell indecated by relid
  void PosInAlice(const Int_t *relid,Float_t &theta,Float_t &phi) const ;
  // Returns an array indicating the Tower/preshower, iz, and iphi for a
  // specific EMCAL indes.
  Bool_t AbsToRelNumbering(Int_t AbsId, Int_t *relid) const;
  /*
  // Returns kTRUE if the two indexs are neighboring towers or preshowers.
  Boot_t AliEMCALGeometry::AreNeighbours(Int_t index1,Int_t index2) const;
  */
  
  void SetNZ(Int_t nz) { fNZ= nz ; Info("SetNZ", "Number of modules in Z set to %d", fNZ) ; }
  void SetNPhi(Int_t nphi) { fNPhi= nphi ; Info("SetNPhi", "Number of modules in Phi set to %d", fNPhi) ; }
  
protected:
  AliEMCALGeometry(const Text_t* name, const Text_t* title="") :
    AliGeometry(name, title) {// ctor only for internal usage (singleton)
    Init();
  };
  void Init(void) ;            // initializes the parameters of EMCAL
  
private:
  static AliEMCALGeometry * fgGeom ; // pointer to the unique instance
  // of the singleton 
  static Bool_t fgInit;// Tells if geometry has been succesfully set up.
  Float_t fAlFrontThick; // Thickness of the front Al face of the support box
  Float_t fPreShowerSintThick; // Thickness of the sintilator for the
  // preshower part of the calorimeter
  Float_t fFullShowerSintThick;// Thickness of the sintilaor for the full
  // shower part of the calorimeter
  Float_t fPbRadThickness; // Thickness of Pb radiators cm.
  Float_t fCuRadThickness; // Thickness of Cu radiators cm.
  Float_t fArm1PhiMin; // Minimum angular position of EMCAL in Phi (degrees)
  Float_t fArm1PhiMax; // Maximum angular position of EMCAL in Phi (degrees)
  Float_t fArm1EtaMin; // Minimum pseudorapidity position of EMCAL in Eta
  Float_t fArm1EtaMax; // Maximum pseudorapidity position of EMCAL in Eta
  
  // It is assumed that Arm1 and Arm2 have the same following parameters
  Float_t fEnvelop[3];     // the GEANT TUB for the detector 
  Float_t fIPDistance;     // Radial Distance of the inner surface of the EMCAL
  Float_t fShellThickness; // Total thickness in (x,y) direction
  Float_t fZLength;        // Total length in z direction
  Float_t fGap2Active;     // Gap between the envelop and the active material
  Int_t   fNECLayers;        // Number of layers of material in the R direction for the electromagnetic calorimeter 
  Int_t   fNPRLayers;      // Number of layers of material in the R direction for the preshower
  Int_t   fNHCLayers;      // Number of layers of material in the R direction for the hadron calorimeter
  Int_t   fNZ;             // Number of Towers in the Z direction
  Int_t   fNPhi;           // Number of Towers in the Phi Direction
  
  ClassDef(AliEMCALGeometry,5) // EMCAL geometry class 
    
    };

#endif // AliEMCALGEOMETRY_H
