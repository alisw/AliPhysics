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

//#include <assert.h> 

// --- ROOT system ---
class TString ;
class TObjArray ;
class TVector3 ;
class TParticle ; 

// --- AliRoot header files ---

#include "AliGeometry.h"

class AliEMCALGeometry : public AliGeometry {
public:
  AliEMCALGeometry() {
    // default ctor,  must be kept public for root persistency purposes,
    // but should never be called by the outside world
  };
  AliEMCALGeometry(const AliEMCALGeometry& geom):AliGeometry(geom) {
    // cpy ctor requested by Coding Convention but not yet needed
    Fatal("Cpy ctor", "Not implemented");
  };
  virtual ~AliEMCALGeometry(void) ; 
  static AliEMCALGeometry * GetInstance(const Text_t* name,
					const Text_t* title="") ; 
  static AliEMCALGeometry * GetInstance() ;
  AliEMCALGeometry & operator = (const AliEMCALGeometry  & /*rvalue*/) const {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *(GetInstance()) ; 
  };

  Bool_t AreInSameTower(Int_t id1, Int_t id2) const ;  
  virtual void GetGlobal(const AliRecPoint *, TVector3 &, TMatrix &) const {}
  virtual void GetGlobal(const AliRecPoint *, TVector3 &) const {}
  virtual Bool_t Impact(const TParticle *) const {return kTRUE;}

  Bool_t IsInEMCAL(Double_t x, Double_t y, Double_t z) const;
  // General
  Bool_t  IsInitialized(void) const { return fgInit ; }
  	// Return EMCA geometrical parameters
  // geometry
  Float_t GetAlFrontThickness() const { return fAlFrontThick;}
  Float_t GetArm1PhiMin() const { return fArm1PhiMin ; }
  Float_t GetArm1PhiMax() const { return fArm1PhiMax ; }
  Float_t GetArm1EtaMin() const { return fArm1EtaMin;}
  Float_t GetArm1EtaMax() const { return fArm1EtaMax;}
  Float_t GetIPDistance() const { return fIPDistance;}   
  Float_t GetIP2ECASection() const { return ( GetIPDistance() + GetAlFrontThickness() + GetGap2Active() ) ; }   
  Float_t GetEnvelop(Int_t index) const { return fEnvelop[index] ; }  
  Float_t GetShellThickness() const { return fShellThickness ; }
  Float_t GetZLength() const { return fZLength ; } 
  Float_t GetGap2Active() const {return  fGap2Active ; }
  Float_t GetDeltaEta() const {return (fArm1EtaMax-fArm1EtaMin)/
				       ((Float_t)fNZ);}
  Float_t GetDeltaPhi() const {return (fArm1PhiMax-fArm1PhiMin)/
				       ((Float_t)fNPhi);}
  Int_t   GetNECLayers() const {return fNECLayers ;}
  Int_t   GetNZ() const {return fNZ ;}
  Int_t   GetNEta() const {return fNZ ;}
  Int_t   GetNPhi() const {return fNPhi ;}
  Int_t   GetNTowers() const {return fNPhi * fNZ ;}
  Float_t GetECPbRadThick()const {return fECPbRadThickness;}
  Float_t GetECScintThick() const {return fECScintThick;}
  Float_t GetSampling() const {return fSampling ; } 
  Bool_t IsInECA(Int_t index) const { if ( (index > 0 && (index <= GetNZ() * GetNPhi()))) return kTRUE; else return kFALSE ;}
 
  Float_t AngleFromEta(Float_t eta){ // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta){ // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }
  Int_t TowerIndex(Int_t iz,Int_t iphi) const; // returns tower index
  	// returns tower indexs iz, iphi.
  void TowerIndexes(Int_t index,Int_t &iz,Int_t &iphi) const;
  	// for a given tower index it returns eta and phi of center of that tower.
  void EtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi) const;
  	// returns x, y, and z (cm) on the inner surface of a given EMCAL Cell specified by relid.
  void XYZFromIndex(const Int_t *relid,Float_t &x,Float_t &y, Float_t &z) const;
  void XYZFromIndex(Int_t absid, TVector3 &v) const;
  	// for a given eta and phi in the EMCAL it returns the tower index.
  Int_t TowerIndexFromEtaPhi(Float_t eta,Float_t phi) const;
  	// for a given eta and phi in the EMCAL it returns the pretower index.
  void PosInAlice(const Int_t *relid, Float_t &theta, Float_t &phi) const ;
  void PosInAlice(Int_t absid, Float_t &theta, Float_t &phi) const ;
  Bool_t AbsToRelNumbering(Int_t AbsId, Int_t *relid) const;
  void SetNZ(Int_t nz) { fNZ= nz ; printf("SetNZ: Number of modules in Z set to %d", fNZ) ; }
  void SetNPhi(Int_t nphi) { fNPhi= nphi ; printf("SetNPhi: Number of modules in Phi set to %d", fNPhi) ; }
  void SetSampling(Float_t samp) { fSampling = samp; printf("SetSampling: Sampling factor set to %f", fSampling) ; }

protected:
  AliEMCALGeometry(const Text_t* name, const Text_t* title="") :
    AliGeometry(name, title) {// ctor only for internal usage (singleton)
    Init();
  };
  void Init(void);     			// initializes the parameters of EMCAL
  
private:
  static AliEMCALGeometry * fgGeom;	// pointer to the unique instance of the singleton
  static Bool_t fgInit;			// Tells if geometry has been succesfully set up.
  Float_t fAlFrontThick;		// Thickness of the front Al face of the support box
  
  Float_t fECPbRadThickness;		// cm, Thickness of the Pb radiators
  Float_t fECScintThick;		// cm, Thickness of the scintillators
  Int_t   fNECLayers;			// number of scintillator layers
  
  Float_t fArm1PhiMin; 			// Minimum angular position of EMCAL in Phi (degrees)
  Float_t fArm1PhiMax;			// Maximum angular position of EMCAL in Phi (degrees)
  Float_t fArm1EtaMin;			// Minimum pseudorapidity position of EMCAL in Eta
  Float_t fArm1EtaMax; 			// Maximum pseudorapidity position of EMCAL in Eta
  
  // Geometry Parameters
  Float_t fEnvelop[3];			// the GEANT TUB for the detector 
  Float_t fIPDistance;			// Radial Distance of the inner surface of the EMCAL
  Float_t fShellThickness;		// Total thickness in (x,y) direction
  Float_t fZLength;			// Total length in z direction
  Float_t fGap2Active;			// Gap between the envelop and the active material
  Int_t   fNZ;				// Number of Towers in the Z direction
  Int_t   fNPhi;			// Number of Towers in the Phi Direction
  Float_t fSampling;			// Sampling factor
  
  ClassDef(AliEMCALGeometry,8) // EMCAL geometry class 
    
    };

#endif // AliEMCALGEOMETRY_H
