#ifndef ALIEMCALGEOMETRY_H
#define ALIEMCALGEOMETRY_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton
// EMCAL consists of a layers of scintillator, and lead.
//                  
//*-- Author: Sahal Yacoob (LBL / UCT)
//*--   and : Yves Schutz (Subatech)
//*--   and : Aleksei Pavlinov (WSU) - shashlyk staff
//*--   and : Gustavo Conesa: Add TRU mapping. TRU parameters still not fixed.

// --- ROOT system ---
class TString ;
class TObjArray;
class TVector3 ;
class TParticle ; 
class TClonesArray ;

// --- AliRoot header files ---
#include "AliGeometry.h"

class AliEMCALGeometry : public AliGeometry {
public:
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

  TClonesArray *  FillTRU(const TClonesArray * digits)  ;
  
  virtual void GetGlobal(const AliRecPoint *, TVector3 &, TMatrixF &) const {}
  virtual void GetGlobal(const AliRecPoint *, TVector3 &) const {}
  virtual Bool_t Impact(const TParticle *) const {return kTRUE;}

  Bool_t IsInEMCAL(Double_t x, Double_t y, Double_t z) const;
  // General
  Bool_t  IsInitialized(void) const { return fgInit ; }
  // Return EMCAL geometrical parameters
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
  Float_t GetGap2Active() const {return  fGap2Active ;}
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

  Int_t   GetNumberOfSuperModules() {return fNumberOfSuperModules;}
  Float_t GetfPhiGapForSuperModules() {return fPhiGapForSM;}
  Float_t GetPhiModuleSize() const  {return fPhiModuleSize;}
  Float_t GetEtaModuleSize() const  {return fEtaModuleSize;}
  Float_t GetFrontSteelStrip() const {return fFrontSteelStrip;}
  Float_t GetLateralSteelStrip() const {return fLateralSteelStrip;}
  Float_t GetPassiveScintThick() const {return fPassiveScintThick;}
  Float_t GetPhiTileSize() const {return fPhiTileSize;}
  Float_t GetEtaTileSize() const {return fEtaTileSize;}
  Int_t   GetNPhiSuperModule() const {return fNPhiSuperModule;}
  Int_t   GetNPHIdiv() const {return fNPHIdiv ;}
  Int_t   GetNETAdiv() const {return fNETAdiv ;}
  Int_t   GetNCells()  const {return fNCells;}

  Int_t   GetNTRU() const    {return fNTRU ; }  
  Int_t   GetNTRUEta() const {return fNTRUEta ; }  
  Int_t   GetNTRUPhi() const {return fNTRUPhi ; }  

  Float_t GetSteelFrontThickness() const { return fSteelFrontThick;}
  Float_t GetLongModuleSize() const {return fLongModuleSize;}

  Float_t GetTrd1Angle() const {return fTrd1Angle;}
  Float_t Get2Trd1Dx2()  const {return f2Trd1Dx2;}
  Float_t GetTrd2AngleY()const {return fTrd2AngleY;}
  Float_t Get2Trd2Dy2()  const {return f2Trd2Dy2;}
  Float_t GetTubsR()     const {return fTubsR;}
  Float_t GetTubsTurnAngle() const {return fTubsTurnAngle;}
  // Dabs id <-> indexes; Shish-kebab case 
  Int_t   GetAbsCellId(Int_t nSupMod, Int_t nTower, Int_t nIphi, Int_t nIeta);
  Bool_t  GetCellIndex(Int_t absId, Int_t &nSupMod, Int_t &nTower, Int_t &nIphi, Int_t &nIeta);
  void    GetTowerPhiEtaIndexInSModule(Int_t nSupMod, Int_t nTower, Int_t &iphit, Int_t &ietat);
  void    GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nTower, Int_t nIphi, Int_t nIeta,
                                      Int_t &iphi, Int_t &ieta);
  Bool_t  CheckAbsCellId(Int_t ind); // replace the IsInECA
  // ---
  Float_t AngleFromEta(Float_t eta){ // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta){ // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }
  // These methods are obsolete but use in AliEMCALRecPoint - keep it now
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
  // --
  void SetNZ(Int_t nz) { fNZ= nz ; printf("SetNZ: Number of modules in Z set to %d", fNZ) ; }
  void SetNPhi(Int_t nphi) { fNPhi= nphi ; printf("SetNPhi: Number of modules in Phi set to %d", fNPhi) ; }

  void SetNTRU(Int_t ntru)    {fNTRU    = ntru; printf("SetNTRU: Number of TRUs per SuperModule set to %d", fNTRU) ; }
  void SetNTRUEta(Int_t ntru) {fNTRUEta = ntru; ; printf("SetNTRU: Number of TRUs per SuperModule in Etaset to %d", fNTRUEta) ;}
  void SetNTRUPhi(Int_t ntru) {fNTRUPhi = ntru; ; printf("SetNTRU: Number of TRUs per SuperModule in Phi set to %d", fNTRUPhi) ;}

  void SetSampling(Float_t samp) { fSampling = samp; printf("SetSampling: Sampling factor set to %f", fSampling) ; }

protected:
  AliEMCALGeometry(const Text_t* name, const Text_t* title="") :
    AliGeometry(name, title) {// ctor only for internal usage (singleton)
    Init();
  };
  AliEMCALGeometry() :
    AliGeometry() {// ctor only for internal usage (singleton)
    Init();
  };
  void Init(void);     			// initializes the parameters of EMCAL
  void CheckAditionalOptions();              //
  
private:
  static AliEMCALGeometry * fgGeom;	// pointer to the unique instance of the singleton
  static Bool_t fgInit;			// Tells if geometry has been succesfully set up.
  TObjArray *fArrayOpts;                //! array of geometry options

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
  Int_t   fNPhi;			// Number of Towers in the PHI direction
  Float_t fSampling;			// Sampling factor

  // Shish-kebab option - 23-aug-04 by PAI; COMPACT, TWIST, TRD1 and TRD2
  Int_t   fNumberOfSuperModules;         // default is 12 = 6 * 2 
  Float_t fSteelFrontThick;		 // Thickness of the front stell face of the support box - 9-sep-04
  Float_t fFrontSteelStrip;              // 13-may-05
  Float_t fLateralSteelStrip;            // 13-may-05
  Float_t fPassiveScintThick;            // 13-may-05
  Float_t fPhiModuleSize;                // Phi -> X 
  Float_t fEtaModuleSize;                // Eta -> Y
  Float_t fPhiTileSize;                  // 
  Float_t fEtaTileSize;                  // 
  Float_t fLongModuleSize;               // 
  Int_t   fNPhiSuperModule;              // 6 - number supermodule in phi direction
  Int_t   fNPHIdiv;                      // number phi divizion of module
  Int_t   fNETAdiv;                      // number eta divizion of module
  //
  Int_t   fNCells;                       // number of cells in calo
  Int_t   fNCellsInSupMod;               // number cell in super module
  Int_t   fNCellsInTower;                // number cell in tower(or module)
  //TRU parameters
  Int_t   fNTRU ;                        //! Number of TRUs per module
  Int_t   fNTRUEta ;                     //! Number of cell rows per Z in one TRU
  Int_t   fNTRUPhi ;                     //! Number of cell rows per Phi in one TRU
  // TRD1 options - 30-sep-04
  Float_t fTrd1Angle;                    // angle in x-z plane (in degree) 
  Float_t f2Trd1Dx2;                     // 2*dx2 for TRD1
  Float_t fPhiGapForSM;                  // Gap betweeen supermodules in phi direction
  Int_t   fKey110DEG;                    // for calculation abs cell id; 19-oct-05 
  // TRD2 options - 27-jan-07
  Float_t fTrd2AngleY;                   // angle in y-z plane (in degree) 
  Float_t f2Trd2Dy2;                     // 2*dy2 for TRD2
  Float_t fEmptySpace;                   // 2mm om fred drawing
  // Super module as TUBS
  Float_t fTubsR;                        // radius of tubs 
  Float_t fTubsTurnAngle;                // turn angle of tubs in degree
  // Service routine 
  static int ParseString(const TString &topt, TObjArray &Opt);

  ClassDef(AliEMCALGeometry,10) // EMCAL geometry class 
  };

#endif // AliEMCALGEOMETRY_H
