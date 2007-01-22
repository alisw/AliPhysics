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
class TVector3;
class TGeoMatrix;
class TParticle ; 
class AliEMCALShishKebabTrd1Module;
class AliEMCALRecPoint;
class TClonesArray ;

// --- AliRoot header files ---
#include <TArrayD.h>
#include <TMath.h>

#include "AliGeometry.h"

class AliEMCALGeometry : public AliGeometry {
public:
  AliEMCALGeometry(const AliEMCALGeometry& geom);
  virtual ~AliEMCALGeometry(void); 

  static AliEMCALGeometry * GetInstance(const Text_t* name,
					const Text_t* title="") ; 
  static AliEMCALGeometry * GetInstance();
  AliEMCALGeometry & operator = (const AliEMCALGeometry  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };
  static Char_t* GetDefaulGeometryName() {return fgDefaultGeometryName;}
  void PrintGeometry();                                           //*MENU*  
  void PrintCellIndexes(Int_t absId=0, int pri=0, char *tit="");  //*MENU*
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const;

  void FillTRU(const TClonesArray * digits, TClonesArray * amptru, TClonesArray * timeRtru)  ; //Fills Trigger Unit matrices with digit amplitudes and time
  void GetCellPhiEtaIndexInSModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru, Int_t &ietaSM, Int_t &iphiSM) const ; // Tranforms Eta-Phi Cell index in TRU into Eta-Phi index in Super Module
  
  // Have to call GetTransformationForSM() before calculation global charachteristics 
  void GetGlobal(const Double_t *loc, Double_t *glob, int ind) const;
  void GetGlobal(const TVector3 &vloc, TVector3 &vglob, int ind) const;
  void GetGlobal(Int_t absId, Double_t glob[3]) const;
  void GetGlobal(Int_t absId, TVector3 &vglob) const;
  // for a given tower index absId returns eta and phi of gravity center of tower.
  void EtaPhiFromIndex(Int_t absId, Double_t &eta, Double_t &phi) const;
  void EtaPhiFromIndex(Int_t absId, Float_t &eta, Float_t &phi) const;
  // 
  Bool_t GetPhiBoundariesOfSM   (Int_t nSupMod, Double_t &phiMin, Double_t &phiMax) const;
  Bool_t GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const;
  Bool_t SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const;

  Bool_t GetAbsCellIdFromEtaPhi(Double_t eta,Double_t phi, Int_t &absId) const;

  //  virtual void GetGlobal(const AliEMCALRecPoint *rp, TVector3 &vglob) const;

  virtual void GetGlobal(const AliRecPoint *rp, TVector3 &vglob) const;
  // Bool_t AreInSameTower(Int_t id1, Int_t id2) const ;  

  virtual void GetGlobal(const AliRecPoint *, TVector3 &, TMatrixF &) const {}

  virtual Bool_t Impact(const TParticle *) const {return kTRUE;}

  Bool_t IsInEMCAL(Double_t x, Double_t y, Double_t z) const;
  // General
  Bool_t  IsInitialized(void) const { return fgInit ; }
  // Return EMCAL geometrical parameters
  // geometry
  Char_t* GetNameOfEMCALEnvelope() const {return "XEN1";}
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
  //  Bool_t IsInECA(Int_t index) const { if ( (index > 0 && (index <= GetNZ() * GetNPhi()))) return kTRUE; else return kFALSE ;}

  Int_t   GetNumberOfSuperModules() const {return fNumberOfSuperModules;}
  Float_t GetfPhiGapForSuperModules() const {return fPhiGapForSM;}
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

  // TRD1 staff
  void    CreateListOfTrd1Modules();
  TList  *GetShishKebabTrd1Modules() const {return fShishKebabTrd1Modules;}
  AliEMCALShishKebabTrd1Module *GetShishKebabModule(Int_t neta);

  void     GetTransformationForSM();
  Float_t *GetSuperModulesPars() {return fParSM;}
  TGeoMatrix *GetTransformationForSM(int i) {
  if(i>=0 && i<GetNumberOfSuperModules()) return fMatrixOfSM[i]; 
                                        else return 0;}
  // May 31, 2006; ALICE numbering scheme: 
  // see ALICE-INT-2003-038: ALICE Coordinate System and Software Numbering Convention
  // All indexes are stared from zero now.
  // 
  // abs id <-> indexes; Shish-kebab case, only TRD1 now.
  // EMCAL -> Super Module -> module -> tower(or cell) - logic tree of EMCAL
  // 
  //**  Usual name of variable - Dec 18,2006 **
  //  nSupMod - index of super module (SM)
  //  nModule - index of module in SM
  //  nIphi   - phi index of tower(cell) in module
  //  nIeta   - eta index of tower(cell) in module
  //  
  //  Inside SM
  //  iphim   - phi index of module in SM  
  //  ietam   - eta index of module in SM  
  //
  //  iphi    - phi index of tower(cell) in SM  
  //  ieta    - eta index of tower(cell) in SM  
  Int_t   GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const;
  Bool_t  CheckAbsCellId(Int_t absId) const;
  Bool_t  GetCellIndex(Int_t absId, Int_t &nSupMod, Int_t &nModule, Int_t &nIphi, Int_t &nIeta) const;
  // Local coordinate of Super Module 
  void    GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t &iphim, Int_t &ietam) const;
  void    GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta,
                                      Int_t &iphi, Int_t &ieta) const ;
  Int_t   GetSuperModuleNumber(Int_t absId)  const;
  Int_t   GetNumberOfModuleInPhiDirection(Int_t nSupMod)  const
  {
    // inline function
    if(fKey110DEG == 1 && nSupMod>=10) return fNPhi/2;
    else                               return fNPhi;
  }
  // From cell indexes to abs cell id
  void    GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
	  Int_t &iphim, Int_t &ietam, Int_t &nModule) const;
  Int_t   GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const;
  // Methods for AliEMCALRecPoint - Feb 19, 2006
  Bool_t   RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t   RelPosCellInSModule(Int_t absId, Double_t loc[3]) const;
  Bool_t   RelPosCellInSModule(Int_t absId, TVector3 &vloc) const;
  // ---
  Float_t AngleFromEta(Float_t eta) const { // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta) const { // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }
  void SetNZ(Int_t nz) { fNZ= nz ; printf("SetNZ: Number of modules in Z set to %d", fNZ) ; }
  void SetNPhi(Int_t nphi) { fNPhi= nphi ; printf("SetNPhi: Number of modules in Phi set to %d", fNPhi) ; }

  void SetNTRU(Int_t ntru)    {fNTRU    = ntru; printf("SetNTRU: Number of TRUs per SuperModule set to %d", fNTRU) ; }
  void SetNTRUEta(Int_t ntru) {fNTRUEta = ntru; ; printf("SetNTRU: Number of TRUs per SuperModule in Etaset to %d", fNTRUEta) ;}
  void SetNTRUPhi(Int_t ntru) {fNTRUPhi = ntru; ; printf("SetNTRU: Number of TRUs per SuperModule in Phi set to %d", fNTRUPhi) ;}

  void SetSampling(Float_t samp) { fSampling = samp; printf("SetSampling: Sampling factor set to %f", fSampling) ; }

  Int_t GetNCellsInSupMod() const {return fNCellsInSupMod;}
  Int_t GetNCellsInModule()  const {return fNCellsInModule; }
  Int_t GetKey110DEG()      const {return fKey110DEG;}
  Int_t GetILOSS() const {return fILOSS;}
  Int_t GetIHADR() const {return fIHADR;}

  AliEMCALGeometry(); // default ctor only for internal usage (singleton)

protected:
  AliEMCALGeometry(const Text_t* name, const Text_t* title);// ctor only for internal usage (singleton)

  void Init(void);     			// initializes the parameters of EMCAL
  void CheckAdditionalOptions();        //
  void DefineSamplingFraction();        // Jun 5, 2006
  
private:
  static AliEMCALGeometry * fgGeom;	// pointer to the unique instance of the singleton
  static Bool_t  fgInit;	        // Tells if geometry has been succesfully set up.
  static Char_t* fgDefaultGeometryName; // Default name of geometry

  TString fGeoName;                     //geometry name

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
  Float_t fPhiTileSize;                  // Size of phi tile
  Float_t fEtaTileSize;                  // Size of eta tile
  Float_t fLongModuleSize;               // Size of long module
  Int_t   fNPhiSuperModule;              // 6 - number supermodule in phi direction
  Int_t   fNPHIdiv;                      // number phi divizion of module
  Int_t   fNETAdiv;                      // number eta divizion of module
  //
  Int_t   fNCells;                       // number of cells in calo
  Int_t   fNCellsInSupMod;               // number cell in super module
  Int_t   fNCellsInModule;               // number cell in module)
  //TRU parameters
  Int_t   fNTRU ;                        //! Number of TRUs per module
  Int_t   fNTRUEta ;                     //! Number of cell rows per Z in one TRU
  Int_t   fNTRUPhi ;                     //! Number of cell rows per Phi in one TRU
  // TRD1 options - 30-sep-04
  Float_t fTrd1Angle;                    // angle in x-z plane (in degree) 
  Float_t f2Trd1Dx2;                     // 2*dx2 for TRD1
  Float_t fPhiGapForSM;                  // Gap betweeen supermodules in phi direction
  Int_t   fKey110DEG;                    // for calculation abs cell id; 19-oct-05 
  TArrayD fPhiBoundariesOfSM;            // phi boundaries of SM in rad; size is fNumberOfSuperModules;
  TArrayD fPhiCentersOfSM;                // phi of centers of SMl size is fNumberOfSuperModules/2
  Float_t fEtaMaxOfTRD1;                 // max eta in case of TRD1 geometry (see AliEMCALShishKebabTrd1Module)
  // TRD2 options - 27-jan-07
  Float_t fTrd2AngleY;                   // angle in y-z plane (in degree) 
  Float_t f2Trd2Dy2;                     // 2*dy2 for TRD2
  Float_t fEmptySpace;                   // 2mm om fred drawing
  // Super module as TUBS
  Float_t fTubsR;                        // radius of tubs 
  Float_t fTubsTurnAngle;                // turn angle of tubs in degree
  // Local Coordinates of SM
  TArrayD  fCentersOfCellsEtaDir;        // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD  fCentersOfCellsXDir;          // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  TArrayD  fCentersOfCellsPhiDir;        // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  //
  TArrayD  fEtaCentersOfCells;           // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  TArrayD  fPhiCentersOfCells;           // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
  // Move from AliEMCALv0 - Feb 19, 2006
  TList *fShishKebabTrd1Modules; //! list of modules
  // Local coordinates of SM for TRD1
  Float_t     fParSM[3];       // SM sizes as in GEANT (TRD1)
  TGeoMatrix* fMatrixOfSM[12]; //![fNumberOfSuperModules]; get from gGeoManager;

  char *fAdditionalOpts[6];  //! some additional options for the geometry type and name
  int  fNAdditionalOpts;     //! size of additional options parameter

  // Options for Geant (MIP business) - will call in AliEMCAL
  Int_t fILOSS;
  Int_t fIHADR;

  ClassDef(AliEMCALGeometry, 11) // EMCAL geometry class 
  };

#endif // AliEMCALGEOMETRY_H
