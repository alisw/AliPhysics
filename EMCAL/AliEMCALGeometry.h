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
//*--   and : Alexei Pavlinov (WSU) - shashlyk staff
//*--   and : Gustavo Conesa: Add TRU mapping. TRU parameters still not fixed.
//*--   and : Magali Estienne : analysis access adaptations

// --- ROOT system ---
#include <TNamed.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TVector3.h>
#include <TGeoMatrix.h> 
class TBrowser ;
class TParticle ;

// --- AliRoot header files ---
#include "AliEMCALEMCGeometry.h"
#include "AliEMCALGeoParams.h"
class AliEMCALShishKebabTrd1Module;
#include "AliLog.h"

class AliEMCALGeometry : public TNamed {

public: 

  AliEMCALGeometry();
  AliEMCALGeometry(const Text_t* name, const Text_t* title="",
                   const Text_t* mcname="", const Text_t* mctitle="");
  AliEMCALGeometry(const AliEMCALGeometry & geom);
  
  virtual ~AliEMCALGeometry(void); 
  AliEMCALGeometry & operator = (const AliEMCALGeometry  & rvalue);
  
  static AliEMCALGeometry * GetInstance(const Text_t* name,      const Text_t* title="",
                                        const Text_t* mcname="TGeant3", const Text_t* mctitle="") ; 
  static AliEMCALGeometry * GetInstance();


  //////////
  // General
  //
  static Bool_t  IsInitialized(void)            {return AliEMCALEMCGeometry::fgInit; }
  static const Char_t* GetDefaultGeometryName() {return AliEMCALEMCGeometry::fgkDefaultGeometryName;}
  
  /////////////
  // TRD1 stuff
  void    CreateListOfTrd1Modules();
  TList  *GetShishKebabTrd1Modules() const {return fShishKebabTrd1Modules;}
  AliEMCALShishKebabTrd1Module *GetShishKebabModule(Int_t neta) const;

  void PrintGeometryGeoUtils();   // *MENU*  
  void PrintCellIndexes(Int_t absId=0, int pri=0, const char *tit="") const ;  //*MENU*
  void PrintLocalTrd1(Int_t pri=0) const;  // *MENU*
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const;

  virtual Bool_t Impact(const TParticle *) const;
  void ImpactOnEmcal(TVector3 vtx, Double_t theta, Double_t phi, Int_t & absId, TVector3 & vimpact) const;
  Bool_t IsInEMCAL(Double_t x, Double_t y, Double_t z) const;

  //////////////////////////////////////
  // Return EMCAL geometrical parameters
  //
  
  AliEMCALEMCGeometry* GetEMCGeometry()       const { return fEMCGeometry                            ; }
  //
  const Char_t*  GetNameOfEMCALEnvelope(void) const { return fEMCGeometry->GetNameOfEMCALEnvelope()  ; }
  Float_t  GetArm1PhiMin(void)                const { return fEMCGeometry->GetArm1PhiMin()           ; }
  Float_t  GetArm1PhiMax(void)                const { return fEMCGeometry->GetArm1PhiMax()           ; }
  Float_t  GetArm1EtaMin(void)                const { return fEMCGeometry->GetArm1EtaMin()           ; }
  Float_t  GetArm1EtaMax(void)                const { return fEMCGeometry->GetArm1EtaMax()           ; }
  Float_t  GetIPDistance(void)                const { return fEMCGeometry->GetIPDistance()           ; }   
  Float_t  GetEnvelop(Int_t index)            const { return fEMCGeometry->GetEnvelop(index)         ; }  
  Float_t  GetShellThickness(void)            const { return fEMCGeometry->GetShellThickness()       ; }
  Float_t  GetZLength(void)                   const { return fEMCGeometry->GetZLength()              ; } 
  Int_t    GetNECLayers(void)                 const { return fEMCGeometry->GetNECLayers()            ; }
  Int_t    GetNZ(void)                        const { return fEMCGeometry->GetNZ()                   ; }
  Int_t    GetNEta(void)                      const { return fEMCGeometry->GetNEta()                 ; }
  Int_t    GetNPhi(void)                      const { return fEMCGeometry->GetNPhi()                 ; }
  Float_t  GetECPbRadThick(void)              const { return fEMCGeometry->GetECPbRadThick()         ; }
  Float_t  GetECScintThick(void)              const { return fEMCGeometry->GetECScintThick()         ; }
  Float_t  GetSampling(void)                  const { return fEMCGeometry->GetSampling()             ; } 
  Int_t    GetNumberOfSuperModules(void)      const { return fEMCGeometry->GetNumberOfSuperModules() ; }
  Float_t  GetPhiGapForSuperModules(void)     const { return fEMCGeometry->GetfPhiGapForSuperModules() ; }
  Float_t  GetPhiModuleSize(void)             const { return fEMCGeometry->GetPhiModuleSize()        ; }
  Float_t  GetEtaModuleSize(void)             const { return fEMCGeometry->GetEtaModuleSize()        ; }
  Float_t  GetFrontSteelStrip(void)           const { return fEMCGeometry->GetFrontSteelStrip()      ; }
  Float_t  GetLateralSteelStrip(void)         const { return fEMCGeometry->GetLateralSteelStrip()    ; }
  Float_t  GetPassiveScintThick(void)         const { return fEMCGeometry->GetPassiveScintThick()    ; }
  Float_t  GetPhiTileSize(void)               const { return fEMCGeometry->GetPhiTileSize()          ; }
  Float_t  GetEtaTileSize(void)               const { return fEMCGeometry->GetEtaTileSize()          ; }
  Int_t    GetNPhiSuperModule(void)           const { return fEMCGeometry->GetNPhiSuperModule()      ; }
  Int_t    GetNPHIdiv(void)                   const { return fEMCGeometry->GetNPHIdiv()              ; }
  Int_t    GetNETAdiv(void)                   const { return fEMCGeometry->GetNETAdiv()              ; }
  Int_t    GetNCells(void)                    const { return fEMCGeometry->GetNCells()               ; }
  Float_t  GetLongModuleSize(void)            const { return fEMCGeometry->GetLongModuleSize()       ; }
  Float_t  GetTrd1Angle(void)                 const { return fEMCGeometry->GetTrd1Angle()            ; }
  Float_t  Get2Trd1Dx2(void)                  const { return fEMCGeometry->Get2Trd1Dx2()             ; }
  Float_t  GetTrd1AlFrontThick()              const { return fEMCGeometry->GetTrd1AlFrontThick()     ; }
  Float_t  GetTrd1BondPaperThick()            const { return fEMCGeometry->GetTrd1BondPaperThick()   ; }
  // --
  Int_t    GetNCellsInSupMod(void)            const { return fEMCGeometry->GetNCellsInSupMod()       ; }
  Int_t    GetNCellsInModule(void)            const { return fEMCGeometry->GetNCellsInModule()       ; }
  Int_t    GetKey110DEG(void)                 const { return fEMCGeometry->GetKey110DEG()            ; }
  Int_t    GetILOSS(void)                     const { return fEMCGeometry->GetILOSS()                ; }
  Int_t    GetIHADR(void)                     const { return fEMCGeometry->GetIHADR()                ; }  
  // --
  Float_t  GetDeltaEta(void)                  const { return fEMCGeometry->GetDeltaEta()             ; }
  Float_t  GetDeltaPhi(void)                  const { return fEMCGeometry->GetDeltaPhi()             ; }
  Int_t    GetNTowers(void)                   const { return fEMCGeometry->GetNTowers()              ; }
  //
  Double_t GetPhiCenterOfSM(Int_t nsupmod)    const { return fEMCGeometry->GetPhiCenterOfSM(nsupmod) ; }
  Float_t  GetSuperModulesPar(Int_t ipar)     const { return fEMCGeometry->GetSuperModulesPar(ipar)  ; }
  //
  Bool_t   GetPhiBoundariesOfSM(Int_t nSupMod, Double_t &phiMin, Double_t &phiMax)    const 
    { return fEMCGeometry->GetPhiBoundariesOfSM(nSupMod, phiMin, phiMax)   ; }
  Bool_t   GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const 
    { return fEMCGeometry->GetPhiBoundariesOfSMGap(nPhiSec, phiMin, phiMax); }
  //
  
  //////////////////////////////////////////////////
  // Obsolete methods to be thrown out when feasible
  Float_t GetGap2Active(void)                 const { return fEMCGeometry->GetGap2Active()           ; }
  Float_t GetSteelFrontThickness(void)        const { return fEMCGeometry->GetSteelFrontThickness()  ; }
  Float_t GetTrd2AngleY(void)                 const { return fEMCGeometry->GetTrd2AngleY()           ; }
  Float_t Get2Trd2Dy2(void)                   const { return fEMCGeometry->Get2Trd2Dy2()             ; }
  Float_t GetTubsR(void)                      const { return fEMCGeometry->GetTubsR()                ; }
  Float_t GetTubsTurnAngle(void)              const { return fEMCGeometry->GetTubsTurnAngle()        ; }
  //Float_t GetAlFrontThickness(void)           const { return fEMCGeometry->GetAlFrontThickness()     ; }
  //Float_t GetIP2ECASection(void)              const { return fEMCGeometry->GetIP2ECASection()        ; }   
  //////////////////////////////////////////////////
  
  ///////////////////////////////
  //Geometry data member setters
  //
  void SetNZ(Int_t nz)           { fEMCGeometry->SetNZ(nz)         ; }
  void SetNPhi(Int_t nphi)       { fEMCGeometry->SetNPhi(nphi)     ; }
  //Trigger
  void SetNTRUEta(Int_t ntru)    { fEMCGeometry->SetNTRUEta(ntru)  ; }
  void SetNTRUPhi(Int_t ntru)    { fEMCGeometry->SetNTRUPhi(ntru)  ; }
  //
  void SetSampling(Float_t samp) { fEMCGeometry->SetSampling(samp) ; }
  //
  void PrintGeometry()           { fEMCGeometry->PrintGeometry()   ; }  //*MENU*  
  
  //////////////////////////
  // Global geometry methods
  //
  void GetGlobal(const Double_t *loc, Double_t *glob, int ind) const;
  void GetGlobal(const TVector3 &vloc, TVector3 &vglob, int ind) const;
  void GetGlobal(Int_t absId, Double_t glob[3]) const;
  void GetGlobal(Int_t absId, TVector3 &vglob) const;

  ////////////////////////////////////////
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
  //
  // for a given tower index absId returns eta and phi of gravity center of tower.
  void    EtaPhiFromIndex(Int_t absId, Double_t &eta, Double_t &phi) const;
  void    EtaPhiFromIndex(Int_t absId, Float_t  &eta, Float_t  &phi) const;

  Bool_t  GetAbsCellIdFromEtaPhi(Double_t eta,Double_t phi, Int_t &absId) const;
  Bool_t  SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const;
  Int_t   GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const;
  Bool_t  CheckAbsCellId(Int_t absId) const;
  Bool_t  GetCellIndex(Int_t absId, Int_t &nSupMod, Int_t &nModule, Int_t &nIphi, 
		       Int_t &nIeta) const;
  // Local coordinate of Super Module 
  void    GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t &iphim, 
					Int_t &ietam) const;
  void    GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta,
                                      Int_t &iphi, Int_t &ieta) const ;
  Int_t   GetSuperModuleNumber(Int_t absId)  const;
  Int_t   GetNumberOfModuleInPhiDirection(Int_t nSupMod)  const
  { 
    if(fKey110DEG == 1 && nSupMod>=10) return fNPhi/2;
    else                               return fNPhi;
  } 
  // From cell indexes to abs cell id
  void    GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
					      Int_t &iphim, Int_t &ietam, Int_t &nModule) const;
  Int_t   GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const;

  // Methods for AliEMCALRecPoint - Feb 19, 2006
  Bool_t  RelPosCellInSModule(Int_t absId, 
                              Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t distEff,
                              Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t loc[3]) const;
  Bool_t  RelPosCellInSModule(Int_t absId, TVector3 &vloc)  const;

  // Local Coordinates of SM
  TArrayD  GetCentersOfCellsEtaDir() const { return fCentersOfCellsEtaDir ; }     // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD  GetCentersOfCellsXDir()   const { return fCentersOfCellsXDir   ; }     // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  TArrayD  GetCentersOfCellsPhiDir() const { return fCentersOfCellsPhiDir ; }     // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  //
  TArrayD  GetEtaCentersOfCells()    const { return fEtaCentersOfCells    ; }     // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  TArrayD  GetPhiCentersOfCells()    const { return fPhiCentersOfCells    ; }     // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)

  
  // For gamma(Jet) trigger simulations *FIXME OLD TO BE REMOVED with AliEMCALTrigger*
  Int_t    GetNTRU()             const { return fEMCGeometry->GetNTRU()             ; }  
  Int_t    GetNTRUEta()          const { return fEMCGeometry->GetNTRUEta()          ; }  
  Int_t    GetNTRUPhi()          const { return fEMCGeometry->GetNTRUPhi()          ; }
  Int_t    GetNEtaSubOfTRU()     const { return fEMCGeometry->GetNEtaSubOfTRU()     ; }
  Int_t    GetNModulesInTRU()    const { return fEMCGeometry->GetNModulesInTRU()    ; }
  Int_t    GetNModulesInTRUEta() const { return fEMCGeometry->GetNModulesInTRUEta() ; }  
  Int_t    GetNModulesInTRUPhi() const { return fEMCGeometry->GetNModulesInTRUPhi() ; }
  // *MEFIX OLD TO BE REMOVED*

  //
  // Tranforms Eta-Phi Module index in TRU into Eta-Phi index in Super Module
  void     GetModulePhiEtaIndexInSModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru, 
                                                     Int_t &ietaSM, Int_t &iphiSM) const;
  Int_t   GetAbsTRUNumberFromNumberInSm(const Int_t row, const Int_t col, const Int_t sm) const ;

	
  void     BuildFastOR2DMap();
  Bool_t   GetAbsFastORIndexFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& id) const;
  Bool_t                    GetAbsFastORIndexFromPositionInTRU(const Int_t iTRU, const Int_t iEta, const Int_t iPhi, Int_t& id) const;	
  Bool_t                    GetAbsFastORIndexFromPositionInSM( const Int_t  iSM, const Int_t iEta, const Int_t iPhi, Int_t& id) const;	
  Bool_t                    GetAbsFastORIndexFromPositionInEMCAL(                const Int_t iEta, const Int_t iPhi, Int_t& id) const;
  Bool_t             GetTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iADC) const;
  Bool_t   GetPositionInTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const;
  Bool_t    GetPositionInSMFromAbsFastORIndex(const Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const;
  Bool_t GetPositionInEMCALFromAbsFastORIndex(const Int_t id, Int_t& iEta, Int_t& iPhi) const;
  Bool_t          GetFastORIndexFromCellIndex(const Int_t id, Int_t& idx) const;
  Bool_t          GetCellIndexFromFastORIndex(const Int_t id, Int_t idx[4]) const;
  Bool_t              GetTRUIndexFromSTUIndex(const Int_t id, Int_t& idx) const;
  Int_t               GetTRUIndexFromSTUIndex(const Int_t id) const;
  Bool_t            GetFastORIndexFromL0Index(const Int_t iTRU, const Int_t id, Int_t idx[], const Int_t size) const;
	
  ///////////////////
  // useful utilities
  //
  Float_t AngleFromEta(Float_t eta)        const { // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta) const { // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }

  //Method to set shift-rotational matrixes from ESDHeader
  void SetMisalMatrix(const TGeoHMatrix * m, Int_t smod) {
    fUseExternalMatrices = kTRUE;
	  if (smod >= 0 && smod < fEMCGeometry->GetNumberOfSuperModules()){
            if(!fkSModuleMatrix[smod]) fkSModuleMatrix[smod] = new TGeoHMatrix(*m) ; //Set only if not set yet
    }
	  else AliFatal(Form("Wrong supermodule index -> %d",smod));
  }
	
  //Alternate geometry that allows to calculate tower position for different particles and different alignments
  void RecalculateTowerPosition(Float_t drow, Float_t dcol, const Int_t sm, const Float_t depth,
                                const Float_t misaligTransShifts[15], const Float_t misaligRotShifts[15],Float_t global[3]) const;
  
  //Returns shift-rotational matrixes for different volumes
  const TGeoHMatrix * GetMatrixForSuperModule(Int_t smod)const ;
	
protected:

  void Init(void);     		     // initializes the parameters of EMCAL
  
  AliEMCALEMCGeometry * fEMCGeometry;// Geometry object for Electromagnetic calorimeter

  TString  fGeoName;                 // geometry name
  Int_t    fKey110DEG;               // for calculation abs cell id; 19-oct-05 
  Int_t    fNCellsInSupMod;          // number cell in super module
  Int_t    fNETAdiv;                 // number eta divizion of module
  Int_t    fNPHIdiv;                 // number phi divizion of module
  Int_t    fNCellsInModule;          // number cell in module
  TArrayD  fPhiBoundariesOfSM;       // phi boundaries of SM in rad; size is fNumberOfSuperModules;
  TArrayD  fPhiCentersOfSM;          // phi of centers of SMl size is fNumberOfSuperModules/2
  // Local Coordinates of SM
  TArrayD  fPhiCentersOfCells;       // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
  TArrayD  fCentersOfCellsEtaDir;    // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD  fCentersOfCellsPhiDir;    // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  TArrayD  fEtaCentersOfCells;       // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  Int_t    fNCells;                  // number of cells in calo
  Int_t    fNPhi;		                 // Number of Towers in the PHI direction
  TArrayD  fCentersOfCellsXDir;      // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  Float_t  fEnvelop[3];              // the GEANT TUB for the detector 
  Float_t  fArm1EtaMin;              // Minimum pseudorapidity position of EMCAL in Eta
  Float_t  fArm1EtaMax;              // Maximum pseudorapidity position of EMCAL in Eta
  Float_t  fArm1PhiMin;              // Minimum angular position of EMCAL in Phi (degrees)
  Float_t  fArm1PhiMax;              // Maximum angular position of EMCAL in Phi (degrees)
  Float_t  fEtaMaxOfTRD1;            // Max eta in case of TRD1 geometry (see AliEMCALShishKebabTrd1Module)
  TList   *fShishKebabTrd1Modules;   // list of modules
  Float_t  fParSM[3];                // SM sizes as in GEANT (TRD1)
  Float_t  fPhiModuleSize;           // Phi -> X 
  Float_t  fEtaModuleSize;           // Eta -> Y 
  Float_t  fPhiTileSize;             // Size of phi tile
  Float_t  fEtaTileSize;             // Size of eta tile
  Int_t    fNZ;                      // Number of Towers in the Z direction
  Float_t  fIPDistance;		     // Radial Distance of the inner surface of the EMCAL
  Float_t  fLongModuleSize;          // Size of long module
  // Geometry Parameters
  Float_t  fShellThickness;	     // Total thickness in (x,y) direction
  Float_t  fZLength;		     // Total length in z direction
  Float_t  fSampling;		     // Sampling factor

  Int_t    fFastOR2DMap[48][64];     // FastOR 2D Map over full EMCal
	
  TGeoHMatrix* fkSModuleMatrix[AliEMCALGeoParams::fgkEMCALModules] ; //Orientations of EMCAL super modules
  Bool_t   fUseExternalMatrices;      // Use the matrices set in fkSModuleMatrix and not those in the geoManager
	
private:
  
  static AliEMCALGeometry *fgGeom;	               // Pointer to the unique instance of the singleton
  static Bool_t            fgInit;	               // Tells if geometry has been succesfully set up.
  static const Char_t     *fgkDefaultGeometryName; // Default name of geometry
  
  
  ClassDef(AliEMCALGeometry,16)       // EMCAL geometry class 

} ;

#endif // AliEMCALGEOUTILS_H

