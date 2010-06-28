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
#include <Riostream.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TMatrixF.h>
class TVector3;

// --- AliRoot header files ---
#include "AliEMCALGeoUtils.h"
#include "AliEMCALEMCGeometry.h"
//class AliRecPoint;
//class AliEMCALRecPoint;

class AliEMCALGeometry : public AliEMCALGeoUtils {

public:

  AliEMCALGeometry(); // default ctor only for internal usage (singleton)
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


  //////////
  // General
  //
  Bool_t  IsInitialized(void) const { return AliEMCALEMCGeometry::fgInit ; }
  static const Char_t* GetDefaultGeometryName() {return AliEMCALEMCGeometry::fgkDefaultGeometryName;}

  //////////////////////////
  // Global geometry methods
  //
  using AliEMCALGeoUtils::GetGlobal;

//  virtual void GetGlobal(const AliRecPoint *rp, TVector3 & gpos, TMatrixF & /* gmat */) 
//    const {GetGlobal(rp,gpos); }
//  virtual void GetGlobalEMCAL(const AliEMCALRecPoint *rp, TVector3 &vglob) const;
//  virtual void GetGlobalEMCAL(const AliEMCALRecPoint *rp, TVector3 & gpos, TMatrixF & /* gmat */) 
//    const {GetGlobalEMCAL(rp,gpos); }

  // Return EMC geometry parameters
  AliEMCALEMCGeometry * GetEMCGeometry()      const {return fEMCGeometry ;}

  //////////////////////////////////////
  // Return EMCAL geometrical parameters
  //
  const Char_t*  GetNameOfEMCALEnvelope(void) const {return fEMCGeometry->GetNameOfEMCALEnvelope();}
  Float_t  GetArm1PhiMin(void) const { return fEMCGeometry->GetArm1PhiMin(); }
  Float_t  GetArm1PhiMax(void) const { return fEMCGeometry->GetArm1PhiMax(); }
  Float_t  GetArm1EtaMin(void) const { return fEMCGeometry->GetArm1EtaMin();}
  Float_t  GetArm1EtaMax(void) const { return fEMCGeometry->GetArm1EtaMax();}
  Float_t  GetIPDistance(void) const { return fEMCGeometry->GetIPDistance();}   
  Float_t  GetEnvelop(Int_t index) const { return fEMCGeometry->GetEnvelop(index); }  
  Float_t  GetShellThickness(void) const { return fEMCGeometry->GetShellThickness(); }
  Float_t  GetZLength(void) const { return fEMCGeometry->GetZLength(); } 
  Int_t    GetNECLayers(void) const {return fEMCGeometry->GetNECLayers();}
  Int_t    GetNZ(void) const {return fEMCGeometry->GetNZ();}
  Int_t    GetNEta(void) const {return fEMCGeometry->GetNEta();}
  Int_t    GetNPhi(void) const {return fEMCGeometry->GetNPhi();}
  Float_t  GetECPbRadThick(void)const {return fEMCGeometry->GetECPbRadThick();}
  Float_t  GetECScintThick(void) const {return fEMCGeometry->GetECScintThick();}
  Float_t  GetSampling(void) const {return fEMCGeometry->GetSampling(); } 
  Int_t    GetNumberOfSuperModules(void) const {return fEMCGeometry->GetNumberOfSuperModules();}
  Float_t  GetfPhiGapForSuperModules(void) const {return fEMCGeometry->GetfPhiGapForSuperModules();}
  Float_t  GetPhiModuleSize(void) const  {return fEMCGeometry->GetPhiModuleSize();}
  Float_t  GetEtaModuleSize(void) const  {return fEMCGeometry->GetEtaModuleSize();}
  Float_t  GetFrontSteelStrip(void) const {return fEMCGeometry->GetFrontSteelStrip();}
  Float_t  GetLateralSteelStrip(void) const {return fEMCGeometry->GetLateralSteelStrip();}
  Float_t  GetPassiveScintThick(void) const {return fEMCGeometry->GetPassiveScintThick();}
  Float_t  GetPhiTileSize(void) const {return fEMCGeometry->GetPhiTileSize();}
  Float_t  GetEtaTileSize(void) const {return fEMCGeometry->GetEtaTileSize();}
  Int_t    GetNPhiSuperModule(void) const {return fEMCGeometry->GetNPhiSuperModule();}
  Int_t    GetNPHIdiv(void) const {return fEMCGeometry->GetNPHIdiv();}
  Int_t    GetNETAdiv(void) const {return fEMCGeometry->GetNETAdiv();}
  Int_t    GetNCells(void)  const {return fEMCGeometry->GetNCells();}
  Float_t  GetLongModuleSize(void) const {return fEMCGeometry->GetLongModuleSize();}
  Float_t  GetTrd1Angle(void) const {return fEMCGeometry->GetTrd1Angle();}
  Float_t  Get2Trd1Dx2(void)  const {return fEMCGeometry->Get2Trd1Dx2();}
  // --
  Int_t    GetNCellsInSupMod(void) const {return fEMCGeometry->GetNCellsInSupMod();}
  Int_t    GetNCellsInModule(void)  const {return fEMCGeometry->GetNCellsInModule(); }
  Int_t    GetKey110DEG(void)      const {return fEMCGeometry->GetKey110DEG();}
  Int_t    GetILOSS(void) const {return fEMCGeometry->GetILOSS();}
  Int_t    GetIHADR(void) const {return fEMCGeometry->GetIHADR();}
  // For gamma(Jet) trigger simulations
  Int_t    GetNTRU() const    {return fEMCGeometry->GetNTRU(); }  
  Int_t    GetNTRUEta() const {return fEMCGeometry->GetNTRUEta(); }  
  Int_t    GetNTRUPhi() const {return fEMCGeometry->GetNTRUPhi(); }
  Int_t    GetNEtaSubOfTRU() const {return fEMCGeometry->GetNEtaSubOfTRU();}
  Int_t    GetNModulesInTRU() const {return fEMCGeometry->GetNModulesInTRU(); }
  Int_t    GetNModulesInTRUEta() const {return fEMCGeometry->GetNModulesInTRUEta(); }  
  Int_t    GetNModulesInTRUPhi() const {return fEMCGeometry->GetNModulesInTRUPhi(); }  

  // --
  Float_t  GetDeltaEta(void) const {return fEMCGeometry->GetDeltaEta();}
  Float_t  GetDeltaPhi(void) const {return fEMCGeometry->GetDeltaPhi();}
  Int_t    GetNTowers(void) const {return fEMCGeometry->GetNTowers();}
  //
  Double_t GetPhiCenterOfSM(Int_t nsupmod) const {return fEMCGeometry->GetPhiCenterOfSM(nsupmod);}
  Float_t *GetSuperModulesPars(void) const {return fEMCGeometry->GetSuperModulesPars();}
  //
  Bool_t   GetPhiBoundariesOfSM(Int_t nSupMod, Double_t &phiMin, Double_t &phiMax) const {return fEMCGeometry->GetPhiBoundariesOfSM(nSupMod, phiMin, phiMax);}
  Bool_t   GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const {return fEMCGeometry->GetPhiBoundariesOfSMGap(nPhiSec, phiMin, phiMax);}
  //
  
  //  Methods for AliEMCALRecPoint with taking into account energy of rec.point - Jul 30. 2007
  using AliEMCALGeoUtils::RelPosCellInSModule;
  Bool_t   RelPosCellInSModule(Int_t absId,Double_t distEff,Double_t &xr,Double_t &yr,
							   Double_t & zr) const;
	
  //Not in use, comment for the moment
  //Bool_t   RelPosCellInSModule(Int_t absId,Int_t maxAbsId,Double_t distEff,Double_t &xr,
  //		       Double_t &yr,Double_t &zr) const;

  ///////////////////////////////
  //Geometry data member setters
  //
  void SetNZ(Int_t nz) { fEMCGeometry->SetNZ(nz);}
  void SetNPhi(Int_t nphi) { fEMCGeometry->SetNPhi(nphi);}

  void SetNTRUEta(Int_t ntru) { fEMCGeometry->SetNTRUEta(ntru);}
  void SetNTRUPhi(Int_t ntru) { fEMCGeometry->SetNTRUPhi(ntru);}
  void SetSampling(Float_t samp) { fEMCGeometry->SetSampling(samp);}
  
/*   /////////////////// */
/*   // useful utilities */
/*   // */
/*   Float_t AngleFromEta(Float_t eta) const { // returns theta in radians for a given pseudorapidity */
/*     return 2.0*TMath::ATan(TMath::Exp(-eta)); */
/*   } */
/*   Float_t ZFromEtaR(Float_t r,Float_t eta) const { // returns z in for a given */
/*     // pseudorapidity and r=sqrt(x*x+y*y). */
/*     return r/TMath::Tan(AngleFromEta(eta)); */
/*   } */

  //////////////////////////////////////////////////
  // Obsolete methods to be thrown out when feasible
  Float_t GetAlFrontThickness(void) const { return fEMCGeometry->GetAlFrontThickness();}
  Float_t GetGap2Active(void) const {return  fEMCGeometry->GetGap2Active();}
  Float_t GetSteelFrontThickness(void) const { return fEMCGeometry->GetSteelFrontThickness();}
  Float_t GetTrd2AngleY(void) const {return fEMCGeometry->GetTrd2AngleY();}
  Float_t Get2Trd2Dy2(void)  const {return fEMCGeometry->Get2Trd2Dy2();}
  Float_t GetTubsR(void)     const {return fEMCGeometry->GetTubsR();}
  Float_t GetTubsTurnAngle(void) const {return fEMCGeometry->GetTubsTurnAngle();}
  Float_t GetIP2ECASection(void) const { return fEMCGeometry->GetIP2ECASection(); }   
  //////////////////////////////////////////////////

protected:

  // ctor only for internal usage (singleton)
  AliEMCALGeometry(const Text_t* name, const Text_t* title);

  void Init(void);     			// initializes the parameters of EMCAL
  
private:

  //Member data
  static AliEMCALGeometry * fgGeom;	// pointer to the unique instance of the singleton
  //  static Bool_t  fgInit;	        // Tells if geometry has been succesfully set up.
  static const Char_t* fgkDefaultGeometryName; // Default name of geometry

  ///////////////////////////////////////////////////////////

  ClassDef(AliEMCALGeometry, 15) // EMCAL geometry class 
};

#endif // AliEMCALGEOMETRY_H
