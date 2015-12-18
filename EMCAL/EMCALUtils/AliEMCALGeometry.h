#ifndef ALIEMCALGEOMETRY_H
#define ALIEMCALGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALGeometry
/// \brief EMCal geometry, singleton
///
/// Geometry class  for EMCAL : singleton
/// EMCAL consists of layers of scintillator and lead
/// with scintillator fiber arranged as "shish-kebab" skewers
/// Places the the Barrel Geometry of The EMCAL at Midrapidity
/// between 80 and 180(or 190) degrees of Phi and
/// -0.7 to 0.7 in eta
///
///    * EMCAL geometry tree:
///    * EMCAL -> superModule -> module -> tower(cell)
///    * Indexes
///    *  absId -> nSupMod     -> nModule -> (nIphi,nIeta)
///
///   Name choices:
///   * EMCAL_PDC06 (geometry used for PDC06 simulations, kept for backward compatibility)
///      * = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG in old notation
///   * EMCAL_COMPLETE (geometry for expected complete detector)
///      * = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG scTh=0.176 pbTh=0.144
///          in old notation
///   * EMCAL_FIRSTYEARV1 - geometry for December 2009 to December 2010 run period;
///      *          fixed bug for positions of modules inside SM
///                 (first module has tilt 0.75 degree);
///                 the sizes updated with last information from production
///                 drawing (end of October 2010).
///
///   * EMCAL_COMPLETEV1: Same fixes as FIRSTYEAR and 10 SM instead of 10 + 2 one_third SM, for 2011 runs
///
///   * EMCAL_COMPLETE12SMV1: contains 12 SM for runs from year 2012 and on
///
///   * EMCAL_COMPLETE12SMV1_DCAL: contains 12 SM and 6 DCAL SM
///
///   * EMCAL_COMPLETE12SMV1_DCAL_8SM: contains 12 SM and 8 DCAL SM including the DCAL extention (2 SM)
///
///   * EMCAL_COMPLETE12SMV1_DCAL_DEV: contains 12 SM shifted and 10 DCAL SM
///
///   * EMCAL_WSUC (Wayne State test stand)
///      * = no definite equivalent in old notation, was only used by
///          Aleksei, but kept for testing purposes
///
///
/// Usage:
///         You can create the AliEMCALGeometry object independently from anything.
///         You have to use just the correct name of geometry. If name is empty string the
///         default name of geometry will be used.
///
///  AliEMCALGeometry* g = AliEMCALGeometry::GetInstance(name,title); // first time
///  ..
///  g = AliEMCALGeometry::GetInstance();                             // after first time
///
///  MC:   If you work with MC data you have to get geometry the next way:
///  ==                                      =============================
///  AliRunLoader    *rl   = AliRunLoader::Instance();
///  AliEMCALGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
///  TGeoManager::Import("geometry.root");
///
/// \author Sahal Yacoob (LBL / UCT)
/// \author Yves Schutz (SUBATECH)
/// \author Jennifer Klay (LBL)
/// \author Alexei Pavlinov (WSU)
///
///  Implementation for analysis usage, before AliEMCALGeometry now (06/2011) merged again
///  in AliEMCALGeometry
///
/// \author Magali Estienne (magali.estienne@subatech.in2p3.fr)
/// \author M.L. Wang CCNU & Subatech Adapted for DCAL Oct-18-2012
///
///
/// Usage:
///        You can create the AliEMCALGeometry object independently from anything.
///        You have to use just the correct name of geometry. If name is empty string the
///        default name of geometry will be used.
///
///  AliEMCALGeometry* geom = new AliEMCALGeometry("EMCAL_COMPLETE12SMV1","EMCAL");
///  TGeoManager::Import("geometry.root");
///
///  MC:   If you work with MC data you have to get geometry the next way:
///  ==                                      =============================
/// !!!!!!!!! This part has to be modified
///  AliRunLoader    *rl   = AliRunLoader::GetRunLoader();
///  AliEMCALEMCGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
///  TGeoManager::Import("geometry.root");
//_________________________________________________________________________

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
#include "AliEMCALTriggerMapping.h"
class AliEMCALShishKebabTrd1Module;
class AliLog;

class AliEMCALGeometry : public TNamed {

public: 
  enum fEMCSMType { kEMCAL_Standard = 0, kEMCAL_Half = 1, kEMCAL_3rd = 2, kDCAL_Standard = 3, kDCAL_Ext= 4 }; // possible SM Type

  AliEMCALGeometry();
  AliEMCALGeometry(const Text_t* name, const Text_t* title="",
                   const Text_t* mcname="", const Text_t* mctitle="");
  AliEMCALGeometry(const AliEMCALGeometry & geom);
  
  virtual ~AliEMCALGeometry(void); 
  AliEMCALGeometry & operator = (const AliEMCALGeometry  & rvalue);
  
  static AliEMCALGeometry * GetInstance();

  static AliEMCALGeometry * GetInstance(const Text_t* name,               const Text_t* title  =  "",
                                        const Text_t* mcname = "TGeant3", const Text_t* mctitle = "" ) ; 
  
  static AliEMCALGeometry * GetInstanceFromRunNumber(Int_t runNumber, 
                                                     TString geoName = "",
                                                     const Text_t* mcname  = "TGeant3", 
                                                     const Text_t* mctitle = ""        ) ;

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

  virtual Bool_t Impact(const TParticle *particle) const;
  void ImpactOnEmcal(TVector3 vtx, Double_t theta, Double_t phi, Int_t & absId, TVector3 & vimpact) const;
  Bool_t IsInEMCAL(Double_t x, Double_t y, Double_t z) const;
  Bool_t IsInDCAL(Double_t x, Double_t y, Double_t z) const;
  Int_t  IsInEMCALOrDCAL(Double_t x, Double_t y, Double_t z) const;

  //////////////////////////////////////
  // Return EMCAL geometrical parameters
  //
  
  AliEMCALEMCGeometry* GetEMCGeometry()       const { return fEMCGeometry                            ; }
  
  AliEMCALTriggerMapping* GetTriggerMapping() const { return fTriggerMapping; }
  
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
  Float_t  GetDCALInnerEdge(void)             const { return fEMCGeometry->GetDCALInnerEdge()        ; }
  Float_t  GetDCALPhiMin(void)                const { return fEMCGeometry->GetDCALPhiMin()           ; }
  Float_t  GetDCALPhiMax(void)                const { return fEMCGeometry->GetDCALPhiMax()           ; }
  Float_t  GetEMCALPhiMax(void)               const { return fEMCGeometry->GetEMCALPhiMax()          ; }
  Int_t    GetNECLayers(void)                 const { return fEMCGeometry->GetNECLayers()            ; }
  Float_t  GetDCALInnerExtandedEta(void)      const { return fEMCGeometry->GetDCALInnerExtandedEta() ; }
  Int_t    GetNZ(void)                        const { return fEMCGeometry->GetNZ()                   ; }
  Int_t    GetNEta(void)                      const { return fEMCGeometry->GetNEta()                 ; }
  Int_t    GetNPhi(void)                      const { return fEMCGeometry->GetNPhi()                 ; }
  Float_t  GetECPbRadThick(void)              const { return fEMCGeometry->GetECPbRadThick()         ; }
  Float_t  GetECScintThick(void)              const { return fEMCGeometry->GetECScintThick()         ; }
  Float_t  GetSampling(void)                  const { return fEMCGeometry->GetSampling()             ; } 
  Int_t    GetNumberOfSuperModules(void)      const { return fEMCGeometry->GetNumberOfSuperModules() ; }
  Float_t  GetPhiGapForSuperModules(void)     const { return fEMCGeometry->GetPhiGapForSuperModules(); }
  Float_t  GetPhiModuleSize(void)             const { return fEMCGeometry->GetPhiModuleSize()        ; }
  Float_t  GetEtaModuleSize(void)             const { return fEMCGeometry->GetEtaModuleSize()        ; }
  Float_t  GetFrontSteelStrip(void)           const { return fEMCGeometry->GetFrontSteelStrip()      ; }
  Float_t  GetLateralSteelStrip(void)         const { return fEMCGeometry->GetLateralSteelStrip()    ; }
  Float_t  GetPassiveScintThick(void)         const { return fEMCGeometry->GetPassiveScintThick()    ; }
  Float_t  GetPhiTileSize(void)               const { return fEMCGeometry->GetPhiTileSize()          ; }
  Float_t  GetEtaTileSize(void)               const { return fEMCGeometry->GetEtaTileSize()          ; }
  Float_t  GetPhiSuperModule(void)            const { return fEMCGeometry->GetPhiSuperModule()       ; }
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
  Int_t    GetnSupModInDCAL(void)             const { return fEMCGeometry->GetnSupModInDCAL()        ; }
  Int_t    GetILOSS(void)                     const { return fEMCGeometry->GetILOSS()                ; }
  Int_t    GetIHADR(void)                     const { return fEMCGeometry->GetIHADR()                ; }  
  // --
  Float_t  GetDeltaEta(void)                  const { return fEMCGeometry->GetDeltaEta()             ; }
  Float_t  GetDeltaPhi(void)                  const { return fEMCGeometry->GetDeltaPhi()             ; }
  Int_t    GetNTowers(void)                   const { return fEMCGeometry->GetNTowers()              ; }
  //
  Double_t GetPhiCenterOfSM(Int_t nsupmod)    const { return fEMCGeometry->GetPhiCenterOfSM(nsupmod) ; }
  Double_t GetPhiCenterOfSMSec(Int_t nsupmod) const { return fEMCGeometry->GetPhiCenterOfSMSec(nsupmod) ; }
  Float_t  GetSuperModulesPar(Int_t ipar)     const { return fEMCGeometry->GetSuperModulesPar(ipar)  ; }
  //
  Int_t    GetSMType(Int_t nSupMod)           const { if( nSupMod > fEMCGeometry->GetNumberOfSuperModules() ) return -1;
                                                      return fEMCSMSystem[nSupMod]		     ; }
  Bool_t   IsDCALSM(Int_t nSupMod) const;
  Bool_t   IsDCALExtSM(Int_t nSupMod) const;
  Bool_t   GetPhiBoundariesOfSM(Int_t nSupMod, Double_t &phiMin, Double_t &phiMax)    const 
    { return fEMCGeometry->GetPhiBoundariesOfSM(nSupMod, phiMin, phiMax)   ; }
  Bool_t   GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const 
    { return fEMCGeometry->GetPhiBoundariesOfSMGap(nPhiSec, phiMin, phiMax); }
  //
  // especially for SM in extension, where center of SM != center of the SM-section.
  // Used in AliEMCALv0 to calculate position.
  
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
    if(     GetSMType(nSupMod) == kEMCAL_Half) return fNPhi/2;
    else if(GetSMType(nSupMod) == kEMCAL_3rd)  return fNPhi/3;
    else if(GetSMType(nSupMod) == kDCAL_Ext)   return fNPhi/3;
    else                                       return fNPhi;
  } 
  // From cell indexes to abs cell id
  void    GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
					      Int_t &iphim, Int_t &ietam, Int_t &nModule) const;
  Int_t   GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const;

  void    ShiftOnlineToOfflineCellIndexes(Int_t sm, Int_t & iphi, Int_t & ieta) const ;
  void    ShiftOfflineToOnlineCellIndexes(Int_t sm, Int_t & iphi, Int_t & ieta) const ;
  
  // Methods for AliEMCALRecPoint - Feb 19, 2006
  Bool_t  RelPosCellInSModule(Int_t absId, 
                              Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t distEff,
                              Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t loc[3]) const;
  Bool_t  RelPosCellInSModule(Int_t absId, TVector3 &vloc)  const;

  Int_t  * GetEMCSystem()            const { return fEMCSMSystem          ; }     //EMC System, SM type list
  // Local Coordinates of SM
  TArrayD  GetCentersOfCellsEtaDir() const { return fCentersOfCellsEtaDir ; }     // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD  GetCentersOfCellsXDir()   const { return fCentersOfCellsXDir   ; }     // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  TArrayD  GetCentersOfCellsPhiDir() const { return fCentersOfCellsPhiDir ; }     // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  //
  TArrayD  GetEtaCentersOfCells()    const { return fEtaCentersOfCells    ; }     // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  TArrayD  GetPhiCentersOfCells()    const { return fPhiCentersOfCells    ; }     // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
	
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
  void SetMisalMatrix(const TGeoHMatrix * m, Int_t smod);
	
  //Alternate geometry that allows to calculate tower position for different particles and different alignments
  void RecalculateTowerPosition(Float_t drow, Float_t dcol, const Int_t sm, const Float_t depth,
                                const Float_t misaligTransShifts[15], const Float_t misaligRotShifts[15],Float_t global[3]) const;
  
  //Returns shift-rotational matrixes for different volumes
  const TGeoHMatrix * GetMatrixForSuperModule(Int_t smod) const ;
  const TGeoHMatrix * GetMatrixForSuperModuleFromGeoManager(Int_t smod) const ;
  const TGeoHMatrix * GetMatrixForSuperModuleFromArray(Int_t smod) const ;

  Bool_t GetAbsFastORIndexFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& id) const { 
    return fTriggerMapping->GetAbsFastORIndexFromTRU(iTRU, iADC, id);
  }
  Bool_t GetAbsFastORIndexFromPositionInTRU(const Int_t iTRU, const Int_t iEta, const Int_t iPhi, Int_t& id) const { 
    return fTriggerMapping->GetAbsFastORIndexFromPositionInTRU(iTRU, iEta, iPhi, id);
  }
  Bool_t GetAbsFastORIndexFromPositionInSM(const Int_t  iSM, const Int_t iEta, const Int_t iPhi, Int_t& id) const { 
    return fTriggerMapping->GetAbsFastORIndexFromPositionInSM( iSM, iEta, iPhi, id);
  }
  Bool_t GetAbsFastORIndexFromPositionInEMCAL(const Int_t iEta, const Int_t iPhi, Int_t& id) const { 
    return fTriggerMapping->GetAbsFastORIndexFromPositionInEMCAL(iEta, iPhi, id);
  }
  Bool_t GetAbsFastORIndexFromPHOSSubregion(const Int_t iPHOS, Int_t& id) const {
    return fTriggerMapping->GetAbsFastORIndexFromPHOSSubregion(iPHOS, id);
  }
  Bool_t GetTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iADC) const { 
    return fTriggerMapping->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
  }
  Bool_t GetPositionInTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const { 
    return fTriggerMapping->GetPositionInTRUFromAbsFastORIndex(id, iTRU, iEta, iPhi);
  }
  Bool_t GetPositionInSMFromAbsFastORIndex(const Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const { 
    return fTriggerMapping->GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi);
  }
  Bool_t GetPositionInEMCALFromAbsFastORIndex(const Int_t id, Int_t& iEta, Int_t& iPhi) const { 
    return fTriggerMapping->GetPositionInEMCALFromAbsFastORIndex(id, iEta, iPhi);
  }
  Bool_t GetFastORIndexFromCellIndex(const Int_t id, Int_t& idx) const { 
    return fTriggerMapping->GetFastORIndexFromCellIndex(id, idx);
  }
  Bool_t GetCellIndexFromFastORIndex(const Int_t id, Int_t idx[4]) const { 
    return fTriggerMapping->GetCellIndexFromFastORIndex(id, idx);
  }
  Bool_t GetTRUIndexFromSTUIndex(const Int_t id, Int_t& idx, Int_t detector) const { 
    return fTriggerMapping->GetTRUIndexFromSTUIndex(id, idx, detector);
  }
  Bool_t GetTRUIndexFromOnlineIndex(const Int_t id, Int_t& idx) const { 
    return fTriggerMapping->GetTRUIndexFromOnlineIndex(id, idx);
  }
  Bool_t GetOnlineIndexFromTRUIndex(const Int_t id, Int_t& idx) const { 
    return fTriggerMapping->GetOnlineIndexFromTRUIndex(id, idx);
  }
  Bool_t GetFastORIndexFromL0Index(const Int_t iTRU, const Int_t id, Int_t idx[], const Int_t size) const { 
    return fTriggerMapping->GetFastORIndexFromL0Index(iTRU, id, idx, size);
  }
  Int_t  GetTRUIndexFromSTUIndex(const Int_t id, Int_t detector) const { 
    return fTriggerMapping->GetTRUIndexFromSTUIndex(id, detector);
  }
  Int_t  GetTRUIndexFromOnlineIndex(const Int_t id) const { 
    return fTriggerMapping->GetTRUIndexFromOnlineIndex(id);
  }
  Int_t  GetOnlineIndexFromTRUIndex(const Int_t id) const {
    return fTriggerMapping->GetOnlineIndexFromTRUIndex(id);
  }
  Int_t  GetNTotalTRU() const { 
    return fTriggerMapping->GetNTRU(); 
  }
  Int_t  GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm)const{
    return fTriggerMapping->GetTRUIndexFromOnlineHwAdd(hwAdd, ddl, sm);
  }
  Bool_t GetSTUIndexFromTRUIndex(    const Int_t id, Int_t& idx                                    ) const { 
    return fTriggerMapping->GetSTUIndexFromTRUIndex(id, idx ); 
  }
  Int_t  GetSTUIndexFromTRUIndex(    const Int_t id                                                ) const { 
    return fTriggerMapping->GetSTUIndexFromTRUIndex(id      ); 
  }
  Bool_t GetTRUFromSTU(const Int_t iTRU, const Int_t iADC, Int_t& oTRU, Int_t& oADC, Int_t detector) const { 
    return fTriggerMapping->GetTRUFromSTU(iTRU, iADC, oTRU, oADC, detector); 
  }
  Bool_t GetSTUFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& oTRU, Int_t& oADC                ) const {
    return fTriggerMapping->GetSTUFromTRU(iTRU, iADC, oTRU, oADC          ); 
  }
  Bool_t GetTRUFromSTU(const Int_t iTRU, const Int_t ieta, const Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi, Int_t detector) const {
    return fTriggerMapping->GetTRUFromSTU(iTRU, ieta, iphi, oTRU, oeta, ophi, detector) ;
  }
  Bool_t GetSTUFromTRU(const Int_t iTRU, const Int_t ieta, const Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi                ) const {
    return fTriggerMapping->GetSTUFromTRU(iTRU, ieta, iphi, oTRU, oeta, ophi          ) ;
  }
  Int_t  GetTriggerMappingVersion() const {
    return fTriggerMapping->GetUniqueID();
  }

protected:

  void Init(void);     		     // initializes the parameters of EMCAL
  
  AliEMCALEMCGeometry * fEMCGeometry;// Geometry object for Electromagnetic calorimeter

  AliEMCALTriggerMapping* fTriggerMapping; // Trigger mapping
  
  TString  fGeoName;                 // geometry name
  Int_t    *fEMCSMSystem;	           // geometry structure
  Int_t    fKey110DEG;               // for calculation abs cell id; 19-oct-05 
  Int_t    fnSupModInDCAL;           // for calculation abs cell id; 06-nov-12
  Int_t    fNCellsInSupMod;          // number cell in super module
  Int_t    fNETAdiv;                 // number eta divizion of module
  Int_t    fNPHIdiv;                 // number phi divizion of module
  Int_t    fNCellsInModule;          // number cell in module
  TArrayD  fPhiBoundariesOfSM;       // phi boundaries of SM in rad; size is fNumberOfSuperModules;
  TArrayD  fPhiCentersOfSM;          // phi of centers of SM; size is fNumberOfSuperModules/2
  TArrayD  fPhiCentersOfSMSec;       // phi of centers of section where SM lies; size is fNumberOfSuperModules/2
  // Local Coordinates of SM
  TArrayD  fPhiCentersOfCells;       // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
  TArrayD  fCentersOfCellsEtaDir;    // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD  fCentersOfCellsPhiDir;    // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  TArrayD  fEtaCentersOfCells;       // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  Int_t    fNCells;                  // number of cells in calo
  Int_t    fNPhi;                    // Number of Towers in the PHI direction
  TArrayD  fCentersOfCellsXDir;      // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  Float_t  fEnvelop[3];              // the GEANT TUB for the detector 
  Float_t  fArm1EtaMin;              // Minimum pseudorapidity position of EMCAL in Eta
  Float_t  fArm1EtaMax;              // Maximum pseudorapidity position of EMCAL in Eta
  Float_t  fArm1PhiMin;              // Minimum angular position of EMCAL in Phi (degrees)
  Float_t  fArm1PhiMax;              // Maximum angular position of EMCAL in Phi (degrees)
  Float_t  fEtaMaxOfTRD1;            // Max eta in case of TRD1 geometry (see AliEMCALShishKebabTrd1Module)
  Float_t  fDCALPhiMin;              // Minimum angular position of DCAL in Phi (degrees)
  Float_t  fDCALPhiMax;              // Maximum angular position of DCAL in Phi (degrees)
  Float_t  fEMCALPhiMax;             // Maximum angular position of EMCAL in Phi (degrees)
  Float_t  fDCALStandardPhiMax;      // special edge for the case that DCAL contian extension
  Float_t  fDCALInnerExtandedEta;    // DCAL inner edge in Eta (with some extension)
  TList   *fShishKebabTrd1Modules;   // list of modules
  Float_t  fParSM[3];                // SM sizes as in GEANT (TRD1)
  Float_t  fPhiModuleSize;           // Phi -> X 
  Float_t  fEtaModuleSize;           // Eta -> Y 
  Float_t  fPhiTileSize;             // Size of phi tile
  Float_t  fEtaTileSize;             // Size of eta tile
  Int_t    fNZ;                      // Number of Towers in the Z direction
  Float_t  fIPDistance;		           // Radial Distance of the inner surface of the EMCAL
  Float_t  fLongModuleSize;          // Size of long module
  // Geometry Parameters
  Float_t  fShellThickness;	         // Total thickness in (x,y) direction
  Float_t  fZLength;		             // Total length in z direction
  Float_t  fSampling;		             // Sampling factor
	
  TGeoHMatrix* fkSModuleMatrix[AliEMCALGeoParams::fgkEMCALModules] ; //Orientations of EMCAL super modules
  Bool_t   fUseExternalMatrices;      // Use the matrices set in fkSModuleMatrix and not those in the geoManager
	
private:
  
  static AliEMCALGeometry *fgGeom;	               // Pointer to the unique instance of the singleton
  static Bool_t            fgInit;	               // Tells if geometry has been succesfully set up.
  static const Char_t     *fgkDefaultGeometryName; // Default name of geometry
  
  
  ClassDef(AliEMCALGeometry,17)       // EMCAL geometry class 

} ;

#endif // AliEMCALGEOUTILS_H

