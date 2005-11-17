#ifndef ALIPHOSGEOMETRY_H
#define ALIPHOSGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : singleton
// PHOS consists of the electromagnetic calorimeter (EMCA)
// and a charged particle veto either in the Subatech's version (PPSD)
// or in the IHEP's one (CPV).
// The EMCA/PPSD/CPV modules are parametrized so that any configuration
// can be easily implemented 
// The title is used to identify the version of CPV used.
// 
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- AliRoot header files ---

#include "AliGeometry.h"
#include "AliPHOSEMCAGeometry.h"
#include "AliPHOSCPVGeometry.h"
#include "AliPHOSSupportGeometry.h"


class AliPHOSGeometry : public AliGeometry {

public: 

  AliPHOSGeometry() ;

  AliPHOSGeometry(const AliPHOSGeometry & geom) : AliGeometry(geom) {
    Fatal("cpy ctor", "not implemented") ; 
  } 
  
  virtual ~AliPHOSGeometry(void) ; 
  static AliPHOSGeometry * GetInstance(const Text_t* name, const Text_t* title="") ; 
  static AliPHOSGeometry * GetInstance() ; 
  virtual void   GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos, TMatrix & gmat) const ;
  virtual void   GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos) const ;
  virtual Bool_t Impact(const TParticle * particle) const ;

  AliPHOSGeometry & operator = (const AliPHOSGeometry  & /*rvalue*/) const {
    Fatal("operator =", "nt implemented") ; return *(GetInstance()) ; }
 
  // General

  static TString Degre(void) { return TString("deg") ; }  // a global for degree (deg)

  static TString Radian(void){ return TString("rad") ; }  // a global for radian (rad)

  Bool_t AbsToRelNumbering(Int_t AbsId, Int_t * RelId) const ; 
                                          // converts the absolute PHOS numbering to a relative 

  void EmcModuleCoverage(Int_t m, Double_t & tm, Double_t & tM, Double_t & pm, 
			                  Double_t & pM, Option_t * opt = Radian() ) const ;
                                         // calculates the angular coverage in theta and phi of a EMC module
  void EmcXtalCoverage(Double_t & theta, Double_t & phi, Option_t * opt = Radian() ) const ; 
                                         // calculates the angular coverage in theta and phi of a  
                                         // single crystal in a EMC module
  void ImpactOnEmc(Double_t theta, Double_t phi, Int_t & ModuleNumber, 
		         Double_t & z, Double_t & x) const ; 
  void ImpactOnEmc(const TVector3& vec, Int_t & ModuleNumber, 
		         Double_t & z, Double_t & x) const ; 
  void ImpactOnEmc(const TParticle& p, Int_t & ModuleNumber, 
		         Double_t & z, Double_t & x) const ; 
                                        // calculates the impact coordinates of a neutral particle  
                                         // emitted in direction theta and phi in ALICE
  Bool_t IsInEMC(Int_t id) const { if (id > GetNModules() *  GetNCristalsInModule() ) return kFALSE; return kTRUE; } 
  void RelPosInModule(const Int_t * RelId, Float_t & y, Float_t & z) const ; 
                                         // gets the position of element (pad or Xtal) relative to 
                                         // center of PHOS module  
  void RelPosInAlice(Int_t AbsId, TVector3 &  pos) const ;             
                                         // gets the position of element (pad or Xtal) relative to Alice
  Bool_t RelToAbsNumbering(const Int_t * RelId, Int_t & AbsId) const ;         
                                         // converts the absolute PHOS numbering to a relative 
  void  RelPosToAbsId(Int_t module, Double_t x, Double_t z, Int_t & AbsId) const; 
                                         // converts local PHOS-module (x, z) coordinates to absId 

  Bool_t IsInitialized(void)                  const { return fgInit ; }  
                                                                       
  // Return general PHOS parameters
  Int_t    GetNModules(void)                    const { return fNModules ; }
  Float_t  GetPHOSAngle(Int_t index)            const { return fPHOSAngle[index-1] ; }
  Float_t* GetPHOSParams(void)                        { return fPHOSParams;}  //Half-sizes of PHOS trapecoid
  Float_t  GetIPtoUpperCPVsurface(void)         const { return fIPtoUpperCPVsurface ; }
  Float_t  GetOuterBoxSize(Int_t index)         const { return 2.*fPHOSParams[index]; }
  Float_t  GetCrystalSize(Int_t index)          const { return fGeometryEMCA->GetCrystalSize(index) ;  }
  Float_t  GetCellStep(void)                    const { return 2*(fGeometryEMCA->GetAirCellHalfSize()[0] + 
							          fGeometryEMCA->GetStripWallWidthOut()) ;}

  // Return EMCA geometry parameters

  AliPHOSEMCAGeometry * GetEMCAGeometry()      const {return fGeometryEMCA ;}
  Float_t   GetIPtoCrystalSurface(void)        const { return fGeometryEMCA->GetIPtoCrystalSurface() ; }
  Float_t   GetIPtoOuterCoverDistance(void)    const { return fGeometryEMCA->GetIPtoOuterCoverDistance() ; }
  Int_t     GetNPhi(void)                      const { return fGeometryEMCA->GetNPhi() ; }
  Int_t     GetNZ(void)                        const { return fGeometryEMCA->GetNZ() ; }
  Int_t     GetNCristalsInModule(void)         const { return fGeometryEMCA->GetNPhi() * fGeometryEMCA->GetNZ() ; }

  // Return CPV geometry parameters
  Int_t   GetNumberOfCPVLayers(void)           const { return fGeometryCPV ->GetNumberOfCPVLayers();      }
  Float_t GetCPVActiveSize(Int_t index)        const { return fGeometryCPV->GetCPVActiveSize(index);      }
  Int_t   GetNumberOfCPVChipsPhi(void)         const { return fGeometryCPV->GetNumberOfCPVChipsPhi();     }
  Int_t   GetNumberOfCPVChipsZ(void)           const { return fGeometryCPV->GetNumberOfCPVChipsZ();       }
  Int_t   GetNumberOfCPVPadsPhi(void)          const { return fGeometryCPV->GetNumberOfCPVPadsPhi();      }
  Int_t   GetNumberOfCPVPadsZ(void)            const { return fGeometryCPV->GetNumberOfCPVPadsZ();        }
  Float_t GetPadSizePhi(void)                  const { return fGeometryCPV->GetCPVPadSizePhi();           }
  Float_t GetPadSizeZ(void)                    const { return fGeometryCPV->GetCPVPadSizeZ();             }
  Float_t GetGassiplexChipSize(Int_t index)    const { return fGeometryCPV->GetGassiplexChipSize(index);  }
  Float_t GetCPVGasThickness(void)             const { return fGeometryCPV->GetCPVGasThickness();         }
  Float_t GetCPVTextoliteThickness(void)       const { return fGeometryCPV->GetCPVTextoliteThickness();   }
  Float_t GetCPVCuNiFoilThickness(void)        const { return fGeometryCPV->GetCPVCuNiFoilThickness();    }
  Float_t GetFTPosition(Int_t index)           const { return fGeometryCPV->GetFTPosition(index);         }
  Float_t GetCPVFrameSize(Int_t index)         const { return fGeometryCPV->GetCPVFrameSize(index);       }
  Float_t GetCPVBoxSize(Int_t index)           const { return fGeometryCPV ->GetCPVBoxSize(index);        } 
  Float_t GetIPtoCPVDistance(void)             const { return  GetIPtoOuterCoverDistance() - 
							       GetCPVBoxSize(1) - 1.0; }
  void GetModuleCenter(TVector3& center, const char *det, Int_t module) const;
  void Global2Local(TVector3& localPosition,
		    const TVector3& globalPosition,
		    Int_t module) const;

  // Return PHOS' support geometry parameters

  Float_t GetRailOuterSize(Int_t index)  const { return fGeometrySUPP->GetRailOuterSize(index); }
  Float_t GetRailPart1    (Int_t index)  const { return fGeometrySUPP->GetRailPart1    (index); }
  Float_t GetRailPart2    (Int_t index)  const { return fGeometrySUPP->GetRailPart2    (index); }
  Float_t GetRailPart3    (Int_t index)  const { return fGeometrySUPP->GetRailPart3    (index); }
  Float_t GetRailPos      (Int_t index)  const { return fGeometrySUPP->GetRailPos      (index); }
  Float_t GetRailLength   (void)         const { return fGeometrySUPP->GetRailLength   ();      }
  Float_t GetDistanceBetwRails(void)     const { return fGeometrySUPP->GetDistanceBetwRails();  }
  Float_t GetRailsDistanceFromIP(void)   const { return fGeometrySUPP->GetRailsDistanceFromIP();}
  Float_t GetRailRoadSize (Int_t index)  const { return fGeometrySUPP->GetRailRoadSize (index); }
  Float_t GetCradleWallThickness(void)   const { return fGeometrySUPP->GetCradleWallThickness();}
  Float_t GetCradleWall   (Int_t index)  const { return fGeometrySUPP->GetCradleWall   (index); }
  Float_t GetCradleWheel  (Int_t index)  const { return fGeometrySUPP->GetCradleWheel  (index); }

protected:

  AliPHOSGeometry(const Text_t* name, const Text_t* title="") : AliGeometry(name, title) { 
    // ctor only for internal usage (singleton)
    Init() ; 
  }
  void Init(void) ;            // steering method for PHOS and PPSD/CPV

private:

  Int_t                    fNModules ;       // Number of modules constituing PHOS
  Float_t                  fAngle ;          // Position angles between modules
  Float_t                 *fPHOSAngle ;      //[fNModules] Position angles of modules
  Float_t                  fPHOSParams[4] ;  // Half-sizes of PHOS trapecoid
  Float_t                  fIPtoUpperCPVsurface; // Minimal distance from IP to PHOS
  TObjArray               *fRotMatrixArray ; // Liste of rotation matrices (one per phos module)
  AliPHOSEMCAGeometry     *fGeometryEMCA ;   // Geometry object for Electromagnetic calorimeter
  AliPHOSCPVGeometry      *fGeometryCPV ;    // Geometry object for CPV  (IHEP)
  AliPHOSSupportGeometry  *fGeometrySUPP ;   // Geometry object for PHOS support

  void                     SetPHOSAngles();  // calculates the PHOS modules PHI angle

  static AliPHOSGeometry * fgGeom ; // pointer to the unique instance of the singleton 
  static Bool_t fgInit ;            // Tells if geometry has been succesfully set up 

  ClassDef(AliPHOSGeometry,1)       // PHOS geometry class 

} ;

#endif // AliPHOSGEOMETRY_H
