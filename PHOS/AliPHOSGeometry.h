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

#include <assert.h> 

// --- ROOT system ---

#include "TNamed.h"
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h" 

// --- AliRoot header files ---

#include "AliGeometry.h"
#include "AliPHOSEMCAGeometry.h"
#include "AliPHOSCPVGeometry.h"
#include "AliPHOSPPSDGeometry.h"
#include "AliPHOSSupportGeometry.h"
#include "AliPHOSRecPoint.h"


class AliPHOSGeometry : public AliGeometry {

public: 

  AliPHOSGeometry() {
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
    fPHOSAngle = 0 ; 
  } ;  

  AliPHOSGeometry(const AliPHOSGeometry & geom) {
    // cpy ctor requested by Coding Convention but not yet needed
    assert(0==1) ;
  } 
  
  virtual ~AliPHOSGeometry(void) ; 
  static AliPHOSGeometry * GetInstance(const Text_t* name, const Text_t* title="") ; 
  static AliPHOSGeometry * GetInstance() ; 
  virtual void  GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos, TMatrix & gmat) const ;
  virtual void  GetGlobal(const AliRecPoint* RecPoint, TVector3 & gpos) const ;

  AliPHOSGeometry & operator = (const AliPHOSGeometry  & rvalue) const {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *(GetInstance()) ; 
  }
 
  // General

  static TString Degre(void) {
    // a global for degree (deg)
    return TString("deg") ; 
  }

  static TString Radian(void) { 
    // a global for radian (rad)
    return TString("rad") ; 
  }

  Bool_t AbsToRelNumbering(const Int_t AbsId, Int_t * RelId) ; // converts the absolute PHOS numbering to a relative 

  void EmcModuleCoverage(const Int_t m, Double_t & tm, Double_t & tM, Double_t & pm, Double_t & pM, Option_t * opt = Radian() );
                                                         // calculates the angular coverage in theta and phi of a EMC module
  void EmcXtalCoverage(Double_t & theta, Double_t & phi, Option_t * opt = Radian() ) ; 
                                                         // calculates the angular coverage in theta and phi of a 
                                                         // single crystal in a EMC module

  void ImpactOnEmc(const Double_t theta, const Double_t phi, Int_t & ModuleNumber, Double_t & x, Double_t & z) ; 
                                                         // calculates the impact coordinates of a neutral particle  
                                                         // emitted in direction theta and phi in ALICE
 
  void   RelPosInModule(const Int_t * RelId, Float_t & y, Float_t & z) ; // gets the position of element (pad or Xtal) relative to 
                                                                         // center of PHOS module  
  void   RelPosInAlice(const Int_t AbsId, TVector3 &  pos) ;             // gets the position of element (pad or Xtal) relative to 
                                                                         // Alice
  Bool_t RelToAbsNumbering(const Int_t * RelId, Int_t & AbsId) ;         // converts the absolute PHOS numbering to a relative 

  Bool_t  IsInitialized(void)                  const { return fgInit ; }  
                                                                       
  // Return general PHOS parameters

  Int_t   GetNModules(void)                    const { return fNModules ; }
  Int_t   GetNPPSDModules(void)                const { return fNPPSDModules ; }
  Int_t   GetNCPVModules(void)                 const { return fNModules - fNPPSDModules ; }
  Float_t GetPHOSAngle(Int_t index)            const { return fPHOSAngle[index-1] ; }

  // Return EMCA geometrical parameters

  Float_t GetOuterBoxSize(Int_t index)         const { return fGeometryEMCA->GetOuterBoxSize(index);            }
  Float_t GetAirFilledBoxSize(Int_t index)     const { return fGeometryEMCA->GetAirFilledBoxSize(index) ;       }
  Float_t GetCrystalHolderThickness(void)      const { return fGeometryEMCA->GetCrystalHolderThickness() ;      }
  Float_t GetCrystalSize(Int_t index)          const { return fGeometryEMCA->GetCrystalSize(index) ;            }
  Float_t GetCrystalSupportHeight(void)        const { return fGeometryEMCA->GetCrystalSupportHeight() ;        }
  Float_t GetCrystalWrapThickness(void)        const { return fGeometryEMCA->GetCrystalWrapThickness() ;        }
  Float_t GetGapBetweenCrystals(void)          const { return fGeometryEMCA->GetGapBetweenCrystals() ;          }
  Float_t GetIPtoCrystalSurface(void)          const { return fGeometryEMCA->GetIPtoCrystalSurface() ;          }
  Float_t GetIPtoOuterCoverDistance(void)      const { return fGeometryEMCA->GetIPtoOuterCoverDistance() ;      }
  Float_t GetLowerThermoPlateThickness(void)   const { return fGeometryEMCA->GetLowerThermoPlateThickness() ;   }
  Float_t GetLowerTextolitPlateThickness(void) const { return fGeometryEMCA->GetLowerTextolitPlateThickness() ; }
  Float_t GetModuleBoxThickness(void)          const { return fGeometryEMCA->GetModuleBoxThickness() ;          }
  Int_t   GetNPhi(void)                        const { return fGeometryEMCA->GetNPhi() ;                        }
  Int_t   GetNZ(void)                          const { return fGeometryEMCA->GetNZ() ;                          }
  Float_t GetOuterBoxThickness(Int_t index)    const { return fGeometryEMCA->GetOuterBoxThickness(index) ;      }
  Float_t GetPinDiodeSize(Int_t index)         const { return fGeometryEMCA->GetPinDiodeSize(index) ;           }
  Float_t GetSecondUpperPlateThickness(void)   const { return fGeometryEMCA->GetSecondUpperPlateThickness() ;   }
  Float_t GetSupportPlateThickness(void)       const { return fGeometryEMCA->GetSupportPlateThickness() ;       }
  Float_t GetTextolitBoxSize(Int_t index)      const { return fGeometryEMCA->GetTextolitBoxSize(index) ;        }
  Float_t GetTextolitBoxThickness(Int_t index) const { return fGeometryEMCA->GetTextolitBoxThickness(index);    }
  Float_t GetUpperPlateThickness(void)         const { return fGeometryEMCA->GetUpperPlateThickness() ;         }
  Float_t GetUpperCoolingPlateThickness(void)  const { return fGeometryEMCA->GetUpperCoolingPlateThickness() ;  }

  // Return PPSD geometrical parameters

  Float_t GetAnodeThickness(void)              const { return fGeometryPPSD->GetAnodeThickness();         }
  Float_t GetAvalancheGap(void)                const { return fGeometryPPSD->GetAvalancheGap();           }
  Float_t GetCathodeThickness(void)            const { return fGeometryPPSD->GetCathodeThickness();       }
  Float_t GetCompositeThickness(void)          const { return fGeometryPPSD->GetCompositeThickness();     }
  Float_t GetConversionGap(void)               const { return fGeometryPPSD->GetConversionGap();          }
  Float_t GetLeadConverterThickness(void)      const { return fGeometryPPSD->GetLeadConverterThickness(); }
  Float_t GetLeadToMicro2Gap(void)             const { return fGeometryPPSD->GetLeadToMicro2Gap();        }
  Float_t GetLidThickness(void)                const { return fGeometryPPSD->GetLidThickness();           }
  Float_t GetMicromegas1Thickness(void)        const { return fGeometryPPSD->GetMicromegas1Thickness();   }
  Float_t GetMicromegas2Thickness(void)        const { return fGeometryPPSD->GetMicromegas2Thickness();   }
  Float_t GetMicromegasWallThickness(void)     const { return fGeometryPPSD->GetMicromegasWallThickness();}
  Float_t GetMicro1ToLeadGap(void)             const { return fGeometryPPSD->GetMicro1ToLeadGap();        }
  Int_t   GetNumberOfModulesPhi(void)          const { return fGeometryPPSD->GetNumberOfModulesPhi();     }
  Int_t   GetNumberOfModulesZ(void)            const { return fGeometryPPSD->GetNumberOfModulesZ();       }
  Int_t   GetNumberOfPadsPhi(void)             const { return fGeometryPPSD->GetNumberOfPadsPhi();        }
  Int_t   GetNumberOfPadsZ(void)               const { return fGeometryPPSD->GetNumberOfPadsZ();          }
  Float_t GetPCThickness(void)                 const { return fGeometryPPSD->GetPCThickness();            }
  Float_t GetPhiDisplacement(void)             const { return fGeometryPPSD->GetPhiDisplacement();        }
  Float_t GetPPSDModuleSize(Int_t index)       const { return fGeometryPPSD->GetPPSDModuleSize(index);    }
  Float_t GetZDisplacement(void)               const { return fGeometryPPSD->GetZDisplacement();          }
  void    SetLeadConverterThickness(Float_t x) const {        fGeometryPPSD->SetLeadConverterThickness(x);}

  // Return CPV geometrical parameters

  Int_t   GetNumberOfCPVLayers(void)           const { return fGeometryCPV ->GetNumberOfCPVLayers();      }
  Bool_t  IsLeadConverterExists(void)          const { return fGeometryCPV->IsLeadConverterExists();      }
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

  // Common PPSD and CPV parameters

  Float_t GetCPVBoxSize(Int_t index)  const { 
    if      (strcmp(fName,"GPS2")==0) return fGeometryPPSD->GetCPVBoxSize(index);
    else if (strcmp(fName,"IHEP")==0) return fGeometryCPV ->GetCPVBoxSize(index);
    else if (strcmp(fName,"MIXT")==0) return TMath::Max(fGeometryCPV ->GetCPVBoxSize(index),
							fGeometryPPSD->GetCPVBoxSize(index));
    else                              return 0;
  }

  // Mixed EMCA and PPSD parameters

  Float_t GetIPtoPpsdUp(void)                  const {
    return (GetIPtoOuterCoverDistance() - GetCPVBoxSize(1) + GetPPSDModuleSize(1)/2 ); } 
  Float_t GetIPtoTopLidDistance(void)          const { 
    return  GetIPtoOuterCoverDistance() - GetCPVBoxSize(1) - 1. ; } 
  Float_t GetIPtoPpsdLow(void)                 const { 
    return (GetIPtoOuterCoverDistance() - GetPPSDModuleSize(1)/2 ); } 

  // Mixed EMCA and CPV parameters

  Float_t GetIPtoCPVDistance(void)             const {
    return  GetIPtoOuterCoverDistance() - GetCPVBoxSize(1) - 1.0; }

  // Return PHOS' support geometrical parameters

  Float_t GetRailOuterSize(Int_t index)  const { return fGeometrySUPP->GetRailOuterSize(index); }
  Float_t GetRailPart1    (Int_t index)  const { return fGeometrySUPP->GetRailPart1    (index); }
  Float_t GetRailPart2    (Int_t index)  const { return fGeometrySUPP->GetRailPart2    (index); }
  Float_t GetRailPart3    (Int_t index)  const { return fGeometrySUPP->GetRailPart3    (index); }
  Float_t GetRailPos      (Int_t index)  const { return fGeometrySUPP->GetRailPos      (index); }
  Float_t GetRailLength   ()             const { return fGeometrySUPP->GetRailLength   ();      }
  Float_t GetDistanceBetwRails()         const { return fGeometrySUPP->GetDistanceBetwRails();  }
  Float_t GetRailsDistanceFromIP()       const { return fGeometrySUPP->GetRailsDistanceFromIP();}
  Float_t GetRailRoadSize (Int_t index)  const { return fGeometrySUPP->GetRailRoadSize (index); }
  Float_t GetCradleWallThickness()       const { return fGeometrySUPP->GetCradleWallThickness();}
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
  Int_t                    fNPPSDModules ;   // Number of PPSD modules
  Float_t                  fAngle ;          // Position angles between modules
  Float_t                 *fPHOSAngle ;      //[fNModules] Position angles of modules
  TObjArray               *fRotMatrixArray ; // Liste of rotation matrices (one per phos module)
  AliPHOSEMCAGeometry     *fGeometryEMCA ;   // Geometry object for Electromagnetic calorimeter
  AliPHOSCPVGeometry      *fGeometryCPV ;    // Geometry object for CPV  (IHEP)
  AliPHOSPPSDGeometry     *fGeometryPPSD ;   // Geometry object for PPSD (GPS2)
  AliPHOSSupportGeometry  *fGeometrySUPP ;   // Geometry object for PHOS support

  void                     SetPHOSAngles();  // calculates the PHOS modules PHI angle

  static AliPHOSGeometry * fgGeom ; // pointer to the unique instance of the singleton 
  static Bool_t fgInit ;            // Tells if geometry has been succesfully set up 

  ClassDef(AliPHOSGeometry,1)       // PHOS geometry class 

} ;

#endif // AliPHOSGEOMETRY_H
