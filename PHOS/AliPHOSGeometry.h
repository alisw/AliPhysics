#ifndef ALIPHOSGEOMETRY_H
#define ALIPHOSGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : singleton  
// The EMC modules are parametrized so that any configuration can be easily implemented 
// The title is used to identify the type of CPV used. PPSD and CPV are implemented
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

  Float_t GetAnodeThickness(void)              const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetAnodeThickness();         }
  Float_t GetAvalancheGap(void)                const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetAvalancheGap();           }
  Float_t GetCathodeThickness(void)            const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetCathodeThickness();       }
  Float_t GetCompositeThickness(void)          const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetCompositeThickness();     }
  Float_t GetConversionGap(void)               const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetConversionGap();          }
  Float_t GetLeadConverterThickness(void)      const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetLeadConverterThickness(); }
  Float_t GetLeadToMicro2Gap(void)             const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetLeadToMicro2Gap();        }
  Float_t GetLidThickness(void)                const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetLidThickness();           }
  Float_t GetMicromegas1Thickness(void)        const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetMicromegas1Thickness();   }
  Float_t GetMicromegas2Thickness(void)        const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetMicromegas2Thickness();   }
  Float_t GetMicromegasWallThickness(void)     const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetMicromegasWallThickness();}
  Float_t GetMicro1ToLeadGap(void)             const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetMicro1ToLeadGap();        }
  Float_t GetPCThickness(void)                 const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetPCThickness();            }
  Float_t GetPhiDisplacement(void)             const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetPhiDisplacement();        }
  Float_t GetPPSDModuleSize(Int_t index)       const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetPPSDModuleSize(index);    }
  Float_t GetZDisplacement(void)               const { return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetZDisplacement();          }

  // Return CPV geometrical parameters

  Bool_t  IsLeadConverterExists(void)          const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->IsLeadConverterExists();      }
  Float_t GetCPVActiveSize(Int_t index)        const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVActiveSize(index);         }
  Int_t   GetNumberOfCPVChipsPhi(void)         const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetNumberOfCPVChipsPhi();     }
  Int_t   GetNumberOfCPVChipsZ(void)           const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetNumberOfCPVChipsZ();       }
  Float_t GetGassiplexChipSize(Int_t index)    const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetGassiplexChipSize(index);  }
  Float_t GetCPVGasThickness(void)             const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVGasThickness();         }
  Float_t GetCPVTextoliteThickness(void)       const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVTextoliteThickness();   }
  Float_t GetCPVCuNiFoilThickness(void)        const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVCuNiFoilThickness();    }
  Float_t GetFTPosition(Int_t index)           const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetFTPosition(index);         }
  Float_t GetCPVFrameSize(Int_t index)         const { return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVFrameSize(index);       }

  // Common PPSD and CPV parameters

  Int_t   GetNumberOfCPVLayers(void)  const {
    if      (strcmp(fName,"GPS2")==0) return 2;
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetNumberOfCPVLayers();
    else                              return 0;
  }

  Float_t GetCPVBoxSize(Int_t index)  const { 
    if      (strcmp(fName,"GPS2")==0) return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetCPVBoxSize(index);
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry* ) fGeometryCPV)->GetCPVBoxSize(index);
    else                              return 0;
  }

  Int_t   GetNumberOfModulesPhi(void) const {
    if      (strcmp(fName,"GPS2")==0) return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetNumberOfModulesPhi();
    else if (strcmp(fName,"IHEP")==0) return 1;
    else                              return 0;
  }

  Int_t   GetNumberOfModulesZ(void)   const {
    if      (strcmp(fName,"GPS2")==0) return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetNumberOfModulesZ();
    else if (strcmp(fName,"IHEP")==0) return 1;
    else                              return 0;
  }

  Int_t   GetNumberOfPadsPhi(void)    const { 
    if      (strcmp(fName,"GPS2")==0) return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetNumberOfPadsPhi();
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry* ) fGeometryCPV)->GetNumberOfCPVPadsPhi();
    else                              return 0;
  }

  Int_t   GetNumberOfPadsZ(void)      const { 
    if      (strcmp(fName,"GPS2")==0) return ((AliPHOSPPSDGeometry*) fGeometryCPV)->GetNumberOfPadsZ();
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry* ) fGeometryCPV)->GetNumberOfCPVPadsZ();
    else                              return 0;
  }

  Float_t GetPadSizePhi(void)         const {
    if      (strcmp(fName,"GPS2")==0) return GetPPSDModuleSize(0) / GetNumberOfPadsPhi();
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVPadSizePhi();
    else                              return 0;
  }

  Float_t GetPadSizeZ(void)           const {
    if      (strcmp(fName,"GPS2")==0) return GetPPSDModuleSize(2) / GetNumberOfPadsZ();
    else if (strcmp(fName,"IHEP")==0) return ((AliPHOSCPVGeometry*) fGeometryCPV)->GetCPVPadSizeZ();
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


protected:

  AliPHOSGeometry(const Text_t* name, const Text_t* title="") : AliGeometry(name, title) { 
    // ctor only for internal usage (singleton)
    Init() ; 
  }
  void Init(void) ;            // steering method for PHOS and PPSD/CPV

private:

  Int_t                    fNModules ;       // Number of modules constituing PHOS
  Float_t                 *fPHOSAngle ;      //[fNModules] Position angles of modules
  TObjArray               *fRotMatrixArray ; // Liste of rotation matrices (one per phos module)
  AliPHOSEMCAGeometry     *fGeometryEMCA ;   // Geometry object for Electromagnetic calorimeter
  AliPHOSCPVBaseGeometry  *fGeometryCPV ;    // Geometry object for CPV (either GPS2 or IHEP)

  void                 SetPHOSAngles(); // calculates the PHOS modules PHI angle

  static AliPHOSGeometry * fgGeom ; // pointer to the unique instance of the singleton 
  static Bool_t fgInit ;            // Tells if geometry has been succesfully set up 

  ClassDef(AliPHOSGeometry,1)       // PHOS geometry class 

} ;

#endif // AliPHOSGEOMETRY_H
