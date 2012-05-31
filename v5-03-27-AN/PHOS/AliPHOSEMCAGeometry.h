#ifndef ALIPHOSEMCAGEOMETRY_H
#define ALIPHOSEMCAGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : EMCA (Electromagnetic Calirometer)
// Its data members provide geometry parametrization of EMCA
// which can be changed in the constructor only.
// Author:   Yves Schutz (Subatech)
// Modified: Yuri Kharlov (IHEP, Protvino)
// 15 September 2000

// --- ROOT system ---
#include "TObject.h"
class TObjArray ;

// --- AliRoot header files ---

class AliPHOSEMCAGeometry : public TObject {

public: 

  AliPHOSEMCAGeometry();
  //Compiler-generated copy ctor and copy-assignment operator _ARE_ OK.

  virtual ~AliPHOSEMCAGeometry(void) {}

  const Float_t * GetStripHalfSize()          const {return fStripHalfSize ;         }
  Float_t   GetStripWallWidthOut()      const {return fStripWallWidthOut ;     }
  const Float_t * GetAirCellHalfSize()        const {return fAirCellHalfSize ;       }
  const Float_t * GetWrappedHalfSize()        const {return fWrappedHalfSize ;       }
  Float_t   GetAirGapLed()              const {return fAirGapLed ;             }
  const Float_t * GetCrystalHalfSize()        const {return fCrystalHalfSize ;       }
  const Float_t * GetSupportPlateHalfSize()   const {return fSupportPlateHalfSize ;  }
  const Float_t * GetSupportPlateInHalfSize() const {return fSupportPlateInHalfSize ;}
  Float_t   GetSupportPlateThickness()  const { return fSupportPlateThickness ;}    

  const Float_t * GetPreampHalfSize()         const {return fPreampHalfSize ;        }
  const Float_t * GetAPDHalfSize()            const {return fPinDiodeHalfSize ;      }
  const Float_t * GetOuterThermoParams()      const {return  fOuterThermoParams ;    }
  const Float_t * GetCoolerHalfSize()         const {return fCoolerHalfSize ;        }
  const Float_t * GetAirGapHalfSize()         const {return fAirGapHalfSize;         }
  const Float_t * GetInnerThermoHalfSize()    const {return  fInnerThermoHalfSize ;  }
  const Float_t * GetAlCoverParams()          const {return fAlCoverParams ;         }
  const Float_t * GetFiberGlassHalfSize()     const {return fFiberGlassHalfSize ;    }
  const Float_t * GetWarmAlCoverHalfSize()    const {return fWarmAlCoverHalfSize ;   }
  const Float_t * GetWarmThermoHalfSize()     const {return fWarmThermoHalfSize ;    }
  const Float_t * GetTSupport1HalfSize()      const {return fTSupport1HalfSize ;     }
  const Float_t * GetTSupport2HalfSize()      const {return fTSupport2HalfSize ;     }
  const Float_t * GetTCables1HalfSize()       const {return fTCables1HalfSize ;      }
  const Float_t * GetTCables2HalfSize()       const {return fTCables2HalfSize ;      }
  Float_t   GetTSupportDist()           const {return fTSupportDist ;          }
  const Float_t * GetFrameXHalfSize()         const {return fFrameXHalfSize ;        }
  const Float_t * GetFrameZHalfSize()         const {return fFrameZHalfSize ;        }
  const Float_t * GetFrameXPosition()         const {return fFrameXPosition ;        }
  const Float_t * GetFrameZPosition()         const {return fFrameZPosition ;        }
  const Float_t * GetFGupXHalfSize()          const {return fFGupXHalfSize ;         }
  const Float_t * GetFGupXPosition()          const {return fFGupXPosition ;         }
  const Float_t * GetFGupZHalfSize()          const {return fFGupZHalfSize ;         }
  const Float_t * GetFGupZPosition()          const {return fFGupZPosition ;         }
  const Float_t * GetFGlowXHalfSize()         const {return fFGlowXHalfSize ;        }
  const Float_t * GetFGlowXPosition()         const {return fFGlowXPosition ;        }
  const Float_t * GetFGlowZHalfSize()         const {return fFGlowZHalfSize ;        }
  const Float_t * GetFGlowZPosition()         const {return fFGlowZPosition ;        }
  const Float_t * GetFEEAirHalfSize()         const {return fFEEAirHalfSize ;        }
  const Float_t * GetFEEAirPosition()         const {return fFEEAirPosition ;        }
  const Float_t * GetEMCParams()              const {return fEMCParams ;             }

  Float_t    GetIPtoCrystalSurface()     const {return fIPtoCrystalSurface ;   }
  Float_t    GetIPtoOuterCoverDistance() const {return fIPtoOuterCoverDistance;}
  Float_t    GetCrystalSize(Int_t index) const {return 2.*fCrystalHalfSize[index] ; }

  Int_t     GetNCellsXInStrip() const { return fNCellsXInStrip;}
  Int_t     GetNCellsZInStrip() const { return fNCellsZInStrip;}
  Int_t     GetNStripX()        const { return fNStripX ;      }
  Int_t     GetNStripZ()        const { return fNStripZ ;      }
  Int_t     GetNTSuppots()      const { return fNTSupports;    }
  Int_t     GetNPhi()           const { return fNPhi ;         }
  Int_t     GetNZ()             const { return fNZ ;           }
 
private:

  Float_t fStripHalfSize[3]   ;        // Strip unit size/2
  Float_t fAirCellHalfSize[3] ;        // geometry parameter
  Float_t fWrappedHalfSize[3] ;        // geometry parameter
  Float_t fSupportPlateHalfSize[3] ;   // geometry parameter
  Float_t fSupportPlateInHalfSize[3] ; // geometry parameter
  Float_t fCrystalHalfSize[3] ;        // crystal size/2
  Float_t fAirGapLed ;                 // geometry parameter
  Float_t fStripWallWidthOut ;         // Side to another strip  
  Float_t fStripWallWidthIn ;          // geometry parameter
  Float_t fTyvecThickness ;            // geometry parameter
  Float_t fTSupport1HalfSize[3] ;      // geometry parameter
  Float_t fTSupport2HalfSize[3] ;      // geometry parameter
  Float_t fPreampHalfSize[3] ;         // geometry parameter
  Float_t fPinDiodeHalfSize[3] ;       // Size of the PIN Diode 

  Float_t fOuterThermoParams[4] ;      // geometry parameter
  Float_t fCoolerHalfSize[3] ;         // geometry parameter
  Float_t fAirGapHalfSize[3] ;         // geometry parameter
  Float_t fInnerThermoHalfSize[3] ;    // geometry parameter
  Float_t fAlCoverParams[4] ;          // geometry parameter
  Float_t fFiberGlassHalfSize[3] ;     // geometry parameter

  Float_t fInnerThermoWidthX ;         // geometry parameter
  Float_t fInnerThermoWidthY ;         // geometry parameter
  Float_t fInnerThermoWidthZ ;         // geometry parameter
  Float_t fAirGapWidthX ;              // geometry parameter
  Float_t fAirGapWidthY ;              // geometry parameter
  Float_t fAirGapWidthZ ;              // geometry parameter
  Float_t fCoolerWidthX ;              // geometry parameter
  Float_t fCoolerWidthY ;              // geometry parameter
  Float_t fCoolerWidthZ ;              // geometry parameter
  Float_t fAlCoverThickness ;          // geometry parameter
  Float_t fOuterThermoWidthXUp ;       // geometry parameter
  Float_t fOuterThermoWidthXLow;       // geometry parameter
  Float_t fOuterThermoWidthY ;         // geometry parameter
  Float_t fOuterThermoWidthZ ;         // geometry parameter
  Float_t fAlFrontCoverX  ;            // geometry parameter
  Float_t fAlFrontCoverZ  ;            // geometry parameter
  Float_t fFiberGlassSup2X ;           // geometry parameter
  Float_t fFiberGlassSup1X ;           // geometry parameter
  Float_t fFrameHeight ;               // geometry parameter
  Float_t fFrameThickness ;            // geometry parameter
  Float_t fAirSpaceFeeX ;              // geometry parameter
  Float_t fAirSpaceFeeZ ;              // geometry parameter
  Float_t fAirSpaceFeeY ;              // geometry parameter
  Float_t fTCables2HalfSize[3] ;       // geometry parameter
  Float_t fTCables1HalfSize[3] ;       // geometry parameter
  Float_t fWarmUpperThickness ;        // geometry parameter
  Float_t fWarmBottomThickness ;       // geometry parameter
  Float_t fWarmAlCoverWidthX ;         // geometry parameter
  Float_t fWarmAlCoverWidthY ;         // geometry parameter
  Float_t fWarmAlCoverWidthZ ;         // geometry parameter
  Float_t fWarmAlCoverHalfSize[3] ;    // geometry parameter
  Float_t fWarmThermoHalfSize[3] ;     // geometry parameter
  Float_t fFiberGlassSup1Y ;           // geometry parameter
  Float_t fFiberGlassSup2Y ;           // geometry parameter
  Float_t fTSupportDist ;              // geometry parameter
  Float_t fTSupport1Thickness ;        // geometry parameter
  Float_t fTSupport2Thickness ;        // geometry parameter
  Float_t fTSupport1Width ;            // geometry parameter
  Float_t fTSupport2Width ;            // geometry parameter
  Float_t fFrameXHalfSize[3] ;         // geometry parameter
  Float_t fFrameZHalfSize[3] ;         // geometry parameter
  Float_t fFrameXPosition[3] ;         // geometry parameter
  Float_t fFrameZPosition[3] ;         // geometry parameter
  Float_t fFGupXHalfSize[3] ;          // geometry parameter
  Float_t fFGupXPosition[3] ;          // geometry parameter
  Float_t fFGupZHalfSize[3] ;          // geometry parameter
  Float_t fFGupZPosition[3] ;          // geometry parameter
  Float_t fFGlowXHalfSize[3] ;         // geometry parameter
  Float_t fFGlowXPosition[3] ;         // geometry parameter
  Float_t fFGlowZHalfSize[3] ;         // geometry parameter
  Float_t fFGlowZPosition[3] ;         // geometry parameter
  Float_t fFEEAirHalfSize[3] ;         // geometry parameter
  Float_t fFEEAirPosition[3] ;         // geometry parameter
  Float_t fEMCParams[4] ;              // geometry parameter
  Float_t fIPtoOuterCoverDistance ;    // Distances from interaction point to outer cover 
  Float_t fIPtoCrystalSurface ;        // Distances from interaction point to Xtal surface

  Float_t fSupportPlateThickness ;     // Thickness of the Aluminium support plate for Strip   

  Int_t  fNCellsXInStrip ;             // Number of cells in a strip unit in X
  Int_t  fNCellsZInStrip ;             // Number of cells in a strip unit in Z
  Int_t  fNStripX ;                    // Number of strip units in X
  Int_t  fNStripZ ;                    // Number of strip units in Z
  Int_t  fNTSupports ;                 // geometry parameter
  Int_t  fNPhi ;                       // Number of crystal units in X (phi) direction
  Int_t  fNZ ;                         // Number of crystal units in Z direction
  ClassDef(AliPHOSEMCAGeometry,1)      // EMCA geometry class 

} ;

#endif // AliPHOSEMCAGEOMETRY_H
