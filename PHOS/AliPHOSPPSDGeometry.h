#ifndef ALIPHOSPPSDGEOMETRY_H
#define ALIPHOSPPSDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry derived class for PHOS:PPSD (PHOS Preshower Detector)
//
//*-- Author  : Yves Schutz
//    Modified: Yuri Kharlov (IHEP, Protvino)
//    14 September 2000

#include <assert.h> 

#include "AliPHOSCPVBaseGeometry.h"

class AliPHOSPPSDGeometry : public AliPHOSCPVBaseGeometry {

public: 

           AliPHOSPPSDGeometry();
  virtual ~AliPHOSPPSDGeometry(void) {};

  // PPSD functions

  Float_t GetAnodeThickness(void)          const { return  fAnodeThickness ;          }
  Float_t GetAvalancheGap(void)            const { return  fAvalancheGap ;            }
  Float_t GetCathodeThickness(void)        const { return  fCathodeThickness ;        }
  Float_t GetCompositeThickness(void)      const { return  fCompositeThickness ;      }
  Float_t GetConversionGap(void)           const { return  fConversionGap ;           }
  Float_t GetLeadConverterThickness(void)  const { return  fLeadConverterThickness ;  }
  Float_t GetLeadToMicro2Gap(void)         const { return  fLeadToMicro2Gap ;         }
  Float_t GetLidThickness(void)            const { return  fLidThickness ;            }
  Float_t GetMicromegas1Thickness(void)    const { return  fMicromegas1Thickness ;    }
  Float_t GetMicromegas2Thickness(void)    const { return  fMicromegas2Thickness ;    }
  Float_t GetMicromegasWallThickness(void) const { return  fMicromegasWallThickness ; }
  Float_t GetMicro1ToLeadGap(void)         const { return  fMicro1ToLeadGap ;         }
  Int_t   GetNumberOfPadsPhi(void)         const { return  fNumberOfPadsPhi ;         }
  Int_t   GetNumberOfPadsZ(void)           const { return  fNumberOfPadsZ ;           }
  Int_t   GetNumberOfModulesPhi(void)      const { return  fNumberOfModulesPhi ;      }
  Int_t   GetNumberOfModulesZ(void)        const { return  fNumberOfModulesZ ;        }
  Float_t GetPCThickness(void)             const { return  fPCThickness ;             }
  Float_t GetPhiDisplacement(void)         const { return  fPhiDisplacement ;         }
  Float_t GetCPVBoxSize(Int_t index)       const { return  fPPSDBoxSize[index] ;      }
  Float_t GetPPSDModuleSize(Int_t index)   const { return  fPPSDModuleSize[index] ;   }
  Float_t GetZDisplacement(void)           const { return  fZDisplacement ;           }
 
  // CPV functions cannot be used for PPSD

  Int_t   GetNumberOfCPVLayers(void)       { AssertCPV(); return 0; }
  Bool_t  IsLeadConverterExists(void)      { AssertCPV(); return 0; }
  Float_t GetCPVActiveSize(Int_t index)    { AssertCPV(); return 0; }
  Int_t   GetNumberOfCPVChipsPhi(void)     { AssertCPV(); return 0; }
  Int_t   GetNumberOfCPVChipsZ(void)       { AssertCPV(); return 0; }
  Float_t GetGassiplexChipSize(Int_t index){ AssertCPV(); return 0; }
  Float_t GetCPVGasThickness(void)         { AssertCPV(); return 0; }
  Float_t GetCPVTextoliteThickness(void)   { AssertCPV(); return 0; }
  Float_t GetCPVCuNiFoilThickness(void)    { AssertCPV(); return 0; }
  Float_t GetFTPosition(Int_t index)       { AssertCPV(); return 0; }
  Float_t GetCPVFrameSize(Int_t index)     { AssertCPV(); return 0; }
  Float_t GetIPtoCPVDistance(void)         { AssertCPV(); return 0; }

private:

  Float_t fAnodeThickness ;          // Thickness of the copper layer which makes the anode 
  Float_t fAvalancheGap ;            // Thickness of the gas in the avalanche stage
  Float_t fCathodeThickness ;        // Thickeness of composite material ensuring rigidity of cathode
  Float_t fCompositeThickness ;      // Thickeness of composite material ensuring rigidity of anode
  Float_t fConversionGap ;           // Thickness of the gas in the conversion stage
  Float_t fLeadConverterThickness ;  // Thickness of the Lead converter 
  Float_t fLeadToMicro2Gap ;         // Thickness of the air gap between the Lead and Micromegas 2        
  Float_t fLidThickness ;            // Thickness of top lid 
  Float_t fMicromegas1Thickness ;    // Thickness of the first downstream Micromegas 
  Float_t fMicromegas2Thickness ;    // Thickness of the second downstream Micromegas 
  Float_t fMicromegasWallThickness ; // Thickness of the Micromegas leak tight box
  Float_t fMicro1ToLeadGap ;         // Thickness of the air gap between Micromegas 1 and the Lead
  Int_t   fNumberOfPadsPhi ;         // Number of pads on a micromegas module ;  
  Int_t   fNumberOfPadsZ ;           // Number of pads on a micromegas module ;  
  Int_t   fNumberOfModulesPhi ;      // Number of micromegas modules in phi
  Int_t   fNumberOfModulesZ ;        // Number of micromegas modules in z
  Float_t fPCThickness ;             // Thickness of the printed circuit board of the anode   
  Float_t fPhiDisplacement ;         // Phi displacement of micromegas1 with respect to micromegas2  
  Float_t fPPSDBoxSize[3] ;          // Size of large box which contains PPSD; matches PHOS module size
  Float_t fPPSDModuleSize[3] ;       // Size of an individual micromegas module
  Float_t fZDisplacement ;           // Z displacement of micromegas1 with respect to micromegas2  

  Float_t fIPtoTopLidDistance ;      // Distance from interaction point to top lid of PPSD

  void    AssertCPV() {
    printf("Function %s should not be called for PPSD geometry\n",__PRETTY_FUNCTION__);
    assert(0==1) ;
  }

  ClassDef(AliPHOSPPSDGeometry,1)        // PPSD geometry class 

} ;

#endif // AliPHOSPPSDGEOMETRY_H
