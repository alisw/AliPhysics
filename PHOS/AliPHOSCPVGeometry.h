#ifndef ALIPHOSCPVGEOMETRY_H
#define ALIPHOSCPVGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry derived class for PHOS:CPV (Charged particle veto, IHEP version)
//
//*-- Author:  Yuri Kharlov (IHEP, Protvino)
//    14 September 2000

#include <assert.h> 

#include "AliPHOSCPVBaseGeometry.h"

class AliPHOSCPVGeometry : public AliPHOSCPVBaseGeometry {

public: 

           AliPHOSCPVGeometry();
  virtual ~AliPHOSCPVGeometry(void) {};

  // CPV functions

  Int_t   GetNumberOfCPVLayers(void)        const { return  fNumberOfCPVLayers;        }
  Bool_t  IsLeadConverterExists(void)       const { return  fLeadConverterExists;      }
  Int_t   GetNumberOfCPVPadsPhi(void)       const { return  fNumberOfCPVPadsPhi ;      }
  Int_t   GetNumberOfCPVPadsZ(void)         const { return  fNumberOfCPVPadsZ ;        }
  Float_t GetCPVPadSizePhi(void)            const { return  fCPVPadSizePhi;            }
  Float_t GetCPVPadSizeZ(void)              const { return  fCPVPadSizeZ;              }
  Float_t GetCPVBoxSize(Int_t index)        const { return  fCPVBoxSize[index];        }
  Float_t GetCPVActiveSize(Int_t index)     const { return  fCPVActiveSize[index];     }
  Int_t   GetNumberOfCPVChipsPhi(void)      const { return  fNumberOfCPVChipsPhi;      }
  Int_t   GetNumberOfCPVChipsZ(void)        const { return  fNumberOfCPVChipsZ;        }
  Float_t GetGassiplexChipSize(Int_t index) const { return  fGassiplexChipSize[index]; }
  Float_t GetCPVGasThickness(void)          const { return  fCPVGasThickness;          }
  Float_t GetCPVTextoliteThickness(void)    const { return  fCPVTextoliteThickness;    }
  Float_t GetCPVCuNiFoilThickness(void)     const { return  fCPVCuNiFoilThickness;     }
  Float_t GetFTPosition(Int_t index)        const { return  fFTPosition[index];        }
  Float_t GetCPVFrameSize(Int_t index)      const { return  fCPVFrameSize[index];      }

  // PPSD functions cannot be used for CPV

  Float_t GetAnodeThickness(void)          { AssertPPSD(); return 0; }
  Float_t GetAvalancheGap(void)            { AssertPPSD(); return 0; }
  Float_t GetCathodeThickness(void)        { AssertPPSD(); return 0; }
  Float_t GetCompositeThickness(void)      { AssertPPSD(); return 0; }
  Float_t GetConversionGap(void)           { AssertPPSD(); return 0; }
  Float_t GetLeadConverterThickness(void)  { AssertPPSD(); return 0; }
  Float_t GetLeadToMicro2Gap(void)         { AssertPPSD(); return 0; }
  Float_t GetLidThickness(void)            { AssertPPSD(); return 0; }
  Float_t GetMicromegas1Thickness(void)    { AssertPPSD(); return 0; }
  Float_t GetMicromegas2Thickness(void)    { AssertPPSD(); return 0; }
  Float_t GetMicromegasWallThickness(void) { AssertPPSD(); return 0; }
  Float_t GetMicro1ToLeadGap(void)         { AssertPPSD(); return 0; }
  Float_t GetPCThickness(void)             { AssertPPSD(); return 0; }
  Float_t GetPhiDisplacement(void)         { AssertPPSD(); return 0; }
  Float_t GetPPSDModuleSize(Int_t index)   { AssertPPSD(); return 0; }
  Float_t GetZDisplacement(void)           { AssertPPSD(); return 0; }
 
private:

  Int_t   fNumberOfCPVLayers;      // Number of CPV identical layers
  Bool_t  fLeadConverterExists;    // kTRUE if the lead converter between CPV layers exists
  Int_t   fNumberOfCPVPadsPhi;     // Number of CPV pads in phi
  Int_t   fNumberOfCPVPadsZ;       // Number of CPV pads in z
  Float_t fCPVPadSizePhi;          // CPV pad size in phi
  Float_t fCPVPadSizeZ;            // CPV pad size in z
  Float_t fCPVBoxSize[3];          // Outer size of CPV box
  Float_t fCPVActiveSize[2];       // Active size of CPV box (x,z)
  Int_t   fNumberOfCPVChipsPhi;    // Number of CPV Gassiplex chips in phi
  Int_t   fNumberOfCPVChipsZ;      // Number of CPV Gassiplex chips in z
  Float_t fGassiplexChipSize[3];   // Size of a Gassiplex chip (0 - in z, 1 - in phi, 2 - thickness (in ALICE radius))
  Float_t fCPVGasThickness;        // Thickness of CPV gas volume
  Float_t fCPVTextoliteThickness;  // Thickness of CPV textolite PCB (without foil)
  Float_t fCPVCuNiFoilThickness;   // Thickness of CPV Copper-Nickel foil of PCB
  Float_t fFTPosition[4];          // Positions of the 4 PCB vs the CPV box center
  Float_t fCPVFrameSize[3];        // CPV frame size (0 - in phi, 1 - in z, 2 - thickness (along ALICE radius))

  void    AssertPPSD() {
    printf("Function %s should not be called for CPV geometry\n",__PRETTY_FUNCTION__);
    assert(0==1) ;
  }

  ClassDef(AliPHOSCPVGeometry,1)       // CPV geometry base class 

} ;

#endif // AliPHOSCPVGEOMETRY_H
