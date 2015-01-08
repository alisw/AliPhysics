#ifndef ALIPHOSCPVGEOMETRY_H
#define ALIPHOSCPVGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry derived class for PHOS:CPV (Charged particle veto, IHEP version)
// Its data members provide geometry parametrization of CPV
// which can be changed in the constructor only.
// Author:  Yuri Kharlov (IHEP, Protvino)
// 7 November 2000

//#include "AliPHOSCPVBaseGeometry.h"
#include "TObject.h"

class AliPHOSCPVGeometry : public TObject {

public: 

           AliPHOSCPVGeometry();
  virtual ~AliPHOSCPVGeometry(void) {};

  // CPV functions

  virtual Int_t   GetNumberOfCPVLayers(void)        const { return  fNumberOfCPVLayers;        }
  virtual Int_t   GetNumberOfCPVPadsPhi(void)       const { return  fNumberOfCPVPadsPhi ;      }
  virtual Int_t   GetNumberOfCPVPadsZ(void)         const { return  fNumberOfCPVPadsZ ;        }
  virtual Float_t GetCPVPadSizePhi(void)            const { return  fCPVPadSizePhi;            }
  virtual Float_t GetCPVPadSizeZ(void)              const { return  fCPVPadSizeZ;              }
  virtual Float_t GetCPVBoxSize(Int_t index)        const { return  fCPVBoxSize[index];        }
  virtual Float_t GetCPVActiveSize(Int_t index)     const { return  fCPVActiveSize[index];     }
  virtual Int_t   GetNumberOfCPVChipsPhi(void)      const { return  fNumberOfCPVChipsPhi;      }
  virtual Int_t   GetNumberOfCPVChipsZ(void)        const { return  fNumberOfCPVChipsZ;        }
  virtual Float_t GetGassiplexChipSize(Int_t index) const { return  fGassiplexChipSize[index]; }
  virtual Float_t GetCPVGasThickness(void)          const { return  fCPVGasThickness;          }
  virtual Float_t GetCPVTextoliteThickness(void)    const { return  fCPVTextoliteThickness;    }
  virtual Float_t GetCPVCuNiFoilThickness(void)     const { return  fCPVCuNiFoilThickness;     }
  virtual Float_t GetFTPosition(Int_t index)        const { return  fFTPosition[index];        }
  virtual Float_t GetCPVFrameSize(Int_t index)      const { return  fCPVFrameSize[index];      }

 
private:

  Int_t   fNumberOfCPVLayers;      // Number of CPV identical layers
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


  ClassDef(AliPHOSCPVGeometry,2)       // CPV geometry base class 

} ;

#endif // AliPHOSCPVGEOMETRY_H
