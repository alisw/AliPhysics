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

#include <assert.h> 

#include "AliPHOSCPVBaseGeometry.h"

class AliPHOSCPVGeometry : public AliPHOSCPVBaseGeometry {

public: 

           AliPHOSCPVGeometry();
  virtual ~AliPHOSCPVGeometry(void) {};

  // CPV functions

  virtual Int_t   GetNumberOfCPVLayers(void)        { return  fNumberOfCPVLayers;        }
  virtual Bool_t  IsLeadConverterExists(void)       { return  fLeadConverterExists;      }
  virtual Int_t   GetNumberOfCPVPadsPhi(void)       { return  fNumberOfCPVPadsPhi ;      }
  virtual Int_t   GetNumberOfCPVPadsZ(void)         { return  fNumberOfCPVPadsZ ;        }
  virtual Float_t GetCPVPadSizePhi(void)            { return  fCPVPadSizePhi;            }
  virtual Float_t GetCPVPadSizeZ(void)              { return  fCPVPadSizeZ;              }
  virtual Float_t GetCPVBoxSize(Int_t index)        { return  fCPVBoxSize[index];        }
  virtual Float_t GetCPVActiveSize(Int_t index)     { return  fCPVActiveSize[index];     }
  virtual Int_t   GetNumberOfCPVChipsPhi(void)      { return  fNumberOfCPVChipsPhi;      }
  virtual Int_t   GetNumberOfCPVChipsZ(void)        { return  fNumberOfCPVChipsZ;        }
  virtual Float_t GetGassiplexChipSize(Int_t index) { return  fGassiplexChipSize[index]; }
  virtual Float_t GetCPVGasThickness(void)          { return  fCPVGasThickness;          }
  virtual Float_t GetCPVTextoliteThickness(void)    { return  fCPVTextoliteThickness;    }
  virtual Float_t GetCPVCuNiFoilThickness(void)     { return  fCPVCuNiFoilThickness;     }
  virtual Float_t GetFTPosition(Int_t index)        { return  fFTPosition[index];        }
  virtual Float_t GetCPVFrameSize(Int_t index)      { return  fCPVFrameSize[index];      }

  // PPSD functions cannot be used for CPV

  virtual Float_t GetAnodeThickness(void)          { AssertPPSD("GetAnodeThickness");          return 0; }
  virtual Float_t GetAvalancheGap(void)            { AssertPPSD("GetAvalancheGap");            return 0; }
  virtual Float_t GetCathodeThickness(void)        { AssertPPSD("GetCathodeThickness");        return 0; }
  virtual Float_t GetCompositeThickness(void)      { AssertPPSD("GetCompositeThickness");      return 0; }
  virtual Float_t GetConversionGap(void)           { AssertPPSD("GetConversionGap");           return 0; }
  virtual Float_t GetLeadConverterThickness(void)  { AssertPPSD("GetLeadConverterThickness");  return 0; }
  virtual Float_t GetLeadToMicro2Gap(void)         { AssertPPSD("GetLeadToMicro2Gap");         return 0; }
  virtual Float_t GetLidThickness(void)            { AssertPPSD("GetLidThickness");            return 0; }
  virtual Float_t GetMicromegas1Thickness(void)    { AssertPPSD("GetMicromegas1Thickness");    return 0; }
  virtual Float_t GetMicromegas2Thickness(void)    { AssertPPSD("GetMicromegas2Thickness");    return 0; }
  virtual Float_t GetMicromegasWallThickness(void) { AssertPPSD("GetMicromegasWallThickness"); return 0; }
  virtual Float_t GetMicro1ToLeadGap(void)         { AssertPPSD("GetMicro1ToLeadGap");         return 0; }
  virtual Float_t GetPCThickness(void)             { AssertPPSD("GetPCThickness");             return 0; }
  virtual Float_t GetPhiDisplacement(void)         { AssertPPSD("GetPhiDisplacement");         return 0; }
  virtual Float_t GetPPSDModuleSize(Int_t index)   { AssertPPSD("GetPPSDModuleSize");          return 0; }
  virtual Float_t GetZDisplacement(void)           { AssertPPSD("GetZDisplacement");           return 0; }
  virtual Int_t   GetNumberOfPadsPhi(void)         { AssertPPSD("GetNumberOfPadsPhi");         return 0; }
  virtual Int_t   GetNumberOfPadsZ(void)           { AssertPPSD("GetNumberOfPadsZ");           return 0; }
  virtual Int_t   GetNumberOfModulesPhi(void)      { AssertPPSD("GetNumberOfModulesPhi");      return 0; }
  virtual Int_t   GetNumberOfModulesZ(void)        { AssertPPSD("GetNumberOfModulesZ");        return 0; }
  virtual void    SetLeadConverterThickness(Float_t x)  { AssertPPSD("SetLeadConverterThickness");       }
 
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

  void    AssertPPSD(char* name) {
    printf("Function AliCPVGeometry::%s should not be called for CPV geometry\n",name);
    assert(0==1) ;
  }

  ClassDef(AliPHOSCPVGeometry,1)       // CPV geometry base class 

} ;

#endif // AliPHOSCPVGEOMETRY_H
