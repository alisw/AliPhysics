#ifndef ALIPHOSPPSDGEOMETRY_H
#define ALIPHOSPPSDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry derived class for PHOS:PPSD (PHOS Preshower Detector)
// Its data members provide geometry parametrization of PPSD
// which can be changed in the constructor only.
// Author  : Yves Schutz
// Modified: Yuri Kharlov (IHEP, Protvino)
// 7 November 2000

#include <assert.h> 

#include "AliPHOSCPVBaseGeometry.h"

class AliPHOSPPSDGeometry : public AliPHOSCPVBaseGeometry {

public: 

           AliPHOSPPSDGeometry();
  virtual ~AliPHOSPPSDGeometry(void) {};

  // PPSD functions

  virtual Float_t GetAnodeThickness(void)          { return  fAnodeThickness ;          }
  virtual Float_t GetAvalancheGap(void)            { return  fAvalancheGap ;            }
  virtual Float_t GetCathodeThickness(void)        { return  fCathodeThickness ;        }
  virtual Float_t GetCompositeThickness(void)      { return  fCompositeThickness ;      }
  virtual Float_t GetConversionGap(void)           { return  fConversionGap ;           }
  virtual Float_t GetLeadConverterThickness(void)  { return  fLeadConverterThickness ;  }
  virtual Float_t GetLeadToMicro2Gap(void)         { return  fLeadToMicro2Gap ;         }
  virtual Float_t GetLidThickness(void)            { return  fLidThickness ;            }
  virtual Float_t GetMicromegas1Thickness(void)    { return  fMicromegas1Thickness ;    }
  virtual Float_t GetMicromegas2Thickness(void)    { return  fMicromegas2Thickness ;    }
  virtual Float_t GetMicromegasWallThickness(void) { return  fMicromegasWallThickness ; }
  virtual Float_t GetMicro1ToLeadGap(void)         { return  fMicro1ToLeadGap ;         }
  virtual Int_t   GetNumberOfPadsPhi(void)         { return  fNumberOfPadsPhi ;         }
  virtual Int_t   GetNumberOfPadsZ(void)           { return  fNumberOfPadsZ ;           }
  virtual Int_t   GetNumberOfModulesPhi(void)      { return  fNumberOfModulesPhi ;      }
  virtual Int_t   GetNumberOfModulesZ(void)        { return  fNumberOfModulesZ ;        }
  virtual Float_t GetPCThickness(void)             { return  fPCThickness ;             }
  virtual Float_t GetPhiDisplacement(void)         { return  fPhiDisplacement ;         }
  virtual Float_t GetCPVBoxSize(Int_t index)       { return  fPPSDBoxSize[index] ;      }
  virtual Float_t GetPPSDModuleSize(Int_t index)   { return  fPPSDModuleSize[index] ;   }
  virtual Float_t GetZDisplacement(void)           { return  fZDisplacement ;           }
  
  // CPV functions cannot be used for PPSD
  
  virtual Int_t   GetNumberOfCPVLayers(void)       { AssertCPV("GetNumberOfCPVLayers");     return 0; }
  virtual Bool_t  IsLeadConverterExists(void)      { AssertCPV("IsLeadConverterExists");    return 0; }
  virtual Float_t GetCPVActiveSize(Int_t index)    { AssertCPV("GetCPVActiveSize");         return 0; }
  virtual Int_t   GetNumberOfCPVChipsPhi(void)     { AssertCPV("GetNumberOfCPVChipsPhi");   return 0; }
  virtual Int_t   GetNumberOfCPVChipsZ(void)       { AssertCPV("GetNumberOfCPVChipsZ");     return 0; }
  virtual Float_t GetGassiplexChipSize(Int_t index){ AssertCPV("GetGassiplexChipSize");     return 0; }
  virtual Float_t GetCPVGasThickness(void)         { AssertCPV("GetCPVGasThickness");       return 0; }
  virtual Float_t GetCPVTextoliteThickness(void)   { AssertCPV("GetCPVTextoliteThickness"); return 0; }
  virtual Float_t GetCPVCuNiFoilThickness(void)    { AssertCPV("GetCPVCuNiFoilThickness");  return 0; }
  virtual Float_t GetFTPosition(Int_t index)       { AssertCPV("GetFTPosition");            return 0; }
  virtual Float_t GetCPVFrameSize(Int_t index)     { AssertCPV("GetCPVFrameSize");          return 0; }
  virtual Float_t GetIPtoCPVDistance(void)         { AssertCPV("GetIPtoCPVDistance");       return 0; }

private:

  Float_t fAnodeThickness ;          // Thickness of the copper layer which makes the anode 
  Float_t fAvalancheGap ;            // Thickness of the gas in the avalanche stage
  Float_t fCathodeThickness ;        // Thickness of composite material ensuring rigidity of cathode
  Float_t fCompositeThickness ;      // Thickness of composite material ensuring rigidity of anode
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

  void    AssertCPV(char* name) {
    printf("Function AliPPSDGeometry::%s should not be called for PPSD geometry\n",name);
    assert(0==1) ;
  }

  ClassDef(AliPHOSPPSDGeometry,1)        // PPSD geometry class 

} ;

#endif // AliPHOSPPSDGEOMETRY_H
