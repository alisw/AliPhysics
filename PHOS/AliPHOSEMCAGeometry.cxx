/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : EMCA (Electromagnetic Calorimeter)  
// Its data members provide geometry parametrization of EMCA
// which can be changed in the constructor only.
// Author   : Yves Schutz (SUBATECH)
// Modified : Yuri Kharlov (IHEP, Protvino)
// 13 September 2000

// --- AliRoot header files ---

#include "AliPHOSEMCAGeometry.h"

ClassImp(AliPHOSEMCAGeometry) ;

//____________________________________________________________________________
AliPHOSEMCAGeometry::AliPHOSEMCAGeometry()
{

  // Initializes the EMC parameters

  fNPhi       = 64 ; 
  fNZ         = 64 ; 

  fXtlSize[0] =  2.2 ;
  fXtlSize[1] = 22.0 ;
  fXtlSize[2] =  2.2 ;

  // all these numbers coming next are subject to changes

  fOuterBoxThickness[0] = 2.5 ;
  fOuterBoxThickness[1] = 5.0 ;      
  fOuterBoxThickness[2] = 5.0 ;
  
  fUpperPlateThickness  = 4.0 ;
  
  fSecondUpperPlateThickness = 5.0 ; 
  
  fCrystalSupportHeight   = 6.95 ; 
  fCrystalWrapThickness   = 0.01 ;
  fCrystalHolderThickness = 0.005 ;
  fModuleBoxThickness     = 2.0 ; 
  fIPtoOuterCoverDistance = 447.0 ;      
  fIPtoCrystalSurface     = 460.0 ;  
  
  fPinDiodeSize[0] = 1.71 ;   //Values given by Odd Harald feb 2000  
  fPinDiodeSize[1] = 0.0280 ; // 0.0280 is the depth of active layer in the silicon     
  fPinDiodeSize[2] = 1.61 ;    
  
  fUpperCoolingPlateThickness   = 0.06 ; 
  fSupportPlateThickness        = 10.0 ;
  fLowerThermoPlateThickness    =  3.0 ; 
  fLowerTextolitPlateThickness  =  1.0 ;
  fGapBetweenCrystals           = 0.03 ;
  
  fTextolitBoxThickness[0] = 1.5 ;  
  fTextolitBoxThickness[1] = 0.0 ;   
  fTextolitBoxThickness[2] = 3.0 ; 
  
  fAirThickness[0] =  0.4   ;
  fAirThickness[1] = 20.5175 ;  
  fAirThickness[2] =  2.48   ;  
  
  Float_t xtalModulePhiSize =  fNPhi * ( fXtlSize[0] + 2 * fGapBetweenCrystals ) ; 
  Float_t xtalModuleZSize   =  fNZ   * ( fXtlSize[2] + 2 * fGapBetweenCrystals ) ;
  
  // The next dimensions are calculated from the above parameters
  
  fOuterBoxSize[0] =  xtalModulePhiSize + 2 * ( fAirThickness[0] + fModuleBoxThickness
						+ fTextolitBoxThickness[0] + fOuterBoxThickness[0] ) ; 
  fOuterBoxSize[1] = ( fXtlSize[1] + fCrystalSupportHeight + fCrystalWrapThickness + fCrystalHolderThickness )
    + 2 * (fAirThickness[1] +  fModuleBoxThickness + fTextolitBoxThickness[1] + fOuterBoxThickness[1] ) ;
  fOuterBoxSize[2] =  xtalModuleZSize   + 2 * ( fAirThickness[2] + fModuleBoxThickness 
						+ fTextolitBoxThickness[2] + fOuterBoxThickness[2] ) ; 
  
  fTextolitBoxSize[0]  = fOuterBoxSize[0] - 2 * fOuterBoxThickness[0] ;
  fTextolitBoxSize[1]  = fOuterBoxSize[1] -  fOuterBoxThickness[1] - fUpperPlateThickness ;
  fTextolitBoxSize[2]  = fOuterBoxSize[2] - 2 * fOuterBoxThickness[2] ;
  
  fAirFilledBoxSize[0] =  fTextolitBoxSize[0] - 2 * fTextolitBoxThickness[0] ; 
  fAirFilledBoxSize[1] =  fTextolitBoxSize[1] - fSecondUpperPlateThickness ; 
  fAirFilledBoxSize[2] =  fTextolitBoxSize[2] - 2 * fTextolitBoxThickness[2] ; 
  
}
//____________________________________________________________________________
