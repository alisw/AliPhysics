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

/*
  $Log$
  Revision 1.1  2000/10/25 08:24:33  schutz
  New CPV(GPS2) geometry class

*/

//_________________________________________________________________________
// Geometry class  for PHOS : PPSD (PHOS Preshower Detector)
//                  
//*-- Author   : Yves Schutz (SUBATECH)
//    Modified : Yuri Kharlov (IHEP, Protvino) 15 September 2000
//
// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPPSDGeometry.h"

ClassImp(AliPPSDGeometry) ;

//____________________________________________________________________________
AliPPSDGeometry::AliPPSDGeometry()
{

  // Initializes the PPSD parameters

  fAnodeThickness           = 0.0009 ; 
  fAvalancheGap             = 0.01 ; 
  fCathodeThickness         = 0.0009 ;
  fCompositeThickness       = 0.3 ; 
  fConversionGap            = 0.6 ; 
  fLeadConverterThickness   = 0.56 ; 
  fLeadToMicro2Gap          = 0.1 ; 
  fLidThickness             = 0.2 ; 
  fMicro1ToLeadGap          = 0.1 ; 
  fMicromegasWallThickness  = 0.6 ; 
  fNumberOfModulesPhi       = 4 ; 
  fNumberOfModulesZ         = 4 ; 
  fNumberOfPadsPhi          = 32 ; 
  fNumberOfPadsZ            = 32 ;   
  fPCThickness              = 0.1 ; 
  fPhiDisplacement          = 0.8 ;  
  fZDisplacement            = 0.8 ;  

  fMicromegas1Thickness     = fLidThickness + 2 * fCompositeThickness + fCathodeThickness 
                            + fPCThickness + fAnodeThickness + fConversionGap + fAvalancheGap ; 
  fMicromegas2Thickness     = fMicromegas1Thickness ; 

  fPPSDModuleSize[0]        = 38.0 ; 
  fPPSDModuleSize[1]        = fMicromegas1Thickness ; 
  fPPSDModuleSize[2]        = 38.0 ; 
 
  fPPSDBoxSize[0]           = fNumberOfModulesPhi * fPPSDModuleSize[0] + 2 * fPhiDisplacement ;  
  fPPSDBoxSize[1]           = fMicromegas2Thickness + fMicromegas2Thickness 
                            + fLeadConverterThickness + fMicro1ToLeadGap + fLeadToMicro2Gap ;    
  fPPSDBoxSize[2]           = fNumberOfModulesZ *  fPPSDModuleSize[2] + 2 * fZDisplacement ;

}

//____________________________________________________________________________
