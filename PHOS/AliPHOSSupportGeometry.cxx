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
// Geometry class  for PHOS : Support which holds all PHOS modules.
// Its data members provide geometry parametrization of
// the PHOS support which can be changed in the constructor only.
// Author:   Yuri Kharlov (IHEP, Protvino)
// 13 November 2000

// --- AliRoot header files ---

#include "AliPHOSSupportGeometry.h"

ClassImp(AliPHOSSupportGeometry) ;

//____________________________________________________________________________
AliPHOSSupportGeometry::AliPHOSSupportGeometry()
{

  // Initializes the PHOS support parameters
  
  fRailLength   = 1200.0;

  fRailPart1[0] =   28.0;
  fRailPart1[1] =    3.0;
  fRailPart1[2] = fRailLength;

  fRailPart2[0] =    1.5;
  fRailPart2[1] =   54.0;
  fRailPart2[2] = fRailLength;

  fRailPart3[0] =    6.0;
  fRailPart3[1] =    5.0;
  fRailPart3[2] = fRailLength;

  fRailOuterSize[0] = fRailPart1[0];
  fRailOuterSize[1] = fRailPart1[1]*2 + fRailPart2[1] + fRailPart3[1];
  fRailOuterSize[2] = fRailLength;

  fDistanceBetwRails = 402.5;
  fRailsDistanceFromIP = 610.;

  fRailRoadSize[0] = fDistanceBetwRails + fRailOuterSize[0];
  fRailRoadSize[1] = fRailOuterSize[1];
  fRailRoadSize[2] = fRailOuterSize[2];

  fCradleWallThickness = 1.0;

  fCradleWall[0] =   0.;  // Inner radius, to be defined from PHOS parameters
  fCradleWall[1] =  65.;  // Diff. between outer and inner radii
  fCradleWall[2] =  18.;
  fCradleWall[3] = 270. - 50.;
  fCradleWall[4] = 270. + 50.;

  fCradleWheel[0] = 30.0;
  fCradleWheel[1] = 80.0;
  fCradleWheel[2] = 30.0;
  
}
//____________________________________________________________________________
