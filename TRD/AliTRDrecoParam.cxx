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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Parameter class for the TRD reconstruction                               //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDrecoParam.h"

ClassImp(AliTRDrecoParam)

//______________________________________________________________
AliTRDrecoParam::AliTRDrecoParam()
  :AliDetectorRecoParam()
  ,fkClusterSharing(0)
  ,fkPIDMethod(1) // LQ PID
  ,fkMaxTheta(1.0)
  ,fkMaxPhi(2.0)
  ,fkRoad0y(6.0)
  ,fkRoad0z(8.5) 
  ,fkRoad1y(2.0)
  ,fkRoad1z(20.0)	
  ,fkRoad2y(3.0)
  ,fkRoad2z(20.0)
  ,fkPlaneQualityThreshold(5.0)// 4.2? under Investigation
  ,fkFindable(.333)
  ,fkChi2Z(30./*14.*//*12.5*/)
  ,fkChi2Y(.25)
  ,fkTrackLikelihood(-15.)
{
  //
  // Default constructor
  //

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetLowFluxParam()
{
  //
  // Parameters for the low flux environment
  //

  return new AliTRDrecoParam();

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetHighFluxParam()
{
  //
  // Parameters for the high flux environment
  //

  return new AliTRDrecoParam();

}
