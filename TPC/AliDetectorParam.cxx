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
Revision 1.1.4.2  2000/04/10 11:36:13  kowal2

New Detector parameters handling class

*/

///////////////////////////////////////////////////////////////////////
//  Paramter class for AliDetector                                   //
//                                                                   //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliDetectorParam.h"




Float_t * AliDetectorParam::GetAnglesAccMomentum(Float_t *x, Int_t * index, Float_t *momentum, Float_t *angle)
{
  //
  //calculates deflection angle of particle with longitudinal
  //longitudinal  momentum[0] and transversal momentum momentum[1]
  //at position (x,y,z) = (x[0],x[1],x[2]) 
  //angle[0] - deep angle
  //angle[1] - magnetic deflection angle 
  if (momentum==0) {
    Float_t rtotal =TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
    if (rtotal==0) angle[0]=0;
    else    
      angle[0] = TMath::ATan(x[2]/rtotal);
    angle[1]=0;
    return angle;
  }
  Float_t mtotal =TMath::Sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]);
  if (mtotal==0) {
    angle[0]= 0;
    angle[1]=0;
    return angle;
  }
  angle[0]= TMath::ATan(momentum[2]/mtotal);
  Float_t radius1 = TMath::Sqrt(x[0]*x[0]+x[1]*x[1]); //axial symetry in z
  Float_t radius2 = 1000*mtotal/(3*fBField);
  if (radius1<radius2)
    angle[1]= TMath::ASin(radius1/radius2);
  else 
    angle[1]=0;
  return angle;
} 





ClassImp(AliDetectorParam)
