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

//-----------------------------------------------------------------------------
/// \class AliMUONSurveyUtil
/// Singleton utility class for the survey processing of the ALICE DiMuon spectrometer 
/// 
/// This class contains various functions to calculate misalignement parameters
/// from survey data and designed positions of survey targets.
/// Macro also includes a method to get the new AliMUONGeometryTranformer.
/// 
/// \author Javier Castillo
//-----------------------------------------------------------------------------

#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TString.h>

#include "AliAlignObjMatrix.h"

#include "AliMUONSurveyUtil.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

/// \cond CLASSIMP
ClassImp(AliMUONSurveyUtil)
/// \endcond

int AliMUONSurveyUtil::fgNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
AliMUONSurveyUtil* AliMUONSurveyUtil::fgInstance(0x0);

AliMUONSurveyUtil::~AliMUONSurveyUtil(){
  printf("WHAT AM I DOING HERE????????????????\n");
  fgInstance = 0x0;
}

AliMUONSurveyUtil* AliMUONSurveyUtil::Instance() {
  ///  Return its instance 
  if (!fgInstance) 
    fgInstance = new AliMUONSurveyUtil();
  
  return fgInstance;
}

Bool_t AliMUONSurveyUtil::MatrixToAngles(const Double_t *rot, Double_t *angles)
{
  /// Calculates the Euler angles in "x y z" notation
  /// using the rotation matrix
  /// Returns false in case the rotation angles can not be

  // extracted from the matrix
  //
  if(TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) {
    printf("Failed to extract roll-pitch-yall angles!");
    return kFALSE;
  }
  //   Double_t raddeg = TMath::RadToDeg();
  angles[0]=TMath::ATan2(-rot[5],rot[8]);
  angles[1]=TMath::ASin(rot[2]);
  angles[2]=TMath::ATan2(-rot[1],rot[0]);
  return kTRUE;
}

void AliMUONSurveyUtil::AnglesToMatrix(const Double_t *angles, Double_t *rot)
{
  /// Calculates the rotation matrix using the 
  /// Euler angles in "x y z" notation
  ///
  //  Double_t degrad = TMath::DegToRad();
  Double_t degrad = 1.;
  Double_t sinpsi = TMath::Sin(degrad*angles[0]);
  Double_t cospsi = TMath::Cos(degrad*angles[0]);
  Double_t sinthe = TMath::Sin(degrad*angles[1]);
  Double_t costhe = TMath::Cos(degrad*angles[1]);
  Double_t sinphi = TMath::Sin(degrad*angles[2]);
  Double_t cosphi = TMath::Cos(degrad*angles[2]);

  rot[0] =  costhe*cosphi;
  rot[1] = -costhe*sinphi;
  rot[2] =  sinthe;
  rot[3] =  sinpsi*sinthe*cosphi + cospsi*sinphi;
  rot[4] = -sinpsi*sinthe*sinphi + cospsi*cosphi;
  rot[5] = -costhe*sinpsi;
  rot[6] = -cospsi*sinthe*cosphi + sinpsi*sinphi;
  rot[7] =  cospsi*sinthe*sinphi + sinpsi*cosphi;
  rot[8] =  costhe*cospsi;
}

Double_t AliMUONSurveyUtil::XpCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center x position using x coord. of 2 button targets. + solution. 
  Double_t lCos2Tht = TMath::Cos(2*par[6]);
  Double_t lSinTht = TMath::Sin(par[6]);

  Double_t inSqrt = TMath::Abs((par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] -x[1])*(x[0] -x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +((par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*lCos2Tht 
			       +4*(x[0] - x[1])*(par[2] - par[5])*lSinTht);

  if (inSqrt<0) return inSqrt*1e10;

  Double_t xD = ((2*(par[0]*par[0]*x[1] 
		     -par[0]*par[3]*(x[0] + x[1]) 
		     +x[1]*par[1]*(par[1] - par[4]) 
		     +x[0]*(par[3]*par[3] - par[1]*par[4] + par[4]*par[4])) 
		  -2*(par[3]*par[3]*par[2] 
		      +par[0]*par[0]*par[5] 
		      -par[0]*par[3]*(par[2] + par[5]) 
		      +(par[1] - par[4])*(-par[4]*par[2] +par[1]*par[5]))*lSinTht 
		  +TMath::Sqrt(2)*(-par[3]*par[1] + par[0]*par[4])
		  *TMath::Sqrt(inSqrt))
		 /(2*((par[0] - par[3])*(par[0] - par[3]) + (par[1] - par[4])*(par[1] - par[4]))));

  return xD;
}

Double_t AliMUONSurveyUtil::XnCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center x position using x coord. of 2 button targets. - solution. 
  Double_t lCos2Tht = TMath::Cos(2*par[6]);
  Double_t lSinTht = TMath::Sin(par[6]);

  Double_t inSqrt = TMath::Abs((par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +((par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*lCos2Tht
			       +4*(x[0] - x[1])*(par[2] - par[5])*lSinTht);

  if (inSqrt<0) return inSqrt*1e10;

  Double_t xD = ((2*(par[0]*par[0]*x[1] 
		     -par[0]*par[3]*(x[0] + x[1]) 
		     +x[1]*par[1]*(par[1] - par[4]) 
		     +x[0]*(par[3]*par[3] - par[1]*par[4] + par[4]*par[4])) 
		  -2*(par[3]*par[3]*par[2] + par[0]*par[0]*par[5] 
		      -par[0]*par[3]*(par[2] + par[5]) 
		      +(par[1] - par[4])*(-par[4]*par[2] + par[1]*par[5]))*lSinTht 
		  +TMath::Sqrt(2)*(par[3]*par[1] - par[0]*par[4])
		  *TMath::Sqrt(inSqrt))
		 /(2*((par[0] - par[3])*(par[0] - par[3]) + (par[1] - par[4])*(par[1] - par[4]))));

  return xD;
}

Double_t AliMUONSurveyUtil::PhiXpn(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using x coord. of 2 button targets. +- solution. 
  Double_t inSqrt = TMath::Abs(((par[0] - par[3])*(par[0] - par[3]) 
				-2*(x[0] - x[1])*(x[0] - x[1]) 
				+(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
				+(+(par[0] - par[3])*(par[0] - par[3]) 
				  +(par[1] - par[4])*(par[1] - par[4]) 
				  +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
				+4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6])));

  if (inSqrt<0) return inSqrt*1e10;
  
  Double_t phix = ((+2*(par[0] - par[3])*(x[0] - x[1]) 
		    -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		    +TMath::Sqrt(2)*(par[1] - par[4])
		    *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = -TMath::ACos(phix);

  return phix;
}

Double_t AliMUONSurveyUtil::PhiXpp(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using x coord. of 2 button targets. ++ solution. 
  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       +4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));

  if (inSqrt<0) return inSqrt*1e10;

  Double_t phix = ((+2*(par[0] - par[3])*(x[0] - x[1]) 
		    -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		    +TMath::Sqrt(2)*(par[1] - par[4])
		    *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = TMath::ACos(phix);

  return phix;
}

Double_t AliMUONSurveyUtil::PhiXnn(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using x coord. of 2 button targets. -- solution. 
  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       + 4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));

  if (inSqrt<0) return inSqrt*1e10;
  
  Double_t phix = (+(+2*(par[0] - par[3])*(x[0] - x[1]) 
		     -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		     +TMath::Sqrt(2)*(-par[1] + par[4])
		     *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = -TMath::ACos(phix);

  return phix;
}

Double_t AliMUONSurveyUtil::PhiXnp(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using x coord. of 2 button targets. +- solution. 
  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       +4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));

  if (inSqrt<0) return inSqrt*1e10;

  Double_t phix = (+(+2*(par[0] - par[3])*(x[0] - x[1]) 
		     -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		     +TMath::Sqrt(2)*(-par[1] + par[4])
		     *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = TMath::ACos(phix);

  return phix;
}

Double_t AliMUONSurveyUtil::YpCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center y position using y coord. of 2 button targets. + solution. 

  // par : x1l, y1l, z1l, x2l, y2l, z2, lpsi, tht,
  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t yD = ((1./((par[0] - par[3])*(par[0] - par[3]) + (par[1] - par[4])*(par[1] - par[4])))
		 *(+par[3]*par[3]*x[0] 
		   +par[0]*par[0]*x[1] 
		   -par[0]*par[3]*(x[0] + x[1]) 
		   +(par[1] - par[4])*(-x[0]*par[4] + par[1]*x[1]) 
		   +(par[3]*par[3]*par[2] 
		     +par[0]*par[0]*par[5] 
		     -par[0]*par[3]*(par[2] + par[5]) 
		     +(par[1] - par[4])*(-par[4]*par[2] + par[1]*par[5]))*lCosTht*lSinPsi 
		   +(-par[3]*par[1] + par[0]*par[4])
		   *TMath::Sqrt(-(x[0] - x[1] 
				  +(par[2] - par[5])*lCosTht*lSinPsi)
				*(x[0] - x[1] 
				  +(par[2] - par[5])*lCosTht*lSinPsi) 
				+ ((par[0] - par[3])*(par[0] - par[3]) 
				   +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
									  +lSinPsi*lSinPsi*lSinTht*lSinTht))));

  return yD;  
}

Double_t AliMUONSurveyUtil::PhiYpn(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using y coord. of 2 button targets. +- solution. 

  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t phiy = ((lCosPsi*((par[1] - par[4])*(x[0] - x[1]) 
			     +(par[1] - par[4])*(par[2] - par[5])*lCosTht*lSinPsi 
			     +(-par[0] + par[3])
			     *TMath::Sqrt(-(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi)
					  *(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi) 
					  +(+(par[0] - par[3])*(par[0] - par[3]) 
					    +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi
										   +lSinPsi*lSinPsi*lSinTht*lSinTht))) 
		    +lSinPsi*lSinTht*((par[0] - par[3])*(x[0] - x[1]) 
				      +(par[0] - par[3])*(par[2] - par[5])*lCosTht*lSinPsi 
				      +(par[1] - par[4])
				      *TMath::Sqrt(-(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi)
						   *(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi) 
						   + ((par[0] - par[3])*(par[0] - par[3]) 
						      +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
											     +lSinPsi*lSinPsi*lSinTht*lSinTht))))
		   /((+(par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
							     +lSinPsi*lSinPsi*lSinTht*lSinTht)));
  
  phiy = -TMath::ACos(phiy);


  return phiy;
}

Double_t AliMUONSurveyUtil::PhiYpp(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using y coord. of 2 button targets. ++ solution. 

  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t phiy = ((lCosPsi*((par[1] - par[4])*(x[0] - x[1]) 
			     +(par[1] - par[4])*(par[2] - par[5])*lCosTht*lSinPsi 
			     +(-par[0] + par[3])
			     *TMath::Sqrt(-(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi)
					  *(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi) 
					  +((par[0] - par[3])*(par[0] - par[3]) 
					    +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
										   +lSinPsi*lSinPsi*lSinTht*lSinTht))) 
		    +lSinPsi*lSinTht*((par[0] - par[3])*(x[0] - x[1]) 
				      +(par[0] - par[3])*(par[2] - par[5])*lCosTht*lSinPsi 
				      +(par[1] - par[4])*TMath::Sqrt(-(x[0] - x[1] 
								       +(par[2] - par[5])*lCosTht*lSinPsi)
								     *(x[0] - x[1] 
								       +(par[2] - par[5])*lCosTht*lSinPsi) 
								     +((par[0] - par[3])*(par[0] - par[3])
								       +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
													      +lSinPsi*lSinPsi*lSinTht*lSinTht))))
		   /(((par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4]))*(lCosPsi*lCosPsi 
							     +lSinPsi*lSinPsi*lSinTht*lSinTht)));
  
  phiy = TMath::ACos(phiy);

  return phiy;
}
 
Double_t AliMUONSurveyUtil::YnCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center y position using y coord. of 2 button targets. - solution. 
  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t yD = ((1./(+(par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4])))
		 *(+par[3]*par[3]*x[0] 
		   +par[0]*par[0]*x[1] 
		   -par[0]*par[3]*(x[0] + x[1]) 
		   +(par[1] - par[4])*(-x[0]*par[4] + par[1]*x[1]) 
		   +(+par[3]*par[3]*par[2] 
		     +par[0]*par[0]*par[5] 
		     -par[0]*par[3]*(par[2] + par[5]) 
		     +(par[1] - par[4])*(-par[4]*par[2] + par[1]*par[5]))*lCosTht*lSinPsi 
		   +(par[3]*par[1] - par[0]*par[4])
		   *TMath::Sqrt(-(+x[0] - x[1] 
				  +(par[2] - par[5])*lCosTht*lSinPsi)
				*(x[0] - x[1] 
				  +(par[2] - par[5])*lCosTht*lSinPsi) 
				+((par[0] - par[3])*(par[0] - par[3]) 
				  +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
									 +lSinPsi*lSinPsi*lSinTht*lSinTht))));
  
  return yD;  
}
 

Double_t AliMUONSurveyUtil::PhiYnn(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using y coord. of 2 button targets. -- solution. 

  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t phiy = ((lCosPsi*(+(par[1] - par[4])*(x[0] - x[1]) 
			     +(par[1] - par[4])*(par[2] - par[5])*lCosTht*lSinPsi 
			     +(par[0] - par[3])
			     *TMath::Sqrt(-(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi)
					  *(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi) 
					  +(+(par[0] - par[3])*(par[0] - par[3]) 
					    +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
										   +lSinPsi*lSinPsi*lSinTht*lSinTht))) 
		    +lSinPsi*lSinTht*(+(par[0] - par[3])*(x[0] - x[1]) 
				      +(par[0] - par[3])*(par[2] - par[5])*lCosTht*lSinPsi 
				      +(-par[1] + par[4])
				      *TMath::Sqrt(-(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi)
						   *(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi) 
						   +(+(par[0] - par[3])*(par[0] - par[3]) 
						     +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
											    +lSinPsi*lSinPsi*lSinTht*lSinTht))))
		   /((+(par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
							     +lSinPsi*lSinPsi*lSinTht*lSinTht)));
  
  phiy = -TMath::ACos(phiy);
  
  return phiy;
}


Double_t AliMUONSurveyUtil::PhiYnp(const Double_t *x, const Double_t *par) const{
  /// Returns phi rot. using y coord. of 2 button targets. -+ solution. 

  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t phiy = ((lCosPsi*(+(par[1] - par[4])*(x[0] - x[1]) 
			     +(par[1] - par[4])*(par[2] - par[5])*lCosTht*lSinPsi 
			     +(par[0] - par[3])
			     *TMath::Sqrt(-(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi)
					  *(x[0] - x[1] 
					    +(par[2] - par[5])*lCosTht*lSinPsi) 
					  +((par[0] - par[3])*(par[0] - par[3]) 
					    +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
										   +lSinPsi*lSinPsi*lSinTht*lSinTht))) 
		    +lSinPsi*lSinTht*(+(par[0] - par[3])*(x[0] - x[1]) 
				      +(par[0] - par[3])*(par[2] - par[5])*lCosTht*lSinPsi 
				      +(-par[1] + par[4])
				      *TMath::Sqrt(-(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi)
						   *(x[0] - x[1] 
						     +(par[2] - par[5])*lCosTht*lSinPsi) 
						   +((par[0] - par[3])*(par[0] - par[3]) 
						     +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
											    +lSinPsi*lSinPsi*lSinTht*lSinTht))))
		   /((+(par[0] - par[3])*(par[0] - par[3])
		      +(par[1] - par[4])*(par[1] - par[4]))*(+lCosPsi*lCosPsi 
							     +lSinPsi*lSinPsi*lSinTht*lSinTht)));
  
  phiy = TMath::ACos(phiy);
  
  return phiy;
}

Double_t AliMUONSurveyUtil::ZnCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center z position using z coord. of 2 button targets. - solution. 

  // par :  x1l, y1l, z1l, x2l, y2l, z2l, psi, tht
  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t inSqrt = ((par[3]*par[1] - par[0]*par[4])*(par[3]*par[1] - par[0]*par[4])
		     *((-(x[0] - x[1])*(x[0] - x[1]))
		       +(((par[0] - par[3])*(par[0] - par[3])
			  +(par[1] - par[4])*(par[1] - par[4])))*lSinPsi*lSinPsi
		       +lCosPsi*((-(par[2] - par[5]))
				 *lCosTht*(-2*x[0]+2*x[1]
					   +(par[2] - par[5])*lCosPsi*lCosTht)
				 +((par[0] - par[3])*(par[0] - par[3])
				   +(par[1] - par[4])*(par[1] - par[4]))*lCosPsi*lSinTht*lSinTht)));

  if (inSqrt<0) return inSqrt*1e10;

  Double_t zD = ((1./((par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4])))
		 *(-par[1]*par[4]*x[0] 
		   +par[4]*par[4]*x[0] 
		   +par[0]*par[0]*x[1] 
		   +par[1]*par[1]*x[1] 
		   -par[1]*par[4]*x[1] 
		   -par[0]*par[3]*(x[0] + x[1]) 
		   +par[3]*par[3]*x[0]
		   +(+par[1]*par[4]*par[2] 
		     -par[4]*par[4]*par[2] 
		     -par[0]*par[0]*par[5] 
		     -par[1]*par[1]*par[5] 
		     +par[1]*par[4]*par[5] 
		     +par[0]*par[3]*(par[2] + par[5])
		     -par[3]*par[3]*par[2])*lCosPsi*lCosTht
		   -TMath::Sqrt(inSqrt)));
  
  return zD;
}

Double_t AliMUONSurveyUtil::ZpCenter(const Double_t *x, const Double_t *par) const{
  /// Returns center z position using z coord. of 2 button targets. + solution. 

  // par :  x1l, y1l, z1l, x2l, y2l, z2l, psi, tht
  Double_t lCosPsi = TMath::Cos(par[6]);
  Double_t lSinPsi = TMath::Sin(par[6]);
  Double_t lCosTht = TMath::Cos(par[7]);
  Double_t lSinTht = TMath::Sin(par[7]);

  Double_t inSqrt = ((par[3]*par[1] - par[0]*par[4])*(par[3]*par[1] - par[0]*par[4])
		     *((-(x[0] - x[1])*(x[0] - x[1]))
		       +(((par[0] - par[3])*(par[0] - par[3])
			  +(par[1] - par[4])*(par[1] - par[4])))*lSinPsi*lSinPsi
		       +lCosPsi*((-(par[2] - par[5]))
				 *lCosTht*(-2*x[0]+2*x[1]
					   +(par[2] - par[5])*lCosPsi*lCosTht)
				 +((par[0] - par[3])*(par[0] - par[3])
				   +(par[1] - par[4])*(par[1] - par[4]))*lCosPsi*lSinTht*lSinTht)));

  if (inSqrt<0) return inSqrt*1e10;  

  Double_t zD = ((1./((par[0] - par[3])*(par[0] - par[3]) 
		      +(par[1] - par[4])*(par[1] - par[4])))
		 *(-par[1]*par[4]*x[0] 
		   +par[4]*par[4]*x[0] 
		   +par[0]*par[0]*x[1] 
		   +par[1]*par[1]*x[1] 
		   -par[1]*par[4]*x[1] 
		   -par[0]*par[3]*(x[0] + x[1]) 
		   +par[3]*par[3]*x[0]
		   +(+par[1]*par[4]*par[2] 
		     -par[4]*par[4]*par[2] 
		     -par[0]*par[0]*par[5] 
		     -par[1]*par[1]*par[5] 
		     +par[1]*par[4]*par[5] 
		     +par[0]*par[3]*(par[2] + par[5])
		     -par[3]*par[3]*par[2])*lCosPsi*lCosTht
		   +TMath::Sqrt(inSqrt)));

  return zD;
}

//______________________________________________________________________
AliMUONGeometryTransformer* AliMUONSurveyUtil::ReAlign(const AliMUONGeometryTransformer * transformer, 
						       int rMod, int rNDetElems, int rDetElemPseudoIdToDetElem[], TGeoCombiTrans deltaDetElemTransf[], Bool_t verbose)
{
  /////////////////////////////////////////////////////////////////////
  ///   Takes the internal geometry module transformers, copies them
  /// and gets the Detection Elements from them.
  /// Takes misalignment parameters and applies these
  /// to the local transform of the Detection Element
  /// Obtains the global transform by multiplying the module transformer
  /// transformation with the local transformation 
  /// Applies the global transform to a new detection element
  /// Adds the new detection element to a new module transformer
  /// Adds the new module transformer to a new geometry transformer
  /// Returns the new geometry transformer

  Int_t iDetElemId = 0;
  Int_t iDetElemPseudoId = 0;

  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer();
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++) {
    // module transformers    
    const AliMUONGeometryModuleTransformer *kModuleTransformer =
      transformer->GetModuleTransformer(iMt, true);
      
    AliMUONGeometryModuleTransformer *newModuleTransformer =
      new AliMUONGeometryModuleTransformer(iMt);
    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    
    TGeoCombiTrans moduleTransform =
      TGeoCombiTrans(*kModuleTransformer->GetTransformation());
    // New module transformation
    TGeoCombiTrans *newModuleTransform;
    if ((rMod<4 && iMt==rMod) || (rMod>=4 && (iMt==4+(rMod-4)*2||iMt==4+(rMod-4)*2+1))) {
      newModuleTransform = new TGeoCombiTrans(moduleTransform*deltaDetElemTransf[rNDetElems]);
    } else {
      newModuleTransform = new TGeoCombiTrans(moduleTransform);
    }
    newModuleTransformer->SetTransformation(*newModuleTransform);
    
    // For the selected chamber add misalign module
    if ((rMod<4 && iMt==rMod) || (rMod>=4 && (iMt==4+(rMod-4)*2||iMt==4+(rMod-4)*2+1))) {
      // Get delta transformation: 
      // Tdelta = Tnew * Told.inverse
      TGeoHMatrix deltaModuleTransform = 
	AliMUONGeometryBuilder::Multiply(*newModuleTransform, 
					 kModuleTransformer->GetTransformation()->Inverse());    
      // Create module mis alignment matrix
      newGeometryTransformer
	->AddMisAlignModule(kModuleTransformer->GetModuleId(), deltaModuleTransform);
    }
    
    AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();
    
    if (verbose)
      printf("%i DEs in old GeometryStore  %i\n",detElements->GetSize(), iMt);
    TGeoCombiTrans *deltaLocalTransform;
    TIter next(detElements->CreateIterator());
    AliMUONGeometryDetElement *detElement;
    while ((detElement = static_cast<AliMUONGeometryDetElement*>(next()))){
      /// make a new detection element
      AliMUONGeometryDetElement *newDetElement =
	new AliMUONGeometryDetElement(detElement->GetId(),
				      detElement->GetVolumePath());
      TString lDetElemName(detElement->GetDEName());
      lDetElemName.ReplaceAll("DE","");
      iDetElemId = lDetElemName.Atoi();
      iDetElemPseudoId = iDetElemId%100;
      if ((rMod<4 && iMt==rMod) || (rMod>=4 && (iMt==4+(rMod-4)*2||iMt==4+(rMod-4)*2+1))) {
	deltaLocalTransform = new TGeoCombiTrans(deltaDetElemTransf[rDetElemPseudoIdToDetElem[iDetElemPseudoId]]);       
      } else {
	deltaLocalTransform = new TGeoCombiTrans(*gGeoIdentity);
      }
      
      // local transformation of this detection element.
      TGeoCombiTrans localTransform
	= TGeoCombiTrans(*detElement->GetLocalTransformation());
      //      TGeoHMatrix newLocalMatrix = localTransform * (*deltaLocalTransform);
      TGeoCombiTrans newLocalTransform 
	= TGeoCombiTrans(localTransform * (*deltaLocalTransform));
      newDetElement->SetLocalTransformation(newLocalTransform);	  
      // global transformation
      TGeoHMatrix newGlobalTransform =
	AliMUONGeometryBuilder::Multiply(*newModuleTransform,
					 newLocalTransform);
      newDetElement->SetGlobalTransformation(newGlobalTransform);
      
      // add this det element to module
      newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
						      newDetElement);
      
      // In the Alice Alignment Framework misalignment objects store
      // global delta transformation
      // Get detection "intermediate" global transformation
      TGeoHMatrix newOldGlobalTransform = (*newModuleTransform) * localTransform;
      // Get detection element global delta transformation: 
      // Tdelta = Tnew * Told.inverse
      TGeoHMatrix  deltaGlobalTransform
	= AliMUONGeometryBuilder::Multiply(newGlobalTransform, 
					   newOldGlobalTransform.Inverse());
      
      // Create mis alignment matrix
      newGeometryTransformer
	->AddMisAlignDetElement(detElement->GetId(), deltaGlobalTransform);
    }
      
    if (verbose)
      printf("Added module transformer %i to the transformer\n", iMt);
    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
  }
  return newGeometryTransformer;
}

void AliMUONSurveyUtil::SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY){
  /// Sets the alignment resolution to the AliAlignObjMatrix correlation matrix

  TMatrixDSym mChCorrMatrix(6);
  mChCorrMatrix[0][0]=chResX*chResX;
  mChCorrMatrix[1][1]=chResY*chResY;
  //  mChCorrMatrix.Print();

  TMatrixDSym mDECorrMatrix(6);
  mDECorrMatrix[0][0]=deResX*deResX;
  mDECorrMatrix[1][1]=deResY*deResY;
  //  mDECorrMatrix.Print();

  AliAlignObjMatrix *alignMat = 0x0;

//  Int_t modId = (chId<4)? chId : 4+(chId-4)*2;   
  TString chName1;
  TString chName2;
  if (chId<4){
    chName1 = Form("GM%d",chId);
    chName2 = Form("GM%d",chId);
  } else {
    chName1 = Form("GM%d",4+(chId-4)*2);
    chName2 = Form("GM%d",4+(chId-4)*2+1);
  }
  
  for (int i=0; i<misAlignArray->GetEntries(); i++) {
    alignMat = (AliAlignObjMatrix*)misAlignArray->At(i);
    TString volName(alignMat->GetSymName());
    if((volName.Contains(chName1)&&
	volName.Last('/')<=volName.Index(chName1)+chName1.Length())||
       (volName.Contains(chName2)&&
	volName.Last('/')<=volName.Index(chName2)+chName2.Length())) {
      volName.Remove(0,volName.Last('/')+1);
      if (volName.Contains("GM")) {
	//	alignMat->Print("NULL");
	alignMat->SetCorrMatrix(mChCorrMatrix);
      } else if (volName.Contains("DE")) {
	//	alignMat->Print("NULL");
	alignMat->SetCorrMatrix(mDECorrMatrix);
      }
    }
  }
}
