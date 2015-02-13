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

/// \ingroup macros
/// \file MUONSurveyUtil.C
/// \brief Utility macro for survey data to alignment transformation.
///  
/// Macro contains various functions to calculate misalignement parameters
/// from survey data and designed positions of survey targets.
/// Macro also includes a method to get the new AliMUONGeometryTransformer.
/// It is intended to be loaded by chamber specific macros.
/// 
/// \author Javier Castillo

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryMisAligner.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>

#include <fstream>

#endif

static int fgNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};

Bool_t MatrixToAngles(const Double_t *rot, Double_t *angles)
{
  // Calculates the Euler angles in "x y z" notation
  // using the rotation matrix
  // Returns false in case the rotation angles can not be

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

Double_t eqPlane(Double_t *x, Double_t *par){
  return (-par[0]*x[0] -par[1]*x[1] -par[2]); 
}

Double_t xpCenter(Double_t *x, Double_t *par){

  Double_t lCos2Tht = TMath::Cos(2*par[6]);
  Double_t lSinTht = TMath::Sin(par[6]);

  Double_t inSqrt = TMath::Abs((par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] -x[1])*(x[0] -x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +((par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*lCos2Tht 
			       +4*(x[0] - x[1])*(par[2] - par[5])*lSinTht);

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

Double_t xnCenter(Double_t *x, Double_t *par){

  Double_t lCos2Tht = TMath::Cos(2*par[6]);
  Double_t lSinTht = TMath::Sin(par[6]);

  Double_t inSqrt = TMath::Abs((par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +((par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*lCos2Tht
			       +4*(x[0] - x[1])*(par[2] - par[5])*lSinTht);

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

Double_t phixpn(Double_t *x, Double_t *par){

  Double_t inSqrt = TMath::Abs(((par[0] - par[3])*(par[0] - par[3]) 
				-2*(x[0] - x[1])*(x[0] - x[1]) 
				+(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
				+(+(par[0] - par[3])*(par[0] - par[3]) 
				  +(par[1] - par[4])*(par[1] - par[4]) 
				  +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
				+4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6])));
  
  Double_t phix = ((+2*(par[0] - par[3])*(x[0] - x[1]) 
		    -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		    +TMath::Sqrt(2)*(par[1] - par[4])
		    *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = -TMath::ACos(phix);

  return phix;
}

Double_t phixpp(Double_t *x, Double_t *par){

  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       +4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));

  Double_t phix = ((+2*(par[0] - par[3])*(x[0] - x[1]) 
		    -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		    +TMath::Sqrt(2)*(par[1] - par[4])
		    *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = TMath::ACos(phix);

  return phix;
}

Double_t phixnn(Double_t *x, Double_t *par){

  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       + 4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));
  
  Double_t phix = (+(+2*(par[0] - par[3])*(x[0] - x[1]) 
		     -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		     +TMath::Sqrt(2)*(-par[1] + par[4])
		     *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = -TMath::ACos(phix);

  return phix;
}

Double_t phixnp(Double_t *x, Double_t *par){

  Double_t inSqrt = TMath::Abs(+(par[0] - par[3])*(par[0] - par[3]) 
			       -2*(x[0] - x[1])*(x[0] - x[1]) 
			       +(par[1] - par[4] + par[2] - par[5])*(par[1] - par[4] - par[2] + par[5]) 
			       +(+(par[0] - par[3])*(par[0] - par[3]) 
				 +(par[1] - par[4])*(par[1] - par[4]) 
				 +(par[2] - par[5])*(par[2] - par[5]))*TMath::Cos(2*par[6]) 
			       +4*(x[0] - x[1])*(par[2] - par[5])*TMath::Sin(par[6]));

  Double_t phix = (+(+2*(par[0] - par[3])*(x[0] - x[1]) 
		     -2*(par[0] - par[3])*(par[2] - par[5])*TMath::Sin(par[6]) 
		     +TMath::Sqrt(2)*(-par[1] + par[4])
		     *TMath::Sqrt(inSqrt))
		   /(2*(+(par[0] - par[3])*(par[0] - par[3]) 
			+(par[1] - par[4])*(par[1] - par[4]))*TMath::Cos(par[6])));

  phix = TMath::ACos(phix);

  return phix;
}

Double_t ypCenter(Double_t *x, Double_t *par){
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

Double_t phiypn(Double_t *x, Double_t *par){

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

Double_t phiypp(Double_t *x, Double_t *par){

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
 
Double_t ynCenter(Double_t *x, Double_t *par){

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
 

Double_t phiynn(Double_t *x, Double_t *par){

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


Double_t phiynp(Double_t *x, Double_t *par){

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

Double_t znCenter(Double_t *x, Double_t *par){
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

Double_t zpCenter(Double_t *x, Double_t *par){
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
AliMUONGeometryTransformer *ReAlign(const AliMUONGeometryTransformer * transformer, 
				    int rMod, TGeoCombiTrans deltaDetElemTransf[], Bool_t verbose)
{
  /////////////////////////////////////////////////////////////////////
  //   Takes the internal geometry module transformers, copies them
  // and gets the Detection Elements from them.
  // Takes misalignment parameters and applies these
  // to the local transform of the Detection Element
  // Obtains the global transform by multiplying the module transformer
  // transformation with the local transformation 
  // Applies the global transform to a new detection element
  // Adds the new detection element to a new module transformer
  // Adds the new module transformer to a new geometry transformer
  // Returns the new geometry transformer


  Int_t iDetElemId = 0;
  Int_t iDetElemNumber = 0;
  Int_t iDetElemIndex = 0;
  Int_t iCh = 0;

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
    if (iMt==rMod) {
      newModuleTransform = new TGeoCombiTrans(moduleTransform);
    } else {
      newModuleTransform = new TGeoCombiTrans(moduleTransform);
    }
    newModuleTransformer->SetTransformation(*newModuleTransform);
    
    // For the selected chamber add misalign module
    if (iMt==rMod) {
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
      iDetElemNumber = iDetElemId%100;
      iCh = iDetElemId/100 -1;
      if(iMt==rMod){
	if (iCh<4) {
	  iDetElemIndex = iDetElemId;
	} else {
	  if ((iDetElemNumber > (fgNDetElemCh[iCh]-2)/4) &&
	      (iDetElemNumber < fgNDetElemCh[iCh]-(fgNDetElemCh[iCh]-2)/4)) {
	    iDetElemIndex = (+fgNDetElemCh[iCh] 
			     -(1+(fgNDetElemCh[iCh]-2)/4) 
			     -iDetElemNumber);
	  } else {
	    iDetElemIndex = (+fgNDetElemCh[iCh] 
			     -fgNDetElemCh[iCh]/2
			     -((1+(fgNDetElemCh[iCh]-2)/4) 
			       -TMath::Min(iDetElemNumber,
					   TMath::Abs(iDetElemNumber-fgNDetElemCh[iCh]))));
	  }
	}
	deltaLocalTransform = new TGeoCombiTrans(deltaDetElemTransf[iDetElemIndex]);       
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
