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
Revision 1.1  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

*/

#include "AliMagFC.h"
#include <stdlib.h>

ClassImp(AliMagFC)

//________________________________________
AliMagFC::AliMagFC(const char *name, const char *title, const Int_t integ, 
		   const Int_t map, const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
{
  // 
  // Standard constructor
  //
  printf("Constant Field %s created: map= %d, factor= %f\n",fName.Data(),map,
	 factor);
  fType = kConst;
}

//________________________________________
void AliMagFC::Field(Float_t *x, Float_t *b)
{
  //
  // Method to return the field in a point
  //
  b[0]=b[1]=b[2]=0;
  if(fMap==1) {
    if(TMath::Abs(x[2])<700 && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
      b[2]=2;
    } else {
      if ( 725 <= x[2] && x[2] <= 1225 ) {
	Float_t dz = TMath::Abs(975-x[2])*0.01;
	b[0]=(1-0.1*dz*dz)*7;
      }
      else {
//This is the ZDC part
	Float_t rad2=x[0]*x[0]+x[1]*x[1];
	if(rad2<kD2RA2) {
	  if(x[2]>kD2BEG) {
	    
//    Separator Dipole D2
	    if(x[2]<kD2END) b[1]=kFDIP;
	  } else if(x[2]>kD1BEG) {
	    
//    Separator Dipole D1
	    if(x[2]<kD1END) b[1]=-kFDIP;
	  }
	  if(rad2<kCORRA2) {

//    First quadrupole of inner triplet de-focussing in x-direction
//    Inner triplet
	    if(x[2]>kZ4BEG) {
	      if(x[2]<kZ4END) {
	      
//    2430 <-> 3060
		b[0]=-kG1*x[1];
		b[1]=-kG1*x[0];
	      }
	    } else if(x[2]>kZ3BEG) {
	      if(x[2]<kZ3END) {

//    1530 <-> 2080
		b[0]=kG1*x[1];
		b[1]=kG1*x[0];
	      }
	    } else if(x[2]>kZ2BEG) {
	      if(x[2]<kZ2END) {
	      
//    890 <-> 1430
		b[0]=kG1*x[1];
		b[1]=kG1*x[0];
	      }
	    } else if(x[2]>kZ1BEG) {
	      if(x[2]<kZ1END) {

//    0 <->  630
		b[0]=-kG1*x[1];
		b[1]=-kG1*x[0];
	      }
	    } else if(x[2]>kCORBEG) {
	      if(x[2]<kCOREND) {
//    Corrector dipole (because of dimuon arm)
		b[0]=kFCORN;
	      }
	    }
	  }
	}
      }
    }
    if(fFactor!=1) {
      b[0]*=fFactor;
      b[1]*=fFactor;
      b[2]*=fFactor;
    }
  } else {
    printf("Invalid field map for constant field %d\n",fMap);
    exit(1);
  }
}

