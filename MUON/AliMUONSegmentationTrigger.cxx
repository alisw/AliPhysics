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
Revision 1.1.2.1  2000/06/09 21:47:24  morsch
Code from AliMUONSegResTrigger.cxx

*/

/*
old Log:
AliMUONSegResTrigger.cxx,v $
Revision 1.1.2.3  2000/04/26 12:32:39  morsch
Mods by P. Crochet:
- adapted to the new Trigger chamber geometry
- method SetZScale removed

Revision 1.1.2.2  2000/02/21 16:13:33  morsch
Full cluster simulation activated by uncommenting corresponding lines in IntXY()

Revision 1.1.2.1  2000/02/17 14:32:40  morsch
Draft version from P. Crochet

*/

#include "AliMUONSegmentationTrigger.h"
#include <TMath.h>
#include <TRandom.h>
#include <TArc.h>
#include "AliMUONChamber.h"
#include <iostream.h>
 
ClassImp(AliMUONSegmentationTrigger)

void AliMUONSegmentationTrigger::Init(AliMUONChamber* Chamber)
{
  // initialize Module geometry
  cout << "Initialize Trigger Chamber Module Geometry " << "\n";    

  Float_t zPos=Chamber->Z();
  Float_t z1Pos=1603.5;
  fZscale = zPos/z1Pos;

  static Int_t nModule=126;
  fgNmodule = nModule; 
// conv : line-column (line : from top to bottom, column : from left to right)
  static Int_t num[126]=
  {11,12,13,14,15,16,17,         // right side of the chamber
   21,22,23,24,25,26,27,
   31,32,33,34,35,36,37,
   41,42,43,44,45,46,47,
   51,52,53,54,55,56,57,
   61,62,63,64,65,66,67,
   71,72,73,74,75,76,77,
   81,82,83,84,85,86,87,
   91,92,93,94,95,96,97,   
   -11,-12,-13,-14,-15,-16,-17,     // left side of the chamber
   -21,-22,-23,-24,-25,-26,-27,
   -31,-32,-33,-34,-35,-36,-37,
   -41,-42,-43,-44,-45,-46,-47,
   -51,-52,-53,-54,-55,-56,-57,
   -61,-62,-63,-64,-65,-66,-67,
   -71,-72,-73,-74,-75,-76,-77,
   -81,-82,-83,-84,-85,-86,-87,
   -91,-92,-93,-94,-95,-96,-97};
  fgNum     = num; 

  static Int_t nStripX[126]=
  {16,16,16,16,16,16,16,  // right side of the chamber 
   32,32,32,32,32,32,16,
   32,32,32,32,32,32,16,
   48,64,64,32,32,32,16,
   0,64,64,32,32,32,16,
   48,64,64,32,32,32,16,
   32,32,32,32,32,32,16,
   32,32,32,32,32,32,16,
   16,16,16,16,16,16,16,  // left side of the chamber
   16,16,16,16,16,16,16,
   32,32,32,32,32,32,16,
   32,32,32,32,32,32,16,
   48,64,64,32,32,32,16,
   0,64,64,32,32,32,16,
   48,64,64,32,32,32,16,
   32,32,32,32,32,32,16,
   32,32,32,32,32,32,16,
   16,16,16,16,16,16,16};
  fgNstripx = nStripX;
  
  static Int_t nStripY[126]=
  { 8, 8, 8, 8, 8, 8,16,  // right side of the chamber
    8, 8, 8, 8, 8, 8,16,
    16,16,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    0, 8,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    8, 8, 8, 8, 8, 8,16,
    8, 8, 8, 8, 8, 8,16,  // left side of the chamber
    8, 8, 8, 8, 8, 8,16,  // right side of the chamber
    8, 8, 8, 8, 8, 8,16,
    16,16,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    0, 8,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    16,16,16,16,16, 8,16,
    8, 8, 8, 8, 8, 8,16,
    8, 8, 8, 8, 8, 8,16};
  fgNstripy = nStripY;
  
  static Float_t xCmin[126]=
  {0.,34.,68.,102.,136.,170.,204., // right
   0.,34.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   0.,51.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   0.,34.,68.,102.,136.,170.,204.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272., //left
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   0.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.,
   -34.,-68.,-102.,-136.,-170.,-204.,-272.};
  fgXcmin   = xCmin;
  
  static Float_t xCmax[126]=
  {34.,68.,102.,136.,170.,204.,272., //right
   34.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272.,
   0.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272.,
   34.,68.,102.,136.,170.,204.,272., 
   0.,-34.,-68.,-102.,-136.,-170.,-204., // left
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-51.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.,
   0.,-34.,-68.,-102.,-136.,-170.,-204.};
  fgXcmax   = xCmax;

  static Float_t yCmin[126];
  static Float_t yCmax[126];
  Float_t y1Cmin[126];
  Float_t y1Cmax[126];

  Float_t dz=7.2;
  Float_t z1PosPlus=z1Pos+dz/2.;
  Float_t z1PosMinus=z1Pos-dz/2.;

  Float_t z1pm=z1PosPlus/z1PosMinus;
  Float_t z1mp=z1PosMinus/z1PosPlus;

  cout << " fZscale = " << fZscale << "\n";

// calculate yCmin and yCmax 
  for (Int_t i=62; i>=0; i--) {
    Int_t j=ModuleNumber(-num[i]);  // i == right, j == left 
    if (Int_t(num[i]/10)==5) {  // start with middle chamber
      if (num[i]==51) {         // special case (empty module)
	yCmin[i]=yCmax[i]=yCmin[j]=yCmax[j]=0.;
      } else {
	y1Cmin[i]=y1Cmin[j]=-34;
	y1Cmax[i]=y1Cmax[j]=34;
	yCmin[i]=yCmin[j]=-34.;
	yCmax[i]=yCmax[j]=34.;
      }
    } else if (Int_t(num[i]/10)==4) { // up
      if (num[i]!=41) {       
	y1Cmin[i]=y1Cmax[i+7]*z1pm;
	y1Cmax[i]=y1Cmin[i]+68.;
	yCmin[i]=y1Cmin[i];
	yCmax[i]=yCmin[i]+68.;

	y1Cmin[j]=y1Cmax[j+7]*z1mp;
	y1Cmax[j]=y1Cmin[j]+68.;
	yCmin[j]=y1Cmin[j];
	yCmax[j]=yCmin[j]+68.;
      } else { 
	y1Cmin[i]=y1Cmin[ModuleNumber(42)]+17;
	y1Cmax[i]=y1Cmin[i]+51.;
	yCmin[i]=y1Cmin[i];
	yCmax[i]=yCmin[i]+51.;

	y1Cmin[j]=y1Cmin[ModuleNumber(-42)]+17;
	y1Cmax[j]=y1Cmin[j]+51.;
	yCmin[j]=y1Cmin[j];
	yCmax[j]=yCmin[j]+51.;
      }
    } else if (Int_t(num[i]/10)==3) { 
      y1Cmin[i]=y1Cmax[i+7]*z1mp;
      y1Cmax[i]=y1Cmin[i]+68.;
      yCmin[i]=y1Cmin[i];
      yCmax[i]=yCmin[i]+68.;

      y1Cmin[j]=y1Cmax[j+7]*z1pm;
      y1Cmax[j]=y1Cmin[j]+68.;
      yCmin[j]=y1Cmin[j];
      yCmax[j]=yCmin[j]+68.;
    } else if (Int_t(num[i]/10)==2) {
      y1Cmin[i]=y1Cmax[i+7]*z1pm;
      y1Cmax[i]=y1Cmin[i]+68.;
      yCmin[i]=y1Cmin[i];
      yCmax[i]=yCmin[i]+68.;

      y1Cmin[j]=y1Cmax[j+7]*z1mp;
      y1Cmax[j]=y1Cmin[j]+68.;
      yCmin[j]=y1Cmin[j];
      yCmax[j]=yCmin[j]+68.;
    } else if (Int_t(num[i]/10)==1) {
      y1Cmin[i]=y1Cmax[i+7]*z1mp;
      y1Cmax[i]=y1Cmin[i]+68.;
      yCmin[i]=y1Cmin[i];
      yCmax[i]=yCmin[i]+68.;

      y1Cmin[j]=y1Cmax[j+7]*z1pm;
      y1Cmax[j]=y1Cmin[j]+68.;
      yCmin[j]=y1Cmin[j];
      yCmax[j]=yCmin[j]+68.;
    }
  }

  for (Int_t i=0; i<63; i++) {      // second loop (fill lower part)
    Int_t j=ModuleNumber(-num[i]);  // i == right, j == left 
    if (TMath::Abs(Int_t(num[i]/10))==6) { 
      yCmin[i]=-yCmax[i-14];
      yCmax[i]=-yCmin[i-14];
      yCmin[j]=-yCmax[j-14];
      yCmax[j]=-yCmin[j-14];
    } else if (TMath::Abs(Int_t(num[i]/10))==7) { 
      yCmin[i]=-yCmax[i-28];
      yCmax[i]=-yCmin[i-28];
      yCmin[j]=-yCmax[j-28];
      yCmax[j]=-yCmin[j-28];
    } else if (TMath::Abs(Int_t(num[i]/10))==8) { 
      yCmin[i]=-yCmax[i-42];
      yCmax[i]=-yCmin[i-42];
      yCmin[j]=-yCmax[j-42];
      yCmax[j]=-yCmin[j-42];
    } else if (TMath::Abs(Int_t(num[i]/10))==9) { 
      yCmin[i]=-yCmax[i-56];
      yCmax[i]=-yCmin[i-56];
      yCmin[j]=-yCmax[j-56];
      yCmax[j]=-yCmin[j-56];
    } 
  }

  fgYcmin   = yCmin;
  fgYcmax   = yCmax;
  
  fNpx=124;
  fNpy=64;  

  cout << "---------------------------------------------------- \n";   

}

//------------------------------------------------------------------
Int_t AliMUONSegmentationTrigger::ModuleNumber(Int_t imodule){
// returns module number (from 0 to 126) corresponding to module imodule
  Int_t imod=0;
  for (Int_t i=0; i<fgNmodule; i++) {
    if (fgNum[i]==imodule) { 
      imod=i;
      break;
    }
  }
  return imod;
}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTrigger::StripSizeX(Int_t imodule){
// Returns x-strip size for given module imodule

  Int_t absimodule=TMath::Abs(imodule); 
  Int_t moduleNum=ModuleNumber(imodule);
  if (fgNum[absimodule]==51) {
    return 0; 
  } else {
    return TMath::Abs((fgYcmax[moduleNum]-fgYcmin[moduleNum])/
		      fgNstripx[moduleNum]);
  }
}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTrigger::StripSizeY(Int_t imodule){
// Returns y-strip size for given module imodule
        
  Int_t absimodule=TMath::Abs(imodule); 
  Int_t moduleNum=ModuleNumber(imodule);
  if (fgNum[absimodule]==51) {
    return 0;
  } else {
    return TMath::Abs((fgXcmax[moduleNum]-fgXcmin[moduleNum])/
		      fgNstripy[moduleNum]);
  }
}

//------------------------------------------------------------------   
void AliMUONSegmentationTrigger::SetHit(Float_t xhit, Float_t yhit)
{
    // Sets virtual hit position, needed for evaluating pad response 
    // outside the tracking program 
    
  fxhit=xhit;
  fyhit=yhit;
}
