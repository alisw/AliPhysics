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

//-------------------------------------------------------------------------
//      Lego generator in x - y - z bins
//    Uses geantino rays to check the material distributions and detector's
//    geometry
//    Author: A.Morsch 
//-------------------------------------------------------------------------

#include "AliLegoGeneratorXYZ.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliLegoGeneratorXYZ)


//___________________________________________
    

AliLegoGeneratorXYZ::AliLegoGeneratorXYZ()
{
// Default Constructor
    fDir1[0]=1.; fDir1[1]=0.; fDir1[2]=0.;
    fDir2[0]=0.; fDir2[1]=1.; fDir2[2]=0.;    
    fDir3[0]=0.; fDir3[1]=0.; fDir3[2]=1.;    
}

AliLegoGeneratorXYZ::AliLegoGeneratorXYZ(char* axis)
{
// Constructor
    if (!strcmp(axis,"x") || !strcmp(axis,"X")) 
    {
	fDir1[0]=0.; fDir1[1]=1.; fDir1[2]=0.;
	fDir2[0]=0.; fDir2[1]=0.; fDir2[2]=1.;    
	fDir3[0]=1.; fDir3[1]=0.; fDir3[2]=0.;    
    }
    else if (!strcmp(axis,"y") || !strcmp(axis,"Y")) 
    {
	fDir1[0]=1.; fDir1[1]=0.; fDir1[2]=0.;
	fDir2[0]=0.; fDir2[1]=0.; fDir2[2]=1;    
	fDir3[0]=0.; fDir3[1]=1.; fDir3[2]=0.;    
    }
    else if (!strcmp(axis,"z") || !strcmp(axis,"Z")) 
    {
	fDir1[0]=1.; fDir1[1]=0.; fDir1[2]=0.;
	fDir2[0]=0.; fDir2[1]=1.; fDir2[2]=0.;    
	fDir3[0]=0.; fDir3[1]=0.; fDir3[2]=1.;    
    }
    else 
    {
	fDir1[0]=1.; fDir1[1]=0.; fDir1[2]=0.;
	fDir2[0]=0.; fDir2[1]=1.; fDir2[2]=0.;    
	fDir3[0]=0.; fDir3[1]=0.; fDir3[2]=1.;    
    }
}


AliLegoGeneratorXYZ::AliLegoGeneratorXYZ(Int_t nc1, Float_t c1min,
					 Float_t c1max, Int_t nc2, 
					 Float_t c2min, Float_t c2max,
					 Float_t rmin, Float_t rmax, Float_t zmax) : 
    AliLegoGenerator(nc1, c1min, c1max, nc2, c2min, c2max,
		     rmin, rmax, zmax)
{
//  Constructor
    fDir1[0]=1.; fDir1[1]=0.; fDir1[2]=0.;
    fDir2[0]=0.; fDir2[1]=1.; fDir2[2]=0.;    
    fDir3[0]=0.; fDir3[1]=0.; fDir3[2]=1.;    
}


//___________________________________________
void AliLegoGeneratorXYZ::Generate()
{
// Create a geantino with kinematics corresponding to the current bins
// Here: Coor1 =  x 
//       Coor2 =  z
   
  //
  // Rootinos are 0
   const Int_t kMpart = 0;
   Float_t orig[3], pmom[3];
   
   // Prepare for next step
   if(fCoor1Bin>=fNCoor1-1)
     if(fCoor2Bin>=fNCoor2-1) {
       AliWarning("End of Lego Generation");
       return;
     } else { 
       fCoor2Bin++;
       AliDebug(1, Form("Generating rays in Coordinate 2 bin:%d",fCoor2Bin));
       fCoor1Bin=0;
     } else fCoor1Bin++;

   
   fCurCoor1 = (fCoor1Min+(fCoor1Bin+0.5)*(fCoor1Max-fCoor1Min)/fNCoor1);
   fCurCoor2 = (fCoor2Min+(fCoor2Bin+0.5)*(fCoor2Max-fCoor2Min)/fNCoor2);

// Origin and direction
   Int_t i;
   for (i=0; i<3; i++) {
       pmom[i]=fDir3[i];
       orig[i]=fCurCoor1*fDir1[i]+fCurCoor2*fDir2[i];
   }
   
   Float_t dalicz = 3000;
   if (fRadMin > 0) {
       Float_t t = PropagateCylinder(orig,pmom,fRadMin,dalicz);
       orig[0] = pmom[0]*t;
       orig[1] = pmom[1]*t;
       orig[2] = pmom[2]*t;
       if (TMath::Abs(orig[2]) > fZMax) return;
   }
   
   Float_t polar[3]={0.,0.,0.};
   Int_t ntr;
   gAlice->GetMCApp()->PushTrack(1, -1, kMpart, pmom, orig, polar, 0, kPPrimary, ntr);
   
}




