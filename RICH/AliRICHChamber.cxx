//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICHChamber.h"
#include <TRotMatrix.h>
#include <AliLog.h>

ClassImp(AliRICHChamber)	
//______________________________________________________________________________
AliRICHChamber::AliRICHChamber(Int_t iChamber):TNamed()
{
//main ctor. Defines all geometry parameters for the given module.
// 7 6     ^      
// 5 4 3   |      
//   2 1   |
// <-----z y   . x     
//  horizontal angle between chambers  19.5 grad
//  vertical angle between chambers    20   grad     
  RotY(90);//rotate around y
  fCenterX3.SetXYZ(490,0,0);fPcX3.SetXYZ(490+8-0.4,0,0);fRadX3.SetXYZ(490-2,0,0); //shift center along x by 490 cm
  
  switch(iChamber){
    case 0:                    //special test beam configuration without rotation.
      break;
    case 1:        
      RotY(19.5); RotZ(-20);   //right and down 
      break;      
    case 2:
      RotZ(-20);              //down
      break;      
    case 3:
      RotY(19.5);             //right 
      break;      
    case 4:          
      break;                  //no rotation
    case 5:
      RotY(-19.5);            //left   
      break;      
    case 6:
      RotZ(20);               //up
      break;      
    case 7:
      RotY(-19.5); RotZ(20);  //left and up 
      break;      
    default:
      Fatal("named ctor","Wrong chamber number %i, check CreateChamber ctor",iChamber);
  }//switch(iModuleN)
  fName=Form("RICHc%i",iChamber);fTitle=Form("RICH chamber %i",iChamber);
  RotZ(30);     //apply common rotation  
  fpRotMatrix=new TRotMatrix("rot"+fName,"rot"+fName, Rot().ThetaX()*TMath::RadToDeg(), Rot().PhiX()*TMath::RadToDeg(),
                                                      Rot().ThetaY()*TMath::RadToDeg(), Rot().PhiY()*TMath::RadToDeg(),
                                                      Rot().ThetaZ()*TMath::RadToDeg(), Rot().PhiZ()*TMath::RadToDeg());
}//main ctor
//__________________________________________________________________________________________________
void AliRICHChamber::Print(Option_t *opt) const
{
// Debug printout
  ToAliInfo(fCenterX3.Print(opt));
}//Print()
//__________________________________________________________________________________________________
