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
#include "AliRICHParam.h"

ClassImp(AliRICHChamber)	
//______________________________________________________________________________
AliRICHChamber::AliRICHChamber(Int_t iChamber):TNamed()
{
//main ctor. Defines all geometry parameters for the given module.
// 7 6     ^      
// 5 4 3   |      
//   2 1   |
// <-----z y   . x     
  if(iChamber!=0){ //if iChamber=0 then special test config: chamber without any shift and rotation
    const Double_t kAngHor=19.5; //  horizontal angle between chambers  19.5 grad
    const Double_t kAngVer=20;   //  vertical angle between chambers    20   grad     
    const Double_t kAngCom=30;   //  common RICH rotation with respect to x axis  20   grad     
  
    RotY(90);//rotate around y since initial position is in XY plane
    fRad   .SetXYZ(490-2    ,0,0);       //position of the entrance to freon 1.5 cm of freon+0.5 cm of quartz window
    fCenter.SetXYZ(490      ,0,0);      //shift center along x by 490 cm
    fAnod  .SetXYZ(490+8-0.2,0,0);  //position of the center of apmlification gap (anod wires plane)
    fPc    .SetXYZ(490+8    ,0,0);        //position of the center of PC 
    switch(iChamber){
      case 1:        
        RotY(kAngHor);  RotZ(-kAngVer);   //right and down 
        break;      
      case 2:
                        RotZ(-kAngVer);   //down
        break;      
      case 3:
        RotY(kAngHor);                   //right 
        break;      
      case 4:          
        break;                           //no rotation
      case 5:
        RotY(-kAngHor);                  //left   
        break;      
      case 6:
                        RotZ(kAngVer);     //up
        break;      
      case 7:
        RotY(-kAngHor); RotZ(kAngVer);  //left and up 
        break;      
      default:
        Fatal("named ctor","Wrong chamber number %i, check CreateChamber ctor",iChamber);
    }//switch(iChamber)
    RotZ(kAngCom);     //apply common rotation  
  }//if(iChamber
  fName=Form("RICHc%i",iChamber);fTitle=Form("RICH chamber %i",iChamber);
  fpRotMatrix=new TRotMatrix("rot"+fName,"rot"+fName, Rot().ThetaX()*TMath::RadToDeg(), Rot().PhiX()*TMath::RadToDeg(),
                                                      Rot().ThetaY()*TMath::RadToDeg(), Rot().PhiY()*TMath::RadToDeg(),
                                                      Rot().ThetaZ()*TMath::RadToDeg(), Rot().PhiZ()*TMath::RadToDeg());
}//main ctor
//__________________________________________________________________________________________________
void AliRICHChamber::Print(Option_t *opt) const
{
// Debug printout
  TNamed::Print(opt);
  fRad.Print(opt);
  fCenter.Print(opt);
  fAnod.Print(opt);
  fPc.Print(opt);
}//Print()
//__________________________________________________________________________________________________
